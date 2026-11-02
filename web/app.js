// @ts-check
// Page logic for Show Clusters. Python (worker.js) does the science; this file handles
// the form, the Gaia requests, the picture sliders and downloads.

const GAVO_TAP = "https://dc.g-vo.org/tap/sync";
const MAX_ROWS = 300000; // keep in step with showclusters.MAX_ROWS
const MAX_UPLOAD_BYTES = 50 * 1024 * 1024;
const GAVO_TIMEOUT_MS = 60000;

/** @template {HTMLElement} T @param {string} id @returns {T} */
function $(id) {
  const value = document.getElementById(id);
  if (!value) throw new Error(`Missing element: ${id}`);
  return /** @type {T} */ (value);
}

const form = /** @type {HTMLFormElement} */ ($("field"));
const objectInput = /** @type {HTMLInputElement} */ ($("object"));
const raInput = /** @type {HTMLInputElement} */ ($("ra"));
const decInput = /** @type {HTMLInputElement} */ ($("dec"));
const fovInput = /** @type {HTMLInputElement} */ ($("fov"));
const maxmagInput = /** @type {HTMLInputElement} */ ($("maxmag"));
const showButton = /** @type {HTMLButtonElement} */ ($("show-stars"));
const status = $("status");
const errorBox = $("error");
const retry = $("retry");
const starsPanel = $("stars-panel");
const hrPanel = $("hr-panel");
const canvas = /** @type {HTMLCanvasElement} */ ($("sky"));
const magtype = /** @type {HTMLSelectElement} */ ($("magtype"));
const hrButton = /** @type {HTMLButtonElement} */ ($("hr-button"));
const hrImage = /** @type {HTMLImageElement} */ ($("hr"));
const hrDownload = /** @type {HTMLAnchorElement} */ ($("hr-download"));
const skyDownload = /** @type {HTMLAnchorElement} */ ($("sky-download"));
const adql = $("adql");
const copyButton = /** @type {HTMLButtonElement} */ ($("copy-adql"));
const fileInput = /** @type {HTMLInputElement} */ ($("file"));
const useFile = /** @type {HTMLButtonElement} */ ($("use-file"));
const uploadDetails = /** @type {HTMLDetailsElement} */ ($("upload"));
const sliders = {
  colsat: /** @type {HTMLInputElement} */ ($("colsat")),
  gamma: /** @type {HTMLInputElement} */ ($("gamma")),
  mini: /** @type {HTMLInputElement} */ ($("mini")),
  maxi: /** @type {HTMLInputElement} */ ($("maxi")),
};

// ---- Worker calls ------------------------------------------------------------------

/** @type {Worker} */
let worker;
let ready = false;
let busy = false;
let nextId = 0;
/** @type {Map<number, {resolve: (v: any) => void, reject: (e: Error) => void}>} */
const pending = new Map();

class StudentError extends Error {}

/** @param {string} action @param {Record<string, unknown>} args @param {Transferable[]} [transfer] */
function call(action, args, transfer = []) {
  const id = ++nextId;
  return new Promise((resolve, reject) => {
    pending.set(id, { resolve, reject });
    worker.postMessage({ id, action, args }, transfer);
  });
}

function startWorker() {
  ready = false;
  retry.hidden = true;
  $("load-details").hidden = true;
  worker = new Worker(new URL("worker.js", import.meta.url), { type: "module" });
  worker.onmessage = ({ data }) => {
    if (data.type === "ready") {
      ready = true;
      setStatus("Ready. Type an object name such as m13 and press Enter.");
      refresh();
      return;
    }
    if (data.type === "fatal") {
      showError("Python could not load. Check your connection and click Reload tool.");
      $("load-error").textContent = `${data.stage}: ${String(data.message).slice(0, 2000)}`;
      $("load-details").hidden = false;
      retry.hidden = false;
      return;
    }
    const waiting = pending.get(data.id);
    if (!waiting) return;
    pending.delete(data.id);
    if (data.ok) waiting.resolve(data.result);
    else waiting.reject(data.expected ? new StudentError(data.message) : new Error(data.message));
  };
  worker.onerror = () => {
    for (const waiting of pending.values()) waiting.reject(new Error("The Python worker stopped."));
    pending.clear();
    ready = false;
    busy = false;
    showError("The Python worker stopped (possibly out of memory). Click Reload tool, then try a smaller field.");
    retry.hidden = false;
    refresh();
  };
}

// ---- Status helpers ----------------------------------------------------------------

/** @param {string} text */
function setStatus(text) {
  status.textContent = text;
}
/** @param {string} text */
function showError(text) {
  errorBox.textContent = text;
  errorBox.hidden = false;
}
function clearError() {
  errorBox.hidden = true;
  errorBox.textContent = "";
}
function refresh() {
  showButton.disabled = !ready || busy;
  useFile.disabled = !ready || busy || !fileInput.files?.length;
  hrButton.disabled = !ready || busy;
  copyButton.disabled = !adql.dataset.query;
}
/** @param {unknown} error @param {string} fallback */
function report(error, fallback) {
  if (error instanceof StudentError) showError(error.message);
  else showError(`${fallback} (${error instanceof Error ? error.message : String(error)})`);
}

function fieldArgs() {
  return {
    ra: raInput.value.trim(),
    dec: decInput.value.trim(),
    fov: Number(fovInput.value),
    maxmag: Number(maxmagInput.value),
  };
}

// ---- Object lookup and query preview ------------------------------------------------

let objectName = "";

async function lookup() {
  const name = objectInput.value.trim();
  if (!name) return true;
  try {
    const found = await call("lookup", { name });
    raInput.value = found.ra;
    decInput.value = found.dec;
    objectName = found.name;
    clearError();
    await preview();
    return true;
  } catch (error) {
    report(error, "Lookup failed");
    return false;
  }
}

async function preview() {
  if (!ready) return;
  try {
    const plan = await call("plan", fieldArgs());
    adql.textContent = plan.esa;
    adql.dataset.query = plan.esa;
  } catch (error) {
    adql.textContent = error instanceof StudentError ? error.message : "The query could not be prepared.";
    delete adql.dataset.query;
  }
  refresh();
}

// ---- Getting star data: GAVO, then a saved copy, or an uploaded file ----------------

/** @type {Promise<Record<string, any>> | null} */
let snapshotIndex = null;

/** @param {string} query */
async function fetchGavo(query) {
  const body = new URLSearchParams({ REQUEST: "doQuery", LANG: "ADQL", FORMAT: "csv", MAXREC: String(MAX_ROWS), QUERY: query });
  const response = await fetch(GAVO_TAP, { method: "POST", body, signal: AbortSignal.timeout(GAVO_TIMEOUT_MS) });
  const type = response.headers.get("content-type") || "";
  if (!response.ok || !type.includes("csv")) throw new Error(`GAVO answered HTTP ${response.status}`);
  return new Uint8Array(await response.arrayBuffer());
}

/** @param {string} name @param {{ra: string, dec: string, fov: number, maxmag: number}} field */
async function savedCopy(name, field) {
  snapshotIndex ??= fetch("data/snapshots/index.json").then((r) => (r.ok ? r.json() : {})).catch(() => ({}));
  const entry = (await snapshotIndex)[name];
  if (!entry || field.ra !== entry.ra || field.dec !== entry.dec || field.fov > entry.fov || field.maxmag > entry.maxmag) {
    return null;
  }
  const response = await fetch(`data/snapshots/${entry.file}`);
  if (!response.ok) return null;
  return { bytes: new Uint8Array(await response.arrayBuffer()), date: entry.date };
}

/**
 * @param {Uint8Array} bytes @param {string} filename @param {string} source @param {number} maxRows
 */
async function showStars(bytes, filename, source, maxRows) {
  const field = fieldArgs();
  const title = objectName || `RA ${field.ra}, Dec ${field.dec}`;
  setStatus("Placing the stars and drawing the picture…");
  const result = await call(
    "showStars",
    { bytes, filename, source, title, maxRows, ...field },
  );
  displayStars(result, source);
  setStatus(`Plotted ${result.count.toLocaleString()} stars. Adjust the sliders, or click HR Diagram.`);
}

async function onShowStars() {
  if (!ready || busy) return;
  busy = true;
  clearError();
  refresh();
  try {
    if (!(await lookup())) return;
    const field = fieldArgs();
    const plan = await call("plan", field);
    setStatus("Downloading star data from Gaia (via GAVO)…");
    let bytes;
    let reason = "";
    try {
      bytes = await fetchGavo(plan.gavo);
    } catch (error) {
      reason = error instanceof Error && error.name === "TimeoutError" ? "it did not answer in time" : String(error instanceof Error ? error.message : error);
    }
    if (bytes) {
      await showStars(bytes, "gavo.csv", "Gaia DR3 via GAVO", MAX_ROWS);
      return;
    }
    const saved = objectName ? await savedCopy(objectName, field) : null;
    if (saved) {
      await showStars(saved.bytes, "saved.csv", `Gaia DR3 saved copy from ESA (${saved.date})`, 0);
      setStatus(`The live data source was unavailable (${reason}), so a saved copy was used.`);
      return;
    }
    uploadDetails.open = true;
    showError(
      `The live Gaia data source is unavailable (${reason}), and no saved copy covers this field. ` +
      "Try again in a few minutes, or use a file from the ESA Gaia archive (below).",
    );
    setStatus("No star data.");
  } catch (error) {
    report(error, "Show Stars failed");
    setStatus("No star data.");
  } finally {
    busy = false;
    refresh();
  }
}

async function onUseFile() {
  const file = fileInput.files?.[0];
  if (!file || !ready || busy) return;
  clearError();
  if (file.size > MAX_UPLOAD_BYTES) {
    showError("The file is larger than 50 MB. Use a smaller field or magnitude limit.");
    return;
  }
  busy = true;
  refresh();
  try {
    if (!(await lookup())) return;
    await showStars(new Uint8Array(await file.arrayBuffer()), file.name, `Gaia DR3 from uploaded ESA file (${file.name})`, 0);
  } catch (error) {
    report(error, "This file could not be used");
  } finally {
    busy = false;
    refresh();
  }
}

// ---- The picture: slidershow from ass.py, in JavaScript ------------------------------

/** @type {{base: Float32Array, size: number} | null} */
let picture = null;

/** @param {any} result @param {string} source */
function displayStars(result, source) {
  $("count").textContent = result.count.toLocaleString();
  $("brightest").textContent = result.brightest === null ? "—" : result.brightest.toFixed(2);
  $("source").textContent = source;
  const notes = $("notes");
  notes.replaceChildren(...result.notes.map((/** @type {string} */ text) => Object.assign(document.createElement("li"), { textContent: text })));
  const bytes = /** @type {Uint8Array} */ (result.image);
  picture = { base: new Float32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4), size: result.size };
  canvas.width = canvas.height = result.size;
  // Slider ranges and starting values as in slidershow (sqrt scale for min/max).
  const low = 0.8 * Math.sqrt(Math.max(result.pct_lo, 0));
  for (const slider of [sliders.mini, sliders.maxi]) {
    slider.min = String(low);
    slider.max = "1";
  }
  sliders.mini.value = String(Math.sqrt(Math.max(result.pct_lo, 0)));
  sliders.maxi.value = String(Math.sqrt(result.pct_hi));
  sliders.colsat.value = "5";
  sliders.gamma.value = "2";
  magtype.value = "Apparent"; // ass.py resets the magnitude type on Show Stars
  hrImage.hidden = true;
  hrDownload.hidden = true;
  starsPanel.hidden = false;
  hrPanel.hidden = false;
  draw();
  skyDownload.download = `${objectName || "field"}.png`;
}

function draw() {
  if (!picture) return;
  const { base, size } = picture;
  const colsat = Number(sliders.colsat.value);
  const gamma = 3 ** Number(sliders.gamma.value) - 2;
  const mini = Number(sliders.mini.value) ** 2;
  const maxi = Number(sliders.maxi.value) ** 2;
  for (const [key, value] of Object.entries({ colsat, gamma, mini, maxi })) {
    $(`${key}-value`).textContent = value.toFixed(key === "gamma" || key === "colsat" ? 2 : 4);
  }
  const pixels = size * size;
  const stretched = new Float32Array(pixels * 3);
  const lum = new Float32Array(pixels);
  const stretch = Math.abs(gamma) < 1e-9
    ? (/** @type {number} */ x) => x
    : (/** @type {number} */ x) => Math.asinh(x * gamma) / Math.asinh(gamma);
  let maxL = 0;
  if (mini <= maxi) {
    const span = maxi - mini || 1e-12;
    for (let p = 0; p < pixels; p++) {
      let sum = 0;
      for (let c = 0; c < 3; c++) {
        const v = stretch(Math.min(Math.max((base[p * 3 + c] - mini) / span, 0), 1));
        stretched[p * 3 + c] = v;
        sum += v;
      }
      lum[p] = sum / 3;
      if (lum[p] > maxL) maxL = lum[p];
    }
  } else {
    const span = mini - maxi;
    const denom = Math.log(1 + gamma);
    for (let i = 0; i < pixels * 3; i++) {
      stretched[i] = Math.log(2 - Math.min(Math.max((base[i] - maxi) / span, 0), 1) * gamma) / denom;
    }
  }
  const image = new ImageData(size, size);
  const out = image.data;
  for (let p = 0; p < pixels; p++) {
    const L = lum[p];
    const weight = mini <= maxi && maxL > 0 ? colsat * (L / maxL) : 0;
    for (let c = 0; c < 3; c++) {
      let v = stretched[p * 3 + c];
      if (mini <= maxi) v = L + weight * (v - L);
      out[p * 4 + c] = Math.round(255 * Math.min(Math.max(Number.isFinite(v) ? v : 0, 0), 1));
    }
    out[p * 4 + 3] = 255;
  }
  /** @type {CanvasRenderingContext2D} */ (canvas.getContext("2d")).putImageData(image, 0, 0);
}

let drawQueued = false;
function queueDraw() {
  if (drawQueued) return;
  drawQueued = true;
  requestAnimationFrame(() => {
    drawQueued = false;
    draw();
  });
}

// ---- HR diagram ----------------------------------------------------------------------

/** @type {string | null} */
let hrUrl = null;

async function onHrDiagram() {
  if (!ready || busy) return;
  busy = true;
  clearError();
  hrImage.hidden = true; // never show a diagram that doesn't match the current choice
  hrDownload.hidden = true;
  refresh();
  setStatus("Drawing the HR diagram…");
  try {
    const png = /** @type {Uint8Array} */ (await call("hrDiagram", { magnitudeType: magtype.value }));
    if (hrUrl) URL.revokeObjectURL(hrUrl);
    hrUrl = URL.createObjectURL(new Blob([png], { type: "image/png" }));
    hrImage.src = hrUrl;
    hrImage.hidden = false;
    hrDownload.href = hrUrl;
    hrDownload.download = `${objectName || "field"}_HR_${magtype.value}.png`;
    hrDownload.hidden = false;
    setStatus("HR diagram ready. Use Download HR diagram to save it.");
  } catch (error) {
    report(error, "The HR diagram could not be drawn");
  } finally {
    busy = false;
    refresh();
  }
}

// ---- Wiring --------------------------------------------------------------------------

objectInput.addEventListener("keydown", (event) => {
  if (event.key === "Enter" || event.key === "Tab") {
    if (event.key === "Enter") event.preventDefault();
    void lookup();
  }
});
objectInput.addEventListener("input", () => {
  objectName = "";
});
for (const input of [raInput, decInput]) {
  input.addEventListener("input", () => {
    objectName = "";
    void preview();
  });
}
for (const input of [fovInput, maxmagInput]) input.addEventListener("change", () => void preview());
form.addEventListener("submit", (event) => {
  event.preventDefault();
  void onShowStars();
});
for (const slider of Object.values(sliders)) slider.addEventListener("input", queueDraw);
hrButton.addEventListener("click", () => void onHrDiagram());
fileInput.addEventListener("change", refresh);
useFile.addEventListener("click", () => void onUseFile());
copyButton.addEventListener("click", async () => {
  const query = adql.dataset.query;
  if (!query) return;
  try {
    await navigator.clipboard.writeText(query);
    copyButton.textContent = "Copied";
  } catch {
    copyButton.textContent = "Select the text and copy it";
  }
  setTimeout(() => (copyButton.textContent = "Copy query"), 2500);
});
skyDownload.addEventListener("click", (event) => {
  if (!picture) return;
  event.preventDefault();
  canvas.toBlob((blob) => {
    if (!blob) return;
    const url = URL.createObjectURL(blob);
    const link = Object.assign(document.createElement("a"), { href: url, download: skyDownload.download });
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
  }, "image/png");
});
retry.addEventListener("click", () => {
  worker.terminate();
  clearError();
  setStatus("Reloading Python…");
  startWorker();
});

startWorker();
refresh();
