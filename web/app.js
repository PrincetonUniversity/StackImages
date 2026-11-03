// @ts-check
// Page logic for Show Clusters. Python (worker.js) does the science; this file handles
// the form, the Gaia requests, the picture sliders and downloads.

const GAVO_TAP = "https://dc.g-vo.org/tap/sync";
const MAX_ROWS = 300000; // keep in step with showclusters.MAX_ROWS
const MAX_UPLOAD_BYTES = 50 * 1024 * 1024;
const GAVO_TIMEOUT_MS = 60000; // per request; GAVO itself gives up after ~20 s
const GAVO_TOTAL_MS = 180000; // whole field, including splitting
const GAVO_PARALLEL = 4;

/** @template {HTMLElement} T @param {string} id @returns {T} */
function $(id) {
  const value = document.getElementById(id);
  if (!value) throw new Error(`Missing element: ${id}`);
  return /** @type {T} */ (value);
}

const form = /** @type {HTMLFormElement} */ ($("field"));
const objectInput = /** @type {HTMLInputElement} */ ($("object"));
const lookupButton = /** @type {HTMLButtonElement} */ ($("lookup"));
const raInput = /** @type {HTMLInputElement} */ ($("ra"));
const decInput = /** @type {HTMLInputElement} */ ($("dec"));
const fovInput = /** @type {HTMLInputElement} */ ($("fov"));
const maxmagInput = /** @type {HTMLInputElement} */ ($("maxmag"));
const showButton = /** @type {HTMLButtonElement} */ ($("show-stars"));
const status = $("status");
const errorBox = $("error");
const retry = $("retry");
const starsContent = $("stars-content");
const starsEmpty = $("stars-empty");
const canvas = /** @type {HTMLCanvasElement} */ ($("sky"));
const magtype = /** @type {HTMLSelectElement} */ ($("magtype"));
const hrButton = /** @type {HTMLButtonElement} */ ($("hr-button"));
const resetButton = /** @type {HTMLButtonElement} */ ($("reset-sliders"));
const fullscreenButton = /** @type {HTMLButtonElement} */ ($("fullscreen"));
const exitButton = /** @type {HTMLButtonElement} */ ($("exit-fullscreen"));
const viewer = /** @type {HTMLElement} */ (document.querySelector(".viewer"));
const hrImage = /** @type {HTMLImageElement} */ ($("hr"));
const hrDownload = /** @type {HTMLAnchorElement} */ ($("hr-download"));
const hrDownloads = $("hr-downloads");
const hrCredit = $("hr-credit");
const dataDownload = /** @type {HTMLAnchorElement} */ ($("data-download"));
const hrCsvDownload = /** @type {HTMLAnchorElement} */ ($("hr-csv-download"));
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
      setStatus("Ready. Type an object name such as m13, then click Look up (or press Enter).");
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
  lookupButton.disabled = !ready || busy;
  useFile.disabled = !ready || busy || !fileInput.files?.length;
  hrButton.disabled = !ready || busy || !picture;
  resetButton.disabled = !picture;
  fullscreenButton.disabled = !picture;
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
    setStatus(`Found ${found.name}: RA ${found.ra}, Dec ${found.dec}.`);
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
  if (response.ok && type.includes("csv")) return new Uint8Array(await response.arrayBuffer());
  // GAVO's quick endpoint stops queries after ~20 s and says so in its error document.
  const text = await response.text().catch(() => "");
  throw new GavoError(`GAVO answered HTTP ${response.status}`, /timed out/i.test(text));
}

class GavoError extends Error {
  /** @param {string} message @param {boolean} timedOut */
  constructor(message, timedOut) {
    super(message);
    this.timedOut = timedOut;
  }
}

class TooDenseError extends Error {}

/**
 * Fetch a field as one or more pieces (large fields are split so each answers quickly);
 * a piece that times out is split again. Returns one result per piece.
 * @param {{query: string, radius: number}[]} tiles
 */
async function fetchGavoTiles(tiles) {
  const queue = [...tiles];
  /** @type {Uint8Array<ArrayBuffer>[]} */
  const parts = [];
  let total = tiles.length;
  const deadline = Date.now() + GAVO_TOTAL_MS;
  const report = () => setStatus(total > 1
    ? `Downloading star data from Gaia (via GAVO)… ${parts.length} of ${total} parts`
    : "Downloading star data from Gaia (via GAVO)…");
  report();
  const fetchNext = async () => {
    while (queue.length) {
      if (Date.now() > deadline) throw new GavoError("the field took too long to download", true);
      const tile = /** @type {{query: string, radius: number}} */ (queue.shift());
      try {
        parts.push(await fetchGavo(tile.query));
      } catch (error) {
        if (!(error instanceof GavoError && error.timedOut)) throw error;
        const smaller = /** @type {{query: string, radius: number}[]} */ (await call("split", { tile }));
        if (!smaller.length) throw new TooDenseError("dense field");
        queue.push(...smaller);
        total += smaller.length - 1;
      }
      report();
    }
  };
  await Promise.all(Array.from({ length: GAVO_PARALLEL }, fetchNext));
  return parts;
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
  return { bytes: new Uint8Array(await response.arrayBuffer()), date: entry.date, from: entry.source || "ESA" };
}

/**
 * @param {Uint8Array | Uint8Array[]} bytes one result, or one per piece of a large field
 * @param {string} filename @param {string} source @param {number} maxRows
 * @param {boolean} [uploaded] a student's file: warn if it is mostly outside the field
 */
async function showStars(bytes, filename, source, maxRows, uploaded = false) {
  const field = fieldArgs();
  const title = objectName || `RA ${field.ra}, Dec ${field.dec}`;
  setStatus("Placing the stars and drawing the picture…");
  const result = await call(
    "showStars",
    { bytes, filename, source, title, maxRows, uploaded, ...field },
  );
  pictureField = { title, source, ...field };
  displayStars(result, source);
  setStatus(result.display
    ? `Plotted ${result.count.toLocaleString()} stars, with the picture settings saved in the file.`
    : `Plotted ${result.count.toLocaleString()} stars. Adjust the sliders, or click HR Diagram.`);
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
    let parts;
    let reason = "";
    let dense = false;
    try {
      parts = await fetchGavoTiles(plan.gavo_tiles);
    } catch (error) {
      dense = error instanceof TooDenseError || (error instanceof GavoError && error.timedOut);
      reason = error instanceof Error && error.name === "TimeoutError" ? "it did not answer in time" : String(error instanceof Error ? error.message : error);
    }
    if (parts) {
      await showStars(parts, "gavo.csv", "Gaia DR3 via GAVO", MAX_ROWS);
      return;
    }
    const saved = objectName ? await savedCopy(objectName, field) : null;
    if (saved) {
      await showStars(saved.bytes, "saved.csv", `Gaia DR3 saved copy from ${saved.from} (${saved.date})`, 0);
      setStatus(`The live data source was unavailable (${reason}), so a saved copy was used.`);
      return;
    }
    uploadDetails.open = true;
    showError(dense
      ? "This field has too many stars for the live data source to send in time (it is wide or in a dense part " +
        "of the Milky Way). Use a smaller Field of View or a lower Maximum Magnitude, or a file from the ESA " +
        "Gaia archive (below)."
      : `The live Gaia data source is unavailable (${reason}), and no saved copy covers this field. ` +
        "Try again in a few minutes, or use a file from the ESA Gaia archive (below).");
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
    await showStars(new Uint8Array(await file.arrayBuffer()), file.name, `Gaia DR3 from uploaded ESA file (${file.name})`, 0, true);
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
  // Remember where the sliders started (as slidershow in ass.py starts them), for Reset.
  initialSliders = Object.fromEntries(Object.entries(sliders).map(([name, slider]) => [name, slider.value]));
  if (result.display) {
    // A star-data file downloaded from this page: put the sliders where they were (Reset still
    // returns to the standard starting positions). Values outside a slider's range are clamped.
    for (const [name, slider] of Object.entries(sliders)) {
      const value = result.display[name];
      if (typeof value === "string" && Number.isFinite(Number(value))) slider.value = value;
    }
    notes.append(Object.assign(document.createElement("li"), { textContent: "Picture settings restored from the file." }));
  }
  magtype.value = "Apparent"; // ass.py resets the magnitude type on Show Stars
  hrImage.hidden = true;
  hrDownloads.hidden = true;
  hrCredit.hidden = true;
  starsContent.hidden = false;
  starsEmpty.hidden = true;
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

/** @type {Record<string, string>} */
let initialSliders = {};

function resetSliders() {
  for (const [name, slider] of Object.entries(sliders)) {
    if (name in initialSliders) slider.value = initialSliders[name];
  }
  draw();
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
  hrDownloads.hidden = true;
  hrCredit.hidden = true;
  refresh();
  setStatus("Drawing the HR diagram…");
  try {
    const png = /** @type {Uint8Array<ArrayBuffer>} */ (await call("hrDiagram", { magnitudeType: magtype.value }));
    if (hrUrl) URL.revokeObjectURL(hrUrl);
    hrUrl = URL.createObjectURL(new Blob([png], { type: "image/png" }));
    hrImage.src = hrUrl;
    hrImage.hidden = false;
    hrDownload.href = hrUrl;
    hrDownload.download = `${objectName || "field"}_HR_${magtype.value}.png`;
    hrDownloads.hidden = false;
    hrCredit.hidden = false;
    hrMagnitudeType = magtype.value;
    setStatus("HR diagram ready. Use Download HR diagram to save it.");
  } catch (error) {
    report(error, "The HR diagram could not be drawn");
  } finally {
    busy = false;
    refresh();
  }
}

// ---- CSV downloads: the star data used, and the HR diagram's points --------------------

/** The diagram's magnitude type, so its CSV matches what is shown even if the menu changed. */
let hrMagnitudeType = "Apparent";

/** Today's date in the student's own time zone, for file names. */
function today() {
  const now = new Date();
  return `${now.getFullYear()}-${String(now.getMonth() + 1).padStart(2, "0")}-${String(now.getDate()).padStart(2, "0")}`;
}

/** @param {Uint8Array<ArrayBuffer>} data @param {string} name */
function saveFile(data, name) {
  const url = URL.createObjectURL(new Blob([data], { type: "text/csv;charset=utf-8" }));
  const link = Object.assign(document.createElement("a"), { href: url, download: name });
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 10000);
}

/** @param {string} action @param {Record<string, unknown>} args @param {string} what */
async function downloadCsv(action, args, what) {
  if (!ready || busy) return;
  clearError();
  try {
    const file = /** @type {{name: string, data: Uint8Array<ArrayBuffer>}} */ (await call(action, { ...args, date: today() }));
    saveFile(file.data, file.name);
    setStatus(`Saved ${what} as ${file.name}.`);
  } catch (error) {
    report(error, `The ${what} could not be saved`);
  }
}

// ---- Recording the picture settings ------------------------------------------------------

/** @type {{title: string, source: string, ra: string, dec: string, fov: number, maxmag: number} | null} */
let pictureField = null;

/**
 * The sliders' exact positions (restored when a star-data file is uploaded again), then the
 * values as shown beside them. Parsed by showclusters.display_settings.
 */
function displaySummary() {
  const exact = Object.entries(sliders).map(([name, slider]) => `${name}=${slider.value}`).join("; ");
  const shown = [["color saturation", "colsat"], ["nonlinear stretch", "gamma"], ["min", "mini"], ["max", "maxi"]]
    .map(([label, name]) => `${label} ${$(`${name}-value`).textContent}`).join(", ");
  return `${exact} (shown as ${shown})`;
}

/** Text fields for the picture PNG (standard PNG keywords; ASCII only). */
function pictureText() {
  const f = pictureField;
  const field = f ? `${f.title}, RA ${f.ra}, Dec ${f.dec}, FOV ${f.fov} deg, G < ${f.maxmag}` : "";
  return [
    ["Title", `Show Clusters picture: ${field}`],
    ["Description", `Display sliders: ${displaySummary()}`],
    ["Source", f ? f.source : ""],
    ["Copyright", "Gaia DR3 data: ESA/Gaia/DPAC, CC BY-NC 3.0 IGO. Bright stars: HYG database v3 (David Nash), CC BY-SA 2.5."],
    ["Software", "Show Clusters (StackImages), https://princetonuniversity.github.io/StackImages/ShowClusters/"],
  ];
}

/** @type {Uint32Array | null} */
let crcTable = null;
/** @param {Uint8Array} bytes */
function crc32(bytes) {
  if (!crcTable) {
    crcTable = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
      let c = n;
      for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
      crcTable[n] = c >>> 0;
    }
  }
  let crc = 0xffffffff;
  for (const byte of bytes) crc = crcTable[(crc ^ byte) & 0xff] ^ (crc >>> 8);
  return (crc ^ 0xffffffff) >>> 0;
}

/**
 * Insert tEXt chunks after a PNG's header chunk (IHDR ends at byte 33).
 * @param {Uint8Array} png @param {string[][]} entries [keyword, text] pairs
 * @returns {Uint8Array<ArrayBuffer>}
 */
function pngWithText(png, entries) {
  const ascii = (/** @type {string} */ text) =>
    Uint8Array.from(text.replace(/[^\x20-\x7e]/g, "?"), (c) => c.charCodeAt(0));
  const chunks = entries.filter(([, text]) => text).map(([keyword, text]) => {
    const body = new Uint8Array([...ascii("tEXt"), ...ascii(keyword), 0, ...ascii(text)]);
    const chunk = new Uint8Array(body.length + 8);
    const view = new DataView(chunk.buffer);
    view.setUint32(0, body.length - 4);
    chunk.set(body, 4);
    view.setUint32(body.length + 4, crc32(body));
    return chunk;
  });
  const size = png.length + chunks.reduce((sum, c) => sum + c.length, 0);
  const out = new Uint8Array(size);
  out.set(png.subarray(0, 33), 0);
  let offset = 33;
  for (const chunk of chunks) {
    out.set(chunk, offset);
    offset += chunk.length;
  }
  out.set(png.subarray(33), offset);
  return out;
}

// ---- Full screen ------------------------------------------------------------------------

function inFullscreen() {
  return viewer.classList.contains("fullscreen");
}

/** Picture and controls fill the window; also the browser's own full screen where allowed. */
function enterFullscreen() {
  if (!picture || inFullscreen()) return;
  viewer.classList.add("fullscreen");
  document.body.classList.add("no-scroll");
  exitButton.hidden = false;
  fullscreenButton.setAttribute("aria-pressed", "true");
  exitButton.focus();
  // Not every browser allows it (e.g. iPhone Safari); the full-window view works regardless.
  viewer.requestFullscreen?.().catch(() => {});
}

function exitFullscreen() {
  if (!inFullscreen()) return;
  viewer.classList.remove("fullscreen");
  document.body.classList.remove("no-scroll");
  exitButton.hidden = true;
  fullscreenButton.setAttribute("aria-pressed", "false");
  if (document.fullscreenElement) document.exitFullscreen().catch(() => {});
  fullscreenButton.focus();
}

// ---- FOV in arcminutes ------------------------------------------------------------------

/** Show the field of view in arcminutes beside its degrees, e.g. "degrees (18 arcminutes)". */
function showFovArcminutes() {
  const degrees = Number(fovInput.value);
  const arcminutes = Math.round(degrees * 60 * 100) / 100; // at most two decimals, no float noise
  $("fov-unit").textContent = fovInput.value.trim() !== "" && Number.isFinite(degrees) && degrees > 0
    ? `degrees (${arcminutes.toLocaleString("en-US", { maximumFractionDigits: 2 })} arcminutes)`
    : "degrees";
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
fovInput.addEventListener("input", showFovArcminutes);
form.addEventListener("submit", (event) => {
  event.preventDefault();
  void onShowStars();
});
for (const slider of Object.values(sliders)) slider.addEventListener("input", queueDraw);
resetButton.addEventListener("click", resetSliders);
fullscreenButton.addEventListener("click", enterFullscreen);
exitButton.addEventListener("click", exitFullscreen);
document.addEventListener("keydown", (event) => {
  if (event.key === "Escape" && inFullscreen()) exitFullscreen();
});
// Leaving the browser's full screen another way (its own Esc handling, a gesture) leaves ours too.
document.addEventListener("fullscreenchange", () => {
  if (!document.fullscreenElement && inFullscreen()) exitFullscreen();
});
hrButton.addEventListener("click", () => void onHrDiagram());
dataDownload.addEventListener("click", (event) => {
  event.preventDefault();
  void downloadCsv("starData", { display: displaySummary() }, "star data");
});
hrCsvDownload.addEventListener("click", (event) => {
  event.preventDefault();
  void downloadCsv("hrPoints", { magnitudeType: hrMagnitudeType }, "HR diagram points");
});
lookupButton.addEventListener("click", () => void lookup());
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
  canvas.toBlob(async (blob) => {
    if (!blob) return;
    const png = pngWithText(new Uint8Array(await blob.arrayBuffer()), pictureText());
    const url = URL.createObjectURL(new Blob([png], { type: "image/png" }));
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
showFovArcminutes();
refresh();
