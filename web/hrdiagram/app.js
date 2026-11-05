// @ts-check
// Page logic for "HR Diagram of Your Own Image". Python (worker.js, hrimage.py) reads the
// image and does the science; this file handles the file, the picture sliders, the two
// diagrams and the downloads. The steps follow the StackImages desktop window.

import { sliderValues, stretchImage } from "../stretch.js"; // the site build puts these beside this file
import { setUpViewer } from "../viewer.js";

const MAX_FILE_BYTES = 250 * 1024 * 1024; // keep in step with hrimage.MAX_FILE_BYTES

/** @param {string} id */
function $(id) {
  const value = document.getElementById(id);
  if (!value) throw new Error(`Missing #${id}`);
  return value;
}

const fileInput = /** @type {HTMLInputElement} */ ($("file"));
const status = $("status");
const errorBox = $("error");
const progress = $("progress");
const progressText = $("progress-text");
const retry = $("retry");
const canvas = /** @type {HTMLCanvasElement} */ ($("picture"));
const brightestInput = /** @type {HTMLInputElement} */ ($("brightest"));
const faintestInput = /** @type {HTMLInputElement} */ ($("faintest"));
const hrButton = /** @type {HTMLButtonElement} */ ($("hr-button"));
const colorButton = /** @type {HTMLButtonElement} */ ($("color-button"));
const resetButton = /** @type {HTMLButtonElement} */ ($("reset-sliders"));
const fullscreenButton = /** @type {HTMLButtonElement} */ ($("fullscreen"));
const diagram = /** @type {HTMLImageElement} */ ($("diagram"));
const diagramSummary = $("diagram-summary");
const diagramDownloads = $("diagram-downloads");
const diagramDownload = /** @type {HTMLAnchorElement} */ ($("diagram-download"));
const pointsDownload = /** @type {HTMLAnchorElement} */ ($("points-download"));
const pictureDownload = /** @type {HTMLAnchorElement} */ ($("picture-download"));
const sliders = {
  colsat: /** @type {HTMLInputElement} */ ($("colsat")),
  gamma: /** @type {HTMLInputElement} */ ($("gamma")),
  mini: /** @type {HTMLInputElement} */ ($("mini")),
  maxi: /** @type {HTMLInputElement} */ ($("maxi")),
};

$("max-mb").textContent = String(MAX_FILE_BYTES / 1024 / 1024);

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
  image = null;
  showPictureEmpty();
  retry.hidden = true;
  $("load-details").hidden = true;
  setStatus("Loading Python… The first visit downloads about 60 MB.");
  worker = new Worker(new URL("worker.js", import.meta.url)); // classic: see worker.js
  worker.onmessage = ({ data }) => {
    if (data.type === "ready") {
      ready = true;
      setStatus("Ready. Choose an RGB FITS image.");
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
    for (const waiting of pending.values()) waiting.reject(new StudentError(
      "The Python worker stopped, probably out of memory. Click Reload tool, then try a smaller crop.",
    ));
    pending.clear();
    ready = false;
    retry.hidden = false;
    refresh();
  };
  refresh();
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
}

function refresh() {
  fileInput.disabled = !ready || busy;
  hrButton.disabled = colorButton.disabled = !ready || busy || !image;
  resetButton.disabled = busy || !image;
  fullscreenButton.disabled = !image;
}

/** @param {unknown} error @param {string} fallback */
function report(error, fallback) {
  if (error instanceof StudentError) showError(error.message);
  else showError(`${fallback}. Details: ${error instanceof Error ? error.message : String(error)}`);
  setStatus(ready ? "Ready." : "Stopped.");
}

/**
 * Run `work` with the page busy and a clock in the status line.
 * @template T @param {string} text @param {() => Promise<T>} work @returns {Promise<T>}
 */
async function working(text, work) {
  busy = true;
  clearError();
  refresh();
  const started = Date.now();
  const tick = () => {
    const seconds = Math.round((Date.now() - started) / 1000);
    progressText.textContent = seconds >= 2 ? `${seconds} s so far` : "";
  };
  setStatus(text);
  progress.hidden = false;
  tick();
  const clock = setInterval(tick, 1000);
  try {
    return await work();
  } finally {
    clearInterval(clock);
    progress.hidden = true;
    busy = false;
    refresh();
  }
}

// ---- 1. Open an RGB fits file -------------------------------------------------------

/** @type {{name: string, filename: string, base: Float32Array, width: number, height: number} | null} */
let image = null;
/** @type {Record<string, string>} */
let initialSliders = {};

function showPictureEmpty() {
  exitFullscreen(); // a new or failed file never leaves the old picture full screen
  $("picture-content").hidden = true;
  $("picture-empty").hidden = false;
  clearDiagram();
}

async function onFile() {
  const file = fileInput.files?.[0];
  if (!file || !ready || busy) return;
  image = null;
  showPictureEmpty();
  if (file.size > MAX_FILE_BYTES) {
    showError(`This file is larger than ${MAX_FILE_BYTES / 1024 / 1024} MB. Use a smaller crop.`);
    refresh();
    return;
  }
  try {
    const result = await working(`Reading ${file.name}…`, async () => {
      const bytes = new Uint8Array(await file.arrayBuffer());
      return call("open", { bytes, filename: file.name }, [bytes.buffer]);
    });
    showPicture(result, file.name);
    setStatus("Picture ready. Adjust it with the sliders, then click HR diagram.");
  } catch (error) {
    report(error, "The image could not be opened");
  } finally {
    refresh();
  }
}

/** @param {any} result @param {string} filename */
function showPicture(result, filename) {
  const bytes = /** @type {Uint8Array} */ (result.image);
  image = {
    name: result.name,
    filename,
    base: new Float32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4),
    width: result.width,
    height: result.height,
  };
  $("object-name").textContent = result.name || "—";
  const [n, m] = result.pixels;
  $("image-size").textContent = `${n.toLocaleString()} × ${m.toLocaleString()} pixels`;
  const hasPosition = result.ra_hours !== null && result.dec_deg !== null;
  $("position-label").hidden = $("position").hidden = !hasPosition;
  if (hasPosition) {
    $("position").textContent = `${formatAngle(result.ra_hours)} h, ${formatAngle(result.dec_deg)}°`;
  }
  $("picture-note").textContent = result.factor > 1
    ? `Shown at 1/${result.factor} size (${result.width} × ${result.height}); the diagrams use every pixel.`
    : "";
  canvas.width = result.width;
  canvas.height = result.height;
  // slidershow's slider ranges and starting values (sqrt scale for min/max).
  const low = 0.8 * Math.sqrt(Math.max(result.pct_lo, 0));
  for (const slider of [sliders.mini, sliders.maxi]) {
    slider.min = String(low);
    slider.max = "1";
  }
  sliders.mini.value = String(Math.sqrt(Math.max(result.pct_lo, 0)));
  sliders.maxi.value = String(Math.sqrt(result.pct_hi));
  sliders.colsat.value = "5";
  sliders.gamma.value = "2";
  initialSliders = Object.fromEntries(Object.entries(sliders).map(([name, slider]) => [name, slider.value]));
  $("picture-content").hidden = false;
  $("picture-empty").hidden = true;
  pictureDownload.download = `${slug(result.name || filename)}.png`;
  draw();
}

/** Sexagesimal, as the desktop app's boxes show it. @param {number} value */
function formatAngle(value) {
  const sign = value < 0 ? "-" : "+";
  let rest = Math.abs(value);
  const whole = Math.floor(rest);
  rest = (rest - whole) * 60;
  const minutes = Math.floor(rest);
  const seconds = Math.round((rest - minutes) * 60);
  return `${sign}${String(whole).padStart(2, "0")} : ${String(minutes).padStart(2, "0")} : ${String(seconds).padStart(2, "0")}`;
}

/** @param {string} text */
function slug(text) {
  return text.replace(/\.fi?ts?$/i, "").replace(/[^A-Za-z0-9]+/g, "_").replace(/^_+|_+$/g, "") || "image";
}

// ---- 2. The picture ------------------------------------------------------------------

function draw() {
  if (!image) return;
  const values = sliderValues({
    colsat: sliders.colsat.value, gamma: sliders.gamma.value, mini: sliders.mini.value, maxi: sliders.maxi.value,
  });
  for (const [key, value] of Object.entries(values)) {
    $(`${key}-value`).textContent = value.toFixed(key === "gamma" || key === "colsat" ? 2 : 4);
  }
  const picture = stretchImage(image.base, image.width, image.height, values);
  /** @type {CanvasRenderingContext2D} */ (canvas.getContext("2d")).putImageData(picture, 0, 0);
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

function resetSliders() {
  for (const [name, slider] of Object.entries(sliders)) {
    if (name in initialSliders) slider.value = initialSliders[name];
  }
  draw();
}

/** @param {Event} event */
function downloadPicture(event) {
  event.preventDefault();
  canvas.toBlob((blob) => {
    if (blob) saveFile(blob, pictureDownload.download || "picture.png");
  }, "image/png");
}

// ---- 3. HR diagram and Color Diagram ------------------------------------------------------

/** @type {string | null} */
let diagramUrl = null;
/** The magnitudes the diagram on screen was drawn with, so its CSV matches it. */
let shownMags = { brightest: "", faintest: "20" };
let measured = false;

function clearDiagram() {
  diagram.hidden = true;
  diagramSummary.hidden = true;
  diagramDownloads.hidden = true;
  measured = false;
}

/** @param {Blob} blob @returns {Promise<string>} */
function dataUrl(blob) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(String(reader.result));
    reader.onerror = () => reject(reader.error);
    reader.readAsDataURL(blob);
  });
}

/** @param {"hrDiagram" | "colorDiagram"} action */
async function onDiagram(action) {
  if (!ready || busy || !image) return;
  const mags = { brightest: brightestInput.value, faintest: faintestInput.value };
  diagram.hidden = true; // never show a diagram that doesn't match the current settings
  diagramSummary.hidden = true;
  diagramDownloads.hidden = true;
  const what = action === "hrDiagram" ? "HR diagram" : "color diagram";
  const text = measured
    ? `Drawing the ${what}…`
    : "Finding and measuring stars… A large image can take a minute or two; please keep this tab open.";
  try {
    const result = await working(text, () => call(action, mags));
    measured = true;
    const blob = new Blob([result.png], { type: "image/png" });
    if (diagramUrl) URL.revokeObjectURL(diagramUrl);
    diagramUrl = URL.createObjectURL(blob);
    // A data: URL for display: strict Content-Security-Policies (vanderbei.princeton.edu)
    // allow data: images but not blob: ones. See ShowClusters' app.js.
    diagram.src = await dataUrl(blob);
    diagram.alt = action === "hrDiagram"
      ? `Hertzsprung–Russell diagram of ${image.name || image.filename}`
      : `Color diagram (stellar locus) of ${image.name || image.filename}`;
    diagram.hidden = false;
    diagramSummary.textContent = `${result.plotted.toLocaleString()} stars plotted, of ${result.detected.toLocaleString()} measured.`;
    diagramSummary.hidden = false;
    diagramDownload.href = diagramUrl;
    diagramDownload.download = `${slug(image.name || image.filename)}_${action === "hrDiagram" ? "HR" : "Color"}.png`;
    diagramDownloads.hidden = false;
    shownMags = mags;
    setStatus(`${action === "hrDiagram" ? "HR diagram" : "Color diagram"} ready.`);
  } catch (error) {
    report(error, `The ${what} could not be drawn`);
  }
}

/** @param {Blob} blob @param {string} name */
function saveFile(blob, name) {
  const url = URL.createObjectURL(blob);
  const link = Object.assign(document.createElement("a"), { href: url, download: name });
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 60000);
}

/** @param {Event} event */
async function downloadPoints(event) {
  event.preventDefault();
  if (!ready || busy) return;
  try {
    const result = await working("Preparing the points…", () => call("points", shownMags));
    saveFile(new Blob([result.data], { type: "text/csv;charset=utf-8" }), result.name);
    setStatus("Points downloaded.");
  } catch (error) {
    report(error, "The points could not be saved");
  }
}

// ---- Wiring ------------------------------------------------------------------------------

fileInput.addEventListener("change", () => void onFile());
for (const slider of Object.values(sliders)) slider.addEventListener("input", queueDraw);
const { exitFullscreen } = setUpViewer({ sliders, initial: () => initialSliders, reset: resetSliders, hasPicture: () => Boolean(image) });
pictureDownload.addEventListener("click", downloadPicture);
hrButton.addEventListener("click", () => void onDiagram("hrDiagram"));
colorButton.addEventListener("click", () => void onDiagram("colorDiagram"));
pointsDownload.addEventListener("click", (event) => void downloadPoints(event));
$("mags").addEventListener("submit", (event) => {
  event.preventDefault();
  void onDiagram("hrDiagram");
});
for (const input of [brightestInput, faintestInput]) {
  input.addEventListener("keydown", (event) => {
    if (event.key === "Enter") {
      event.preventDefault();
      void onDiagram("hrDiagram");
    }
  });
}
retry.addEventListener("click", () => {
  worker.terminate();
  clearError();
  startWorker();
});
startWorker();
