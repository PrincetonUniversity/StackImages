// @ts-check
// Runs Pyodide and web/showclusters.py off the main thread. All files come from this
// site (vendor/ and data/); the page itself does every network request to Gaia services.
/** @type {any} */
const context = self;
/** @type {any} */
let runtime;
/** @type {any} */
let session;
let stage = "Python runtime";

/** @param {string} path */
async function fetchBytes(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`${path}: HTTP ${response.status}`);
  return new Uint8Array(await response.arrayBuffer());
}

try {
  const runtimeModule = new URL("./vendor/pyodide.mjs", import.meta.url).href;
  const { loadPyodide } = await import(runtimeModule);
  runtime = await loadPyodide({ indexURL: new URL("./vendor/", import.meta.url).href });
  stage = "Scientific packages";
  await runtime.loadPackage(["numpy", "matplotlib", "astropy"]);
  stage = "Star catalogs";
  runtime.FS.mkdirTree("/data/csv");
  for (const name of ["NI2021_MessierObjects.csv", "NI2021_NGCObjects.csv", "ICObj.csv"]) {
    runtime.FS.writeFile(`/data/csv/${name}`, await fetchBytes(`data/csv/${name}`));
  }
  runtime.FS.writeFile("/data/bright_stars.npy", await fetchBytes("data/bright_stars.npy"));
  runtime.FS.writeFile("/data/psf.npy", await fetchBytes("data/psf.npy"));
  stage = "Analysis code";
  runtime.FS.writeFile("/home/pyodide/showclusters.py", await fetchBytes("showclusters.py"));
  session = runtime.runPython(`
import numpy as np
from pathlib import Path
import showclusters as sc
sc.Session(Path("/data/csv"), np.load("/data/bright_stars.npy"), np.load("/data/psf.npy"))
`);
  context.postMessage({ type: "ready" });
} catch (error) {
  context.postMessage({ type: "fatal", stage, message: String(error) });
}

/** Convert a Python result to plain JS, keeping bytes as Uint8Array. @param {any} value */
function toJs(value) {
  if (value === null || value === undefined || typeof value !== "object") return value;
  const converted = value.toJs({ dict_converter: Object.fromEntries, create_pyproxies: false });
  value.destroy();
  return converted;
}

/** Python InputError messages are meant for students; anything else is a bug. @param {unknown} error */
function explain(error) {
  const text = String(error);
  const match = text.match(/InputError: ([^\n]+)/);
  return match ? { message: match[1], expected: true } : { message: text.slice(-1500), expected: false };
}

context.onmessage = (/** @type {MessageEvent} */ { data }) => {
  if (!session) return;
  const { id, action, args } = data;
  try {
    let result;
    if (action === "lookup") {
      result = toJs(session.lookup(args.name));
    } else if (action === "plan") {
      result = toJs(session.plan(args.ra, args.dec, args.fov, args.maxmag));
    } else if (action === "split") {
      result = toJs(session.split(args.tile));
    } else if (action === "showStars") {
      const py = session.show_stars(
        args.bytes, args.filename, args.source, args.ra, args.dec, args.fov, args.maxmag,
        args.title, args.maxRows, Boolean(args.uploaded),
      );
      result = toJs(py);
      context.postMessage({ id, ok: true, result }, [result.image.buffer]);
      return;
    } else if (action === "starData" || action === "hrPoints") {
      const py = action === "starData"
        ? session.star_data(args.date, args.display || "")
        : session.hr_points(args.magnitudeType, args.date);
      result = toJs(py);
      context.postMessage({ id, ok: true, result }, [result.data.buffer]);
      return;
    } else if (action === "hrDiagram") {
      result = toJs(session.hr_diagram(args.magnitudeType));
      context.postMessage({ id, ok: true, result }, [result.buffer]);
      return;
    } else {
      throw new Error(`Unknown action ${action}`);
    }
    context.postMessage({ id, ok: true, result });
  } catch (error) {
    context.postMessage({ id, ok: false, ...explain(error) });
  }
};
