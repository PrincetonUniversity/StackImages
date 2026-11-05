// @ts-check
// Runs Pyodide and web/hrimage.py off the main thread. Everything comes from this site;
// Pyodide itself is shared with Show Clusters (../ShowClusters/vendor/). The student's
// image is never sent anywhere.
// Wrapped in a function: tsc checks both classic workers in one global scope.
(() => {
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

  // A classic worker, and .wasm compiled from bytes when its type is wrong: the same reasons
  // as ShowClusters' worker.js (static hosts that don't know these file types).
  const streaming = WebAssembly.instantiateStreaming;
  if (typeof streaming === "function") {
    WebAssembly.instantiateStreaming = async (source, imports) => {
      const response = await source;
      if (/^application\/wasm\b/i.test(response.headers.get("content-type") ?? "")) {
        return streaming(response, imports);
      }
      return WebAssembly.instantiate(await response.arrayBuffer(), imports);
    };
  }

  async function start() {
    const vendor = new URL("../ShowClusters/vendor/", self.location.href).href;
    context.importScripts(`${vendor}pyodide.js`);
    runtime = await context.loadPyodide({ indexURL: vendor });
    stage = "Scientific packages";
    await runtime.loadPackage(["numpy", "scipy", "matplotlib", "astropy"]);
    stage = "Analysis code";
    for (const name of ["showclusters.py", "hrimage.py"]) {
      runtime.FS.writeFile(`/home/pyodide/${name}`, await fetchBytes(name));
    }
    session = runtime.runPython("import hrimage\nhrimage.Session()");
  }

  start().then(
    () => context.postMessage({ type: "ready" }),
    (error) => context.postMessage({ type: "fatal", stage, message: String(error) }),
  );

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
    if (match) return { message: match[1], expected: true };
    if (/MemoryError|memory access out of bounds|Cannot allocate|out of memory/i.test(text)) {
      return { message: "This image needs more memory than the browser allows. Use a smaller crop.", expected: true };
    }
    return { message: text.slice(-1500), expected: false };
  }

  context.onmessage = (/** @type {MessageEvent} */ { data }) => {
    if (!session) return;
    const { id, action, args } = data;
    try {
      if (action === "open") {
        const result = toJs(session.open_image(args.bytes, args.filename));
        context.postMessage({ id, ok: true, result }, [result.image.buffer]);
      } else if (action === "hrDiagram" || action === "colorDiagram") {
        const py = action === "hrDiagram"
          ? session.hr_diagram(args.brightest, args.faintest)
          : session.color_diagram(args.brightest, args.faintest);
        const result = toJs(py);
        context.postMessage({ id, ok: true, result }, [result.png.buffer]);
      } else if (action === "points") {
        const result = toJs(session.points(args.brightest, args.faintest));
        context.postMessage({ id, ok: true, result }, [result.data.buffer]);
      } else {
        throw new Error(`Unknown action ${action}`);
      }
    } catch (error) {
      context.postMessage({ id, ok: false, ...explain(error) });
    }
  };
})();
