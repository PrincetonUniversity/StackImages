// @ts-check
// slidershow's display stretch (ass.py), shared by the Show Clusters and HR diagram pages.
// `base` is display_base's output (web/showclusters.py): width × height × 3 floats.

/**
 * Slider positions as slidershow uses them: the page's sliders hold the raw values.
 * @param {{colsat: string, gamma: string, mini: string, maxi: string}} raw
 */
export function sliderValues(raw) {
  return {
    colsat: Number(raw.colsat),
    gamma: 3 ** Number(raw.gamma) - 2,
    mini: Number(raw.mini) ** 2,
    maxi: Number(raw.maxi) ** 2,
  };
}

/**
 * The stretched picture for the given slider values.
 * @param {Float32Array} base @param {number} width @param {number} height
 * @param {{colsat: number, gamma: number, mini: number, maxi: number}} values
 * @returns {ImageData}
 */
export function stretchImage(base, width, height, { colsat, gamma, mini, maxi }) {
  const pixels = width * height;
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
  const image = new ImageData(width, height);
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
  return image;
}
