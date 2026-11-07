// @ts-check
// The picture viewer's full screen and its "Reset sliders?" question, shared by the Show
// Clusters and HR diagram pages. Both pages use the same markup: a .viewer region holding
// #fullscreen (the icon on the picture), #exit-fullscreen, #reset-sliders and the
// #reset-confirm dialog with #reset-ok and #reset-cancel.

/** @param {string} id */
function $(id) {
  const value = document.getElementById(id);
  if (!value) throw new Error(`Missing #${id}`);
  return value;
}

/**
 * Wire up full screen and the reset question.
 * @param {{
 *   sliders: Record<string, HTMLInputElement>,
 *   initial: () => Record<string, string>,
 *   reset: () => void,
 *   hasPicture: () => boolean,
 * }} options
 */
export function setUpViewer({ sliders, initial, reset, hasPicture }) {
  const viewer = /** @type {HTMLElement} */ (document.querySelector(".viewer"));
  const fullscreenButton = /** @type {HTMLButtonElement} */ ($("fullscreen"));
  const exitButton = /** @type {HTMLButtonElement} */ ($("exit-fullscreen"));
  const resetButton = /** @type {HTMLButtonElement} */ ($("reset-sliders"));
  const resetConfirm = /** @type {HTMLDialogElement} */ ($("reset-confirm"));

  function inFullscreen() {
    return viewer.classList.contains("fullscreen");
  }

  /** Picture and controls fill the window; also the browser's own full screen where allowed. */
  function enterFullscreen() {
    if (!hasPicture() || inFullscreen()) return;
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

  /** Reset loses the student's tuning, so ask first, in the page's own dialog. Nothing moved: just reset. */
  function askToResetSliders() {
    const start = initial();
    const moved = Object.entries(sliders).some(([name, slider]) => name in start && slider.value !== start[name]);
    if (!moved) {
      reset();
      return;
    }
    resetConfirm.showModal(); // focuses "Keep my settings" (autofocus); Esc keeps them too
  }

  resetButton.addEventListener("click", askToResetSliders);
  $("reset-ok").addEventListener("click", () => {
    reset();
    resetConfirm.close();
  });
  $("reset-cancel").addEventListener("click", () => resetConfirm.close());
  // The close event arrives a moment after the dialog closes: return focus to Reset only if it
  // hasn't gone elsewhere since (e.g. Esc again, leaving full screen, focuses that button).
  resetConfirm.addEventListener("close", () => {
    const at = document.activeElement;
    if (!at || at === document.body || resetConfirm.contains(at)) resetButton.focus();
  });
  fullscreenButton.addEventListener("click", enterFullscreen);
  exitButton.addEventListener("click", exitFullscreen);
  document.addEventListener("keydown", (event) => {
    // Esc in the reset question closes the question only.
    if (event.key === "Escape" && inFullscreen() && !resetConfirm.open) exitFullscreen();
  });
  // Leaving the browser's full screen another way (its own Esc handling, a gesture) leaves ours too.
  document.addEventListener("fullscreenchange", () => {
    if (!document.fullscreenElement && inFullscreen()) exitFullscreen();
  });

  return { exitFullscreen };
}
