import json
import subprocess
import sys
import textwrap


def _run_and_capture_title(script_path, expected_title):
    # We run Python and pre-patch wx.Frame to capture titles created by the script.
    # The inline Python imports wx, wraps Frame, runs the target script via runpy,
    # and prints a JSON with the captured titles.
    runner = textwrap.dedent(
        """
    import runpy, sys, json
    import wx as _wx

    Captured = []

    _OrigFrame = getattr(_wx, 'Frame', None)

    class FrameWrapper(object):
        def __init__(self, *args, **kwargs):
            # title is sometimes passed positionally or as a kwarg
            if len(args) >= 2:
                title = args[1]
            else:
                title = kwargs.get('title', '')
            Captured.append(title)
            # create a thin surrogate object with Show(), SetSize(), Centre(), SetSizerAndFit()
            class Dummy:
                def Show(self): pass
                def SetSize(self, *a, **k): pass
                def Centre(self): pass
                def SetSizerAndFit(self, *a, **k): pass
            self._d = Dummy()

        def __getattr__(self, name):
            return getattr(self._d, name)

    _wx.Frame = FrameWrapper

    # run the script
    try:
        runpy.run_path(r'%s', run_name='__main__')
    except SystemExit:
        pass
    except Exception as e:
        # If the script errors it's okay; we just want the title capture
        pass

    print(json.dumps(Captured))
"""
        % script_path
    )

    cmd = [sys.executable, "-c", runner]
    # A GUI that reaches MainLoop() would block forever; fail loudly instead.
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=120)
    out = proc.stdout.strip()
    # Try to parse JSON from the last line of stdout
    if not out:
        return []
    try:
        captured = json.loads(out.splitlines()[-1])
    except Exception:
        # Fallback: try to parse entire stdout
        try:
            captured = json.loads(out)
        except Exception:
            captured = []
    return captured


def test_getfit_title():
    script = "GetFit.py"
    titles = _run_and_capture_title(script, "Download Astro Fits Files")
    assert "Download Astro Fits Files" in titles, f"Expected title not found in {titles}"


def test_stackimages_title():
    script = "StackImages.py"
    titles = _run_and_capture_title(script, "Align, Stack, Show")
    assert "Align, Stack, Show" in titles, f"Expected title not found in {titles}"
