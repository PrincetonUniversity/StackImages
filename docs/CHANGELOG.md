# Changelog: StackImages web pages

Versions of <https://princetonuniversity.github.io/StackImages/>, shown at the bottom of every
page and released as `web-1.N`. Each deploy that changes the site gets the next number; add an
entry under the number CI names in the pull request, in words a student would follow (it becomes
the release notes). Raise the major number in `web/VERSION` for a milestone. See
[WEB.md](WEB.md#versions).

## 1.0 (2026-11-07)

The first numbered version. What the pages do now:

- **Show Clusters** draws a star field and its HR diagram from Gaia data, the same way as
  the desktop app's ShowClusters. Enter a Messier, NGC or IC object, or type RA and Dec.
- **Fields up to 90°**, for whole constellations (for example Orion at 60°). Wide fields up
  to magnitude 8.5 come back in seconds. Wider and deeper ones (for example 45° at
  magnitude 10, or 30° at magnitude 12) take 2 to 5 minutes: keep the tab open while the
  clock runs. Above 5°, the HR diagram starts on absolute magnitudes, as in the desktop app.
- **Your own Gaia file:** a file downloaded from the ESA Gaia archive can be used instead. A
  file that ESA cut off at 50,000 stars (its limit for users who aren't logged in) is
  refused, because the field would be incomplete.
- **HR Diagram of Your Own Image** makes an HR diagram and a color diagram from your own
  FITS image, entirely in your browser.
- Picture sliders, full screen, and downloads of the picture, the star data and the HR
  diagram (PNG and CSV).
