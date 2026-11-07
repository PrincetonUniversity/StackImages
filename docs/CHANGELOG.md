# Changelog: StackImages web pages

The web pages at <https://princetonuniversity.github.io/StackImages/>: Show Clusters and HR
Diagram of Your Own Image. Each version is shown at the bottom of every page ("Version 1.0
(2026-11-07)") and published as a `web-1.N` release with the site as a zip.

How versions work: each deploy that changes the site gets the next number (1.0, 1.1, …).
Add an entry here, under the number CI names in the pull request, in words a student would
follow; the release notes are this entry. Raise the first number in `web/VERSION` for a
milestone (2.0 starts again at .0). The desktop app has its own versions (v21.0, …).

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
