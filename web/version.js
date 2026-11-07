// @ts-check
// The site's version, shown in each page's footer ([data-site-version]). When CI deploys,
// scripts/web_version.py replaces the line below with the real values; a build that was
// never stamped (a local or test build) says "Development build". Loaded with defer.
const SITE_VERSION = { version: "", date: "" };
for (const label of document.querySelectorAll("[data-site-version]")) {
  label.textContent = SITE_VERSION.version
    ? `Version ${SITE_VERSION.version} (${SITE_VERSION.date})`
    : "Development build";
}
