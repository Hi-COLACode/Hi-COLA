// Resolve the version directory from this script's URL so archived pages need no edits.
const versionRoot = new URL("../", document.currentScript.src);
const siteRoot = new URL("../", versionRoot);
document.addEventListener("DOMContentLoaded", async () => {
  const selector = document.getElementById("release-switcher");
  if (!selector) return;
  try {
    const response = await fetch(new URL("versions.json", siteRoot));
    if (!response.ok) return;
    const versions = await response.json();
    const page = window.location.pathname.slice(versionRoot.pathname.length);
    for (const version of versions) {
      const root = new URL(`${encodeURIComponent(version.version)}/`, siteRoot);
      const option = new Option(version.name, new URL(page || "index.html", root).href);
      option.dataset.home = root.href;
      option.selected = root.pathname === versionRoot.pathname;
      selector.add(option);
    }
    selector.hidden = false;
    selector.addEventListener("change", async () => {
      const home = selector.selectedOptions[0].dataset.home;
      try {
        const page = await fetch(selector.value, { method: "HEAD" });
        window.location.assign(page.ok ? selector.value : home);
      } catch (error) {
        window.location.assign(home);
      }
    });
  } catch (error) {
    console.debug("Version list unavailable", error);
  }
});
