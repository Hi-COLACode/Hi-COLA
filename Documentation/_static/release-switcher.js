// Fill the version selector and the version banner from the site's shared
// versions.json. Each version directory loads this script, so archived pages
// pick up new releases without being rebuilt.
const versionRoot = new URL("../", document.currentScript.src);
const siteRoot = new URL("../", versionRoot);

// Return the URL of the same page in another version, or that version's home
// page when the page does not exist there.
async function pageIn(version, page) {
  const root = new URL(`${encodeURIComponent(version)}/`, siteRoot);
  const target = new URL(page || "index.html", root).href;
  try {
    const response = await fetch(target, { method: "HEAD" });
    return response.ok ? target : root.href;
  } catch (error) {
    return root.href;
  }
}

// Show a banner on the development version and on releases other than the latest.
async function showBanner(versions, current, page) {
  const latest = versions.find((entry) => entry.latest) ||
    versions.find((entry) => entry.version !== "current");
  if (latest && current === latest.version) return;
  const message = current === "current" ? "Development version." : "This is an old version.";
  const banner = document.createElement("div");
  banner.className = "hicola-version-banner";
  banner.textContent = `${message} `;
  if (latest) {
    const link = document.createElement("a");
    link.href = await pageIn(latest.version, page);
    link.textContent = "Click here";
    banner.append(link, ` for the latest release, ${latest.name}.`);
  }
  document.querySelector(".bd-article")?.prepend(banner);
}

document.addEventListener("DOMContentLoaded", async () => {
  const selector = document.getElementById("release-switcher");
  const page = window.location.pathname.slice(versionRoot.pathname.length);
  const current = decodeURIComponent(versionRoot.pathname.split("/").filter(Boolean).pop() || "");
  try {
    const response = await fetch(new URL("versions.json", siteRoot));
    if (!response.ok) return;
    const versions = await response.json();
    if (selector) {
      for (const entry of versions) {
        const label = entry.latest ? `${entry.name} (latest)` : entry.name;
        const option = new Option(label, entry.version);
        option.selected = entry.version === current;
        selector.add(option);
      }
      selector.hidden = false;
      selector.addEventListener("change", async () => {
        window.location.assign(await pageIn(selector.value, page));
      });
    }
    await showBanner(versions, current, page);
  } catch (error) {
    console.debug("Version list unavailable", error);
  }
});
