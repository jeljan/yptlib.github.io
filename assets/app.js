const DATA_ROOTS = window.YPTLIB_DATA_ROOT
  ? [window.YPTLIB_DATA_ROOT]
  : [
      "assets/data/",
      "https://cdn.jsdelivr.net/gh/jeljan/yptlib@ba1864e98dcc4ceba93893650276211d0b182bd7/assets/data/",
      "https://cdn.statically.io/gh/jeljan/yptlib@ba1864e98dcc4ceba93893650276211d0b182bd7/assets/data/",
      "https://raw.githubusercontent.com/jeljan/yptlib/ba1864e98dcc4ceba93893650276211d0b182bd7/assets/data/",
    ];
const RELEASE_DATA_ARCHIVE_URL = window.YPTLIB_DATA_ARCHIVE_URL || "";
const RELEASE_DATA_ARCHIVE_ENABLED = Boolean(RELEASE_DATA_ARCHIVE_URL);
let releaseArchivePromise = null;
const ASSET_VERSION = "pages-data-v29";
const FETCH_RETRIES_PER_ROOT = 2;
const FETCH_TIMEOUT_MS = 12000;
const COMPOUND_CACHE_WARMUP_LIMIT = 24;
const SUMMARY_RENDER_DEFER_MS = 5000;
const SEARCHABLE_OPTION_LIMIT = 80;
const CONTACT_TYPES = ["PPI", "Dna", "Rna", "Metal", "Ligand", "Cofactor"];
const CIRCULAR_WARHEAD_CLASS_ORDER = ["FS", "SuFEx", "SuTEx"];
const MIN_CLUSTER_SITE_ROWS = 3;
const MIN_CLUSTER_SHARED_COMPOUNDS = 10;
const MIN_CLUSTER_VALID_PAIR_FRACTION = 0.5;
const MIN_COMPLETE_CLUSTER_COMPOUNDS = 3;
const MIN_VISIBLE_DENDROGRAM_SITES = 8;
const MIN_DENDROGRAM_STROKE_PX = 0.55;
const CIRCULAR_RING_WIDTH = 0.0775;
const CIRCULAR_RING_GAP = 0.0175;
const CIRCULAR_BANDS = {
  superfamily: [0.47 - CIRCULAR_RING_WIDTH, 0.47],
  family: [0.47 - CIRCULAR_RING_WIDTH * 2 - CIRCULAR_RING_GAP, 0.47 - CIRCULAR_RING_WIDTH - CIRCULAR_RING_GAP],
  domain: [0.47 - CIRCULAR_RING_WIDTH * 3 - CIRCULAR_RING_GAP * 2, 0.47 - CIRCULAR_RING_WIDTH * 2 - CIRCULAR_RING_GAP * 2],
  site: [0.47 - CIRCULAR_RING_WIDTH * 4 - CIRCULAR_RING_GAP * 3, 0.47 - CIRCULAR_RING_WIDTH * 3 - CIRCULAR_RING_GAP * 3],
};
const CIRCULAR_LAYER_LABELS = [
  { level: "superfamily", label: "Superfamily" },
  { level: "family", label: "Family" },
  { level: "domain", label: "Domain" },
  { level: "site", label: "Site" },
];
const CIRCULAR_LAYER_LABEL_ANGLE = Math.PI;
const CONTACT_LABELS = {
  PPI: "PPI",
  Dna: "DNA",
  Rna: "RNA",
  Metal: "Metal",
  Ligand: "Ligand",
  Cofactor: "Cofactor",
};
const TARGET_LIST_OPTIONS = [
  { value: "all", label: "Whole proteome" },
  { value: "cancer", label: "Cancer-driver genes" },
  { value: "PPI", label: "Sites at PPI interface" },
  { value: "Metal", label: "Sites near metal" },
  { value: "Cofactor", label: "Sites near cofactor" },
  { value: "Dna", label: "Sites near DNA" },
  { value: "Rna", label: "Sites near RNA" },
  { value: "Ligand", label: "Sites near ligand" },
  { value: "custom", label: "Custom gene list" },
];
const SIG_NEGLOG10 = -Math.log10(0.05);
const MORANDI = {
  blueDark: "#2f5f8f",
  blue: "#5f7f9a",
  blueMid: "#8ea9c1",
  blueLight: "#d8e3ed",
  gray: "#d7d3ca",
  green: "#6f8372",
  red: "#a65f58",
  amber: "#b49a6a",
  mauve: "#9c8795",
};
const PROTEOME_MUTED_GREY = "rgba(200,200,200,0.3)";
const PVALUE_BLUE_SCALE = [
  [0, "#1f4f82"],
  [0.02, "#2f6ea3"],
  [0.1, "#5f91bd"],
  [0.5, "#a9bfd1"],
  [1, "#dbe4eb"],
];
const SPECTRAL_HIT_BLUE = "#2b83ba";
const SPECTRAL_CANCER_RED = "#d7191c";
const CUSTOM_GENE_PURPLE = "#7b3294";
const ANGSTROM = "\u00c5";

const state = {
  manifest: null,
  datasetKey: "os",
  dataset: null,
  catalog: [],
  rowById: new Map(),
  compoundCache: new Map(),
  compoundLoadPromises: new Map(),
  rawHoverCache: new Map(),
  rawHoverAliasCache: new Map(),
  geneSiteCache: new Map(),
  contactFileIndex: null,
  contactFileIndexPromise: null,
  contactTargetRowCache: new Map(),
  contactTargetRowSetCache: new Map(),
  activeGeneFilterIndexCache: new Map(),
  contactCache: new Map(),
  filteredSitesByGene: new Map(),
  compoundChoiceCountCache: new Map(),
  compoundChoiceRefreshToken: 0,
  compoundChoiceRefreshTimer: null,
  compoundWarmupToken: 0,
  searchableControls: new Map(),
  persistentHitTooltip: null,
  copyToastTimer: null,
  summaryDatasetRenderToken: 0,
  summaryDatasetRenderTimer: null,
  summaryRenderTimer: null,
  currentSiteRow: null,
  currentSiteSummary: null,
  currentPdbEntries: [],
  activeSummaryDrug: null,
  activeSiteDrug: null,
  siteRefreshSeq: 0,
  globalProteome: null,
  globalProteomePromise: null,
  circularAtlas: null,
  circularSectors: [],
  circularCellMap: new Map(),
  circularVisibleSiteDendrogram: null,
  circularVisibleSiteDendrogramKey: "",
  circularVisibleSiteDendrogramPromise: null,
  circularConstrainedOrderMap: new Map(),
  circularConstrainedOrderKey: "",
  circularConstrainedOrderPromise: null,
  circularConstrainedOrderStatus: "",
  circularZoom: { scale: 1, offsetX: 0, offsetY: 0, dragging: false, lastX: 0, lastY: 0 },
};
const BIO_FAMILY_META = {
  catalytic: { label: "Catalytic/active-site", color: MORANDI.red, className: "bio-family-catalytic" },
  regulatory: { label: "PTM/regulatory", color: MORANDI.mauve, className: "bio-family-regulatory" },
  variant: { label: "Variant evidence", color: MORANDI.blueDark, className: "bio-family-variant" },
  structure: { label: "Structure/contact", color: MORANDI.green, className: "bio-family-structure" },
  constraint: { label: "Constraint/conservation", color: MORANDI.amber, className: "bio-family-constraint" },
};
const BIO_FAMILIES = Object.keys(BIO_FAMILY_META);

const el = (id) => document.getElementById(id);

function setStatus(text) {
  el("statusText").textContent = text;
}

function activeDatasetDrugs() {
  const noReplicateDrugs = new Set(state.dataset?.noReplicateDrugs || []);
  const tier2OnlyDrugs = new Set(state.dataset?.tier2OnlyDrugs || []);
  return (state.dataset?.rawDrugs || []).filter(
    (drug) => !checkboxChecked("dropNoReplicates") || (!noReplicateDrugs.has(drug) && !tier2OnlyDrugs.has(drug))
  );
}

function updateDatasetStatus() {
  if (!state.dataset || !state.catalog) return;
  const active = activeDatasetDrugs().length;
  setStatus(`${state.dataset.label}: ${state.catalog.length.toLocaleString()} sites, ${active.toLocaleString()} compounds`);
}

function normalizeArchivePath(path) {
  return path.split("?")[0].replace(/^\/+/, "");
}

async function getReleaseArchive() {
  if (!RELEASE_DATA_ARCHIVE_ENABLED) return null;
  if (typeof zip === "undefined") {
    throw new Error("zip.js did not load; cannot read release data archive");
  }
  if (!releaseArchivePromise) {
    releaseArchivePromise = (async () => {
      const response = await fetch(RELEASE_DATA_ARCHIVE_URL, { cache: "force-cache" });
      if (!response.ok) {
        throw new Error(`Unable to load release data archive: ${response.status}`);
      }
      const blob = await response.blob();
      const reader = new zip.ZipReader(new zip.BlobReader(blob));
      const entries = await reader.getEntries();
      const byName = new Map(entries.filter((entry) => !entry.directory).map((entry) => [entry.filename.replace(/^assets\/data\//, ""), entry]));
      return { reader, byName };
    })();
  }
  return releaseArchivePromise;
}

async function fetchJson(path) {
  if (RELEASE_DATA_ARCHIVE_ENABLED) {
    const archive = await getReleaseArchive();
    const normalized = normalizeArchivePath(path);
    const entry = archive.byName.get(normalized);
    if (!entry) {
      throw new Error(`Unable to load ${path}: missing from release archive`);
    }
    const text = await entry.getData(new zip.TextWriter());
    return JSON.parse(text);
  }
  const separator = path.includes("?") ? "&" : "?";
  let lastError = null;
  for (const root of DATA_ROOTS) {
    for (let attempt = 0; attempt < FETCH_RETRIES_PER_ROOT; attempt += 1) {
      let timeout = null;
      try {
        const controller = new AbortController();
        timeout = setTimeout(() => controller.abort(), FETCH_TIMEOUT_MS);
        const response = await fetch(`${root}${path}${separator}v=${ASSET_VERSION}`, { cache: "no-store", signal: controller.signal });
        if (!response.ok) {
          throw new Error(`${response.status}`);
        }
        return response.json();
      } catch (err) {
        lastError = err;
        if (attempt + 1 < FETCH_RETRIES_PER_ROOT) {
          await new Promise((resolve) => setTimeout(resolve, 250 * (attempt + 1)));
        }
      } finally {
        if (timeout) clearTimeout(timeout);
      }
    }
  }
  throw new Error(`Unable to load ${path}: ${lastError?.message || "network error"}`);
}

function setOptions(select, options, selected = null) {
  select.innerHTML = "";
  for (const option of options) {
    const node = document.createElement("option");
    node.value = option.value;
    node.textContent = option.label;
    node.label = option.label;
    select.appendChild(node);
  }
  if (selected != null) {
    select.value = selected;
  }
  syncSearchableControl(select.id);
}

function searchableLabel(option) {
  return option?.label || option?.textContent || option?.value || "";
}

function searchableOptions(select) {
  return [...(select?.options || [])].filter((option) => option.value || searchableLabel(option));
}

function searchableElements(selectId) {
  const select = el(selectId);
  const root = document.querySelector(`[data-searchable-for="${selectId}"]`);
  const input = root?.querySelector(".searchable-input");
  const list = root?.querySelector(".searchable-list");
  return { select, root, input, list };
}

function syncSearchableControl(selectId) {
  const { select, input, list } = searchableElements(selectId);
  if (!select || !input) return;
  const control = state.searchableControls.get(selectId);
  if (control?.open) {
    openSearchableList(selectId, input.value);
    return;
  }
  const selected = [...select.options].find((option) => option.value === select.value) || select.options[select.selectedIndex];
  input.value = searchableLabel(selected);
  input.title = input.value;
  input.setAttribute("aria-expanded", "false");
  if (list) list.hidden = true;
}

function closeSearchableList(selectId, { restore = true } = {}) {
  const { input, list } = searchableElements(selectId);
  const control = state.searchableControls.get(selectId);
  if (control) {
    control.open = false;
    control.activeIndex = -1;
  }
  if (list) {
    list.hidden = true;
    list.innerHTML = "";
  }
  if (input) input.setAttribute("aria-expanded", "false");
  if (restore) syncSearchableControl(selectId);
}

function searchableMatches(option, query) {
  if (!query) return true;
  const text = `${option.value} ${searchableLabel(option)}`.toLowerCase();
  return text.includes(query);
}

function searchableHighlightedLabel(option, query) {
  const label = searchableLabel(option);
  const needle = String(query || "").trim().toLowerCase();
  if (!needle) return escapeHtml(label);
  const haystack = label.toLowerCase();
  let cursor = 0;
  let html = "";
  while (cursor < label.length) {
    const index = haystack.indexOf(needle, cursor);
    if (index < 0) {
      html += escapeHtml(label.slice(cursor));
      break;
    }
    html += escapeHtml(label.slice(cursor, index));
    const end = index + needle.length;
    html += `<span class="searchable-match">${escapeHtml(label.slice(index, end))}</span>`;
    cursor = end;
  }
  return html;
}

function chooseSearchableOption(selectId, value) {
  const { select } = searchableElements(selectId);
  if (!select || select.value === value && value === "") return;
  select.value = value;
  closeSearchableList(selectId, { restore: true });
  select.dispatchEvent(new Event("input", { bubbles: true }));
  select.dispatchEvent(new Event("change", { bubbles: true }));
}

function setSearchableActiveOption(selectId, nextIndex) {
  const { list } = searchableElements(selectId);
  const control = state.searchableControls.get(selectId);
  if (!list || !control) return;
  const options = [...list.querySelectorAll(".searchable-option")];
  if (!options.length) {
    control.activeIndex = -1;
    return;
  }
  control.activeIndex = Math.max(0, Math.min(nextIndex, options.length - 1));
  options.forEach((option, idx) => option.classList.toggle("active", idx === control.activeIndex));
  options[control.activeIndex]?.scrollIntoView({ block: "nearest" });
}

function openSearchableList(selectId, rawQuery = "") {
  const { select, input, list } = searchableElements(selectId);
  if (!select || !input || !list) return;
  const control = state.searchableControls.get(selectId) || { open: false, activeIndex: -1 };
  state.searchableControls.set(selectId, control);
  const query = String(rawQuery || "").trim().toLowerCase();
  const matches = searchableOptions(select)
    .filter((option) => searchableMatches(option, query))
    .slice(0, SEARCHABLE_OPTION_LIMIT);
  control.open = true;
  control.activeIndex = matches.length ? 0 : -1;
  input.setAttribute("aria-expanded", "true");
  list.hidden = false;
  list.innerHTML = matches.length
    ? matches.map((option, idx) => `<button type="button" class="searchable-option${idx === 0 ? " active" : ""}" role="option" data-value="${escapeHtml(option.value)}">${searchableHighlightedLabel(option, query)}</button>`).join("")
    : `<div class="searchable-empty">No matches</div>`;
  list.querySelectorAll(".searchable-option").forEach((option) => {
    option.addEventListener("mousedown", (event) => event.preventDefault());
    option.addEventListener("click", () => chooseSearchableOption(selectId, option.dataset.value || ""));
  });
}

function initSearchableControl(selectId) {
  const { input, list } = searchableElements(selectId);
  if (!input || !list || input.dataset.searchableReady === "true") return;
  input.dataset.searchableReady = "true";
  state.searchableControls.set(selectId, { open: false, activeIndex: -1 });
  input.addEventListener("focus", () => {
    const control = state.searchableControls.get(selectId);
    if (!control?.open) input.value = "";
    openSearchableList(selectId, input.value);
  });
  input.addEventListener("click", () => {
    const control = state.searchableControls.get(selectId);
    if (!control?.open) input.value = "";
    openSearchableList(selectId, input.value);
  });
  input.addEventListener("input", () => openSearchableList(selectId, input.value));
  input.addEventListener("keydown", (event) => {
    const control = state.searchableControls.get(selectId);
    const options = [...list.querySelectorAll(".searchable-option")];
    if (event.key === "ArrowDown") {
      event.preventDefault();
      if (!control?.open) openSearchableList(selectId, input.value);
      else setSearchableActiveOption(selectId, (control.activeIndex ?? -1) + 1);
    } else if (event.key === "ArrowUp") {
      event.preventDefault();
      setSearchableActiveOption(selectId, (control?.activeIndex ?? 0) - 1);
    } else if (event.key === "Enter") {
      if (control?.open && options[control.activeIndex]) {
        event.preventDefault();
        chooseSearchableOption(selectId, options[control.activeIndex].dataset.value || "");
      }
    } else if (event.key === "Escape") {
      closeSearchableList(selectId, { restore: true });
      input.blur();
    }
  });
  input.addEventListener("blur", () => {
    window.setTimeout(() => closeSearchableList(selectId, { restore: true }), 120);
  });
  syncSearchableControl(selectId);
}

function customGeneSet(value) {
  return new Set(
    String(value || "")
      .split(",")
      .map((g) => g.trim().toUpperCase())
      .filter(Boolean)
  );
}

function isCancerGene(gene) {
  return state.cancerSet.has(String(gene || "").trim().toUpperCase());
}

function compoundTypeLabel(rawType, drug = "") {
  if (String(drug).toUpperCase().includes("SRPKIN-1")) return "Sulfonyl Fluorides";
  return (
    {
      FS: "Fluorosulfates",
      SuFEx: "Sulfonyl Fluorides",
      SuTEx: "Sulfonyl Triazoles",
      "Sulfonyl Fluorides": "Fluorosulfates",
      Fluorosulfates: "Sulfonyl Fluorides",
      "Sulfonyl Triazoles": "Sulfonyl Triazoles",
    }[rawType] || "Sulfonyl Fluorides"
  );
}

function compoundTypeLabelSingular(rawType, drug = "") {
  return (
    {
      "Sulfonyl Fluorides": "Sulfonyl fluoride",
      Fluorosulfates: "Fluorosulfate",
      "Sulfonyl Triazoles": "Sulfonyl triazole",
    }[compoundTypeLabel(rawType, drug)] || compoundTypeLabel(rawType, drug)
  );
}

function drugHitCount(drug) {
  return Number(state.dataset?.drugHitCounts?.[drug] ?? 0);
}

function compoundLabel(drug, hitCountOverride = null) {
  const choice = state.dataset.drugChoices[drug];
  const rawType = state.manifest?.compoundTypes?.[drug] || choice?.match(/\(([^,]+),/)?.[1] || "";
  const hitCount = hitCountOverride ?? drugHitCount(drug);
  return `${drug} (${compoundTypeLabelSingular(rawType, drug)}, ${hitCount} site${hitCount === 1 ? "" : "s"} hit)`;
}

function numericValue(id, fallback = 0) {
  const value = Number(el(id).value);
  return Number.isFinite(value) ? value : fallback;
}

function roundLabel(value, digits = 2) {
  return value == null || !Number.isFinite(Number(value)) ? "N/A" : Number(value).toFixed(digits);
}

function hitCountLimit(id) {
  const raw = el(id)?.value;
  return raw === "all" ? Infinity : Number(raw);
}

function formatPValue(value) {
  const num = Number(value);
  if (!Number.isFinite(num)) return "N/A";
  if (num === 0) return "0.00e+0";
  return num < 0.00001 ? num.toExponential(2) : num.toFixed(5);
}

function cleanDescription(value) {
  return String(value || "").split(/\s+OS=/)[0].trim();
}

function compactPValue(value) {
  return formatPValue(value);
}

function contactLabel(type) {
  return CONTACT_LABELS[type] || type;
}

function siteFilterIds() {
  return {
    sigOnly: "siteSigOnly",
    hideVariance: "siteHideVariance",
    minSn: "siteMinSn",
  };
}

function viewFilterIds(view) {
  return {
    sigOnly: `${view}SigOnly`,
    hideVariance: `${view}HideVariance`,
    minSn: `${view}MinSn`,
  };
}

function copySharedFilterValues(sourceView) {
  const source = viewFilterIds(sourceView);
  for (const view of ["compound", "site", "bio"]) {
    if (view === sourceView) continue;
    const target = viewFilterIds(view);
    if (el(target.sigOnly)) el(target.sigOnly).checked = el(source.sigOnly).checked;
    if (el(target.hideVariance)) el(target.hideVariance).checked = el(source.hideVariance).checked;
    if (el(target.minSn)) el(target.minSn).value = el(source.minSn).value;
  }
}

function copySelectivityValue(sourceId, targetId) {
  if (el(sourceId) && el(targetId)) el(targetId).value = el(sourceId).value;
}

function checkboxChecked(id) {
  return Boolean(id && el(id)?.checked);
}

function compoundQualityFiltersActive() {
  return checkboxChecked("compoundSigOnly") || checkboxChecked("compoundHideVariance") || numericValue("compoundMinSn", 0) > 0;
}

function invalidateCompoundChoiceCounts() {
  state.compoundChoiceCountCache.clear();
}

function yieldToBrowser() {
  return new Promise((resolve) => {
    if (typeof window !== "undefined" && typeof window.requestIdleCallback === "function") {
      window.requestIdleCallback(resolve, { timeout: 50 });
      return;
    }
    setTimeout(resolve, 0);
  });
}

function qualityPassFromCompound(row, ids) {
  const minSn = numericValue(ids.minSn, 0);
  const sigOnly = checkboxChecked(ids.sigOnly);
  const hideVariance = checkboxChecked(ids.hideVariance);
  const hideDmso = hideVariance || checkboxChecked(ids.hideDmso);
  const hideCpd = hideVariance || checkboxChecked(ids.hideCpd);
  if (sigOnly && row[4] <= SIG_NEGLOG10) return false;
  if (hideDmso && row[6]) return false;
  if (hideCpd && row[7]) return false;
  if ((row[5] ?? 0) < minSn) return false;
  return true;
}

function qualityPassFromCompoundSnapshot(row, filters) {
  const hideDmso = filters.hideVariance;
  const hideCpd = filters.hideVariance;
  if (filters.sigOnly && row[4] <= SIG_NEGLOG10) return false;
  if (hideDmso && row[6]) return false;
  if (hideCpd && row[7]) return false;
  if ((row[5] ?? 0) < filters.minSn) return false;
  return true;
}

function qualityPassFromHit(hit, ids) {
  const minSn = numericValue(ids.minSn, 0);
  const sigOnly = checkboxChecked(ids.sigOnly);
  const hideVariance = checkboxChecked(ids.hideVariance);
  const hideDmso = hideVariance || checkboxChecked(ids.hideDmso);
  const hideCpd = hideVariance || checkboxChecked(ids.hideCpd);
  if (sigOnly && hit[3] <= SIG_NEGLOG10) return false;
  if (hideDmso && hit[6]) return false;
  if (hideCpd && hit[7]) return false;
  if ((hit[5] ?? 0) < minSn) return false;
  return true;
}

function plotMessage(target, message) {
  const node = el(target);
  node.innerHTML = `<div class="muted">${message}</div>`;
}

function safePlot(target, traces, layout, config = {}) {
  const node = typeof target === "string" ? el(target) : target;
  if (!window.Plotly) {
    plotMessage(node?.id || target, "Plotly failed to load. Check network access to the CDN.");
    return false;
  }
  Plotly.react(node, traces, { autosize: true, ...layout }, {
    responsive: true,
    displaylogo: false,
    ...config,
  });
  schedulePlotResize(node);
  return true;
}

function schedulePlotResize(target) {
  const node = typeof target === "string" ? el(target) : target;
  if (!node || !window.Plotly?.Plots?.resize) return;
  [0, 80, 240].forEach((delay) => {
    window.setTimeout(() => Plotly.Plots.resize(node), delay);
  });
}

function baseLayout(extra = {}) {
  return {
    margin: { l: 48, r: 18, t: 18, b: 48 },
    paper_bgcolor: "white",
    plot_bgcolor: "white",
    hovermode: "closest",
    font: { family: "Helvetica, Arial, sans-serif", size: 12 },
    xaxis: { showgrid: false, zeroline: false, linecolor: "#cfcac1" },
    yaxis: { showgrid: false, zeroline: false, linecolor: "#cfcac1" },
    ...extra,
  };
}

function mergeAxes(layout, xaxis = {}, yaxis = {}) {
  return {
    ...layout,
    xaxis: { ...(layout.xaxis || {}), ...xaxis },
    yaxis: { ...(layout.yaxis || {}), ...yaxis },
  };
}

function setMolecule(targetId, drug) {
  const target = el(targetId);
  if (!drug) {
    target.innerHTML = `<p class="muted">Hover over a bar to view the compound structure.</p>`;
    return;
  }
  const path = `structures/${drug}.png`;
  target.innerHTML = `
    <strong>${drug}</strong>
    <img src="${path}" alt="${drug} structure" onerror="this.style.display='none'; this.nextElementSibling.style.display='block'">
    <p class="muted" style="display:none">Structure image not found.</p>
  `;
}

async function loadRawHover(drug) {
  const key = `${state.datasetKey}:${drug}`;
  if (state.rawHoverCache.has(key)) return state.rawHoverCache.get(key);
  const path = state.dataset?.rawHoverParts?.[drug];
  if (!path) {
    state.rawHoverCache.set(key, null);
    return null;
  }
  try {
    const payload = await fetchJson(path);
    state.rawHoverCache.set(key, payload?.rows || {});
    return payload?.rows || {};
  } catch {
    state.rawHoverCache.set(key, null);
    return null;
  }
}

function currentTooltipRowId(fallback = null) {
  if (Number.isFinite(Number(fallback))) return Number(fallback);
  if (el("site")?.classList.contains("active") && state.currentSiteRow?.i != null) return Number(state.currentSiteRow.i);
  return null;
}

function rawHoverChartHtml(values) {
  if (!values) return `<div class="raw-hover-empty muted">Raw intensities unavailable.</div>`;
  const normalized = Array.isArray(values?.[0])
    ? { dmso: values[0], compound: values[1] }
    : { dmso: [values[0], values[1]], compound: [values[2], values[3]] };
  const dmsoVals = (normalized.dmso || []).map((v) => Number(v)).filter(Number.isFinite);
  const compoundVals = (normalized.compound || []).map((v) => Number(v)).filter(Number.isFinite);
  if (!dmsoVals.length || !compoundVals.length) return `<div class="raw-hover-empty muted">Raw intensities unavailable.</div>`;
  const groups = [
    { label: "DMSO", vals: dmsoVals, color: MORANDI.gray },
    { label: "Compound", vals: compoundVals, color: MORANDI.blueDark },
  ];
  const maxVal = Math.max(...groups.flatMap((g) => g.vals).filter((v) => Number.isFinite(v)), 1);
  const axisMax = Math.max(10, Math.ceil(maxVal / 10) * 10);
  const axisMid = axisMax / 2;
  return `
    <div class="raw-hover">
      <div class="raw-hover-title">TMT SN</div>
      <div class="raw-horizontal">
        ${groups
          .map((group) => {
            const mean = group.vals.reduce((a, b) => a + b, 0) / group.vals.length;
            const meanPct = Math.max(3, Math.min(100, (mean / axisMax) * 100));
            const dots = group.vals
              .map((value) => {
                const left = Math.max(0, Math.min(100, (value / axisMax) * 100));
                return `<span class="raw-dot" style="left:${left}%; background:${group.color}" title="${roundLabel(value, 2)}"></span>`;
              })
              .join("");
            return `
              <div class="raw-row">
                <div class="raw-label">${group.label}</div>
                <div class="raw-track">
                  <div class="raw-bar" style="width:${meanPct}%; background:${group.color}"></div>
                  ${dots}
                </div>
              </div>
            `;
          })
          .join("")}
        <div class="raw-axis"><span>0</span><span>${axisMid}</span><span>${axisMax}</span></div>
      </div>
    </div>
  `;
}

function rawHoverUnavailableHtml() {
  return `<div class="raw-hover-empty muted">Raw DMSO/compound SN unavailable for this site.</div>`;
}

function rawHoverCandidateKeys(rowId) {
  const numericRowId = Number(rowId);
  if (!Number.isFinite(numericRowId)) return [];
  if (state.rawHoverAliasCache.has(numericRowId)) return state.rawHoverAliasCache.get(numericRowId);

  const catalogRow = state.rowById.get(numericRowId);
  const keys = new Set([String(numericRowId)]);
  if (catalogRow) {
    const positions = new Set((catalogRow.positions || []).map((pos) => String(pos)));
    for (const row of state.catalog) {
      if (row.i === numericRowId) continue;
      if (row.gene !== catalogRow.gene || row.proteinId !== catalogRow.proteinId || row.sequence !== catalogRow.sequence) continue;
      if ((row.positions || []).some((pos) => positions.has(String(pos)))) {
        keys.add(String(row.i));
      }
    }
  }

  const candidates = [...keys];
  state.rawHoverAliasCache.set(numericRowId, candidates);
  return candidates;
}

async function hydrateHitTooltip(drug, rowId) {
  const tooltip = el("moleculeTooltip");
  const key = `${drug}|${rowId}`;
  tooltip.dataset.key = key;
  const rows = await loadRawHover(drug);
  if (tooltip.dataset.key !== key || tooltip.style.display === "none") return;
  const slot = tooltip.querySelector(".raw-hover-slot");
  if (!slot) return;
  const candidateKeys = rawHoverCandidateKeys(rowId);
  const values = candidateKeys.map((candidate) => rows?.[candidate]).find(Boolean);
  slot.innerHTML = values ? rawHoverChartHtml(values) : rawHoverUnavailableHtml();
  repositionTooltipFromDataset(tooltip);
}

async function hydratePointTooltip(drug, rowId) {
  const tooltip = el("moleculeTooltip");
  const key = `${drug}|point|${rowId}`;
  tooltip.dataset.key = key;
  const rows = await loadRawHover(drug);
  if (tooltip.dataset.key !== key || tooltip.style.display === "none") return;
  const slot = tooltip.querySelector(".raw-hover-slot");
  if (!slot) return;
  const candidateKeys = rawHoverCandidateKeys(rowId);
  const values = candidateKeys.map((candidate) => rows?.[candidate]).find(Boolean);
  slot.innerHTML = values ? rawHoverChartHtml(values) : rawHoverUnavailableHtml();
  repositionTooltipFromDataset(tooltip);
}

function showMoleculeTooltip(drug, event) {
  const tooltip = el("moleculeTooltip");
  tooltip.classList.remove("compact-tooltip", "hit-tooltip");
  if (!drug || !event?.event) {
    tooltip.style.display = "none";
    return;
  }
  tooltip.innerHTML = `
    <strong>${drug}</strong>
    <img src="structures/${drug}.png" alt="${drug} structure" onerror="this.style.display='none'; this.nextElementSibling.style.display='block'">
    <p class="muted" style="display:none">Structure image not found.</p>
  `;
  tooltip.style.display = "block";
  positionTooltip(tooltip, event.event, 250, 250);
}

function computeTooltipPosition(tooltip, pointerEvent, fallbackWidth = 286, fallbackHeight = 250) {
  const rect = tooltip.getBoundingClientRect();
  const width = Math.max(rect.width || fallbackWidth, fallbackWidth);
  const height = Math.max(tooltip.scrollHeight || rect.height || fallbackHeight, fallbackHeight);
  const edgePad = 32;
  const pointerPad = 24;
  const bottomReserve = 92;
  let x = pointerEvent.clientX + pointerPad;
  let y = pointerEvent.clientY + pointerPad;
  if (x + width + edgePad > window.innerWidth) x = pointerEvent.clientX - width - pointerPad;
  if (y + height + bottomReserve > window.innerHeight) y = pointerEvent.clientY - height - pointerPad;
  return {
    x: Math.max(edgePad, Math.min(x, window.innerWidth - width - edgePad)),
    y: Math.max(edgePad, Math.min(y, window.innerHeight - height - bottomReserve)),
  };
}

function positionTooltip(tooltip, pointerEvent, fallbackWidth = 286, fallbackHeight = 250) {
  tooltip.dataset.pointerX = String(pointerEvent.clientX);
  tooltip.dataset.pointerY = String(pointerEvent.clientY);
  const pos = computeTooltipPosition(tooltip, pointerEvent, fallbackWidth, fallbackHeight);
  tooltip.style.left = `${pos.x}px`;
  tooltip.style.top = `${pos.y}px`;
}

function repositionTooltipFromDataset(tooltip) {
  const clientX = Number(tooltip.dataset.pointerX);
  const clientY = Number(tooltip.dataset.pointerY);
  if (!Number.isFinite(clientX) || !Number.isFinite(clientY)) return;
  positionTooltip(tooltip, { clientX, clientY }, 286, 300);
}

function showHitTooltip(payload, event) {
  if (!payload) {
    hideMoleculeTooltip();
    return;
  }
  const tooltip = el("moleculeTooltip");
  tooltip.classList.remove("compact-tooltip");
  tooltip.classList.add("hit-tooltip");
  const hit = payload.hit || payload;
  const drug = hit[0];
  const rowId = currentTooltipRowId(payload.rowId);
  tooltip.innerHTML = `
    <strong>${drug}</strong>
    <img src="structures/${drug}.png" alt="${drug} structure" onerror="this.style.display='none'; this.nextElementSibling.style.display='block'">
    <p class="muted" style="display:none">Structure image not found.</p>
    <div class="tooltip-stats one-row">
      <div><strong>R:</strong><span>${roundLabel(hit[1])}</span></div>
      <div><strong>P-value:</strong><span>${formatPValue(hit[2])}</span></div>
      <div><strong>Sites Hit:</strong><span>${drugHitCount(drug)}</span></div>
    </div>
    <div class="raw-hover-slot"><div class="raw-hover-empty muted">Loading raw intensities...</div></div>
  `;
  tooltip.style.display = "block";
  positionTooltip(tooltip, event.event, 286, 320);
  hydrateHitTooltip(drug, rowId);
}

function showTransientHitTooltip(payload, event) {
  if (state.persistentHitTooltip) return;
  showHitTooltip(payload, event);
}

function hideTransientHitTooltip() {
  if (state.persistentHitTooltip) return;
  hideMoleculeTooltip();
}

function ensureCopyToast() {
  let toast = document.querySelector(".copy-toast");
  if (toast) return toast;
  toast = document.createElement("div");
  toast.className = "copy-toast";
  toast.setAttribute("role", "status");
  toast.setAttribute("aria-live", "polite");
  document.body.appendChild(toast);
  return toast;
}

function showCopyToast(message = "Compound ID copied to clipboard") {
  const toast = ensureCopyToast();
  toast.textContent = message;
  toast.classList.add("show");
  if (state.copyToastTimer) window.clearTimeout(state.copyToastTimer);
  state.copyToastTimer = window.setTimeout(() => {
    toast.classList.remove("show");
    state.copyToastTimer = null;
  }, 1600);
}

function fallbackCopyText(text) {
  const textarea = document.createElement("textarea");
  textarea.value = text;
  textarea.setAttribute("readonly", "");
  textarea.style.position = "fixed";
  textarea.style.left = "-9999px";
  textarea.style.top = "0";
  document.body.appendChild(textarea);
  textarea.focus();
  textarea.select();
  textarea.setSelectionRange(0, textarea.value.length);
  let copied = false;
  try {
    copied = document.execCommand("copy");
  } catch {
    copied = false;
  }
  textarea.remove();
  return copied;
}

function copyPersistentHitCompoundId(payload) {
  const hit = payload?.hit || payload;
  const drug = hit?.[0] ? String(hit[0]) : "";
  if (!drug) return;
  const fallbackCopied = fallbackCopyText(drug);
  if (fallbackCopied) showCopyToast();
  const write = window.navigator?.clipboard?.writeText?.(drug);
  if (!write) {
    if (!fallbackCopied) showCopyToast("Unable to copy compound ID");
    return;
  }
  write
    .then(() => {
      if (!fallbackCopied) showCopyToast();
    })
    .catch(() => {
      if (!fallbackCopied) showCopyToast("Unable to copy compound ID");
    });
}

function pinHitTooltip(payload, event) {
  if (!payload) return;
  event?.event?.stopPropagation?.();
  state.persistentHitTooltip = { payload };
  copyPersistentHitCompoundId(payload);
  showHitTooltip(payload, event);
}

function clearPersistentHitTooltip({ hide = true } = {}) {
  if (!state.persistentHitTooltip) return;
  state.persistentHitTooltip = null;
  if (hide) hideMoleculeTooltip();
}

function showCompoundPointTooltip(payload, event) {
  if (!payload || !event?.event) {
    hideMoleculeTooltip();
    return;
  }
  const drug = el("compoundSelect").value;
  const tooltip = el("moleculeTooltip");
  tooltip.classList.remove("compact-tooltip", "hit-tooltip");
  const description = cleanDescription(payload.description);
  tooltip.innerHTML = `
    <strong>${payload.label}</strong>
    <div class="tooltip-uniprot">${payload.uniprot || "N/A"}</div>
    <div class="tooltip-detail">${description || "No description available."}</div>
    <div class="tooltip-stats one-row">
      <div><strong>R:</strong><span>${roundLabel(payload.r)}</span></div>
      <div><strong>P-value:</strong><span>${compactPValue(payload.p)}</span></div>
      <div><strong>Site rank:</strong><span>${payload.rRank ?? "N/A"}/${payload.totalSites ?? "N/A"}</span></div>
    </div>
    <div class="raw-hover-slot"><div class="raw-hover-empty muted">Loading raw intensities...</div></div>
  `;
  tooltip.style.display = "block";
  positionTooltip(tooltip, event.event, 286, 320);
  hydratePointTooltip(drug, payload.rowId);
}

function showSiteDistributionTooltip(payload, event) {
  if (!payload || !event?.event) {
    hideMoleculeTooltip();
    return;
  }
  const tooltip = el("moleculeTooltip");
  tooltip.classList.remove("hit-tooltip");
  tooltip.classList.add("compact-tooltip");
  tooltip.innerHTML = `
    <div class="tooltip-compact-value">${escapeHtml(payload.title)} (${Number(payload.sites?.length || 0).toLocaleString()} sites)</div>
    <div class="tooltip-detail tooltip-click-hint"><em>Click to see all sites</em></div>
  `;
  tooltip.style.display = "block";
  positionTooltip(tooltip, event.event, 0, 92);
}

function showSiteHitsTooltip(payload, event) {
  if (!payload || !event?.event) {
    hideMoleculeTooltip();
    return;
  }
  const tooltip = el("moleculeTooltip");
  tooltip.classList.remove("hit-tooltip");
  tooltip.classList.add("compact-tooltip");
  tooltip.innerHTML = `
    <div class="tooltip-compact-value">${Number(payload.count || 0).toLocaleString()} sites (${roundLabel(payload.percent, 1)}%)</div>
  `;
  tooltip.style.display = "block";
  positionTooltip(tooltip, event.event, 0, 72);
}

function showCompoundSummaryTooltip(payload, event) {
  if (!payload || !event?.event) {
    hideMoleculeTooltip();
    return;
  }
  const tooltip = el("moleculeTooltip");
  tooltip.classList.remove("hit-tooltip");
  tooltip.classList.add("compact-tooltip");
  tooltip.innerHTML = `
    <div class="tooltip-compact-value">${Number(payload.drugs?.length || 0).toLocaleString()} compounds</div>
    <div class="tooltip-detail tooltip-click-hint"><em>Click to see all compounds</em></div>
  `;
  tooltip.style.display = "block";
  positionTooltip(tooltip, event.event, 0, 92);
}

function hideMoleculeTooltip() {
  const tooltip = el("moleculeTooltip");
  tooltip.style.display = "none";
  tooltip.classList.remove("compact-tooltip", "hit-tooltip");
}

function closeSummaryDetailModal() {
  const modal = el("summaryDetailModal");
  if (!modal) return;
  modal.classList.remove("open");
  modal.classList.remove("circular-heatmap-modal");
  modal.setAttribute("aria-hidden", "true");
}

function openSummaryDetailModal(title, meta, bodyHtml) {
  hideMoleculeTooltip();
  el("summaryDetailTitle").textContent = title;
  el("summaryDetailMeta").textContent = meta || "";
  el("summaryDetailBody").innerHTML = bodyHtml || '<p class="muted">No matching records.</p>';
  const modal = el("summaryDetailModal");
  modal.classList.add("open");
  modal.setAttribute("aria-hidden", "false");
}

async function openCircularHeatmapModal(title, meta, rows, compounds, cells) {
  const clustered = await clusterCircularHeatmapSites(rows, cells);
  const clusterMeta = clustered.status ? `${meta || ""}${meta ? " · " : ""}${clustered.status}` : meta;
  openSummaryDetailModal(title, meta, '<canvas id="circularFlatHeatmap" aria-label="Flat compound by biology heatmap"></canvas>');
  el("summaryDetailMeta").textContent = clusterMeta || "";
  el("summaryDetailModal").classList.add("circular-heatmap-modal");
  drawFlatHeatmap(el("circularFlatHeatmap"), clustered.rows, compounds, clustered.cells, { dendrogram: clustered.dendrogram });
}

async function circularSiteHit(drug, rowId) {
  const numericRowId = Number(rowId);
  if (!Number.isFinite(numericRowId)) return null;
  const row = state.rowById.get(numericRowId);
  const summary = await siteSummaryForRow(row);
  return (summary?.hits || []).find((hit) => String(hit[0]) === String(drug)) || null;
}

async function hydrateCircularSiteTooltip(drug, rowId) {
  const tooltip = el("moleculeTooltip");
  const key = `${drug}|circular-site|${rowId}`;
  tooltip.dataset.key = key;
  const [rows, siteHit] = await Promise.all([loadRawHover(drug), circularSiteHit(drug, rowId)]);
  if (tooltip.dataset.key !== key || tooltip.style.display === "none") return;
  const rValue = tooltip.querySelector("[data-circular-r]");
  const pValue = tooltip.querySelector("[data-circular-pvalue]");
  if (siteHit) {
    if (rValue) rValue.textContent = roundLabel(siteHit[1]);
    if (pValue) pValue.textContent = formatPValue(siteHit[2]);
  }
  const slot = tooltip.querySelector(".raw-hover-slot");
  if (!slot) return;
  const candidateKeys = rawHoverCandidateKeys(rowId);
  const values = candidateKeys.map((candidate) => rows?.[candidate]).find(Boolean);
  slot.innerHTML = values ? rawHoverChartHtml(values) : rawHoverUnavailableHtml();
  repositionTooltipFromDataset(tooltip);
}

function showCircularAggregateTooltip(payload, event) {
  if (!payload?.hit || !event?.event) {
    hideMoleculeTooltip();
    return;
  }
  const tooltip = el("moleculeTooltip");
  const hit = payload.hit;
  const drug = hit[0];
  tooltip.classList.remove("compact-tooltip");
  tooltip.classList.add("hit-tooltip");
  tooltip.innerHTML = `
    <strong>${escapeHtml(drug)}</strong>
    <img src="structures/${escapeHtml(drug)}.png" alt="${escapeHtml(drug)} structure" onerror="this.style.display='none'; this.nextElementSibling.style.display='block'">
    <p class="muted" style="display:none">Structure image not found.</p>
    <div class="tooltip-stats one-row">
      <div><strong>Max R:</strong><span>${roundLabel(payload.maxR ?? hit[1])}</span></div>
    </div>
  `;
  tooltip.style.display = "block";
  positionTooltip(tooltip, event.event, 286, 250);
}

function showCircularSiteTooltip(payload, event) {
  if (!payload?.hit || !event?.event) {
    hideMoleculeTooltip();
    return;
  }
  const tooltip = el("moleculeTooltip");
  const hit = payload.hit;
  const drug = hit[0];
  tooltip.classList.remove("compact-tooltip");
  tooltip.classList.add("hit-tooltip");
  tooltip.innerHTML = `
    <strong>${escapeHtml(drug)}</strong>
    <img src="structures/${escapeHtml(drug)}.png" alt="${escapeHtml(drug)} structure" onerror="this.style.display='none'; this.nextElementSibling.style.display='block'">
    <p class="muted" style="display:none">Structure image not found.</p>
    <div class="tooltip-stats one-row">
      <div><strong>R:</strong><span data-circular-r>${roundLabel(hit[1])}</span></div>
      <div><strong>P-value:</strong><span data-circular-pvalue>${formatPValue(hit[2])}</span></div>
      <div><strong>Sites Hit:</strong><span>${drugHitCount(drug)}</span></div>
    </div>
    <div class="raw-hover-slot"><div class="raw-hover-empty muted">Loading raw intensities...</div></div>
  `;
  tooltip.style.display = "block";
  positionTooltip(tooltip, event.event, 286, 320);
  hydrateCircularSiteTooltip(drug, payload.rowId);
}

function showCircularHeatmapTooltip(payload, event) {
  if (payload?.mode === "site") {
    showCircularSiteTooltip(payload, event);
    return;
  }
  showCircularAggregateTooltip(payload, event);
}

async function openCircularSiteInSiteView(row) {
  const siteRow = Number(row?.siteRow ?? row?.rowId);
  if (!Number.isFinite(siteRow) || siteRow < 0) return;
  const catalogRow = state.rowById.get(siteRow);
  if (!catalogRow?.gene) return;
  closeSummaryDetailModal();
  switchTab("site");
  if (el("siteGeneSelect").value !== catalogRow.gene) {
    await refreshSiteControls(catalogRow.gene, { allowFallback: true });
  }
  el("siteSelect").value = String(siteRow);
  await renderSite();
}

function compoundCardHtml(drug) {
  return `
    <article class="compound-card">
      <img src="structures/${escapeHtml(drug)}.png" alt="${escapeHtml(drug)} structure" onerror="this.style.display='none'">
      <strong>${escapeHtml(drug)}</strong>
      <span class="muted">${escapeHtml(compoundTypeLabelSingular(state.manifest?.compoundTypes?.[drug], drug))}</span>
    </article>
  `;
}

function showCompoundBinModal(payload) {
  const drugs = payload?.drugs || [];
  const title = payload?.title || "Compounds";
  const meta = `${drugs.length.toLocaleString()} compound${drugs.length === 1 ? "" : "s"}`;
  const body = drugs.length ? `<div class="compound-card-grid">${drugs.map(compoundCardHtml).join("")}</div>` : "";
  openSummaryDetailModal(title, meta, body);
}

function showSiteBinModal(payload) {
  const sites = payload?.sites || [];
  const rows = sites
    .map((site) => {
      const hitCount = Number(site.hitCount || 0);
      const observedCount = Number(site.observedCount || 0);
      return `<li><strong>${escapeHtml(site.label)}</strong><span>${roundLabel(site.value)}%, ${hitCount.toLocaleString()}/${observedCount.toLocaleString()}</span></li>`;
    })
    .join("");
  openSummaryDetailModal(
    payload?.title || "Reactive Sites",
    `${sites.length.toLocaleString()} site${sites.length === 1 ? "" : "s"}`,
    rows ? `<ul class="site-bin-list">${rows}</ul>` : ""
  );
}

function bindSiteDistributionHover() {
  const plot = el("siteReactivityDist");
  if (typeof plot.on !== "function") return;
  plot.removeAllListeners?.("plotly_hover");
  plot.removeAllListeners?.("plotly_unhover");
  plot.removeAllListeners?.("plotly_click");
  plot.on("plotly_hover", (event) => showSiteDistributionTooltip(event.points[0].customdata, event));
  plot.on("plotly_unhover", hideMoleculeTooltip);
  plot.on("plotly_click", (event) => showSiteBinModal(event.points[0].customdata));
}

function bindSiteHitsHover() {
  const plot = el("siteReactivityBinary");
  const payloads = plot._siteHitPayloads || [];
  const slices = [...plot.querySelectorAll(".slice")];
  slices.forEach((slice, idx) => {
    const payload = payloads[idx];
    const target = slice.querySelector(".surface") || slice;
    target.onmouseenter = (event) => showSiteHitsTooltip(payload, { event });
    target.onmousemove = (event) => showSiteHitsTooltip(payload, { event });
    target.onmouseleave = hideMoleculeTooltip;
  });
}

function bindCompoundSummaryClick(targetId) {
  const plot = el(targetId);
  if (typeof plot.on !== "function") return;
  plot.removeAllListeners?.("plotly_hover");
  plot.removeAllListeners?.("plotly_unhover");
  plot.removeAllListeners?.("plotly_click");
  plot.on("plotly_hover", (event) => showCompoundSummaryTooltip(event.points[0].customdata, event));
  plot.on("plotly_unhover", hideMoleculeTooltip);
  plot.on("plotly_click", (event) => showCompoundBinModal(event.points[0].customdata));
}

function bindBarHover(targetId, payloads) {
  const plot = el(targetId);
  const nodes = [...plot.querySelectorAll(".barlayer .bars .point")];
  nodes.forEach((node, idx) => {
    const payload = payloads[idx];
    node.onmouseenter = (event) => showTransientHitTooltip(payload, { event });
    node.onmousemove = (event) => showTransientHitTooltip(payload, { event });
    node.onclick = (event) => pinHitTooltip(payload, { event });
    node.onmouseleave = hideTransientHitTooltip;
  });
  const dragLayer = plot.querySelector(".nsewdrag");
  if (dragLayer) {
    dragLayer.onmouseleave = hideTransientHitTooltip;
  }
}

function scheduleSummaryRender() {
  if (state.summaryRenderTimer) {
    window.clearTimeout(state.summaryRenderTimer);
    state.summaryRenderTimer = null;
  }
  state.summaryRenderTimer = window.setTimeout(() => {
    state.summaryRenderTimer = null;
    if (!el("summary")?.classList.contains("active")) return;
    const render = () => renderGlobalProteome();
    if (typeof window.requestIdleCallback === "function") {
      window.requestIdleCallback(render, { timeout: 3000 });
    } else {
      window.setTimeout(render, 0);
    }
  }, SUMMARY_RENDER_DEFER_MS);
}

function cancelScheduledSummaryRenders() {
  state.summaryDatasetRenderToken++;
  if (state.summaryDatasetRenderTimer) {
    window.clearTimeout(state.summaryDatasetRenderTimer);
    state.summaryDatasetRenderTimer = null;
  }
  if (state.summaryRenderTimer) {
    window.clearTimeout(state.summaryRenderTimer);
    state.summaryRenderTimer = null;
  }
}

function scheduleSummaryDatasetRender() {
  if (state.summaryDatasetRenderTimer) {
    window.clearTimeout(state.summaryDatasetRenderTimer);
    state.summaryDatasetRenderTimer = null;
  }
  const datasetKey = state.datasetKey;
  const renderToken = ++state.summaryDatasetRenderToken;
  state.summaryDatasetRenderTimer = window.setTimeout(async () => {
    state.summaryDatasetRenderTimer = null;
    if (!el("summary")?.classList.contains("active")) return;
    if (renderToken !== state.summaryDatasetRenderToken || datasetKey !== state.datasetKey) return;
    renderDatasetStaticPlots();
    await loadCircularCompoundAtlas();
    if (renderToken !== state.summaryDatasetRenderToken || datasetKey !== state.datasetKey) return;
    renderCircularCompoundAtlas();
    await renderCompound();
    if (renderToken !== state.summaryDatasetRenderToken || datasetKey !== state.datasetKey) return;
    scheduleSummaryRender();
  }, SUMMARY_RENDER_DEFER_MS);
}

async function loadDataset(datasetKey) {
  cancelScheduledSummaryRenders();
  setDatasetLoading(true, datasetKey);
  state.datasetKey = datasetKey;
  state.dataset = state.manifest.datasets[datasetKey];
  el("datasetSwitch").checked = datasetKey === "frac";
  try {
    state.catalog = await fetchJson(state.dataset.catalog);
    state.rowById = new Map(state.catalog.map((row) => [row.i, row]));
    state.compoundCache.clear();
    state.compoundLoadPromises.clear();
    state.compoundWarmupToken++;
    invalidateCompoundChoiceCounts();
    state.geneSiteCache.clear();
    state.contactTargetRowCache.clear();
    state.contactTargetRowSetCache.clear();
    state.filteredSitesByGene.clear();
    state.rawHoverAliasCache.clear();
    state.currentSiteRow = null;
    state.currentSiteSummary = null;

    await populateDatasetControls();
    scheduleCompoundCacheWarmup();
    await refreshSiteControls();
    if (el("summary")?.classList.contains("active")) scheduleSummaryDatasetRender();
    if (el("compound")?.classList.contains("active")) await renderCompound();
    updateDatasetStatus();
  } finally {
    setDatasetLoading(false, datasetKey);
  }
}

function setDatasetLoading(isLoading, datasetKey = state.datasetKey) {
  const toggle = document.querySelector(".dataset-toggle");
  const checkbox = el("datasetSwitch");
  toggle?.classList.toggle("loading", isLoading);
  toggle?.setAttribute("aria-busy", isLoading ? "true" : "false");
  if (checkbox) checkbox.disabled = isLoading;
  if (isLoading) {
    const datasetLabel = datasetKey === "frac" ? "Fractionated" : "One-Shot";
    setStatus(`Loading ${datasetLabel} data...`);
  }
}

async function populateDatasetControls() {
  await populateCompoundChoices();

  const geneOptions = state.dataset.geneChoices.map((gene) => ({ value: gene, label: gene }));
  setOptions(el("siteGeneSelect"), geneOptions, state.dataset.geneChoices[0] || state.dataset.defaultGene);

  const hitCounts = [...new Set(Object.values(state.dataset.drugHitCounts || {}).map((v) => Number(v)).filter((v) => Number.isFinite(v) && v > 0))].sort(
    (a, b) => a - b
  );
  setOptions(
    el("siteMaxHits"),
    [{ value: "all", label: "Any" }, ...hitCounts.map((v) => ({ value: String(v), label: `≤ ${v} site${v === 1 ? "" : "s"} hit` }))],
    "all"
  );
  setOptions(el("siteTargetList"), TARGET_LIST_OPTIONS, "all");

  const highlightModes = [
    { value: "threshold", label: "Above threshold" },
    { value: "pvalue", label: "P-value gradient" },
    { value: "cancer", label: "Cancer-driver genes" },
    { value: "contact:PPI", label: "Sites at PPI interface" },
    { value: "contact:Metal", label: "Sites near metal" },
    { value: "contact:Cofactor", label: "Sites near cofactor" },
    { value: "contact:Dna", label: "Sites near DNA" },
    { value: "contact:Rna", label: "Sites near RNA" },
    { value: "contact:Ligand", label: "Sites near ligand" },
    { value: "custom", label: "Custom gene list" },
  ];
  setOptions(el("compoundColorMode"), highlightModes);
  updateConditionalFields();
}

function compoundChoiceFilterKey() {
  return [
    numericValue("compoundThreshold", 2),
    checkboxChecked("compoundSigOnly") ? 1 : 0,
    checkboxChecked("compoundHideVariance") ? 1 : 0,
    numericValue("compoundMinSn", 0),
    checkboxChecked("dropNoReplicates") ? 1 : 0,
  ].join("|");
}

async function filteredCompoundHitCount(drug) {
  const key = compoundChoiceFilterKey();
  const threshold = numericValue("compoundThreshold", 2);
  const filters = {
    sigOnly: checkboxChecked("compoundSigOnly"),
    hideVariance: checkboxChecked("compoundHideVariance"),
    minSn: numericValue("compoundMinSn", 0),
  };
  if (!state.compoundChoiceCountCache.has(key)) {
    state.compoundChoiceCountCache.set(key, new Map());
  }
  const cached = state.compoundChoiceCountCache.get(key);
  if (cached.has(drug)) return cached.get(drug);
  const rows = await loadCompound(drug);
  const count = rows.filter(
    (row) =>
      row[1] > threshold &&
      qualityPassFromCompoundSnapshot(row, filters)
  ).length;
  cached.set(drug, count);
  return count;
}

async function updateSelectedCompoundChoiceLabel(drug = el("compoundSelect").value) {
  if (!drug) return;
  const count = await filteredCompoundHitCount(drug);
  const select = el("compoundSelect");
  if (select.value !== drug) return;
  const option = [...select.options].find((candidate) => candidate.value === drug);
  if (!option) return;
  const label = compoundLabel(drug, count);
  option.textContent = label;
  option.label = label;
  syncSearchableControl("compoundSelect");
}

function setCompoundListUpdating(isUpdating) {
  const status = el("compoundListStatus");
  if (!status) return;
  status.hidden = !isUpdating;
}

async function refreshCompoundChoiceLabelsInPlace(refreshToken = state.compoundChoiceRefreshToken) {
  const select = el("compoundSelect");
  const options = [...select.options].map((option) => option.value).filter(Boolean);
  await Promise.all(options.map(async (drug) => {
    const count = await filteredCompoundHitCount(drug);
    if (refreshToken !== state.compoundChoiceRefreshToken) return;
    const option = [...select.options].find((candidate) => candidate.value === drug);
    if (!option) return;
    const label = compoundLabel(drug, count);
    option.textContent = label;
    option.label = label;
  }));
  syncSearchableControl("compoundSelect");
}

async function dynamicCompoundHitCounts(drugs, refreshToken = null) {
  const counts = new Map();
  const batchSize = 16;
  for (let i = 0; i < drugs.length; i += batchSize) {
    if (refreshToken != null && refreshToken !== state.compoundChoiceRefreshToken) return null;
    const batch = drugs.slice(i, i + batchSize);
    const entries = await Promise.all(batch.map(async (drug) => [drug, await filteredCompoundHitCount(drug)]));
    for (const [drug, count] of entries) counts.set(drug, count);
    await yieldToBrowser();
  }
  return counts;
}

function scheduleCompoundCacheWarmup() {
  const datasetKey = state.datasetKey;
  const warmupToken = ++state.compoundWarmupToken;
  setTimeout(() => warmCompoundCache(datasetKey, warmupToken), 100);
}

async function warmCompoundCache(datasetKey, warmupToken) {
  const dataset = state.manifest?.datasets?.[datasetKey];
  if (!dataset) return;
  const selected = el("compoundSelect")?.value || dataset.defaultDrug;
  const active = activeDatasetDrugs();
  const drugs = [
    selected,
    dataset.defaultDrug,
    ...active,
  ].filter(Boolean).filter((drug, idx, arr) => arr.indexOf(drug) === idx).slice(0, COMPOUND_CACHE_WARMUP_LIMIT);
  const batchSize = 12;
  for (let i = 0; i < drugs.length; i += batchSize) {
    if (warmupToken !== state.compoundWarmupToken || datasetKey !== state.datasetKey) return;
    await Promise.allSettled(drugs.slice(i, i + batchSize).map((drug) => loadCompoundFromDataset(datasetKey, dataset, drug)));
    await yieldToBrowser();
  }
}

async function populateCompoundChoices(refreshToken = null) {
  const activeOnly = checkboxChecked("compoundActiveOnly");
  const dropNoReplicates = checkboxChecked("dropNoReplicates");
  const noReplicateDrugs = new Set(state.dataset.noReplicateDrugs || []);
  const tier2OnlyDrugs = new Set(state.dataset.tier2OnlyDrugs || []);
  const current = el("compoundSelect").value || state.dataset.defaultDrug;
  let dynamicCounts = null;
  if (compoundQualityFiltersActive() || numericValue("compoundThreshold", 2) !== 2) {
    dynamicCounts = await dynamicCompoundHitCounts(state.dataset.rawDrugs, refreshToken);
    if (!dynamicCounts) return false;
  }
  const shownHitCount = (drug) => dynamicCounts?.get(drug) ?? drugHitCount(drug);
  const drugOptions = state.dataset.rawDrugs
    .filter((drug) => !dropNoReplicates || (!noReplicateDrugs.has(drug) && !tier2OnlyDrugs.has(drug)))
    .filter((drug) => !activeOnly || shownHitCount(drug) > 0)
    .sort((a, b) => shownHitCount(a) - shownHitCount(b) || compoundLabel(a, shownHitCount(a)).localeCompare(compoundLabel(b, shownHitCount(b))))
    .map((drug) => ({
      value: drug,
      label: compoundLabel(drug, shownHitCount(drug)),
    }));
  const selected = drugOptions.some((option) => option.value === current)
    ? current
    : drugOptions[0]?.value || "";
  if (refreshToken != null && refreshToken !== state.compoundChoiceRefreshToken) return false;
  setOptions(el("compoundSelect"), drugOptions.length ? drugOptions : [{ value: "", label: "No compounds found" }], selected);
  return true;
}

function scheduleCompoundChoiceListRefresh({ invalidate = false } = {}) {
  const refreshToken = ++state.compoundChoiceRefreshToken;
  setCompoundListUpdating(true);
  if (state.compoundChoiceRefreshTimer) clearTimeout(state.compoundChoiceRefreshTimer);
  state.compoundChoiceRefreshTimer = setTimeout(() => {
    state.compoundChoiceRefreshTimer = null;
    refreshCompoundChoiceList({ refreshToken, invalidate });
  }, 0);
}

async function refreshCompoundChoiceList({ refreshToken = null, invalidate = true } = {}) {
  const activeRefreshToken = refreshToken ?? ++state.compoundChoiceRefreshToken;
  setCompoundListUpdating(true);
  if (invalidate) invalidateCompoundChoiceCounts();
  try {
    const updated = await populateCompoundChoices(activeRefreshToken);
    if (updated) await renderCompound();
  } catch (err) {
    console.warn("Unable to refresh compound choice list", err);
  } finally {
    if (activeRefreshToken === state.compoundChoiceRefreshToken) setCompoundListUpdating(false);
  }
}

function updateConditionalFields() {
  el("siteCustomGenesWrap").hidden = el("siteTargetList").value !== "custom";
  el("compoundCustomGenesWrap").hidden = el("compoundColorMode").value !== "custom";
  el("globalCustomGenesWrap").hidden = !checkboxChecked("globalShowCustom");
  el("globalThresholdWrap").hidden = !checkboxChecked("globalShowHits");
}

function updateStructureLinks(structureHref = "", contactHref = "") {
  const structureLink = el("structureSourceLink");
  const contactLink = el("contactSourceLink");
  structureLink.hidden = !structureHref;
  contactLink.hidden = !contactHref;
  structureLink.href = structureHref || "#";
  contactLink.href = contactHref || "#";
}

function renderDatasetStaticPlots() {
  const summaryDrugs = activeDatasetDrugs();
  const sitePromiscuity = (state.dataset.sitePromiscuity || []).map(Number);
  const siteHits = sitePromiscuity
    .filter((v) => Number.isFinite(v))
    .map((v) => Math.round((v / 100) * summaryDrugs.length));
  const siteZero = siteHits.filter((v) => v === 0).length;
  const siteActive = siteHits.length - siteZero;
  const siteHitPayloads = [
    { label: "\u22651 Hit", count: siteActive, percent: (100 * siteActive) / Math.max(siteHits.length, 1) },
    { label: "0 Hits", count: siteZero, percent: (100 * siteZero) / Math.max(siteHits.length, 1) },
  ];
  el("siteReactivityBinary")._siteHitPayloads = siteHitPayloads;
  safePlot(
    "siteReactivityBinary",
    [
      {
        type: "pie",
        labels: ["\u22651 Hit", "0 Hits"],
        values: [siteActive, siteZero],
        customdata: siteHitPayloads,
        domain: { x: [0, 0.82], y: [0.02, 1] },
        hole: 0.45,
        sort: false,
        direction: "clockwise",
        marker: { colors: [MORANDI.red, MORANDI.gray] },
        textinfo: "none",
        hoverinfo: "none",
      },
    ],
    baseLayout({
      margin: { l: 2, r: 36, t: 4, b: 4 },
      showlegend: true,
      legend: { orientation: "v", traceorder: "normal", x: 0.82, xanchor: "left", y: 0.9, yanchor: "top" },
    })
  );
  bindSiteHitsHover();
  window.setTimeout(bindSiteHitsHover, 80);
  window.setTimeout(bindSiteHitsHover, 240);
  window.setTimeout(bindSiteHitsHover, 600);

  const siteBins = Array.from({ length: 50 }, (_, idx) => ({
    low: idx * 2,
    high: (idx + 1) * 2,
    mid: idx * 2 + 1,
    sites: [],
  }));
  sitePromiscuity.forEach((value, idx) => {
    if (!Number.isFinite(value) || value <= 0) return;
    const binIndex = Math.min(49, Math.max(0, Math.floor(value / 2)));
    const row = state.catalog[idx];
    const hitCount = Number(state.dataset.siteHitCounts?.[idx] ?? row?.hitCount ?? 0);
    const observedCount = Number(state.dataset.siteObservedCounts?.[idx] ?? row?.observedCount ?? 0);
    siteBins[binIndex].sites.push({
      label: row?.label || `Site ${idx + 1}`,
      gene: row?.gene || "",
      value,
      hitCount,
      observedCount,
    });
  });
  const activeBins = siteBins.filter((bin) => bin.sites.length);
  const plottedSiteDist = safePlot(
    "siteReactivityDist",
    [{
      type: "bar",
      x: activeBins.map((bin) => bin.mid),
      y: activeBins.map((bin) => bin.sites.length),
      width: activeBins.map(() => 1.8),
      marker: { color: MORANDI.blue },
      opacity: 0.82,
      customdata: activeBins.map((bin) => ({
        title: `${bin.low}-${bin.high}%`,
        sites: bin.sites.sort((a, b) => b.value - a.value || a.label.localeCompare(b.label)),
      })),
      hoverinfo: "none",
    }],
    mergeAxes(baseLayout({ margin: { l: 44, r: 2, t: 8, b: 54 } }), { title: "Reactivity (% compounds hit at R > 2)", range: [0, 100], automargin: true }, { title: "Number of Sites" })
  );
  if (plottedSiteDist) bindSiteDistributionHover();

  const promiscuityRecords = state.dataset.compoundPromiscuityRecords || [];
  const typeByDrug = state.manifest.compoundTypes || {};
  const normalizeType = (drug, record) => {
    return compoundTypeLabel(typeByDrug[drug] || record?.Type, drug);
  };
  const compoundRecords = summaryDrugs.map((drug) => {
    const idx = state.dataset.rawDrugs.indexOf(drug);
    const record = promiscuityRecords[idx] || {};
    const promiscuity = Number(state.dataset.drugPromiscuity?.[drug] ?? record.Promiscuity ?? 0);
    return {
      Drug: drug,
      Promiscuity: promiscuity,
      Hits: drugHitCount(drug),
      Type: normalizeType(drug, record),
    };
  });
  const statusCounts = new Map();
  const statusDrugs = new Map();
  for (const row of compoundRecords) {
    const status = row.Hits === 0 ? "0 Hits" : "\u22651 Hit";
    const key = `${row.Type}|||${status}`;
    statusCounts.set(key, (statusCounts.get(key) || 0) + 1);
    if (!statusDrugs.has(key)) statusDrugs.set(key, []);
    statusDrugs.get(key).push(row.Drug);
  }
  const preferredTypeOrder = ["Sulfonyl Fluorides", "Fluorosulfates", "Sulfonyl Triazoles"];
  const types = [...new Set(compoundRecords.map((row) => row.Type))].sort((a, b) => {
    const ai = preferredTypeOrder.indexOf(a);
    const bi = preferredTypeOrder.indexOf(b);
    if (ai !== -1 || bi !== -1) return (ai === -1 ? 99 : ai) - (bi === -1 ? 99 : bi);
    return a.localeCompare(b);
  });
  safePlot(
    "compoundHitsBinary",
    ["0 Hits", "\u22651 Hit"].map((status) => ({
      type: "bar",
      name: status,
      x: types,
      y: types.map((type) => statusCounts.get(`${type}|||${status}`) || 0),
      customdata: types.map((type) => ({
        title: `${status} - ${type}`,
        drugs: (statusDrugs.get(`${type}|||${status}`) || []).sort(),
      })),
      marker: { color: status === "0 Hits" ? MORANDI.gray : MORANDI.red },
      hoverinfo: "none",
    })),
    mergeAxes(baseLayout({ barmode: "stack", margin: { l: 52, r: 2, t: 4, b: 58 } }), { title: "", tickangle: 28, automargin: true }, { title: "Number of Compounds" })
  );
  bindCompoundSummaryClick("compoundHitsBinary");

  const bins = ["1 Hit", "2-4 Hits", "5-10 Hits", ">10 Hits"];
  const binColor = {
    "1 Hit": MORANDI.blueLight,
    "2-4 Hits": MORANDI.blueMid,
    "5-10 Hits": MORANDI.blue,
    ">10 Hits": MORANDI.blueDark,
  };
  const binCounts = new Map();
  const binDrugs = new Map();
  for (const row of compoundRecords.filter((r) => r.Hits > 0)) {
    const bin = row.Hits === 1 ? "1 Hit" : row.Hits <= 4 ? "2-4 Hits" : row.Hits <= 10 ? "5-10 Hits" : ">10 Hits";
    const key = `${row.Type}|||${bin}`;
    binCounts.set(key, (binCounts.get(key) || 0) + 1);
    if (!binDrugs.has(key)) binDrugs.set(key, []);
    binDrugs.get(key).push(row.Drug);
  }
  safePlot(
    "compoundHitsBinned",
    bins.map((bin) => ({
      type: "bar",
      name: bin,
      x: types,
      y: types.map((type) => binCounts.get(`${type}|||${bin}`) || 0),
      customdata: types.map((type) => ({
        title: `${bin} - ${type}`,
        drugs: (binDrugs.get(`${type}|||${bin}`) || []).sort(),
      })),
      marker: { color: binColor[bin] },
      hoverinfo: "none",
    })),
    mergeAxes(baseLayout({ barmode: "group", margin: { l: 52, r: 2, t: 4, b: 58 } }), { title: "", tickangle: 28, automargin: true }, { title: "Number of Compounds" })
  );
  bindCompoundSummaryClick("compoundHitsBinned");

}

async function loadCircularCompoundAtlas() {
  const rel = state.manifest?.atlas?.circularCompoundDomain;
  if (!rel) return null;
  if (!state.circularAtlas) state.circularAtlas = await fetchJson(`atlas/${rel}`);
  return state.circularAtlas;
}

function circularColor(value) {
  const ratio = Math.max(0, Number(value || 0) / 1000);
  if (ratio <= 0) return "rgba(255,255,255,0)";
  if (ratio < 2.5) return "#d9edf4";
  if (ratio < 5) return "#8dc6df";
  if (ratio < 10) return "#3e88bf";
  return "#1f3136";
}

function circularDataset() {
  return state.circularAtlas?.datasets?.[state.datasetKey] || null;
}

function circularCellKey(rowIndex, compoundIndex) {
  return `${rowIndex}:${compoundIndex}`;
}

function renderCircularLegend() {
  const legend = el("circularAtlasLegend");
  if (!legend) return;
  legend.innerHTML = `
    <span style="background:${circularColor(2200)}"></span><small>R > 2</small>
    <span style="background:${circularColor(5000)}"></span><small>R > 5</small>
    <span style="background:${circularColor(10000)}"></span><small>R > 10</small>
    <span style="background:#fff"></span><small>No hit</small>`;
}

function renderCircularLayerLabels() {
  const canvas = el("circularCompoundAtlas");
  const stage = canvas?.parentElement;
  if (!stage) return;
  let labels = stage.querySelector(".circular-layer-labels");
  if (!labels) {
    labels = document.createElement("div");
    labels.className = "circular-layer-labels";
    stage.appendChild(labels);
  }
  const canvasRect = canvas.getBoundingClientRect();
  const stageRect = stage.getBoundingClientRect();
  const size = canvasRect.width || 760;
  labels.style.left = `${canvasRect.left - stageRect.left}px`;
  labels.style.top = `${canvasRect.top - stageRect.top}px`;
  labels.style.width = `${size}px`;
  labels.style.height = `${size}px`;
  const centerX = size / 2;
  const centerY = size / 2;
  const zoom = state.circularZoom;
  labels.innerHTML = "";
  CIRCULAR_LAYER_LABELS.forEach(({ level, label }) => {
    const band = CIRCULAR_BANDS[level];
    const radius = size * ((band[0] + band[1]) / 2);
    const node = document.createElement("span");
    node.textContent = label;
    node.style.left = `${centerX + zoom.offsetX + Math.cos(CIRCULAR_LAYER_LABEL_ANGLE) * radius * zoom.scale}px`;
    node.style.top = `${centerY + zoom.offsetY + Math.sin(CIRCULAR_LAYER_LABEL_ANGLE) * radius * zoom.scale}px`;
    labels.appendChild(node);
  });
}

function drawCircularCell(ctx, cx, cy, startAngle, endAngle, innerRadius, outerRadius, color) {
  ctx.beginPath();
  ctx.arc(cx, cy, outerRadius, startAngle, endAngle);
  ctx.arc(cx, cy, innerRadius, endAngle, startAngle, true);
  ctx.closePath();
  ctx.fillStyle = color;
  ctx.fill();
}

function maxByRow(cells = []) {
  const byRow = new Map();
  for (const [rowIndex, , maxR, hitCount] of cells) {
    const current = byRow.get(rowIndex) || { maxR: 0, hitCount: 0 };
    current.maxR = Math.max(current.maxR, Number(maxR || 0));
    current.hitCount += Number(hitCount || 0);
    byRow.set(rowIndex, current);
  }
  return byRow;
}

function circularWarheadClasses(data) {
  const present = new Set((data?.compounds || []).map((compound) => compound.type || "Unknown"));
  const ordered = [...CIRCULAR_WARHEAD_CLASS_ORDER];
  const extras = [...present].filter((className) => !ordered.includes(className)).sort((a, b) => a.localeCompare(b));
  return [...ordered, ...extras];
}

function circularWarheadLabel(className) {
  return className === "Unknown" ? "Unknown" : compoundTypeLabel(className);
}

function maxByRowAndWarhead(data, cells = []) {
  const classes = circularWarheadClasses(data);
  const classIndex = new Map(classes.map((name, index) => [name, index]));
  const byRow = new Map();
  for (const [rowIndex, compoundIndex, maxR, hitCount] of cells) {
    const compound = data?.compounds?.[compoundIndex] || {};
    const className = compound.type || "Unknown";
    const index = classIndex.get(className) ?? 0;
    if (!byRow.has(rowIndex)) byRow.set(rowIndex, classes.map((name) => ({ className: name, maxR: 0, hitCount: 0 })));
    const current = byRow.get(rowIndex)[index];
    current.maxR = Math.max(current.maxR, Number(maxR || 0));
    current.hitCount += Number(hitCount || 0);
  }
  return byRow;
}

function hierarchyChildren(data, parentId) {
  const drilldown = data?.drilldowns?.[parentId];
  return drilldown ? { level: drilldown.level || "", rows: drilldown.rows || [], cells: drilldown.cells || [] } : { level: "", rows: [], cells: [] };
}

function addCircularAggregateCell(target, compoundIndex, maxR, hitCount) {
  const cIndex = Number(compoundIndex);
  if (!Number.isFinite(cIndex)) return;
  if (Number(hitCount || 0) <= 0) return;
  const current = target.get(cIndex) || { maxR: 0, hitCount: 0 };
  current.maxR = Math.max(current.maxR, Number(maxR || 0));
  current.hitCount += Number(hitCount || 0);
  target.set(cIndex, current);
}

function mergeCircularAggregateCells(target, source) {
  for (const [compoundIndex, cell] of source.entries()) {
    addCircularAggregateCell(target, compoundIndex, cell.maxR, cell.hitCount);
  }
  return target;
}

function circularCellsForRow(cells = [], rowIndex = -1) {
  const aggregate = new Map();
  for (const [cellRowIndex, compoundIndex, maxR, hitCount] of cells || []) {
    if (Number(cellRowIndex) !== Number(rowIndex)) continue;
    addCircularAggregateCell(aggregate, compoundIndex, maxR, hitCount);
  }
  return aggregate;
}

function circularInheritedSingleSiteCells(data, row, parentCells = null, parentRowIndex = -1) {
  const drilldown = hierarchyChildren(data, row?.id);
  if (!drilldown.rows?.length || drilldown.cells?.length) return null;
  const directAggregate = circularCellsForRow(parentCells || [], parentRowIndex);
  if (!circularAggregateHasHits(directAggregate)) return null;
  const descendantSites = descendantCircularSiteRows(data, row);
  if (descendantSites.length === 1) {
    const siteRow = Number(descendantSites[0]?.siteRow ?? descendantSites[0]?.rowId);
    const childIndex = drilldown.rows.findIndex((childRow) => (
      childRow.level === "site" &&
      Number(childRow.siteRow ?? childRow.rowId) === siteRow
    ));
    if (childIndex >= 0) return circularCellsFromAggregate(directAggregate, childIndex);
  }
  return null;
}

function circularDrilldownCellsForRow(data, row, parentCells = null, parentRowIndex = -1) {
  const drilldown = hierarchyChildren(data, row?.id);
  if (drilldown.cells?.length) return drilldown.cells;
  return circularInheritedSingleSiteCells(data, row, parentCells, parentRowIndex) || drilldown.cells || [];
}

function circularAggregateCellsForRow(data, row, parentCells = null, parentRowIndex = -1, memo = null) {
  if (!row) return new Map();
  if (row.level === "site") return circularCellsForRow(parentCells || [], parentRowIndex);
  const cache = memo || new Map();
  if (cache.has(row.id)) return new Map(cache.get(row.id));
  const drilldown = hierarchyChildren(data, row.id);
  const drilldownCells = circularDrilldownCellsForRow(data, row, parentCells, parentRowIndex);
  const aggregate = new Map();
  (drilldown.rows || []).forEach((childRow, childIndex) => {
    const childAggregate =
      childRow.level === "site"
        ? circularCellsForRow(drilldownCells, childIndex)
        : circularAggregateCellsForRow(data, childRow, drilldownCells, childIndex, cache);
    mergeCircularAggregateCells(aggregate, childAggregate);
  });
  cache.set(row.id, new Map(aggregate));
  return aggregate;
}

function circularAggregateSummary(aggregate) {
  const summary = { maxR: 0, hitCount: 0 };
  for (const cell of aggregate.values()) {
    summary.maxR = Math.max(summary.maxR, Number(cell.maxR || 0));
    summary.hitCount += Number(cell.hitCount || 0);
  }
  return summary;
}

function circularAggregateHasHits(aggregate) {
  for (const cell of aggregate.values()) {
    if (Number(cell.hitCount || 0) > 0) return true;
  }
  return false;
}

function circularClassStatsFromAggregate(data, aggregate) {
  const classes = circularWarheadClasses(data);
  const classIndex = new Map(classes.map((name, index) => [name, index]));
  const stats = classes.map((name) => ({ className: circularWarheadLabel(name), maxR: 0, hitCount: 0, siteCount: 0 }));
  for (const [compoundIndex, cell] of aggregate.entries()) {
    const compound = data?.compounds?.[compoundIndex] || {};
    const className = compound.type || "Unknown";
    const index = classIndex.get(className) ?? 0;
    const current = stats[index];
    current.maxR = Math.max(current.maxR, Number(cell.maxR || 0));
    current.hitCount += Number(cell.hitCount || 0);
  }
  return stats;
}

function circularClassSiteCountsForRow(data, row, parentCells = null, parentRowIndex = -1) {
  const classes = circularWarheadClasses(data);
  const classIndex = new Map(classes.map((name, index) => [name, index]));
  const sitesByClass = classes.map(() => new Set());
  const addSiteCells = (siteRow, cells = [], rowIndex = -1) => {
    const siteKey = String(siteRow?.siteRow ?? siteRow?.rowId ?? siteRow?.id ?? rowIndex);
    for (const [cellRowIndex, compoundIndex, , hitCount] of cells || []) {
      if (Number(cellRowIndex) !== Number(rowIndex) || Number(hitCount || 0) <= 0) continue;
      const compound = data?.compounds?.[compoundIndex] || {};
      const index = classIndex.get(compound.type || "Unknown") ?? 0;
      sitesByClass[index].add(siteKey);
    }
  };
  const visit = (currentRow, cells, rowIndex) => {
    if (!currentRow) return;
    if (currentRow.level === "site") {
      addSiteCells(currentRow, cells, rowIndex);
      return;
    }
    const drilldown = hierarchyChildren(data, currentRow.id);
    const drilldownCells = circularDrilldownCellsForRow(data, currentRow, cells, rowIndex);
    (drilldown.rows || []).forEach((childRow, childIndex) => visit(childRow, drilldownCells, childIndex));
  };
  visit(row, parentCells || [], parentRowIndex);
  return sitesByClass.map((sites) => sites.size);
}

function circularCellsFromAggregate(aggregate, rowIndex) {
  return [...aggregate.entries()]
    .filter(([, cell]) => Number(cell.hitCount || 0) > 0)
    .map(([compoundIndex, cell]) => [rowIndex, compoundIndex, cell.maxR, cell.hitCount])
    .sort((a, b) => a[1] - b[1]);
}

function circularHitSiteSpanCount(data, row, parentCells = null, parentRowIndex = -1, memo = null) {
  if (!row) return 0;
  const aggregate = circularAggregateCellsForRow(data, row, parentCells, parentRowIndex, memo);
  if (!circularAggregateHasHits(aggregate)) return 0;
  if (row.level === "site") return 1;
  const drilldown = hierarchyChildren(data, row.id);
  if (!drilldown.level || !drilldown.rows.length) return 0;
  const drilldownCells = circularDrilldownCellsForRow(data, row, parentCells, parentRowIndex);
  return drilldown.rows.reduce(
    (sum, childRow, childIndex) => sum + circularHitSiteSpanCount(data, childRow, drilldownCells, childIndex, memo),
    0,
  );
}

const CIRCULAR_ROOT_ORDER_KEY = "__root__";

function circularOrderKey(parentId) {
  return parentId || CIRCULAR_ROOT_ORDER_KEY;
}

function orderedCircularChildren(parentId, rows) {
  return orderedCircularChildEntries(parentId, rows).map((item) => item.row);
}

function orderedCircularChildEntries(parentId, rows) {
  const order = state.circularConstrainedOrderMap.get(circularOrderKey(parentId));
  const entries = (rows || []).map((row, index) => ({ row, index }));
  if (!order?.length) return entries;
  const byId = new Map(order.map((id, index) => [id, index]));
  return entries
    .sort((a, b) => {
      const aOrder = byId.has(a.row.id) ? byId.get(a.row.id) : Number.MAX_SAFE_INTEGER;
      const bOrder = byId.has(b.row.id) ? byId.get(b.row.id) : Number.MAX_SAFE_INTEGER;
      return aOrder - bOrder || a.index - b.index;
    });
}

function descendantCircularSiteRows(data, row) {
  if (!row) return [];
  if (row.level === "site" && Number.isFinite(Number(row.siteRow ?? row.rowId))) return [row];
  const drilldown = hierarchyChildren(data, row.id);
  return (drilldown.rows || []).flatMap((child) => descendantCircularSiteRows(data, child));
}

function circularPathKey(pathIds, row, rowIndex) {
  return [...pathIds, `${row?.id || "row"}#${rowIndex}`].join(">");
}

function circularEntrySpecificity(row) {
  if (!row) return 0;
  if (row.level === "domain" || row.level === "repeat") return 3;
  if (row.level === "family") return 2;
  if (row.level === "superfamily") return 1;
  return 0;
}

function circularBetterSiteCandidate(candidate, current) {
  if (!current) return true;
  if (candidate.specificity !== current.specificity) return candidate.specificity > current.specificity;
  const candidateInformative = candidate.parentHitSites > 1;
  const currentInformative = current.parentHitSites > 1;
  if (candidateInformative !== currentInformative) return candidateInformative;
  if (candidate.parentHitSites !== current.parentHitSites) return candidate.parentHitSites < current.parentHitSites;
  if (candidate.pathDepth !== current.pathDepth) return candidate.pathDepth > current.pathDepth;
  return candidate.pathKey.localeCompare(current.pathKey) < 0;
}

function buildCircularCanonicalSitePathMap(data) {
  const candidates = new Map();
  const aggregateMemo = new Map();
  const parentSizeMemo = new Map();
  const parentHitSiteCount = (row, parentCells, parentRowIndex, pathKey) => {
    if (!row) return 0;
    if (parentSizeMemo.has(pathKey)) return parentSizeMemo.get(pathKey);
    const count = circularHitSiteSpanCount(data, row, parentCells, parentRowIndex, aggregateMemo);
    parentSizeMemo.set(pathKey, count);
    return count;
  };
  const visit = (row, parentCells = [], rowIndex = -1, path = [], pathIds = []) => {
    const aggregate = circularAggregateCellsForRow(data, row, parentCells, rowIndex, aggregateMemo);
    if (!circularAggregateHasHits(aggregate)) return;
    const pathKey = circularPathKey(pathIds, row, rowIndex);
    if (row.level === "site") {
      const siteRow = Number(row.siteRow ?? row.rowId);
      if (!Number.isFinite(siteRow)) return;
      const terminal = path[path.length - 1] || null;
      const parentHitSites = terminal
        ? parentHitSiteCount(terminal.row, terminal.parentCells, terminal.rowIndex, terminal.pathKey)
        : 1;
      const candidate = {
        pathKey,
        specificity: circularEntrySpecificity(terminal?.row),
        parentHitSites,
        pathDepth: path.length,
      };
      if (circularBetterSiteCandidate(candidate, candidates.get(siteRow))) candidates.set(siteRow, candidate);
      return;
    }
    const drilldown = hierarchyChildren(data, row.id);
    const drilldownCells = circularDrilldownCellsForRow(data, row, parentCells, rowIndex);
    const nextPath = [...path, { row, parentCells, rowIndex, pathKey }];
    const nextPathIds = [...pathIds, `${row.id}#${rowIndex}`];
    (drilldown.rows || []).forEach((childRow, childIndex) => visit(childRow, drilldownCells, childIndex, nextPath, nextPathIds));
  };
  (data.rows || []).forEach((row, index) => visit(row, data.cells || [], index, [], []));
  return new Map([...candidates.entries()].map(([siteRow, candidate]) => [siteRow, candidate.pathKey]));
}

async function circularRepresentativeFeature(data, row) {
  const siteRows = descendantCircularSiteRows(data, row);
  const siteFeatures = await Promise.all(siteRows.map(circularSiteFeature));
  const valuesByCompound = new Map();
  for (const feature of siteFeatures) {
    for (const [compound, value] of feature.values.entries()) {
      if (!Number.isFinite(value) || value <= 0) continue;
      const current = valuesByCompound.get(compound) || { sum: 0, count: 0 };
      current.sum += Math.log2(value);
      current.count += 1;
      valuesByCompound.set(compound, current);
    }
  }
  const values = new Map();
  for (const [compound, current] of valuesByCompound.entries()) {
    if (current.count > 0) values.set(compound, 2 ** (current.sum / current.count));
  }
  return { row, values };
}

async function clusteredCircularChildOrder(data, rows) {
  if (!rows?.length || rows.length < 3) return rows.map((row) => row.id);
  const features = await Promise.all(rows.map((row) => circularRepresentativeFeature(data, row)));
  const usable = features
    .map((feature, index) => ({ feature, index }))
    .filter((item) => item.feature.values.size >= MIN_CLUSTER_SHARED_COMPOUNDS);
  if (usable.length < 3) return rows.map((row) => row.id);
  const usableFeatures = usable.map((item) => item.feature);
  const root = buildMissingAwareAgglomerativeDendrogram(usableFeatures, { allowBridges: true }).root;
  if (!root) return rows.map((row) => row.id);
  const clustered = rotateOrderAtLargestGap(optimizeDendrogramLeafOrder(root, usableFeatures), usableFeatures)
    .map((localIndex) => usable[localIndex].index);
  const included = new Set(clustered);
  const remainder = rows.map((_, index) => index).filter((index) => !included.has(index));
  return [...clustered, ...remainder].map((index) => rows[index].id);
}

async function buildCircularConstrainedOrderMap(data) {
  const orderMap = new Map();
  const aggregateMemo = new Map();
  const visitGroup = async (parentId, rows, cells = null) => {
    const visibleRows = (rows || []).filter((row, index) => circularHitSiteSpanCount(data, row, cells || [], index, aggregateMemo) > 0);
    if (visibleRows.length) {
      orderMap.set(circularOrderKey(parentId), await clusteredCircularChildOrder(data, visibleRows));
    }
    await mapWithConcurrency(visibleRows, 8, async (row) => {
      const drilldown = hierarchyChildren(data, row.id);
      if (drilldown.rows?.length) await visitGroup(row.id, drilldown.rows, circularDrilldownCellsForRow(data, row, cells || [], rows.indexOf(row)));
    });
  };
  await visitGroup("", data.rows || [], data.cells || []);
  return orderMap;
}

function circularConstrainedOrderKey(data) {
  const rows = data?.rows || [];
  return `${state.datasetKey}:${rows.map((row) => row.id).join("|")}`;
}

function prepareCircularConstrainedOrders(data) {
  const key = circularConstrainedOrderKey(data);
  if (state.circularConstrainedOrderKey === key && state.circularConstrainedOrderMap.size) return state.circularConstrainedOrderStatus;
  if (state.circularConstrainedOrderKey === key && state.circularConstrainedOrderPromise) return "ordering sectors by parent-constrained ligandability clustering...";
  state.circularConstrainedOrderKey = key;
  state.circularConstrainedOrderMap = new Map();
  state.circularConstrainedOrderStatus = "ordering sectors by parent-constrained ligandability clustering...";
  state.circularConstrainedOrderPromise = buildCircularConstrainedOrderMap(data)
    .then((orderMap) => {
      if (state.circularConstrainedOrderKey === key) {
        state.circularConstrainedOrderMap = orderMap;
        state.circularConstrainedOrderStatus = "parent-constrained ligandability sector ordering active";
        state.circularConstrainedOrderPromise = null;
        renderCircularCompoundAtlas();
      }
      return orderMap;
    })
    .catch((err) => {
      console.warn(`Unable to order circular sectors by ligandability: ${err.message}`);
      if (state.circularConstrainedOrderKey === key) {
        state.circularConstrainedOrderStatus = "parent-constrained ligandability sector ordering unavailable";
        state.circularConstrainedOrderPromise = null;
      }
    });
  return state.circularConstrainedOrderStatus;
}

function collectCircularHierarchy(data) {
  const segments = [];
  const topRows = data.rows || [];
  const full = Math.PI * 2;
  const classTemplate = circularClassStatsFromAggregate(data, new Map());
  const canonicalSitePathMap = buildCircularCanonicalSitePathMap(data);
  const addSegment = (level, row, start, end, stats, band, classCells, parentId = "", parentCells = null, parentRowIndex = -1) => {
    if (!row || end <= start) return;
    segments.push({
      level,
      row,
      id: row.id,
      parentId,
      parentCells,
      parentRowIndex,
      label: row.label || row.id,
      start,
      end,
      inner: band[0],
      outer: band[1],
      maxR: stats?.maxR || 0,
      hitCount: stats?.hitCount || 0,
      classCells: classCells || classTemplate,
    });
  };
  const bandForLevel = (level) => {
    if (level === "superfamily") return CIRCULAR_BANDS.superfamily;
    if (level === "family") return CIRCULAR_BANDS.family;
    if (level === "domain" || level === "repeat") return CIRCULAR_BANDS.domain;
    if (level === "site" || level === "site engagement") return CIRCULAR_BANDS.site;
    return CIRCULAR_BANDS.domain;
  };
  const displayLevel = (level) => (level === "site" ? "site engagement" : level);
  const aggregateMemo = new Map();
  const localAggregate = (row, parentCells = null, parentRowIndex = -1, pathIds = []) => {
    if (!row) return new Map();
    const pathKey = circularPathKey(pathIds, row, parentRowIndex);
    if (row.level === "site") {
      const siteRow = Number(row.siteRow ?? row.rowId);
      return canonicalSitePathMap.get(siteRow) === pathKey ? circularCellsForRow(parentCells || [], parentRowIndex) : new Map();
    }
    const drilldown = hierarchyChildren(data, row.id);
    const drilldownCells = circularDrilldownCellsForRow(data, row, parentCells, parentRowIndex);
    const nextPathIds = [...pathIds, `${row.id}#${parentRowIndex}`];
    const aggregate = new Map();
    (drilldown.rows || []).forEach((childRow, childIndex) => {
      mergeCircularAggregateCells(aggregate, localAggregate(childRow, drilldownCells, childIndex, nextPathIds));
    });
    return aggregate;
  };
  const localClassSiteCounts = (row, parentCells = null, parentRowIndex = -1, pathIds = []) => {
    const classes = circularWarheadClasses(data);
    const classIndex = new Map(classes.map((name, index) => [name, index]));
    const sitesByClass = classes.map(() => new Set());
    const visit = (currentRow, cells, rowIndex, currentPathIds) => {
      if (!currentRow) return;
      const pathKey = circularPathKey(currentPathIds, currentRow, rowIndex);
      if (currentRow.level === "site") {
        const siteRow = Number(currentRow.siteRow ?? currentRow.rowId);
        if (canonicalSitePathMap.get(siteRow) !== pathKey) return;
        for (const [cellRowIndex, compoundIndex, , hitCount] of cells || []) {
          if (Number(cellRowIndex) !== Number(rowIndex) || Number(hitCount || 0) <= 0) continue;
          const compound = data?.compounds?.[compoundIndex] || {};
          const index = classIndex.get(compound.type || "Unknown") ?? 0;
          sitesByClass[index].add(String(siteRow));
        }
        return;
      }
      const drilldown = hierarchyChildren(data, currentRow.id);
      const drilldownCells = circularDrilldownCellsForRow(data, currentRow, cells, rowIndex);
      const nextPathIds = [...currentPathIds, `${currentRow.id}#${rowIndex}`];
      (drilldown.rows || []).forEach((childRow, childIndex) => visit(childRow, drilldownCells, childIndex, nextPathIds));
    };
    visit(row, parentCells || [], parentRowIndex, pathIds);
    return sitesByClass.map((sites) => sites.size);
  };
  const siteSpanCount = (row, parentCells = null, parentRowIndex = -1, pathIds = []) => {
    const aggregate = localAggregate(row, parentCells, parentRowIndex, pathIds);
    if (!circularAggregateHasHits(aggregate)) return 0;
    if (row.level === "site") return 1;
    const drilldown = hierarchyChildren(data, row.id);
    const drilldownCells = circularDrilldownCellsForRow(data, row, parentCells, parentRowIndex);
    const nextPathIds = [...pathIds, `${row.id}#${parentRowIndex}`];
    return (drilldown.rows || []).reduce(
      (sum, childRow, childIndex) => sum + siteSpanCount(childRow, drilldownCells, childIndex, nextPathIds),
      0,
    );
  };
  const retainedTopRows = orderedCircularChildEntries("", topRows)
    .map(({ row, index }) => ({ row, index, count: siteSpanCount(row, data.cells || [], index, []) }))
    .filter((item) => item.count > 0);
  const totalSites = retainedTopRows.reduce((sum, item) => sum + item.count, 0);
  const anglePerSite = full / Math.max(1, totalSites);

  const addHierarchyNode = (row, start, end, parentId = "", parentCells = null, parentRowIndex = -1, pathIds = []) => {
    const band = bandForLevel(row.level);
    const aggregate = localAggregate(row, parentCells, parentRowIndex, pathIds);
    if (!circularAggregateHasHits(aggregate)) return;
    const classCells = circularClassStatsFromAggregate(data, aggregate);
    const siteCounts = localClassSiteCounts(row, parentCells, parentRowIndex, pathIds);
    classCells.forEach((cell, index) => {
      cell.siteCount = siteCounts[index] || 0;
    });
    addSegment(displayLevel(row.level), row, start, end, circularAggregateSummary(aggregate), band, classCells, parentId, parentCells, parentRowIndex);
    const drilldown = hierarchyChildren(data, row.id);
    if (!drilldown.level || !drilldown.rows.length) return;
    const drilldownCells = circularDrilldownCellsForRow(data, row, parentCells, parentRowIndex);
    const nextPathIds = [...pathIds, `${row.id}#${parentRowIndex}`];
    let localStart = start;
    orderedCircularChildEntries(row.id, drilldown.rows).forEach(({ row: childRow, index: childIndex }) => {
      const count = siteSpanCount(childRow, drilldownCells, childIndex, nextPathIds);
      if (!count) return;
      const childStart = localStart;
      const childEnd = childStart + count * anglePerSite;
      localStart = childEnd;
      addHierarchyNode(childRow, childStart, childEnd, row.id, drilldownCells, childIndex, nextPathIds);
    });
  };

  let cursor = -Math.PI / 2;
  retainedTopRows.forEach(({ row, index, count }) => {
    const start = cursor;
    const end = start + count * anglePerSite;
    cursor = end;
    addHierarchyNode(row, start, end, "", data.cells || [], index, []);
  });
  return { segments };
}

function descendantCircularSiteCount(data, row) {
  let count = 0;
  const stack = [hierarchyChildren(data, row.id)];
  while (stack.length) {
    const drilldown = stack.pop();
    if (!drilldown?.level) continue;
    if (drilldown.level === "site") {
      count += drilldown.rows.length;
      continue;
    }
    for (const child of drilldown.rows) stack.push(hierarchyChildren(data, child.id));
  }
  return count;
}

function applyCircularTransform(ctx, size) {
  const zoom = state.circularZoom;
  ctx.translate(size / 2 + zoom.offsetX, size / 2 + zoom.offsetY);
  ctx.scale(zoom.scale, zoom.scale);
  ctx.translate(-size / 2, -size / 2);
}

function circularCanvasPoint(event) {
  const canvas = el("circularCompoundAtlas");
  const rect = canvas.getBoundingClientRect();
  const size = rect.width || 760;
  const zoom = state.circularZoom;
  return {
    size,
    x: (event.clientX - rect.left - zoom.offsetX - size / 2) / zoom.scale + size / 2,
    y: (event.clientY - rect.top - zoom.offsetY - size / 2) / zoom.scale + size / 2,
  };
}

function renderCircularCompoundAtlas() {
  const canvas = el("circularCompoundAtlas");
  const status = el("circularAtlasStatus");
  if (!canvas || !status) return;
  const data = circularDataset();
  if (!data) {
    status.textContent = "Compound-domain matrix unavailable.";
    return;
  }
  const rows = data.rows || [];
  const compounds = data.compounds || [];
  const rect = canvas.getBoundingClientRect();
  const size = Math.max(420, Math.floor(rect.width || 760));
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.round(size * dpr);
  canvas.height = Math.round(size * dpr);
  canvas.style.height = `${size}px`;
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, size, size);
  state.circularSectors = [];
  state.circularCellMap = new Map();
  if (!rows.length || !compounds.length) {
    status.textContent = "No compound-domain hit matrix for the current dataset.";
    return;
  }
  const cx = size / 2;
  const cy = size / 2;
  ctx.save();
  applyCircularTransform(ctx, size);
  const constrainedOrderStatus = prepareCircularConstrainedOrders(data);
  const { segments } = collectCircularHierarchy(data);

  for (const segment of segments) {
    const innerRadius = size * segment.inner;
    const outerRadius = size * segment.outer;
    const angularGap = segment.level === "site engagement" ? 0.0005 : 0.0012;
    const start = segment.start + angularGap;
    const end = segment.end - angularGap;
    drawCircularCell(ctx, cx, cy, start, end, innerRadius, outerRadius, "#f8fbfc");
    const classCells = segment.classCells || [];
    const classWidth = (outerRadius - innerRadius) / Math.max(1, classCells.length);
    classCells.forEach((cell, index) => {
      drawCircularCell(
        ctx,
        cx,
        cy,
        start,
        end,
        innerRadius + index * classWidth + 0.2,
        innerRadius + (index + 1) * classWidth - 0.2,
        circularColor(cell.maxR),
      );
    });
    state.circularSectors.push(segment);
  }
  const ringGaps = [
    CIRCULAR_BANDS.site[0],
    CIRCULAR_BANDS.site[1],
    CIRCULAR_BANDS.domain[0],
    CIRCULAR_BANDS.domain[1],
    CIRCULAR_BANDS.family[0],
    CIRCULAR_BANDS.family[1],
    CIRCULAR_BANDS.superfamily[0],
  ];
  ringGaps.forEach((radius) => {
    ctx.strokeStyle = "#fff";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(cx, cy, size * radius, 0, Math.PI * 2);
    ctx.stroke();
  });
  ctx.beginPath();
  ctx.arc(cx, cy, size * Math.max(0, CIRCULAR_BANDS.site[0] - 0.005), 0, Math.PI * 2);
  ctx.fillStyle = "#fff";
  ctx.fill();
  const visibleSiteSegments = segments.filter((segment) => segment.level === "site engagement" && Number.isFinite(Number(segment.row?.siteRow ?? segment.row?.rowId)));
  const centerDendrogram = prepareCircularVisibleSiteDendrogram(visibleSiteSegments);
  if (centerDendrogram?.root) drawCircularCenterDendrogram(ctx, centerDendrogram, size);
  ctx.restore();
  const siteSegments = segments.filter((segment) => segment.level === "site engagement").length;
  const familySegments = segments.filter((segment) => segment.level === "family").length;
  const domainSegments = segments.filter((segment) => segment.level === "domain" || segment.level === "repeat").length;
  const dendrogramStatus =
    centerDendrogram?.status ||
    (state.circularVisibleSiteDendrogramPromise ? "center dendrogram: clustering visible sites..." : "");
  status.textContent = `${familySegments.toLocaleString()} family, ${domainSegments.toLocaleString()} domain/repeat, ${siteSegments.toLocaleString()} site sectors. ${compounds.length.toLocaleString()} hit compounds. Wheel to zoom, drag to pan, double-click to reset.${constrainedOrderStatus ? ` ${constrainedOrderStatus}.` : ""}${dendrogramStatus ? ` ${dendrogramStatus}.` : ""}`;
  renderCircularLegend();
  renderCircularLayerLabels();
}

function circularHitAtEvent(event) {
  const data = circularDataset();
  if (!el("circularCompoundAtlas") || !data) return null;
  const point = circularCanvasPoint(event);
  const size = point.size;
  const x = point.x - size / 2;
  const y = point.y - size / 2;
  const radius = Math.sqrt(x * x + y * y);
  let angle = Math.atan2(y, x) + Math.PI / 2;
  if (angle < 0) angle += Math.PI * 2;
  const worldAngle = angle - Math.PI / 2;
  const segment = [...state.circularSectors].reverse().find((item) => {
    const inner = size * item.inner;
    const outer = size * item.outer;
    return radius >= inner && radius <= outer && worldAngle >= item.start && worldAngle <= item.end;
  });
  if (!segment) return null;
  const classCells = segment.classCells || [];
  const classIndex = classCells.length
    ? Math.max(0, Math.min(classCells.length - 1, Math.floor(((radius / size - segment.inner) / Math.max(0.0001, segment.outer - segment.inner)) * classCells.length)))
    : -1;
  const classCell = classIndex >= 0 ? classCells[classIndex] : null;
  return { row: segment.row, segment, classIndex, cell: classCell || { maxR: segment.maxR, hitCount: segment.hitCount } };
}

function circularLayerLabel(level) {
  if (level === "site engagement" || level === "site") return "Site";
  if (level === "domain" || level === "repeat") return "Domain";
  if (level === "family") return "Family";
  if (level === "superfamily") return "Superfamily";
  return String(level || "biology").replace(/\b\w/g, (letter) => letter.toUpperCase());
}

function showCircularTooltip(event, hit) {
  const tooltip = el("circularAtlasTooltip");
  if (!tooltip) return;
  if (!hit?.row) {
    tooltip.hidden = true;
    return;
  }
  const className = hit.cell?.className ? `${escapeHtml(hit.cell.className)} · ` : "";
  tooltip.innerHTML = `<strong>${escapeHtml(hit.row.label || hit.row.id)}</strong><br>${escapeHtml(circularLayerLabel(hit.segment?.level))}<br>${
    hit.cell ? `${className}Max R ${(hit.cell.maxR / 1000).toFixed(2)} · ${Number(hit.cell.siteCount || 0).toLocaleString()} hit site${hit.cell.siteCount === 1 ? "" : "s"}` : "No hit"
  }`;
  const stage = el("circularCompoundAtlas").parentElement.getBoundingClientRect();
  tooltip.style.left = `${event.clientX - stage.left + 12}px`;
  tooltip.style.top = `${event.clientY - stage.top + 12}px`;
  tooltip.hidden = false;
}

function rowOrderFromDendrogram(node) {
  if (!node) return [];
  if (node.leaf != null) return [node.leaf];
  return [...rowOrderFromDendrogram(node.left), ...rowOrderFromDendrogram(node.right)];
}

function sharedCompoundDistance(a, b) {
  return sharedCompoundDistanceWithMinimum(a, b, MIN_CLUSTER_SHARED_COMPOUNDS);
}

function sharedCompoundDistanceWithMinimum(a, b, minShared) {
  let sum = 0;
  let count = 0;
  for (const [compound, aValue] of a.values) {
    if (!b.values.has(compound)) continue;
    const bValue = b.values.get(compound);
    if (!Number.isFinite(aValue) || !Number.isFinite(bValue) || aValue <= 0 || bValue <= 0) continue;
    const delta = Math.log2(aValue) - Math.log2(bValue);
    sum += delta * delta;
    count += 1;
  }
  return count >= minShared ? { distance: Math.sqrt(sum / count), shared: count } : null;
}

function buildMissingAwareAgglomerativeDendrogram(features, options = {}) {
  const minValidPairFraction = options.minValidPairFraction ?? MIN_CLUSTER_VALID_PAIR_FRACTION;
  const allowBridges = options.allowBridges ?? false;
  const pairDistances = new Map();
  let validPairs = 0;
  const pairKey = (a, b) => `${Math.min(a, b)}:${Math.max(a, b)}`;
  for (let i = 0; i < features.length; i += 1) {
    for (let j = i + 1; j < features.length; j += 1) {
      const result = sharedCompoundDistance(features[i], features[j]);
      if (!result) continue;
      pairDistances.set(pairKey(i, j), result.distance);
      validPairs += 1;
    }
  }
  const totalPairs = (features.length * (features.length - 1)) / 2;
  if (!totalPairs || validPairs / totalPairs < minValidPairFraction) {
    return { root: null, validPairs, totalPairs };
  }
  let nextId = features.length;
  let clusters = features.map((_, index) => ({ id: index, members: [index], leaf: index, height: 0 }));
  let maxMergeHeight = 0;
  let bridgedMerges = 0;
  const clusterDistance = (a, b) => {
    let sum = 0;
    let count = 0;
    for (const ai of a.members) {
      for (const bi of b.members) {
        const distance = pairDistances.get(pairKey(ai, bi));
        if (!Number.isFinite(distance)) continue;
        sum += distance;
        count += 1;
      }
    }
    return count ? sum / count : Infinity;
  };
  const relaxedClusterDistance = (a, b) => {
    let sum = 0;
    let count = 0;
    for (const ai of a.members) {
      for (const bi of b.members) {
        const result = sharedCompoundDistanceWithMinimum(features[ai], features[bi], 1);
        if (!result) continue;
        sum += result.distance;
        count += 1;
      }
    }
    return count ? sum / count : Infinity;
  };
  while (clusters.length > 1) {
    let best = { i: -1, j: -1, distance: Infinity };
    for (let i = 0; i < clusters.length; i += 1) {
      for (let j = i + 1; j < clusters.length; j += 1) {
        const distance = clusterDistance(clusters[i], clusters[j]);
        if (distance < best.distance) best = { i, j, distance };
      }
    }
    let bridged = false;
    if (!Number.isFinite(best.distance)) {
      if (!allowBridges) break;
      for (let i = 0; i < clusters.length; i += 1) {
        for (let j = i + 1; j < clusters.length; j += 1) {
          const distance = relaxedClusterDistance(clusters[i], clusters[j]);
          if (distance < best.distance) best = { i, j, distance };
        }
      }
      if (!Number.isFinite(best.distance)) best = { i: 0, j: 1, distance: maxMergeHeight || 1 };
      bridged = true;
      bridgedMerges += 1;
    }
    const right = clusters.splice(best.j, 1)[0];
    const left = clusters.splice(best.i, 1)[0];
    const height = bridged ? Math.max(best.distance, maxMergeHeight * 1.04, maxMergeHeight + 0.02) : best.distance;
    maxMergeHeight = Math.max(maxMergeHeight, height);
    clusters.push({
      id: nextId,
      left,
      right,
      members: [...left.members, ...right.members],
      height,
      bridged,
    });
    nextId += 1;
  }
  return { root: clusters.length === 1 ? clusters[0] : null, validPairs, totalPairs, bridgedMerges };
}

function buildMissingAwareDendrogram(features) {
  return buildMissingAwareAgglomerativeDendrogram(features);
}

function completeCompoundCount(features) {
  if (!features.length) return 0;
  let shared = new Set(features[0].values.keys());
  for (const feature of features.slice(1)) {
    shared = new Set([...shared].filter((compound) => feature.values.has(compound)));
  }
  return shared.size;
}

async function circularSiteFeature(row) {
  const siteRow = Number(row?.siteRow ?? row?.rowId);
  const catalogRow = Number.isFinite(siteRow) ? state.rowById.get(siteRow) : null;
  const summary = catalogRow ? await siteSummaryForRow(catalogRow) : null;
  const values = new Map();
  for (const hit of summary?.hits || []) {
    const value = Number(hit[1]);
    if (hit[0] && Number.isFinite(value)) values.set(String(hit[0]), value);
  }
  return { row, values };
}

function reorderCircularHeatmap(rows, cells, order) {
  const rowMap = new Map(order.map((oldIndex, newIndex) => [oldIndex, newIndex]));
  return {
    rows: order.map((index) => rows[index]),
    cells: (cells || [])
      .filter(([rowIndex]) => rowMap.has(rowIndex))
      .map(([rowIndex, ...rest]) => [rowMap.get(rowIndex), ...rest]),
  };
}

async function clusterCircularHeatmapSites(rows, cells) {
  if (!rows?.length || rows.length < MIN_CLUSTER_SITE_ROWS) return { rows, cells, dendrogram: null, status: "" };
  if (!rows.every((row) => row.level === "site" && Number.isFinite(Number(row.siteRow ?? row.rowId)))) {
    return { rows, cells, dendrogram: null, status: "" };
  }
  const features = await Promise.all(rows.map(circularSiteFeature));
  if (features.some((feature) => feature.values.size < MIN_CLUSTER_SHARED_COMPOUNDS)) {
    return { rows, cells, dendrogram: null, status: `Not clustered: at least one site has fewer than ${MIN_CLUSTER_SHARED_COMPOUNDS} observed compounds` };
  }
  const completeCount = completeCompoundCount(features);
  const dendrogram = buildMissingAwareDendrogram(features);
  if (!dendrogram.root) {
    const pairSummary = dendrogram.totalPairs ? `${dendrogram.validPairs}/${dendrogram.totalPairs}` : "0/0";
    return { rows, cells, dendrogram: null, status: `Not clustered: only ${pairSummary} site pairs have >= ${MIN_CLUSTER_SHARED_COMPOUNDS} shared compounds` };
  }
  const order = rowOrderFromDendrogram(dendrogram.root);
  const reordered = reorderCircularHeatmap(rows, cells, order);
  const completeStatus =
    completeCount < MIN_COMPLETE_CLUSTER_COMPOUNDS
      ? `Complete-compound clustering skipped (${completeCount} compounds shared by all sites); missing-aware clustering used`
      : `Missing-aware clustering used (${completeCount} complete compounds also available)`;
  return {
    rows: reordered.rows,
    cells: reordered.cells,
    dendrogram: { root: dendrogram.root, order },
    status: `${completeStatus}; min shared pair overlap ${MIN_CLUSTER_SHARED_COMPOUNDS}`,
  };
}

function drawCircularDendrogram(ctx, dendrogram, orderedRows, x, y, width, cellSize) {
  if (!dendrogram) return;
  const root = dendrogram.root || dendrogram;
  const order = dendrogram.order || orderedRows.map((_, index) => index);
  const leafY = new Map(order.map((leafIndex, displayIndex) => [leafIndex, y + displayIndex * cellSize + cellSize / 2]));
  const maxHeight = Math.max(0.0001, root.height || 0.0001);
  const xForHeight = (height) => x + width - 4 - (height / maxHeight) * (width - 12);
  const drawNode = (node) => {
    if (node.leaf != null) return { x: x + width - 4, y: leafY.get(node.leaf) };
    const left = drawNode(node.left);
    const right = drawNode(node.right);
    const nodeX = xForHeight(node.height || 0);
    ctx.beginPath();
    ctx.moveTo(left.x, left.y);
    ctx.lineTo(nodeX, left.y);
    ctx.lineTo(nodeX, right.y);
    ctx.lineTo(right.x, right.y);
    ctx.stroke();
    return { x: nodeX, y: (left.y + right.y) / 2 };
  };
  ctx.save();
  ctx.strokeStyle = "rgba(47, 49, 54, 0.62)";
  ctx.lineWidth = 1;
  drawNode(root);
  ctx.restore();
}

function circularVisibleSiteDendrogramKey(siteSegments) {
  return `${state.datasetKey}:${siteSegments.map((segment) => segment.row?.siteRow ?? segment.row?.rowId ?? segment.id).join("|")}`;
}

function compoundVarianceForGroup(features, indexes, compound) {
  const values = [];
  for (const index of indexes) {
    const value = features[index].values.get(compound);
    if (Number.isFinite(value) && value > 0) values.push(Math.log2(value));
  }
  if (values.length < Math.max(3, Math.min(MIN_CLUSTER_SHARED_COMPOUNDS, Math.ceil(indexes.length * 0.2)))) return null;
  const mean = values.reduce((sum, value) => sum + value, 0) / values.length;
  const variance = values.reduce((sum, value) => sum + (value - mean) ** 2, 0) / values.length;
  return { variance, values };
}

function bestSplitCompound(features, indexes) {
  const compounds = new Set();
  indexes.forEach((index) => features[index].values.forEach((_, compound) => compounds.add(compound)));
  let best = null;
  compounds.forEach((compound) => {
    const stats = compoundVarianceForGroup(features, indexes, compound);
    if (!stats) return;
    if (!best || stats.variance > best.variance) best = { compound, ...stats };
  });
  return best;
}

function meanDistanceToCluster(feature, cluster, features) {
  let sum = 0;
  let count = 0;
  cluster.forEach((index) => {
    const result = sharedCompoundDistance(feature, features[index]);
    if (!result) return;
    sum += result.distance;
    count += 1;
  });
  return count ? sum / count : Infinity;
}

function buildMissingAwareDivisiveDendrogram(features, indexes = features.map((_, index) => index), depth = 0) {
  if (indexes.length === 1) return { leaf: indexes[0], height: 0 };
  const split = bestSplitCompound(features, indexes);
  if (!split) {
    const mid = Math.ceil(indexes.length / 2);
    return {
      left: buildMissingAwareDivisiveDendrogram(features, indexes.slice(0, mid), depth + 1),
      right: buildMissingAwareDivisiveDendrogram(features, indexes.slice(mid), depth + 1),
      height: indexes.length + depth,
    };
  }
  const observed = indexes
    .map((index) => ({ index, value: features[index].values.get(split.compound) }))
    .filter((item) => Number.isFinite(item.value) && item.value > 0)
    .sort((a, b) => a.value - b.value);
  const median = observed[Math.floor(observed.length / 2)]?.value ?? 0;
  const left = [];
  const right = [];
  const missing = [];
  indexes.forEach((index) => {
    const value = features[index].values.get(split.compound);
    if (!Number.isFinite(value) || value <= 0) missing.push(index);
    else if (value <= median) left.push(index);
    else right.push(index);
  });
  missing.forEach((index) => {
    const leftDistance = left.length ? meanDistanceToCluster(features[index], left, features) : Infinity;
    const rightDistance = right.length ? meanDistanceToCluster(features[index], right, features) : Infinity;
    if (leftDistance === Infinity && rightDistance === Infinity) (left.length <= right.length ? left : right).push(index);
    else (leftDistance <= rightDistance ? left : right).push(index);
  });
  if (!left.length || !right.length) {
    const mid = Math.ceil(indexes.length / 2);
    left.splice(0, left.length, ...indexes.slice(0, mid));
    right.splice(0, right.length, ...indexes.slice(mid));
  }
  return {
    left: buildMissingAwareDivisiveDendrogram(features, left, depth + 1),
    right: buildMissingAwareDivisiveDendrogram(features, right, depth + 1),
    height: Math.max(split.variance, 0.0001) * indexes.length + depth,
  };
}

function dendrogramBoundaryDistance(features, leftIndex, rightIndex) {
  const result = sharedCompoundDistance(features[leftIndex], features[rightIndex]);
  return result?.distance ?? Infinity;
}

function optimizeDendrogramLeafOrder(root, features) {
  if (!root) return [];
  if (root.leaf != null) return [root.leaf];
  const leftOrder = optimizeDendrogramLeafOrder(root.left, features);
  const rightOrder = optimizeDendrogramLeafOrder(root.right, features);
  if (!leftOrder.length) return rightOrder;
  if (!rightOrder.length) return leftOrder;
  const currentDistance = dendrogramBoundaryDistance(features, leftOrder[leftOrder.length - 1], rightOrder[0]);
  const swappedDistance = dendrogramBoundaryDistance(features, rightOrder[rightOrder.length - 1], leftOrder[0]);
  if (swappedDistance < currentDistance) return [...rightOrder, ...leftOrder];
  return [...leftOrder, ...rightOrder];
}

function rotateOrderAtLargestGap(order, features) {
  if (order.length < 3) return order;
  let bestIndex = -1;
  let bestDistance = -Infinity;
  for (let index = 0; index < order.length; index += 1) {
    const next = (index + 1) % order.length;
    const distance = dendrogramBoundaryDistance(features, order[index], order[next]);
    if (distance > bestDistance) {
      bestDistance = distance;
      bestIndex = index;
    }
  }
  if (bestIndex < 0 || !Number.isFinite(bestDistance)) return order;
  return [...order.slice(bestIndex + 1), ...order.slice(0, bestIndex + 1)];
}

function remapDendrogramLeaves(node, mapLeaf) {
  if (!node) return null;
  if (node.leaf != null) return { ...node, leaf: mapLeaf(node.leaf) };
  return {
    ...node,
    left: remapDendrogramLeaves(node.left, mapLeaf),
    right: remapDendrogramLeaves(node.right, mapLeaf),
  };
}

function orientDendrogramToOrder(node, positionByLeaf) {
  if (!node || node.leaf != null) return node;
  const left = orientDendrogramToOrder(node.left, positionByLeaf);
  const right = orientDendrogramToOrder(node.right, positionByLeaf);
  const leftOrder = rowOrderFromDendrogram(left);
  const rightOrder = rowOrderFromDendrogram(right);
  const leftMin = Math.min(...leftOrder.map((leaf) => positionByLeaf.get(leaf) ?? Number.MAX_SAFE_INTEGER));
  const rightMin = Math.min(...rightOrder.map((leaf) => positionByLeaf.get(leaf) ?? Number.MAX_SAFE_INTEGER));
  return rightMin < leftMin ? { ...node, left: right, right: left } : { ...node, left, right };
}

function mergeOrderedDendrogramBlocks(blocks, heightBase) {
  if (!blocks.length) return null;
  if (blocks.length === 1) return blocks[0];
  const mid = Math.ceil(blocks.length / 2);
  const left = mergeOrderedDendrogramBlocks(blocks.slice(0, mid), heightBase);
  const right = mergeOrderedDendrogramBlocks(blocks.slice(mid), heightBase);
  const leftMembers = left?.members || rowOrderFromDendrogram(left);
  const rightMembers = right?.members || rowOrderFromDendrogram(right);
  return {
    left,
    right,
    members: [...leftMembers, ...rightMembers],
    height: heightBase + leftMembers.length + rightMembers.length,
    bridged: true,
  };
}

function buildParentBlockedCircularDendrogram(siteSegments, usable, usableFeatures) {
  const blocks = [];
  const byOriginalIndex = new Map(usable.map((item, localIndex) => [item.index, { localIndex, feature: item.feature }]));
  for (const segment of siteSegments) {
    const originalIndex = siteSegments.indexOf(segment);
    const usableItem = byOriginalIndex.get(originalIndex);
    if (!usableItem) continue;
    const parentBlockId = segment.parentId || segment.row?.parentId || CIRCULAR_ROOT_ORDER_KEY;
    const last = blocks[blocks.length - 1];
    if (last?.parentBlockId === parentBlockId) {
      last.items.push(usableItem);
    } else {
      blocks.push({ parentBlockId, items: [usableItem] });
    }
  }
  const order = blocks.flatMap((block) => block.items.map((item) => item.localIndex));
  const positionByLeaf = new Map(order.map((leaf, position) => [leaf, position]));
  let bridgedMerges = 0;
  let maxHeight = 0;
  const roots = blocks
    .map((block) => {
      if (block.items.length === 1) return { leaf: block.items[0].localIndex, members: [block.items[0].localIndex], height: 0, parentBlockId: block.parentBlockId };
      const localFeatures = block.items.map((item) => item.feature);
      const dendrogram = buildMissingAwareAgglomerativeDendrogram(localFeatures, { allowBridges: true, minValidPairFraction: 0 });
      bridgedMerges += dendrogram.bridgedMerges || 0;
      maxHeight = Math.max(maxHeight, dendrogram.root?.height || 0);
      const root = remapDendrogramLeaves(dendrogram.root, (localLeaf) => block.items[localLeaf].localIndex);
      return orientDendrogramToOrder(root, positionByLeaf);
    })
    .filter(Boolean);
  const root = mergeOrderedDendrogramBlocks(roots, Math.max(1, maxHeight));
  const orientedRoot = orientDendrogramToOrder(root, positionByLeaf);
  const treeOrder = rowOrderFromDendrogram(orientedRoot);
  return { root: orientedRoot, order: treeOrder, bridgedMerges: bridgedMerges + Math.max(0, roots.length - 1), blockCount: blocks.length };
}

async function computeCircularVisibleSiteDendrogram(siteSegments, key) {
  const rows = siteSegments.map((segment) => segment.row).filter(Boolean);
  if (rows.length < MIN_VISIBLE_DENDROGRAM_SITES) return { key, root: null, order: [], status: "too few visible sites" };
  const features = await Promise.all(rows.map(circularSiteFeature));
  const usable = features
    .map((feature, index) => ({ feature, index }))
    .filter((item) => item.feature.values.size >= MIN_CLUSTER_SHARED_COMPOUNDS);
  if (usable.length < MIN_VISIBLE_DENDROGRAM_SITES) {
    return { key, root: null, order: [], status: `only ${usable.length} visible sites have >= ${MIN_CLUSTER_SHARED_COMPOUNDS} observed compounds` };
  }
  const usableFeatures = usable.map((item) => item.feature);
  const dendrogram = buildParentBlockedCircularDendrogram(siteSegments, usable, usableFeatures);
  const root = dendrogram.root;
  const order = dendrogram.order.map((localIndex) => usable[localIndex].index);
  const completeCount = completeCompoundCount(usableFeatures);
  return {
    key,
    root,
    order,
    status: completeCount < MIN_COMPLETE_CLUSTER_COMPOUNDS
      ? `center dendrogram: parent-blocked agglomerative visible-site clustering; complete compounds shared by all visible sites=${completeCount}; parent blocks=${dendrogram.blockCount}; bridged joins=${dendrogram.bridgedMerges}`
      : `center dendrogram: parent-blocked agglomerative visible-site clustering; complete compounds shared by all visible sites=${completeCount}; parent blocks=${dendrogram.blockCount}; bridged joins=${dendrogram.bridgedMerges}`,
  };
}

function prepareCircularVisibleSiteDendrogram(siteSegments) {
  const key = circularVisibleSiteDendrogramKey(siteSegments);
  if (state.circularVisibleSiteDendrogramKey === key && state.circularVisibleSiteDendrogram) return state.circularVisibleSiteDendrogram;
  if (state.circularVisibleSiteDendrogramKey === key && state.circularVisibleSiteDendrogramPromise) return null;
  state.circularVisibleSiteDendrogramKey = key;
  state.circularVisibleSiteDendrogram = null;
  state.circularVisibleSiteDendrogramPromise = computeCircularVisibleSiteDendrogram(siteSegments, key)
    .then((result) => {
      if (state.circularVisibleSiteDendrogramKey === key) {
        state.circularVisibleSiteDendrogram = result;
        state.circularVisibleSiteDendrogramPromise = null;
        renderCircularCompoundAtlas();
      }
      return result;
    })
    .catch((err) => {
      console.warn(`Unable to cluster visible ligandability sites: ${err.message}`);
      if (state.circularVisibleSiteDendrogramKey === key) {
        state.circularVisibleSiteDendrogram = { key, root: null, order: [], status: "center dendrogram unavailable" };
        state.circularVisibleSiteDendrogramPromise = null;
      }
    });
  return null;
}

function drawCircularCenterDendrogram(ctx, dendrogram, size) {
  if (!dendrogram?.root || !dendrogram.order?.length) return;
  const cx = size / 2;
  const cy = size / 2;
  const outerRadius = size * Math.max(0.02, CIRCULAR_BANDS.site[0] - 0.017);
  const rootRadius = 3;
  const dendrogramLeafDistance = (node) => {
    if (!node || node.leaf != null) return 0;
    return 1 + Math.max(dendrogramLeafDistance(node.left), dendrogramLeafDistance(node.right));
  };
  const maxLeafDistance = Math.max(1, dendrogramLeafDistance(dendrogram.root));
  const radiusForLeafDistance = (distance) => rootRadius + (1 - distance / maxLeafDistance) * (outerRadius - rootRadius);
  const angleByLeaf = new Map(
    dendrogram.order.map((leafIndex, position) => [
      leafIndex,
      -Math.PI / 2 + ((position + 0.5) / Math.max(1, dendrogram.order.length)) * Math.PI * 2,
    ]),
  );
  const normalizeRadialAngle = (angle, reference = 0) => {
    let adjusted = angle;
    while (adjusted - reference > Math.PI) adjusted -= Math.PI * 2;
    while (adjusted - reference < -Math.PI) adjusted += Math.PI * 2;
    return adjusted;
  };
  const radialPoint = (angle, radius) => ({ x: cx + Math.cos(angle) * radius, y: cy + Math.sin(angle) * radius });
  const radialDendrogramPath = (from, to) => {
    if (Math.hypot(to.x - from.x, to.y - from.y) < MIN_DENDROGRAM_STROKE_PX) return;
    ctx.beginPath();
    ctx.moveTo(from.x, from.y);
    ctx.lineTo(to.x, to.y);
    ctx.stroke();
  };
  const drawArc = (radius, startAngle, endAngle) => {
    if (Math.abs(endAngle - startAngle) * radius < MIN_DENDROGRAM_STROKE_PX) return;
    ctx.beginPath();
    ctx.arc(cx, cy, radius, startAngle, endAngle, endAngle < startAngle);
    ctx.stroke();
  };
  const drawRadialTree = (node) => {
    if (node.leaf != null) {
      const angle = angleByLeaf.get(node.leaf) ?? -Math.PI / 2;
      return { ...radialPoint(angle, outerRadius), angle, radius: outerRadius, startAngle: angle, endAngle: angle, leafDistance: 0 };
    }
    const left = drawRadialTree(node.left);
    const right = drawRadialTree(node.right);
    const rightStart = normalizeRadialAngle(right.startAngle, left.startAngle);
    const rightEnd = normalizeRadialAngle(right.endAngle, left.startAngle);
    const startAngle = Math.min(left.startAngle, rightStart);
    const endAngle = Math.max(left.endAngle, rightEnd);
    const angle = (startAngle + endAngle) / 2;
    const radius = radiusForLeafDistance(dendrogramLeafDistance(node));
    const here = radialPoint(angle, radius);
    const leftJoin = radialPoint(left.angle, radius);
    const rightAngle = normalizeRadialAngle(right.angle, left.angle);
    const rightJoin = radialPoint(rightAngle, radius);
    const leafDistance = dendrogramLeafDistance(node);
    radialDendrogramPath(left, leftJoin);
    radialDendrogramPath(right, rightJoin);
    drawArc(radius, left.angle, rightAngle);
    return { ...here, angle, radius, startAngle, endAngle, leafDistance };
  };
  ctx.save();
  ctx.strokeStyle = "#2f3136";
  ctx.lineWidth = 0.35;
  const treeRoot = drawRadialTree(dendrogram.root);
  ctx.beginPath();
  ctx.moveTo(cx, cy);
  ctx.lineTo(treeRoot.x, treeRoot.y);
  ctx.stroke();
  ctx.beginPath();
  ctx.arc(cx, cy, 2.2, 0, Math.PI * 2);
  ctx.fillStyle = "rgba(47, 49, 54, 0.5)";
  ctx.fill();
  ctx.restore();
}

function drawFlatHeatmap(canvas, rows, compounds, cells, options = {}) {
  const dendrogramWidth = options.dendrogram ? 72 : 0;
  const labelWidth = 180;
  const headerHeight = 96;
  const cellSize = 12;
  const heatmapX = dendrogramWidth + labelWidth;
  const width = heatmapX + Math.max(1, compounds.length) * cellSize;
  const height = headerHeight + Math.max(1, rows.length) * cellSize;
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(height * dpr);
  canvas.style.width = `${width}px`;
  canvas.style.height = `${height}px`;
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#fff";
  ctx.fillRect(0, 0, width, height);
  ctx.font = "11px Helvetica, Arial, sans-serif";
  if (options.dendrogram) drawCircularDendrogram(ctx, options.dendrogram, rows, 4, headerHeight, dendrogramWidth - 8, cellSize);
  const labelTargets = [];
  rows.forEach((row, index) => {
    const label = String(row.label || row.id).slice(0, 28);
    const rowId = row.level === "site" ? Number(row.siteRow ?? row.rowId) : NaN;
    const isLinkedSite = row.level === "site" && Number.isFinite(rowId) && rowId >= 0;
    ctx.fillStyle = isLinkedSite ? MORANDI.blueDark : "#1f3136";
    ctx.fillText(label, dendrogramWidth + 6, headerHeight + index * cellSize + 10);
    if (isLinkedSite) {
      labelTargets.push({
        x: dendrogramWidth,
        y: headerHeight + index * cellSize,
        width: labelWidth,
        height: cellSize,
        row,
        rowId,
      });
    }
  });
  ctx.fillStyle = "#1f3136";
  compounds.forEach((compound, index) => {
    if (index % Math.max(1, Math.ceil(compounds.length / 80)) !== 0) return;
    ctx.save();
    ctx.translate(heatmapX + index * cellSize + 8, headerHeight - 8);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(compound.id, 0, 0);
    ctx.restore();
  });
  const hitTargets = [];
  for (const [rowIndex, compoundIndex, maxR, hitCount, pValue] of cells) {
    ctx.fillStyle = circularColor(maxR);
    ctx.fillRect(heatmapX + compoundIndex * cellSize, headerHeight + rowIndex * cellSize, cellSize, cellSize);
    const row = rows[rowIndex] || {};
    const compound = compounds[compoundIndex] || {};
    const siteMode = row.level === "site";
    hitTargets.push({
      x: heatmapX + compoundIndex * cellSize,
      y: headerHeight + rowIndex * cellSize,
      width: cellSize,
      height: cellSize,
      payload: {
        mode: siteMode ? "site" : "aggregate",
        hit: [compound.id, Number(maxR || 0) / 1000, pValue],
        maxR: Number(maxR || 0) / 1000,
        hitCount: Number(hitCount || 0),
        rowId: siteMode ? (row.siteRow ?? row.rowId) : null,
      },
    });
  }
  canvas._circularHeatmapTargets = hitTargets;
  canvas._circularHeatmapLabelTargets = labelTargets;
  canvas.onmousemove = (event) => {
    const rect = canvas.getBoundingClientRect();
    const x = ((event.clientX - rect.left) / Math.max(1, rect.width)) * width;
    const y = ((event.clientY - rect.top) / Math.max(1, rect.height)) * height;
    const target = hitTargets.find((item) => x >= item.x && x <= item.x + item.width && y >= item.y && y <= item.y + item.height);
    const labelTarget = labelTargets.find((item) => x >= item.x && x <= item.x + item.width && y >= item.y && y <= item.y + item.height);
    canvas.style.cursor = labelTarget ? "pointer" : "crosshair";
    if (!target) {
      hideMoleculeTooltip();
      return;
    }
    showCircularHeatmapTooltip(target.payload, { event });
  };
  canvas.onclick = (event) => {
    const rect = canvas.getBoundingClientRect();
    const x = ((event.clientX - rect.left) / Math.max(1, rect.width)) * width;
    const y = ((event.clientY - rect.top) / Math.max(1, rect.height)) * height;
    const labelTarget = labelTargets.find((item) => x >= item.x && x <= item.x + item.width && y >= item.y && y <= item.y + item.height);
    if (labelTarget) {
      openCircularSiteInSiteView(labelTarget.row);
      return;
    }
    const target = hitTargets.find((item) => x >= item.x && x <= item.x + item.width && y >= item.y && y <= item.y + item.height);
    if (target) showCircularHeatmapTooltip(target.payload, { event });
  };
  canvas.onmouseleave = () => {
    canvas.style.cursor = "crosshair";
    hideMoleculeTooltip();
  };
}

function findCircularRowContext(data, nodeId) {
  const topIndex = (data.rows || []).findIndex((row) => row.id === nodeId);
  if (topIndex >= 0) {
    const row = data.rows[topIndex];
    return {
      row,
      rowIndex: topIndex,
      parentRows: data.rows || [],
      parentCells: data.cells || [],
      drilldown: data.drilldowns?.[nodeId] || null,
    };
  }
  for (const drilldown of Object.values(data.drilldowns || {})) {
    const rows = drilldown.rows || [];
    const rowIndex = rows.findIndex((row) => row.id === nodeId);
    if (rowIndex < 0) continue;
    const row = rows[rowIndex];
    return {
      row,
      rowIndex,
      parentRows: rows,
      parentCells: drilldown.cells || [],
      drilldown: data.drilldowns?.[nodeId] || null,
    };
  }
  return null;
}

function circularHeatmapPayload(data, nodeId, segmentContext = null) {
  const context = segmentContext?.row
    ? {
        row: segmentContext.row,
        rowIndex: segmentContext.parentRowIndex ?? -1,
        parentRows: [],
        parentCells: segmentContext.parentCells || [],
        drilldown: data.drilldowns?.[segmentContext.row.id] || null,
      }
    : findCircularRowContext(data, nodeId);
  if (!context) return null;
  const allCompounds = data.compounds || [];
  let rows = [];
  let cells = [];
  if (context.drilldown?.rows?.length) {
    const aggregateMemo = new Map();
    const drilldownCells = circularDrilldownCellsForRow(data, context.row, context.parentCells || [], context.rowIndex);
    (context.drilldown.rows || []).forEach((row, rowIndex) => {
      const aggregate = circularAggregateCellsForRow(data, row, drilldownCells, rowIndex, aggregateMemo);
      if (!circularAggregateHasHits(aggregate)) return;
      const filteredRowIndex = rows.length;
      rows.push(row);
      cells.push(...circularCellsFromAggregate(aggregate, filteredRowIndex));
    });
  } else {
    const aggregate = circularAggregateCellsForRow(data, context.row, context.parentCells || [], context.rowIndex);
    if (circularAggregateHasHits(aggregate)) {
      rows = [context.row];
      cells = circularCellsFromAggregate(aggregate, 0);
    }
  }
  const hitCompoundIndexes = [
    ...new Set(
      cells
        .filter(([, , , hitCount]) => Number(hitCount || 0) > 0)
        .map(([, compoundIndex]) => compoundIndex)
    ),
  ].sort((a, b) => a - b);
  const indexMap = new Map(hitCompoundIndexes.map((compoundIndex, index) => [compoundIndex, index]));
  const hitCompounds = hitCompoundIndexes.map((compoundIndex) => allCompounds[compoundIndex]).filter(Boolean);
  const filteredCells = cells
    .filter(([, compoundIndex, , hitCount]) => indexMap.has(compoundIndex) && Number(hitCount || 0) > 0)
    .map(([rowIndex, compoundIndex, maxR, hitCount]) => [rowIndex, indexMap.get(compoundIndex), maxR, hitCount]);
  return {
    title: context.row.label || context.row.id || "Detail",
    meta: `${hitCompounds.length.toLocaleString()}/${allCompounds.length.toLocaleString()} compounds`,
    rows,
    compounds: hitCompounds,
    cells: filteredCells,
  };
}

async function renderCircularDrilldown(nodeId, segmentContext = null) {
  const data = circularDataset();
  if (!data) return;
  const payload = circularHeatmapPayload(data, nodeId, segmentContext);
  if (!payload) {
    openSummaryDetailModal("No lower-level heatmap", "", '<p class="muted">No lower-level heatmap is available for this sector.</p>');
    return;
  }
  await openCircularHeatmapModal(
    payload.title,
    payload.meta,
    payload.rows,
    payload.compounds,
    payload.cells,
  );
}

function rowMatchesTarget(row, targetList, customGenes) {
  if (!row) return false;
  if (targetList === "all") return true;
  if (targetList === "cancer") return isCancerGene(row.gene);
  if (CONTACT_TYPES.includes(targetList)) return row.contactTypes.includes(targetList);
  if (targetList === "custom") return customGenes.has(String(row.gene).toUpperCase());
  return true;
}

async function mapWithConcurrency(items, concurrency, worker) {
  const results = [];
  let index = 0;
  const workers = Array.from({ length: Math.min(concurrency, items.length) }, async () => {
    while (index < items.length) {
      const currentIndex = index;
      index += 1;
      results[currentIndex] = await worker(items[currentIndex]);
    }
  });
  await Promise.all(workers);
  return results;
}

function siteTargetListValue() {
  return el("siteTargetList")?.value || "all";
}

function siteTargetCustomGenes() {
  return customGeneSet(el("siteCustomGenes")?.value);
}

function isContactTargetList(targetList) {
  return CONTACT_TYPES.includes(targetList);
}

function rowHasContactType(row, contacts, contactType) {
  if (!row || !contacts || !contactType) return false;
  return (row.positions || []).some((position) =>
    (contacts[String(position)] || []).some((entry) => (entry.contacts || []).includes(contactType))
  );
}

async function contactCandidateUniprots(contactType) {
  const payload = await loadGlobalProteome();
  const points = payload?.points || [];
  if (!points.length) return null;
  return new Set(points.filter((point) => (point.contactTypes || []).includes(contactType)).map((point) => point.uniprot));
}

async function loadContactFileIndex() {
  if (state.contactFileIndex) return state.contactFileIndex;
  if (state.contactFileIndexPromise) return state.contactFileIndexPromise;
  state.contactFileIndexPromise = fetchJson("contacts/index.json")
    .then((items) => {
      state.contactFileIndex = new Set(items || []);
      return state.contactFileIndex;
    })
    .catch(() => null)
    .finally(() => {
      state.contactFileIndexPromise = null;
    });
  return state.contactFileIndexPromise;
}

async function buildContactTargetRows(contactType) {
  const candidateUniprots = await contactCandidateUniprots(contactType);
  const contactFileIndex = await loadContactFileIndex();
  const rowsByUniprot = new Map();
  for (const row of state.catalog) {
    if (Number(row.observedCount || 0) <= 0) continue;
    if (!Number.isFinite(Number(row.maxR))) continue;
    if (candidateUniprots && !candidateUniprots.has(row.uniprot)) continue;
    if (contactFileIndex && !contactFileIndex.has(row.uniprot)) continue;
    if (!rowsByUniprot.has(row.uniprot)) rowsByUniprot.set(row.uniprot, []);
    rowsByUniprot.get(row.uniprot).push(row);
  }
  const rowGroups = [...rowsByUniprot.values()];
  const matchedGroups = await mapWithConcurrency(rowGroups, 8, async (rows) => {
    const contacts = await loadContactsForRow(rows[0]);
    return rows.filter((row) => rowHasContactType(row, contacts, contactType));
  });
  return matchedGroups.flat();
}

async function contactTargetRows(contactType) {
  const key = `${state.datasetKey}:${contactType}`;
  if (!state.contactTargetRowCache.has(key)) {
    state.contactTargetRowCache.set(key, buildContactTargetRows(contactType));
  }
  return state.contactTargetRowCache.get(key);
}

async function contactTargetRowSet(contactType) {
  const key = `${state.datasetKey}:${contactType}`;
  if (!state.contactTargetRowSetCache.has(key)) {
    state.contactTargetRowSetCache.set(key, contactTargetRows(contactType).then((rows) => new Set(rows.map((row) => Number(row.i)))));
  }
  return state.contactTargetRowSetCache.get(key);
}

async function siteTargetRows(targetList, customGenes) {
  if (isContactTargetList(targetList)) return contactTargetRows(targetList);
  return state.catalog.filter((row) => rowMatchesTarget(row, targetList, customGenes));
}

async function loadGeneSites(gene) {
  if (!gene) return null;
  if (state.geneSiteCache.has(gene)) return state.geneSiteCache.get(gene);
  const path = state.dataset.geneFiles[gene];
  if (!path) return null;
  try {
    const payload = await fetchJson(path.replace(/%/g, "%25"));
    state.geneSiteCache.set(gene, payload);
    return payload;
  } catch (err) {
    state.geneSiteCache.set(gene, null);
    console.warn(`Unable to load target data for ${gene}: ${err.message}`);
    return null;
  }
}

async function siteSummaryForRow(row) {
  if (!row) return null;
  const payload = await loadGeneSites(row.gene);
  return payload?.sites.find((site) => site.row === row.i || site.label === row.label) || null;
}

function filteredHits(summary, ids, maxHitCount = Infinity) {
  if (!summary) return [];
  return summary.hits
    .filter((hit) => drugHitCount(hit[0]) <= maxHitCount)
    .filter((hit) => qualityPassFromHit(hit, ids))
    .sort((a, b) => b[1] - a[1]);
}

function siteFilterSignature(activeOnly = checkboxChecked("siteActiveOnly")) {
  return [
    state.datasetKey,
    siteTargetListValue(),
    [...siteTargetCustomGenes()].sort().join(","),
    activeOnly ? 1 : 0,
    el("siteMaxHits")?.value || "all",
    checkboxChecked("siteSigOnly") ? 1 : 0,
    checkboxChecked("siteHideVariance") ? 1 : 0,
    numericValue("siteMinSn", 0),
  ].join("|");
}

function renderHitBar(targetId, hits, threshold, rowId = null) {
  if (targetId === "siteBar") clearPersistentHitTooltip();
  const isSiteBar = targetId === "siteBar";
  const topHits = hits.slice(0, 20);
  if (!topHits.length) {
    plotMessage(targetId, "No compounds meet the selected filters.");
    hideMoleculeTooltip();
    return;
  }
  const hoverPayload = topHits.map((hit) => ({ hit, rowId }));
  const trace = {
    type: "bar",
    x: topHits.map((hit) => hit[0]),
    y: topHits.map((hit) => hit[1]),
    marker: {
      color: topHits.map((hit) => hit[2] ?? 1),
      colorscale: PVALUE_BLUE_SCALE,
      cmin: 0,
      cmax: 1,
      colorbar: {
        title: { text: "P-value", side: "right" },
        len: 0.82,
        ...(isSiteBar ? { x: 0.995, xanchor: "left" } : {}),
      },
    },
    customdata: hoverPayload,
    hoverinfo: "none",
  };
  const syncedSize = isSiteBar ? syncSiteBarHeight() : null;
  const plotted = safePlot(
    targetId,
    [trace],
    mergeAxes(baseLayout({
      ...(syncedSize ? { height: syncedSize.height, width: syncedSize.width } : {}),
      ...(isSiteBar ? { margin: { l: 48, r: 66, t: 18, b: 88, autoexpand: false } } : {}),
      shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: threshold, y1: threshold, line: { dash: "dash", color: "#666" } }],
    }), { title: `Top ${topHits.length} Compounds (/${hits.length})`, tickangle: 42, automargin: true }, { title: "R" })
  );
  if (!plotted) return;
  const plot = el(targetId);
  if (typeof plot.on === "function") {
    plot.removeAllListeners?.("plotly_hover");
    plot.removeAllListeners?.("plotly_unhover");
    plot.removeAllListeners?.("plotly_click");
    plot.on("plotly_hover", (event) => showTransientHitTooltip(event.points[0].customdata, event));
    plot.on("plotly_unhover", hideTransientHitTooltip);
    plot.on("plotly_click", (event) => pinHitTooltip(event.points[0].customdata, event));
  }
  requestAnimationFrame(() => bindBarHover(targetId, hoverPayload));
}

function syncSiteBarHeight() {
  const layout = document.querySelector("#site .site-analysis-layout");
  const settingsPanel = layout?.querySelector("aside.panel");
  const chartPanel = layout?.querySelector("section.panel");
  const bar = el("siteBar");
  if (!layout || !settingsPanel || !chartPanel || !bar) return null;
  const settingsHeight = settingsPanel.getBoundingClientRect().height;
  if (!Number.isFinite(settingsHeight) || settingsHeight <= 0) return null;
  const panelStyles = getComputedStyle(chartPanel);
  const verticalInset =
    Number.parseFloat(panelStyles.paddingTop || "0") +
    Number.parseFloat(panelStyles.paddingBottom || "0") +
    Number.parseFloat(panelStyles.borderTopWidth || "0") +
    Number.parseFloat(panelStyles.borderBottomWidth || "0");
  const horizontalInset =
    Number.parseFloat(panelStyles.paddingLeft || "0") +
    Number.parseFloat(panelStyles.paddingRight || "0") +
    Number.parseFloat(panelStyles.borderLeftWidth || "0") +
    Number.parseFloat(panelStyles.borderRightWidth || "0");
  const plotHeight = Math.max(320, Math.round(settingsHeight - verticalInset));
  const contentWidth = Math.max(320, Math.round(chartPanel.getBoundingClientRect().width - horizontalInset));
  chartPanel.style.height = `${Math.round(settingsHeight)}px`;
  bar.style.setProperty("--site-bar-height", `${plotHeight}px`);
  bar.style.setProperty("--site-bar-width", `${contentWidth}px`);
  bar.style.height = `${plotHeight}px`;
  bar.style.width = `${contentWidth}px`;
  return { height: plotHeight, width: contentWidth };
}

function resizeSyncedSiteBar() {
  const size = syncSiteBarHeight();
  if (!size) return;
  if (window.Plotly?.relayout) Plotly.relayout(el("siteBar"), size);
  schedulePlotResize("siteBar");
}

async function loadCompoundFromDataset(datasetKey, dataset, drug) {
  const key = `${datasetKey}:${drug}`;
  if (state.compoundCache.has(key)) return state.compoundCache.get(key);
  if (state.compoundLoadPromises.has(key)) return state.compoundLoadPromises.get(key);
  const promise = Promise.all((dataset.compoundParts[drug] || []).map((path) => fetchJson(path)))
    .then((payloads) => {
      const rows = payloads.flatMap((payload) => payload.rows || []);
      state.compoundCache.set(key, rows);
      return rows;
    })
    .finally(() => {
      state.compoundLoadPromises.delete(key);
    });
  state.compoundLoadPromises.set(key, promise);
  return promise;
}

async function loadCompound(drug) {
  return loadCompoundFromDataset(state.datasetKey, state.dataset, drug);
}

function colorForCompound(row, catalogRow, mode, threshold, customGenes) {
  if (mode === "custom") return customGenes.has(String(catalogRow.gene).toUpperCase()) ? "highlight" : "non";
  if (mode === "cancer") return isCancerGene(catalogRow.gene) ? "highlight" : "non";
  if (mode === "contacts") return catalogRow.contactTypes.length ? "highlight" : "non";
  if (mode.startsWith("contact:")) return catalogRow.contactTypes.includes(mode.split(":")[1]) ? "highlight" : "non";
  if (mode === "pvalue") return "pvalue";
  return row[1] > threshold ? "high" : "non";
}

function compoundPlotRows(rows, applyQualityFilters = true) {
  const mode = el("compoundColorMode").value;
  const threshold = numericValue("compoundThreshold", 2);
  const customGenes = customGeneSet(el("compoundCustomGenes").value);
  return rows
    .filter((row) =>
      !applyQualityFilters ||
      qualityPassFromCompound(row, {
        sigOnly: "compoundSigOnly",
        hideVariance: "compoundHideVariance",
        minSn: "compoundMinSn",
      })
    )
    .map((row, siteIndex) => {
      const catalogRow = state.rowById.get(row[0]);
      return { row, catalogRow, siteIndex, color: colorForCompound(row, catalogRow, mode, threshold, customGenes) };
    })
    .filter((item) => item.catalogRow)
    .map((item, _, filtered) => item);
}

function rankCompoundItems(items) {
  const ranked = [...items].sort((a, b) => Number(b.row[1] ?? -Infinity) - Number(a.row[1] ?? -Infinity));
  ranked.forEach((item, idx) => {
    item.rRank = idx + 1;
    item.totalSites = items.length;
  });
  return items;
}

function compoundPointPayload(item) {
  return {
    rowId: item.row[0],
    label: item.catalogRow.label,
    uniprot: item.catalogRow.uniprot,
    description: item.catalogRow.description,
    r: item.row[1],
    p: item.row[2],
    sn: item.row[5],
    rRank: item.rRank,
    totalSites: item.totalSites,
  };
}

function traceForCompoundItems(items, color, name) {
  const filtered = items.filter((item) => item.color === color);
  const isMuted = color === "non";
  return {
    type: "scattergl",
    mode: "markers",
    name,
    x: filtered.map((item) => item.siteIndex),
    y: filtered.map((item) => item.row[1]),
    text: filtered.map((item) => item.catalogRow.label),
    customdata: filtered.map(compoundPointPayload),
    marker: { color: isMuted ? PROTEOME_MUTED_GREY : color === "highlight" ? MORANDI.red : MORANDI.blueDark, opacity: isMuted ? 1 : 0.9, size: 5 },
    hoverinfo: "none",
  };
}

function pValueGradientTrace(items) {
  return {
    type: "scattergl",
    mode: "markers",
    name: "P-value",
    x: items.map((item) => item.siteIndex),
    y: items.map((item) => item.row[1]),
    text: items.map((item) => item.catalogRow.label),
    customdata: items.map(compoundPointPayload),
    marker: {
      color: items.map((item) => item.row[2] ?? 1),
      colorscale: PVALUE_BLUE_SCALE,
      cmin: 0,
      cmax: 1,
      showscale: true,
      colorbar: { title: { text: "P-value", side: "right" }, len: 0.82 },
      opacity: 0.82,
      size: 5,
    },
    hoverinfo: "none",
  };
}

function pValueVolcanoTrace(items) {
  return {
    type: "scattergl",
    mode: "markers",
    name: "P-value",
    x: items.map((item) => item.row[3]),
    y: items.map((item) => item.row[4]),
    text: items.map((item) => item.catalogRow.label),
    customdata: items.map(compoundPointPayload),
    marker: {
      color: items.map((item) => item.row[2] ?? 1),
      colorscale: PVALUE_BLUE_SCALE,
      cmin: 0,
      cmax: 1,
      showscale: true,
      colorbar: { title: { text: "P-value", side: "right" }, len: 0.82 },
      opacity: 0.82,
      size: 5,
    },
    hoverinfo: "none",
  };
}

async function renderCompound() {
  updateConditionalFields();
  const drug = el("compoundSelect").value;
  setMolecule("compoundMolecule", drug);
  if (!drug) {
    plotMessage("compoundScatter", "No compounds meet the selected filters.");
    plotMessage("compoundVolcano", "No compounds meet the selected filters.");
    return;
  }
  const rows = await loadCompound(drug);
  const items = rankCompoundItems(compoundPlotRows(rows, true));
  const volcanoItemsBase = items;
  const threshold = numericValue("compoundThreshold", 2);
  const labelCount = numericValue("compoundLabels", 5);
  const colorMode = el("compoundColorMode").value;
  const topLabels = items
    .sort((a, b) => b.row[1] - a.row[1])
    .slice(0, labelCount);
  const scatterX = items.map((item) => item.siteIndex).filter((value) => Number.isFinite(Number(value)));
  const scatterX0 = scatterX.length ? Math.min(...scatterX) : 0;
  const scatterX1 = scatterX.length ? Math.max(...scatterX) : 1;

  safePlot(
    "compoundScatter",
    colorMode === "pvalue"
      ? [pValueGradientTrace(items)]
      : [traceForCompoundItems(items, "non", "Other"), traceForCompoundItems(items, "high", "Above threshold"), traceForCompoundItems(items, "highlight", "Highlighted")],
    baseLayout({
      height: 380,
      margin: { l: 48, r: 8, t: 18, b: 48 },
      xaxis: { title: "Tyrosine Sites", range: [scatterX0, scatterX1], showgrid: false, zeroline: false },
      yaxis: { title: "R", showgrid: false, zeroline: false },
      showlegend: false,
      shapes: [{ type: "line", xref: "x", x0: scatterX0, x1: scatterX1, y0: threshold, y1: threshold, line: { dash: "dash", color: "#666" } }],
      annotations: topLabels.map((item) => ({
        x: item.siteIndex,
        y: item.row[1],
        text: item.catalogRow.label,
        showarrow: true,
        arrowhead: 0,
        arrowwidth: 1,
        arrowsize: 0.5,
        ax: 10,
        ay: -12,
        font: { size: 10 },
      })),
    }),
    { editable: true, edits: { annotationPosition: true } }
  );

  const volcanoItems = volcanoItemsBase.map((item) => ({ ...item, sig: item.row[1] > threshold && item.row[4] > SIG_NEGLOG10 }));
  const volcanoTraces = colorMode === "pvalue" ? [pValueVolcanoTrace(volcanoItems)] : ["non", "high", "highlight"].map((color) => {
    const filtered = volcanoItems.filter((item) => item.color === color);
    const isMuted = color === "non";
    return {
      type: "scattergl",
      mode: "markers",
      name: color,
      x: filtered.map((item) => item.row[3]),
      y: filtered.map((item) => item.row[4]),
      text: filtered.map((item) => item.catalogRow.label),
      customdata: filtered.map(compoundPointPayload),
      marker: { color: isMuted ? PROTEOME_MUTED_GREY : color === "highlight" ? MORANDI.red : MORANDI.blueDark, opacity: isMuted ? 1 : 0.88, size: 5 },
      hoverinfo: "none",
    };
  });
  const volcanoLabels = volcanoItems
    .sort((a, b) => b.row[1] - a.row[1])
    .slice(0, labelCount);
  safePlot(
    "compoundVolcano",
    volcanoTraces,
    baseLayout({
      height: 380,
      xaxis: { title: "log\u2082R", range: [-1, 2], showgrid: false, zeroline: false },
      yaxis: { title: "-log\u2081\u2080P", showgrid: false, zeroline: false },
      showlegend: false,
      shapes: [
        { type: "line", x0: Math.log2(Math.max(threshold, 0.00001)), x1: Math.log2(Math.max(threshold, 0.00001)), yref: "paper", y0: 0, y1: 1, line: { dash: "dash", color: "#666" } },
        { type: "line", xref: "paper", x0: 0, x1: 1, y0: SIG_NEGLOG10, y1: SIG_NEGLOG10, line: { dash: "dash", color: "#666" } },
      ],
      annotations: volcanoLabels.map((item) => ({
        x: item.row[3],
        y: item.row[4],
        text: item.catalogRow.label,
        showarrow: true,
        arrowhead: 0,
        arrowwidth: 1,
        arrowsize: 0.45,
        ax: 8,
        ay: -10,
        font: { size: 10 },
      })),
    }),
    { editable: true, edits: { annotationPosition: true } }
  );
  bindCompoundPointHover("compoundScatter");
  bindCompoundPointHover("compoundVolcano");
}

function bindCompoundPointHover(targetId) {
  const plot = el(targetId);
  if (typeof plot.on !== "function") return;
  plot.removeAllListeners?.("plotly_hover");
  plot.removeAllListeners?.("plotly_unhover");
  plot.on("plotly_hover", (event) => showCompoundPointTooltip(event.points[0].customdata, event));
  plot.on("plotly_unhover", hideMoleculeTooltip);
}

async function computeFilteredSitesForGene(gene, activeOnly) {
  const maxHits = hitCountLimit("siteMaxHits");
  const targetList = siteTargetListValue();
  const customGenes = siteTargetCustomGenes();
  const targetRowSet = isContactTargetList(targetList) ? await contactTargetRowSet(targetList) : null;
  const payload = await loadGeneSites(gene);
  return (payload?.sites || [])
    .filter((site) => {
      const row = state.rowById.get(Number(site.row));
      return targetRowSet ? targetRowSet.has(Number(site.row)) : rowMatchesTarget(row, targetList, customGenes);
    })
    .map((site) => {
      if (!activeOnly) return { ...site, filteredMaxR: Number(site.maxR) };
      const hits = filteredHits(site, siteFilterIds(), maxHits);
      const filteredMaxR = hits[0]?.[1] ?? null;
      return { ...site, filteredMaxR };
    })
    .filter((site) => Number.isFinite(Number(site.filteredMaxR)))
    .filter((site) => !activeOnly || Number(site.filteredMaxR) > 2)
    .sort((a, b) => (b.filteredMaxR ?? -1) - (a.filteredMaxR ?? -1) || a.label.localeCompare(b.label));
}

async function filteredSitesForGene(gene, activeOnly) {
  if (!gene) return [];
  const key = `gene:${siteFilterSignature(activeOnly)}:${gene}`;
  if (!state.filteredSitesByGene.has(key)) {
    state.filteredSitesByGene.set(key, computeFilteredSitesForGene(gene, activeOnly));
  }
  return state.filteredSitesByGene.get(key);
}

async function loadActiveGeneFilterIndex() {
  const path = state.dataset?.activeGeneFilterIndex;
  if (!path) return null;
  const key = `${state.datasetKey}:${path}`;
  if (!state.activeGeneFilterIndexCache.has(key)) {
    state.activeGeneFilterIndexCache.set(key, fetchJson(path));
  }
  return state.activeGeneFilterIndexCache.get(key);
}

function activeIndexEntryAtMinSn(entries, minSn, binSize) {
  if (!entries?.length) return null;
  const queryBin = Math.max(0, Math.ceil(Number(minSn || 0) / binSize) * binSize);
  if (queryBin > entries[entries.length - 1][0]) return null;
  let lo = 0;
  let hi = entries.length - 1;
  let best = -1;
  while (lo <= hi) {
    const mid = Math.floor((lo + hi) / 2);
    if (entries[mid][0] <= queryBin) {
      best = mid;
      lo = mid + 1;
    } else {
      hi = mid - 1;
    }
  }
  return best >= 0 ? entries[best] : null;
}

async function activeGeneStatsForCurrentFilters() {
  const index = await loadActiveGeneFilterIndex();
  if (!index) return null;
  const maxHits = el("siteMaxHits")?.value || "all";
  const comboKey = `${maxHits}|${checkboxChecked("siteSigOnly") ? 1 : 0}|${checkboxChecked("siteHideVariance") ? 1 : 0}`;
  const comboRows = index.stats?.[comboKey] || [];
  const stats = new Map();
  for (const [geneIdx, entries] of comboRows) {
    const entry = activeIndexEntryAtMinSn(entries, numericValue("siteMinSn", 0), index.binSize || 10);
    if (!entry) continue;
    stats.set(index.genes[geneIdx], { count: entry[1], maxR: entry[2] });
  }
  return stats;
}

async function refreshSiteControls(preferredGene = null, { allowFallback = false } = {}) {
  const seq = ++state.siteRefreshSeq;
  const activeOnly = checkboxChecked("siteActiveOnly");
  const targetList = siteTargetListValue();
  const customGenes = siteTargetCustomGenes();
  const currentGene = preferredGene || el("siteGeneSelect").value || state.dataset.defaultGene;
  const currentSite = el("siteSelect").value;
  const targetRows = await siteTargetRows(targetList, customGenes);
  let geneStats = activeOnly && targetList === "all" ? await activeGeneStatsForCurrentFilters() : null;
  if (seq !== state.siteRefreshSeq) return;

  if (!geneStats) {
    geneStats = new Map();
    for (const row of targetRows) {
      if (Number(row.observedCount || 0) <= 0) continue;
      if (activeOnly && Number(row.hitCount || 0) <= 0) continue;
      if (!Number.isFinite(Number(row.maxR))) continue;
      const stats = geneStats.get(row.gene) || { count: 0, maxR: -Infinity };
      stats.count += 1;
      stats.maxR = Math.max(stats.maxR, Number(row.maxR));
      geneStats.set(row.gene, stats);
    }
  }

  const geneValues = [...geneStats.keys()].sort((a, b) => String(a).localeCompare(String(b)));
  let selectedGene = geneValues.includes(currentGene) ? currentGene : geneValues[0] || "";
  let sites = selectedGene ? await filteredSitesForGene(selectedGene, activeOnly) : [];
  if (seq !== state.siteRefreshSeq) return;

  if (allowFallback && !sites.length) {
    for (const gene of geneValues) {
      if (gene === selectedGene) continue;
      const fallbackSites = await filteredSitesForGene(gene, activeOnly);
      if (seq !== state.siteRefreshSeq) return;
      if (fallbackSites.length) {
        selectedGene = gene;
        sites = fallbackSites;
        break;
      }
    }
  }
  if (seq !== state.siteRefreshSeq) return;

  const geneOptions = geneValues.map((gene) => {
    const stats = geneStats.get(gene);
    return { value: gene, label: `${gene} (Sites Seen: ${stats.count}, Max R: ${roundLabel(stats.maxR)})` };
  });
  if (seq !== state.siteRefreshSeq) return;
  setOptions(el("siteGeneSelect"), geneOptions.length ? geneOptions : [{ value: "", label: "No genes found" }], selectedGene);
  const selectedSite = sites.some((site) => String(site.row) === currentSite)
    ? currentSite
    : sites[0]
      ? String(sites[0].row)
      : "";
  setOptions(
    el("siteSelect"),
    sites.length
      ? sites.map((site) => ({ value: String(site.row), label: `${site.site} (Max R: ${roundLabel(site.filteredMaxR)})` }))
      : [{ value: "", label: "No sites found" }],
    selectedSite
  );
  await renderSite({ reloadStructures: selectedSite !== currentSite });
}

async function renderSite({ reloadStructures = true } = {}) {
  const row = state.rowById.get(Number(el("siteSelect").value));
  state.currentSiteRow = row || null;
  state.currentSiteSummary = row ? await siteSummaryForRow(row) : null;
  if (!row) {
    plotMessage("siteBar", "No sites meet the selected filters.");
    el("structureViewer").innerHTML = "";
    el("contactViewer").innerHTML = "";
    el("structureMessage").textContent = "No structure is available for the current selection.";
    el("contactMessage").textContent = "";
    updateStructureLinks("", "");
    return;
  }
  const threshold = 2;
  const maxHits = hitCountLimit("siteMaxHits");
  const hits = filteredHits(state.currentSiteSummary, {
    sigOnly: "siteSigOnly",
    hideVariance: "siteHideVariance",
    minSn: "siteMinSn",
  }, maxHits);
  renderHitBar("siteBar", hits, threshold, row?.i ?? null);
  if (reloadStructures) {
    await populatePdbChoices();
  }
  if (reloadStructures && el("site").classList.contains("active")) {
    autoLoadStructures();
  }
}

async function loadContactsForRow(row) {
  if (!row?.uniprot) return {};
  if (state.contactCache.has(row.uniprot)) return state.contactCache.get(row.uniprot);
  try {
    const payload = await fetchJson(state.manifest.contactsPath.replace("{uniprot}", row.uniprot));
    state.contactCache.set(row.uniprot, payload);
    return payload;
  } catch {
    state.contactCache.set(row.uniprot, {});
    return {};
  }
}

async function populatePdbChoices() {
  const row = state.currentSiteRow;
  const contacts = await loadContactsForRow(row);
  const entries = [];
  for (const pos of row?.positions || []) {
    for (const entry of contacts[pos] || []) {
      entries.push({ ...entry, position: pos });
    }
  }
  const seen = new Set();
  state.currentPdbEntries = entries.filter((entry) => {
    const key = `${entry.pdb}:${entry.chain}:${entry.authSeq}`;
    if (seen.has(key)) return false;
    seen.add(key);
    return true;
  }).sort((a, b) => {
    const methodRank = (entry) => {
      const method = String(entry.method || "").toLowerCase();
      if (method.includes("electron") || method.includes("cryo")) return 0;
      if (method.includes("x-ray")) return 1;
      if (method.includes("nmr")) return 2;
      return 3;
    };
    const resA = Number.isFinite(Number(a.resolution)) ? Number(a.resolution) : 999;
    const resB = Number.isFinite(Number(b.resolution)) ? Number(b.resolution) : 999;
    return methodRank(a) - methodRank(b) || resA - resB || Number(a.minDist || 999) - Number(b.minDist || 999) || String(a.pdb).localeCompare(String(b.pdb));
  });
  const structureEntries = state.currentPdbEntries;
  const contactEntries = state.currentPdbEntries.filter((entry) => Number.isFinite(Number(entry.minDist)) && (entry.contacts || []).length);
  const formatEntryLabel = (entry, includeContact = false) => {
    const resolution = Number(entry.resolution);
    const methodResolution = Number.isFinite(resolution)
      ? `${entry.method || "Method"}, ${roundLabel(resolution, 1)} ${ANGSTROM}`
      : `${entry.method || "Method"}`;
    const hasContact = Number.isFinite(Number(entry.minDist)) && (entry.contacts || []).length;
    if (!hasContact || !includeContact) return `${entry.pdb} | ${methodResolution}`;
    const contactText = (entry.contacts || []).join(", ");
    return `${entry.pdb} | ${methodResolution} | ${contactText}, ${roundLabel(entry.minDist, 1)} ${ANGSTROM}`;
  };
  const structureOptions = structureEntries.length
    ? [
        ...structureEntries.map((entry) => ({
          value: `pdb:${state.currentPdbEntries.indexOf(entry)}`,
          label: formatEntryLabel(entry, true),
        })),
        { value: "af", label: "AlphaFold" },
      ]
    : [{ value: "af", label: "AlphaFold" }];
  const contactOptions = contactEntries.length
    ? contactEntries.map((entry) => ({
        value: String(state.currentPdbEntries.indexOf(entry)),
        label: formatEntryLabel(entry, true),
      }))
      : [{ value: "", label: "No interactors found from PDB structures" }];
  const currentStructure = el("structureSelect").value;
  const currentStructureValue = currentStructure.startsWith("pdb:") && structureOptions.some((option) => option.value === currentStructure)
    ? currentStructure
    : structureOptions[0]?.value || "af";
  setOptions(el("structureSelect"), structureOptions, currentStructureValue);
  setOptions(el("contactPdbSelect"), contactOptions, contactOptions[0]?.value || "");
  el("structureMessage").textContent = structureEntries.length ? "Structure will load inline." : "AlphaFold structure will load inline.";
  el("contactMessage").textContent = contactEntries.length ? "Indexed contact residues will load inline." : "";
  updateConditionalFields();
}

function selectedStructureEntry() {
  const value = el("structureSelect").value || "af";
  if (!value.startsWith("pdb:")) return null;
  return state.currentPdbEntries[Number(value.split(":")[1])] || null;
}

function selectedContactEntry() {
  return state.currentPdbEntries[Number(el("contactPdbSelect").value)] || null;
}

function renderMolstar(containerId, messageId, options, selectionData) {
  const container = el(containerId);
  const message = el(messageId);
  container.innerHTML = "";
  if (!window.PDBeMolstarPlugin) {
    message.textContent = "Mol* failed to load. Check network access to the PDBe Mol* CDN.";
    return;
  }
  try {
    message.textContent = "Loading structure...";
    const viewer = new PDBeMolstarPlugin();
    viewer.render(container, {
      hideStructure: ["water"],
      bgColor: { r: 255, g: 255, b: 255 },
      hideControls: true,
      sequencePanel: false,
      pdbeLink: false,
      expanded: false,
      ...options,
    });
    viewer.events.loadComplete.subscribe(() => {
      message.textContent = "";
      try {
        const canvas3d = viewer.plugin.canvas3d;
        if (canvas3d) canvas3d.setProps({ cameraClipping: { far: false } });
        if (selectionData?.length) {
          viewer.visual.select({ data: selectionData });
        }
      } catch (err) {
        message.textContent = "Structure loaded, but residue highlighting was unavailable.";
      }
    });
  } catch (err) {
    message.textContent = "Unable to initialize Mol* for this structure.";
  }
}

function loadPdbViewer(includeNeighbors = false, target = "structure", chosenEntry = null) {
  const entry = chosenEntry || (target === "contact" ? selectedContactEntry() : selectedStructureEntry());
  const messageId = target === "contact" ? "contactMessage" : "structureMessage";
  const viewerId = target === "contact" ? "contactViewer" : "structureViewer";
  if (!entry) {
    el(messageId).textContent = target === "contact" ? "No contact structure is selected." : "No experimental PDB entry is selected.";
    el(viewerId).innerHTML = "";
    return;
  }
  const selection = [];
  if (entry.chain && entry.authSeq != null) {
    selection.push({
      auth_asym_id: entry.chain,
      auth_residue_number: entry.authSeq,
      representationColor: { r: 255, g: 0, b: 0 },
      representation: "ball-and-stick",
      sideChain: true,
      focus: true,
    });
  }
  if (includeNeighbors) {
    for (const n of entry.neighbors || []) {
      selection.push({
        auth_asym_id: n.chain,
        auth_residue_number: n.residue,
        representationColor: { r: 52, g: 103, b: 177 },
        representation: "ball-and-stick",
        sideChain: true,
        focus: false,
      });
    }
  }
  renderMolstar(
    viewerId,
    messageId,
    { moleculeId: String(entry.pdb).toLowerCase(), visualStyle: "cartoon" },
    selection
  );
  if (target === "structure") {
    updateStructureLinks(`https://www.rcsb.org/structure/${entry.pdb}`, el("contactSourceLink").hidden ? "" : el("contactSourceLink").href);
  } else {
    updateStructureLinks(el("structureSourceLink").hidden ? "" : el("structureSourceLink").href, `https://www.rcsb.org/structure/${entry.pdb}`);
  }
}

function loadAlphaFoldViewer(target = "structure") {
  const messageId = target === "contact" ? "contactMessage" : "structureMessage";
  const viewerId = target === "contact" ? "contactViewer" : "structureViewer";
  const row = state.currentSiteRow;
  if (!row?.uniprot) {
    el(messageId).textContent = "No UniProt ID is available for this site.";
    el(viewerId).innerHTML = "";
    return;
  }
  const numericPos = row.positions?.[0] ? Number(row.positions[0]) : null;
  const selection = Number.isFinite(numericPos)
    ? [
        {
          struct_asym_id: "A",
          residue_number: numericPos,
          representationColor: { r: 255, g: 0, b: 0 },
          representation: "ball-and-stick",
          sideChain: true,
          focus: true,
        },
      ]
    : [];
  renderMolstar(
    viewerId,
    messageId,
    {
      alphafoldView: true,
      customData: { url: `https://alphafold.ebi.ac.uk/files/AF-${row.uniprot}-F1-model_v6.cif`, format: "cif" },
    },
    selection
  );
  if (target === "structure") {
    updateStructureLinks(`https://alphafold.ebi.ac.uk/entry/${row.uniprot}`, el("contactSourceLink").hidden ? "" : el("contactSourceLink").href);
  }
}

function renderStructureViewer() {
  if (selectedStructureEntry()) {
    loadPdbViewer(false, "structure");
    return;
  }
  loadAlphaFoldViewer("structure");
}

function renderContactViewer() {
  if (state.currentPdbEntries.length) {
    loadPdbViewer(true, "contact");
    return;
  }
  el("contactViewer").innerHTML = "";
  el("contactMessage").textContent = "";
  updateStructureLinks(el("structureSourceLink").hidden ? "" : el("structureSourceLink").href, "");
}

function autoLoadStructures() {
  if (!state.currentSiteRow) return;
  renderStructureViewer();
  renderContactViewer();
}

async function loadGlobalProteome() {
  if (state.globalProteome) return state.globalProteome;
  if (state.globalProteomePromise) return state.globalProteomePromise;
  const path = state.manifest?.globalProteome;
  if (!path) {
    plotMessage("globalScatter", "Run python build_global_proteome_data.py to generate the global proteome view.");
    return null;
  }
  state.globalProteomePromise = fetchJson(path)
    .then((payload) => {
      state.globalProteome = payload;
      return payload;
    })
    .catch((err) => {
      state.globalProteomePromise = null;
      plotMessage("globalScatter", `Unable to load global proteome data: ${err.message}`);
      return null;
    });
  return state.globalProteomePromise;
}

function pointSeen(point) {
  if (typeof point.seen === "boolean") return point.seen;
  return point.source === "assay" || point.source === "both";
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function pointMaxR(point, datasetKey = state.datasetKey) {
  const raw = point.maxRByDataset?.[datasetKey];
  if (raw == null) return null;
  const value = Number(raw);
  return Number.isFinite(value) ? value : null;
}

function pointHitCount(point, datasetKey = state.datasetKey) {
  const value = Number(point.hitCountByDataset?.[datasetKey]);
  return Number.isFinite(value) ? value : 0;
}

function globalCustomData(points, clusterMap) {
  return points.map((point) => {
    const cluster = clusterMap.get(String(point.cluster));
    return {
      gene: point.gene || point.uniprot,
      uniprot: point.uniprot,
      description: point.description || "",
      annotations: point.annotations || (point.annotation ? [point.annotation] : []),
      maxR: pointMaxR(point),
      hitCount: pointHitCount(point),
      seen: pointSeen(point),
      clusterId: String(point.cluster),
      clusterLabel: cluster?.label || `Cluster ${point.cluster}`,
      clusterAnnotation: cluster?.annotation || "N/A",
      clusterAnnotationPercent: cluster?.annotationPercent,
      clusterTopAnnotations: cluster?.topAnnotations || [],
    };
  });
}

function uniqueGlobalPoints(points) {
  const seen = new Set();
  return points.filter((point) => {
    const key = point.uniprot || point.gene || `${point.x}:${point.y}`;
    if (seen.has(key)) return false;
    seen.add(key);
    return true;
  });
}

function globalOverlayTrace(points, clusterMap, name, color = SPECTRAL_HIT_BLUE, glowColor = "rgba(43,131,186,0.22)") {
  return [
    {
      type: "scattergl",
      mode: "markers",
      x: points.map((point) => point.x),
      y: points.map((point) => point.y),
      marker: {
        size: 30,
        color: glowColor,
        line: { width: 0 },
      },
      hoverinfo: "skip",
      showlegend: false,
    },
    {
      type: "scattergl",
      mode: "markers",
      x: points.map((point) => point.x),
      y: points.map((point) => point.y),
      customdata: globalCustomData(points, clusterMap),
      marker: {
        size: 9,
        color,
        opacity: 0.92,
        line: { color: "black", width: 1.1 },
      },
      name,
      hoverinfo: "none",
    },
    {
      type: "scattergl",
      mode: "markers",
      x: points.map((point) => point.x),
      y: points.map((point) => point.y),
      customdata: globalCustomData(points, clusterMap),
      marker: {
        size: 46,
        color: "rgba(255,255,255,0.01)",
        opacity: 0.01,
        line: { width: 0 },
      },
      name: `${name} hover target`,
      hoverinfo: "none",
      showlegend: false,
    },
  ];
}

function formatClusterAnnotation(payload) {
  if (payload.clusterId === "-1") return null;
  const topAnnotations = (payload.clusterTopAnnotations || []).slice(0, 2);
  if (!topAnnotations.length && payload.clusterAnnotation) {
    topAnnotations.push({ annotation: payload.clusterAnnotation, percent: payload.clusterAnnotationPercent });
  }
  if (!topAnnotations.length) return null;
  return topAnnotations
    .map((item) => {
      const percent = Number(item.percent);
      const percentLabel = Number.isFinite(percent) ? `${percent.toFixed(1)}%` : "N/A";
      return `${escapeHtml(item.annotation || "N/A")} (${percentLabel})`;
    })
    .join("; ");
}

function formatProteinAnnotations(payload) {
  const clusterTerms =
    payload.clusterId === "-1"
      ? new Set()
      : new Set((payload.clusterTopAnnotations || []).slice(0, 2).map((item) => String(item.annotation || "").toLowerCase()));
  const annotations = payload.annotations || [];
  if (!annotations.length) return "N/A";
  return annotations
    .map((annotation) => {
      const escaped = escapeHtml(annotation);
      return clusterTerms.has(String(annotation).toLowerCase()) ? `<em>${escaped}</em>` : escaped;
    })
    .join("; ");
}

function showGlobalProteinTooltip(payload, event) {
  if (!payload || !event?.event) {
    hideMoleculeTooltip();
    return;
  }
  const tooltip = el("moleculeTooltip");
  const description = cleanDescription(payload.description);
  const annotations = formatProteinAnnotations(payload);
  const clusterAnnotation = formatClusterAnnotation(payload);
  tooltip.innerHTML = `
    <strong>${escapeHtml(payload.gene || payload.uniprot)}</strong>
    <div class="tooltip-uniprot">${escapeHtml(payload.uniprot || "N/A")}</div>
    <div class="tooltip-detail">${escapeHtml(description || "No description available.")}</div>
    <div class="tooltip-stats one-row">
      <div><strong>Max R</strong><span>${roundLabel(payload.maxR)}</span></div>
      <div><strong>Sites hit</strong><span>${payload.hitCount ?? 0}</span></div>
      <div><strong>Seen</strong><span>${payload.seen ? "True" : "False"}</span></div>
    </div>
    <div class="tooltip-cluster">${escapeHtml(payload.clusterLabel || "Cluster N/A")}</div>
    ${clusterAnnotation ? `<div class="tooltip-detail"><strong>Cluster annotations:</strong> ${clusterAnnotation}</div>` : ""}
    <div class="tooltip-detail"><strong>Annotations:</strong> ${annotations}</div>
  `;
  tooltip.style.display = "block";
  positionTooltip(tooltip, event.event, 360, 300);
}

function bindGlobalProteinHover() {
  const plot = el("globalScatter");
  if (typeof plot.on !== "function") return;
  plot.removeAllListeners?.("plotly_hover");
  plot.removeAllListeners?.("plotly_unhover");
  plot.on("plotly_hover", (event) => showGlobalProteinTooltip(event.points[0].customdata, event));
  plot.on("plotly_unhover", hideMoleculeTooltip);
}

function globalSquareRange(points) {
  let minX = Infinity;
  let maxX = -Infinity;
  let minY = Infinity;
  let maxY = -Infinity;
  for (const point of points) {
    if (!Number.isFinite(point.x) || !Number.isFinite(point.y)) continue;
    minX = Math.min(minX, point.x);
    maxX = Math.max(maxX, point.x);
    minY = Math.min(minY, point.y);
    maxY = Math.max(maxY, point.y);
  }
  if (![minX, maxX, minY, maxY].every(Number.isFinite)) return null;
  const centerX = (minX + maxX) / 2;
  const centerY = (minY + maxY) / 2;
  const span = Math.max(maxX - minX, maxY - minY) * 1.06;
  return {
    x: [centerX - span / 2, centerX + span / 2],
    y: [centerY - span / 2, centerY + span / 2],
  };
}

async function renderGlobalProteome() {
  if (!state.globalProteome) plotMessage("globalScatter", "Loading global proteome embedding...");
  const payload = await loadGlobalProteome();
  if (!payload) return;
  try {
    const threshold = numericValue("globalThreshold", 2);
    const showHits = checkboxChecked("globalShowHits");
    const seenOnly = checkboxChecked("globalSeenOnly");
    const showCancer = checkboxChecked("globalShowCancer");
    const showCustom = checkboxChecked("globalShowCustom");
    const customGenes = customGeneSet(el("globalCustomGenes").value);
    const points = payload.points || [];
    const fixedRange = globalSquareRange(points);
    const clusterMap = new Map((payload.clusters || []).map((cluster) => [String(cluster.id), cluster]));
    const traces = [];

    const backgroundColors = points.map((point) =>
      seenOnly && !pointSeen(point) ? PROTEOME_MUTED_GREY : clusterMap.get(String(point.cluster))?.color || MORANDI.blueLight
    );
    traces.push({
      type: "scattergl",
      mode: "markers",
      x: points.map((point) => point.x),
      y: points.map((point) => point.y),
      customdata: globalCustomData(points, clusterMap),
      marker: { color: backgroundColors, size: 5.6, opacity: 0.56, line: { width: 0 } },
      name: "Protein clusters",
      hoverinfo: "none",
    });

    const cancerHighlighted = points.filter((point) => showCancer && isCancerGene(point.gene));
    traces.push(...globalOverlayTrace(cancerHighlighted, clusterMap, "Cancer-driver genes", SPECTRAL_CANCER_RED, "rgba(215,25,28,0.2)"));

    const customHighlighted = uniqueGlobalPoints(
      points.filter((point) => showCustom && customGenes.has(String(point.gene || "").toUpperCase()))
    );
    traces.push(...globalOverlayTrace(customHighlighted, clusterMap, "Custom genes", CUSTOM_GENE_PURPLE, "rgba(123,50,148,0.22)"));

    const hitHighlighted = points
      .filter((point) => showHits && (pointMaxR(point) ?? 0) >= threshold)
      .sort((a, b) => (pointMaxR(a) ?? 0) - (pointMaxR(b) ?? 0));
    traces.push(...globalOverlayTrace(hitHighlighted, clusterMap, "Hits"));

    safePlot(
      "globalScatter",
      traces,
      baseLayout({
        margin: { l: 18, r: 18, t: 18, b: 18 },
        xaxis: { visible: false, showgrid: false, zeroline: false, constrain: "domain", range: fixedRange?.x },
        yaxis: { visible: false, showgrid: false, zeroline: false, scaleanchor: "x", scaleratio: 1, range: fixedRange?.y },
        showlegend: false,
        uirevision: "global-proteome-fixed-range",
        hoverdistance: 40,
      })
    );
    bindGlobalProteinHover();
  } catch (err) {
    plotMessage("globalScatter", `Unable to render global proteome data: ${err.message}`);
  }
}

function switchTab(tabName) {
  document.querySelectorAll(".tab-button").forEach((button) => {
    button.classList.toggle("active", button.dataset.tab === tabName);
  });
  document.querySelectorAll(".tab-panel").forEach((panel) => {
    panel.classList.toggle("active", panel.id === tabName);
  });
  updateConditionalFields();
  if (tabName === "summary") {
    scheduleSummaryDatasetRender();
  }
  if (tabName === "compound") renderCompound();
  if (tabName === "site") {
    resizeSyncedSiteBar();
    autoLoadStructures();
  }
  document.querySelectorAll(`#${tabName} .plot`).forEach((plot) => schedulePlotResize(plot));
}

window.switchTab = switchTab;

async function refreshAllForSharedFilters(sourceView) {
  copySharedFilterValues(sourceView);
  invalidateCompoundChoiceCounts();
  const refreshCompoundChoices =
    sourceView === "compound" || el("compound")?.classList.contains("active")
      ? populateCompoundChoices()
      : Promise.resolve();
  await Promise.all([refreshCompoundChoices, renderCompound(), refreshSiteControls(el("siteGeneSelect").value, { allowFallback: true })]);
}

async function refreshAllForSelectivity() {
  await refreshSiteControls(el("siteGeneSelect").value, { allowFallback: true });
}

function bindEvents() {
  initSearchableControl("compoundSelect");
  initSearchableControl("siteGeneSelect");
  document.addEventListener("pointerdown", (event) => {
    if (!state.persistentHitTooltip) return;
    const target = event.target;
    if (target?.closest?.("#siteBar .barlayer .bars .point") || target?.closest?.("#moleculeTooltip")) return;
    clearPersistentHitTooltip();
  });
  el("datasetSwitch").addEventListener("change", (event) => loadDataset(event.target.checked ? "frac" : "os"));
  el("summaryDetailClose").addEventListener("click", closeSummaryDetailModal);
  el("summaryDetailModal").addEventListener("click", (event) => {
    if (event.target === el("summaryDetailModal")) closeSummaryDetailModal();
  });
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape") closeSummaryDetailModal();
  });
  window.addEventListener("resize", () => {
    if (el("site")?.classList.contains("active")) resizeSyncedSiteBar();
  });
  el("dropNoReplicates").addEventListener("input", async () => {
    invalidateCompoundChoiceCounts();
    await populateCompoundChoices();
    await renderDatasetStaticPlots();
    await renderCompound();
    updateDatasetStatus();
  });
  document.querySelectorAll(".tab-button").forEach((button) => {
    button.addEventListener("click", () => {
      switchTab(button.dataset.tab);
    });
  });

  ["compoundSelect", "compoundLabels", "compoundColorMode", "compoundCustomGenes"].forEach((id) =>
    el(id).addEventListener("input", () => renderCompound())
  );
  el("compoundThreshold").addEventListener("input", async () => {
    invalidateCompoundChoiceCounts();
    await renderCompound();
    await updateSelectedCompoundChoiceLabel();
    scheduleCompoundChoiceListRefresh({ invalidate: false });
  });
  el("compoundActiveOnly").addEventListener("input", async () => {
    await refreshCompoundChoiceList();
  });
  ["compoundSigOnly", "compoundHideVariance", "compoundMinSn"].forEach((id) => el(id).addEventListener("input", () => refreshAllForSharedFilters("compound")));

  ["siteTargetList", "siteCustomGenes"].forEach((id) =>
    el(id).addEventListener("input", () => {
      updateConditionalFields();
      if (isContactTargetList(siteTargetListValue())) {
        setOptions(el("siteGeneSelect"), [{ value: "", label: "Loading targets..." }], "");
        setOptions(el("siteSelect"), [{ value: "", label: "Loading sites..." }], "");
      }
      refreshSiteControls(null, { allowFallback: true });
    })
  );
  el("siteGeneSelect").addEventListener("change", (event) => refreshSiteControls(event.target.value));
  el("siteMaxHits").addEventListener("input", () => refreshAllForSelectivity("siteMaxHits"));
  ["siteSigOnly", "siteHideVariance", "siteMinSn"].forEach((id) => el(id).addEventListener("input", () => refreshAllForSharedFilters("site")));
  el("siteActiveOnly").addEventListener("input", () => refreshSiteControls(el("siteGeneSelect").value, { allowFallback: true }));
  el("siteSelect").addEventListener("input", () => renderSite());
  el("structureSelect").addEventListener("input", () => {
    if (el("site").classList.contains("active")) renderStructureViewer();
  });
  el("contactPdbSelect").addEventListener("input", () => {
    if (el("site").classList.contains("active")) renderContactViewer();
  });

  ["globalThreshold", "globalCustomGenes"].forEach((id) => el(id).addEventListener("input", () => renderGlobalProteome()));
  ["globalShowHits", "globalSeenOnly", "globalShowCancer", "globalShowCustom"].forEach((id) => el(id).addEventListener("input", () => {
    updateConditionalFields();
    renderGlobalProteome();
  }));
  el("circularCompoundAtlas").addEventListener("mousemove", (event) => showCircularTooltip(event, circularHitAtEvent(event)));
  el("circularCompoundAtlas").addEventListener("mouseleave", () => {
    const tooltip = el("circularAtlasTooltip");
    if (tooltip) tooltip.hidden = true;
  });
  el("circularCompoundAtlas").addEventListener("click", (event) => {
    const hit = circularHitAtEvent(event);
    if (hit?.row?.id) void renderCircularDrilldown(hit.row.id, hit.segment);
  });
  el("circularCompoundAtlas").addEventListener("wheel", (event) => {
    event.preventDefault();
    const zoom = state.circularZoom;
    const factor = event.deltaY < 0 ? 1.16 : 0.86;
    zoom.scale = Math.max(1, Math.min(10, zoom.scale * factor));
    if (zoom.scale === 1) {
      zoom.offsetX = 0;
      zoom.offsetY = 0;
    }
    renderCircularCompoundAtlas();
  }, { passive: false });
  el("circularCompoundAtlas").addEventListener("mousedown", (event) => {
    const zoom = state.circularZoom;
    zoom.dragging = true;
    zoom.lastX = event.clientX;
    zoom.lastY = event.clientY;
  });
  window.addEventListener("mousemove", (event) => {
    const zoom = state.circularZoom;
    if (!zoom.dragging) return;
    zoom.offsetX += event.clientX - zoom.lastX;
    zoom.offsetY += event.clientY - zoom.lastY;
    zoom.lastX = event.clientX;
    zoom.lastY = event.clientY;
    renderCircularCompoundAtlas();
  });
  window.addEventListener("mouseup", () => {
    state.circularZoom.dragging = false;
  });
  el("circularCompoundAtlas").addEventListener("dblclick", () => {
    state.circularZoom = { scale: 1, offsetX: 0, offsetY: 0, dragging: false, lastX: 0, lastY: 0 };
    renderCircularCompoundAtlas();
  });
  window.addEventListener("resize", () => {
    if (el("summary")?.classList.contains("active")) renderCircularCompoundAtlas();
  });
}

async function init() {
  try {
    state.manifest = await fetchJson("manifest.json");
    state.cancerSet = new Set((state.manifest.cancerGenes || []).map((gene) => gene.toUpperCase()));
    bindEvents();
    await loadDataset("os");
  } catch (err) {
    setStatus(`${err.message}. Run "python build_static_data.py" and serve with "python -m http.server".`);
  }
}

init();

