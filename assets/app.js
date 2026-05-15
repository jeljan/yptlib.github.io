const DATA_ROOTS = window.YPTLIB_DATA_ROOT
  ? [window.YPTLIB_DATA_ROOT]
  : [
      "assets/data/",
      "https://cdn.jsdelivr.net/gh/jeljan/yptlib@a7e92171d1e6127153f516a19d37a9dceddb5cad/assets/data/",
      "https://cdn.statically.io/gh/jeljan/yptlib/a7e92171d1e6127153f516a19d37a9dceddb5cad/assets/data/",
      "https://raw.githubusercontent.com/jeljan/yptlib/a7e92171d1e6127153f516a19d37a9dceddb5cad/assets/data/",
    ];
const RELEASE_DATA_ARCHIVE_URL = window.YPTLIB_DATA_ARCHIVE_URL || "";
const RELEASE_DATA_ARCHIVE_ENABLED = Boolean(RELEASE_DATA_ARCHIVE_URL);
let releaseArchivePromise = null;
const ASSET_VERSION = "pages-data-v12";
const FETCH_RETRIES_PER_ROOT = 2;
const FETCH_TIMEOUT_MS = 12000;
const SUMMARY_GENE_FETCH_CONCURRENCY = 8;
const CONTACT_TYPES = ["PPI", "Dna", "Rna", "Metal", "Ligand", "Cofactor"];
const CONTACT_LABELS = {
  PPI: "PPI",
  Dna: "DNA",
  Rna: "RNA",
  Metal: "Metal",
  Ligand: "Ligand",
  Cofactor: "Cofactor",
};
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
  activeGeneFilterIndexCache: new Map(),
  contactCache: new Map(),
  filteredSitesByGene: new Map(),
  compoundChoiceCountCache: new Map(),
  compoundChoiceRefreshToken: 0,
  compoundChoiceRefreshTimer: null,
  compoundWarmupToken: 0,
  currentSiteRow: null,
  currentSiteSummary: null,
  currentPdbEntries: [],
  activeSummaryDrug: null,
  activeSiteDrug: null,
  siteRefreshSeq: 0,
  globalProteome: null,
  globalProteomePromise: null,
};

const el = (id) => document.getElementById(id);

function setStatus(text) {
  el("statusText").textContent = text;
}

function activeDatasetDrugs() {
  const noReplicateDrugs = new Set(state.dataset?.noReplicateDrugs || []);
  return (state.dataset?.rawDrugs || []).filter((drug) => !checkboxChecked("dropNoReplicates") || !noReplicateDrugs.has(drug));
}

function updateDatasetStatus() {
  if (!state.dataset || !state.catalog) return;
  const total = state.dataset.rawDrugs.length;
  const active = activeDatasetDrugs().length;
  const suffix = checkboxChecked("dropNoReplicates")
    ? `${active.toLocaleString()} compounds after no-replicate filter (${(total - active).toLocaleString()} dropped)`
    : `${total.toLocaleString()} compounds`;
  setStatus(`${state.dataset.label}: ${state.catalog.length.toLocaleString()} sites, ${suffix}`);
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
    select.appendChild(node);
  }
  if (selected != null) {
    select.value = selected;
  }
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
      FS: "Sulfonyl Fluorides",
      SuFEx: "Fluorosulfates",
      SuTEx: "Sulfonyl Triazoles",
      "Sulfonyl Fluorides": "Sulfonyl Fluorides",
      Fluorosulfates: "Fluorosulfates",
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
  if (num === 0) return "0.00000";
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

function summaryFilterIds() {
  return {
    sigOnly: "summarySigOnly",
    hideVariance: "summaryHideVariance",
    minSn: "summaryMinSn",
  };
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
  for (const view of ["summary", "compound", "site"]) {
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
    font: { family: "system-ui, sans-serif", size: 12 },
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
  const summaryRow = Number(el("summarySiteSelect")?.value);
  return Number.isFinite(summaryRow) ? summaryRow : null;
}

function rawHoverChartHtml(values) {
  if (!values || values.length < 4) return `<div class="raw-hover-empty muted">Raw intensities unavailable.</div>`;
  const [d1, d2, c1, c2] = values.map((v) => Number(v));
  const groups = [
    { label: "DMSO", vals: [d1, d2], color: MORANDI.gray },
    { label: "Compound", vals: [c1, c2], color: MORANDI.blueDark },
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
    .map((site) => `<li><strong>${escapeHtml(site.label)}</strong><span>${roundLabel(site.value)}%</span></li>`)
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
    node.onmouseenter = (event) => showHitTooltip(payload, { event });
    node.onmousemove = (event) => showHitTooltip(payload, { event });
    node.onclick = (event) => showHitTooltip(payload, { event });
    node.onmouseleave = hideMoleculeTooltip;
  });
  const dragLayer = plot.querySelector(".nsewdrag");
  if (dragLayer) {
    dragLayer.onmouseleave = hideMoleculeTooltip;
  }
}

async function loadDataset(datasetKey) {
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
    state.filteredSitesByGene.clear();
    state.rawHoverAliasCache.clear();
    state.currentSiteRow = null;
    state.currentSiteSummary = null;

    await populateDatasetControls();
    renderDatasetStaticPlots();
    await renderSummary();
    await renderCompound();
    scheduleCompoundCacheWarmup();
    await refreshSiteControls();
    if (el("summary")?.classList.contains("active")) await renderGlobalProteome();
    renderBioTable();
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
  setOptions(el("siteGeneSelect"), geneOptions, state.dataset.defaultGene);

  const hitCounts = [...new Set(Object.values(state.dataset.drugHitCounts || {}).map((v) => Number(v)).filter((v) => Number.isFinite(v) && v > 0))].sort(
    (a, b) => a - b
  );
  setOptions(
    el("summaryMaxHits"),
    [{ value: "all", label: "Any" }, ...hitCounts.map((v) => ({ value: String(v), label: `≤ ${v} site${v === 1 ? "" : "s"} hit` }))],
    "all"
  );
  setOptions(
    el("siteMaxHits"),
    [{ value: "all", label: "Any" }, ...hitCounts.map((v) => ({ value: String(v), label: `≤ ${v} site${v === 1 ? "" : "s"} hit` }))],
    el("summaryMaxHits").value || "all"
  );

  setOptions(el("summaryTargetList"), [
    { value: "cancer", label: "Cancer-driver genes" },
    { value: "PPI", label: "Sites at PPI interface" },
    { value: "Metal", label: "Sites near metal" },
    { value: "Cofactor", label: "Sites near cofactor" },
    { value: "Dna", label: "Sites near DNA" },
    { value: "Rna", label: "Sites near RNA" },
    { value: "Ligand", label: "Sites near ligand" },
    { value: "custom", label: "Custom gene list" },
  ]);

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
  const drugs = [...(dataset.rawDrugs || [])];
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
  const current = el("compoundSelect").value || state.dataset.defaultDrug;
  let dynamicCounts = null;
  if (compoundQualityFiltersActive() || numericValue("compoundThreshold", 2) !== 2) {
    dynamicCounts = await dynamicCompoundHitCounts(state.dataset.rawDrugs, refreshToken);
    if (!dynamicCounts) return false;
  }
  const shownHitCount = (drug) => dynamicCounts?.get(drug) ?? drugHitCount(drug);
  const drugOptions = state.dataset.rawDrugs
    .filter((drug) => !dropNoReplicates || !noReplicateDrugs.has(drug))
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
  if (state.compoundChoiceRefreshTimer) clearTimeout(state.compoundChoiceRefreshTimer);
  state.compoundChoiceRefreshTimer = setTimeout(() => {
    state.compoundChoiceRefreshTimer = null;
    refreshCompoundChoiceList({ refreshToken, invalidate });
  }, 0);
}

async function refreshCompoundChoiceList({ refreshToken = null, invalidate = true } = {}) {
  const activeRefreshToken = refreshToken ?? ++state.compoundChoiceRefreshToken;
  if (invalidate) invalidateCompoundChoiceCounts();
  try {
    const updated = await populateCompoundChoices(activeRefreshToken);
    if (updated) await renderCompound();
  } catch (err) {
    console.warn("Unable to refresh compound choice list", err);
  }
}

function updateConditionalFields() {
  el("summaryCustomGenesWrap").hidden = el("summaryTargetList").value !== "custom";
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
    siteBins[binIndex].sites.push({
      label: row?.label || `Site ${idx + 1}`,
      gene: row?.gene || "",
      value,
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

function rowMatchesTarget(row, targetList, customGenes) {
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

async function renderSummary() {
  updateConditionalFields();
  const targetList = el("summaryTargetList").value;
  const customGenes = customGeneSet(el("summaryCustomGenes").value);
  const maxHits = hitCountLimit("summaryMaxHits");
  const candidates = state.catalog.filter((row) => rowMatchesTarget(row, targetList, customGenes));
  const uniqueGenes = [...new Set(candidates.map((row) => row.gene))];
  await mapWithConcurrency(uniqueGenes, SUMMARY_GENE_FETCH_CONCURRENCY, (gene) => loadGeneSites(gene));
  const scoredRows = candidates
    .map((row) => {
      const payload = state.geneSiteCache.get(row.gene);
      const summary = payload?.sites.find((site) => site.row === row.i || site.label === row.label);
      const hits = filteredHits(summary, summaryFilterIds(), maxHits);
      const filteredMaxR = hits[0]?.[1] ?? null;
      return { ...row, filteredMaxR };
    })
    .filter((row) => Number.isFinite(Number(row.filteredMaxR)));
  const dedupedRows = new Map();
  for (const row of scoredRows) {
    const key = row.label;
    const prev = dedupedRows.get(key);
    if (!prev || row.filteredMaxR > prev.filteredMaxR || (row.filteredMaxR === prev.filteredMaxR && row.i < prev.i)) {
      dedupedRows.set(key, row);
    }
  }
  const rows = [...dedupedRows.values()]
    .sort((a, b) => (b.filteredMaxR ?? -1) - (a.filteredMaxR ?? -1) || a.label.localeCompare(b.label))
    .slice(0, 2000);
  const current = el("summarySiteSelect").value;
  const selected = rows.some((row) => String(row.i) === current)
    ? current
    : rows[0]
      ? String(rows[0].i)
      : "";

  setOptions(
    el("summarySiteSelect"),
    rows.length
      ? rows.map((row) => ({ value: String(row.i), label: `${row.label} (Max R: ${roundLabel(row.filteredMaxR)})` }))
      : [{ value: "", label: "No sites found" }],
    selected
  );
  await renderSummaryBar();
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
    activeOnly ? 1 : 0,
    el("siteMaxHits")?.value || "all",
    checkboxChecked("siteSigOnly") ? 1 : 0,
    checkboxChecked("siteHideVariance") ? 1 : 0,
    numericValue("siteMinSn", 0),
  ].join("|");
}

function renderHitBar(targetId, hits, threshold, rowId = null) {
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
      colorbar: { title: { text: "P-value", side: "right" }, len: 0.82 },
    },
    customdata: hoverPayload,
    hoverinfo: "none",
  };
  const plotted = safePlot(
    targetId,
    [trace],
    mergeAxes(baseLayout({
      shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: threshold, y1: threshold, line: { dash: "dash", color: "#666" } }],
    }), { title: `Top ${topHits.length} Compounds (/${hits.length})`, tickangle: 42, automargin: true }, { title: "R" })
  );
  if (!plotted) return;
  const plot = el(targetId);
  if (typeof plot.on === "function") {
    plot.removeAllListeners?.("plotly_hover");
    plot.removeAllListeners?.("plotly_unhover");
    plot.removeAllListeners?.("plotly_click");
    plot.on("plotly_hover", (event) => showHitTooltip(event.points[0].customdata, event));
    plot.on("plotly_unhover", hideMoleculeTooltip);
    plot.on("plotly_click", (event) => showHitTooltip(event.points[0].customdata, event));
  }
  requestAnimationFrame(() => bindBarHover(targetId, hoverPayload));
}

async function renderSummaryBar() {
  const row = state.rowById.get(Number(el("summarySiteSelect").value));
  const summary = await siteSummaryForRow(row);
  const maxProm = hitCountLimit("summaryMaxHits");
  const hits = filteredHits(summary, {
    sigOnly: "summarySigOnly",
    hideVariance: "summaryHideVariance",
    minSn: "summaryMinSn",
  }, maxProm);
  renderHitBar("summaryBar", hits, 2.0, row?.i ?? null);
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
  const payload = await loadGeneSites(gene);
  return (payload?.sites || [])
    .map((site) => {
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
  const currentGene = preferredGene || el("siteGeneSelect").value || state.dataset.defaultGene;
  const currentSite = el("siteSelect").value;
  let geneStats = activeOnly ? await activeGeneStatsForCurrentFilters() : null;
  if (seq !== state.siteRefreshSeq) return;

  if (!geneStats) {
    geneStats = new Map();
    for (const row of state.catalog) {
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

function renderBioTable() {
  const filter = String(el("bioGeneFilter").value || "").trim().toUpperCase();
  const rows = state.catalog
    .filter((row) => !filter || String(row.gene).toUpperCase().includes(filter))
    .sort((a, b) => (b.maxR ?? -1) - (a.maxR ?? -1))
    .slice(0, 250);
  el("bioTable").querySelector("tbody").innerHTML = rows
    .map(
      (row) =>
        `<tr><td>${row.gene}</td><td>${row.site}</td><td>${roundLabel(row.maxR)}</td><td>Pending data</td><td>${row.contactTypes.length ? row.contactTypes.join(", ") : "No mapped contact"}</td></tr>`
    )
    .join("");
}

function switchTab(tabName) {
  document.querySelectorAll(".tab-button").forEach((button) => {
    button.classList.toggle("active", button.dataset.tab === tabName);
  });
  document.querySelectorAll(".tab-panel").forEach((panel) => {
    panel.classList.toggle("active", panel.id === tabName);
  });
  updateConditionalFields();
  if (tabName === "summary") renderGlobalProteome();
  if (tabName === "bio") renderBioTable();
  if (tabName === "site") autoLoadStructures();
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
  await Promise.all([refreshCompoundChoices, renderSummary(), renderCompound(), refreshSiteControls(el("siteGeneSelect").value, { allowFallback: true })]);
}

async function refreshAllForSelectivity(sourceId) {
  const targetId = sourceId === "summaryMaxHits" ? "siteMaxHits" : "summaryMaxHits";
  copySelectivityValue(sourceId, targetId);
  await Promise.all([renderSummary(), refreshSiteControls(el("siteGeneSelect").value, { allowFallback: true })]);
}

function bindEvents() {
  el("datasetSwitch").addEventListener("change", (event) => loadDataset(event.target.checked ? "frac" : "os"));
  el("summaryDetailClose").addEventListener("click", closeSummaryDetailModal);
  el("summaryDetailModal").addEventListener("click", (event) => {
    if (event.target === el("summaryDetailModal")) closeSummaryDetailModal();
  });
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape") closeSummaryDetailModal();
  });
  el("dropNoReplicates").addEventListener("input", async () => {
    invalidateCompoundChoiceCounts();
    await populateCompoundChoices();
    await renderDatasetStaticPlots();
    await renderSummary();
    await renderCompound();
    updateDatasetStatus();
  });
  document.querySelectorAll(".tab-button").forEach((button) => {
    button.addEventListener("click", () => {
      switchTab(button.dataset.tab);
    });
  });

  ["summaryTargetList", "summaryCustomGenes"].forEach((id) => el(id).addEventListener("input", () => renderSummary()));
  el("summaryMaxHits").addEventListener("input", () => refreshAllForSelectivity("summaryMaxHits"));
  ["summarySigOnly", "summaryHideVariance", "summaryMinSn"].forEach((id) => el(id).addEventListener("input", () => refreshAllForSharedFilters("summary")));
  el("summarySiteSelect").addEventListener("input", () => renderSummaryBar());

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
  el("bioGeneFilter").addEventListener("input", () => renderBioTable());
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
