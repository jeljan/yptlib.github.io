const ASSET_VERSION = "pages-data-v26";
const BIO_FAMILIES = ["catalytic", "regulatory", "variant", "structure", "constraint"];

const state = {
  manifest: null,
  manifestBase: null,
  dictionaries: {},
  rowPartCache: new Map(),
  indexCache: new Map(),
  proteinIndexCache: new Map(),
  proteinDomainCache: new Map(),
  domainIndexCache: new Map(),
  positiveSetCache: new Map(),
  resultCache: new Map(),
  domainEntries: null,
};

async function fetchJson(url) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Unable to fetch ${url}: ${response.status}`);
  return response.json();
}

function atlasUrl(path) {
  return new URL(path, state.manifestBase).href;
}

async function loadAtlas(manifestPath) {
  if (state.manifest) return;
  const manifestUrl = new URL(`data/${manifestPath || "atlas/manifest.json"}?v=${ASSET_VERSION}`, self.location.href);
  state.manifest = await fetchJson(manifestUrl.href);
  state.manifestBase = new URL(".", manifestUrl.href);
  await Promise.all(
    Object.entries(state.manifest.dictionaries || {}).map(async ([name, rel]) => {
      state.dictionaries[name] = await fetchJson(atlasUrl(`${rel}?v=${ASSET_VERSION}`));
    })
  );
}

async function loadRowPart(partIndex) {
  if (state.rowPartCache.has(partIndex)) return state.rowPartCache.get(partIndex);
  const rel = state.manifest.rowParts[partIndex];
  const payload = await fetchJson(atlasUrl(`${rel}?v=${ASSET_VERSION}`));
  state.rowPartCache.set(partIndex, payload.rows);
  return payload.rows;
}

async function rawRow(rowId) {
  const partIndex = Math.floor(rowId / state.manifest.partSize);
  const rows = await loadRowPart(partIndex);
  return rows[rowId - partIndex * state.manifest.partSize];
}

function sortKey(column, datasetKey) {
  if (column === "seen") return `seen_${datasetKey}`;
  if (column === "maxR") return `maxR_${datasetKey}`;
  if (column === "hits") return `hits_${datasetKey}`;
  return column || "score";
}

async function loadIndexParts(indexSpec) {
  const cacheKey = indexSpec.parts.join("|");
  if (state.indexCache.has(cacheKey)) return state.indexCache.get(cacheKey);
  const chunks = await Promise.all(indexSpec.parts.map((rel) => fetchJson(atlasUrl(`${rel}?v=${ASSET_VERSION}`))));
  const ids = chunks.flatMap((chunk) => chunk.rowIds || []);
  state.indexCache.set(cacheKey, ids);
  return ids;
}

async function orderedSections(sort, datasetKey) {
  const key = sortKey(sort.column, datasetKey);
  const direction = sort.direction || "desc";
  const spec = state.manifest.sortIndexes[key]?.[direction] || state.manifest.sortIndexes.score.desc;
  const valid = await loadIndexParts(spec);
  if (spec.remainder !== "default") return [{ ids: valid, skip: null }];
  return [
    { ids: valid, skip: null },
    { ids: await loadIndexParts(state.manifest.defaultIndex), skip: new Set(valid) },
  ];
}

async function positiveSet(key) {
  if (state.positiveSetCache.has(key)) return state.positiveSetCache.get(key);
  const spec = state.manifest.sortIndexes[key]?.desc;
  if (!spec) return new Set();
  const ids = await loadIndexParts(spec);
  const set = new Set(ids);
  state.positiveSetCache.set(key, set);
  return set;
}

async function evidencePositiveSet(family) {
  if (!family || family === "all") return null;
  if (family === "catalytic") return positiveSet("catalytic");
  if (family === "regulatory") return positiveSet("regulatory");
  if (family === "variant") return positiveSet("variant");
  if (family === "structure") return positiveSet("structure");
  if (family === "constraint") return positiveSet("constraint");
  if (family === "any") {
    const sets = await Promise.all(["catalytic", "regulatory", "variant", "structure", "constraint"].map((key) => positiveSet(key)));
    return new Set(sets.flatMap((set) => [...set]));
  }
  return null;
}

function intersectCandidateSets(sets) {
  const active = sets.filter(Boolean);
  if (!active.length) return null;
  active.sort((a, b) => a.size - b.size);
  const result = new Set();
  for (const rowId of active[0]) {
    if (active.every((set) => set.has(rowId))) result.add(rowId);
  }
  return result;
}

function dict(name, id) {
  return state.dictionaries[name]?.[id] ?? "";
}

function scale(value, scaleName) {
  const scaleValue = state.manifest.scales?.[scaleName] || 1;
  return Number(value || 0) / scaleValue;
}

function decodeRow(rowId, row) {
  const i = state.manifest.rowIndex;
  return {
    rowId,
    gene: dict("genes", row[i.gene]),
    site: row[i.site],
    uniprot: dict("uniprots", row[i.uniprot]),
    score: scale(row[i.score], "score"),
    tier: dict("tiers", row[i.tier]),
    catalytic: dict("catalytic", row[i.catalytic]),
    regulatory: dict("regulatory", row[i.regulatory]),
    variant: row[i.variant],
    structure: row[i.structure],
    constraint: dict("constraint", row[i.constraint]),
    amMax: scale(row[i.amMax], "float"),
    gerp: scale(row[i.gerp], "float"),
    phyloP: scale(row[i.phyloP], "float"),
    phastCons: scale(row[i.phastCons], "float"),
    plddt: scale(row[i.plddt], "float"),
    rsa: scale(row[i.rsa], "float"),
    pdbCount: row[i.pdbCount],
    essentiality: dict("essentiality", row[i.essentiality]),
    otMax: scale(row[i.otMax], "float"),
    systemsDegree: row[i.systemsDegree],
    expressionMax: scale(row[i.expressionMax], "float"),
    activeSiteDistance: scale(row[i.activeSiteDistance], "float"),
    gnomadMaxAf: scale(row[i.gnomadMaxAf], "score"),
    cosmicSamples: row[i.cosmicSamples],
    clinvarPathogenic: row[i.clinvarPathogenic],
    mc3Variants: row[i.mc3Variants],
    interproCount: row[i.interproCount],
    epsdCount: row[i.epsdCount],
    intrinsicRegion: dict("intrinsicRegion", row[i.intrinsicRegion]),
    prolineContext: row[i.prolineContext],
    modResDistance: scale(row[i.modResDistance], "float"),
    bindingDistance: scale(row[i.bindingDistance], "float"),
    siteDistance: scale(row[i.siteDistance], "float"),
    evidenceBasis: dict("evidenceBasis", row[i.evidenceBasis]),
    tierConfidence: dict("tierConfidence", row[i.tierConfidence]),
    dataset: {
      os: {
        seen: Boolean(row[i.osSeen]),
        maxR: scale(row[i.osMaxR], "ratio"),
        hits: row[i.osHits],
        observed: row[i.osObserved],
        siteRow: row[i.osSiteRow],
        qualitySummary: row[i.osQuality] || [],
      },
      frac: {
        seen: Boolean(row[i.fracSeen]),
        maxR: scale(row[i.fracMaxR], "ratio"),
        hits: row[i.fracHits],
        observed: row[i.fracObserved],
        siteRow: row[i.fracSiteRow],
        qualitySummary: row[i.fracQuality] || [],
      },
    },
  };
}

function dataFiltersActive(filters) {
  return Boolean(
    filters.sigOnly ||
      filters.hideVariance ||
      Number(filters.minSn || 0) > 0 ||
      (filters.maxHits && filters.maxHits !== "all")
  );
}

function qualitySummary(row, filters) {
  const datasetKey = filters.datasetKey || "os";
  const dataset = row.dataset?.[datasetKey] || {};
  const buckets = dataset.qualitySummary || [];
  const maxHits = filters.maxHits && filters.maxHits !== "all" ? Number(filters.maxHits) : Infinity;
  const minSn = Number(filters.minSn || 0);
  let observed = 0;
  let hits = 0;
  let maxR = 0;
  for (const bucket of buckets) {
    const compoundHitCount = Number(bucket[0] || 0);
    const meanSnFloor = Number(bucket[1] || 0);
    const flags = Number(bucket[2] || 0);
    if (compoundHitCount > maxHits) continue;
    if (meanSnFloor < minSn) continue;
    if (filters.hideVariance && !(flags & 1)) continue;
    if (filters.sigOnly && !(flags & 2)) continue;
    observed += Number(bucket[4] || 0);
    hits += Number(bucket[5] || 0);
    maxR = Math.max(maxR, scale(bucket[3], "ratio"));
  }
  return {
    ...dataset,
    seen: observed > 0,
    maxR,
    hits,
    observed,
  };
}

function filteredDatasetForRow(row, filters) {
  const dataset = row.dataset?.[filters.datasetKey || "os"] || {};
  return dataFiltersActive(filters) ? qualitySummary(row, filters) : dataset;
}

function applyDataFilters(row, filters) {
  if (!dataFiltersActive(filters)) return row;
  const datasetKey = filters.datasetKey || "os";
  return {
    ...row,
    dataset: {
      ...row.dataset,
      [datasetKey]: filteredDatasetForRow(row, filters),
    },
  };
}

function evidenceMatches(row, family) {
  if (!family || family === "all") return true;
  if (family === "any") return row.catalytic !== "NONE" || row.regulatory !== "NONE" || row.variant > 0 || row.structure > 0 || row.constraint !== "NONE";
  if (family === "catalytic") return row.catalytic !== "NONE";
  if (family === "regulatory") return row.regulatory !== "NONE";
  if (family === "variant") return row.variant > 0;
  if (family === "structure") return row.structure > 0;
  if (family === "constraint") return row.constraint !== "NONE";
  return true;
}

function familyPositive(row, family) {
  if (family === "catalytic") return row.catalytic !== "NONE";
  if (family === "regulatory") return row.regulatory !== "NONE";
  if (family === "variant") return row.variant > 0;
  if (family === "structure") return row.structure > 0;
  if (family === "constraint") return row.constraint !== "NONE";
  return false;
}

function dominantFamily(row) {
  return BIO_FAMILIES.find((family) => familyPositive(row, family)) || "none";
}

function contextMatches(row, context) {
  if (!context || context === "all") return true;
  if (context === "interpro") return Number(row.interproCount || 0) > 0;
  if (context === "epsd") return Number(row.epsdCount || 0) > 0;
  if (context === "mc3") return Number(row.mc3Variants || 0) > 0;
  if (context === "modResProximity") return Number(row.modResDistance || 0) > 0;
  if (context === "bindingProximity") return Number(row.bindingDistance || 0) > 0;
  if (context === "siteProximity") return Number(row.siteDistance || 0) > 0;
  if (context === "prolineContext") return Number(row.prolineContext || 0) > 0;
  return true;
}

function rowMatchesDecoded(row, filters) {
  const dataset = filteredDatasetForRow(row, filters);
  if (filters.targetList === "cancer" && !(filters.cancerGenes || []).includes(row.gene.toUpperCase())) return false;
  if (filters.targetList === "custom" && !(filters.customGenes || []).includes(row.gene.toUpperCase())) return false;
  if (filters.tier && filters.tier !== "all" && row.tier !== filters.tier) return false;
  if (Number(filters.minScore || 0) > 0 && row.score < Number(filters.minScore || 0)) return false;
  if (!evidenceMatches(row, filters.evidenceFamily)) return false;
  if (!contextMatches(row, filters.contextFilter)) return false;
  const search = String(filters.search || "").trim().toUpperCase();
  if (search) {
    const label = `${row.gene} ${row.uniprot} Y${row.site} ${row.site}`.toUpperCase();
    if (!label.includes(search)) return false;
  }
  if (filters.seenOnly && !dataset.seen) return false;
  if (dataFiltersActive(filters) && !dataset.seen) return false;
  if (filters.hitOnly && !(Number(dataset.hits || 0) > 0)) return false;
  if (filters.bin) {
    const maxR = Number(dataset.maxR || 0);
    const score = Number(row.score || 0);
    if (maxR < Number(filters.bin.maxRMin || 0) || maxR >= Number(filters.bin.maxRMax || Infinity)) return false;
    if (score < Number(filters.bin.scoreMin || 0) || score >= Number(filters.bin.scoreMax || Infinity)) return false;
  }
  return true;
}

function defaultCompare(a, b) {
  return (b.score - a.score) || a.gene.localeCompare(b.gene) || (a.site - b.site) || a.uniprot.localeCompare(b.uniprot);
}

function sortValue(row, column, datasetKey, filters = { datasetKey }) {
  if (column === "gene") return row.gene;
  if (column === "site") return row.site;
  if (column === "uniprot") return row.uniprot;
  const dataset = filteredDatasetForRow(row, filters);
  if (column === "seen") return dataset.seen ? 1 : 0;
  if (column === "maxR") return dataset.maxR || 0;
  if (column === "hits") return dataset.hits || 0;
  return row[column || "score"];
}

function sortValueIsValid(row, column, datasetKey, filters = { datasetKey }) {
  const value = sortValue(row, column, datasetKey, filters);
  if (column === "gene" || column === "site" || column === "uniprot") return true;
  if (column === "tier") return value && value !== "LOW";
  if (column === "catalytic" || column === "regulatory" || column === "constraint") return value && value !== "NONE";
  if (typeof value === "boolean") return value;
  if (typeof value === "number") return value !== 0;
  return Boolean(value);
}

function compareRowsForSort(sort, datasetKey, filters = { datasetKey }) {
  const column = sort.column || "score";
  const direction = sort.direction || "desc";
  const factor = direction === "asc" ? 1 : -1;
  return (a, b) => {
    const aValid = sortValueIsValid(a, column, datasetKey, filters);
    const bValid = sortValueIsValid(b, column, datasetKey, filters);
    if (aValid !== bValid) return aValid ? -1 : 1;
    if (!aValid && !bValid) return defaultCompare(a, b);
    const aValue = sortValue(a, column, datasetKey, filters);
    const bValue = sortValue(b, column, datasetKey, filters);
    if (typeof aValue === "string" || typeof bValue === "string") {
      const textCompare = String(aValue).localeCompare(String(bValue));
      return textCompare ? textCompare * factor : defaultCompare(a, b);
    }
    const numericCompare = (Number(aValue) - Number(bValue)) * factor;
    return numericCompare || defaultCompare(a, b);
  };
}

function needsDecodedFilter(filters) {
  return Boolean(
    String(filters.search || "").trim() ||
      (filters.targetList && filters.targetList !== "all") ||
      (filters.tier && filters.tier !== "all") ||
      (filters.contextFilter && filters.contextFilter !== "all") ||
      Number(filters.minScore || 0) > 0 ||
      Boolean(filters.bin) ||
      dataFiltersActive(filters)
  );
}

async function candidateSets(filters) {
  const datasetKey = filters.datasetKey || "os";
  return {
    seenSet: filters.seenOnly || dataFiltersActive(filters) ? await positiveSet(`seen_${datasetKey}`) : null,
    hitSet: filters.hitOnly && !dataFiltersActive(filters) ? await positiveSet(`hits_${datasetKey}`) : null,
    evidenceSet: await evidencePositiveSet(filters.evidenceFamily),
  };
}

async function query(message) {
  await loadAtlas(message.manifestPath);
  const filters = message.filters || {};
  const pageSize = state.manifest.pageSize || 100;
  const page = Math.max(0, Number(message.page || 0));
  const cacheKey = JSON.stringify({ filters, sort: message.sort || {}, page });
  if (state.resultCache.has(cacheKey)) return state.resultCache.get(cacheKey);

  const sections = await orderedSections(message.sort || {}, filters.datasetKey || "os");
  const { seenSet, hitSet, evidenceSet } = await candidateSets(filters);
  const decodeForFilter = needsDecodedFilter(filters);
  const start = page * pageSize;
  const end = start + pageSize;
  let total = 0;
  const rows = [];

  const limitedSet = intersectCandidateSets([seenSet, hitSet, evidenceSet]);
  const hasDatasetLimitedSet = Boolean(seenSet || hitSet);
  if (hasDatasetLimitedSet && limitedSet && limitedSet.size <= 400_000) {
    const decodedRows = [];
    for (const rowId of limitedSet) {
      if (seenSet && !seenSet.has(rowId)) continue;
      if (hitSet && !hitSet.has(rowId)) continue;
      const decoded = applyDataFilters(decodeRow(rowId, await rawRow(rowId)), filters);
      if (!rowMatchesDecoded(decoded, filters)) continue;
      decodedRows.push(decoded);
    }
      decodedRows.sort(compareRowsForSort(message.sort || {}, filters.datasetKey || "os", filters));
    const result = {
      type: "result",
      seq: message.seq,
      rows: decodedRows.slice(start, end),
      total: decodedRows.length,
      page,
      pageSize,
    };
    state.resultCache.set(cacheKey, result);
    return result;
  }

  const knownTotal = !decodeForFilter && limitedSet ? limitedSet.size : null;
  outer:
  for (const section of sections) {
    for (const rowId of section.ids) {
      if (section.skip?.has(rowId)) continue;
      if (seenSet && !seenSet.has(rowId)) continue;
      if (hitSet && !hitSet.has(rowId)) continue;
      if (evidenceSet && !evidenceSet.has(rowId)) continue;
      let decoded = null;
      if (decodeForFilter) {
        decoded = applyDataFilters(decodeRow(rowId, await rawRow(rowId)), filters);
        if (!rowMatchesDecoded(decoded, filters)) continue;
      }
      if (total >= start && total < end) {
        rows.push(decoded || applyDataFilters(decodeRow(rowId, await rawRow(rowId)), filters));
      }
      total++;
      if (knownTotal != null && total >= end) break outer;
    }
  }
  const result = { type: "result", seq: message.seq, rows, total: knownTotal ?? total, page, pageSize };
  state.resultCache.set(cacheKey, result);
  return result;
}

function blankFamilyCounts() {
  return Object.fromEntries(BIO_FAMILIES.map((family) => [family, 0]));
}

function addFamilyCounts(target, row) {
  for (const family of BIO_FAMILIES) {
    if (familyPositive(row, family)) target[family] += 1;
  }
}

function scoreBin(score) {
  const edges = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.000001];
  for (let index = 0; index < edges.length - 1; index += 1) {
    if (score >= edges[index] && score < edges[index + 1]) return { index, min: edges[index], max: edges[index + 1] };
  }
  return { index: edges.length - 2, min: 0.8, max: 1.000001 };
}

function maxRBin(maxR) {
  const edges = [0, 1, 2, 5, 10, 25, Number.POSITIVE_INFINITY];
  for (let index = 0; index < edges.length - 1; index += 1) {
    if (maxR >= edges[index] && maxR < edges[index + 1]) {
      const max = Number.isFinite(edges[index + 1]) ? edges[index + 1] : Math.max(25.01, maxR + 0.01);
      return { index, min: edges[index], max };
    }
  }
  return { index: edges.length - 2, min: 25, max: Math.max(25.01, maxR + 0.01) };
}

function binCenter(min, max) {
  return (min + max) / 2;
}

async function matchingDecodedRows(filters) {
  const { seenSet, hitSet, evidenceSet } = await candidateSets(filters);
  const limitedSet = intersectCandidateSets([seenSet, hitSet, evidenceSet]);
  const ids = limitedSet ? [...limitedSet] : Array.from({ length: state.manifest.rowCount }, (_, rowId) => rowId);
  const rows = [];
  for (const rowId of ids) {
    if (seenSet && !seenSet.has(rowId)) continue;
    if (hitSet && !hitSet.has(rowId)) continue;
    if (evidenceSet && !evidenceSet.has(rowId)) continue;
    const decoded = applyDataFilters(decodeRow(rowId, await rawRow(rowId)), filters);
    if (!rowMatchesDecoded(decoded, filters)) continue;
    rows.push(decoded);
  }
  return rows;
}

async function summary(message) {
  await loadAtlas(message.manifestPath);
  const filters = message.filters || {};
  const datasetKey = filters.datasetKey || "os";
  const rows = await matchingDecodedRows(filters);
  const groupMap = {
    top: { value: "top", label: "Top score", count: 0, families: blankFamilyCounts() },
    hits: { value: "hits", label: "YPT hits", count: 0, families: blankFamilyCounts() },
    seenNoHit: { value: "seenNoHit", label: "Seen no hit", count: 0, families: blankFamilyCounts() },
    all: { value: "all", label: "All matching", count: 0, families: blankFamilyCounts() },
  };
  const overlap = {
    matching: rows.length,
    osSeen: 0,
    fracSeen: 0,
    families: blankFamilyCounts(),
  };
  const binMap = new Map();

  for (const row of rows) {
    const current = row.dataset[datasetKey] || {};
    overlap.osSeen += row.dataset.os?.seen ? 1 : 0;
    overlap.fracSeen += row.dataset.frac?.seen ? 1 : 0;
    addFamilyCounts(overlap.families, row);

    const activeGroups = [groupMap.all];
    if (row.score >= 0.2) activeGroups.push(groupMap.top);
    if (Number(current.hits || 0) > 0) activeGroups.push(groupMap.hits);
    if (current.seen && !(Number(current.hits || 0) > 0)) activeGroups.push(groupMap.seenNoHit);
    for (const group of activeGroups) {
      group.count += 1;
      addFamilyCounts(group.families, row);
    }

    const xBin = maxRBin(Number(current.maxR || 0));
    const yBin = scoreBin(Number(row.score || 0));
    const family = dominantFamily(row);
    const key = `${xBin.index}:${yBin.index}:${family}`;
    const existing = binMap.get(key) || {
      maxRMin: xBin.min,
      maxRMax: xBin.max,
      scoreMin: yBin.min,
      scoreMax: yBin.max,
      maxRCenter: binCenter(xBin.min, xBin.max),
      scoreCenter: binCenter(yBin.min, yBin.max),
      family,
      count: 0,
    };
    existing.count += 1;
    binMap.set(key, existing);
  }

  return {
    type: "summary",
    seq: message.seq,
    landscape: [groupMap.top, groupMap.hits, groupMap.seenNoHit, groupMap.all],
    overlap,
    quadrantBins: [...binMap.values()],
  };
}

async function loadProteinIndexPart(partIndex) {
  if (state.proteinIndexCache.has(partIndex)) return state.proteinIndexCache.get(partIndex);
  const rel = state.manifest.proteinIndexParts?.[partIndex];
  if (!rel) return null;
  const payload = await fetchJson(atlasUrl(`${rel}?v=${ASSET_VERSION}`));
  state.proteinIndexCache.set(partIndex, payload);
  return payload;
}

async function loadDomainEntries() {
  if (state.domainEntries) return state.domainEntries;
  state.domainEntries = await fetchJson(atlasUrl(`${state.manifest.interproEntries}?v=${ASSET_VERSION}`));
  return state.domainEntries;
}

async function loadProteinDomainPart(partIndex) {
  if (state.proteinDomainCache.has(partIndex)) return state.proteinDomainCache.get(partIndex);
  const rel = state.manifest.proteinDomainParts?.[partIndex];
  if (!rel) return null;
  const payload = await fetchJson(atlasUrl(`${rel}?v=${ASSET_VERSION}`));
  state.proteinDomainCache.set(partIndex, payload);
  return payload;
}

async function loadDomainIndexPart(partIndex) {
  if (state.domainIndexCache.has(partIndex)) return state.domainIndexCache.get(partIndex);
  const rel = state.manifest.domainIndexParts?.[partIndex];
  if (!rel) return null;
  const payload = await fetchJson(atlasUrl(`${rel}?v=${ASSET_VERSION}`));
  state.domainIndexCache.set(partIndex, payload);
  return payload;
}

function decodeDomainInterval(interval, entries) {
  const entry = entries[interval[0]] || {};
  return {
    entryId: interval[0],
    interproId: entry.interproId || "",
    name: entry.name || entry.interproId || "Domain",
    memberDatabases: entry.memberDatabases || "",
    memberSignatures: entry.memberSignatures || "",
    members: entry.members || [],
    suppressedMembers: entry.suppressedMembers || [],
    start: interval[1],
    end: interval[2],
  };
}

function entrySearchText(entry) {
  const members = entry.members || [];
  return [
    entry.interproId,
    entry.name,
    entry.memberDatabases,
    entry.memberSignatures,
    ...members.flatMap((member) => [member.interproId, member.name, member.memberDatabases, member.memberSignatures]),
  ]
    .join(" ")
    .toUpperCase();
}

async function protein(message) {
  await loadAtlas(message.manifestPath);
  const entries = await loadDomainEntries();
  const uniprotId = state.dictionaries.uniprots?.indexOf(message.uniprot);
  if (uniprotId == null || uniprotId < 0) {
    return { type: "protein", seq: message.seq, rows: [], domains: [] };
  }
  const partSize = state.manifest.proteinIndexPartSize || 5000;
  const part = await loadProteinIndexPart(Math.floor(uniprotId / partSize));
  const domainPart = await loadProteinDomainPart(Math.floor(uniprotId / partSize));
  const rowIds = part?.proteins?.[String(uniprotId)] || [];
  const domainIntervals = domainPart?.proteins?.[String(uniprotId)] || [];
  const rows = [];
  for (const rowId of rowIds) {
    rows.push(decodeRow(rowId, await rawRow(rowId)));
  }
  return {
    type: "protein",
    seq: message.seq,
    rows,
    domains: domainIntervals.map((interval) => decodeDomainInterval(interval, entries)),
  };
}

function domainRowPasses(row, filters) {
  const dataset = filteredDatasetForRow(row, filters);
  if (filters.seenOnly && !dataset.seen) return false;
  if (dataFiltersActive(filters) && !dataset.seen) return false;
  if (filters.hitOnly && !(Number(dataset.hits || 0) > 0)) return false;
  return true;
}

async function domainFilterSets(filters) {
  const datasetKey = filters.datasetKey || "os";
  return {
    seenSet: filters.seenOnly || dataFiltersActive(filters) ? await positiveSet(`seen_${datasetKey}`) : null,
    hitSet: filters.hitOnly && !dataFiltersActive(filters) ? await positiveSet(`hits_${datasetKey}`) : null,
  };
}

function domainRowIdPassesSets(rowId, sets) {
  if (sets.seenSet && !sets.seenSet.has(rowId)) return false;
  if (sets.hitSet && !sets.hitSet.has(rowId)) return false;
  return true;
}

function scoreDomainRows(rows, filters) {
  const datasetKey = filters.datasetKey || "os";
  return rows
    .map((row) => {
      const current = filteredDatasetForRow(row, filters);
      return { ...row, domainMaxR: Number(current.maxR || 0), domainHits: Number(current.hits || 0), domainObserved: Number(current.observed || 0) };
    })
    .sort((a, b) => (b.domainMaxR - a.domainMaxR) || (b.domainHits - a.domainHits) || a.gene.localeCompare(b.gene) || (a.site - b.site));
}

async function domainSignalStats(entryId, filters, sets = null) {
  const partSize = state.manifest.domainIndexPartSize || 5000;
  const part = await loadDomainIndexPart(Math.floor(entryId / partSize));
  const rowIds = part?.domains?.[String(entryId)] || [];
  const activeSets = sets || await domainFilterSets(filters);
  const datasetKey = filters.datasetKey || "os";
  const proteins = new Set();
  let filteredSites = 0;
  let seenSites = 0;
  let hitSites = 0;
  let domainMaxR = 0;
  let domainTopSite = "";
  for (const rowId of rowIds) {
    if (!domainRowIdPassesSets(rowId, activeSets)) continue;
    const row = applyDataFilters(decodeRow(rowId, await rawRow(rowId)), filters);
    if (!domainRowPasses(row, filters)) continue;
    const current = filteredDatasetForRow(row, filters);
    const maxR = Number(current.maxR || 0);
    filteredSites += 1;
    proteins.add(row.uniprot);
    if (current.seen) seenSites += 1;
    if (Number(current.hits || 0) > 0) hitSites += 1;
    if (maxR > domainMaxR) {
      domainMaxR = maxR;
      domainTopSite = `${row.gene} Y${row.site}`;
    }
  }
  return { filteredSites, seenSites, hitSites, proteinCount: proteins.size, domainMaxR, domainTopSite };
}

async function domains(message) {
  await loadAtlas(message.manifestPath);
  const entries = await loadDomainEntries();
  const search = String(message.search || "").trim().toUpperCase();
  const filters = message.filters || {};
  const sets = await domainFilterSets(filters);
  const candidates = entries
    .map((entry, entryId) => ({ ...entry, entryId }))
    .filter((entry) => {
      if (!search) return true;
      return entrySearchText(entry).includes(search);
    })
    .sort((a, b) => (b.siteCount - a.siteCount) || (b.proteinCount - a.proteinCount) || a.name.localeCompare(b.name))
    .slice(0, search ? 600 : 400);
  const rows = [];
  for (const entry of candidates) {
    const stats = await domainSignalStats(entry.entryId, filters, sets);
    if (stats.filteredSites > 0) rows.push({ ...entry, ...stats });
  }
  rows.sort((a, b) => (b.domainMaxR - a.domainMaxR) || (b.hitSites - a.hitSites) || (b.filteredSites - a.filteredSites) || a.name.localeCompare(b.name));
  return { type: "domains", seq: message.seq, rows: rows.slice(0, 250) };
}

async function selectedDomainIntervals(entryId, proteinIdsByKey) {
  const entries = await loadDomainEntries();
  const partSize = state.manifest.proteinIndexPartSize || 5000;
  const result = {};
  for (const [key, uniprotId] of proteinIdsByKey.entries()) {
    if (uniprotId == null || uniprotId < 0) continue;
    const part = await loadProteinDomainPart(Math.floor(uniprotId / partSize));
    const intervals = (part?.proteins?.[String(uniprotId)] || [])
      .filter((interval) => Number(interval[0]) === entryId)
      .map((interval) => decodeDomainInterval(interval, entries));
    if (intervals.length) result[key] = intervals;
  }
  return result;
}

async function domain(message) {
  await loadAtlas(message.manifestPath);
  const entries = await loadDomainEntries();
  const entryId = Number(message.entryId || 0);
  const partSize = state.manifest.domainIndexPartSize || 5000;
  const part = await loadDomainIndexPart(Math.floor(entryId / partSize));
  const rowIds = part?.domains?.[String(entryId)] || [];
  const filters = message.filters || {};
  const sets = await domainFilterSets(filters);
  const decoded = [];
  const families = blankFamilyCounts();
  const proteins = new Set();
  const proteinIdsByKey = new Map();
  let hitSites = 0;
  let seenSites = 0;
  let domainMaxR = 0;
  let domainTopSite = "";
  for (const rowId of rowIds) {
    if (!domainRowIdPassesSets(rowId, sets)) continue;
    const raw = await rawRow(rowId);
    const row = applyDataFilters(decodeRow(rowId, raw), filters);
    if (!domainRowPasses(row, filters)) continue;
    const key = `${row.gene}|${row.uniprot}`;
    const current = filteredDatasetForRow(row, filters);
    const maxR = Number(current.maxR || 0);
    proteins.add(row.uniprot);
    proteinIdsByKey.set(key, raw[state.manifest.rowIndex.uniprot]);
    if (current.seen) seenSites += 1;
    if (Number(current.hits || 0) > 0) hitSites += 1;
    if (maxR > domainMaxR) {
      domainMaxR = maxR;
      domainTopSite = `${row.gene} Y${row.site}`;
    }
    addFamilyCounts(families, row);
    decoded.push(row);
  }
  const scoredRows = scoreDomainRows(decoded, filters);
  return {
    type: "domain",
    seq: message.seq,
    entry: entries[entryId] ? { ...entries[entryId], entryId } : null,
    rows: scoredRows.slice(0, 500),
    total: decoded.length,
    proteinCount: proteins.size,
    seenSites,
    hitSites,
    domainMaxR,
    domainTopSite,
    families,
    domainIntervals: await selectedDomainIntervals(entryId, proteinIdsByKey),
  };
}

self.onmessage = async (event) => {
  try {
    const message = event.data || {};
    if (message.type === "query") {
      self.postMessage(await query(message));
      return;
    }
    if (message.type === "summary") {
      self.postMessage(await summary(message));
      return;
    }
    if (message.type === "protein") {
      self.postMessage(await protein(message));
      return;
    }
    if (message.type === "domains") {
      self.postMessage(await domains(message));
      return;
    }
    if (message.type === "domain") {
      self.postMessage(await domain(message));
    }
  } catch (err) {
    self.postMessage({ type: "error", seq: event.data?.seq, message: err.message || String(err) });
  }
};

