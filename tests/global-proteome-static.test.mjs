import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import path from "node:path";

const root = path.resolve(".");

async function readJson(relativePath) {
  return JSON.parse(await readFile(path.join(root, relativePath), "utf8"));
}

const manifest = await readJson("assets/data/manifest.json");
assert.equal(
  manifest.globalProteome,
  "global_proteome.json",
  "manifest should expose the global proteome asset path"
);

const payload = await readJson("assets/data/global_proteome.json");
assert.equal(payload.version, 1, "global proteome payload should be versioned");
assert.ok(Array.isArray(payload.points), "payload.points should be an array");
assert.ok(payload.points.length > 1000, "payload should contain a proteome-scale point set");
assert.ok(Array.isArray(payload.clusters), "payload.clusters should be an array");
assert.ok(payload.clusters.length > 0, "payload should include cluster metadata");
assert.ok(
  payload.clusters.some((cluster) => cluster.annotation && Number.isFinite(Number(cluster.annotationPercent))),
  "clusters should include a dominant annotation and prevalence percentage"
);
assert.ok(
  payload.clusters.some((cluster) => Array.isArray(cluster.topAnnotations) && cluster.topAnnotations.length > 1),
  "clusters should include the top two dominant annotations when available"
);
assert.equal(payload.defaults?.threshold, 2, "payload should preserve the default hit threshold");

const seen = new Set();
let activeOverlayCount = 0;
let finiteCoordinateCount = 0;
let annotationCount = 0;
let seenCount = 0;
const clusterCounts = new Map();
const termCountsByCluster = new Map();
for (const point of payload.points) {
  assert.ok(point.uniprot, "each point should include a UniProt identifier");
  assert.ok(!seen.has(point.uniprot), `duplicate UniProt point found: ${point.uniprot}`);
  seen.add(point.uniprot);
  assert.equal(Number.isFinite(point.x), true, `point ${point.uniprot} should have finite x`);
  assert.equal(Number.isFinite(point.y), true, `point ${point.uniprot} should have finite y`);
  finiteCoordinateCount += 1;
  assert.ok(Object.hasOwn(point, "maxRByDataset"), "point should include dataset max-R overlays");
  assert.ok(Object.hasOwn(point, "hitCountByDataset"), "point should include dataset hit-count overlays");
  assert.equal(typeof point.seen, "boolean", "point should expose whether it was seen in the assay data");
  assert.equal(
    point.seen,
    (point.maxRByDataset?.os != null && Number.isFinite(Number(point.maxRByDataset.os))) ||
      (point.maxRByDataset?.frac != null && Number.isFinite(Number(point.maxRByDataset.frac))),
    "seen should require finite assay max-R overlay data"
  );
  if (point.seen) seenCount += 1;
  assert.ok(Array.isArray(point.annotations), "point should expose all compressed InterPro annotations for hover");
  assert.ok(Array.isArray(point.contactTypes), "point.contactTypes should be an array");
  const cluster = String(point.cluster);
  clusterCounts.set(cluster, (clusterCounts.get(cluster) || 0) + 1);
  if (!termCountsByCluster.has(cluster)) termCountsByCluster.set(cluster, new Map());
  const termCounts = termCountsByCluster.get(cluster);
  for (const term of new Set(point.annotations)) {
    termCounts.set(term, (termCounts.get(term) || 0) + 1);
  }
  if ((point.maxRByDataset?.os ?? 0) >= 2 || (point.maxRByDataset?.frac ?? 0) >= 2) {
    activeOverlayCount += 1;
  }
  if (point.annotation) {
    annotationCount += 1;
  }
}

assert.equal(finiteCoordinateCount, payload.points.length, "all points should have finite coordinates");
assert.ok(activeOverlayCount > 0, "at least one protein should have active assay overlay data");
assert.ok(seenCount > 0, "some proteins should be marked as seen in the assay data");
assert.ok(annotationCount > 1000, "InterPro annotations should be present for a substantial share of proteins");

let orderedAnnotationChecks = 0;
for (const point of payload.points) {
  const annotations = point.annotations || [];
  if (annotations.length < 2) continue;
  const termCounts = termCountsByCluster.get(String(point.cluster));
  for (let i = 1; i < annotations.length; i += 1) {
    const prev = termCounts.get(annotations[i - 1]) / clusterCounts.get(String(point.cluster));
    const current = termCounts.get(annotations[i]) / clusterCounts.get(String(point.cluster));
    assert.ok(prev >= current, "annotations should be sorted by within-cluster prevalence");
    orderedAnnotationChecks += 1;
  }
}
assert.ok(orderedAnnotationChecks > 100, "annotation ordering should be validated across many proteins");

const appJs = await readFile(path.join(root, "assets/app.js"), "utf8");
const html = await readFile(path.join(root, "index.html"), "utf8");
assert.match(appJs, /state\.globalProteome/, "runtime state should cache global proteome data");
assert.match(appJs, /loadGlobalProteome/, "runtime should have a global proteome loader");
assert.doesNotMatch(appJs, /Assay activity/, "global highlight modes should not use the old assay activity label");
assert.doesNotMatch(appJs, /value: "cluster"/, "cluster highlight mode should be removed");
assert.doesNotMatch(appJs, /globalHighlightMode/, "global proteome controls should not use the old dropdown");
assert.doesNotMatch(appJs, /Contacts=%\{customdata/, "global hover should not show contact labels");
assert.doesNotMatch(appJs, /Source=%\{customdata/, "global hover should not show source labels");
assert.match(appJs, /showGlobalProteinTooltip/, "global hover should use the app tooltip");
assert.match(appJs, /Cluster annotations:/, "global hover should show dominant cluster annotations");
assert.match(appJs, /Annotations:/, "global hover should show all protein annotations");
assert.match(appJs, /formatProteinAnnotations/, "global hover should format matching annotations specially");
assert.match(appJs, /globalSeenOnly/, "global proteome should support greying unseen background points");
assert.match(appJs, /if \(raw == null\) return null/, "missing max-R values should not display as zero");
assert.match(appJs, /globalThresholdWrap/, "hit threshold control should be conditional");

assert.doesNotMatch(html, /data-tab="global"/, "global proteome should no longer be a top-level tab");
assert.doesNotMatch(html, /id="global" class="tab-panel"/, "global proteome should no longer be a standalone panel");
assert.match(html, /Proteome Explorer/, "global proteome heading should use the requested title");
assert.match(html, /Target Explorer/, "summary target section should use the requested title");
assert.match(html, /Hits by Warhead/, "summary warhead chart should use the requested title");
assert.match(appJs, /Cancer-driver genes/, "target list should use the requested cancer-driver spelling");
assert.doesNotMatch(html, /globalHighlightMode/, "global proteome controls should not use a dropdown");
assert.match(html, /> Hits</, "global proteome should label the hits checkbox");
assert.match(html, /id="globalShowHits"/, "global proteome should expose a hits checkbox");
assert.match(html, /id="globalSeenOnly"/, "global proteome should expose a seen-only checkbox");
assert.match(html, /id="globalShowCancer"/, "global proteome should expose a cancer-driver checkbox");
assert.match(html, /id="globalShowCustom"/, "global proteome should expose a custom genes checkbox");
assert.match(html, /id="globalScatter" class="plot square"/, "global proteome plot should render in a square container");
assert.match(html, /class="global-proteome-layout"/, "global proteome plot should have a side-control layout");
assert.match(html, /id="globalThresholdWrap"/, "hit threshold label should be hideable");
assert.ok(
  html.indexOf('id="globalThresholdWrap"') < html.indexOf('id="globalShowHits"'),
  "hit threshold should appear above the hits checkbox when shown"
);
assert.ok(
  html.indexOf('id="globalScatter"') < html.indexOf('id="globalProteomeControls"'),
  "global proteome controls should sit below the plot"
);
assert.ok(
  html.indexOf('id="compoundHitsBinned"') < html.indexOf('id="globalScatter"') &&
    html.indexOf('id="globalScatter"') < html.indexOf("Target Explorer"),
  "global proteome plot should sit below the four summary plots and above Target Explorer"
);
assert.match(appJs, /el\("summary"\)\?\.classList\.contains\("active"\)[\s\S]*renderGlobalProteome/, "summary dataset refresh should render the proteome view");
assert.match(appJs, /tabName === "summary"[\s\S]*renderGlobalProteome/, "summary tab activation should render the proteome view");
