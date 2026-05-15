import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

const appJs = await readFile(new URL("../assets/app.js", import.meta.url), "utf8");
const indexHtml = await readFile(new URL("../index.html", import.meta.url), "utf8");

assert.match(
  appJs,
  /const DATA_ROOTS = window\.YPTLIB_DATA_ROOT/,
  "deployed builds should support an overrideable external static JSON data root"
);

assert.match(
  appJs,
  /"assets\/data\/"/,
  "deployed builds should try the co-located static JSON data first"
);

assert.match(
  appJs,
  /https:\/\/cdn\.jsdelivr\.net\/gh\/jeljan\/yptlib@a7e92171d1e6127153f516a19d37a9dceddb5cad\/assets\/data\//,
  "deployed builds should default to a CORS-enabled CDN data root"
);

assert.match(
  appJs,
  /https:\/\/cdn\.statically\.io\/gh\/jeljan\/yptlib\/a7e92171d1e6127153f516a19d37a9dceddb5cad\/assets\/data\//,
  "deployed builds should include a second CORS-enabled CDN fallback"
);

assert.match(
  appJs,
  /for \(const root of DATA_ROOTS\)/,
  "data fetches should try fallback roots before failing"
);

assert.match(
  appJs,
  /const RELEASE_DATA_ARCHIVE_URL = window\.YPTLIB_DATA_ARCHIVE_URL \|\| "";/,
  "release ZIP loading should be opt-in because GitHub release downloads are not CORS-safe for Pages"
);

assert.match(
  appJs,
  /const ASSET_VERSION = "pages-data-v14";/,
  "data fetches should use the current deployment cache key"
);

assert.match(
  indexHtml,
  /<script src="assets\/app\.js\?v=pages-data-v14"><\/script>/,
  "index.html should load app.js with the current cache-busting key"
);

assert.match(
  appJs,
  /const FETCH_RETRIES_PER_ROOT = 2;/,
  "data fetches should retry each root before falling back"
);

assert.match(
  appJs,
  /const FETCH_TIMEOUT_MS = 12000;/,
  "data fetches should time out instead of hanging on one CDN"
);

assert.match(
  appJs,
  /const SUMMARY_GENE_FETCH_CONCURRENCY = 8;/,
  "Target Engagement should limit concurrent gene data requests"
);

assert.match(
  appJs,
  /await mapWithConcurrency\(uniqueGenes, SUMMARY_GENE_FETCH_CONCURRENCY, \(gene\) => loadGeneSites\(gene\)\);/,
  "Target Engagement should not flood the browser with gene data requests"
);
