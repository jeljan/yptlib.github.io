import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

const appJs = await readFile(new URL("../assets/app.js", import.meta.url), "utf8");
const indexHtml = await readFile(new URL("../index.html", import.meta.url), "utf8");

assert.match(
  appJs,
  /const DATA_ROOT = window\.YPTLIB_DATA_ROOT \|\| "https:\/\/raw\.githubusercontent\.com\/jeljan\/yptlib\/a7e92171d1e6127153f516a19d37a9dceddb5cad\/assets\/data\/";/,
  "deployed builds should default to a CORS-enabled static JSON data root"
);

assert.match(
  appJs,
  /const RELEASE_DATA_ARCHIVE_URL = window\.YPTLIB_DATA_ARCHIVE_URL \|\| "";/,
  "release ZIP loading should be opt-in because GitHub release downloads are not CORS-safe for Pages"
);

assert.match(
  appJs,
  /const ASSET_VERSION = "pages-data-v2";/,
  "data fetches should use the current deployment cache key"
);

assert.match(
  indexHtml,
  /<script src="assets\/app\.js\?v=pages-data-v2"><\/script>/,
  "index.html should load app.js with the current cache-busting key"
);
