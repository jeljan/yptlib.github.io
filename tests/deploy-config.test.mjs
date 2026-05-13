import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

const appJs = await readFile(new URL("../assets/app.js", import.meta.url), "utf8");

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
