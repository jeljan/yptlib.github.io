import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

const appJs = await readFile(new URL("../assets/app.js", import.meta.url), "utf8");

assert.match(appJs, /FS:\s*"Fluorosulfates"/, "FS raw type should display as Fluorosulfates");
assert.match(appJs, /SuFEx:\s*"Sulfonyl Fluorides"/, "SuFEx raw type should display as Sulfonyl Fluorides");
assert.match(
  appJs,
  /"Sulfonyl Fluorides":\s*"Fluorosulfates"/,
  "Already-expanded Sulfonyl Fluorides labels should be swapped to Fluorosulfates"
);
assert.match(
  appJs,
  /Fluorosulfates:\s*"Sulfonyl Fluorides"/,
  "Already-expanded Fluorosulfates labels should be swapped to Sulfonyl Fluorides"
);
