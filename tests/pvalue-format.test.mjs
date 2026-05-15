import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

const appJs = await readFile(new URL("../assets/app.js", import.meta.url), "utf8");

assert.match(appJs, /if \(num === 0\) return "0\.00e\+0";/, "zero p-values should render in scientific notation");
assert.match(appJs, /num < 0\.00001 \? num\.toExponential\(2\)/, "tiny p-values should render in scientific notation");
