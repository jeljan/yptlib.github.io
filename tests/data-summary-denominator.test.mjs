import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";

const manifest = JSON.parse(await readFile(new URL("../assets/data/manifest.json", import.meta.url), "utf8"));

for (const key of ["os", "frac"]) {
  const dataset = manifest.datasets[key];
  const maxSitePromiscuity = Math.max(...dataset.sitePromiscuity);
  assert.equal(maxSitePromiscuity, 100, `${key} site reactivity should use observed compounds as denominator`);

  const minHitCount = Math.min(...Object.values(dataset.drugHitCounts));
  assert.equal(dataset.drugHitCounts[dataset.defaultDrug], minHitCount, `${key} should default to a least-promiscuous compound`);
}
