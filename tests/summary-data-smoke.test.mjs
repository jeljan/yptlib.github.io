import assert from "node:assert/strict";

const roots = [
  "https://cdn.jsdelivr.net/gh/jeljan/yptlib@a7e92171d1e6127153f516a19d37a9dceddb5cad/assets/data/",
  "https://cdn.statically.io/gh/jeljan/yptlib/a7e92171d1e6127153f516a19d37a9dceddb5cad/assets/data/",
  "https://raw.githubusercontent.com/jeljan/yptlib/a7e92171d1e6127153f516a19d37a9dceddb5cad/assets/data/",
];
const version = "pages-data-v5";

async function fetchJson(path) {
  let lastError = null;
  for (const root of roots) {
    for (let attempt = 0; attempt < 2; attempt += 1) {
      let timeout = null;
      try {
        const controller = new AbortController();
        timeout = setTimeout(() => controller.abort(), 12000);
        const response = await fetch(`${root}${path}?v=${version}`, { signal: controller.signal });
        if (!response.ok) throw new Error(`${response.status} ${response.statusText}`);
        return response.json();
      } catch (err) {
        lastError = err;
        if (attempt === 0) await new Promise((resolve) => setTimeout(resolve, 250));
      } finally {
        if (timeout) clearTimeout(timeout);
      }
    }
  }
  throw lastError;
}

const manifest = await fetchJson("manifest.json");
const dataset = manifest.datasets.os;
const targetGenes = ["STAT2", dataset.defaultGene].filter(Boolean);

const failures = [];
let loaded = 0;
async function checkGene(gene) {
  const path = dataset.geneFiles[gene];
  if (!path) {
    return `${gene}: missing path`;
  }
  try {
    const payload = await fetchJson(path.replace(/%/g, "%25"));
    assert.ok(Array.isArray(payload.sites), `${gene} should have a sites array`);
    return null;
  } catch (err) {
    return `${gene}: ${err.message}`;
  }
}

const concurrency = 24;
for (let i = 0; i < targetGenes.length; i += concurrency) {
  const chunk = targetGenes.slice(i, i + concurrency);
  const results = await Promise.all(chunk.map(checkGene));
  for (const failure of results) {
    if (failure) failures.push(failure);
    else loaded += 1;
  }
  console.log(`Checked ${Math.min(i + concurrency, targetGenes.length)} / ${targetGenes.length} target genes`);
}

assert.deepEqual(failures, []);
assert.ok(loaded > 0, "should load at least one target gene file");
console.log(`Loaded ${loaded} Summary target gene files successfully.`);
