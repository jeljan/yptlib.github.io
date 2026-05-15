# TASKS

## Production Ship - 2026-05-14

Completed and shipped-review-ready in this workspace:

- Fixed the initial Compound `Hit threshold` dropdown freeze by updating the selected compound label immediately and rebuilding full dropdown counts in cancellable, browser-yielded batches.
- Added background compound-row cache warming after dataset load so the first threshold change can use hot compound data instead of cold-loading every compound.
- Removed the Compound `Update list` fallback control.
- Bumped the shipped static asset cache key to `pages-data-v11`.
- Changed Compound tab inline fallback action from `Refresh list` to `Update list`.
- Tightened Site tab structure panel spacing so the structure box sits one standard 16px gap below the site analysis row.
- Synced the local review copy in `output/local-review/`.

Verification run:

```powershell
node --check assets/app.js
node --check output/local-review/assets/app.js
node tests\local-review-static.test.mjs
node tests\deploy-config.test.mjs
node tests\global-proteome-static.test.mjs
```

Browser smoke verified `assets/app.js?v=pages-data-v11`, no `Update list` control, a 0 ms warmed threshold input dispatch, immediate selected-label update to the new count, and all 740 options retained.

## Immediate Next Tasks - 2026-05-14

The user asked to update this file for immediate continuation in a new window. The following items are the active unfinished request.

### Browser Diff Comments To Implement

- Move the Summary `Site Hits` donut legend to the right/top to match `Hits by Warhead`, while keeping it clear of the donut.
- Make the `Reactive Site Distribution` bin details easier to read than the current scroll-heavy hover. Preferred implementation: click a bar to open a text/list modal of all sites in that bin.
- Remove the legend title `Selectivity` from `Reactive Compound Selectivity`.
- Add click popups/modals for the two Summary compound charts:
  - `Hits by Warhead`
  - `Reactive Compound Selectivity`
  Each popup should show all compound structures in the clicked bar segment/bin.
- Make Proteome Explorer overlaid points easier to hover at default zoom without making visible markers larger. Preferred implementation: add a larger transparent hover-target trace.
- Fix Compound tab `Hit threshold` so threshold changes update the plot and compound dropdown labels even when `Significant only` is not selected.
- In Compound profile/volcano hovers, replace `Sites hit` with `Site rank`, formatted as `R-rank/total sites`.
- Change Compound volcano axis labels to literal `log_2R` and `-log_10P`.
- Normalize whitespace around the Site tab `Source` link and `No structural contacts found.` message so both have matching top/bottom spacing.

### Data / Manifest Tasks

- Remove `SRPKIN-1-1` and `SRPKIN-1-2` from compound lists and generated assets.
- Known generated SRPKIN assets found:
  - `assets/data/frac/compounds/SRPKIN-1-1_part_0.json`
  - `assets/data/frac/compounds/SRPKIN-1-1_part_1.json`
  - `assets/data/frac/compounds/SRPKIN-1-1_part_2.json`
  - `assets/data/frac/hover/SRPKIN-1-1.json`
  - `assets/data/os/compounds/SRPKIN-1-1_part_0.json`
  - `assets/data/os/compounds/SRPKIN-1-1_part_1.json`
  - `assets/data/os/compounds/SRPKIN-1-1_part_2.json`
  - `assets/data/os/compounds/SRPKIN-1-1_part_3.json`
  - `assets/data/os/compounds/SRPKIN-1-2_part_0.json`
  - `assets/data/os/compounds/SRPKIN-1-2_part_1.json`
  - `assets/data/os/compounds/SRPKIN-1-2_part_2.json`
  - `assets/data/os/compounds/SRPKIN-1-2_part_3.json`
  - `assets/data/os/hover/SRPKIN-1-1.json`
  - `assets/data/os/hover/SRPKIN-1-2.json`
  - `structures/SRPKIN-1.png`
- Also remove SRPKIN entries from `assets/data/manifest.json`, including `rawDrugs`, `drugChoices`, `drugPromiscuity`, `drugHitCounts`, `compoundTypes`, `compoundParts`, `rawHoverParts`, and the parallel `compoundPromiscuityRecords` array.
- Add a checkbox below the one-shot/fractionated toggle to drop compounds with no replicates.
- Derive the no-replicate compound list from the pipeline pair-QC metadata:
  - `C:\Users\ethan\.codex\worktrees\49ce\yptlib\data\raw_hover\ypt_functions.py`
  - `C:\Users\ethan\.codex\worktrees\49ce\yptlib\data\raw_hover\rc2.csv`
  - `C:\Users\ethan\.codex\worktrees\49ce\yptlib\data\raw_hover\ypt_samples.csv`
  - `C:\Users\ethan\.codex\worktrees\49ce\yptlib\data\raw_hover\data_list.csv`
- Store the result in the manifest, for example `noReplicateDrugs` per dataset, then filter compound dropdowns and compound summary charts when the checkbox is selected.

### Files Likely To Touch

- `index.html`
- `assets/app.js`
- `assets/styles.css`
- `assets/data/manifest.json`
- SRPKIN files under `assets/data/**`
- `structures/SRPKIN-1.png`

### Verification For This Batch

Run at minimum:

```powershell
node --check assets/app.js
node tests\deploy-config.test.mjs
node tests\summary-data-smoke.test.mjs
```

If a local server is used:

```powershell
python -m http.server 8000
```

Then smoke the Summary, Compound, Site, and Proteome Explorer views. Confirm the modal bar-click behavior, threshold-driven dropdown labels, and Proteome hover targets.

### Shipping Notes

- Bump `ASSET_VERSION` in `assets/app.js`, likely to `pages-data-v8`, to avoid stale browser cache.
- Current repo may require a Git safe-directory fix in the sandbox before git commands work.
- Push final production changes to `main` for GitHub Pages after verification.

## Completed From Static-App Worktree

- Added static app files: `index.html`, `assets/app.js`, `assets/styles.css`.
- Added `build_static_data.py` for compact browser assets.
- Added `build_raw_hover_data.py` for raw DMSO/compound duplicate SN hover assets.
- Kept `app.py` as the legacy/reference Shiny implementation.
- Added Summary, Compound, Site, Global Proteome, and Bioinformatics tabs.
- Added one-shot/fractionated dataset switching.
- Added significance, high-variance, minimum-SN, and selectivity/hit-count filters.
- Replaced percent-promiscuity filters with number-of-sites-hit selectivity controls where requested.
- Added Summary four-plot layout based on the earlier Shiny code.
- Normalized compound type labels to `Sulfonyl Fluorides`, `Fluorosulfates`, and `Sulfonyl Triazoles`.
- Treated SRPKIN-1-like compounds as `Sulfonyl Fluorides`.
- Updated Summary and Site bar hovers to include compound structures and raw SN duplicate values.
- Fixed Compound hover Mean SN to use the correct mean compound SN field.
- Updated Compound hover to show UniProt rather than repeating gene symbol.
- Added p-value gradient color mode and contact subtype color modes.
- Added direct links for PDB and AlphaFold structures.
- Consolidated Site structure controls into Structure and Contact neighborhood dropdowns.
- Kept Mol* viewers inline and automatically loaded for the selected Site view.
- Reduced viewer whitespace and prevented filter-only changes from unnecessarily refreshing structures.
- Deduplicated Summary site dropdown entries after filters are applied.
- Added `Sites Seen` and `Max R` annotations to Site gene dropdown labels.
- Limited the Compound profile threshold line to the span where points exist.
- Added lightweight Volcano labels for top threshold-crossing sites.
- Made Structure default to a PDB entry when mapped PDB structures are present.
- Filtered Contact neighborhood dropdown entries to those with real contact distances.
- Added synchronized shared filters across Summary, Compound, and Site.
- Added Site selectivity filter synchronized with Summary selectivity.
- Restored Summary/Site bars to show Top 20 when more than 20 compounds remain, or Top N when fewer remain.
- Restored the compound structure card in the Compound sidebar.
- Added styled custom hover cards with Site SN raw-hover bars for Compound profile and volcano points.
- Removed Plotly editable title/subtitle placeholders with CSS while keeping draggable annotation editing.
- Updated Compound color modes to sentence case and reordered contact subtype options.
- Reordered Summary contact targets to match Compound contact subtype order.
- Removed misleading mean Site SN fallback bars from raw-intensity hover cards.
- Fixed raw-hover lookup for related catalog rows and fixed raw-hover asset generation by preserving cleaned peptide sequences that contain `K`.
- Regenerated raw hover assets under `assets/data/{os,frac}/hover/`.
- Added `activeGeneFilterIndex` assets and wired Site View to use them when `>= 1 Compound hit only` is enabled.
- Updated Compound volcano to use the same quality-filtered rows as Compound profile.
- Added runtime cache busting for generated `assets/data/*` JSON and `assets/app.js`.
- Added a cache-busted local review copy under `output/local-review` so the in-app browser can review the full static app while the main checkout remains pruned for Pages.
- Updated review copy and `49ce` static-app worktree for latest browser comments:
  - Renamed `Site Reactivity Distribution` to `Reactive Site Distribution`.
  - Renamed `Compound Selectivity Profiles` to `Reactive Compound Selectivity`.
  - Moved the Site Hits donut legend below the plot to prevent legend/plot overlap.
  - Forced Plotly containers to fill panel width and resize after render to reduce large right-side whitespace in Compound profile/volcano panels.
  - Added loading state and disabled toggle behavior for slow dataset switches.
  - Increased hover collision spacing and repositioned hovers after raw values load.
  - Changed raw hover title from `Site SN` to `TMT SN`.
  - Put `R`, `P-value`, and `Sites Hit` in a one-row hover stat layout where applicable.
- Added `tests/local-review-static.test.mjs` to catch the local-review title/loading/hover-label hooks.

## Completed From Pages Deployment Worktree

- Fixed GitHub Pages deployment for the small static site.
- Kept large `assets/data/` out of the Pages branch to stay below GitHub Pages limits.
- Disabled unreliable browser loading from GitHub Release ZIP assets.
- Added external JSON data roots with fallback order: jsDelivr, Statically, raw.githubusercontent.com.
- Hardened external data loading with fallback, retries, timeout, and bounded Target Engagement concurrency.
- Bumped the app script/data cache key to `pages-data-v5`.
- Added deployment/config smoke tests:
  - `tests/deploy-config.test.mjs`
  - `tests/summary-data-smoke.test.mjs`
- Verified Target Engagement plots on the live site, including the reported `STAT2` path.
- Confirmed GitHub Pages build status reported `built`.

## Current Known Limitations

- Full Shiny parity is not complete; current scope is core workflow parity plus requested UI improvements.
- Contact neighborhood highlighting is limited by `pdb_interactions.parquet`, which stores nearest contact residues, not all residues within 5A.
- If true 5A neighborhood highlighting is required, generate a richer coordinate-derived asset from PDB/mmCIF coordinates.
- PDBe Mol* is asked to focus and show selected residues as side chains, but exact sidechain-facing camera orientation is not exposed by the current simple viewer code.
- Bioinformatics tab has placeholder UI only; conservation/trait logic is deferred until data is provided.
- Static app depends on CDN/network access for Plotly, PDBe Mol*, and deployed JSON data.
- Python Playwright is available; `npx`/Playwright CLI may not be installed globally.
- The generated `activeGeneFilterIndex` currently supports the active-hit Site gene dropdown path. If non-active Site gene labels must become fully filter-aware too, add a separate compact index for all-site counts.
- `output/local-review` is a disposable review surface, not the canonical deployment checkout. Canonical static-app edits are mirrored into `C:\Users\ethan\.codex\worktrees\49ce\yptlib`.

## High-Priority Next Tasks

- Monitor whether users still report data fetch failures after a hard refresh on `pages-data-v5`.
- If GitHub-hosted CDN fallbacks remain flaky, move generated JSON data to a durable object store/CDN.
- Add a browser automation smoke test that can be run locally before each Pages deploy.
- Consider reducing first-load pressure further by precomputing compact Summary Target Engagement indexes for each target list.
- Decide whether to generate full 5A residue-neighborhood assets for contact viewers.
- If full 5A neighborhoods are needed, add an offline builder that parses PDB/mmCIF coordinates and writes all residues within the chosen radius for each mapped site/structure.
- Continue UI review in browser and respond to any new direct diff comments.
- Spot-check static outputs against Shiny calculations for top compounds, Max R, significance, SN filters, and contact filters.
- Check Summary/Site/Compound hover cards in multiple viewport sizes after the latest tooltip collision changes, especially near the bottom edge.
- Verify Mol* selection behavior on several PDB and AlphaFold examples, including no-PDB and no-contact cases.

## Useful Verification Commands

Static app checks:

```powershell
python build_static_data.py
python build_raw_hover_data.py
node --check assets/app.js
python -m py_compile build_static_data.py build_raw_hover_data.py
python -m http.server 8765
```

Local review checks:

```powershell
node tests\local-review-static.test.mjs
node --check output\local-review\assets\app.js
node --check C:\Users\ethan\.codex\worktrees\49ce\yptlib\assets\app.js
```

Current in-app browser review URL:

```text
http://127.0.0.1:8765/index.html?review=raw-hover-v4
```

Deployment checks:

```powershell
node tests\deploy-config.test.mjs
node tests\summary-data-smoke.test.mjs
& 'C:\Program Files\GitHub CLI\gh.exe' api repos/jeljan/yptlib/pages/builds/latest
```

Python Playwright smoke-test pattern:

```powershell
@'
from playwright.sync_api import sync_playwright
with sync_playwright() as p:
    browser = p.chromium.launch(channel='chrome', headless=True)
    page = browser.new_page(viewport={'width': 1355, 'height': 897})
    page.goto('http://127.0.0.1:8765/index.html', wait_until='domcontentloaded', timeout=30000)
    page.wait_for_timeout(3500)
    print(page.locator('h1').inner_text())
    print(page.locator('#summaryTargetList option').evaluate_all('(opts) => opts.map(o => o.textContent).join("|")'))
    browser.close()
'@ | python -
```

## Data Questions To Preserve

- Raw hover data is generated from the current `data/raw_hover/os` and `data/raw_hover/frac` files.
- Each plex of 32 `sn_sum` channels maps DMSO duplicates and compound duplicates as described by the user:
  - Channels 0 and 1 are DMSO duplicates for the first seven compounds.
  - Compound duplicate pairs follow as 2/3, 4/5, and so on.
  - Channels 30 and 31 are DMSO duplicates for the last seven compounds.
  - Compound identities per plex are referenced through `ypt_samples`.
- If raw-hover values look wrong, inspect `build_raw_hover_data.py` first and compare against `ypt_functions.py` / pipeline logic.
- The previous raw-hover miss was caused by copying the pipeline's `clean_seq contains K` exclusion into `build_raw_hover_data.py`.
- After removing that exclusion and regenerating hover assets, audited coverage was 100% for tested OS/frac combinations, including strict Min SN/significance/variance/selectivity cases.
- Specific checked case: `Z4262428218`, significant only, hide high variance, Min SN 10 had 296 filtered Compound View points and zero missing raw-hover rows.

## Worktree With No Actionable Notes

`C:\Users\ethan\.codex\worktrees\d233\yptlib` only contains a minimal `README.md` title (`# yptlib.github.io`). No `HANDOFF.md`, `TASKS.md`, TODOs, or follow-up notes were found there.
