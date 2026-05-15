# HANDOFF

## Latest Production Ship - 2026-05-14

Latest completed polish before shipping:

- Compound `Hit threshold` changes no longer synchronously rebuild every dropdown count on the input event. The selected option label updates immediately, then the full list refresh runs in cancellable batches that yield back to the browser.
- Compound rows now warm in the background after dataset load, with load de-duplication between warmup and foreground refreshes, so the first threshold edit gets the same hot-cache behavior as later edits.
- Removed the Compound `Update list` fallback control.
- Static asset cache key bumped to `pages-data-v11`.
- Compound tab inline fallback action now reads `Update list`.
- Site tab structure/contact panel now has a 16px top gap from the site analysis row, with `.site-structure-panel` margin reset to `0`.
- Local review mirror refreshed under `output/local-review/`.

Verification run:

```powershell
node --check assets/app.js
node --check output/local-review/assets/app.js
node tests\local-review-static.test.mjs
node tests\deploy-config.test.mjs
node tests\global-proteome-static.test.mjs
```

Browser smoke verified `assets/app.js?v=pages-data-v11`, no `Update list` control, immediate selected-label threshold count update, and a warmed first threshold dispatch.

## Immediate Continuation - 2026-05-14

Current workspace:

- Canonical checkout: `C:\Users\ethan\Documents\GitHub\yptlib`
- Shell: PowerShell
- Writable root: `C:\Users\ethan\Documents\GitHub\yptlib`
- Current target branch from prior work: `codex/yptlib`
- Live site under review: `https://jeljan.github.io/yptlib/`
- Production branch is `main`; the previous pushed production commit was `041e453f Refine compound and site views`.
- `git status` currently failed in the sandbox because Git marks the repo as dubious ownership for `CodexSandboxOffline`. Use the normal repo-aware path or add a safe.directory exception only if needed.

## GitHub CLI Auth Note - 2026-05-14

`gh` is authenticated for the user in normal PowerShell, but the Codex sandbox may not see that auth directly.

Confirmed details:

- `gh.exe` path: `C:\Program Files\GitHub CLI\gh.exe`
- Interactive PowerShell shows `gh auth status` as logged in to `github.com` account `jeljan` using `(keyring)`.
- Token scopes shown there: `gist`, `read:org`, `repo`, `workflow`.
- The config file `C:\Users\ethan\AppData\Roaming\GitHub CLI\hosts.yml` does not store the token, which is expected for Windows keyring auth.
- In the sandboxed Codex shell, `gh auth token` reported `no oauth token found for github.com`, and `gh api user` failed with a socket-permission error.
- Running `gh auth status -h github.com` outside the sandbox succeeded and showed the keyring login.
- Running `gh api user --jq .login` outside the sandbox succeeded and returned `jeljan`.

Practical instruction for future Codex windows: if a normal sandboxed `gh` command claims auth is missing or network is blocked, do not assume the user is logged out. Request escalation for the specific `gh` command that needs GitHub/keyring/network access. Avoid `--show-token` in commands, because it exposes credentials.

## Recurrent Sandbox Networking Note - 2026-05-14

Sandbox networking is a known recurrent issue in this workspace. Do not spend time debugging the app, CDN URLs, GitHub auth, or localhost server when the only symptom is a sandbox/network denial.

Observed examples:

- `node tests\summary-data-smoke.test.mjs` failed in the sandbox with `fetch failed` / `connect EACCES ...:443`, then passed when rerun with network escalation.
- In-app Browser navigation to `http://127.0.0.1:8766/...` and `http://localhost:8766/...` was blocked with `net::ERR_BLOCKED_BY_CLIENT`.
- Local Chrome/Playwright without network escalation loaded the local shell and data but reported CDN assets such as Plotly/Mol* as `net::ERR_NETWORK_ACCESS_DENIED`.

Future instruction: for networked verification or browser checks that need CDN/GitHub/local browser access, first run the cheap static checks. If a network/browser check fails with one of the above denial patterns, skip further diagnosis and immediately rerun the specific command with `sandbox_permissions="require_escalated"` and a narrow justification. Treat a passing escalated run as authoritative.

The user asked to pause and update `TASKS.md` / `HANDOFF.md` for a new window. Do not continue implementation before rereading the latest browser comments below.

Latest requested work, still pending:

- Summary `Site Hits` donut: move the `0 Hits` / `>=1 Hit` legend to the right/top like `Hits by Warhead`, with enough chart domain/margin so it does not overlap the donut.
- `Reactive Site Distribution`: current hover content is acceptable but hard to scroll. Prefer a text-based full list on bar click, or another reliable way to view all sites in the selected bin.
- `Reactive Compound Selectivity`: remove the legend title `Selectivity`.
- Summary compound charts: for both `Hits by Warhead` and `Reactive Compound Selectivity`, add a click popup/modal that shows all compound structures belonging to the clicked bar segment/bin.
- Proteome Explorer: overlaid points are visually large but still hard to hover at default zoom. Add a larger invisible/near-invisible hover target trace or equivalent without increasing visible marker size.
- Compound tab `Hit threshold`: changing the threshold must update both the plot and compound dropdown labels regardless of whether `Significant only` is selected.
- Compound plot hovers: replace `Sites hit` with `Site rank`, formatted as `R-rank/total sites`.
- Compound volcano axes: use literal `log_2R` for x-axis and `-log_10P` for y-axis.
- Site tab structure/contact area: make `Source` link whitespace and `No structural contacts found.` whitespace consistent above and below.
- Remove `SRPKIN-1-1` and `SRPKIN-1-2` from compound lists and generated assets. Known files include:
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
- Add a checkbox below the one-shot/fractionated toggle to drop compounds with no replicates. This should be based on the `run_pipeline` pair-QC readout where a replicate is completely dropped due to low correlation. The mapping still needs to be finalized from raw pipeline metadata.

Useful discoveries already made for the replicate checkbox:

- Full raw-hover/pipeline worktree: `C:\Users\ethan\.codex\worktrees\49ce\yptlib\data\raw_hover`
- `ypt_functions.py` has `run_pipeline()` and `_summarize()`; it prints:
  - `PAIR QC (rho >= 0.95)`
  - statuses like `full`, `rep0`, `rep1`, and `dead`
- `rc2.csv` is read with `index_col=0`, rows are TMT channels / sample rows, columns are plex numbers.
- `ypt_samples.csv` columns include `Sample ID`, `Sample Name`, and `Plex Number`. The first rows show duplicate channel assignments per compound, for example DMSO and `SRPKIN-1` in plex 1.
- `data_list.csv` maps `Plex Number` to one-shot/fractionated run IDs.
- Next step: inspect `ypt_functions.py` around sample/channel mapping and derive which compound duplicate pairs are `dead` / no usable replicate per dataset. Then store those compound IDs in the manifest, e.g. `manifest.datasets.os.noReplicateDrugs` and `manifest.datasets.frac.noReplicateDrugs`, and wire the new checkbox to filter compound dropdowns/compound summary charts.

Commands already used for discovery:

```powershell
rg -n "SRPKIN|replicate|correlation|low correlation|drop|dropped|run_pipeline|drugHitCounts|compoundPromiscuity|rawDrugs|compoundHitsBinned|compoundHitsBinary|globalOverlay|overlay|hovertemplate|showCompoundPointTooltip|compoundPromiscuityRecords" -S .
Get-Content C:\Users\ethan\.codex\worktrees\49ce\yptlib\data\raw_hover\ypt_samples.csv -TotalCount 8
Get-Content C:\Users\ethan\.codex\worktrees\49ce\yptlib\data\raw_hover\data_list.csv -TotalCount 8
rg -n "exp_names|name_list|control|rep|rc2|PAIR QC|sample" C:\Users\ethan\.codex\worktrees\49ce\yptlib\data\raw_hover\ypt_functions.py
```

Likely files to edit:

- `index.html`
- `assets/app.js`
- `assets/styles.css`
- `assets/data/manifest.json`
- possibly delete generated SRPKIN JSON/PNG assets after verifying resolved paths are inside the workspace.

Implementation notes:

- Bump `ASSET_VERSION` in `assets/app.js` after changing shipped data/UI, likely from `pages-data-v7` to `pages-data-v8`.
- For summary bar popups, add a shared modal in `index.html` and CSS. Use clicked Plotly bar `customdata` to populate compound cards with `structures/{drug}.png`.
- For site-distribution bar detail, reuse the same modal and show a text/list of all sites in the bin.
- For Proteome Explorer hoverability, keep visible glow/core sizes as requested earlier, but add an invisible hover target trace with larger marker size and the hovertemplate/customdata.
- For compound threshold updates, avoid stale caches. Clear compound choice count cache on threshold/filter changes and ensure dropdown labels always derive from the active threshold/filter state.
- For SRPKIN removal, update manifest structured arrays/objects and remove matching generated files. `compoundPromiscuityRecords` is parallel to `rawDrugs`, so remove matching indexes carefully.
- Preserve existing user-requested titles:
  - `Proteome Explorer`
  - `Target Explorer`
  - `Hits by Warhead`
  - `Hit threshold`

Verify before shipping:

```powershell
node --check assets/app.js
node tests\deploy-config.test.mjs
node tests\summary-data-smoke.test.mjs
python -m http.server 8000
```

If shipping to production:

```powershell
git add index.html assets/app.js assets/styles.css assets/data/manifest.json assets/data structures
git commit -m "Refine summary and compound interactions"
git checkout main
git merge codex/yptlib
git push origin main
& 'C:\Program Files\GitHub CLI\gh.exe' api repos/jeljan/yptlib/pages/builds/latest
```

## Worktrees Checked

Relevant project worktrees found with `git worktree list`:

- `C:\Users\ethan\.codex\worktrees\49ce\yptlib`
  - Detached at `171b5aaf` (`Updates`).
  - Contains the full static-app migration work and generated data/build scripts.
  - Has explicit `HANDOFF.md` and `TASKS.md`.
- `C:\Users\ethan\.codex\worktrees\c187\yptlib`
  - Detached at `ef301d17` (`Bump target engagement cache key`), matching `origin/main`.
  - Contains the GitHub Pages deployment fix.
  - Has explicit `HANDOFF.md` and `TASKS.md`.
- `C:\Users\ethan\.codex\worktrees\d233\yptlib`
  - Detached at `171b5aaf` (`Updates`).
  - Only contains `README.md` with `# yptlib.github.io`.
  - No `HANDOFF.md`, `TASKS.md`, or inline task notes were found.

The main checkout at `C:\Users\ethan\Documents\GitHub\yptlib` also had older untracked handoff/task notes about release-zip hosting. Those notes are superseded by the `c187` Pages fix.

## Current Deployment State

- GitHub Pages is enabled for `jeljan/yptlib` from `main` at `/`.
- Live site: `https://jeljan.github.io/yptlib/`
- Latest deployed commit verified in the deployment-fix worktree: `ef301d1733ab328a29caf9bb6c6c7350d6e0594d`.
- The Pages site intentionally does not include `assets/data/` because the full static data tree exceeds GitHub Pages size limits.
- The current Pages app does not rely on the older GitHub Release ZIP approach. Browser release-asset fetches were unreliable, so data loading was changed to external JSON roots with fallback order:
  1. jsDelivr
  2. Statically
  3. raw.githubusercontent.com
- JSON fetches include retry and timeout handling.
- Summary Target Engagement gene-data fetch concurrency is bounded to avoid overwhelming the browser/CDN.
- App script/data cache key is `pages-data-v5`.

Users with an older cached script should hard refresh or open:

```text
https://jeljan.github.io/yptlib/?fresh=target-v5
```

First load can still be slower than an all-same-origin deployment because the app depends on external CDN JSON files. If CDN reliability remains an issue, the next structural fix is to host `assets/data` on a dedicated object/CDN service.

## Static App Shape

The full static-app migration lives in the `49ce` worktree.

Primary files:

- `index.html`: static app entry point.
- `assets/app.js`: runtime UI, filtering, plotting, hover, data loading, and Mol* logic.
- `assets/styles.css`: styling.
- `build_static_data.py`: generates compact static browser data under `assets/data/`.
- `build_raw_hover_data.py`: generates per-compound raw hover assets under `assets/data/{os,frac}/hover/`.
- `app.py`: legacy/reference Shiny implementation.

Serve locally from the full static worktree with:

```powershell
python -m http.server 8765
```

Then open:

```text
http://127.0.0.1:8765/index.html
```

Current local review setup:

- Port `8765` is currently served from `C:\Users\ethan\Documents\GitHub\yptlib\output\local-review`.
- `output/local-review` is a disposable cache-busted review copy of the full static app.
- It links `assets/data` and `structures` back to the full generated assets in `C:\Users\ethan\.codex\worktrees\49ce\yptlib`.
- Use this URL in the in-app browser for current review:

```text
http://127.0.0.1:8765/index.html?review=raw-hover-v4
```

Why this exists: the main checkout is the pruned Pages checkout and does not contain full `assets/data`. Serving the main checkout directly loads the shell but fails app data fetches. The review copy avoids stale browser cache by loading `assets/app.js?v=raw-hover-v4`.

## Data Assets

Important generated assets from the static migration:

- `assets/data/manifest.json`: dataset metadata, choices, counts, and asset paths.
- `assets/data/os/` and `assets/data/frac/`: one-shot and fractionated data.
- `assets/data/*/compound/`: lazy-loaded per-compound data.
- `assets/data/*/genes/`: lazy-loaded per-gene site summaries.
- `assets/data/*/contacts/`: contact/PDB indices keyed by UniProt.
- `assets/data/*/hover/`: raw DMSO/compound duplicate SN values for hover cards.
- `assets/data/*/gene_filter_index_active.json`: compact active-hit gene index used by Site View when `>= 1 Compound hit only` is enabled.

Raw hover generation expects source files in `data/raw_hover/os` and `data/raw_hover/frac`. These correspond to the original `os_0.7` and `frac_0.7` folders.

Important raw-hover caveat: `build_raw_hover_data.py` must not drop peptides merely because the cleaned sequence contains `K`. That old pipeline filter caused real raw SN rows, for example `ANXA5_Y297 / Z4607531286`, to be omitted from hover assets.

## Implemented App Behavior

- Tabs: Summary, Compound, Site, Global Proteome, Bioinformatics.
- Summary tab includes four plots, compact filters, human-readable target labels, synchronized quality filters, deduplicated site dropdown labels, and raw duplicate TMT SN hover cards.
- Summary panel titles currently include `Site Hits`, `Reactive Site Distribution`, `Compound Hits By Warhead`, and `Reactive Compound Selectivity`.
- Summary `Site Hits` donut legend is horizontal below the plot to avoid legend/plot collision.
- Compound tab includes quality-filtered profile/volcano plots, p-value gradient coloring, contact subtype coloring, custom gene list mode, singular warhead labels, draggable Plotly annotations, structure cards, and raw hover cards without mean-SN fallback.
- Compound profile/volcano Plotly containers are forced to fill panel width and resize after render to avoid large right-side whitespace.
- Dataset switching disables the toggle while loading and changes the subtitle to `Loading One-Shot data...` or `Loading Fractionated data...` so the old dataset label is not shown during slow switches.
- Summary, Site, Compound profile, and Compound volcano hover cards use stronger viewport collision spacing; hovers are repositioned after raw TMT values load.
- Hover stats show `R`, `P-value`, and `Sites Hit` on one row where applicable, and raw duplicate bars are labeled `TMT SN`.
- Site tab includes alphabetic gene sorting, `Sites Seen` and `Max R` labels, selectivity filters synchronized with Summary, active-hit gene filtering, inline Mol* viewers, one Structure dropdown, one Contact neighborhood dropdown, PDB-first defaults, and AlphaFold fallback.
- Bioinformatics tab remains a placeholder for future conservation/trait data.

Mol* behavior:

- Experimental PDB viewer highlights the selected mapped residue in red.
- Contact neighborhood viewer highlights the selected residue in red and indexed nearby/contact residues in blue.
- AlphaFold viewer highlights the selected site when UniProt and numeric site position are available.
- Viewers reload only when the selected site or relevant structure/contact control changes.

Important limitation: `pdb_interactions.parquet` stores nearest-residue fields per contact class, not a full list of every residue within 5A. The app can highlight indexed contact-nearest residues, but true all-residues-within-5A highlighting requires a richer coordinate-derived asset.

## Verification Already Completed

Static app worktree checks that passed:

```powershell
node --check assets/app.js
python -m py_compile build_static_data.py build_raw_hover_data.py
```

Latest local-review checks that passed after the browser diff comments:

```powershell
node tests\local-review-static.test.mjs
node --check output\local-review\assets\app.js
node --check C:\Users\ethan\.codex\worktrees\49ce\yptlib\assets\app.js
```

HTTP checks returned `200` for:

```text
http://127.0.0.1:8765/index.html?verify=raw-hover-v4
http://127.0.0.1:8765/assets/data/manifest.json
```

In-app browser verification on `http://127.0.0.1:8765/index.html?review=raw-hover-v4` confirmed:

- Summary loaded with `One-Shot (Complete): 31,821 sites, 742 compounds`.
- `Reactive Site Distribution` and `Reactive Compound Selectivity` are visible; old titles are absent.
- Dataset switch showed `Loading Fractionated data...` and then resolved to `Fractionated (Incomplete): 58,181 sites, 126 compounds`.
- No browser warnings/errors were reported during the verification pass.

Python Playwright smoke checks verified:

- Page loads.
- Summary target options include `Sites at PPI interface`.
- Summary site dropdown deduplicates labels.
- Site gene labels include `Sites Seen` and `Max R`.
- Strict Site filter combinations only list genes with matching active sites.
- Contact-bearing examples default to PDB structures and show only contact entries with indexed contact distances.
- Raw hover coverage was 100% for tested one-shot/fractionated filter combinations.
- `Z4262428218` with significant only, hide high variance, and Min SN 10 had 296 filtered Compound View points and zero missing raw-hover rows.

Deployment-fix checks that passed:

```powershell
node tests/deploy-config.test.mjs
node tests/summary-data-smoke.test.mjs
```

Live browser verification on `pages-data-v5`:

- Initial Summary Target Engagement rendered 20 bars.
- Custom gene `STAT2` rendered 3 site options and 20 bars.
- Built-in target checks:
  - `cancer`: 2000 site options, 20 bars.
  - `PPI`: 2000 site options, 20 bars.
  - `Metal`: 116 site options, 14 bars.
- No visible `NetworkError` text.
- No current-version browser warnings/errors.

## User Preferences To Preserve

- Mol* structure visualization is a core feature and should not be removed.
- Prefer inline structure widgets, not popups or new windows.
- Keep filters compact, aligned, and clear about scope.
- Use human-readable contact and compound type labels.
- Keep molecule structure hover in the same hover card as plot data.
- Raw hover cards should show true DMSO/compound duplicate SN values; do not substitute mean SN in raw-intensity hovers.
- Use direct `structures/*.png` images rather than repeated base64 generation.
- Prefer Morandi-style colors where possible.
- Use matplotlib-Blues-like gradients for p-value encoding, with low p-values more intense.
- Avoid large annotation arrows; use small labels or lightweight lines.
- Preserve deterministic sorting and empty-result behavior.
