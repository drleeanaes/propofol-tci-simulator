# Propofol TCI Simulator — Marsh (effect-site)

**Live:** https://drleeanaes.github.io/propofol-tci-simulator/

An educational, single-file web app that simulates propofol effect-site
concentration using the **modified Marsh** pharmacokinetic model (3 compartments
+ effect site, ke0 = 1.2 min⁻¹). It computes a loading bolus to reach a target
effect-site concentration (Ce), then prompts top-up boluses to maintain it, and
plots Ce and plasma Cp live over time.

> ⚠️ **Educational simulation only.** This reproduces published PK math for
> teaching. It is **not** a medical device and must **not** be used for real
> patient care or dosing decisions.

## Run / host

It's a single self-contained `index.html` — no build, no dependencies.

- **Locally:** open `index.html` in a browser, or serve the folder:
  `python3 -m http.server 8777` → http://localhost:8777
- **Host it:** drop `index.html` on any static host (GitHub Pages, Netlify,
  Cloudflare Pages, S3, etc.). Nothing server-side is required.

This repo auto-deploys to GitHub Pages via `.github/workflows/deploy.yml` on every
push to `main`. To redeploy, just push; to deploy manually, run the "Deploy static
site to GitHub Pages" workflow from the Actions tab.

## How it works

1. Enter **age**, **body weight (kg)**, and a **target Ce (µg/mL)**, then Start.
   (Marsh does not use age — it is recorded only; all volumes/clearances scale
   with weight.)
2. The app prompts a **loading bolus** sized so predicted Ce *peaks* exactly at
   the target, then asks you to administer it.
3. As Ce decays to **90% of target**, the app prompts a **top-up bolus** sized to
   bring Ce back to target. You administer each prompted dose (nothing is given
   automatically) — producing the characteristic maintenance sawtooth.
4. **Give bolus now** lets you inject any dose (mg) at any time; the model
   responds immediately.
5. A live chart plots **Ce** (effect site), **Cp** (plasma), the **target** line,
   and bolus markers. Speed control (1×–10×) fast-forwards the simulation.

## Marsh parameters (adult)

| Parameter | Value |
|---|---|
| V1 | 0.228 L/kg × weight |
| k10 | 0.119 min⁻¹ |
| k12 / k21 | 0.112 / 0.055 min⁻¹ |
| k13 / k31 | 0.042 / 0.0033 min⁻¹ |
| ke0 (modified Marsh) | 1.2 min⁻¹ |

Effect site: `dCe/dt = ke0 · (Cp − Ce)`, with `Cp = A1 / V1`.

## Files

- `index.html` — the app (Marsh model inlined; single self-contained file).
- `marsh.js` — the Marsh model as a standalone module (source of truth; the same
  code is inlined into `index.html`).
- `marsh.test.js` — Node test harness proving the PK math. Run: `node marsh.test.js`.

## Development notes

`index.html` exposes `window.__app` and `window.__test` (with `advance(minutes)`,
`administer()`, `manual()`, `snapshot()`) for deterministic testing/inspection.
These are harmless debug hooks and do not affect the UI.
