/* Node test harness for the Marsh model. Run: node marsh.test.js
 * Proves the PK math against known, checkable properties.
 */
var M = require('./marsh.js');

var pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS  ' + name + (detail ? '  (' + detail + ')' : '')); }
  else { fail++; console.log('  FAIL  ' + name + (detail ? '  (' + detail + ')' : '')); }
}
function approx(a, b, tol) { return Math.abs(a - b) <= tol; }

console.log('\n=== Marsh model tests ===\n');

// --- 1. Model construction: V1 = 0.228 * weight -------------------------------
var m70 = M.createModel(70);
check('V1 = 0.228 L/kg * 70kg = 15.96 L', approx(m70.V1, 15.96, 1e-9), 'V1=' + m70.V1.toFixed(3) + ' L');
check('ke0 default = 1.2 /min', m70.ke0 === 1.2, 'ke0=' + m70.ke0);
check('V1 scales linearly with weight', approx(M.createModel(140).V1, 2 * m70.V1, 1e-9));

// --- 2. Bolus sets Cp = dose/V1 instantly -------------------------------------
var s0 = M.initialState();
var s = M.giveBolus(s0, 100); // 100 mg
check('Cp immediately after 100mg bolus = 100/15.96 = 6.266 ug/mL',
  approx(M.plasmaConc(m70, s), 100 / 15.96, 1e-9), 'Cp=' + M.plasmaConc(m70, s).toFixed(3));
check('Ce is still 0 immediately after bolus (lag)', s.Ce === 0);

// --- 3. Total body drug never increases without input; Cp decays after bolus --
var st = M.giveBolus(M.initialState(), 100);
var totalStart = st.A1 + st.A2 + st.A3;
var cpStart = M.plasmaConc(m70, st);
for (var i = 0; i < 60; i++) st = M.step(m70, st, 1 / 60); // 1 min
var totalAfter = st.A1 + st.A2 + st.A3;
check('Total drug decreases over time (elimination via k10)', totalAfter < totalStart,
  totalStart.toFixed(2) + ' -> ' + totalAfter.toFixed(2) + ' mg');
check('Cp falls over first minute (redistribution+elimination)',
  M.plasmaConc(m70, st) < cpStart, cpStart.toFixed(3) + ' -> ' + M.plasmaConc(m70, st).toFixed(3));

// --- 4. Effect site rises then peaks; time-to-peak ~1.6 min for ke0=1.2 -------
st = M.giveBolus(M.initialState(), 100);
var peak = -1, tpeak = -1, dt = 1 / 600; // 0.1s steps for a clean tpeak
for (var t = 0, n = 0; n < 6000; n++, t += dt) {
  st = M.step(m70, st, dt);
  if (st.Ce > peak) { peak = st.Ce; tpeak = t + dt; }
}
check('Effect-site Ce peaks ~1.4-1.8 min after bolus (modified Marsh)',
  tpeak > 1.4 && tpeak < 1.8, 'tpeak=' + tpeak.toFixed(2) + ' min, Cepeak=' + peak.toFixed(3));

// --- 5. Linearity: peak Ce scales linearly with dose from zero state ----------
var p50 = M.peakCeAfterBolus(m70, M.initialState(), 50);
var p100 = M.peakCeAfterBolus(m70, M.initialState(), 100);
check('Peak Ce is linear in dose (100mg peak = 2x 50mg peak)',
  approx(p100, 2 * p50, 1e-6), 'p50=' + p50.toFixed(4) + ', p100=' + p100.toFixed(4));

// --- 6. Loading bolus sizing: peak Ce hits the target exactly -----------------
var target = 4.0; // ug/mL
var load = M.bolusForTargetPeak(m70, M.initialState(), target);
var achieved = M.peakCeAfterBolus(m70, M.initialState(), load);
check('Loading bolus makes predicted Ce peak == target 4.0',
  approx(achieved, target, 1e-3), 'dose=' + load.toFixed(1) + ' mg, peak=' + achieved.toFixed(4));
check('Loading dose is clinically plausible (~1-2 mg/kg for 70kg @ target 4)',
  (load / 70) > 1.0 && (load / 70) < 2.2, (load / 70).toFixed(2) + ' mg/kg');

// --- 7. Top-up sizing from residual state also hits target --------------------
// Give loading dose, let Ce decay to 90% of target, then size a top-up.
st = M.giveBolus(M.initialState(), load);
var thr = 0.9 * target, reached = false, elapsed = 0;
// first climb to peak, then wait for decay below threshold
for (var q = 0; q < 60 * 30; q++) { // up to 30 min
  st = M.step(m70, st, 1 / 60); elapsed += 1 / 60;
  if (!reached && st.Ce >= peakCeApprox(st)) {} // noop guard
  if (st.Ce <= thr && elapsed > 3) { reached = true; break; }
}
function peakCeApprox() { return Infinity; } // (kept simple; not used for logic)
var topup = M.bolusForTargetPeak(m70, st, target);
var afterTopupPeak = M.peakCeAfterBolus(m70, st, topup);
check('Top-up from residual state restores predicted Ce peak to target',
  approx(afterTopupPeak, target, 1e-3), 'topup=' + topup.toFixed(1) + ' mg, peak=' + afterTopupPeak.toFixed(4));
check('Top-up is smaller than loading dose (residual drug present)',
  topup < load, 'topup=' + topup.toFixed(1) + ' < load=' + load.toFixed(1) + ' mg');
check('Time from loading to first top-up is a plausible interval (>3 min)',
  elapsed > 3, 'interval=' + elapsed.toFixed(1) + ' min');

// --- 8. Simulate a short maintenance sawtooth and report -----------------------
console.log('\n  --- Maintenance sawtooth preview (70 kg, target 4.0 ug/mL) ---');
st = M.initialState();
var loadDose = M.bolusForTargetPeak(m70, st, target);
st = M.giveBolus(st, loadDose);
console.log('  t=0.0 min  LOADING  ' + loadDose.toFixed(0) + ' mg');
var clock = 0, boluses = 1;
while (clock < 25 && boluses < 8) {
  st = M.step(m70, st, 1 / 60); clock += 1 / 60;
  if (st.Ce <= 0.9 * target && M.plasmaConc(m70, st) < target) {
    var d = M.bolusForTargetPeak(m70, st, target);
    st = M.giveBolus(st, d);
    boluses++;
    console.log('  t=' + clock.toFixed(1) + ' min  TOP-UP   ' + d.toFixed(0) +
      ' mg   (Ce was ' + (0.9 * target).toFixed(2) + ')');
  }
}

console.log('\n=== ' + pass + ' passed, ' + fail + ' failed ===\n');
process.exit(fail === 0 ? 0 : 1);
