/*
 * Marsh propofol pharmacokinetic model (modified, effect-site).
 *
 * 3-compartment PK + effect-site compartment.
 * - Weight-proportional (adult Marsh does NOT use age).
 * - Effect site: dCe/dt = ke0 * (Cp - Ce), ke0 = 1.2/min ("modified Marsh").
 *
 * Units:
 *   weight  : kg
 *   A1,A2,A3: mg (drug amount in each compartment)
 *   V1      : L
 *   Cp, Ce  : mg/L  (numerically equal to microg/mL for propofol)
 *   time    : minutes internally
 *
 * Educational simulation of published PK math. NOT for clinical patient care.
 */
(function (global) {
  'use strict';

  var MARSH_RATES = {
    k10: 0.119,   // /min  elimination from central
    k12: 0.112,   // /min  central -> fast peripheral
    k13: 0.042,   // /min  central -> slow peripheral
    k21: 0.055,   // /min  fast peripheral -> central
    k31: 0.0033,  // /min  slow peripheral -> central
    V1_per_kg: 0.228 // L/kg central volume
  };

  var DEFAULT_KE0 = 1.2; // /min modified Marsh effect-site rate

  function createModel(weightKg, options) {
    options = options || {};
    var ke0 = (options.ke0 != null) ? options.ke0 : DEFAULT_KE0;
    return {
      weight: weightKg,
      V1: MARSH_RATES.V1_per_kg * weightKg,
      k10: MARSH_RATES.k10,
      k12: MARSH_RATES.k12,
      k13: MARSH_RATES.k13,
      k21: MARSH_RATES.k21,
      k31: MARSH_RATES.k31,
      ke0: ke0
    };
  }

  function initialState() {
    return { A1: 0, A2: 0, A3: 0, Ce: 0 };
  }

  function cloneState(s) {
    return { A1: s.A1, A2: s.A2, A3: s.A3, Ce: s.Ce };
  }

  function plasmaConc(m, s) {
    return s.A1 / m.V1;
  }

  // Time derivatives of the state vector. infusionMgPerMin optional (continuous input).
  function derivs(m, s, infusionMgPerMin) {
    var inf = infusionMgPerMin || 0;
    var Cp = s.A1 / m.V1;
    return {
      A1: -(m.k10 + m.k12 + m.k13) * s.A1 + m.k21 * s.A2 + m.k31 * s.A3 + inf,
      A2: m.k12 * s.A1 - m.k21 * s.A2,
      A3: m.k13 * s.A1 - m.k31 * s.A3,
      Ce: m.ke0 * (Cp - s.Ce)
    };
  }

  // One RK4 integration step. dtMin in minutes. Returns new state (does not mutate).
  function step(m, s, dtMin, infusionMgPerMin) {
    var k1 = derivs(m, s, infusionMgPerMin);
    var s2 = _advance(s, k1, 0.5 * dtMin);
    var k2 = derivs(m, s2, infusionMgPerMin);
    var s3 = _advance(s, k2, 0.5 * dtMin);
    var k3 = derivs(m, s3, infusionMgPerMin);
    var s4 = _advance(s, k3, dtMin);
    var k4 = derivs(m, s4, infusionMgPerMin);
    return {
      A1: s.A1 + dtMin / 6 * (k1.A1 + 2 * k2.A1 + 2 * k3.A1 + k4.A1),
      A2: s.A2 + dtMin / 6 * (k1.A2 + 2 * k2.A2 + 2 * k3.A2 + k4.A2),
      A3: s.A3 + dtMin / 6 * (k1.A3 + 2 * k2.A3 + 2 * k3.A3 + k4.A3),
      Ce: s.Ce + dtMin / 6 * (k1.Ce + 2 * k2.Ce + 2 * k3.Ce + k4.Ce)
    };
  }

  function _advance(s, k, h) {
    return {
      A1: s.A1 + h * k.A1,
      A2: s.A2 + h * k.A2,
      A3: s.A3 + h * k.A3,
      Ce: s.Ce + h * k.Ce
    };
  }

  // Inject an instantaneous bolus (mg) into the central compartment. Returns new state.
  function giveBolus(s, doseMg) {
    var ns = cloneState(s);
    ns.A1 += doseMg;
    return ns;
  }

  // Predicted peak Ce reached after adding a bolus of doseMg to the current state,
  // simulated over horizonMin. Used to size boluses to a target effect-site conc.
  function peakCeAfterBolus(m, s, doseMg, horizonMin, dtMin) {
    horizonMin = horizonMin || 8;
    dtMin = dtMin || (1 / 60);
    var st = giveBolus(s, doseMg);
    var peak = st.Ce;
    var n = Math.round(horizonMin / dtMin);
    for (var i = 0; i < n; i++) {
      st = step(m, st, dtMin);
      if (st.Ce > peak) peak = st.Ce;
    }
    return peak;
  }

  // Smallest bolus (mg) whose predicted peak Ce reaches targetCe, given current state.
  // Binary search (peak Ce is monotonic in dose). Works for the loading dose (zero
  // state) and for top-ups (residual state).
  function bolusForTargetPeak(m, s, targetCe, horizonMin, dtMin) {
    if (peakCeAfterBolus(m, s, 0, horizonMin, dtMin) >= targetCe) {
      return 0; // already at/above target without more drug
    }
    var lo = 0, hi = 1;
    while (peakCeAfterBolus(m, s, hi, horizonMin, dtMin) < targetCe && hi < 1e6) {
      hi *= 2;
    }
    for (var i = 0; i < 50; i++) {
      var mid = 0.5 * (lo + hi);
      if (peakCeAfterBolus(m, s, mid, horizonMin, dtMin) < targetCe) lo = mid;
      else hi = mid;
    }
    return 0.5 * (lo + hi);
  }

  var api = {
    MARSH_RATES: MARSH_RATES,
    DEFAULT_KE0: DEFAULT_KE0,
    createModel: createModel,
    initialState: initialState,
    cloneState: cloneState,
    plasmaConc: plasmaConc,
    derivs: derivs,
    step: step,
    giveBolus: giveBolus,
    peakCeAfterBolus: peakCeAfterBolus,
    bolusForTargetPeak: bolusForTargetPeak
  };

  if (typeof module !== 'undefined' && module.exports) {
    module.exports = api;
  } else {
    global.Marsh = api;
  }
})(typeof globalThis !== 'undefined' ? globalThis : this);
