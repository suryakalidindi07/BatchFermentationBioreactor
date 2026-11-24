(function (root) {
'use strict';

const PRESET = Object.freeze({
  organism: 'E. coli',
  mode: 'Batch',
  muOptPerHour: 0.60,
  ksGramsPerLiter: 0.050,
  yieldBiomassPerSubstrate: 0.45,
  substrateInitialGramsPerLiter: 20.0,
  biomassInitialGramsPerLiter: 0.05,
  volumeLiters: 20.0,
  substrateName: 'Glucose',
  productName: 'Acetate (illustrative estimate)',
  productAlphaGramsPerGramBiomass: 0.10,
  productBetaGramsPerGramBiomassHour: 0.0,
  temperatureMinimumC: 3.06,
  temperatureOptimumC: 41.10,
  temperatureMaximumC: 45.06,
  phMinimum: 3.88,
  phOptimum: 7.20,
  phMaximum: 12.17,
  oxygenMonodKmmolPerLiter: 0.010,
  oxygenMmolPerGramBiomass: 8.0,
  refRpm: 480,
  refVvm: 1.0,
  refPowerPerVolumeWm3: 2521,
  refOxygenSaturationMgPerLiter: 7.0,
  refTemperatureC: 37.0,
  tankHeightToDiameter: 1.16,
  vanTRietCoefficient: 0.026,
  vanTRietPowerExponent: 0.4,
  vanTRietGasExponent: 0.5,
  defaultDurationHours: 24.0,
  integrationStepHours: 0.002,
  sampleEveryHours: 0.1
});

const DEFAULT_CONTROLS = Object.freeze({
  rpm: 480, vvm: 1.0, temperatureC: 37, ph: 7.0
});

const INPUT_RANGES = Object.freeze({
  rpm: Object.freeze([100, 1400]),
  vvm: Object.freeze([0, 2.5]),
  temperatureC: Object.freeze([15, 50]),
  ph: Object.freeze([2, 11])
});

function validateControls(input = DEFAULT_CONTROLS) {
  const values = {};
  for (const key of Object.keys(INPUT_RANGES)) {
    const n = Number(input[key]);
    const [min, max] = INPUT_RANGES[key];
    if (!Number.isFinite(n) || n < min || n > max) {
      throw new RangeError(`${key} must be a finite number from ${min} to ${max}.`);
    }
    values[key] = n;
  }
  return values;
}

function clamp(value, min, max) {
  return Math.min(max, Math.max(min, value));
}

function temperatureFactor(tempC) {
  const a = PRESET.temperatureMinimumC;
  const b = PRESET.temperatureOptimumC;
  const c = PRESET.temperatureMaximumC;
  if (tempC <= a || tempC >= c) return 0;
  const numer = (tempC - c) * (tempC - a) ** 2;
  const denom = (b - a) * ((b - a) * (tempC - b) -
    (b - c) * (b + a - 2 * tempC));
  return clamp(numer / denom, 0, 1);
}

function phFactor(ph) {
  const a = PRESET.phMinimum;
  const b = PRESET.phOptimum;
  const c = PRESET.phMaximum;
  if (ph <= a || ph >= c) return 0;
  const p = (ph - a) * (ph - c);
  return clamp(p / (p - (ph - b) ** 2), 0, 1);
}

function oxygenSaturationMgPerLiter(tempC) {
  return PRESET.refOxygenSaturationMgPerLiter *
    Math.exp(-0.022 * (tempC - PRESET.refTemperatureC));
}

function tankDiameterMeters() {
  const volumeM3 = PRESET.volumeLiters / 1000;
  return Math.cbrt((4 * volumeM3) / (Math.PI * PRESET.tankHeightToDiameter));
}

function transferConditions(controls) {
  const { rpm, vvm, temperatureC } = validateControls(controls);
  const d = tankDiameterMeters();
  const areaM2 = Math.PI * d * d / 4;
  const gasFlowM3s = (vvm * (PRESET.volumeLiters / 1000)) / 60;
  const superficialVelocityMs = gasFlowM3s / areaM2;
  const powerPerVolumeWm3 = PRESET.refPowerPerVolumeWm3 *
    (rpm / PRESET.refRpm) ** 3;
  const klaPerHour = 3600 * PRESET.vanTRietCoefficient *
    powerPerVolumeWm3 ** PRESET.vanTRietPowerExponent *
    superficialVelocityMs ** PRESET.vanTRietGasExponent;
  const oxygenSatMgL = oxygenSaturationMgPerLiter(temperatureC);
  const oxygenSatMmolL = oxygenSatMgL / 32;
  return {
    klaPerHour,
    oxygenSatMgL,
    oxygenSatMmolL,
    maximumOtrMmolLh: klaPerHour * oxygenSatMmolL,
    powerPerVolumeWm3,
    superficialVelocityMs,
    tankDiameterM: d
  };
}

function evaluatedRates(state, muEnvironment, conditions) {
  const x = Math.max(0, state.x);
  const s = Math.max(0, state.s);
  const c = clamp(state.c, 0, conditions.oxygenSatMmolL);
  const substrateFactor = s <= 0 ? 0 : s / (PRESET.ksGramsPerLiter + s);
  const oxygenFactor = c <= 0 ? 0 :
    c / (PRESET.oxygenMonodKmmolPerLiter + c);
  const muPerHour = muEnvironment * substrateFactor * oxygenFactor;
  const biomassRate = muPerHour * x;
  const uptake = PRESET.oxygenMmolPerGramBiomass * biomassRate;
  const transfer = conditions.klaPerHour * (conditions.oxygenSatMmolL - c);
  return {
    dx: biomassRate,
    ds: -biomassRate / PRESET.yieldBiomassPerSubstrate,
    dc: transfer - uptake,
    muPerHour,
    uptakeMmolLh: uptake,
    transferMmolLh: transfer
  };
}

function rk4Step(state, dt, muEnvironment, conditions, maxBiomass) {
  const rates = (st) => evaluatedRates(st, muEnvironment, conditions);
  const increment = (st, k, multiplier) => ({
    x: st.x + multiplier * k.dx,
    s: st.s + multiplier * k.ds,
    c: st.c + multiplier * k.dc
  });
  const k1 = rates(state);
  const k2 = rates(increment(state, k1, dt / 2));
  const k3 = rates(increment(state, k2, dt / 2));
  const k4 = rates(increment(state, k3, dt));
  const nextX = clamp(state.x + dt / 6 *
    (k1.dx + 2 * k2.dx + 2 * k3.dx + k4.dx), state.x, maxBiomass);
  return {
    x: nextX,
    s: Math.max(0, (maxBiomass - nextX) / PRESET.yieldBiomassPerSubstrate),
    c: clamp(state.c + dt / 6 *
      (k1.dc + 2 * k2.dc + 2 * k3.dc + k4.dc),
      0, conditions.oxygenSatMmolL)
  };
}


function simulate(input = DEFAULT_CONTROLS, durationHours = PRESET.defaultDurationHours) {
  const controls = validateControls(input);
  if (!(Number.isFinite(durationHours) && durationHours > 0 && durationHours <= 168)) {
    throw new RangeError('durationHours must be within (0,168].');
  }
  const conditions = transferConditions(controls);
  const gammaTemperature = temperatureFactor(controls.temperatureC);
  const gammaPh = phFactor(controls.ph);
  const muEnvironment = PRESET.muOptPerHour * gammaTemperature * gammaPh;
  const maxBiomass = PRESET.biomassInitialGramsPerLiter +
    PRESET.yieldBiomassPerSubstrate * PRESET.substrateInitialGramsPerLiter;
  const steps = Math.ceil(durationHours / PRESET.integrationStepHours);
  const dt = durationHours / steps;
  const sampleEverySteps = Math.max(1, Math.round(PRESET.sampleEveryHours / dt));
  let state = {
    x: PRESET.biomassInitialGramsPerLiter,
    s: PRESET.substrateInitialGramsPerLiter,
    c: conditions.oxygenSatMmolL
  };
  const samples = [];
  function addSample(hour) {
    const rates = evaluatedRates(state, muEnvironment, conditions);
    const productGramsPerLiter = PRESET.productAlphaGramsPerGramBiomass *
      Math.max(0, state.x - PRESET.biomassInitialGramsPerLiter);
    samples.push({
      productGramsPerLiter,
      productGrams: productGramsPerLiter * PRESET.volumeLiters,
      productFormationRateGramsPerLiterHour: PRESET.productAlphaGramsPerGramBiomass * rates.dx,
      timeHours: hour,
      biomassGramsPerLiter: state.x,
      substrateGramsPerLiter: state.s,
      dissolvedOxygenMgPerLiter: state.c * 32,
      dissolvedOxygenPercent: conditions.oxygenSatMmolL > 0 ?
        100 * state.c / conditions.oxygenSatMmolL : 0,
      muPerHour: rates.muPerHour,
      otrMmolLh: rates.transferMmolLh,
      ourMmolLh: rates.uptakeMmolLh
    });
  }
  addSample(0);
  for (let i = 1; i <= steps; i++) {
    state = rk4Step(state, dt, muEnvironment, conditions, maxBiomass);
    if (i % sampleEverySteps === 0 || i === steps) {
      addSample(i * dt);
    }
  }
  const final = samples[samples.length - 1];
  let firstDoBelow20Hours = null;
  for (const sample of samples) {
    if (sample.dissolvedOxygenPercent < 20) {
      firstDoBelow20Hours = sample.timeHours;
      break;
    }
  }
  return {
    controls,
    preset: PRESET,
    conditions: {
      ...conditions, gammaTemperature, gammaPh, muEnvironmentPerHour: muEnvironment
    },
    samples,
    summary: {
      finalProductGrams: final.productGrams,
      finalProductGramsPerLiter: final.productGramsPerLiter,
      averageProductivityGramsPerLiterHour: final.productGramsPerLiter / durationHours,
      finalBiomassGramsPerLiter: final.biomassGramsPerLiter,
      finalSubstrateGramsPerLiter: final.substrateGramsPerLiter,
      finalOxygenPercent: final.dissolvedOxygenPercent,
      maximumBiomassFromSubstrateGramsPerLiter: maxBiomass,
      minimumOxygenPercent: Math.min(...samples.map((row) => row.dissolvedOxygenPercent)),
      firstDoBelow20Hours,
      durationHours
    }
  };
}

function csvFromResult(result) {
  const fields = [
    'timeHours', 'productGrams', 'productGramsPerLiter',
    'productFormationRateGramsPerLiterHour', 'biomassGramsPerLiter', 'substrateGramsPerLiter',
    'dissolvedOxygenMgPerLiter', 'dissolvedOxygenPercent', 'muPerHour',
    'otrMmolLh', 'ourMmolLh'
  ];
  const inputFields = ['rpm', 'vvm', 'temperatureC', 'ph'];
  const header = [...inputFields, ...fields].join(',');
  const rows = result.samples.map((s) => [
    ...inputFields.map((k) => Number(result.controls[k]).toFixed(7)),
    ...fields.map((k) => Number(s[k]).toFixed(7))
  ].join(','));
  return [header, ...rows].join('\n') + '\n';
}

const EcoliBatchModel = Object.freeze({
  PRESET,
  DEFAULT_CONTROLS,
  INPUT_RANGES,
  validateControls,
  temperatureFactor,
  phFactor,
  oxygenSaturationMgPerLiter,
  tankDiameterMeters,
  transferConditions,
  simulate,
  csvFromResult
});
root.EcoliBatchModel = EcoliBatchModel;
if (typeof module === 'object' && module && module.exports) {
  module.exports = EcoliBatchModel;
}
})(typeof globalThis !== 'undefined' ? globalThis : this);
