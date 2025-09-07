function extractGaugeDataFromCSV(csvText) {
  // Split lines, skip header
  const lines = csvText.trim().split(/\r?\n/).slice(1);
  // For each line, extract columns 1-3 (index 0,1,2)
  return lines.map(line => {
    const cols = line.split(',');
    return [parseFloat(cols[0]), parseFloat(cols[1]), parseFloat(cols[2])];
  }).filter(row => row.every(v => !isNaN(v)));
}
// Wave App - Full Functionality (migrated from previous website)

// --- Complex algebra helpers for reflection algorithm ---
function wavelen(h, T) {
  // Accurate wavelength calculation using iterative solution for dispersion relation
  // h: water depth, T: wave period
  const g = 9.81;
  const omega = 2 * Math.PI / T;
  let L = (g * T * T) / (2 * Math.PI); // deep water guess
  for (let iter = 0; iter < 20; iter++) {
    const k = 2 * Math.PI / L;
    L = (g * T * T) / (2 * Math.PI) * Math.tanh(k * h);
  }
  return L;
}

function zeroCrossing(data, dt) {
  // Accurate zero-crossing wave statistics
  // Remove mean
  const mean = data.reduce((a, b) => a + b, 0) / data.length;
  const detrended = data.map(v => v - mean);
  // Find upward zero crossings
  let crossings = [];
  for (let i = 1; i < detrended.length; i++) {
    if (detrended[i - 1] < 0 && detrended[i] >= 0) {
      crossings.push(i);
    }
  }
  // Remove first crossing if data starts above zero
  if (detrended[0] > 0 && crossings.length) crossings = crossings.slice(1);
  // Only keep upward crossings
  if (crossings.length < 2) {
    return { nWaves: 0, Hs: NaN, Hmean: NaN, H10: NaN, Tmean: NaN, Ts: NaN, waves: [] };
  }
  // Calculate crest, trough, period for each wave
  let wave = [];
  for (let i = 0; i < crossings.length - 1; i++) {
    let seg = detrended.slice(crossings[i], crossings[i + 1]);
    let crest = Math.max(...seg);
    let trough = -Math.min(...seg);
    let period = (crossings[i + 1] - crossings[i]) * dt;
    wave.push({ height: crest + trough, crest, trough, period });
  }
  // Threshold: 1% of max wave height
  let threshold = 0.01 * Math.max(...wave.map(w => w.height));
  // Merge/remove small waves
  let i = 0;
  while (i < wave.length) {
    if (wave[i].crest < threshold) {
      if (i > 0) {
        wave[i - 1].crest = Math.max(wave[i - 1].crest, wave[i].crest);
        wave[i - 1].trough = Math.max(wave[i - 1].trough, wave[i].trough);
        wave[i - 1].period += wave[i].period;
        wave.splice(i, 1);
        i--;
      } else {
        wave.splice(i, 1);
        i--;
      }
    } else if (wave[i].trough < threshold) {
      if (i < wave.length - 1) {
        wave[i].crest = Math.max(wave[i].crest, wave[i + 1].crest);
        wave[i].trough = Math.max(wave[i].trough, wave[i + 1].trough);
        wave[i].period += wave[i + 1].period;
        wave.splice(i + 1, 1);
      } else {
        wave.splice(i, 1);
        i--;
      }
    }
    i++;
  }
  // Compute total wave height
  wave.forEach(w => { w.height = w.crest + w.trough; });
  let nb = wave.length;
  if (nb === 0) {
    return { nWaves: 0, Hs: NaN, Hmean: NaN, H10: NaN, Tmean: NaN, Ts: NaN, waves: [] };
  }
  // Sort by height descending
  let sorted = wave.slice().sort((a, b) => b.height - a.height);
  // Statistics
  let Hs = sorted.slice(0, Math.round(nb / 3)).reduce((a, w) => a + w.height, 0) / Math.round(nb / 3);
  let Hmean = sorted.reduce((a, w) => a + w.height, 0) / nb;
  let H10 = sorted.slice(0, Math.max(1, Math.round(nb / 10))).reduce((a, w) => a + w.height, 0) / Math.max(1, Math.round(nb / 10));
  let Hmax = sorted[0].height;
  let Tmean = sorted.reduce((a, w) => a + w.period, 0) / nb;
  let Ts = sorted.slice(0, Math.round(nb / 3)).reduce((a, w) => a + w.period, 0) / Math.round(nb / 3);
  return {
    nWaves: nb,
    Hs,
    Hmean,
    H10,
    Hmax,
    Tmean,
    Ts,
    waves: wave.map(w => ({ height: w.height, period: w.period }))
  };
}
function threeArray(h, dt, z, gpos) {
  // Direct conversion from MATLAB three_array function
  // Calculates the incident and reflected waves for random wave experiments
  // conducted in laboratory flumes using three gauge array.
  // Based on Goda and Suzuki (1976) expanded to three gauges by Kobayashi and Cox
  
  console.log('threeArray called with:');
  console.log('h (depth):', h);
  console.log('dt (time step):', dt);
  console.log('z type:', typeof z, 'is array:', Array.isArray(z));
  
  // Handle different input formats
  let gaugeData;
  if (Array.isArray(z) && Array.isArray(z[0])) {
    // z is already in format [[gauge1_data], [gauge2_data], [gauge3_data]]
    gaugeData = z;
    console.log('z data lengths:', z.map(arr => arr.length));
  } else if (Array.isArray(z)) {
    // z might be in MATLAB format (array of rows, each row has 3 columns)
    // Convert to [[gauge1], [gauge2], [gauge3]] format
    gaugeData = [
      z.map(row => row[0]),
      z.map(row => row[1]), 
      z.map(row => row[2])
    ];
    console.log('Converted z data lengths:', gaugeData.map(arr => arr.length));
  } else {
    throw new Error('Invalid data format for z parameter');
  }
  
  console.log('gpos (positions):', gpos);
  
  // Detrend data (remove mean)
  const zDetrended = gaugeData.map(gauge => {
    const mean = gauge.reduce((a, b) => a + b, 0) / gauge.length;
    return gauge.map(val => val - mean);
  });
  
  const g = 9.81;
  const nfft = gaugeData[0].length;  // Use actual length, not next power of 2
  const df = 1 / (nfft * dt);  // frequency resolution
  const half = Math.round(nfft / 2);
  
  console.log('nfft:', nfft, 'df:', df, 'half:', half);
  
  // Preallocation
  const An = Array(half - 1).fill(0).map(() => Array(3).fill(0));
  const Bn = Array(half - 1).fill(0).map(() => Array(3).fill(0));
  const Sn = Array(half - 1).fill(0).map(() => Array(3).fill(0));
  const k = Array(half - 1).fill(0);
  const Ainc = Array(half - 1).fill(0).map(() => Array(3).fill(0));
  const Binc = Array(half - 1).fill(0).map(() => Array(3).fill(0));
  const Aref = Array(half - 1).fill(0).map(() => Array(3).fill(0));
  const Bref = Array(half - 1).fill(0).map(() => Array(3).fill(0));
  const nmin = Array(3).fill(0);
  const nmax = Array(3).fill(0);
  
  // Solving Fourier coefficients for all three gauges
  for (let j = 0; j < 3; j++) {
    // Apply Hann window (equivalent to MATLAB hann function)
    const win = Array(nfft).fill(0).map((_, i) => 0.5 * (1 - Math.cos(2 * Math.PI * i / (nfft - 1))));
    const windowed = zDetrended[j].map((val, i) => val * win[i]);
    
    const fn = FFT.fft(windowed);
    
    // Extract Fourier coefficients
    for (let i = 0; i < half - 1; i++) {
      const idx = i + 1; // MATLAB index 2:half corresponds to JS index 1:half-1
      An[i][j] = 2 * fn[idx].re / nfft;  // Real component
      Bn[i][j] = -2 * fn[idx].im / nfft; // Imaginary component
      
      // Spectral density
      const fnSquared = fn[idx].re * fn[idx].re + fn[idx].im * fn[idx].im;
      Sn[i][j] = dt * fnSquared * 2 / nfft;
    }
  }
  
  // Frequency array
  const f = Array(An.length).fill(0).map((_, i) => df * (i + 1));
  
  console.log('Frequency array length:', f.length, 'sample f[0:5]:', f.slice(0, 5));
  
  // Calculate wavenumber at each frequency
  for (let j = 0; j < f.length; j++) {
    const L = wavelen(h, 1 / f[j]);
    k[j] = (2 * Math.PI) / L;
  }
  
  // Three Gauge Array setup (MATLAB 1-based to JS 0-based)
  const g1 = [0, 0, 1]; // [1 1 2] in MATLAB becomes [0 0 1] in JS
  const g2 = [1, 2, 2]; // [2 3 3] in MATLAB becomes [1 2 2] in JS
  
  // Solve for incident and reflected wave amplitudes
  for (let j = 0; j < 3; j++) {
    const A1 = An.map(row => row[g1[j]]);
    const A2 = An.map(row => row[g2[j]]);
    const B1 = Bn.map(row => row[g1[j]]);
    const B2 = Bn.map(row => row[g2[j]]);
    const pos1 = gpos[g1[j]];
    const pos2 = gpos[g2[j]];
    
    for (let i = 0; i < k.length; i++) {
      const term1 = -A2[i] * Math.sin(k[i] * pos1) + A1[i] * Math.sin(k[i] * pos2) + 
                    B2[i] * Math.cos(k[i] * pos1) - B1[i] * Math.cos(k[i] * pos2);
      const term2 = A2[i] * Math.cos(k[i] * pos1) - A1[i] * Math.cos(k[i] * pos2) + 
                    B2[i] * Math.sin(k[i] * pos1) - B1[i] * Math.sin(k[i] * pos2);
      const term3 = -A2[i] * Math.sin(k[i] * pos1) + A1[i] * Math.sin(k[i] * pos2) - 
                    B2[i] * Math.cos(k[i] * pos1) + B1[i] * Math.cos(k[i] * pos2);
      const term4 = A2[i] * Math.cos(k[i] * pos1) - A1[i] * Math.cos(k[i] * pos2) - 
                    B2[i] * Math.sin(k[i] * pos1) + B1[i] * Math.sin(k[i] * pos2);
      
      const denominator = 2 * Math.sin(k[i] * Math.abs(pos2 - pos1));
      
      Ainc[i][j] = term1 / denominator;
      Binc[i][j] = term2 / denominator;
      Aref[i][j] = term3 / denominator;
      Bref[i][j] = term4 / denominator;
    }
    
    // Upper and lower limits of significant spectra
    const Lmin = Math.abs((pos2 - pos1) / 0.45);
    const Lmax = Math.abs((pos2 - pos1) / 0.05);
    const kmax = 2 * Math.PI / Lmin;
    const kmin = 2 * Math.PI / Lmax;
    
    const validIndices = k.map((kval, idx) => (kval <= kmax && kval >= kmin) ? idx : -1)
                          .filter(idx => idx !== -1);
    
    nmin[j] = Math.min(...validIndices);
    nmax[j] = Math.max(...validIndices);
  }
  
  const rangeStart = Math.min(...nmin);
  const rangeEnd = Math.max(...nmax);
  const range = [];
  for (let i = rangeStart; i <= rangeEnd; i++) {
    range.push(i);
  }
  
  console.log('Valid range:', range.length, 'frequencies from', rangeStart, 'to', rangeEnd);
  
  // Set invalid values to NaN
  for (let j = 0; j < 3; j++) {
    for (let i = 0; i < Ainc.length; i++) {
      if (i < nmin[j] || i > nmax[j]) {
        Ainc[i][j] = NaN;
        Binc[i][j] = NaN;
        Aref[i][j] = NaN;
        Bref[i][j] = NaN;
      }
    }
  }
  
  // Helper function for nanmean
  const nanmean = (arr) => {
    const validValues = arr.filter(val => !isNaN(val) && isFinite(val));
    return validValues.length > 0 ? validValues.reduce((a, b) => a + b, 0) / validValues.length : NaN;
  };
  
  // Average overlapped amplitudes across gauges
  const Aincav = Array(Ainc.length).fill(NaN);
  const Bincav = Array(Binc.length).fill(NaN);
  const Arefav = Array(Aref.length).fill(NaN);
  const Brefav = Array(Bref.length).fill(NaN);
  
  for (let i of range) {
    Aincav[i] = nanmean([Ainc[i][0], Ainc[i][1], Ainc[i][2]]);
    Bincav[i] = nanmean([Binc[i][0], Binc[i][1], Binc[i][2]]);
    Arefav[i] = nanmean([Aref[i][0], Aref[i][1], Aref[i][2]]);
    Brefav[i] = nanmean([Bref[i][0], Bref[i][1], Bref[i][2]]);
  }
  
  // Calculate spectra
  const Si = range.map(i => (Aincav[i] * Aincav[i] + Bincav[i] * Bincav[i]) / (2 * df));
  const Sr = range.map(i => (Arefav[i] * Arefav[i] + Brefav[i] * Brefav[i]) / (2 * df));
  
  console.log('Calculated spectra lengths - Si:', Si.length, 'Sr:', Sr.length);
  console.log('Sample Si values:', Si.slice(0, 5));
  console.log('Sample Sr values:', Sr.slice(0, 5));
  
  // Evaluate energies
  const Ei = Si.reduce((sum, val) => sum + val, 0) * df;
  const Er = Sr.reduce((sum, val) => sum + val, 0) * df;
  
  // Reflection coefficient
  const refco = Math.sqrt(Er / Ei);
  
  // Calculate wave heights
  const mo = range.map(i => Sn[i][0]).reduce((sum, val) => sum + val, 0) * df;
  const Htot = 4.004 * Math.sqrt(mo);
  const Hi = 4.004 * Math.sqrt(Ei);
  const Hr = 4.004 * Math.sqrt(Er);
  
  console.log('Final results: Ei=', Ei, 'Er=', Er, 'refco=', refco);
  console.log('Wave heights: Hi=', Hi, 'Hr=', Hr, 'Htot=', Htot);
  
  // Return results matching the expected structure
  return {
    Hi: Hi,
    Hr: Hr,
    Kr: refco,
    validFrequencies: range.length,
    details: range.map((i, idx) => ({
      freq: f[i],
      Si: Si[idx],
      Sr: Sr[idx],
      Sf: Sn[i][0]  // Use first gauge composite spectrum
    }))
  };
}
function parseSingleColumnSeries(txt) {
  // Parse single column of numbers from text
  return txt.split(/\r?\n/)
    .map(line => parseFloat(line.trim()))
    .filter(v => !isNaN(v));
}
function parseThreeGaugeMatrix(txt) {
  // Parse 3-column CSV/tab/space delimited data
  return txt.split(/\r?\n/)
    .map(line => line.trim().split(/,|\s+/).map(Number))
    .filter(row => row.length >= 3 && row.every(v => !isNaN(v)));
}

// Example dataset generator
const REAL_WAVE_DATA={ depth:0.30, waveFrequency:0.75, samplingFrequency:100, gaugeSpacing:[0,0.3,0.9], generateData:function(duration=20){ const dt=1/this.samplingFrequency; const n=Math.floor(duration/dt); const omega=2*Math.PI*this.waveFrequency; const T=1/this.waveFrequency; const L=wavelen(this.depth,T); const k=2*Math.PI/L; const Ai=0.025; const Ar=0.008; const phase=0.2; const data=[]; for (let i=0;i<n;i++){ const t=i*dt; const row=[]; for (let j=0;j<3;j++){ const x=this.gaugeSpacing[j]; const incident=Ai*Math.sin(k*x - omega*t); const reflected=Ar*Math.sin(-k*x - omega*t + phase); const secondary=0.003*Math.sin(1.5*omega*t + 0.5*j); const elev=incident+reflected+secondary; row.push(elev.toFixed(6)); } data.push(row);} return data; }};
const CACHED_EXAMPLE_DATA = REAL_WAVE_DATA.generateData(20, 20250208);

// State
let currentWaveData=null, currentWaveResults=null, currentReflectionData=null;

function byId(id){ return document.getElementById(id); }

// Wavelength example & calculation
window.loadWavelengthExample=function(showAlert=true){
  byId('water-depth').value=REAL_WAVE_DATA.depth;
  byId('wave-period').value=(1/REAL_WAVE_DATA.waveFrequency).toFixed(2);
  if(showAlert) alert('Example loaded.');
};
window.calculateWavelength=function(){
  try {
    const h=parseFloat(byId('water-depth').value);
    const T=parseFloat(byId('wave-period').value);
    if(!(h>0)&&!(T>0)) return alert('Enter positive values.');
    const L=wavelen(h,T);
    const k=2*Math.PI/L;
    const c=L/T;
    const omega=2*Math.PI/T;
    const n=0.5*(1+2*k*h/Math.sinh(2*k*h));
    const cg=n*c;
    let regime=h/L>0.5?'Deep Water': h/L<0.05?'Shallow Water':'Intermediate Water';
    byId('wavelength-value').textContent=L.toFixed(3);
    byId('wavelength-details').innerHTML=`<strong>Additional Parameters:</strong><br>• Wave Number (k): ${k.toFixed(4)} rad/m<br>• Wave Celerity (c): ${c.toFixed(2)} m/s<br>• Group Velocity (cg): ${cg.toFixed(2)} m/s<br>• Relative Depth (h/L): ${(h/L).toFixed(3)}<br>• Water Regime: ${regime}<br>• Angular Frequency: ${omega.toFixed(3)} rad/s`;
    byId('wavelength-result').style.display='block';
  } catch(e){
    alert('Error calculating wavelength.');
  }
};

// Wave stats
window.loadWaveStatsExample=function(showAlert=true){
  const waveData=CACHED_EXAMPLE_DATA;
  const first=waveData.map(r=>parseFloat(r[0]));
  byId('wave-data').value=first.map(v=>v.toFixed(6)).join('\n');
  byId('sampling-frequency').value=REAL_WAVE_DATA.samplingFrequency;
  if(showAlert) alert('Wave stats example loaded.');
};
window.calculateWaveStats=function(){
  try {
    const txt=byId('wave-data').value.trim();
    const fs=parseFloat(byId('sampling-frequency').value);
    if(!txt||!(fs>0)) return alert('Enter data and valid frequency');
    const data=parseSingleColumnSeries(txt);
    if(!data.length){ alert('No numeric values parsed. Ensure you pasted only one column.'); return; }
    if(data.length<10) return alert('Need at least 10 points (parsed '+data.length+')');
    const dt=1/fs;
    const res=zeroCrossing(data, dt);
    if(!res.nWaves) return alert('No valid waves found');
    currentWaveData=data;
    currentWaveResults=res;
    const totalDuration=data.length*dt;
    const heights=res.waves.map(w=>w.height);
    const maxH=Math.max(...heights);
    const minH=Math.min(...heights);
    const stdH=Math.sqrt(heights.reduce((s,h)=>s+(h-res.Hmean)**2,0)/res.nWaves);
    byId('hs-value').textContent=res.Hs.toFixed(3);
    byId('hmean-value').textContent=res.Hmean.toFixed(3);
    byId('tmean-value').textContent=res.Tmean.toFixed(2);
    byId('h10-value').textContent=res.H10.toFixed(3);
    byId('wave-count-value').textContent=res.nWaves;
    const container=byId('wave-stats-result');
    let ext=container.querySelector('.additional-results');
    if(!ext){ ext=document.createElement('div'); ext.className='additional-results'; container.appendChild(ext);}
    ext.innerHTML=`<strong>Extended Statistics:</strong><br>• Maximum Wave Height: ${maxH.toFixed(3)} m<br>• Minimum Wave Height: ${minH.toFixed(3)} m<br>• Standard Deviation: ${stdH.toFixed(3)} m<br>• Significant Period (Ts): ${res.Ts.toFixed(2)} s<br>• Data Duration: ${totalDuration.toFixed(1)} s<br>• Sampling Rate: ${fs} Hz<br>• Data Points: ${data.length}<br>• Parsing: Excel/CSV compatible (first numeric per row used)`;
    container.style.display='block';
    byId('wave-stats-plot-controls').style.display='block';
  } catch(e){
    alert('Error analyzing wave data.');
  }
};

// Reflection
window.loadReflectionExample=function(showAlert=true){
  // Set correct parameters for user's CSV data
  byId('gauge-positions').value = '0,0.3,0.6';  // User's requested positions
  byId('reflection-depth').value = '0.3';
  byId('time-step').value = '0.01';
  
  // Always load from myfile.csv
  fetch('myfile.csv').then(resp => resp.text()).then(txt => {
    const gauge = extractGaugeDataFromCSV(txt);
    if (gauge && gauge.length > 0) {
      byId('gauge-data').value = gauge.map(r => r.join(',')).join('\n');
      if(showAlert) alert('Reflection example loaded from myfile.csv with ' + gauge.length + ' data points.');
    } else {
      throw new Error('No valid data extracted from myfile.csv');
    }
  }).catch(e => {
    console.warn('Could not load myfile.csv:', e);
    // Fallback to generated data with correct positions
    const gauge = CACHED_EXAMPLE_DATA;
    byId('gauge-data').value = gauge.map(r => r.join(',')).join('\n');
    if(showAlert) alert('Reflection example loaded (using fallback data).');
  });
};
window.calculateReflection=function(){
  try {
    const posStr = byId('gauge-positions').value.trim();
    const depthStr = byId('reflection-depth').value.trim();
    const dtStr = byId('time-step').value.trim();
    const txt = byId('gauge-data').value.trim();
    
    console.log('Input validation:');
    console.log('Positions string:', posStr);
    console.log('Depth string:', depthStr);
    console.log('Time step string:', dtStr);
    console.log('Data length:', txt.length);
    
    // Validate inputs
    if (!posStr || !depthStr || !dtStr || !txt) {
      alert('Please fill in all fields: gauge positions, water depth, time step, and gauge data.');
      return;
    }
    
    const pos = posStr.split(',').map(x => parseFloat(x.trim())).filter(x => !isNaN(x));
    const h = parseFloat(depthStr);
    const dt = parseFloat(dtStr);
    
    console.log('Parsed values:');
    console.log('Positions:', pos);
    console.log('Depth:', h);
    console.log('Time step:', dt);
    
    if (pos.length !== 3) {
      alert('Please enter exactly 3 gauge positions separated by commas (e.g., 0,0.3,0.6)');
      return;
    }
    
    if (isNaN(h) || h <= 0) {
      alert('Please enter a valid positive water depth');
      return;
    }
    
    if (isNaN(dt) || dt <= 0) {
      alert('Please enter a valid positive time step');
      return;
    }
    
    const data = parseThreeGaugeMatrix(txt);
    console.log('Parsed data rows:', data.length);
    
    if (data.length < 10) {
      alert(`Need at least 10 rows of data with 3 columns each. Found ${data.length} valid rows.`);
      return;
    }
    
    currentReflectionData = data;
    const eta = [data.map(r => r[0]), data.map(r => r[1]), data.map(r => r[2])];
    
    console.log('Starting threeArray analysis...');
    console.log('eta structure:', eta.length, 'arrays of length:', eta.map(arr => arr.length));
    console.log('eta sample data:', eta.map(arr => arr.slice(0, 3)));
    
    // MATLAB function signature: three_array(h, dt, z, gpos)
    // h = water depth, dt = time step, z = gauge data matrix, gpos = gauge positions
    const res = threeArray(h, dt, eta, pos);
    console.log('Analysis complete:', res);
    
    // Check for required elements before accessing
    const hiValueEl = byId('hi-value');
    const hrValueEl = byId('hr-value');
    const krValueEl = byId('kr-value');
    const freqCountEl = byId('freq-count-value');
    const containerEl = byId('reflection-result');
    const plotControlsEl = byId('reflection-plot-controls');
    
    if (!hiValueEl || !hrValueEl || !krValueEl || !freqCountEl || !containerEl || !plotControlsEl) {
      console.error('Missing DOM elements:', {
        'hi-value': !!hiValueEl,
        'hr-value': !!hrValueEl,
        'kr-value': !!krValueEl,
        'freq-count-value': !!freqCountEl,
        'reflection-result': !!containerEl,
        'reflection-plot-controls': !!plotControlsEl
      });
      alert('Error: Required page elements not found. Please refresh the page.');
      return;
    }
    
    hiValueEl.textContent = res.Hi.toFixed(3);
    hrValueEl.textContent = res.Hr.toFixed(3);
    krValueEl.textContent = res.Kr.toFixed(3);
    freqCountEl.textContent = res.validFrequencies;
    
    let ext = containerEl.querySelector('.additional-results');
    if (!ext) {
      ext = document.createElement('div');
      ext.className = 'additional-results';
      containerEl.appendChild(ext);
    }
    
    const kr = res.Kr;
    let perf = kr < 0.3 ? 'Low reflection - Good absorption' : 
               kr < 0.7 ? 'Moderate reflection - Partial absorption' : 
               'High reflection - Strong reflection';
    const totalDuration = data.length * dt;
    const spacing = Math.abs(pos[1] - pos[0]);
    const energyRefPct = (kr * kr * 100).toFixed(1);
    
    ext.innerHTML = `<strong>Analysis Details:</strong><br>• Gauge Spacing: ${spacing.toFixed(2)} m<br>• Data Duration: ${totalDuration.toFixed(1)} s<br>• Performance: ${perf}<br>• Energy Reflection: ${energyRefPct}%<br>• Data Points: ${data.length}<br>• Method: 3-gauge LSQ (freq-domain)<br>• Parsing: Excel/CSV/tab delimited accepted (first 3 numerics per row)`;
    
    containerEl.style.display = 'block';
    containerEl.classList.add('success');
    plotControlsEl.style.display = 'block';
    
    // Automatically plot results
    plotGaugeData();
    plotReflectionSpectra();
    
  } catch (e) {
    console.error('Error in calculateReflection:', e);
    alert('Error analyzing reflection: ' + e.message);
    const containerEl = byId('reflection-result');
    if (containerEl) {
      containerEl.classList.add('error');
    }
  }
};

// Plotting helpers (Plotly.js)
window.plotReflectionSpectra=function(){
  if(!currentReflectionData) return alert('No gauge data');
  const pos=byId('gauge-positions').value.split(',').map(x=>parseFloat(x.trim()));
  const h=parseFloat(byId('reflection-depth').value);
  const dt=parseFloat(byId('time-step').value);
  const data=currentReflectionData;
  const eta=[data.map(r=>r[0]), data.map(r=>r[1]), data.map(r=>r[2])];
  const res=threeArray(pos, eta, dt, h);
  // Extract spectra and frequency arrays
  const Si=res.details.map(d=>d.Si);
  const Sr=res.details.map(d=>d.Sr);
  const Sf=res.details.map(d=>d.Sf);
  const f=res.details.map(d=>d.freq);
  // Band averaging (5 bands)
  const noBands=5;
  const bandSize=Math.floor(Si.length/noBands);
  const flimBand=[], SiBand=[], SrBand=[], SfBand=[];
  for(let j=0;j<noBands;j++){
    const start=j*bandSize;
    const end=(j+1)*bandSize;
    flimBand.push(mean(f.slice(start,end)));
    SiBand.push(mean(Si.slice(start,end)));
    SrBand.push(mean(Sr.slice(start,end)));
    SfBand.push(mean(Sf.slice(start,end)));
  }
  function mean(arr){
    if(!arr.length) return NaN;
    return arr.reduce((a,b)=>a+b,0)/arr.length;
  }
  // Plot 1: Incident, Reflected, Composite Spectra
  Plotly.newPlot('reflection-spectra-plot', [
    {
      x: f, 
      y: Si, 
      mode: 'lines', 
      name: 'Incident', 
      line: {color: '#0ea5e9', width: 3},
      hovertemplate: 'Frequency: %{x:.2f} Hz<br>Spectral Density: %{y:.4f} m²·s<extra></extra>'
    },
    {
      x: f, 
      y: Sr, 
      mode: 'lines', 
      name: 'Reflected', 
      line: {color: '#dc2626', width: 3, dash: 'dash'},
      hovertemplate: 'Frequency: %{x:.2f} Hz<br>Spectral Density: %{y:.4f} m²·s<extra></extra>'
    },
    {
      x: f, 
      y: Sf, 
      mode: 'lines', 
      name: 'Composite', 
      line: {color: '#059669', width: 3},
      hovertemplate: 'Frequency: %{x:.2f} Hz<br>Spectral Density: %{y:.4f} m²·s<extra></extra>'
    }
  ], {
    title: {
      text: 'Wave Spectra Analysis',
      font: {size: 20, family: 'Inter, Arial', color: '#0c4a6e'},
      x: 0.5
    },
    xaxis: {
      title: 'Frequency [Hz]',
      gridcolor: 'rgba(14, 165, 233, 0.1)',
      tickfont: {size: 12}
    },
    yaxis: {
      title: 'Spectral Density [m²·s]',
      gridcolor: 'rgba(14, 165, 233, 0.1)',
      tickfont: {size: 12}
    },
    legend: {
      orientation: 'h',
      x: 0.5,
      xanchor: 'center',
      y: -0.25,
      font: {size: 12}
    },
    margin: {t: 60, l: 60, r: 40, b: 100},
    autosize: true,
    plot_bgcolor: 'white',
    paper_bgcolor: 'white',
    font: {family: 'Inter, Arial'},
    hoverlabel: {
      bgcolor: 'white',
      bordercolor: '#0c4a6e',
      font: {color: '#0c4a6e', size: 14, family: 'Inter, Arial'}
    }
  }, {
    responsive: true, 
    displayModeBar: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'png',
      filename: 'wave_spectra_analysis',
      scale: 2
    }
  });
  byId('reflection-spectra-plot-container').style.display='block';
};

window.plotReflectionSpectra=function(){
  if(!currentReflectionData) return alert('No gauge data');
  
  const pos=byId('gauge-positions').value.split(',').map(x=>parseFloat(x.trim()));
  const h=parseFloat(byId('reflection-depth').value);
  const dt=parseFloat(byId('time-step').value);
  const data=currentReflectionData;
  const eta=[data.map(r=>r[0]), data.map(r=>r[1]), data.map(r=>r[2])];
  
  console.log('DEBUG: Starting reflection spectra plotting...');
  console.log('Gauge positions:', pos);
  console.log('Water depth:', h);
  console.log('Time step:', dt);
  console.log('Data points:', data.length);
  console.log('Sample data rows:', data.slice(0, 3));
  
  console.log('eta structure:', eta.length, 'arrays of length:', eta.map(arr => arr.length));
  console.log('eta sample data:', eta.map(arr => arr.slice(0, 3)));
  
  // MATLAB function signature: three_array(h, dt, z, gpos)
  const res=threeArray(h, dt, eta, pos);
  console.log('Analysis result:', res);
  
  // Extract spectra and frequency arrays
  const Si=res.details.map(d=>d.Si);
  const Sr=res.details.map(d=>d.Sr);
  const Sf=res.details.map(d=>d.Sf);
  const f=res.details.map(d=>d.freq);
  
  console.log('Extracted arrays:');
  console.log('Si (Incident):', Si.slice(0, 5), '... length:', Si.length);
  console.log('Sr (Reflected):', Sr.slice(0, 5), '... length:', Sr.length);
  console.log('Sf (Composite):', Sf.slice(0, 5), '... length:', Sf.length);
  console.log('f (Frequency):', f.slice(0, 5), '... length:', f.length);
  
  // NaN check and warning
  let errorMsg = '';
  if (Si.every(isNaN) || Sr.every(isNaN)) {
    errorMsg += 'Incident and/or Reflected spectra are all NaN.\n';
    if (eta[0].length < 32) {
      errorMsg += 'Gauge data is very short. FFT and reflection analysis require at least 32 points.\n';
    }
    if (pos.length !== 3 || new Set(pos).size !== 3) {
      errorMsg += 'Gauge positions must be three distinct values.\n';
    }
    if (!(h > 0) || !(dt > 0)) {
      errorMsg += 'Water depth and time step must be positive numbers.\n';
    }
    if (errorMsg) {
      alert(errorMsg.trim());
      return;
    }
  }
  
  // Check for plot container
  const plotEl = byId('reflection-spectra-plot');
  const containerEl = byId('reflection-spectra-plot-container');
  
  if (!plotEl || !containerEl) {
    console.error('Plot elements not found:', {
      'reflection-spectra-plot': !!plotEl,
      'reflection-spectra-plot-container': !!containerEl
    });
    alert('Plot container not found. Please refresh the page.');
    return;
  }
  
  // Plot 1: Incident, Reflected, Composite Spectra
  Plotly.newPlot('reflection-spectra-plot', [
    {
      x: f, 
      y: Si, 
      mode: 'lines+markers', 
      name: 'Incident Wave', 
      line: {color: '#0ea5e9', width: 4, dash: 'dot', shape: 'spline'}, 
      marker: {size: 7, color: '#0ea5e9', symbol: 'circle'},
      hovertemplate: '<b>Incident Wave</b><br>Frequency: %{x:.3f} Hz<br>Spectral Density: %{y:.6f} m²·s<extra></extra>'
    },
    {
      x: f, 
      y: Sr, 
      mode: 'lines+markers', 
      name: 'Reflected Wave', 
      line: {color: '#dc2626', width: 4, dash: 'dashdot', shape: 'spline'}, 
      marker: {size: 7, color: '#dc2626', symbol: 'diamond'},
      hovertemplate: '<b>Reflected Wave</b><br>Frequency: %{x:.3f} Hz<br>Spectral Density: %{y:.6f} m²·s<extra></extra>'
    },
    {
      x: f, 
      y: Sf, 
      mode: 'lines+markers', 
      name: 'Composite Spectrum', 
      line: {color: '#059669', width: 5, shape: 'spline'}, 
      marker: {size: 7, color: '#059669', symbol: 'square'},
      hovertemplate: '<b>Composite Spectrum</b><br>Frequency: %{x:.3f} Hz<br>Spectral Density: %{y:.6f} m²·s<extra></extra>'
    }
  ], {
    title: {
      text: 'Three-Gauge Reflection Analysis: Wave Spectra Decomposition',
      font: {size: 20, family: 'Inter, Arial', color: '#0c4a6e'},
      x: 0.5
    },
    xaxis: {
      title: {
        text: 'Frequency f [Hz]',
        font: {size: 20, color: '#0369a1', weight: 'bold'}
      },
      gridcolor: 'rgba(14, 165, 233, 0.2)',
      zerolinecolor: 'rgba(14, 165, 233, 0.4)',
      tickfont: {size: 16, color: '#0c4a6e'},
      showline: true,
      linewidth: 3,
      linecolor: 'rgba(14, 165, 233, 0.4)',
      mirror: true
    },
    yaxis: {
      title: {
        text: 'Spectral Density S<sub>f</sub> [m²·s]',
        font: {size: 20, color: '#0369a1', weight: 'bold'}
      },
      gridcolor: 'rgba(14, 165, 233, 0.2)',
      zerolinecolor: 'rgba(14, 165, 233, 0.4)',
      tickfont: {size: 16, color: '#0c4a6e'},
      showline: true,
      linewidth: 3,
      linecolor: 'rgba(14, 165, 233, 0.4)',
      mirror: true
    },
    legend: {
      orientation: 'h',
      x: 0.5,
      xanchor: 'center',
      y: -0.25,
      bgcolor: 'rgba(255, 255, 255, 0.95)',
      bordercolor: 'rgba(14, 165, 233, 0.4)',
      borderwidth: 3,
      font: {size: 18, color: '#0c4a6e', weight: 'bold'}
    },
    margin: {t: 120, l: 100, r: 60, b: 140},
    autosize: true,
    plot_bgcolor: 'rgba(240, 249, 255, 0.6)',
    paper_bgcolor: 'rgba(255, 255, 255, 0.95)',
    font: {family: 'Inter, Arial', color: '#0c4a6e'},
    hoverlabel: {
      bgcolor: 'white',
      bordercolor: '#0c4a6e',
      font: {color: '#0c4a6e', size: 14, family: 'Inter, Arial'}
    },
    annotations: [{
      text: 'Goda & Suzuki (1976) Method',
      x: 1,
      y: 1,
      xref: 'paper',
      yref: 'paper',
      showarrow: false,
      font: {size: 14, color: '#64748b', style: 'italic'},
      xanchor: 'right',
      yanchor: 'top'
    }]
  }, {
    responsive: true, 
    displayModeBar: true,
    modeBarButtonsToAdd: ['drawline', 'drawopenpath', 'drawclosedpath', 'drawcircle', 'drawrect', 'eraseshape'],
    displaylogo: false,
    toImageButtonOptions: {
      format: 'png',
      filename: 'reflection_analysis_spectra',
      height: 700,
      width: 1200,
      scale: 2
    }
  });
  
  containerEl.style.display='block';
};
window.plotGaugeData=function(){
  if(!currentReflectionData) return alert('No gauge data');
  const dt=parseFloat(byId('time-step').value);
  const time=currentReflectionData.map((_,i)=>(i*dt).toFixed(3));
  const pos=byId('gauge-positions').value.split(',').map(x=>parseFloat(x.trim()));
  const colors=['#0ea5e9','#dc2626','#059669'];
  const symbols=['circle','diamond','square'];
  const labels = [`Gauge 1 (x=${pos[0]}m)`, `Gauge 2 (x=${pos[1]}m)`, `Gauge 3 (x=${pos[2]}m)`];
  const traces = [0,1,2].map(i => ({
    x: time,
    y: currentReflectionData.map(r => r[i]),
    mode: 'lines',
    name: labels[i],
    line: {color: colors[i], width: 2},
    hovertemplate: `${labels[i]}<br>Time: %{x:.2f} s<br>Elevation: %{y:.3f} m<extra></extra>`
  }));
  
  Plotly.newPlot('reflection-plot', traces, {
    title: {
      text: 'Three-Gauge Wave Elevation Time Series',
      font: {size: 20, family: 'Inter, Arial', color: '#0c4a6e'},
      x: 0.5
    },
    xaxis: {
      title: 'Time [s]',
      gridcolor: 'rgba(14, 165, 233, 0.1)',
      tickfont: {size: 12}
    },
    yaxis: {
      title: 'Wave Elevation [m]',
      gridcolor: 'rgba(14, 165, 233, 0.1)',
      tickfont: {size: 12}
    },
    legend: {
      orientation: 'h',
      x: 0.5,
      xanchor: 'center',
      y: -0.25,
      font: {size: 12}
    },
    margin: {t: 60, l: 60, r: 40, b: 100},
    autosize: true,
    plot_bgcolor: 'white',
    paper_bgcolor: 'white',
    font: {family: 'Inter, Arial'},
    hoverlabel: {
      bgcolor: 'white',
      bordercolor: '#0c4a6e',
      font: {color: '#0c4a6e', size: 14, family: 'Inter, Arial'}
    }
  }, {
    paper_bgcolor: 'white',
    font: {family: 'Inter, Arial'}
  }, {
    responsive: true,
    displayModeBar: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'png',
      filename: 'gauge_time_series',
      scale: 2
    }
  });
    byId('reflection-plot-container').style.display='block';
  };
  

  window.plotWaveTimeSeries=function(){
    if(!currentWaveData) return alert('No wave data');
    const fs=parseFloat(byId('sampling-frequency').value);
    const time=currentWaveData.map((_,i)=>(i/fs).toFixed(3));
    Plotly.newPlot('wave-plot', [{
      x: time,
      y: currentWaveData,
      type: 'scatter',
      mode: 'lines',
      name: 'Wave Elevation',
      line: {color: '#0ea5e9', width: 2},
      hovertemplate: 'Time: %{x:.2f} s<br>Elevation: %{y:.3f} m<extra></extra>'
    }], {
      title: {
        text: 'Wave Elevation Time Series',
        font: {size: 20, family: 'Inter, Arial', color: '#0c4a6e'},
        x: 0.5
      },
      xaxis: {
        title: {
          text: 'Time t [s]',
          font: {size: 20, color: '#0369a1', weight: 'bold'}
        },
        gridcolor: 'rgba(14, 165, 233, 0.2)',
        zerolinecolor: 'rgba(14, 165, 233, 0.4)',
        tickfont: {size: 16, color: '#0c4a6e'},
        showline: true,
        linewidth: 3,
        linecolor: 'rgba(14, 165, 233, 0.4)',
        mirror: true
      },
      yaxis: {
        title: {
          text: 'Wave Elevation η [m]',
          font: {size: 20, color: '#0369a1', weight: 'bold'}
        },
        gridcolor: 'rgba(14, 165, 233, 0.2)',
        zerolinecolor: 'rgba(14, 165, 233, 0.4)',
        tickfont: {size: 16, color: '#0c4a6e'},
        showline: true,
        linewidth: 3,
        linecolor: 'rgba(14, 165, 233, 0.4)',
        mirror: true
      },
      legend: {
        orientation: 'h',
        x: 0.5,
        xanchor: 'center',
        y: -0.2,
        bgcolor: 'rgba(255, 255, 255, 0.95)',
        bordercolor: 'rgba(14, 165, 233, 0.4)',
        borderwidth: 3,
        font: {size: 18, color: '#0c4a6e', weight: 'bold'}
      },
      margin: {t: 120, l: 100, r: 60, b: 120},
      autosize: true,
      plot_bgcolor: 'rgba(240, 249, 255, 0.6)',
      paper_bgcolor: 'rgba(255, 255, 255, 0.95)',
      font: {family: 'Inter, Arial', color: '#0c4a6e'},
      hoverlabel: {
        bgcolor: 'white',
        bordercolor: '#0c4a6e',
        font: {color: '#0c4a6e', size: 14, family: 'Inter, Arial'}
      },
      annotations: [{
        text: `Sampling Frequency: ${fs} Hz`,
        x: 1,
        y: 1,
        xref: 'paper',
        yref: 'paper',
        showarrow: false,
        font: {size: 14, color: '#64748b', style: 'italic'},
        xanchor: 'right',
        yanchor: 'top'
      }]
    }, {
      responsive: true, 
      displayModeBar: true,
      modeBarButtonsToAdd: ['drawline', 'drawopenpath', 'drawclosedpath', 'drawcircle', 'drawrect', 'eraseshape'],
      displaylogo: false,
      toImageButtonOptions: {
        format: 'png',
        filename: 'wave_time_series',
        height: 700,
        width: 1200,
        scale: 2
      }
    });
    byId('wave-plot-container').style.display='block';
  };

  window.plotWaveHeights = function() {
    if(!currentWaveResults || !currentWaveResults.waves) return alert('No wave stats');
    const heights = currentWaveResults.waves.map(w => w.height);
    
    Plotly.newPlot('wave-heights-plot', [{
      x: heights.map((_, i) => i + 1),
      y: heights,
      type: 'bar',
      name: 'Wave Heights',
      marker: {
        color: 'rgba(14, 165, 233, 0.7)',
        line: {width: 1, color: '#0ea5e9'}
      },
      hovertemplate: 'Wave %{x}<br>Height: %{y:.3f} m<extra></extra>'
    }], {
      title: {
        text: 'Individual Wave Heights Distribution',
        font: {size: 20, family: 'Inter, Arial', color: '#0c4a6e'},
        x: 0.5
      },
      xaxis: {
        title: 'Wave Number',
        gridcolor: 'rgba(14, 165, 233, 0.1)',
        tickfont: {size: 12}
      },
      yaxis: {
        title: 'Wave Height [m]',
        gridcolor: 'rgba(14, 165, 233, 0.1)',
        tickfont: {size: 12}
      },
      margin: {t: 60, l: 60, r: 40, b: 60},
      autosize: true,
      plot_bgcolor: 'white',
      paper_bgcolor: 'white',
      font: {family: 'Inter, Arial'},
      showlegend: false
    }, {
      responsive: true,
      displayModeBar: true,
      displaylogo: false,
      toImageButtonOptions: {
        format: 'png',
        filename: 'wave_heights_distribution',
        scale: 2
      }
    });
    byId('wave-heights-plot-container').style.display = 'block';
  };

window.plotGaugeData=function(){
  if(!currentReflectionData) return alert('No gauge data');
  
  // Check for required elements
  const plotEl = byId('reflection-plot');
  const containerEl = byId('reflection-plot-container');
  const spectraContainerEl = byId('reflection-spectra-plot-container');
  
  if (!plotEl || !containerEl) {
    console.error('Plot elements not found:', {
      'reflection-plot': !!plotEl,
      'reflection-plot-container': !!containerEl
    });
    alert('Plot container not found. Please refresh the page.');
    return;
  }
  
  const dt=parseFloat(byId('time-step').value);
  const time=currentReflectionData.map((_,i)=>(i*dt).toFixed(3));
  const pos=byId('gauge-positions').value.split(',').map(x=>parseFloat(x.trim()));
  const colors=['#667eea','#f093fb','#4facfe'];
  const labels=[`Gauge 1 (x=${pos[0]}m)`,`Gauge 2 (x=${pos[1]}m)`,`Gauge 3 (x=${pos[2]}m)`];
  const traces=[0,1,2].map(i=>({
    x: time,
    y: currentReflectionData.map(r=>r[i]),
    mode: 'lines',
    name: labels[i],
    line: {color: colors[i], width: 3, shape: 'spline'}
  }));
  Plotly.newPlot('reflection-plot', traces, {
    title: {text:'Three-Gauge Time Series Data', font:{size:24, family:'Inter, Arial'}},
    xaxis: {title: 'Time t (s)', gridcolor:'#e0e7ff', zerolinecolor:'#cbd5e1'},
    yaxis: {title: 'Wave Elevation η (m)', gridcolor:'#e0e7ff', zerolinecolor:'#cbd5e1'},
    legend: {orientation: 'h', font:{size:16}},
    margin: {t: 60, l: 60, r: 30, b: 50},
    autosize: true,
    plot_bgcolor:'#f8fafc',
    paper_bgcolor:'#f8fafc',
  }, {responsive: true, displayModeBar: false});
  
  containerEl.style.display='block';
};

window.plotFrequencySpectrum=function(){
    if(!currentReflectionData) return alert('No gauge data');
    const dt=parseFloat(byId('time-step').value);
    const fs=1/dt;
    const fft=FFT.fft(currentReflectionData.map(r=>r[0]));
    const n=fft.length;
    const freqs=[];
    const energy=[];
    for (let i=0;i<n/2;i++){
      freqs.push(i*fs/n);
      const mag=FFT.magnitude(fft[i])/n;
      energy.push((mag*mag)/fs);
    }
    Plotly.newPlot('reflection-spectrum-plot', [{
      x: freqs,
      y: energy,
      type: 'scatter',
      mode: 'lines',
      name: 'Energy Spectrum',
      line: {color: '#667eea', width: 3, shape: 'spline'},
      fill: 'tozeroy',
      fillcolor: 'rgba(102,126,234,0.1)'
    }], {
      title: {text:'Frequency Spectrum Analysis (Gauge 1)', font:{size:24, family:'Inter, Arial'}},
      xaxis: {title: 'Frequency f (Hz)', gridcolor:'#e0e7ff', zerolinecolor:'#cbd5e1'},
      yaxis: {title: 'Energy Spectrum (m²/Hz)', gridcolor:'#e0e7ff', zerolinecolor:'#cbd5e1'},
      margin: {t: 60, l: 60, r: 30, b: 50},
      autosize: true,
      plot_bgcolor:'#f8fafc',
      paper_bgcolor:'#f8fafc',
    }, {responsive: true, displayModeBar: false});
    byId('reflection-spectra-plot-container').style.display='block';
  };
// Removed duplicate and misplaced code. All plotting functions now have correct closing brackets.

// DOM Ready check for essential elements & auto-load examples
window.addEventListener('DOMContentLoaded', ()=>{
  console.log('DOM loaded, initializing Wave App...');
  
  // Add a small delay to ensure all elements are rendered
  setTimeout(() => {
    try {
      // Check for essential elements
      const requiredElements = [
        'water-depth', 'wave-period', 'wavelength-value',
        'wave-data', 'sampling-frequency', 'hs-value',
        'gauge-positions', 'reflection-depth', 'time-step', 'gauge-data',
        'hi-value', 'hr-value', 'kr-value', 'freq-count-value',
        'reflection-result', 'reflection-plot-controls',
        'reflection-plot', 'reflection-plot-container',
        'reflection-spectra-plot', 'reflection-spectra-plot-container'
      ];
      
      const missingElements = requiredElements.filter(id => !byId(id));
      if (missingElements.length > 0) {
        console.error('Missing required elements:', missingElements);
        alert('Page not fully loaded. Please refresh the browser.');
        return;
      }
      
      console.log('All DOM elements found, proceeding with initialization...');
      
      // Auto-load examples and initialize UI
      console.log('Loading wavelength example...');
      loadWavelengthExample(false);
      calculateWavelength();
      
      console.log('Loading wave stats example...');
      loadWaveStatsExample(false);
      calculateWaveStats();
      
      console.log('Loading reflection example...');
      loadReflectionExample(false);
      
      // Don't auto-run reflection analysis - let user click the button
      console.log('Wave App initialization complete');
    } catch(e){
      console.error('Initialization error:', e);
      alert('Error initializing app: ' + e.message);
    }
  }, 100); // 100ms delay to ensure DOM is fully rendered
});
