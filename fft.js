// Minimal FFT implementation for Wave App
// Source: https://github.com/dntj/jsfft (MIT License, simplified)
// Only real input supported, returns array of {re, im}
function fft(input) {
  // Pad input to next power of 2 if needed
  let N = input.length;
  let pow2 = 1 << Math.ceil(Math.log2(N));
  let arr = input.slice();
  if (N < pow2) {
    arr = arr.concat(Array(pow2 - N).fill(0));
    N = pow2;
  }
  if (N <= 1) return [{re: arr[0], im: 0}];
  // ...existing code...
  // Bit reversal
  let output = new Array(N);
  for (let i = 0; i < N; i++) {
    let j = 0;
    for (let k = 0; k < Math.log2(N); k++) {
      j |= ((i >> k) & 1) << (Math.log2(N) - 1 - k);
    }
    output[j] = {re: arr[i], im: 0};
  }
  // Cooley-Tukey
  for (let s = 1; s <= Math.log2(N); s++) {
    const m = 1 << s;
    const m2 = m >> 1;
    const wPhaseStep = -2 * Math.PI / m;
    for (let k = 0; k < N; k += m) {
      for (let j = 0; j < m2; j++) {
        const wPhase = wPhaseStep * j;
        const w = {re: Math.cos(wPhase), im: Math.sin(wPhase)};
        const t = {
          re: w.re * output[k + j + m2].re - w.im * output[k + j + m2].im,
          im: w.re * output[k + j + m2].im + w.im * output[k + j + m2].re
        };
        const u = output[k + j];
        output[k + j] = {re: u.re + t.re, im: u.im + t.im};
        output[k + j + m2] = {re: u.re - t.re, im: u.im - t.im};
      }
    }
  }
  return output;
}
function magnitude(c) { return Math.sqrt(c.re * c.re + c.im * c.im); }
window.FFT = { fft, magnitude };
