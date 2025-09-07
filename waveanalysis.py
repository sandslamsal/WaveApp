import numpy as np
from scipy.fftpack import fft, ifft
from scipy.signal import detrend
import matplotlib.pyplot as plt
import pandas as pd


def wavelen(h, T):
    tol=50
    g=9.81
    L0 = (g * T**2) / (2 * np.pi)
    L = L0
    L_next = L0
    for _ in range(tol):
        L_next = (g * T**2) / (2 * np.pi) * np.tanh(2 * np.pi * h / L)
        if np.abs(L_next - L) < 1e-5:
            break
        L = L_next
    return L



def three_array(h, dt, z, gpos, plot_flag=1):
    """
     Calculates the incident and reflected waves for random wave experiments
     conducted in laboratory flumes using a three-gauge array.

     Based on Goda and Suzuki (1976), expanded to three gauges by Kobayashi and Cox.
     Original code called IRSORT by Dan Cox, created by Mary Bryant on July 19, 2019 in MATLAB.
     Revised and adapted to Python by Sandesh Lamsal in August 14, 2023.

    Parameters:
      h        : water depth
      dt       : time step
      z        : (N x 3) matrix of surface elevations from the three gauges,
                 where z[:,0] is the first gauge, z[:,1] the second, and z[:,2] the third.
      gpos     : positions of the gauges (e.g., [0, 0.3, 0.9])
      plot_flag: flag to produce plots (default is 1)

    Returns:
      refanalysis: dictionary containing reflection analysis parameters and spectra.
    """
    
    # Remove the mean (detrend with a constant) from the gauge data
    z = detrend(z, type='constant')
    g = 9.81

    # nfft: number of points for the FFT (here, simply the length of z)
    # (In IRSORT, nfft=2^nextpow2(length(z)); here we use nfft=length(z))
    nfft = len(z)
    df = 1 / (nfft * dt)  # frequency resolution
    half = round(nfft / 2)
    
    # Preallocation
    An   = np.zeros((half - 1, 3))
    Bn   = np.zeros((half - 1, 3))
    Sn   = np.zeros((half - 1, 3))
    k    = np.zeros(half - 1)       # wavenumber (1D array)
    Ainc = np.zeros((half - 1, 3))
    Binc = np.zeros((half - 1, 3))
    Aref = np.zeros((half - 1, 3))
    Bref = np.zeros((half - 1, 3))
    nmin = np.zeros(3, dtype=int)
    nmax = np.zeros(3, dtype=int)
    
    # Solving Fourier coefficients (A1-B2) for all three gauges - eq.2
    # eta = A*cos(omega*t) + B*sin(omega*t)
    for j in range(3):
        fn = np.fft.fft(z[:, j], nfft)
        a0 = fn[0] / nfft
        An[:, j] = 2 * np.real(fn[1:half]) / nfft   # Real component
        Bn[:, j] = -2 * np.imag(fn[1:half]) / nfft   # Imaginary component
        
        fn_squared = fn * np.conj(fn)
        fn_fold = fn_squared[1:half] * 2
        Sn[:, j] = dt * np.real(fn_fold) / nfft       # Spectral/energy density
    
    # Frequency vector corresponding to the FFT coefficients
    f = df * np.arange(1, Sn.shape[0] + 1)
    
    # Required to solve for the wavenumber at each frequency - eq.6
    # w^2 = g*k*tanh(k*h)
    for j in range(len(f)):
        L = wavelen(h, 1 / f[j])
        k[j] = (2 * np.pi) / L
    
    # Solve for amplitudes of incident (ai) and reflected waves (ar) - eq.5
    # Three Gauge Array: gauge pairs defined by
    g1 = [1, 1, 2]
    g2 = [2, 3, 3]
    # (gpos example: [0, 0.3, 0.9] defines the distance from the first gauge)
    for j in range(3):
        A1 = An[:, g1[j] - 1]
        A2 = An[:, g2[j] - 1]
        B1 = Bn[:, g1[j] - 1]
        B2 = Bn[:, g2[j] - 1]
        pos1 = gpos[g1[j] - 1]
        pos2 = gpos[g2[j] - 1]
        
        term1 = -A2 * np.sin(k * pos1) + A1 * np.sin(k * pos2) + B2 * np.cos(k * pos1) - B1 * np.cos(k * pos2)
        term2 =  A2 * np.cos(k * pos1) - A1 * np.cos(k * pos2) + B2 * np.sin(k * pos1) - B1 * np.sin(k * pos2)
        term3 = -A2 * np.sin(k * pos1) + A1 * np.sin(k * pos2) - B2 * np.cos(k * pos1) + B1 * np.cos(k * pos2)
        term4 =  A2 * np.cos(k * pos1) - A1 * np.cos(k * pos2) - B2 * np.sin(k * pos1) + B1 * np.sin(k * pos2)
        
        Ainc[:, j] = term1 / (2 * np.sin(k * np.abs(pos2 - pos1)))
        Binc[:, j] = term2 / (2 * np.sin(k * np.abs(pos2 - pos1)))
        Aref[:, j] = term3 / (2 * np.sin(k * np.abs(pos2 - pos1)))
        Bref[:, j] = term4 / (2 * np.sin(k * np.abs(pos2 - pos1)))
        
        # Upper and lower limits of significant spectra - eq.7
        Lmin = np.abs((pos2 - pos1) / 0.45)  # Ranges suggested by Goda and Suzuki (1976)
        Lmax = np.abs((pos2 - pos1) / 0.05)
        kmax = 2 * np.pi / Lmin
        kmin = 2 * np.pi / Lmax
        ind = np.where((k <= kmax) & (k >= kmin))[0]
        nmin[j] = np.min(ind)
        nmax[j] = np.max(ind)
    
    # Set amplitudes outside the significant index range to NaN (averaging only valid data)
    for j in range(3):
        for i in range(Ainc.shape[0]):
            if i < nmin[j] or i > nmax[j]:
                Ainc[i, j] = np.nan
                Binc[i, j] = np.nan
                Aref[i, j] = np.nan
                Bref[i, j] = np.nan

    # Determine the overall index range from the three gauge pairs
    rng = np.arange(np.min(nmin), np.max(nmax) + 1, dtype=int)
    
    # Averaging overlapped amplitudes across gauges
    # Here we create full-length arrays (size = half-1) and only fill the 'rng' indices
    Aincav_full = np.zeros(Ainc.shape[0])
    Bincav_full = np.zeros(Binc.shape[0])
    Arefav_full = np.zeros(Aref.shape[0])
    Brefav_full = np.zeros(Bref.shape[0])
    
    Aincav_full[rng] = np.nanmean(Ainc[rng, :], axis=1)
    Bincav_full[rng] = np.nanmean(Binc[rng, :], axis=1)
    Arefav_full[rng] = np.nanmean(Aref[rng, :], axis=1)
    Brefav_full[rng] = np.nanmean(Bref[rng, :], axis=1)
    
    # Backing out spectra
    Si = (Aincav_full[rng]**2 + Bincav_full[rng]**2) / (2 * df)
    Sr = (Arefav_full[rng]**2 + Brefav_full[rng]**2) / (2 * df)
    Sfcheck = (An**2 + Bn**2) / (2 * df)  # should match Sn

    # Evaluate energies of resolved incident and reflected waves - eq.8
    Ei = np.sum(Si) * df
    Er = np.sum(Sr) * df

    # Reflection coefficient - eq.9
    refco = np.sqrt(Er / Ei)

    # Calculating incident, reflected, and total Hmo wave height
    mo = np.sum(Sn[rng, 0]) * df
    Htot = 4.004 * np.sqrt(mo)  # Hs = 4*sqrt(mo), Hs = Htot
    Hi = 4.004 * np.sqrt(Ei)
    Hr = 4.004 * np.sqrt(Er)
    Hicheck = np.mean(Htot) / (np.sqrt(1 + refco**2))
    Hrcheck = refco * np.mean(Htot) / (np.sqrt(1 + refco**2))
    
    # Band averaging spectra for plotting
    no_bands = 5
    bands = int(np.floor(len(Si) / no_bands))
    flim_band = np.zeros(bands)
    Si_band   = np.zeros(bands)
    Sr_band   = np.zeros(bands)
    Sf_band   = np.zeros(bands)
    
    flim = f[rng]
    Sf   = Sn[rng, 0]
    
    for j in range(bands):
        flim_band[j] = np.mean(flim[j * no_bands:(j + 1) * no_bands])
        Si_band[j]   = np.mean(Si[j * no_bands:(j + 1) * no_bands])
        Sr_band[j]   = np.mean(Sr[j * no_bands:(j + 1) * no_bands])
        Sf_band[j]   = np.mean(Sf[j * no_bands:(j + 1) * no_bands])
    
    # Only plot if plot_flag is set to 1
    if plot_flag == 1:
        plt.figure()
        fontName = 'Arial'
        fontSize = 12

        # Plot the incident, reflected, and composite spectra
        plt.subplot(3, 1, 1)
        plt.plot(f[rng], Si, ':b', f[rng], Sr, 'r-.', f[rng], Sn[rng, 0], 'k',
                 linewidth=2.5, marker='o')
        plt.xlim([f[rng[0]], f[rng[-1]]])
        plt.title('Incident, Reflected, and Composite Spectra')
        plt.legend(['Incident', 'Reflected', 'Composite'])
        plt.xlabel('Frequency [Hz]', fontname=fontName, fontsize=fontSize)
        plt.ylabel('S_f [m^2*s]')
        plt.ylim([0, np.max(Sn[rng, 0]) + np.max(Sn[rng, 0]) * 0.05])
        plt.grid(True)

        # Plot the band-averaged spectra
        plt.subplot(3, 1, 2)
        plt.plot(flim_band, Si_band, ':b', flim_band, Sr_band, 'r-.', flim_band, Sf_band, 'k',
                 linewidth=2.5, marker='o')
        plt.xlim([flim[0], flim[-1]])
        plt.title('Band-Averaged Spectra')
        plt.legend(['Incident', 'Reflected', 'Composite'])
        plt.xlabel('Frequency [Hz]', fontname=fontName, fontsize=fontSize)
        plt.ylabel('S_f [m^2*s]')
        plt.ylim([0, np.max(Sf_band) + np.max(Sf_band) * 0.05])
        plt.grid(True)

        # Plot the incident spectrum
        plt.subplot(3, 1, 3)
        plt.plot(flim_band, Si_band, '-b', linewidth=2.5, marker='o', markersize=0.5)
        plt.xlim([flim[0], flim[-1]])
        plt.title('Incident Spectrum', fontname=fontName, fontsize=fontSize)
        plt.xlabel('Frequency [Hz]', fontname=fontName, fontsize=fontSize)
        plt.ylabel('S_f [m^2*s]')
        plt.ylim([0, np.max(Sf_band) + np.max(Sf_band) * 0.05])
        plt.grid(True)

        # Reconstruct time series using Inverse Fourier Transform
        N = len(z)  # Number of time points
        eta_inc = np.zeros(N)
        eta_ref = np.zeros(N)
        t = dt * np.arange(N)
        # Loop over all frequency bins (only nonzero in 'rng' will contribute)
        for i in range(len(k)):
            omega = 2 * np.pi * f[i]
            eta_inc += Aincav_full[i] * np.cos(omega * t) + Bincav_full[i] * np.sin(omega * t)
            eta_ref += Arefav_full[i] * np.cos(omega * t) + Brefav_full[i] * np.sin(omega * t)
        
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(t, eta_inc, 'b', linewidth=2.0, marker='o', markersize=0.5)
        plt.title('Time Series of the Incident Wave', fontweight='bold',
                  fontname=fontName, fontsize=fontSize)
        plt.xlabel('Time [s]', fontname=fontName, fontsize=fontSize)
        plt.ylabel('Surface Elevation (η) [m]', fontname=fontName, fontsize=fontSize)
        plt.grid(True)
        
        plt.subplot(2, 1, 2)
        plt.plot(t, eta_ref, 'r', linewidth=2.5, marker='o', markersize=0.5)
        plt.title('Time Series of the Reflected Wave', fontweight='bold',
                  fontname=fontName, fontsize=fontSize)
        plt.xlabel('Time [s]', fontname=fontName, fontsize=fontSize)
        plt.ylabel('Surface Elevation (η) [m]', fontname=fontName, fontsize=fontSize)
        plt.grid(True)
        
        plt.show()
    
    # Output the reflection analysis structure as a dictionary
    refanalysis = {
        # 'refco': refco,
        # 'Htot': Htot,
        'Hi': Hi,
        'Hr': Hr,
        # 'Si': Si,
        # 'Sr': Sr,
        # 'Sn': Sn,
        # 'f': f,
        # 'flim': flim,
        # 'indrange': rng,
        # 'Hicheck': Hicheck,
        # 'Hrcheck': Hrcheck
    }
    
    return refanalysis


import numpy as np
import pandas as pd
from scipy.signal import detrend

def zero_crossing(data, frequency, threshold=None):

    print(f"Running zero_crossing from: {__file__}")
    

    """
    Zero crossing analysis of wave data - Python version adapted from original MATLAB code.
    Original MATLAB code by Urs Neumeier (modified by Sandesh Lamsal).
    Python version updated for full compatibility with MATLAB logic.

    Parameters
    ----------
    data : array-like, pandas Series, pandas DataFrame, or list of arrays
        Water surface elevation data. Can be a single array or list of arrays.
    frequency : float
        Sampling frequency in Hz.
    threshold : float, optional
        Minimum crest or trough height to be considered a valid wave. If None, 1% of Hmax is used.

    Returns
    -------
    res : dict
        Dictionary containing wave statistics and wave data (height, period).
    names : list
        List of output parameter names.
    """

    if frequency <= 0:
        raise ValueError('Frequency must be greater than zero')

    if isinstance(data, list):
        results = []
        names = None
        for sub_data in data:
            res, names = zero_crossing(sub_data, frequency, threshold)
            results.append(res)
        return results, names

    if isinstance(data, pd.DataFrame):
        data = data.values
    elif isinstance(data, pd.Series):
        data = data.values[:, np.newaxis]

    if data.ndim == 1:
        data = data[:, np.newaxis]

    if data.shape[1] > 1:
        if np.all(data[:, 0] > 720000) and np.all(data[:, 1] < 740000):
            data = data[:, 1:]
        if data.shape[1] > 1:
            raise ValueError('Data must have only one column')

    # Convert to column vector and invert (for downward crossing like in MATLAB)
    data = -data.flatten()

    names = ['Hs', 'Hmean', 'H1_10', 'Hmax', 'Tmean', 'Ts']
    res = [np.nan] * 6

    # Remove mean or trend
    data = detrend(data)

    nonzero_indices = np.where(data != 0)[0]
    d0 = data[nonzero_indices]

    # Find zero crossings
    crossings = nonzero_indices[np.where(d0[:-1] * d0[1:] < 0)[0]]

    # Remove first downward crossing if required
    if data[0] > 0:
        crossings = crossings[1:]

    crossings = crossings[::2]  # Keep only upward zero crossings

    if len(crossings) < 2:
        # No valid waves found
        res = {names[i]: res[i] for i in range(len(names))}
        res['wave'] = np.empty((0, 2))  # Empty wave array
        return res, names

    # Initialize wave matrix (height, crest, trough, period)
    wave = np.zeros((len(crossings) - 1, 4))

    for i in range(len(crossings) - 1):
        wave[i, 1] = np.max(data[crossings[i]:crossings[i+1]])  # crest
        wave[i, 2] = -np.min(data[crossings[i]:crossings[i+1]])  # trough
        wave[i, 3] = (crossings[i+1] - crossings[i]) / frequency  # period

    # Threshold: 1% of max wave height if not provided
    if threshold is None:
        threshold = 0.01 * np.max(wave[:, 1] + wave[:, 2])
    elif threshold < 0:
        raise ValueError('Wave threshold must not be negative')

    # Remove small waves by merging with neighbors
    i = 0
    while i < len(wave):
        if wave[i, 1] < threshold:
            if i > 0:
                # Merge with previous wave
                wave[i - 1, 1] = np.max([wave[i - 1, 1], wave[i, 1]])
                wave[i - 1, 2] = np.max([wave[i - 1, 2], wave[i, 2]])
                wave[i - 1, 3] += wave[i, 3]
                wave = np.delete(wave, i, axis=0)
                i -= 1
            else:
                # Remove first wave if it's too small
                wave = np.delete(wave, i, axis=0)
                i -= 1
        elif wave[i, 2] < threshold:
            if i < len(wave) - 1:
                # Merge with next wave
                wave[i, 1] = np.max([wave[i, 1], wave[i + 1, 1]])
                wave[i, 2] = np.max([wave[i, 2], wave[i + 1, 2]])
                wave[i, 3] += wave[i + 1, 3]
                wave = np.delete(wave, i + 1, axis=0)
            else:
                # Remove last wave if it's too small
                wave = np.delete(wave, i, axis=0)
                i -= 1
        i += 1

    # Compute total wave height = crest + trough
    wave[:, 0] = wave[:, 1] + wave[:, 2]

    # Number of valid waves
    nb = len(wave)

    if nb == 0:
        res = {names[i]: res[i] for i in range(len(names))}
        res['wave'] = np.empty((0, 2))
        return res, names

    # Sort waves by height (descending order)
    wave_unsorted = wave.copy()
    wave = wave[np.argsort(wave[:, 0])[::-1]]

    # Compute statistics
    res[0] = np.mean(wave[:round(nb / 3), 0])  # Hs (significant wave height)
    res[1] = np.mean(wave[:, 0])               # Hmean
    res[2] = np.mean(wave[:round(nb / 10), 0])  # H1/10
    res[3] = np.max(wave[:, 0])                 # Hmax
    res[4] = np.mean(wave[:, 3])                 # Tmean
    res[5] = np.mean(wave[:round(nb / 3), 3])   # Ts (significant period)

    # Convert result to dictionary
    res = {names[i]: res[i] for i in range(len(names))}
    res['wave'] = wave_unsorted[:, [0, 3]]  # [Height, Period]

    return res, names


