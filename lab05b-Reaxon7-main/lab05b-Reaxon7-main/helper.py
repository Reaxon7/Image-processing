import numpy as np


def DTFT(x, M):
    """
    Parameters:
    ---
    x: a signal which is assumed to start at time n = 0
    M: the number of output points of the DTFT
    
    Returns:
    ---
    X: the samples of the DTFT
    w: corresponding frequencies of these samples
    """
    N = max(M, len(x))
    N = int(np.power(2, np.ceil(np.log(N) / np.log(2))))
    X = np.fft.fft(x, N)
    w = np.arange(N) / N * 2 * np.pi
    w = w - 2 * np.pi * (w >= np.pi).astype(int)
    X = np.fft.fftshift(X)
    w = np.fft.fftshift(w)
    return X, w


def hanning(N):
    """
    Parameters:
    ---
    N: the length of the Hanning window

    Returns:
    ---
    w: the Hanning window of length N
    """
    n = np.arange(N)
    w = 0.5 - 0.5 * np.cos(2 * np.pi * n / (N - 1))
    return w


def hamming(N):
    """
    Parameters:
    ---
    N: the length of the Hamming window

    Returns:
    ---
    w: the Hamming window of length N
    """
    n = np.arange(N)
    w = 0.54 - 0.46 * np.cos(2 * np.pi * n / (N - 1))
    return w


def blackman(N):
    """
    Parameters:
    ---
    N: the length of the Blackman window

    Returns:
    ---
    w: the Blackman window of length N
    """
    n = np.arange(N)
    w = 0.42 - 0.5 * np.cos(2 * np.pi * n / (N - 1)) + 0.08 * np.cos(4 * np.pi * n / (N - 1))
    return w


def kaiser(N, beta):
    """
    Parameters:
    ---
    N: the length of the Kaiser window
    beta: 

    Returns:
    ---
    w: the Kaiser window of length N and beta
    """
    n = np.arange(N)
    w = np.i0(beta * np.sqrt(1 - ((n - (N - 1) / 2)/((N - 1) / 2)) ** 2)) / np.i0(beta)
    return w


def remlpord(freq1, freq2, delta1, delta2):
    AA = np.array([[-4.278e-01, -4.761e-01, 0], [-5.941e-01, 7.114e-02, 0], [-2.660e-03, 5.309e-03, 0]])
    d1 = np.log10(delta1.item())
    d2 = np.log10(delta2.item())
    D = np.array([[1, d1, d1 * d1]]) @ AA @ np.array([[1], [d2], [d2 * d2]])
    bb = np.array([[11.01217], [0.51244]])
    fK = np.array([[1.0, d1 - d2]]) @ bb
    df = abs(freq2 - freq1)
    L = D / df - fK * df + 1
    return L


def firpmord(fcuts, mags, devs, fsamp):
    fcuts = np.array(fcuts)
    mags = np.array(mags)
    devs = np.array(devs)
    fcuts = fcuts / fsamp
    if np.any(fcuts > 1 / 2):
        raise ValueError("Invalid range.")
    if np.any(fcuts < 0):
        raise ValueError("The frequency must be positive.")
    
    # turn vectors into column vectors
    fcuts = fcuts.reshape(-1, 1)
    mags = mags.reshape(-1, 1)
    devs = devs.reshape(-1, 1)

    mf = fcuts.shape[0]
    mm = mags.shape[0]
    nbands = mm

    if len(mags) != len(devs):
        raise ValueError("Mismatched vector length.")

    if mf != 2 * (nbands - 1):
        raise ValueError("Invalid length.")
    
    zz = (mags == 0).astype(float)
    devs = devs / (zz + mags)

    f1 = fcuts[:mf - 1:2]
    f2 = fcuts[1:mf:2]

    n = np.argmin(f2 - f1)
    if nbands == 2:
        L = remlpord(f1[n], f2[n], devs[0], devs[1])
    else:
        L = 0
        for i in range(1, nbands - 1):
            L1 = remlpord(f1[i - 1], f2[i - 1], devs[i], devs[i - 1])
            L2 = remlpord(f1[i], f2[i], devs[i], devs[i + 1])
            L = max(L, max(L1, L2))
    
    N =  int(np.ceil(L) - 1)

    # ff = np.array([[0], [2 * fcuts], [1]])
    ff = np.concatenate((np.array([[0]]), 2 * fcuts, np.array([[1]])), axis=0)
    # aa = np.repeat(mags, 2, axis=0)
    aa = mags
    wts = np.ones_like(devs) * np.max(devs) / devs
    if aa[-1] != 0 and N % 2 != 0:
        N += 1
    
    return N, (ff * fsamp / 2).reshape(-1), aa.reshape(-1), wts.reshape(-1)
    

