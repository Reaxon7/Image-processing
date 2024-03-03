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