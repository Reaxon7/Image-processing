import numpy as np


def xcorr(x, y, maxlags):
    if len(x) != len(y):
        raise ValueError(f"The lengths of x and y are not equal!")
    r = np.correlate(x, y, "full")
    start_idx = np.argmax(r) - maxlags
    lags = np.arange(-maxlags, maxlags + 1)
    r = r[start_idx:start_idx + 2 * maxlags + 1]
    r = r / (len(x) - np.abs(lags))
    return lags, r


def quantiz(sig, partition, codebook):
    indx = np.zeros(len(sig)).astype(int)
    for i in range(len(partition)):
        indx += sig > partition[i]
    quantv = codebook[indx]
    distor = 0
    for i in range(len(codebook)):
        distor += np.sum(np.square(sig[indx == i] - codebook[i]))
    distor /= len(sig)
    return indx, quantv, distor



def lloyds(training_set, initial_codebook, tol=10e-7):
    """
    The code is adapted from the source code of lloyds.m from Matlab.

    Parameters
    ---
    training_set: the set that is used to estimate the probability density function
    initial_codebook: it contains an initial guess of the optimal quantization levels

    Returns
    ---
    partition:
    codebook:
    """

    if len(initial_codebook) < 2:
        raise ValueError("the length of the initial_book is less than 2")

    min_training = min(training_set)
    max_training = max(training_set)
    codebook = np.sort(initial_codebook)
    initial_codebook = len(codebook)
    partition = (codebook[1:] + codebook[:-1]) / 2

    index, _, distor = quantiz(training_set, partition, codebook)

    last_distor = 0

    ter_cond2 = np.finfo(float).eps * max_training
    if distor > ter_cond2:
        rel_distor = abs(distor - last_distor) / distor
    else:
        rel_distor = distor

    while rel_distor > tol:
        for i in range(initial_codebook):
            waste1 = np.argwhere(index == i).flatten()
            if len(waste1) > 0:
                codebook[i] = np.mean(training_set[waste1])
            else:
                if i == 0:
                    tmp = training_set[training_set <= partition[0]]
                    if len(tmp) == 0:
                        codebook[0] = (partition[0] + min_training) / 2
                    else:
                        codebook[0] = np.mean(tmp)
                elif i == initial_codebook - 1:
                    tmp = training_set[training_set >= partition[i - 1]]
                    if len(tmp) == 0:
                        codebook[i] = (max_training + partition[i - 1]) / 2
                    else:
                        codebook[i] = np.mean(tmp)
                else:
                    tmp = training_set[training_set >= partition[i - 1]]
                    tmp = tmp[tmp <= partition[i]]
                    if len(tmp) == 0:
                        codebook[i] = (partition[i] + partition[i - 1]) / 2
                    else:
                        codebook[i] = np.mean(tmp)

        partition = np.sort((codebook[1:] + codebook[:-1]) / 2)

        last_distor = distor
        index, _, distor = quantiz(training_set, partition, codebook)
        if distor > ter_cond2:
            rel_distor = abs(distor - last_distor) / distor
        else:
            rel_distor = distor
    return partition, codebook
