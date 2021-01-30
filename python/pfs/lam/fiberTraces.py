from pfs.lam.fileHandling import *
from pfs.lam.instrModelAlign import *
import matplotlib.pyplot as plt
from matplotlib import style
import os
from astropy.stats import sigma_clip
from collections import OrderedDict
import pfs.imageAnalysis as imeas
from scipy.interpolate import interp1d

shape = (4176, 4096)
ypix = np.arange(shape[0])

def paramGauss1D(y, x=None):
    x = np.arange(len(y)) if x is None else x
    offset = np.median(y)
    fy = y - offset
    amp = np.max(fy)
    mean = x[np.argmin(np.abs(fy - amp))]
    hmean = x[np.argmin(np.abs(fy - amp / 2))]
    sig = np.abs(mean - hmean) / (np.sqrt(2 * np.log(2)))

    popt1, pcov = imeas.curve_fit(imeas.oneD_Gaussian, x, y, p0=[amp, mean, sig, offset], maxfev=10000)

    newx = np.linspace(np.min(x), np.max(x), 1000)
    data = np.zeros((len(newx), 2))
    data[:, 0] = newx
    data[:, 1] = imeas.oneD_Gaussian(newx, *popt1)

    return popt1[1]


def estimatePeak(integ, peak, doGaussian):
    if doGaussian:
        try:
            peak = int(peak - 2) + paramGauss1D(integ[int(peak - 2):int(peak + 3)])
        except RuntimeError:
            pass

    return peak


def findThr(integ, fibers, perc=5):
    thr = np.percentile(integ, perc)
    peaks = peakdet(integ, thr)

    if len(peaks[0]) != len(fibers):
        return findThr(integ, fibers, perc=perc + 1)
    print(thr, perc)
    return peaks[0]


def distFibers(fibers, cam=None):
    camOffsets = dict(r1=70, b1=79)
    gap = 70 if cam is None else camOffsets[cam] 
    res = []
    for i in (range(len(fibers) - 1)):
        dist = (fibers[i] - fibers[i + 1]) * 6.171
        offset = gap if fibers[i] == 337 else 0
        res.append(dist - offset)

    return np.array(res)


def meanPeaks(data, thr, doPlot=False):
    mid = int(12 * data.shape[0] / 16)
    rows = data[mid - 10:mid + 10, :]
    integ = np.sum(rows, axis=0)

    if doPlot:
        plt.figure(figsize=(12, 6))
        plt.plot(integ)
        plt.grid()

    return peakdet(integ, thr)[0]


def missingFibers(peaks, fibers, cam, doPlot=True):
    print(len(peaks), len(fibers))
    missing = []
    dist = np.array([peaks[i + 1, 0] - peaks[i, 0] for i in range(peaks.shape[0] - 1)])
    disth = distFibers(fibers, cam)

    for i in range(len(dist)):
        delta = np.abs(dist[i] - disth[i])
        if delta > 3:
            print(delta)
            rem = fibers[i + 1]
            print(f'missing fiber {rem}')
            missing.append(rem)
            fibers.remove(rem)
            disth = distFibers(fibers, cam)

    if doPlot:
        fig = plt.figure(figsize=(12, 6))
        plt.plot(dist - disth, 'o-', label='dst')
        plt.grid()
        plt.legend()

    return missing, fibers


def fiberPeaks(data, fibers, sumRows=10, offPix=50, stepPix=100):
    xfibers = []
    yfibers = np.arange(offPix, data.shape[0] - offPix, stepPix)

    for i in yfibers:
        rows = data[i - sumRows:i + sumRows:]
        integ = np.sum(rows, axis=0)
        try:
            peaks = findThr(integ, fibers)
            xfibers.append([estimatePeak(integ, peak, doGaussian=True) for peak in peaks[:, 0]])

        except ValueError:
            print(f'row:{i} failed')
            xfibers.append([np.nan] * len(fibers))

        except RuntimeError:
            print(peak)

    return np.array(xfibers), yfibers


def fitFiberPeaks(xfibers, yfibers, doPlot=True):
    poly = []
    if doPlot:
        fig = plt.figure(figsize=(12, 6))
        plt.grid()
    for i in range(xfibers.shape[1]):
        mask = np.isnan(xfibers[:, i])
        xfiber, yfiber = xfibers[:, i][~mask], yfibers[~mask]
        p = np.polyfit(yfiber, xfiber, deg=2)
        resid = xfiber - np.polyval(p, yfiber)
        mask2 = sigma_clip(resid)

        p = np.polyfit(yfiber[~mask2.mask], xfiber[~mask2.mask], deg=4)
        poly.append(p)

        if doPlot:
            resid = xfiber - np.polyval(p, yfiber)
            plt.plot(xfiber[~mask2.mask], resid[~mask2.mask])

    return poly


def neighbourPoly(poly, missingFiber, doPlot=True):
    minFiber, maxFiber = (339, 650) if missingFiber >= 339 else (2, 315)
    minNeighBour, maxNeighBour = missingFiber - 10, missingFiber + 10
    minNeighBour = minFiber if minNeighBour < minFiber else minNeighBour
    maxNeighBour = maxFiber if maxNeighBour > maxFiber else maxNeighBour
    x = []
    coeff = []
    fibs = np.arange(minNeighBour, maxNeighBour)

    for fib in fibs:
        if fib in poly:
            x.append(fib)
            coeff.append(poly[fib])

    coeff = np.array(coeff)
    res = []
    for i in range(coeff.shape[1]):
        f = interp1d(x, coeff[:, i], kind='cubic')
        res.append(f(missingFiber))
        if doPlot:
            fig = plt.figure()
            plt.plot(x, coeff[:, i], 'o-')
            plt.plot(fibs, f(fibs), 'o-')

    return np.array(res)
