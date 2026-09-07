"""
metrics.py -- turning an aerial image into numbers: edges, CD, contrast, NILS
"""

import numpy as np
import matplotlib.pyplot as plt

from mask import make_grid, single_line, line_space, N, PIXEL
from optics import aerial_image


#extracting 1D profile
def center_cut(I, X):
    """Intensity profile and x coordinates along the middle row"""
    row = I.shape[0] // 2
    return X[row, :], I[row, :]


#subpixel edges
def find_edges(x, I, threshold):
    """Subpixel x positions where the profile crosses 'threshold'"""
    above = I >= threshold
    idx = np.nonzero(np.diff(above))[0]     #where the state flips

    edges = []
    for i in idx:
        I1, I2 = I[i], I[i + 1]
        x1, x2 = x[i], x[i + 1]
        frac = (threshold - I1) / (I2 - I1)     #linear interpolation
        edges.append(x1 + frac * (x2 - x1))
    return np.array(edges)


#metrics
def measure_cd(x, I, threshold, center=0.0):
    """Width of the bright feature containing 'center' in nm"""
    edges = find_edges(x, I, threshold)
    left = edges[edges < center]
    right = edges[edges > center]
    if len(left) == 0 or len(right) == 0:
        return np.nan                       #nothing printed at this dose
    return right.min() - left.max()

def contrast(I, margin=0.25):
    """Michelson contrast over central region and avoiding wrap around edges"""
    n = I.shape[0]
    lo, hi = int(n * margin), int(n * (1 - margin))
    core = I[lo:hi, lo:hi]
    Imax, Imin = core.max(), core.min()
    if Imax + Imin == 0:
        return 0.0
    return (Imax - Imin) / (Imax + Imin)

def nils(x, I, threshold, cd, center=0.0):
    """Normalized image log slope at the left edge of the feature"""
    edges = find_edges(x, I, threshold)
    left = edges[edges < center]
    if len(left) == 0 or not np.isfinite(cd):
        return np.nan
    xe = left.max()
    dIdx = np.gradient(I, x)                #numerical derivative
    slope = np.interp(xe, x, dIdx)          #its value at subpixel edge
    return cd * abs(slope) / threshold


#display and sweeps
def show_cd(mask, X, threshold=0.30, focus=0.0, title="", zoom=300):
    I = aerial_image(mask, focus=focus)
    x, prof = center_cut(I, X)
    edges = find_edges(x, prof, threshold)
    cd = measure_cd(x, prof, threshold)

    fig, ax = plt.subplots()
    ax.plot(x, mask[N // 2, :], "k--", lw=1, alpha=0.4, label="mask")
    ax.plot(x, prof, label="aerial image")
    ax.axhline(threshold, color="red", ls=":", label=f"threshold = {threshold}")
    for e in edges:
        ax.axvline(e, color="green", ls=":", lw=0.8)
    ax.set_xlim(-zoom, zoom)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("intensity")
    ax.set_title(f"{title}  CD = {cd:.1f} nm")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    plt.show


def dose_sweep(mask, X, thresholds=(0.15, 0.20, 0.25, 0.30, 0.35), focus=0.0):
    """Print CD and NILS as the threshold (i.e. the dose) is varied"""
    x, prof = center_cut(aerial_image(mask, focus=focus), X)
    print(f"{'threshold':>10} {'CD (nm)':>10} {'NILS':>8}")
    for t in thresholds:
        cd = measure_cd(x, prof, t)
        print(f"{t:10.2f} {cd:10.1f} {nils(x, prof, t, cd):8.2f}")

if __name__ == "__main__":
    X, Y = make_grid()

    m = single_line(X, 90)
    show_cd(m, X, threshold=0.30, title="90 nm line at best focus")

    print("\n90 nm isolated line, best focus:")
    dose_sweep(m, X)

    print("\n90 nm isolated line, 150 nm defocus")
    dose_sweep(m, X, focus=150.0)

    g = line_space(X, 128)
    print("\ncontrast, 128 nm half pitch:")
    print("  best focus:  ", round(contrast(aerial_image(g)), 3))
    print("  140.7  nm defocus:", round(contrast(aerial_image(g, focus=140.7)), 3))