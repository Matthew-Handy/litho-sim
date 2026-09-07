"""
process_window.py -- focus exposure matrix, Bossung curves, and the process window
"""

import numpy as np
import matplotlib.pyplot as plt

from mask import make_grid, single_line
from optics import aerial_image
from metrics import center_cut, measure_cd

TARGET_CD = 90.0        #nm 
THRESHOLD = 0.339       #resist clearing thresholf at nominal dose
TOLERANCE = 0.10        # +/- 10 % CD spec

def cd_matrix(mask, X, focus_list, dose_list):
    """CD at evry (focus, dose) pair. Rows are dose, columns are focus"""
    cd = np.full((len(dose_list), len(focus_list)), np.nan)
    for j, z in enumerate(focus_list):
        x, prof = center_cut(aerial_image(mask, focus=z), X)        #one fft pair
        for i, d in enumerate(dose_list):
            cd[i, j] = measure_cd(x, d * prof, THRESHOLD)           #easy inner loop
    return cd

def bossung(cd, focus_list, dose_list):
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, d in enumerate(dose_list):
        ax.plot(focus_list, cd[i, :], marker="o", ms=3, label=f"dose {d:.2f}")
    ax.axhline(TARGET_CD, color="k", ls="--", lw=1)
    ax.axhline(TARGET_CD * (1 - TOLERANCE), color="k", ls=":", lw=0.8)
    ax.axhline(TARGET_CD * (1 + TOLERANCE), color="k", ls=":", lw=0.8)
    ax.set_xlabel("focus (nm)")
    ax.set_ylabel("printed CD (nm)")
    ax.set_title("Bossung curves: CD through focus at several doses")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    plt.show()


def show_window(cd, focus_list, dose_list):
    lo = TARGET_CD * (1 - TOLERANCE)
    hi = TARGET_CD * (1 + TOLERANCE)
    ok = ((cd >= lo) & (cd <= hi)).astype(float)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.contourf(focus_list, dose_list, ok, levels=[0.5, 1.5],
                colors=["#4CAF50"], alpha=0.5)
    cs = ax.contour(focus_list, dose_list, cd,
                    levels=[lo, TARGET_CD, hi], colors="k", linewidths=1)
    ax.clabel(cs, fmt="%.0f nm")
    ax.set_xlabel("focus (nm)")
    ax.set_ylabel("relative dose")
    ax.set_title(f"Process window: CD = {TARGET_CD:.0f} nm "
                 f"+/- {TOLERANCE*100:.0f}%")
    ax.grid(alpha=0.3)
    plt.show()


def longest_in_spec(row_ok, focus_list):
    """Longest CONTIGUOUS in-spec focus range in this dose row, in nm."""
    best = 0.0
    start = None
    for k, good in enumerate(row_ok):
        if good and start is None:
            start = k                                        # a run begins
        if not good and start is not None:
            best = max(best, focus_list[k - 1] - focus_list[start])
            start = None                                     # a run ends
    if start is not None:                                    # run hits edge
        best = max(best, focus_list[-1] - focus_list[start])
    return best


def report(cd, focus_list, dose_list):
    lo = TARGET_CD * (1 - TOLERANCE)
    hi = TARGET_CD * (1 + TOLERANCE)
    ok = (cd >= lo) & (cd <= hi)

    best_span, best_dose = 0.0, np.nan
    for i, d in enumerate(dose_list):
        span = longest_in_spec(ok[i, :], focus_list)
        if span > best_span:
            best_span, best_dose = span, d

    print(f"target CD: {TARGET_CD:.0f} nm  +/- {TOLERANCE*100:.0f}%")
    print(f"largest depth of focus: {best_span:.0f} nm at dose {best_dose:.2f}")
    print(f"Rayleigh estimate k2*lambda/NA^2: {193/0.93**2:.0f} nm")

if __name__ == "__main__":
    X, Y = make_grid()
    m = single_line(X, 90)

    focus_list = np.arange(-200.0, 201.0, 20.0)
    dose_list = np.arange(0.85, 1.21, 0.025)

    x, prof = center_cut(aerial_image(m, focus=0.0), X)
    print("THRESHOLD =", THRESHOLD)
    print("peak intensity =", prof.max())
    print("x range =", x.min(), x.max())
    print("CD at dose 1.0 =", measure_cd(x, 1.0 * prof, THRESHOLD))

    cd = cd_matrix(m, X, focus_list, dose_list)
    bossung(cd, focus_list, dose_list)
    show_window(cd, focus_list, dose_list)
    report(cd, focus_list, dose_list) 