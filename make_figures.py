"""
make_figures.py -- regen every figure used in README

Writes PNGs into figures/.  Run with:  python make_figures.py
"""

import os
import matplotlib
matplotlib.use("Agg")            # render to files instead of opening windows
import matplotlib.pyplot as plt
import numpy as np

from mask import make_grid, single_line, line_space, N
from fourier import freq_grid, spectrum
from optics import aerial_image, WAVELENGTH, NA
from resolution import min_half_pitch
from process_window import (cd_matrix, longest_in_spec,
                            TARGET_CD, TOLERANCE)

OUT = "figures"
os.makedirs(OUT, exist_ok=True)

X, Y = make_grid()
row = N // 2


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def fig_mask_vs_aerial():
    m = single_line(X, 90)
    I = aerial_image(m)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(X[row, :], m[row, :], label="mask (90 nm line)")
    ax.plot(X[row, :], I[row, :], label="aerial image on wafer")
    ax.set_xlim(-400, 400)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("intensity")
    ax.set_title(f"193 nm, NA {NA}: sharp mask edges cannot survive the lens")
    ax.legend()
    ax.grid(alpha=0.3)
    save(fig, "mask_vs_aerial.png")


def fig_sinc_validation():
    width = 90.0
    FX, FY = freq_grid()
    fx = FX[row, :]
    sim = np.abs(spectrum(single_line(X, width))[row, :])
    sim = sim / sim.max()
    ana = np.abs(np.sinc(width * fx))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(fx, sim, label="FFT of mask")
    ax.plot(fx, ana, "--", label="analytic sinc")
    ax.set_xlim(-0.05, 0.05)
    ax.set_xlabel("spatial frequency (cycles/nm)")
    ax.set_ylabel("normalized |E|")
    ax.set_title("Validation 1: single-slit diffraction vs theory")
    ax.legend()
    ax.grid(alpha=0.3)
    save(fig, "sinc_validation.png")


def fig_focus_sweep():
    m = line_space(X, 120)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(X[row, :], m[row, :], "k--", lw=1, alpha=0.4, label="mask")
    for z in (0, 60, 119, 180, 238):
        ax.plot(X[row, :], aerial_image(m, focus=z)[row, :],
                label=f"focus = {z} nm")
    ax.set_xlim(-500, 500)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("intensity")
    ax.set_title("Contrast reversal through focus (240 nm pitch)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    save(fig, "focus_sweep.png")


def fig_process():
    m = single_line(X, 90)
    focus_list = np.arange(-200.0, 201.0, 20.0)
    dose_list = np.arange(0.85, 1.21, 0.025)
    cd = cd_matrix(m, X, focus_list, dose_list)

    lo = TARGET_CD * (1 - TOLERANCE)
    hi = TARGET_CD * (1 + TOLERANCE)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for i, d in enumerate(dose_list):
        ax.plot(focus_list, cd[i, :], marker="o", ms=3, label=f"{d:.2f}")
    ax.axhline(TARGET_CD, color="k", ls="--", lw=1)
    ax.axhline(lo, color="k", ls=":", lw=0.8)
    ax.axhline(hi, color="k", ls=":", lw=0.8)
    ax.set_xlabel("focus (nm)")
    ax.set_ylabel("printed CD (nm)")
    ax.set_title("Bossung curves: CD through focus at several doses")
    ax.legend(fontsize=6, ncol=3, title="dose")
    ax.grid(alpha=0.3)
    save(fig, "bossung.png")

    ok = ((cd >= lo) & (cd <= hi)).astype(float)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.contourf(focus_list, dose_list, ok, levels=[0.5, 1.5],
                colors=["#4CAF50"], alpha=0.5)
    cs = ax.contour(focus_list, dose_list, cd,
                    levels=[lo, TARGET_CD, hi], colors="k", linewidths=1)
    ax.clabel(cs, fmt="%.0f nm")
    ax.set_xlabel("focus (nm)")
    ax.set_ylabel("relative dose")
    ax.set_title(f"Process window: {TARGET_CD:.0f} nm +/- {TOLERANCE*100:.0f}%")
    ax.grid(alpha=0.3)
    save(fig, "process_window.png")

    best = max(longest_in_spec(ok[i, :] > 0.5, focus_list)
               for i in range(len(dose_list)))
    print(f"  depth of focus = {best:.0f} nm, "
          f"Rayleigh estimate = {WAVELENGTH / NA**2:.0f} nm")


def fig_resolution_vs_na():
    na_list = np.arange(0.40, 1.40, 0.05)
    measured = np.array([min_half_pitch(X, WAVELENGTH, na) for na in na_list])
    theory = WAVELENGTH / (2.0 * na_list)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(na_list, measured, "o", label="simulated")
    ax.plot(na_list, theory, "-", label=r"theory:  $\lambda / 2\,\mathrm{NA}$")
    ax.set_xlabel("numerical aperture")
    ax.set_ylabel("minimum half-pitch (nm)")
    ax.set_title("Validation 2: resolution limit at 193 nm")
    ax.legend()
    ax.grid(alpha=0.3)
    save(fig, "resolution_vs_na.png")


if __name__ == "__main__":
    fig_mask_vs_aerial()
    fig_sinc_validation()
    fig_focus_sweep()
    fig_process()
    fig_resolution_vs_na()      # kept for last since slow