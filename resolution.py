"""
resolution.py -- measuring the resolution limit rather than assuming it
"""

import numpy as np
import matplotlib.pyplot as plt

from mask import make_grid, line_space, PIXEL
from optics import aerial_image
from metrics import contrast


def min_half_pitch(X, wavelength, na, criterion=0.99, hp_min=10.0, hp_max=400.0):
    """Smallest half pitch that still images with contrast >= criterion"""
    for hp in np.arange(hp_min, hp_max, PIXEL):
        I = aerial_image(line_space(X, hp), wavelength=wavelength, na=na)
        if contrast(I) >= criterion:
            return hp
    return np.nan

def sweep_na(X, na_list=np.arange(0.40, 1.40, 0.05), wavelength=193.0):
    measured = np.array([min_half_pitch(X, wavelength, na) for na in na_list])
    theory = wavelength / (2.0 * na_list)

    fig, ax = plt.subplots()
    ax.plot(na_list, measured, "o", label="simulated")
    ax.plot(na_list, theory, "-", label=r"theory:   $\lambda / 2\,\mathrm{NA}$")
    ax.set_xlabel("numerical aperture")
    ax.set_ylabel("minimum half pitch (nm)")
    ax.set_title(f"Resolution vs NA at {wavelength:.0f} nm")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.show

def sweep_wavelength(X, wl_list=(436, 365, 248, 193, 157, 13.5), na=0.93):
    measured = np.array([min_half_pitch(X, wl, na) for wl in wl_list])
    wl = np.array(wl_list, dtype=float)
    theory = wl / (2.0 * na)

    fig, ax = plt.subplots()
    ax.plot(wl, measured, "o", label="simulated")
    ax.plot(wl, theory, "-", label=r"theory:    $\lambda / 2\,\mathrm{NA}$")
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("minimum half pitch (nm)")
    ax.set_title(f"Resolution vs wavelength at NA {na}")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.show()


if __name__ == "__main__":
    X, Y = make_grid()
    sweep_na(X)
    sweep_wavelength(X)


