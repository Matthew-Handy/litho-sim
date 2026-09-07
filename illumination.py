"""
illumination.py -- off-axis illumination: beating the coherent resolution limit
"""

import numpy as np
import matplotlib.pyplot as plt

from mask import make_grid, line_space, N
from optics import aerial_image, WAVELENGTH
from metrics import contrast

IMMERSION_NA = 1.35     #193nm waterimmersion scanner
IMMERSION_N = 1.44

def best_sigma(pitch, wavelength=WAVELENGTH, na=IMMERSION_NA):
    """Illumination offset that puts the 0th and the +1st orders at opposite pupil edges. Returns a fraction of the pupil radius"""
    return -wavelength / (2.0 * pitch * na)

def compare(X, half_pitch=45.0, na=IMMERSION_NA):
    """Same grating, same lens, on-axis versus optimally tilted illumination"""
    m = line_space(X, half_pitch)
    s = best_sigma(2.0 * half_pitch, na=na)

    on = aerial_image(m, na=na, n_medium=IMMERSION_N)
    off = aerial_image(m, na=na, n_medium=IMMERSION_N, sigma_x=s)

    row = N // 2
    x = X[row, :]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x, m[row, :], "k--", lw=1, alpha=0.4, label="mask")
    ax.plot(x, on[row, :], label=f"on-axis  (contrast {contrast(on):.3f})")
    ax.plot(x, off[row, :], label=f"off-axis, sigma = {s:.2f}   (contrast {contrast(off):.3f})")
    ax.set_xlim(-300, 300)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("intensity")
    ax.set_title(f"{half_pitch:.0f} nm half-pitch at {WAVELENGTH:.0f} nm, " f"NA {na}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    plt.show()

def table(X, na=IMMERSION_NA, wavelength=WAVELENGTH):
    """Contrast versus half-pitch for both illumination schemes"""
    print(f"lambda = {wavelength:.0f} nm, NA = {na}")
    print(f"    on-axis limit   lambda/(2 NA) = {wavelength / (2*na):.1f} nm")
    print(f"    off-axis limit  lambda/(4 NA) = {wavelength / (4*na):.1f} nm\n")
    print(f"{'half-pitch':>11} {'sigma':>7} {'on-axis':>9} {'off-axis':>9}")

    for hp in (30, 35, 40, 45, 50, 55, 70, 75, 90):
        m = line_space(X, float(hp))
        s = best_sigma(2.0 * hp, wavelength, na)
        on = contrast(aerial_image(m, wavelength=wavelength, na=na, n_medium=IMMERSION_N))

        if abs(s) > 1.0:
            print(f"{hp:11.0f} {s:7.2f} {on:9.3f}    tilt exceeds pupil")
        else:
            off = contrast(aerial_image(m, wavelength=wavelength, na=na, n_medium=IMMERSION_N, sigma_x=s))
            print(f"{hp:11.0f} {s:7.2f} {on:9.3f} {off:9.3f}")


if __name__ == "__main__":
    X, Y = make_grid()
    table(X)
    compare(X, half_pitch=45.0)
