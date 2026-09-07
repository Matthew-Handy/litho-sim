"""
optics.py -- the projection lens: pupil filtering and aerial image
"""

import numpy as np
import matplotlib.pyplot as plt

from mask import make_grid, single_line, line_space, N, PIXEL
from fourier import freq_grid, spectrum

#defaults model to a 193nm ArF dry scanner
WAVELENGTH = 193.0
NA = 0.93


#the lens
def defocus_phase(FX, FY, focus=0.0, wavelength=WAVELENGTH, n_medium=1.0):
    """Phase error from placing the wafer `focus` nm away from best focus."""
    sin2 = (FX**2 + FY**2) * (wavelength / n_medium)**2
    sin2 = np.clip(sin2, 0.0, 1.0)          # keeps sqrt real outside pupil
    cos_theta = np.sqrt(1.0 - sin2)
    return (2 * np.pi * n_medium * focus / wavelength) * (1.0 - cos_theta)

def pupil(FX, FY, wavelength=WAVELENGTH, na=NA, focus=0.0, n_medium=1.0):
    """Complex pupil: a hard circular aperture times a defocus phase."""
    f_max = na / wavelength
    aperture = (np.sqrt(FX**2 + FY**2) <= f_max).astype(float)
    phase = defocus_phase(FX, FY, focus, wavelength, n_medium)
    return aperture * np.exp(1j * phase)


#image chain
def aerial_image(mask, wavelength=WAVELENGTH, na=NA, focus=0.0, n=N, pixel=PIXEL, n_medium=1.0):
    """Coherent aerial image at a given focus offset (nm)."""
    FX, FY = freq_grid(n, pixel)
    F = spectrum(mask)
    P = pupil(FX, FY, wavelength, na, focus, n_medium)
    E = np.fft.ifft2(np.fft.ifftshift(F * P))
    return np.abs(E) ** 2


#Display
def show_image(I, X, title=""):
    extent = [X.min(), X.max(), X.min(), X.max()]
    fig, ax = plt.subplots()
    im = ax.imshow(I, cmap="inferno", origin="lower", extent=extent)
    ax.set_title(title)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    fig.colorbar(im, ax=ax, label="intensity")
    plt.show()

def compare_cut(mask, I, X, title="", zoom=400):
    """Mask vs aerial image along the middle row"""
    row = N // 2
    x = X[row, :]
    fig, ax = plt.subplots()
    ax.plot(x, mask[row, :], label="mask (what you drew)")
    ax.plot(x, I[row, :], label="aerial image (what lands)")
    ax.set_xlim(-zoom, zoom)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("intensity")
    ax.set_title(title)
    ax.legend()
    ax.grid(alpha=0.3)
    plt.show()

def focus_sweep(mask, X, focus_list=(0, 60, 119, 180, 238), title="", zoom=500):
    """Overlay aerial images taken at several focus offsets."""
    row = N // 2
    x = X[row, :]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, mask[row, :], "k--", lw=1, alpha=0.5, label="mask")
    for z in focus_list:
        I = aerial_image(mask, focus=z)
        ax.plot(x, I[row, :], label=f"focus = {z:.0f} nm")
    ax.set_xlim(-zoom, zoom)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("intensity")
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    plt.show()

if __name__ == "__main__":
    X, Y = make_grid()

    #sanity check. image = mask since absurd NA passes everything
    m = single_line(X, 90)
    I_perfect = aerial_image(m, na=100.0)
    print("perfect-lens max error:", np.abs(I_perfect - m).max())
    
    #193nm dry scanner
    I = aerial_image(m)
    show_image(I, X, "Aerial image: 90 nm line, 193 nm, NA 0.93")
    compare_cut(m, I, X, "90 nm Isolated line")

    # grating too fine for lens
    g = line_space(X, 40)
    compare_cut(g, aerial_image(g), X, "45 nm half-pitch grating", zoom=300)


    #grating coarse enough to live
    g2 = line_space(X, 120)
    compare_cut(g2, aerial_image(g2), X, "120 nm half-pitch grating", zoom=700)

    focus_sweep(line_space(X, 120), X, title="120 nm half-pitch through focus", zoom=500)
    focus_sweep(single_line(X, 90), X, title="90 nm isolated line through focus", zoom=400)


