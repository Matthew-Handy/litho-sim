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
def pupil(FX, FY, wavelength=WAVELENGTH, na=NA):
    """Circular low-pass filter : 1 inside the lens acceptance, 0 outside"""
    f_max = na / wavelength
    return (np.sqrt(FX**2 + FY**2) <= f_max).astype(float)


#image chain
def aerial_image(mask, wavelength=WAVELENGTH, na=NA, n=N, pixel=PIXEL):
    """Coherent aerial image: mask -> spectrum -> pupil -> back -> intensity"""
    FX, FY = freq_grid(n, pixel)
    F = spectrum(mask)      #centered spectrum
    P = pupil(FX, FY, wavelength, na)
    E = np.fft.ifft2(np.fft.ifftshift(F * P))       #undo shift then invert
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


