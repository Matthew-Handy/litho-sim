"""
fourier.py -- moves patterns between real space and spatial frequency space.
"""

import numpy as np
import matplotlib.pyplot as plt

from mask import make_grid, single_line, line_space, N, PIXEL


#Freq cords
def freq_grid(n=N, pixel=PIXEL):
    """Spatial frequency cords in cycles per nm, centered on 0"""
    f= np.fft.fftshift(np.fft.fftfreq(n, d=pixel))
    FX, FY = np.meshgrid(f, f)
    return FX, FY


#Fourier transform
def spectrum(mask):
    """2D Fourier transform of a mask, zero freq moved to center"""
    return np.fft.fftshift(np.fft.fft2(mask))


#Display
def show_spectrum(F, FX, title="", zoom=0.05):
    """Draw spectrum on a log scale so weak orders stay visible"""
    extent = [FX.min(), FX.max(), FX.min(), FX.max()]
    plt.imshow(np.log10(1 + np.abs(F)), cmap="inferno", origin="lower", extent=extent)
    plt.title(title)
    plt.xlabel("fx (cycles/nm)")
    plt.ylabel("fy (cycles/nm)")
    plt.xlim(-zoom, zoom)
    plt.ylim(-zoom, zoom)
    plt.colorbar(label="log10(1 + |E|)")
    plt.show()


#validate against analytic theory
def sinc_check(width=90.0):
    """Compare the simulated slit spectrum to the analytic sinc func"""
    X, Y = make_grid()
    FX, FY = freq_grid()
    F = spectrum(single_line(X, width))

    row = N // 2        #fy = 0 cut through center
    fx = FX[row, :]
    simulated = np.abs(F[row, :])
    simulated = simulated / simulated.max()

    analytic = np.abs(np.sinc(width * fx))      # np.sinc(u) = sin(pi*u)/(pi*u)

    plt.plot(fx, simulated, label="FFT of mask")
    plt.plot(fx, analytic, "--", label="analytic sinc")
    plt.xlim(-0.05, 0.05)
    plt.xlabel("fx (cycles/nm)")
    plt.ylabel("normalized |E|")
    plt.title(f"Single {width: .0f} nm slit: simulation vs theory")
    plt.legend()
    plt.show()

def order_check(half_pitch=45.0):
    """Show a grating's diffraction orders against theory"""
    X, Y = make_grid()
    FX, FY = freq_grid()
    F = spectrum(line_space(X, half_pitch))

    row = N // 2
    fx = FX[row, :]
    amp = np.abs(F[row, :])
    amp = amp / amp.max()

    pitch = 2.0 * half_pitch
    fig, ax = plt.subplots()
    ax.semilogy(fx,amp + 1e-12, label="FFT of grating")
    for m in range(-4, 5):
        ax.axvline(m / pitch, color="red", ls=":", lw=0.8)
    ax.set_xlim(-0.05, 0.05)
    ax.set_ylim(1e-6, 2)
    ax.set_xlabel("fx (cycles/nm)")
    ax.set_ylabel("normalized |E| (log scale)")
    ax.set_title(f"Grating orders, predicted at multiples of 1/{pitch:.0f} nm")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    X, Y = make_grid()
    FX, FY = freq_grid()

    show_spectrum(spectrum(single_line(X, 90)), FX, "Spectrum: isolated 90nm line")
    show_spectrum(spectrum(line_space(X, 45)), FX, "Spectrum: 45nm half-pitch grating")
    sinc_check(90.0)
    order_check(45.0)


