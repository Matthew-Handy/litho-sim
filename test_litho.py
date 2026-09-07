"""
test_litho.py -- automated checks, lock in every result validated by hand

Run with:   pytest -v
"""

import numpy as np
import pytest

from mask import (make_grid, single_line, line_space, contact_hole, measured_width, N, PIXEL)
from fourier import freq_grid, spectrum
from optics import aerial_image, pupil, WAVELENGTH, NA
from metrics import center_cut, measure_cd, contrast

X, Y = make_grid()


#masks
def test_mask_widths_are_exact():
    """The off-by-one fix from Step 1 stays fixed"""
    assert measured_width(single_line(X, 90)) == pytest.approx(90.0)
    assert measured_width(contact_hole(X, Y, 80)) == pytest.approx(80.0)


def test_line_space_is_fifty_percent_duty():
    assert line_space(X, 40).mean() == pytest.approx(0.5, abs=0.02)


#Fourier machinery
def test_perfect_lens_returns_the_mask():
    """A pupil that filters nothing has to reproduce the mask exactly"""
    m = single_line(X, 90)
    assert np.abs(aerial_image(m, na=100.0) - m).max() < 1e-20


def test_slit_spectrum_matches_analytic_sinc():
    width = 90.0
    FX, FY = freq_grid()
    fx = FX[N // 2, :]
    sim = np.abs(spectrum(single_line(X, width))[N // 2, :])
    sim = sim / sim.max()
    ana = np.abs(np.sinc(width * fx))
    band = np.abs(fx) < 0.03
    assert np.abs(sim[band] - ana[band]).max() < 0.02


def test_even_grating_orders_cancel():
    """A 50% duty grating suppresses its even diffraction orders."""
    half_pitch = 40.0                 # 80 nm pitch tiles the window 
    pitch = 2 * half_pitch
    FX, FY = freq_grid()
    fx = FX[N // 2, :]
    amp = np.abs(spectrum(line_space(X, half_pitch))[N // 2, :])
    amp = amp / amp.max()

    assert np.interp(1 / pitch, fx, amp) > 0.1      # 1st order present
    assert np.interp(2 / pitch, fx, amp) < 1e-6     # 2nd order cancelled
    assert np.interp(3 / pitch, fx, amp) > 0.01     # 3rd order present


#the lens
def test_pupil_cuts_off_at_na_over_lambda():
    FX, FY = freq_grid()
    P = np.abs(pupil(FX, FY, WAVELENGTH, NA))
    f = np.sqrt(FX**2 + FY**2)
    f_max = NA / WAVELENGTH
    assert P[f <= f_max * 0.99].min() == pytest.approx(1.0)
    assert P[f >= f_max * 1.01].max() == pytest.approx(0.0)


def test_defocus_is_symmetric():
    """With no aberrations, +z and -z must give identical images."""
    m = line_space(X, 120)
    assert np.abs(aerial_image(m, focus=+150.0)
                  - aerial_image(m, focus=-150.0)).max() < 1e-12


def test_defocus_reduces_contrast():
    m = line_space(X, 160)
    assert contrast(aerial_image(m, focus=0.0)) > \
           contrast(aerial_image(m, focus=200.0))


#metrics
def test_cd_shrinks_as_threshold_rises():
    x, prof = center_cut(aerial_image(single_line(X, 90)), X)
    cds = [measure_cd(x, prof, t) for t in (0.20, 0.25, 0.30, 0.35)]
    assert all(cds[k] > cds[k + 1] for k in range(len(cds) - 1))


#physics that cant be violated
def test_resolution_never_beats_the_diffraction_limit():
    from resolution import min_half_pitch
    for na in (0.60, 0.93, 1.35):
        measured = min_half_pitch(X, WAVELENGTH, na)
        assert measured >= WAVELENGTH / (2 * na) - PIXEL


def test_off_axis_beats_the_coherent_limit():
    """45 nm half-pitch fails on-axis at NA 1.35 but prints with a tilt."""
    from illumination import best_sigma, IMMERSION_NA, IMMERSION_N
    hp = 45.0
    m = line_space(X, hp)
    s = best_sigma(2 * hp, WAVELENGTH, IMMERSION_NA)
    on = contrast(aerial_image(m, na=IMMERSION_NA, n_medium=IMMERSION_N))
    off = contrast(aerial_image(m, na=IMMERSION_NA, n_medium=IMMERSION_N,
                                sigma_x=s))
    assert on < 0.2
    assert off > 0.8