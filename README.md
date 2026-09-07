# Computational Lithography Simulator

A physically motivated model of optical projection lithography, built from scratch in
Python. It propagates light from a photomask through a projection lens using Fourier
optics, computes the aerial image that lands on the wafer, and extracts the process
metrics a fab actually cares about — critical dimension, image contrast, NILS, depth
of focus, and the focus-exposure process window.

The goal is not to reproduce an industrial scanner. It is to build a model whose
behavior is validated against analytic optical theory, and to use it to show *why*
semiconductor lithography gets harder as features shrink.

**Validated against theory in three independent places:** single-slit sinc
diffraction, the $\lambda/2\mathrm{NA}$ resolution limit, and the Rayleigh depth of
focus (200 nm measured vs 223 nm predicted).

---

## The physics

A photomask is a stencil. Light passes through it, a lens collects that light, and it
is projected onto a wafer. If lenses were perfect the wafer would receive a perfect
copy — but when light passes through small openings it diffracts into a fan of angles,
and a lens is a finite piece of glass that can only collect angles up to some maximum.
Light at steeper angles misses the lens entirely and is lost. The image is rebuilt from
only what survived.

The whole simulation is that sentence, in four steps:

```
mask  --FFT-->  angular spectrum  --x pupil-->  filtered  --inverse FFT-->  E field  --|E|^2-->  intensity
```

Each spatial frequency in the mask corresponds to light leaving at a physical angle,

$$\sin\theta = \lambda f$$

so the lens, with $\mathrm{NA} = n\sin\theta_{max}$, acts as a **low-pass filter** with
a hard cutoff at

$$f_{max} = \frac{\mathrm{NA}}{\lambda}$$

Everything finer than that never reaches the wafer. Defocus is modeled as a phase term
on the pupil — light at different angles travels different path lengths to a displaced
plane, so it arrives out of step:

$$\Phi(f) = \frac{2\pi z}{\lambda}\left(1 - \sqrt{1 - (\lambda f)^2}\right)$$

This is the exact scalar form, not the paraxial approximation $\Phi \approx \pi\lambda z f^2$,
which is off by nearly a factor of two at NA 0.93.

Photoresist is modeled as a constant threshold: material clears where intensity exceeds
$I_{th}$. Exposure dose scales the intensity, so dose and threshold are the same knob.

**Assumptions:** scalar diffraction, coherent on-axis illumination, aberration-free
lens, binary chrome-on-glass masks, no resist diffusion or development kinetics.

---

## Results

### Sharp mask edges cannot survive the lens

![mask vs aerial image](figures/mask_vs_aerial.png)

A perfectly rectangular 90 nm mask feature becomes a rounded hump peaking at 0.63 with
sloped sides and faint side lobes. Reproducing a vertical wall requires infinitely high
spatial frequencies; the pupil discards everything past 0.0048 cycles/nm. The ringing is
the Gibbs phenomenon, and in lithography those side lobes are a real defect mechanism.

### Validation 1 — single-slit diffraction

![sinc validation](figures/sinc_validation.png)

The FFT of a 90 nm slit against the closed-form result $w\,\mathrm{sinc}(wf)$. Building
this comparison first caught a real bug: an inclusive-endpoint comparison was making the
mask 95 nm wide instead of 90, which showed up as progressively drifting side lobes.

### Contrast reversal through focus

![focus sweep](figures/focus_sweep.png)

A 240 nm-pitch grating imaged at five focus offsets. Only three diffraction orders
survive the pupil, so the intensity is

$$I(x) = 0.453 + 0.637\cos\Phi\cos(kx) + 0.203\cos(2kx)$$

At $z = 119$ nm the defocus phase reaches $\pi/2$, the $\cos\Phi$ term vanishes, and the
pattern prints at **twice the intended pitch**. At $z = 238$ nm it returns at full
contrast, completely **inverted**. Both effects fall out of the model; neither was
programmed in.

### Validation 2 — the resolution limit

![resolution vs NA](figures/resolution_vs_na.png)

Minimum resolvable half-pitch measured by shrinking the pitch until image contrast
collapses, swept across NA and compared against $\lambda/2\mathrm{NA}$. Every point sits
on or above the theory line — a point below would mean the simulator was beating the
diffraction limit, which is the validation criterion for this figure.

This confirms $k_1 = 0.5$, the hard floor for coherent on-axis illumination. Production
193 nm immersion tools reach $k_1 \approx 0.28$, which is only possible because they do
not use on-axis coherent light (see Roadmap).

### Bossung curves and the process window

![bossung curves](figures/bossung.png)

![process window](figures/process_window.png)

Printed CD across a focus-exposure matrix. The curves are symmetric parabolas about best
focus, as they must be for an aberration-free lens — asymmetry in a real fab indicates
spherical aberration or coma. The green region marks where CD stays within ±10% of the
90 nm target.

### Validation 3 — depth of focus

Largest contiguous in-spec focus range: **200 nm**, against the Rayleigh estimate
$\lambda/\mathrm{NA}^2 = 223$ nm. The focus axis is sampled every 20 nm, so the measured
value is quantized to multiples of 20; and $k_2$ is a convention rather than a derived
constant. Agreement within 10% is about as close as this comparison can meaningfully get.

---

## Known limitations

Deviations from theory in this model are understood and quantified rather than ignored.

**Spectral leakage.** The DFT treats the field as periodic. When the simulation window
does not contain a whole number of grating periods, diffraction orders smear across bins
and energy that should be blocked leaks through the pupil. This is why the contrast
criterion in `resolution.py` is 0.99 rather than 0.5: leakage can produce up to ~0.976
contrast, while a genuinely resolved three-beam image is pinned at exactly 1.0 because
the field crosses zero. Pitches that tile the window exactly (80, 160, 320 nm) are free
of this.

**Duty-cycle quantization.** A pitch that is not an integer number of pixels cannot have
exactly 50% duty cycle, which shifts the harmonic amplitudes slightly.

**Residual scatter in the resolution curve.** Frequency bins are spaced $1/2560$
cycles/nm, and half-pitch relates to bin index as $p = 1280/b$, so one bin of ambiguity
costs $\Delta p \approx 1280/b^2$ nm — about ±36 nm at NA 0.45, ±4 nm at NA 1.35. This
fully accounts for the scatter, and explains why the fit tightens at high NA. The
principled fix is to size the simulation window to an integer number of periods for each
pitch tested.

**Nyquist floor.** At 5 nm pixels the finest representable period is 10 nm, so EUV at
13.5 nm cannot be simulated meaningfully on this grid.

**Physics not yet modeled:** partial coherence, off-axis illumination, lens aberrations,
phase-shifting masks, vector effects at high NA, resist diffusion and development
kinetics, and stochastic photon shot noise.

---

## Repository layout

| File | Contents |
|---|---|
| `mask.py` | Simulation grid and mask transmission functions |
| `fourier.py` | FFT propagation, frequency coordinates, analytic validation |
| `optics.py` | Pupil function, NA, wavelength, defocus, aerial image |
| `metrics.py` | Threshold resist, sub-pixel edge detection, CD, contrast, NILS |
| `resolution.py` | Resolution vs NA and vs wavelength sweeps |
| `process_window.py` | Focus-exposure matrix, Bossung curves, process window |
| `test_litho.py` | pytest suite validating the optics against theory |
| `make_figures.py` | Regenerates every figure in this README |

## Running it

```bash
git clone https://github.com/<your-username>/litho-sim.git
cd litho-sim
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy matplotlib pytest

pytest -v                    # 10 checks against analytic theory
python make_figures.py       # regenerate all figures
python process_window.py     # focus-exposure matrix and process window
```

## Roadmap

- **Off-axis illumination** — tilt the source so the 0th and 1st orders straddle the
  pupil, doubling the resolvable frequency and pushing $k_1$ toward 0.25. This is how
  193 nm tools printed 38 nm half-pitch, and it should let the simulator break the
  $k_1 = 0.5$ limit it currently proves.
- **Partial coherence** — sum over an extended source rather than a single plane wave.
- **Phase-shifting masks** — allow complex mask transmission with 0 and $\pi$ regions.
- **Lens aberrations** — Zernike terms on the pupil; spherical aberration should tilt
  the Bossung curves asymmetrically.
- **Stochastic effects** — photon shot noise and Monte Carlo CD distributions.
- **Interactive dashboard** — live parameter exploration.

## Built with

Python, NumPy, SciPy, Matplotlib, pytest.