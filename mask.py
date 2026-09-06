"""
mask.py -- builds 2D photomask patterns for lithography simulator.

Every mask is a transmission function: 1.0 = light passes, 0.0 = blocked.
"""

import numpy as np
import matplotlib.pyplot as plt


#Sim grid
N = 512  # Number of pixels as N x N
PIXEL = 5.0 #nm rep by 1 pixel

def make_grid(n=N , pixel=PIXEL):
   """Return X,Y cord array in nm, centered at zero """
   coords = (np.arange(n) - n // 2) * pixel
   X, Y = np.meshgrid(coords, coords)
   return X, Y


#mask shapes
def single_line(X,width):
  """clear line of given width (nm) running along y"""
  return (np.abs(X) <= width / 2).astype(float)

def line_space(X, half_pitch):
  """repeating equal lines && spaces. half_pitch = line width = space width"""
  pitch = 2.0 * half_pitch
  return (np.mod(X + half_pitch / 2, pitch) < half_pitch).astype(float)

def contact_hole(X, Y, size):
    """Square opening of the given side length (nm)."""
    return ((np.abs(X) <= size / 2) & (np.abs(Y) <= size / 2)).astype(float)


#display helper
def show(mask, X, title=""):
   """Draw mask with real nm axes"""
   extent = [X.min(), X.max(), X.min(), X.max()]
   plt.imshow(mask, cmap="gray", origin="lower", extent=extent)
   plt.title(title)
   plt.xlabel("x (nm)")
   plt.ylabel("y (nm)")
   plt.colorbar(label="transmission")
   plt.show()


if __name__ == "__main__":
   X, Y = make_grid()
   show(single_line(X, 90), X, "Isolated 90nm line")
   show(line_space(X,45), X, "45nm half-pitch lines/spaces")
   show(contact_hole(X, Y, 80), X, "80 nm contact hole")

