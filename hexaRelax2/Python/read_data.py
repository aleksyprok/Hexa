import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

nx = 128
ny = 128
nz = 128
dx = 6 / nx
dy = 6 / ny
dz = 6 / nz

filename = 'run1/relax_00000'
file = FortranFile(filename, 'r')
bbx = file.read_reals('float32').reshape((nz + 2, ny + 2, nx + 1), order = "C")
bby = file.read_reals('float32').reshape((nz + 2, ny + 1, nx + 2), order = "C")
bbz = file.read_reals('float32').reshape((nz + 1, ny + 2, nx + 2), order = "C")

filename = 'run1/poten_00003p'
file = FortranFile(filename, 'r')
opt = file.read_ints()
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
aaz = file.read_reals('float32').reshape((nz, ny + 1, nx + 1), order = "C")

bbx_p = (aaz[ :, 1:, :] - aaz[  :, :-1, :]) / dy \
      - (aay[1:,  :, :] - aay[:-1,   :, :]) / dz
bby_p = (aax[1:, :,  :] - aax[:-1, :,   :]) / dz \
      - (aaz[ :, :, 1:] - aaz[  :, :, :-1]) / dx
bbz_p = (aay[:,  :, 1:] - aay[0,   :, :-1]) / dx \
      - (aax[:, 1:,  :] - aax[0, :-1,   :]) / dy

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 3
fig.set_size_inches(fig_size)

ax = fig.add_subplot(331)
im = ax.imshow(bbx[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(332)
im = ax.imshow(bby[1, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(333)
im = ax.imshow(bbz[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(334)
im = ax.imshow(bbx_p[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(335)
im = ax.imshow(bby_p[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(336)
im = ax.imshow(bbz_p[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(337)
im = ax.imshow(bbx[1, 1:-1, :] - bbx_p[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(338)
im = ax.imshow(bby[1, :, 1:-1] - bby_p[0, :, :], origin='lower')
cb = fig.colorbar(im)

ax = fig.add_subplot(339)
im = ax.imshow(bbz[0, 1:-1, 1:-1] - bbz_p[0, :, :], origin='lower')
cb = fig.colorbar(im)


plt.show()
