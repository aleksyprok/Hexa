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
# aax = file.read_record(np.dtype(('float32', (nz + 1, ny + 1, nx))))
aax = file.read_reals('float32').reshape((nz + 1, ny + 1, nx), order = "C")
aay = file.read_reals('float32').reshape((nz + 1, ny, nx + 1), order = "C")
bbz0 = (aay[0,  :, 1:] - aay[0,   :, :-1]) / dx \
     - (aax[0, 1:,  :] - aax[0, :-1,   :]) / dy

fig = plt.figure()

ax = fig.add_subplot(121)
ax.imshow(bbz[0, :, :])

ax = fig.add_subplot(122)
ax.imshow(bbz0[:, :])

plt.show()
