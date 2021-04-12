import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import os
import time

start_time = time.time()

output_dir = 'Figures/Video/x_y_plane'
os.makedirs(output_dir, exist_ok = True)

nx = 128
ny = 128
nz = 128
dx = np.float32(6 / nx)
dy = np.float32(6 / ny)
dz = np.float32(6 / nz)

fig = plt.figure()
fig_size = fig.get_size_inches()
fig_size[0] = fig_size[0] * 3
fig_size[1] = fig_size[1] * 2
fig.set_size_inches(fig_size)

n = 0
filename = 'run1/relax_' + '{:05d}'.format(n)
file = FortranFile(filename, 'r')
bbx = file.read_reals('float32').reshape((nz + 2, ny + 2, nx + 1), order = "C")
bby = file.read_reals('float32').reshape((nz + 2, ny + 1, nx + 2), order = "C")
bbz = file.read_reals('float32').reshape((nz + 1, ny + 2, nx + 2), order = "C")

for n in range(149):

    print(n)

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    bbx = file.read_reals('float32').reshape((nz + 2, ny + 2, nx + 1), order = "C")
    bby = file.read_reals('float32').reshape((nz + 2, ny + 1, nx + 2), order = "C")
    bbz = file.read_reals('float32').reshape((nz + 1, ny + 2, nx + 2), order = "C")

    ax = fig.add_subplot(331)
    im = ax.imshow(bbx[0, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('bbx, z=-1/2')

    ax = fig.add_subplot(332)
    im = ax.imshow(bby[0, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('bby, z=-1/2')

    ax = fig.add_subplot(333)
    im = ax.imshow(bbz[0, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('bbz, z=0')

    ax = fig.add_subplot(334)
    im = ax.imshow(bbx[1, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('bbx, z=1/2')

    ax = fig.add_subplot(335)
    im = ax.imshow(bby[1, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('bby, z=1/2')

    ax = fig.add_subplot(336)
    im = ax.imshow(bbz[1, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('bbz, z=1')

    ax = fig.add_subplot(337)
    im = ax.imshow(bbx[2, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('bbx, z=3/2')

    ax = fig.add_subplot(338)
    im = ax.imshow(bby[2, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('bby, z=3/2')

    ax = fig.add_subplot(339)
    im = ax.imshow(bbz[2, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('bbz, z=2')

    fig.savefig(output_dir + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()

output_dir = 'Figures/Video/current_x_y_plane'
os.makedirs(output_dir, exist_ok = True)

for n in range(149):

    print(n)

    filename = 'run1/relax_' + '{:05d}'.format(n)
    file = FortranFile(filename, 'r')
    bbx = file.read_reals('float32').reshape((nz + 2, ny + 2, nx + 1), order = "C")
    bby = file.read_reals('float32').reshape((nz + 2, ny + 1, nx + 2), order = "C")
    bbz = file.read_reals('float32').reshape((nz + 1, ny + 2, nx + 2), order = "C")

    ccx = (bbz[:, 1:, :] - bbz[:, :-1, :]) / dy \
        - (bby[1:, :, :] - bby[:-1, :, :]) / dz
    ccy = (bbx[1:, :, :] - bbx[:-1, :, :]) / dz \
        - (bbz[:, :, 1:] - bbz[:, :, :-1]) / dx
    ccz = (bby[:, :, 1:] - bby[:, :, :-1]) / dx \
        - (bbx[:, 1:, :] - bbx[:, :-1, :]) / dy

    ax = fig.add_subplot(331)
    im = ax.imshow(ccx[0, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('ccx, z=0')

    ax = fig.add_subplot(332)
    im = ax.imshow(ccy[0, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('ccy, z=0')

    ax = fig.add_subplot(333)
    im = ax.imshow(ccz[0, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('ccz, z=0.5')

    ax = fig.add_subplot(334)
    im = ax.imshow(ccx[1, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('ccx, z=1')

    ax = fig.add_subplot(335)
    im = ax.imshow(ccy[1, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('ccy, z=1')

    ax = fig.add_subplot(336)
    im = ax.imshow(ccz[1, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('ccz, z=1.5')

    ax = fig.add_subplot(337)
    im = ax.imshow(ccx[2, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('ccx, z=2')

    ax = fig.add_subplot(338)
    im = ax.imshow(ccy[2, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('ccy, z=2')

    ax = fig.add_subplot(339)
    im = ax.imshow(ccz[2, :, :], origin='lower')
    cb = fig.colorbar(im)
    ax.set_title('ccz, z=2.5')

    fig.savefig(output_dir + '/' + '{:05d}'.format(n) + '.png',  bbox_inches='tight')
    fig.clf()

print("--- %s seconds ---" % (time.time() - start_time))
