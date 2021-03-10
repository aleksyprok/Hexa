import numpy as np
from scipy.io import FortranFile

dir = '../hexaf90/run1/'

# Get nx, ny, nz
f = open(dir + 'param1')
num_hex_cells = int(f.readline())
nx_ny_nz_string = f.readline()
nx = int(nx_ny_nz_string.split()[0])
ny = int(nx_ny_nz_string.split()[1])
nz = int(nx_ny_nz_string.split()[2])
print(nx)
print(ny)
print(nz)

# Get aax, aay, aaz
root = 'run1_00001p'
f = FortranFile( dir + root)
opt = f.read_ints()[0]
aax = f.read_reals(dtype = 'float32').reshape((nz + 1, ny + 1, nx    ))
aay = f.read_reals(dtype = 'float32').reshape((nz + 1, ny    , nx + 1))
aaz = f.read_reals(dtype = 'float32').reshape((nz    , ny + 1, nx + 1))
