import ImagineModels as im
import numpy as np


# load model

jf12 = im.JF12RegularField()
jaffe = im.JaffeMagneticField()
helix = im.HelixMagneticField()

# define a position in cartesian Galactic coordinates (kpc)

x, y, z = [1., -3., 0.1]

# evaluate the models at the position

b_at_pos_jf12 = jf12.at_position(x, y, z)
b_at_pos_jaffe = jaffe.at_position(x, y, z)
b_at_pos_helix = helix.at_position(x, y, z)

# define grid axes in cartesian Galactic coordinates (kpc)

grid_x = np.asarray([-3., 0., 1.])
grid_y = np.asarray([-1., 1.5, 2.])
grid_z = np.asarray([.1, 0., -1.])

# evaluate grid
b_grid = jf12.on_grid(grid_x, grid_y, grid_z)


# change a parameter

# jf12.disk_amp = 24

# evaluate model at position with updated parameter
b_at_pos = jf12.at_position(x, y, z)

print(b_at_pos)
for i in range(5):
    xx = np.random.uniform(-4, 4, i + 3)
    yy = np.random.uniform(-4, 4, 2*i + 2)
    zz = np.random.uniform(-4, 4, 3*i)
    print(jf12.on_grid(xx, yy, zz))
