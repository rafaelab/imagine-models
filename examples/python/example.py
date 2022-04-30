import ImagineModels as im
import numpy as np


# load model

jf12 = im.JF12MagneticField()

# define a position in cartesian Galactic coordinates (kpc)

pos = [1., -3., 0.1]

# evaluate the model at the position

b_at_pos = jf12.evaluate_pos(pos)

# define grid axes in cartesian Galactic coordinates (kpc)

grid_x = [-3., 0., 1.]
grid_y = [-1., 1.5, 2.]
grid_z = [.1, 0., -1.]

# evaluate grid
b_grid = jf12.evaluate_grid(grid_x, grid_y, grid_z)

# NOTE: for further processing, the evalauted grid should be converted to a numpy array, this may become unnecessary in
# the future

b_grid = np.asarray(b_grid)

# change a parameter

jf12.disk_amp = 24

# evaluate model at position with updated parameter
b_at_pos = jf12.evaluate_pos(pos)
