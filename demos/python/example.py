import ImagineModels as im
import numpy as np

x, y, z = [1., -3., 0.1]

# irregular
grid_x = np.asarray([-3., 0., 1.])
grid_y = np.asarray([-1., 1.5, 2.])
grid_z = np.asarray([.1, 0., -1.])

#regular

shape = [2, 2, 2]
increment = [2.1, 0.3, 1.]
zeropoint = [-4., 0.1, -0.3]


# load model

jf12regular = im.JF12RegularField(grid_x, grid_y, grid_z)
jf12random = im.JF12RandomField()
GaussRandom = im.GaussianScalarField(shape, zeropoint, increment)
# jf12random2 = im.JF12RandomField(shape, zeropoint, increment)
jaffe = im.JaffeMagneticField(grid_x, grid_y, grid_z)
# helix = im.HelixMagneticField()

# define a position in cartesian Galactic coordinates (kpc)



# evaluate the models at the position

# b_at_pos_jf12 = jf12regular.at_position(x, y, z)
# b_at_pos_jaffe = jaffe.at_position(x, y, z)
# b_at_pos_helix = helix.at_position(x, y, z)

# define grid axes in cartesian Galactic coordinates (kpc)



# jf12random2 = im.JF12RandomField(shape, zeropoint, increment)
# evaluate grid
# b_grid = jf12regular.on_grid(grid_x, grid_y, grid_z)




# change a parameter

#jf12regular.w_disk = 24

# print(jf12regular2.on_grid())

# evaluate model at position with updated parameter
#b_at_pos = jf12regular.at_position(x, y, z)

#print(jf12regular.on_grid(shape, zeropoint, increment))
#print(jf12regular2.on_grid()) #TODO: this line leads to bad malloc
# print(b_at_pos)
for i in range(1):
    xx = np.random.uniform(-4, 4, i + 3)
    yy = np.random.uniform(-4, 4, 2*i + 2)
    zz = np.random.uniform(-4, 4, 3*i + 1)
 #   print(xx, yy, zz)
    
 #   print(jf12regular.on_grid(xx, yy, zz))
    
# print(jf12random2.on_grid(23))
# print(jf12random.on_grid(shape, zeropoint, increment, 23))
gr = GaussRandom.on_grid(42)
print(gr)
print(gr.shape)
print(jaffe.on_grid())
print(jf12regular.on_grid())
