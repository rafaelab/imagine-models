import ImagineModels as im
import numpy as np

jf12 = im.JF12MagneticField()
s = jf12.evaluate_grid([-1.*3.0856775806e+21, 0.*3.0856775806e+21, 1.*3.0856775806e+21],
                         [1.*3.0856775806e+21, 0.*3.0856775806e+21, -.1*3.0856775806e+21],
                         [.1*3.0856775806e+21, 0.*3.0856775806e+21, -.1*3.0856775806e+21])
print(np.asarray(s))
print(np.asarray(s).shape)
