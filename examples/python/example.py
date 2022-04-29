import ImagineModels as im
import numpy as np

jf12 = im.JF12MagneticField()
s = jf12.evaluate_grid([-1., 0., 1.], [1., 0., -.1], [.1, 0., -.1])
print(np.asarray(s))
print(np.asarray(s).shape)
