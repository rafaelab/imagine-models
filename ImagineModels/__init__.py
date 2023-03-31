from _ImagineModels import VectorFieldBase, ScalarFieldBase, RegularVectorField, RegularScalarField, \
    JF12RegularField,  JaffeMagneticField, HelixMagneticField, UniformMagneticField, YMW16, Sun2008MagneticField
    
try: 
    from _ImagineModels import JF12RandomField, ESRandomField, GaussianScalarField, LogNormalScalarField
    __has_random_fields__ = True
except ImportError:
    __has_random_fields__ = False
    
        
from .HelperFunctions.CoordinateConversions import cyl2cart

from .MagneticFields.RegularMagneticFields import AxiSymmetricSpiral
