import os

from _ImagineModels import VectorFieldBase, ScalarFieldBase, RegularVectorField, RegularScalarField, \
    JF12RegularField,  JaffeMagneticField, HelixMagneticField, UniformMagneticField, YMW16, Sun2008MagneticField
    
if os.environ['FFTW_FOUND']:
    from _ImagineModels import JF12RandomField, ESRandomField, GaussianScalarField, LogNormalScalarField
        
from .HelperFunctions.CoordinateConversions import cyl2cart

from .MagneticFields.RegularMagneticFields import AxiSymmetricSpiral
