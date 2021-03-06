r"""

**Forward modeling of wave equation**

1D, 2D modeling using implicit and explicit algorithms

For example:.

* :mod:`~WaveFd.Wavelet`: energy source functions

"""

from Wavelet import RickerSource, LinearSin, Sinc, Triangle, Ricker
from Imp1DLuWave import Imp1DLuWave
#from Exp1DLaxWave import LaxWand1DWave
from Imp2DLuWave import Imp2DLuWave as _Imp2DLuWave
from Imp2DLuWave import Imp2DLuSparseWave
from WaveAnim import Wave1DAnim, Wave2DAnim, Wave1DPlot, Wave2DShow
from Exp2DFourierWave import Exp2DFourierWave
import Exp2DWave
import Utils
