import Wavelet
import Fourier2D
import Filters

"""
to test the derivatives based on polynomial interpolations...
"""

def main():
    wv = Wavelet.SourceWavelet(20.0)
    wf = Wave2D.WaveField2D(N=256, M=256, Wavelet=wv)
    wf.Sj=50
    wf.Si=50
    wf.SetVe(200.0)
    wf.SetRh(210.0)
    wf.t = 1

def smallscaleTest():
    """
works with the Wave2DinterpZero but with a very strong dispersion
    """
    times = zeros([200, 50, 50])
    wv = Wavelet.SourceWavelet(30.0)
    wf = Wave2D.WaveField2D(N=50, M=50, Wavelet=wv)
    wf.Sj=4
    wf.Si=4
    wf.SetVe(2000.0)
    wf.SetRh(2100.0)
    wf.Ds=20
    wf.Dt=0.005
    wf.t = 1
    wf._GetWavelet()

if __name__ == '__main__':
    main()
