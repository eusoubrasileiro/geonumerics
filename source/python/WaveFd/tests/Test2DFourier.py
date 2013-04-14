import Wavelet
import Fourier2D


def main():
    wv = Wavelet.SourceWavelet(20.0)
    wf = Fourier2D.FourierWaveField2D(N=256, M=256, Wavelet=wv)
    wf.SetVe(1000.0)
    wf.SetRh(2.0)
    wf._GetDeltaTime()
    wf._GetWavelet()
    wf.Wavelet = Wavelet.SincWavelet(513, 20.0, wf.Dt)
    wf.TotalEstimatedTime()    

if __name__ == '__main__':
    main()
