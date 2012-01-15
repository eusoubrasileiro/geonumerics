import pyplot
import numpy
from numpy.random import rand


def Periodic(N=200, dt=0.05, plot=False):
    """
    periodic signal defined by
    y=2*sin(2*pi*x)+3*cos(2*pi*3*x)+sin(2*pi*3*x)+1.5*cos(2*pi*6*x)
    1*1Hz + 3*3Hz+ 1.5*6Hz
    N = number of samples
    dt = sample rate
    plot = show espectrum after
    """
    x0=0 # inicio da amostrag comeca em 0
    x1=dt*N # fim da amostragem
    # k usados para amostrar a funcao
    k=numpy.arange(x0,x1,dt)
    x=k[range(1,N)] # soh N amostras??
    fn = 1/(2*dt)
    print "nyquist = %f" % fn
    y=2*numpy.sin(2*numpy.pi*x)+3*numpy.cos(2*numpy.pi*3*x)
    y = y + numpy.sin(2*numpy.pi*3*x)+1.5*numpy.cos(2*numpy.pi*6*x)
    if(plot==True):
        pyplot.plotfftNabs_phase(y, dt, True);
    return y;


def Periodic_Noise(N=200, dt=0.05, plot=False):
    """
    periodic signal defined by
    y=2*sin(2*pi*x)+0.4*cos(2*pi*3*x)+sin(2*pi*3*x)+0.4*cos(2*pi*5*x)
    2*2Hz + 0.4*3Hz+ 1*3Hz+0.4*5Hz
    plus random noise
    N = number of samples
    dt = sample rate
    plot = show espectrum after
    """
    #uma maneira de montar a escala do espectro de frequencia baseado
    #na frequencia de nyquist
    x0=0 # inicio da amostrag comeca em 0
    x1=dt*N # fim da amostragem
    # x usados para amostrar a funca
    x=numpy.arange(x0,x1,dt) # soh N amostras
    fn=1/(2*dt)
    print "frequencia de nyquest = %f " % fn

    y=2*numpy.sin(2*numpy.pi*x)+0.4*numpy.cos(2*numpy.pi*3*x)+numpy.sin(2*numpy.pi*3*x)+0.4*numpy.cos(2*numpy.pi*5*x)
    y += -2+4*rand(N);

    if(plot==True):
        pyplot.plotfftNabs_phase(y, dt, True);

    return y;


# no dominio continuo, a serie de fourier torna-se transformada de fourier
# no limite
#
# frequencia fundamental fp= 1/T com T = (N*dt)
# no maximo N/2*fp


def NonLinear(N=200, dt=0.005, fc=50, ct=500, plot=False):
    """
    intervalo de amostragem, amostras, fq. central e constante de decaimento exponencial
    y=cos(2*pi*fc*x).*exp(-ct*x.*x)
    na frequencia de nyquist
    N = 200;
    dt = 0.01; 50Hz Nyquist
    fc = 25;
    ct = 500;
    espc bonito por exemplo
    """
    i = numpy.arange((-N*dt/2),N*dt/2,dt) # divide o mesmo numero de amostras pros dois lados
    x = i[range(1,N)]
    fn=1/(2*dt)
    print "frequencia de nyquest = %f" % fn
    y=numpy.cos(2*numpy.pi*fc*x)*numpy.exp(-ct*x*x)
    #multiplicacao no tempo por um coseno(2*pi*fc*x)
    # igual a convolucao com deltas de frequencia fc
    # igual a deslocar o centro do espectro para a frequencia fc
    if(plot==True):
        pyplot.plotfftNabs_phase(y, dt, True);

    return y

