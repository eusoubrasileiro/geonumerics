#!/usr/python

import numpy


# filtro caixa passa baixa (frequencia corte fc)
# infinite impulse response? e filtro nao causal
# 1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
# logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
# mesmo que convolver com o sinc dessa caixa na frequencia
# que acarreta gibs no modulo e fase?
# 2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
# amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia
# com deltas espacados de 1/dt e amplitude 1/dt
def SincLowPass(n, fc, dt):
    """
    n is number of samples for the filter
    fc is the F cut-off
    dt is the sample rate
    Sample the Box Sinc filter simetrically around the zero
    """

    if(1/(2*dt) < fc):
        print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
        return

    # Amostra simetricamento o operador do filtro em torno do zero
    # nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
    x = y = numpy.arange(0, 1, 1)

    # caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
    #um workaround eh utilizado para evitar divsao por 0 e utilizasse o limite em 0 para setar o valor no 0
    if(n%2!=0):
        #Nota: arange omite o ultimo valor, por isso o + dt, para o intervalo ficar [  , ]  e nao [ , )
        x = numpy.arange(-dt*(n-1)/2,(dt*(n-1)/2)+dt, dt)
        # seta o zero como Nan para evitar, excessao de divisao por 0
        x[numpy.size(x)/2]=numpy.NaN
        y = dt*numpy.sin(2*numpy.pi*fc*x)/(numpy.pi*x) # sinc da frequencia de corte, o termo
        # dt multiplicando serve para garantir a o espc. amplitude em 1, devido a amostragem do filtro no tempo! convolucao com deltas!!
        # inverso da caixa na frequencia
        # dt multiplicando serve para garantir a o espc. amplitude em 1
        # set o valor de y no 0, baseado no limite x->0 do operador
        # fazendo o limite (derivando em baixo e em cima) chega-se para Fc => cos(2*pi*Fc*t)*2*Fc
        # com t =0 => 2*Fc, do not forget also the dt to normalize to 1
        y[numpy.size(y)/2]= 2*fc*dt
    else:
        # nao faz muito sentido no mundo real e ninguem usa um filtro nao simetrico mas vai ai
        #amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero
        #Nota: arange omite o ultimo valor, por isso o + dt para o intervalo ficar [  , ]  e nao [ , )
        dispX = dt*(n.__float__()-1)/2
        x = numpy.arange(-dispX,dispX+dt, dt)
        dt*numpy.sin(2*numpy.pi*fc*x)/(numpy.pi*x); # sinc da frequencia de corte, o termo
        # # dt multiplicando serve para garantir a o espc. amplitude em 1, devido a amostragem do filtro no tempo! convolucao com deltas!!

    print " Filter number of samples %d" % y.__len__()

    return y



def SincBandPass(n=200, dt=0.01, f1=10, f2=20):
    """
    filtro caixa passa banda (frequencias inicias e finais f1 e f2)
    infinite impulse response? e filtro nao causal
    1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
    logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
    mesmo que convolver com o sinc dessa caixa na frequencia
    que acarreta gibs no modulo e fase?
    2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
    amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia
    com deltas espacados de 1/dt e amplitude 1/dt
    To create a band pass filter based on a box, we have to multiply the sinc in time
    for the frequency we want to be the center of our band pass filter.
    Multiplying in time displaces our sinc in frequency, convolving with one delta of the
    that frequency in the frequency domain
    plot = to show or not the response for the filter
    """

    if(f2 < f1):
        print "Dont be stupid f2 must be bigger than f1"
        return

    if(1/(2*dt) < max(f2,f1)):
        print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
        return

    Fdelta = (f1+f2)/2 # frequencia central
    Fbox = (f2-f1)/2 # box translated to the central frequency

    # Amostra simetricamento o operador do filtro em torno do zero
    # nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
    x = y = numpy.arange(0, 1, 1)

    # caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
    #um workaround eh utilizado para evitar divsao por 0 e utilizasse o limite em 0 para setar o valor no 0
    if(n%2!=0):
        #Nota: arange omite o ultimo valor, por isso o + dt, para o intervalo ficar [  , ]  e nao [ , )
        x = numpy.arange(-dt*(n-1)/2,(dt*(n-1)/2)+dt, dt)
        # seta o zero como Nan para evitar, excessao de divisao por 0
        x[numpy.size(x)/2]=numpy.NaN
        y = 2*dt*numpy.sin(2*numpy.pi*Fbox*x)*numpy.cos(2*numpy.pi*Fdelta*x)/(numpy.pi*x); # inverso da caixa (soh parte real), o termo
        # dt multiplicando serve para garantir a o espc. amplitude em 1
        # nao sei o pq do 2? matematicamente
        # set o valor de y no 0, baseado no limite x->0 do operador
        # fazendo o limite (derivando em baixo e em cima) chega-se para a => cos(2*pi*Fbox*t)*2*Fbox
        # com t =0 => 2*Fbox
        # limite da multplicacao eh a multiplicacao dos limites
        # limite soh do sinc utilizando * 2 anterior
        y[numpy.size(y)/2]= 2*dt*2*Fbox#2*numpy.pi*Fdelta

    else:
        #amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero
        #Nota: arange omite o ultimo valor, por isso o + dt para o intervalo ficar [  , ]  e nao [ , )
        dispX = dt*(n.__float__()-1)/2
        x = numpy.arange(-dispX,dispX+dt, dt)
        y = 2*dt*numpy.sin(2*numpy.pi*Fbox*x)*numpy.cos(2*numpy.pi*Fdelta*x)/(numpy.pi*x); # inverso da caixa (soh parte real), o termo

    print " Filter number of samples %f" % y.__len__()
    print " Filter Nyquest frequency %f" % (1/(2*dt))
    print " Filter F1: %f F2: %f" % (f1, f2)
    
    return y


def SincTrapezoidalLowPass(N, dt, Ramp, Fc):
    """
    filtro trapezoidal passa baixa ( rampa = Ramp, Frequencia de corte = Fc )
    infinite impulse response cortada ficando finita (fir) e filtro nao causal
    1) filtro eh truncado no tempo pois nao podemos ter infinitas amostras...
    logo o filtro eh multiplicado por uma funcao caixa no tempo que eh o
    mesmo que convolver com o sinc dessa caixa na frequencia
    que acarreta gibs no modulo e fase
    2) o filtro eh amostrado no tempo, ou seja, multiplicado por deltas de
    amplitude 1 e espacados de dt. Que eh o mesmo que convolver na frequencia
    com deltas espacados de 1/dt e amplitude 1/dt
    dt = dt  Deve ser  tal que  === 1/2*dt > Fc
    """
    a = Ramp/2
    b = Fc - (Ramp/2)

    # when sampling the filter operator simetrical to 0
    # the number of samples of the filter is :  N*2+1 where N is the order of the filter
    # should I change here? the implementation

    print " Filter Nyquest frequency %f" % (1/(2*dt))
    print " Filter cut-off frequency Fc: %f Ramp: %f" % (Fc, Ramp)

    if(Ramp > Fc):
        print "Dont be stupid Ramp bigger than Fc?"
        return

    if(1/(2*dt) < Fc):
        print "Dont be stupid you have to sample the filter at least at the nyquest frequency"
        return

    # this doesnt work, ramp = 0
    if(a==0):
        print "Doesnt make sense a = 0"
        return

    # Amostra simetricamento o operador do filtro em torno do zero
    # nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
    x = y = numpy.arange(0, 1, 1)

    # caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
    #um workaround eh utilizado para evitar divsao por 0 e utilizasse o limite em 0 para setar o valor no 0
    if(N%2!=0):
        #Nota: arange omite o ultimo valor, por isso o + dt, para o intervalo ficar [  , ]  e nao [ , )
        x = numpy.arange(-dt*(N-1)/2,(dt*(N-1)/2)+dt, dt)
        # seta o zero como Nan para evitar, excessao de divisao por 0
        x[numpy.size(x)/2]=numpy.NaN
        y = dt*numpy.sin(2*numpy.pi*a*x)*numpy.sin(2*numpy.pi*b*x)/(numpy.pi*numpy.pi*x*x*2*a);
        # inverso da convolucao de
        # duas caixas na frequencia
        # dt multiplicando serve para garantir a o espc. amplitude em 1
        # a divisao por 2*a, serve para garantir o spec. amplitude em 1,
        # pois o resultado da convolucao de duas caixas de amplitude 1 na frequencia (-a, a), (-b, b) com b > a
        #resulta no trapezio com amplitude maxima igual 2*a
        # set o valor de y no 0, baseado no limite x->0 do operador
        # fazendo o limite (derivando em baixo e em cima) chega-se para a => cos(2*pi*a*t)*2*a
        # com t =0 => 2*a, ou seja para a e b => 2*a*2*b
        # limite da multplicacao eh a multiplicacao dos limites
        y[numpy.size(y)/2]= 2*a*2*b*dt/(2*a)
    else:
        #amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero
        #Nota: arange omite o ultimo valor, por isso o + dt para o intervalo ficar [  , ]  e nao [ , )
        dispX = dt*(N.__float__()-1)/2
        x = numpy.arange(-dispX,dispX+dt, dt)
        y = dt*numpy.sin(2*numpy.pi*a*x)*numpy.sin(2*numpy.pi*b*x)/(numpy.pi*numpy.pi*x*x*2*a);

    print " Filter number of samples %f" % y.__len__()


    return y
