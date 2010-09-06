using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace LogFilterModule
{
    public static class Filtering
    {

        public static double[] LowPass(double[] Signal, double Dt, double Fc, double Ramp)
        {   
            // filter operator in time
            double[] Filter_kernel;

            // Size must be always ODD!! ALWAYS!! 2n+1 where n is the order of the filter            
            // Previous problem... if the sample rate is for example 1E-5 using 
            // a filter with 101 samples, cant represent the spectrum properly            
            // and transition bandwith, the length of the transition part,
            // for the two sinc's the transition btw can be aproximated by
            // BTW = (4/Nf)*(1/(2*dt)) -> e.g dt=2E-5 Nf=101 715 Hz uhuuu considering
            // that I want to filter 125 hz...ahhaa Horrible!!
            // relative band width transition achived by dividing by the cutoff frequency            
            // 20% has an accepted ripple/transition
            int Nf = FilterSize(0.20, Fc, Dt);
            Nf = (Nf < 101) ? 101 : Nf; // minimum size 101            
            Filter_kernel = TrapezoidalLowPass(Nf, Dt, Fc, Ramp);           

            // the best option found to contourn the border effects in the convolution
            return FilterConvolutionCopying(Signal, Filter_kernel);
        }

        /// <summary>
        ///"""
        ///Gives the number of samples (odd number) needed to sample our filter operator
        ///based on:
        ///1) the relative transition bandwidth RTbtw WE WANT 
        ///(Transition band size divide by cut-off frequency)
        ///2) Filter cut-off frequency 
        ///3) Sample rate that's gonna be used to sample the filter
        ///This aproximation is valid for Sinc Filters.
        ///a RTbtw = 20% is a reasonable value for a non distorded frequency response	
        ///"""
        /// </summary>
        /// <param name="RTbtw">Relative bandwidth transition desired percentage</param>
        /// <param name="Fc">cut-off frequency</param>
        /// <param name="dt">sample rate which the filter will be sampled</param>
        /// <returns>the number of samples needed</returns>
        public static int FilterSize(double RTbtw, double Fc, double dt){            
	        int NFs =(int) (2/(RTbtw*Fc*dt));
            NFs = (NFs % 2 == 0) ? NFs + 1 : NFs; // get the next ODD
	        return NFs;
        }
      
        /// <summary>
        /// Effectively makes the convolution and clips the result to
        /// the input signal size
        /// </summary>
        /// <param name="Signal"></param>
        /// <param name="Filter"></param>
        /// <returns></returns>
        public static double[] FilterConvolutionCopying(double[] Signal, double[] Filter)
        {
            double[] InputSignal = (double[])Signal.Clone(); // just to avoid modifications in the original one

            // unaceptable number of samples for a filter
            if (Filter.Length % 2 == 0)
                return null;

            //para garantir tamanho final identico ao inicial
            // sem modificacoes na fase, todo sinal de entrada deve ser impar	
            // to avoid the border effect due filtering
            // after trying padding with zeros the input signal/ the filter
            // after put the filter in the wrap-around form (numerical recipies tip (convolution))
            // after padding with random values between the min and max values of the input signal
            // after adding one period ... at the beginning and at the end add Period/2 parts
            // the actual solution
            // add samples add the end and beginning with copying the last and first sample
            // those are signal extensions at the beginning and end,
            // they contain begin part (first sample copies) and end part (end sample copies)
            double[] SBegin;
            double[] SEnd;
            
            
            if (InputSignal.Length % 2 == 1)
            {
                // In this case we allways add an EVEN number of samples
                // at the begin and at the end
                // Because EVEN+(signal)ODD+EVEN = ODD what we want
                //NS is ODD, so the begin part is equal the end part
                int NS = InputSignal.Length;
                // e.g. size = (7-1)/2 +1 = 4 ou (5-1)2 + 1 = 3
                // so it's need to make one if condition to get one even number of samples                
                int NSAdd = NS / 2;
                NSAdd = (NSAdd % 2 == 1) ? NSAdd + 1 : NSAdd;

                SBegin = new double[NSAdd];
                for (int i = 0; i < SBegin.Length; i++)
                    SBegin[i] = InputSignal[0];
                
                SEnd = new double[NSAdd];
                for (int i = 0; i < SEnd.Length; i++)
                    SEnd[i] = InputSignal[NS-1];                
            }
            else
            {// even case, 
                //NS is EVEN, so the begin part is smaller than the end part, because 
                // of one more sample added at the end
                // Because EVEN+(signal)EVEN+ODD = ODD what we want    
                // added to make the signal odd
                int NS = InputSignal.Length;
                // e.g. size = 6/2 = 3 or 8/2 = 4                
                // so it's need to make one if condition to get one even number of samples                
                int NSAdd = NS / 2;
                NSAdd = (NSAdd % 2 == 1) ? NSAdd + 1 : NSAdd;

                SBegin = new double[NSAdd];
                for (int i = 0; i < SBegin.Length; i++)
                    SBegin[i] = InputSignal[0];

                SEnd = new double[NSAdd];
                for (int i = 0; i < SEnd.Length; i++)
                    SEnd[i] = InputSignal[NS - 1]; 

                // we also add one more sample due the input signal not being odd, 
                // in the end part that will be putted at the end... 
                SEnd = SEnd.Concat(new double[] { InputSignal[NS - 1] }).ToArray();
            }


            // add a begin part with a lot of copies of the first sample      
            // and a end part with a lot of copies of the last sample
            InputSignal = SBegin.Concat<double>(InputSignal.Concat<double>(SEnd)).ToArray();

            // Convolve the Filter with the signal
            // both have to be odd so N+M-1 is odd also
            // the convolution creates negative values for the ILD log, when the input
            // log doesnt have any value in this range, it's a log values scaled data
            // the convolution algoritmo it's wrong problably wrong
            double[] FilterResult = Convolve(Filter, InputSignal);

            // After filtering get what is interesting for us
            // the central part of the convolution

            //RESULTADO DA FILTRAGEM CUTED, 
            //remove o numero total de amostras do filtro
            //deixa o sinal do tamanho do sinal inicial
            //Ns  numero de amostras do sinal inicial
            //Nf numero de amostras do filtro	

            //simetrical samples around Nrs/2 for Nrs odd, all Integer operations
            int Nrs = FilterResult.Length;
            int Ns = InputSignal.Length;

            // samples the we are interested
            int beg = (Nrs / 2) - (Ns - 1) / 2;
            int size = Ns;

            // corta, removendo o numero equivalente as amostras do filtro somente
            double[] Result = new double[size];
            Array.Copy(FilterResult, beg, Result, 0, size);

            // Remove also the undesired added period, begin part and final part            
            double[] Final = new double[Signal.Length];
            Array.ConstrainedCopy(Result, SBegin.Length, Final, 0, Final.Length);

            return Final;
        }
        

        public static double[] TrapezoidalLowPass(int N, double Dt, double Fc, double Ramp)
        {
            // this filter is the result of convolution of two box in frequency
            // the box A with its side/2 = a
            double a = Ramp/2;
            // th box B with its side/2 = b
            double b = Fc - (Ramp/2);
	
	        // this doesnt work, ramp = 0, error 
	        if(a==0)
		        throw new Exception("Impossible size of ramp");
	
            double[] x = null; // imput x/t values to sample the filter operator 
            double[] filter = new double[N]; // filter kernel sampled

	        // Amostra simetricamento o operador do filtro em torno do zero
	        // nao importa se o numero de pontos for par ou impar, exitem opcoes pros dois casos
    	    
        	
	        // caso seja impar amostra incluindo o 0,  amostra perfeitamente simetrico em torno do zero
	        // um workaround eh utilizado para evitar divsao por 0 e utiliza-se o limite em 0 para setar o valor no 0
            if (N % 2 != 0)
            {
                x = Range(-Dt * (N - 1) / 2, (Dt * (N - 1) / 2), Dt);
                // seta o zero como Nan para evitar, excessao de divisao por 0		
                x[x.Length / 2] = Double.NaN;
                // sample the function
                for (int i = 0; i < x.Length; i++)
                    filter[i] = Dt * Math.Sin(2 * Math.PI * a * x[i]) * Math.Sin(2 * Math.PI * b * x[i]) / (Math.PI * Math.PI * x[i] * x[i] * 2 * a);

                // inverso da convolucao de
                // duas caixas na frequencia
                // Dt multiplicando serve para garantir a o espc. amplitude em 1	
                // a divisao por 2*a, serve para garantir o spec. amplitude em 1, 
                // pois o resultado da convolucao de duas caixas de amplitude 1 na frequencia [-a, a], [-b, b] com b > a 
                // resulta no trapezio com amplitude maxima igual 2*a		
                // set o valor de y no 0, baseado no limite x->0 do operador
                // fazendo o limite (derivando em baixo e em cima) chega-se para a => cos(2*pi*a*t)*2*a
                // com t =0 => 2*a, ou seja para a e b => 2*a*2*b
                // limite da multplicacao eh a multiplicacao dos limites
                filter[filter.Length / 2] = 2 * a * 2 * b * Dt / (2 * a);
            }
            else
            {
                //amostra  simetricamente em torno do zero sem passar pelo 0, mais facil, nao precisa evitar divisao por zero		        
                double Xdisp = Dt * (double)(N - 1) / 2;
                x = Range(-Xdisp, Xdisp, Dt);
                // sample the function
                for (int i = 0; i < x.Length; i++)
                    filter[i] = Dt * Math.Sin(2 * Math.PI * a * x[i]) * Math.Sin(2 * Math.PI * b * x[i]) / (Math.PI * Math.PI * x[i] * x[i] * 2 * a);
            }

            // uhuu the filter done
            return filter;
        }

        private static double[] Range(double Begin, double End, double Step)
        {
            int N = (int) ((double)(End-Begin)/Step);
            // one sample more because its a closed interval
 	        double[] range = new double[++N];

            for(int i=0; i<N; i++)                
                range[i] = Begin+i*Step;

            return range;
        }


}
