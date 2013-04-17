#!/usr/python

sys.path.append('..\\source\\python');

import ImplicitAcousticWave as wv

def Demo():
    """
    works but due imshow palete color
    oscilation is not very clear
    """
    field = wv.WaveField()
    
    field.Loop()
    

if __name__ == '__main__':
    Demo()
    





