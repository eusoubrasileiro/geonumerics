#!/usr/bin

import numpy as np
import pylab as py
from os import system


def _Clear(name, imfmt):
    system("rm "+name+"*."+imfmt)
    
    return

def _GifFromPng(name, delay=100):
    """
    requires imagemgick installed for gif conversion
    something like convert example.gif -delay 100 example*.png -loop 0
    """ 
    imgconvert = "convert " + "-delay " + str(delay) 
    imgconvert += " -dispose None " + name + "*.png -loop 0 " + name + ".gif"
    
    system(imgconvert)
    print imgconvert
    
    return


def _AviFromPng(name):
    """
    requires ffmpeg for avi conversion
    something like 
    aviconvert = "ffmpeg -i +" name + "*.png " + name ".avi" 
    """
    aviconvert = "ffmpeg -i " + name + "%03d.png " + name + ".avi" 
    
    system(aviconvert)
    print aviconvert
    
    return


def Wave1DAnim(snapshots, name='Wave1DAnim', bckg=None, anim="gif", Ds=1, Dt=0.01, frameint=10):
    """
    snapshot is a 2d matrix [time][frames]
    anim = animation tipe (gif or avi)
    name is the file name of the animation image
    bckg is something else to plot as background (like velocity)  
    Ds space in x axis (1 meter)
    Dt time increment between simulation steps (0.01 second)
    frameint is in miliseconds  = the interval between frames
    """
    #clear before creating new files
    _Clear(name,"png")
    _Clear(name, anim)
    
    py.ion()
    
    #max index time and max index space
    tm = np.shape(snapshots)[0]
    sm = np.shape(snapshots)[1]
    # time axis and space axis
    timeaxis = np.arange(0,tm*Dt, Dt)
    spaceaxis = np.arange(0, sm*Ds, Ds)
    
    # get the maximum and minimum y value to not blow the scale
    mx = mn = 0.0
    for i in range(tm):
        mxi = max(snapshots[i])
        mni = min(snapshots[i]) 
        if( mxi > mx ):
            mx = mxi
        if( mni < mn ):
            mn = mni
    
    for t in range(tm):
        py.plot(spaceaxis, snapshots[t])
        
        if bckg != None: # plot background
            py.plot(spaceaxis, bckg)
        
        py.axis(ymax=mx, ymin=mn) # set axis ranges
        py.text(0.7*sm*Ds, 0.7*mx, "{0:1.5f}".format(timeaxis[t])) # draw time
        # since its just math objects would be perfect
        # will be something like Wave1DAnim001.png
        py.savefig(name+"{0:03d}".format(t)+'.png', dpi=50)
        py.clf()
    
    py.close()
    
    if ( anim == "gif"):
        _GifFromPng(name, frameint)
        _Clear(name, "png")
        
    if ( anim == "avi"):
        _AviFromPng(name)
    
    #clear after creating new files
    #_Clear(name,"png")
