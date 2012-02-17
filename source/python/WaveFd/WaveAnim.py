#!/usr/bin

import numpy as np
import pylab as py
from os import system
from matplotlib import cm
from pylab import Axes


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
    aviconvert = "ffmpeg -i " + name+"%03d.png " + name + ".avi" 
    
    system(aviconvert)
    print aviconvert
    
    return
    
    
def Wave1DField(ufield, vel, ds=None, spaceaxis=None, extent=None):
    """
    Draw a 1D wavefield at some instant of time
    As background is used the vel field (same dimension as ufield)
    ufield - displacement field at an instant of time
    vel - velocity field
    ds space increment (optional if spaceaxis is set) 
    spaceaxis (distance for each index of ufield)
    extent (xmin, xmax, ymin, ymax) - to use to limits the plots
    """
    
    if ds == None and spaceaxis == None:
        print "nothing to do!"
        return
    
    # max index space
    sm = np.size(ufield)
    
    if spaceaxis == None:
        # time axis and space axis
        spaceaxis = np.arange(0, sm*ds, ds)
        
    if extent == None:
        # put the limits and position of the image to be plotted
        xmin, xmax, ymin, ymax = min(spaceaxis), max(spaceaxis), min(ufield), max(ufield)
        extent= xmin, xmax, ymin, ymax
    
    # plot initial curve
    py.plot(spaceaxis, ufield,'w-',linewidth=2)
    py.xlabel("distance (m)")
    
    # create velocity overlay to plot behind
    velbar = np.zeros((2,sm))
    velbar[:] = vel
    
    # optional cmap=cm.jet, apect='auto' adjust aspect to the previous plot
    py.imshow(velbar, interpolation='bilinear', extent=extent, alpha=0.7, origin='lower', aspect='auto')
    
    cb = py.colorbar(orientation='horizontal')
    cb.set_label("velocity (m/s)")
    
    # just to certify
    py.axis(extent)
    
    py.show()
    
    
    

def Wave1DAnim(snapshots, name='Wave1DAnim', vel=None, anim="gif", ds=1, dt=0.01, frameint=10):
    """
    snapshot is a 2d matrix [time][frames]
    anim = animation tipe (gif or avi)
    name is the file name of the animation image
    vel is the background velocity  
    ds space in x axis (1 meter)
    dt time increment between simulation steps (0.01 second)
    frameint is in miliseconds  = the interval between frames
    """
    #clear before creating new files
    _Clear(name,"png")
    
    py.ion()
    
    #max index time and max index space
    tm = np.shape(snapshots)[0]
    sm = np.shape(snapshots)[1]
    # time axis and space axis
    timeaxis = np.arange(0,tm*dt, dt)
    spaceaxis = np.arange(0, sm*ds, ds)
    
    # get the maximum and minimum y value to not blow the scale
    # during the movie
    ymax = ymin = 0.0
    for i in range(tm):
        ymaxi = max(snapshots[i])
        ymini = min(snapshots[i]) 
        if( ymaxi > ymax ):
            ymax = ymaxi
        if( ymini < ymin ):
            ymin = ymini
    
    # extents of the picture
    xmin, xmax = min(spaceaxis), max(spaceaxis) 
    extent= xmin, xmax, ymin, ymax
    
    # font position
    width = xmax-xmin
    height = ymax-ymin
    posx = 0.8*width+xmin
    posy = 0.8*height+ymin
    # not working?
    # verticalalignment='top', 
    # horizontalalignment='right',
    
    for t in range(tm):
        Wave1DField(snapshots[t], vel, spaceaxis=spaceaxis, extent=extent)
        # set axis ranges
        py.axis(extent) 
        # draw time
        py.text(posx, posy, 
                "{0:1.5f}".format(timeaxis[t]),
                backgroundcolor='w',
                style='italic',
                color='b') 
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
        _Clear(name, "png")
        
    #clear after creating new files
    #_Clear(name,"png")
    