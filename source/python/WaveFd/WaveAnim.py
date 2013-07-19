import numpy as np
import pylab as py
import sys
from os import system
from matplotlib import cm

def _ClearTempImages(name, imfmt):
    try:
        system("rm "+name+"*."+imfmt)
    except:
        pass

def AnimFromPng(name, gif=True, fps=15):
    r"""
    Converts a collection of image files .png format into an avi or gif file.
    Requires ffmpeg .avi.
    Requires ImageMagick .gif.

    * name  :  file name begin of the image files .png
      also the name of the output .avi file
    * gif   :  if False create an .avi file instead of a .gif
    * fps   :  frames per second
    """
    if(gif):
        imgconvert = "convert " + "-delay " + str(int(1000/fps))
        imgconvert += " -dispose None " + name + "*.png -loop 0 " + name + ".gif"
        system(imgconvert)
        print imgconvert
    else:
        aviconvert = "ffmpeg -i " + name + "%03d.png -b:v 2048k -r " + str(fps) + " " + name + ".avi"
        system(aviconvert)
        print aviconvert


def Wave1DPlot(ufield, extent, vel=None):
    r"""
    Draw a 1D pressure function at some instant of time.
    As background is shown velocity field.
    Same dimension as ufield.

    * ufield    : pressure function at an instant of time
    * extent    : (xmin, xmax, ymin, ymax) imshow limits
    * vel       : 1d background velocity function
    """
    # plot curve
    py.hold(True)
    py.plot(ufield,'w-',linewidth=2)
    # create velocity overlay to plot behind
    velbar = np.zeros((2,np.size(ufield)))
    velbar[:] = vel
    # optional cmap=cm.jet, apect='auto' adjust aspect to the previous plot
    py.imshow(velbar, interpolation='bilinear', extent=extent, alpha=0.7, origin='lower', aspect='auto')
    cb = py.colorbar(orientation='horizontal')
    cb.set_label("velocity (m/s)")
    py.xlabel("distance (m)")
    py.hold(False)
    py.show()


def Wave1DAnim(snapshots, ds, dt, vel=None, filename='wave1danim', anim="gif", fps=10):
    r"""
    Create an animation file from a matrix resulting from a simulation of a 1d wave field.
    Creates many intermediate files to achieve that, uses ImageMagick

    * snapshots : is a 2d matrix [time][nx] - pressure field
    * ds        : space equal in x axis
    * dt        : time increment between simulation steps
    * vel       : 1d background velocity field
    * filename  : file name for the animation file
    * anim      : animation type (gif or avi)
    * fps       : frames per second
    """

    py.ion()
    #max index time and max index space
    maxt = np.shape(snapshots)[0]
    maxk = np.shape(snapshots)[1]
    # get the maximum and minimum u value to not blow the scale
    # during the movie
    ymax = ymin = snapshots[0][0]
    for snapshot in snapshots:
        ymaxi = max(snapshot)
        ymini = min(snapshot)
        if( ymaxi > ymax ):
            ymax = ymaxi
        if( ymini < ymin ):
            ymin = ymini
    # extents of the picture x starts at 0
    xmin, xmax = 0, maxk*ds
    extent= xmin, xmax, ymin, ymax
    # font position
    width = xmax-xmin
    height = ymax-ymin
    posx = 0.8*width+xmin
    posy = 0.8*height+ymin
    # not working?
    # verticalalignment='top',
    # horizontalalignment='right',
    _ClearTempImages(filename, "png")   # clear any previous existing
    for t in range(maxt):
        Wave1DPlot(snapshots[t], extent, vel)
        # set axis ranges
        py.hold(True)
        # draw time
        py.text(posx, posy, "{0:1.5f}".format(t*dt), alpha=0.8, color='b')
        # since its just math objects would be perfect
        # will be something like Wave1DAnim001.png
        py.savefig(filename+"{0:03d}".format(t)+'.png', dpi=150)
        sys.stdout.write("\r progressing .. %.1f%%" %(100.0*float(t)/maxt))
        sys.stdout.flush()
        py.clf()
    sys.stdout.write(" done! \n")
    py.ioff()
    py.hold(False)
    py.close()
    if ( anim == "gif"):
        AnimFromPng(filename, fps=fps)
    else :
        AnimFromPng(filename, False, fps)
    _ClearTempImages(filename, "png")
    #clear after creating new files


def Wave2DShow(ufield, ds, vel=None, vmin=None, vmax=None):
    r"""
    Show a 2D pressure field at some instant of time.
    As background is shown velocity field.
    Same dimension as ufield.

    * ufield    : 2d pressure field at an instant of time
    * ds        : space discretization
    * vel       : 2d background velocity field
    * vmin/vmax : vmin/vmax of imshow
    """
    #max index time and max index space
    maxt = np.shape(snapshots)[0]
    maxk = np.shape(snapshots)[1]
    maxi = np.shape(snapshots)[2]    
    if vmin == None or vmax == None :
        # get the maximum and minimum values of the last 5% 
        # snapshots to not blow the scale during the animation
        snaptmp = snapshots[-int(0.05*maxt):]
        vmax = vmin = snaptmp[0][0][0]
        for snapshot in snaptmp:
            for line in snapshot:
                linemax = max(line)
                linemin = min(line)
                if(linemax > vmax):
                    vmax = linemax
                if(linemin < vmin):
                    vmin = linemin

    print "vmin : ", vmin, "vmax : ", vmax
    # space axis starting at 0 in x and z (using y coz' plotting)
    # extents of the picture,
    xmin, xmax = 0, ds*maxi
    ymin, ymax = 0, ds*maxk
    extent= xmin, xmax, ymax, ymin
    py.hold(True)
    if not vel == None:
        py.imshow(vel, interpolation='bilinear', cmap=cm.jet, extent=extent,  origin='upper', aspect='auto')

    py.imshow(ufield, interpolation='bilinear', cmap=cm.Greys_r, alpha=0.8, extent=extent, origin='upper', aspect='auto', vmin=vmin, vmax=vmax)
    py.hold(False)
    # optional cmap=cm.jet, apect='auto' adjust aspect to the previous plot
    py.show()
    #py.xlabel("distance (m)")
    #cb = py.colorbar(orientation='horizontal')
    #cb.set_label("velocity (m/s)")


def Wave2DAnim(snapshots, ds, dt, vel, filename='wave2danim', norm=True, vmin=None, vmax=None, anim="avi", fps=15):
    r"""
    Create an animation file from a matrix resulting from a simulation of a 2d wave field.
    Creates many intermediate files to achieve that, uses ImageMagick.
    Z is downward.

    * snapshots : is a 3d matrix [time][nz][nx] - pressure field
    * ds        : space equal in x/y axis
    * dt        : time increment between simulation steps
    * vel       : 2d background velocity field
    * filename  : file name for the animation file
    * anim      : animation type (gif or avi)
    * fps       : frames per second
    * norm      : scale the values getting the general max and min (vmax/vmin)
    * vmin      : global minimum of snapshots
    * vmax      : global maximum of snapshots
    """
    py.ion()
    #max index time and max index space
    maxt = np.shape(snapshots)[0]
    maxk = np.shape(snapshots)[1]
    maxi = np.shape(snapshots)[2]    
    if norm :
        # get the maximum and minimum values of the last 5% 
        # snapshots to not blow the scale during the animation
        snaptmp = snapshots[-int(0.05*maxt):]
        vmax = vmin = snaptmp[0][0][0]
        for snapshot in snaptmp:
            for line in snapshot:
                linemax = max(line)
                linemin = min(line)
                if(linemax > vmax):
                    vmax = linemax
                if(linemin < vmin):
                    vmin = linemin

        print "vmin : ", vmin, "vmax : ", vmax
    # space axis starting at 0 in x and z (using y coz' plotting)
    # extents of the picture,
    xmin, xmax = 0, ds*maxi
    ymin, ymax = 0, ds*maxk
    extent= xmin, xmax, ymax, ymin
    # font position
    width = xmax-xmin
    height = ymax-ymin
    posx = 0.8*width+xmin
    posz = 0.8*height+ymin
    # not working?
    # verticalalignment='top',
    # horizontalalignment='right'
    _ClearTempImages(filename, "png") # clear any previous existing
    for t in range(maxt):
        py.hold(True)
        py.imshow(vel, interpolation='bilinear', cmap=cm.jet, extent=extent,  origin='upper', aspect='auto')
        py.imshow(snapshots[t], interpolation='bilinear', cmap=cm.Greys_r, alpha=0.8, extent=extent, origin='upper', aspect='auto', vmin=vmin, vmax=vmax)
        # optional cmap=cm.jet, apect='auto' adjust aspect to the previous plot
        py.show()            
        # draw time
        py.text(posx, posz,
                "{0:1.5f}".format(t*dt),
                alpha=0.8,
                style='italic',
                color='b')
        # since its just math objects would be perfect
        # will be something like Wave1DAnim001.png
        py.savefig(filename+"{0:03d}".format(t)+'.png', dpi=150)
        sys.stderr.write("\r progressing .. %.1f%%" %(100.0*float(t)/maxt))
        sys.stderr.flush()
        py.clf()
    sys.stdout.write(" done! \n")
    py.hold(False)
    py.ioff()
    py.close()
    if ( anim == "gif"):
        AnimFromPng(filename, fps=fps)
    else :
        AnimFromPng(filename, False, fps)
    _ClearTempImages(filename, "png")
    #clear after creating new files


