# This Python file uses the following encoding: utf-8
"""
Seismic: 2D finite difference simulation of elastic P and SV wave propagation
Simulates Parana Basin geometry of basalt cover with horizontal layers.
Acquires just z displacement (default seismic acquisition).

4 layers :
bauru ~ 100 meters ~ 6x scaled to 600 m
serra geral flow 1 ~ 20 meters ~ 6x scaled to 120 m
serra geral flow 2 ~ 20 meters ~ 6x scaled to 120 m
botucatu ~ till end

Example based on pg 281/283 - Análise do Sinal Sísmico Rosa, André L.
plus personal insight of vibroseis data of Parana Basin acquisition 0314 . scaled x 6

rock properties

bauru : vp = 2181 vs = 1166 rho = 2.12
serra geral flow 1 : vp = 4800 vs = 2611 rho = 2.75
serra geral flow 2 : vp = 3540 vs = 1906 rho = 2.54 (vulgular)
botucatu : vp = 3200 vs = 1600 rho = 2.65 (~botucatu sandstone)

"""
import numpy as np
from matplotlib import animation
from fatiando.seismic import wavefd
from fatiando.vis import mpl
from fatiando import utils

# Set the parameters of the finite difference grid
shape = (1, 400, 400)  # z, y, x shape
area = [0, 8000, 0, 2000, 0, 1]  # x min and xmax etc...
dx = area[1]/shape[2]
dy = area[3]/shape[1]
dz = area[5]/shape[0]
# due the third dimension the simulation is held in the x, y plane (x horizontal and y vertical)

# Geologic model make a density and S wave velocity model
density = 2650*np.ones(shape)  # botucatu
density[:, 0:int(600/dy), :] = 2120 # bauru
density[:, int(600/dy):int(600/dy)+int(120/dy), :] = 2750  # basalt flow 1
density[:, int(600/dy)+int(120/dy):int(600/dy)+int(240/dy), :] = 2540  # basalt flow 2
pvel = 3200*np.ones(shape)  # botucatu
pvel[:, 0:int(600/dy), :] = 2181  # bauru
pvel[:, int(600/dy):int(600/dy)+int(120/dy), :] = 4800  # basalt flow 1
pvel[:, int(600/dy)+int(120/dy):int(600/dy)+int(240/dy), :] = 3540  # basalt flow 2

# Make a wave source from a gauss derivative that vibrates in the z direction only
sources = [wavefd.GaussSource((2000, 0, 0), area, shape, 100.0, 5.)]  # x source or shear source

# Get the iterator for the simulation
dt = wavefd.maxdt(area, shape, np.max(pvel))
duration = 3.0
maxit = int(duration/dt)
stations = [[2200+i*dx, 0, 0] for i in range(220)]  # x, z coordinate of the seismometers
snapshot = 1 # int(0.004/dt)  # Plot a snapshot of the simulation every 4 miliseconds
print "dt for simulation is ", dt
print "max iteration for simulation is ", maxit
print "duration for simulation is ", duration
# TODO implement 2D esg version with abs condition
simulation = wavefd.scalar3_esg(pvel, density, area, dt, maxit, sources,
        stations, snapshot, padding=50, taper=0.01)

# This part makes an animation using matplotlibs animation API
fig = mpl.figure(figsize=(12, 5))
# Start with everything zero and grab the plot so that it can be updated later
mpl.imshow(pvel[0][::-1], extent=area[:4], alpha=0.25)
wavefield = mpl.imshow(np.zeros(shape[1:]),  extent=area[:4], vmin=-10**-7, vmax=10**-7, alpha=0.75, cmap=mpl.cm.gray_r)
mpl.points([ stations[i][:2] for i in range(len(stations))], '.k')
mpl.ylim(area[2:4][::-1])
mpl.xlabel('x (km)')
mpl.ylabel('z (km)')
mpl.m2km()
times = np.linspace(0, maxit*dt, maxit)
# the only way for getting the feedback response of the seimograms
global seismograms

# This function updates the plot every few timesteps
def animate(i):
    # t, pressure, seismograms
    t, u, zcomp = simulation.next()
    mpl.title('time: %0.3f s' % (times[t]))
    # wavefield.set_array((p + s)[::-1])
    wavefield.set_data(u[0][::-1])
    global seismograms
    seismograms = zcomp
    return wavefield, seismograms
anim = animation.FuncAnimation(fig, animate, frames=maxit/snapshot, interval=1)
mpl.show()
# turn the seismogram data in a numpy array
# traces = np.array(seismograms)
# parana_shot = utils.matrix2stream(traces)
# mpl.seismic_wiggle(parana_shot, ranges=[200, 130*dx], scale=2, color='b', normalize=True)
# mpl.xlabel('offset (m)')
# mpl.ylabel('time (ms)')
# mpl.show()