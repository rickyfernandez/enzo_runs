import yt
import h5py
import numpy as np

from yt.utilities.physical_constants import mh, kboltz

def pressure_k(field, data):
    return data['gas', 'pressure']/kboltz

pc_to_cm = 3.086E18  # parsec
G  = 6.67259E-8      # gravity
#k  = 1.380658E-16    # boltzmann
#mh = 1.67262171E-24  # mass of hydrogen
mu = 1.3
gamma = 1.67

plot_range = range(0,24,2) 

file_names = [] 
for i in plot_range:
    file_names.append("../../DD" + `i`.zfill(4) + "/" + "cloud_collision_" + `i`.zfill(4))

f = h5py.File("equilibrium_curve_0_01.hdf5", "r")
dens = f['/denstity'][:]
temp = f['/temperature'][:]
pres_k = dens * temp


for file in file_names:

    ds = yt.load(file)
    ds.add_field(('gas', 'pressure_k'), units="", function=pressure_k,
            display_name='P/k')

    my_sphere = ds.sphere("c", (75, "pc"))
    p = yt.PhasePlot(my_sphere, "temperature", "pressure", ["cell_mass"],
            weight_field=None)

    plot = p.plots[('gas', 'cell_mass')]
    plot.axes.plot(temp, dens*kboltz*temp/(mu*mh), "k--")

    plot.axes.set_xlim(1e2, 2e4)
    plot.axes.set_ylim(1e-12, 1e-8)
    p.save()
