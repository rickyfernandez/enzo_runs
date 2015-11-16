import yt
import numpy as np

from yt.utilities.physical_constants import mh, kboltz

pc_to_cm = 3.086E18  # parsec
G  = 6.67259E-8      # gravity
mu = 1.3
gamma = 1.67

plot_range = range(0,14,2) 
plot_range.append(13)

file_names = [] 
for i in plot_range:
    file_names.append("../../DD" + `i`.zfill(4) + "/" + "cloud_collision_" + `i`.zfill(4))


for file in file_names:

    ds = yt.load(file)
    my_sphere = ds.sphere("c", (75, "pc"))
    p = yt.PhasePlot(my_sphere, "temperature", "pressure", ["cell_mass"],
            weight_field=None)

    plot = p.plots[('gas', 'cell_mass')]

    plot.axes.set_xlim(1e2, 2e4)
    plot.axes.set_ylim(1e-16, 1e-8)
    p.save()
