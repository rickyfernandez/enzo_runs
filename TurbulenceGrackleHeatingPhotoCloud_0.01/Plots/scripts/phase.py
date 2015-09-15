import yt
import numpy as np

pc_to_cm = 3.086E18  # parsec
G  = 6.67259E-8      # gravity
k  = 1.380658E-16    # boltzmann
mh = 1.67262171E-24  # mass of hydrogen
mu = 1.3
gamma = 1.67

fac = gamma * k /(mu*mh*G)

jeans_length = 4 * 150.0 * pc_to_cm/(128 * 2**4)

plot_range = range(0,24,3) 
plot_range.append(23)

file_names = [] 
for i in plot_range:
    file_names.append("../../DD" + `i`.zfill(4) + "/" + "cloud_collision_" + `i`.zfill(4))

for file in file_names:

    ds = yt.load(file)
    my_sphere = ds.sphere("c", (75, "pc"))
    p = yt.PhasePlot(my_sphere, "density", "temperature", ["cell_mass"],
            weight_field=None)

    plot = p.plots[('gas', 'cell_mass')]

    x = np.arange(-22,-14)
    rho = np.power(10.0, x)
    tem = rho*jeans_length**2/fac

    plot.axes.plot(rho,tem, "k--")

    plot.axes.set_xlim(1e-24, 1e-15)
    plot.axes.set_ylim(1, 1e5)
    p.save()
