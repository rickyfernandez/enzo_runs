import yt
import os
import numpy as np
import matplotlib.pyplot as plt

pc_to_cm = 3.0857E18 # 1pc in cm
sec_to_Myr = 60.0*60.0*24.0*365*1.0E6
plot_range = [0, 3, 6, 9, 12, 15]

num_plots = len(plot_range)
colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0,0.9, num_plots)])

file_names = []
for i in plot_range:
    file_names.append("DD" + `i`.zfill(4) + "/" + "cloud_collision_" + `i`.zfill(4))

labels = []
for file in file_names:

    ds = yt.load(file)
    sphere = ds.sphere("c", (75, "pc"))
    plot = yt.ProfilePlot(sphere, "radius", ["Cooling_Time"],
            weight_field="cell_mass")
    profile = plot.profiles[0]
    plt.loglog(profile.x/pc_to_cm, profile["Cooling_Time"]/sec_to_Myr)
    labels.append(r"%0.2f Myr" % ds.current_time.value)

plt.xlabel(r"Radius $(\mathrm{pc})$")
plt.ylabel(r"Cooling Time $(Myr)$")
plt.xlim(1,75)
plt.legend(labels, loc="lower right", frameon=False, prop={'size':11})
plt.savefig("cooling_series")
