import yt
import numpy as np
import matplotlib.pyplot as plt

plot_range = [0, 10, 20, 30, 40, 50, 60, 70, 80]

num_plots = len(plot_range)
colormap = plt.cm.gist_ncar
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0,0.9, num_plots)])

file_names = []
for i in plot_range:
    file_names.append("../../DD" + `i`.zfill(4) + "/" + "cloud_collision_" + `i`.zfill(4))

labels = []
for file in file_names:

    ds = yt.load(file)
    sphere = ds.sphere('max', (75, 'pc'))

    # calculate and stor the bulk velocity of the sphere
    bulk_velocity = sphere.quantities.bulk_velocity()
    sphere.set_field_parameter('bulk_velocity', bulk_velocity)

    # create 1d profile object for profiles over radius
    # and add a velocity profile.
    prof = yt.create_profile(sphere, 'radius', ('gas', 'velocity_magnitude'),
        units = {'radius': 'pc'},
        extrema = {'radius': ((0.1, 'pc'), (75, 'pc'))},
        weight_field ='cell_mass')

    # create arrays to plot
    radius = prof.x.value
    mean = prof['gas', 'velocity_magnitude'].value
    variance = prof.variance['gas', 'velocity_magnitude'].value

    # plot the variance of the velocity magnitude
    plt.loglog(radius, variance*1.0E-5, label="Standarad Deviation")
    labels.append(r"%0.2f Myr" % ds.current_time.value)

plt.xlabel(r"Radius $(\mathrm{pc})$")
plt.ylabel('v [km/s]')
plt.xlim(1,75)
plt.legend(labels, loc="lower left", frameon=False)
plt.savefig('velocity_disperssion_profiles.png')
