import yt

ds = yt.load("")
my_sphere = ds.sphere("c", )
plot = yt.PhasePlot(my_sphere, "density", "temperature", ["cell_mass"],
        weight_field=None)
plot.save()
