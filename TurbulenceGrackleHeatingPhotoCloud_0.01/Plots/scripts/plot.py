import yt

plot_range = range(0,24,3) 
plot_range.append(23)

file_names = [] 
for i in plot_range:
    file_names.append("../../DD" + `i`.zfill(4) + "/" + "cloud_collision_" + `i`.zfill(4))

for file in file_names:
    ds = yt.load(file)
    #p = yt.SlicePlot(ds, "z", "velocity_magnitude")
    p = yt.SlicePlot(ds, "z", "temperature")
    #p.set_zlim("number_density", 5.0E-2, 1.0E3)
    #p.set_zlim("pressure", 3.0E-13, 4.0E-10)
    p.set_zlim("temperature", 1.0E1, 2.0E6)
    #p.annotate_grids()
    p.save()
