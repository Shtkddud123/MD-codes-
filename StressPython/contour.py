import matplotlib, matplotlib.pyplot, numpy
import MDAnalysis, MDAnalysis.visualization.streamlines

u1, v1, average_displacement,standard_deviation_of_displacement = MDAnalysis.visualization.streamlines.generate_streamlines('testing.gro','testing_filtered.xtc',grid_spacing = 20, MDA_selection = 'name PO4',start_frame=2,end_frame=3,xmin=-8.73000049591,xmax= 1225.96008301,ymin= -12.5799999237, ymax=1224.34008789,maximum_delta_magnitude = 1.0,num_cores=16)
x = numpy.linspace(0,1200,61)
y = numpy.linspace(0,1200,61)
speed = numpy.sqrt(u1*u1 + v1*v1)
fig = matplotlib.pyplot.figure()
ax = fig.add_subplot(111,aspect='equal')
ax.set_xlabel('x ($\AA$)')
ax.set_ylabel('y ($\AA$)')
ax.streamplot(x,y,u1,v1,density=(10,10),color=speed,linewidth=3*speed/speed.max())
fig.savefig('testing_streamline.png',dpi=300)
