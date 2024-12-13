import os
import get_data
import get_plot

path = os.getcwd()

x, y = get_data.get_free_energy( path )
get_plot.plot_barrier( x, y )
