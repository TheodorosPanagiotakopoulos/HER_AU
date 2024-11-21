import get_data
import get_plot

path = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/tutorial_2/"

x, y = get_data.get_free_energy( path )
get_plot.plot_barrier( x, y )
