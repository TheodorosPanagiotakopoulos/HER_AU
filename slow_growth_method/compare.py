import matplotlib.pyplot as plt

def plot_barrier( list_y, x_labels, yticks, name ):
	path = "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/Na/pics_Na/"
	state_l = 0.9
	state_h = 3

	m_left = 0.18
	m_right = 0.98
	m_bottom = 0.17
	m_top = 0.99

	plt_h = 2.5
	plt_w = 3.5


	x_values = [ 1, 11 ]
	fig = plt.figure(figsize=(plt_w,plt_h))
	ax2 = fig.add_subplot(1, 1, 1 )
	plt.bar( x_values, list_y)
	x_filtered = [x for x, y in zip(x_values, list_y) if y != 0]
	y_filtered = [y for y in list_y if y != 0]
	plt.plot( x_filtered, y_filtered, marker='o', color='red', linestyle='-', linewidth=2)
	xmax = max( x_values ) + 4
	xmin = min( x_values ) - 2
	ax2.set_xlabel(r'Categories', fontsize = 12, labelpad = 0.5)
	ax2.set_ylabel( 'Relative free energy (eV)', fontsize = 12, labelpad = 2 )
	ymin = -0.026
	ymax = max( list_y ) + 0.1 #1.3
	ax2.set_ylim( ymin, ymax)
	#yticks =  [0.0, 0.4, 0.8, 1.2 ] 
	ax2.set_yticks( yticks )
	ax2.set_yticklabels( [ str(x) for x in yticks ] )
	ax2.set_xlim( xmin, xmax)
	xticks = x_values
	ax2.set_xticks( xticks )
	ax2.set_xticklabels( [ str( x ) for x in x_labels ], rotation = 0, fontsize = 10, ha = "center"  )
	plt.subplots_adjust(left=m_left, right=m_right, top=m_top, bottom=m_bottom, wspace=0.00, hspace= 0.0 )
	plt.savefig( path + name, dpi = 600, transparent = True )
	plt.show()
