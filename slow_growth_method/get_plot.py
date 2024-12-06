import os
import matplotlib.pyplot as plt

def plot_barrier( list_x, list_y ):
	
	diff = round( max( list_y ) - min( list_y ), 2 )
	max_index = list_y.index( max( list_y ) )

	state_l = 0.9
	state_h = 3

	m_left = 0.18
	m_right = 0.98
	m_bottom = 0.17
	m_top = 0.99

	plt_h = 2.5
	plt_w = 3.5


	fig = plt.figure(figsize=(plt_w,plt_h))
	ax2 = fig.add_subplot(1, 1, 1 )
	ax2.plot( list_x, list_y, "o", markersize = 2)

	xmin = min( list_x )
	xmax = max( list_x )

	ax2.set_xlabel(r'O-H distance ($\mathrm{\AA}$)', fontsize = 12, labelpad = 0.5)
	ax2.set_ylabel( 'Relative free energy (eV)', fontsize = 12, labelpad = 2 )
	ymin = -0.1
	ymax = 1.1
	plt.vlines( x = list_x[ max_index ], ymin = 0, ymax = max( list_y ), color = 'red', linestyle='-', label=f'x = { list_x[ max_index ] }, y = { max( list_y ) }')
	mid_y = ( max(list_y ) + -0.2 ) / 2 
	plt.text( list_x[ max_index ] - 0.02, mid_y, f'{ diff } eV', verticalalignment = "center", horizontalalignment = "right", color = "black", fontsize = 12 )
	plt.axhline( 0, color = 'black', linestyle='--', linewidth = 2)
	#ax2.set_ylim( ymin, ymax)
	#yticks =  [0, 0.5, 1.0 ] 
	#ax2.set_yticks( yticks )
	#ax2.set_yticklabels( [ str(x) for x in yticks ] )
	#ax2.set_xlim( xmin, xmax)
	#xticks =  [ 1.0, 1.2, 1.4, 1.6, 1.8 ] 
	#ax2.set_xticks( xticks )
	#ax2.set_xticklabels( [ str( x ) for x in xticks ] )
	os.chdir( "/home/theodoros/PROJ_ElectroCat/theodoros/HER/Au/HER_Au/slow_grow_method/Na/pics" )
	plt.subplots_adjust(left=m_left, right=m_right, top=m_top, bottom=m_bottom, wspace=0.00, hspace= 0.0 )
	#plt.savefig( '3_NH4_hydration.png', dpi = 600, transparent = True )
	plt.show()
