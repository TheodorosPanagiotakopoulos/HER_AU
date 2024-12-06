import matplotlib.pyplot as plt

def plot_barrier( list_x, list_y ):
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
	ax2.plot( list_x, list_y, "o", markersize = 4)

	xmax = -0.155
	xmin = -0.195
	ax2.set_xlabel(r'Na-H distance ($\mathrm{\AA}$)', fontsize = 12, labelpad = 0.5)
	ax2.set_ylabel( 'Relative free energy (eV)', fontsize = 12, labelpad = 2 )
	ymin = -0.026
	ymax = 0.16
	#ax2.set_ylim( ymin, ymax)
	#yticks =  [0, 0.05, 0.1, 0.15 ] 
	#ax2.set_yticks( yticks )
	#ax2.set_yticklabels( [ str(x) for x in yticks ] )
	#ax2.set_xlim( xmin, xmax)
	#xticks =  [ -0.2, -0.18, -0.16] 
	#ax2.set_xticks( xticks )
	#ax2.set_xticklabels( [ str( x ) for x in xticks ] )
	#ax2.hlines(0, -0.2, xmax, lw = 0.5, color = 'gray' )
	#ax2.text(-0.175,  0.004 , '$\Phi$ = -1.4 (V vs RHE)', fontsize = 9)
	plt.subplots_adjust(left=m_left, right=m_right, top=m_top, bottom=m_bottom, wspace=0.00, hspace= 0.0 )
	#plt.savefig( '1_Na_free_energy.png', dpi = 600 )
	plt.show()

