import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import gzip as gz
import matplotlib.image as mpimg
from tqdm import tqdm 
img = mpimg.imread(  'frame_0.png' )
ny = img.shape[ 0 ]
nx = img.shape[ 1 ]
lx = 4
ly = lx/nx * ny

def add_time_stamp( image_name, outdir = 'outdir' ):
    t = int( image_name.replace( 'frame_','' ).replace( '.png','' ) )
    fig = plt.figure( figsize = ( lx, ly ) )
    ax  = fig.add_subplot( 1, 1, 1  )
    img =  mpimg.imread(  image_name )
    ax.imshow( img )
    time_text = str( round( t/1000, 3 ) )
    if len( time_text ) == 3:
        time_text += '00'
    elif len( time_text ) == 4:
        time_text += '0'
    ax.text( 1, 1, 't = ' + time_text  + ' ps', color = 'k', fontsize = 12, fontweight = 'bold', va = 'top', ha = 'left', bbox=dict(facecolor='white', edgecolor='none', alpha = 0.3, pad=1.0) )
    ax.axis( 'off' )
    plt.subplots_adjust(left=0.0, right=1., top=1, bottom=0.0, wspace=0.00, hspace= 0.00)
    if not os.path.isdir( outdir ):
        os.mkdir( outdir )
    plt.savefig( outdir + '/t' + image_name, dpi = 200 )
    plt.clf()
    plt.close()
frames = [ x for x in os.listdir( '.' ) if '.png' in x and 'frame' in x ]
for image_name in tqdm( frames ):
    if ( not os.path.isfile( 'outdir/t' + image_name)):
        add_time_stamp( image_name )
