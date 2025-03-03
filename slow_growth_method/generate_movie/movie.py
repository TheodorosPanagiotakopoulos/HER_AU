import os

image_folder='./outdir/'
#image_files = sorted( [img for img in os.listdir(image_folder) if img.endswith(".png")] )
image_files = []
for i in range( 30000 ):
    if os.path.isfile( image_folder + 'tframe_' + str(i) +'.png' ):
        image_files.append( 'tframe_' + str(i) +'.png' )
image_list  = [ os.path.join(image_folder,img) for img in image_files]
print ( image_list)

fps = 10
output = "MD_"
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_list, fps=fps)
clip.write_videofile(output + '.mp4')
