#!/bin/bash


myArray=( "$@" )
directory="${myArray[0]}"
cd $directory
#convert -quality 100 *.png outvideo.mpeg
#ffmpeg -f image2 -r 1/4 -i tmp_%04d.png -vcodec mpeg4 -y movie.mp4
#file_name="${myArray[0]:0:(-1)}"
##echo $file_name
#mencoder "mf://*.png" -mf fps=2 -o animation.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=500
mencoder "mf://*.png" -mf fps=2 -o animation.ogg -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=500
mencoder "mf://*.png" -mf fps=2 -o animation.mp4 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=500

