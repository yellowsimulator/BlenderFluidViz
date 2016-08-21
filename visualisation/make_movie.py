
# read from file and sort

from natsort import natsorted
import sys
import os
path = sys.argv[1]
import shutil
import Image
im = Image.open(os.path.isdir(path))
if os.path.isdir(path):
    images = natsorted([f for f in os.listdir(path)])
    if not 'sorted_{}'.format(path):
        new_path = os.mkdir('sorted_{}'.format(path))
        shutil.copytree(images, dst)
else:
    print path, " is not a directory "
    

print im


