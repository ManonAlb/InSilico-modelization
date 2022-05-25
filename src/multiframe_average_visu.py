# Visualization of the multiframe averaging

import matplotlib.pyplot as plt
import numpy as np
import os
from PIL import Image


#Create an array of the concerned images

path = '/home/manon/Documents/M1/internship/Ubiquitin/Ubiquitin_vacuum/'

# list files in img directory

files = os.listdir(path)
img_path = []
for file in files:

    if file.endswith(('.png')):
        #img_path = path + file
        img_path.append(os.path.join(path, file))


N = len(img_path)
s = 0
for i in range(0,N):
    s+= np.asarray(Image.open(img_path[i]),dtype=np.float)

mean = s/N
out = Image.fromarray(np.uint8(mean))

out.show()
