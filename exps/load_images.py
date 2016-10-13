from PIL import Image
import os, os.path
import numpy as np


path = "./"

print "Creating imgs (dat)..."

command = "r mf_image_gen_rmp3D.R"
os.system(command)

print "Converting imgs dat -> tif..."

valid_images = [".dat"]

for f in os.listdir(path):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_images:
        continue
    command = "python convert_img_ascii-tif.py " + os.path.join(path,f)
    os.system(command)

print "tif -> npy..."

imgs = np.zeros((128,128,128))
path = "./"
valid_images = [".tif"]
i = 0
for f in os.listdir(path):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_images:
        continue
    imgs[i] = np.asarray(Image.open(os.path.join(path,f)).convert("L"))
    i = i+1

print "Deleting imgs (dat) and (tif)..."

np.save('img3d.npy', imgs)

command = "rm *.dat"
os.system(command)

command = "rm *.tif"
os.system(command)
