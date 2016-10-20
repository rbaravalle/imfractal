import vtk
import numpy as np
import os
import matplotlib.pyplot as plt
import PIL
import Image

DEBUG =True
#directory="splitted_mri/"

#l = []

#k=0 #add the next picture in a differente level of depth/z-positions
#for file in os.listdir(directory):
#    img = directory + file
#    if DEBUG : print img
#    l.append(img)
# the os.listdir function do not give the files in the right order 
#so we need to sort them
#l=sorted(l)

#temp = Image.open(l[0])
#h, w = temp.size
#d = len(l)*5 #with our sample each images will be displayed 5times to get a better view
#if DEBUG : print 'width, height, depth : ',w,h,d

#stack = np.zeros((w,d,h),dtype=np.uint8)
stack = np.load("../exps/data/bone_sample.npy")

#for i in l:
#    im = Image.open(i)
#    temp = np.asarray(im, dtype=int)
#    for i in range(5):
#        stack[:,k+i,:]= temp
#    k+=5
    #~ stack[:,k,:]= temp
    #~ k+=1

if DEBUG :
    res = np.amax(stack)
    print 'max value',res
    res1 = np.amin(stack)
    print 'min value',res1

#convert the stack in the right dtype
stack = np.min(stack) + stack
#stack = np.clip(stack, np.amin(stack),1500)


if DEBUG :
    res = np.amax(stack)
    print 'max value',res
    res1 = np.amin(stack)
    print 'min value',res1

stack = 255*(stack / np.amax(stack))


if DEBUG :
    res = np.amax(stack)
    print 'max value',res
    res1 = np.amin(stack)
    print 'min value',res1


stack = np.require(stack,dtype=np.uint8)

if DEBUG :
    res = np.amax(stack)
    print 'max value',res
    res1 = np.amin(stack)
    print 'min value',res1

if DEBUG :#check if the image have not been modified
    test = stack [:,:,0]
    plt.imshow(test,cmap='gray')
    plt.show()

if DEBUG : print 'stack shape & dtype' ,stack.shape,',',stack.dtype

dataImporter = vtk.vtkImageImport()
data_string = stack.tostring()

dataImporter.CopyImportVoidPointer(data_string, len(data_string))
dataImporter.SetDataScalarTypeToUnsignedChar()
dataImporter.SetNumberOfScalarComponents(1)

#vtk uses an array in the order : height, depth, width which is 
#different of numpy (w,h,d) 
w, d, h = stack.shape
dataImporter.SetDataExtent(0, h-1, 0, d-1, 0, w-1)
dataImporter.SetWholeExtent(0, h-1, 0, d-1, 0, w-1)

alphaChannelFunc = vtk.vtkPiecewiseFunction()
colorFunc = vtk.vtkColorTransferFunction()

l = 0.1
for i in range(256):
    #alphaChannelFunc.AddPoint(i, i/127.0)
    if i > 1:
        #print ""
    	colorFunc.AddRGBPoint(i,l*0.5,l*0.5,l*0.4)
        alphaChannelFunc.AddPoint(i, i/127.0)
    else:
        alphaChannelFunc.AddPoint(i, 0.05)
        colorFunc.AddRGBPoint(i,0.8,0.8,0.8)

# for our test sample, we set the black opacity to 0 (transparent) so as
#to see the sample  
alphaChannelFunc.AddPoint(255, 0.0)
colorFunc.AddRGBPoint(255,1.0,1.0,1.0)

volumeProperty = vtk.vtkVolumeProperty()
volumeProperty.SetColor(colorFunc)
#volumeProperty.ShadeOn()
volumeProperty.SetScalarOpacity(alphaChannelFunc)

# This class describes how the volume is rendered (through ray tracing).
compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()
# We can finally create our volume. We also have to specify the data for
# it, as well as how the data will be rendered.
volumeMapper = vtk.vtkVolumeRayCastMapper()
# function to reduce the spacing between each image
volumeMapper.SetMaximumImageSampleDistance(0.01)

volumeMapper.SetVolumeRayCastFunction(compositeFunction)
volumeMapper.SetInputConnection(dataImporter.GetOutputPort())

# The class vtkVolume is used to pair the preaviusly declared volume as 
#well as the properties to be used when rendering that volume.
volume = vtk.vtkVolume()
volume.SetMapper(volumeMapper)
volume.SetProperty(volumeProperty)

# With almost everything else ready, its time to initialize the renderer and window,
# as well as creating a method for exiting the application
renderer = vtk.vtkRenderer()
renderWin = vtk.vtkRenderWindow()
renderWin.AddRenderer(renderer)
renderInteractor = vtk.vtkRenderWindowInteractor()
renderInteractor.SetRenderWindow(renderWin)

# We add the volume to the renderer ...
renderer.AddVolume(volume)
# ... set background color to white ...
renderer.SetBackground(0.2,0.2,0.2)
# ... and set window size.
renderWin.SetSize(550, 550)
renderWin.SetMultiSamples(4)
# A simple function to be called when the user decides to quit the application.
def exitCheck(obj, event):
    if obj.GetEventPending() != 0:
        obj.SetAbortRender(1)

# Tell the application to use the function as an exit check.
renderWin.AddObserver("AbortCheckEvent", exitCheck)

#to auit, press q
renderInteractor.Initialize()
# Because nothing will be rendered without any input, we order the first
# render manually before control is handed over to the main-loop.
renderWin.Render()
renderInteractor.Start()
