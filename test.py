import Image
import time
import src.Algorithm.Singularity as Singularity
import src.Algorithm.Sandbox as Sandbox

filename = 'images/baguette2.tif'

i = Singularity.Singularity(20)

print "Calculating Singularity Multifractal Spectrum..."
t =  time.clock()
fds = i.getFDs(filename)
t =  time.clock()-t
print "Time Singularity: ", t
print fds

i = Sandbox.Sandbox(14)
i.setDef(40,1.15)

print "Calculating Sandbox Multifractal Spectrum..."
t =  time.clock()
fds2 = i.getFDs(filename)
t =  time.clock()-t
print "Time Sandbox: ", t
print fds2

