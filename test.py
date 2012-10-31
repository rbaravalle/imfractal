import Image

import src.Algorithm.Singularity as Singularity
import src.Algorithm.Sandbox as Sandbox

filename = 'images/baguette2.tif'

i = Singularity.Singularity()
i.setCuantas(20)
fds = i.getFDs(filename)

i = Sandbox.Sandbox()
i.setDef(40,1.15)
fds2 = i.getFDs(filename)

print fds,fds2

