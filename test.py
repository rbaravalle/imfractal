import Image

import src.Algorithm.Singularity as Singularity
import src.Algorithm.Sandbox as Sandbox

filename = '/media/5f4f962d-ed66-4166-bfb2-da773f9a77cb/rodrigo/mecom2012/mecom/imagenes/scanner/baguette/baguette1.tif'

i = Singularity.Singularity()
i.setCuantas(20)
fds = i.getFDs(filename)

i = Sandbox.Sandbox()
i.setDef(40,1.15)
fds2 = i.getFDs(filename)

print fds,fds2

