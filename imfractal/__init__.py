"""
Copyright (c) 2016 Rodrigo Baravalle
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

#from Path import Path
import Algorithm.Sandbox

cython = False

if(cython):
    import Algorithm.CSandbox
    import Algorithm.CSandbox3D

#import Algorithm.SandboxCL
import Algorithm.MFS
import Algorithm.Singularity
import Algorithm.Boxdimension
import Algorithm.MFS_3D
import Algorithm.Local_MFS_3D
import Algorithm.Local_MFS_Pyramid_3D
import Algorithm.MFS_3D_Slices
import Algorithm.Stats_MFS_3D


#SandboxCL = Algorithm.SandboxCL.SandboxCL
Sandbox = Algorithm.Sandbox.Sandbox

if(cython):
    CSandbox = Algorithm.CSandbox.CSandbox
    CSandbox3D = Algorithm.CSandbox3D.CSandbox3D

MFS = Algorithm.MFS.MFS
MFS_3D = Algorithm.MFS_3D.MFS_3D
MFS_3D_Slices = Algorithm.MFS_3D_Slices.MFS_3D_Slices
Local_MFS_3D = Algorithm.Local_MFS_3D.Local_MFS_3D
Stats_MFS_3D = Algorithm.Stats_MFS_3D.Stats_MFS_3D
Local_MFS_Pyramid_3D = Algorithm.Local_MFS_Pyramid_3D.Local_MFS_Pyramid_3D
Singularity = Algorithm.Singularity.Singularity
Boxdimension = Algorithm.Boxdimension.Boxdimension

# Global Variables

MFS_HOLDER = True
LOCAL = True
# only one of the following equals true, or none
APPLY_LAPLACIAN = False
APPLY_GRADIENT = False
# 2.5 D
SLICES_MFS = False
Stats_MFS = False



TRANSFORMED_INPUT_STR = ''
MFS_STR = ''
ADAPTIVE_STR = ''

if APPLY_LAPLACIAN:
    TRANSFORMED_INPUT_STR = '_laplacian'
else:
    if APPLY_GRADIENT:
        TRANSFORMED_INPUT_STR = '_gradient'

MFS_STR = ''

if MFS_HOLDER :
    MFS_STR = '_holder'
    ADAPTIVE_STR = ''

else:
    ADAPTIVE_STR = '_adaptive_0.75'
    ADAPTIVE_STR = '_normalized_input'

if LOCAL:
    MFS_STR = MFS_STR+'_local'

data_path = "exps/data/"

BASE_NAME = 'mfs' + MFS_STR + TRANSFORMED_INPUT_STR + '_BioAsset' + ADAPTIVE_STR


#REIL_Path = lambda trace_filename,first,last: Path(trace_filename, ReilParser, first, last).
