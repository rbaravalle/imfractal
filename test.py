"""
Copyright (c) 2013 Rodrigo Baravalle
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

if(False):
    import os
    arr = [ "imfractal/Algorithm/qs3D", "imfractal/Algorithm/qs"]

    print "WTF"

    for i in range(len(arr)):

        command1 = "cython "+arr[i]+".pyx "
        command2 = "gcc -c -fPIC -I/usr/include/python2.7/ "+arr[i]+".c"+" -o "+arr[i]+".o"
        command3 = "gcc -shared "+arr[i]+".o -o "+arr[i]+".so"

        print command1
        os.system(command1)
        print command2
        os.system(command2)
        print command3
        os.system(command3)
    exit()

#import tests.test_MFS as tmfs
#except ImportError:
#    print "Error: the module tests/tests_MFS is missing"

if(False):
    import tests.test_bones as ts
    ts.do_test()
    exit()

    import tests.test_bonesFelix as ts
    ts.do_test()

    #exit()

if(False):
    import tests.test_bonesR2 as ts
    ts.do_test()

    import tests.test_bonesBVTV as tb
    tb.do_test()

#import tests.test_bonesBVTVspearman as tb
#tb.do_test()

import tests.test_MFS as tmfs
tmfs.do_test()
exit()

import tests.test_bones_BioAsset as tbba
tbba.do_test("/home/rodrigo/members.imaglabs.org/felix.thomsen/Rodrigo/BioAsset/mats/")

#import tests.test_boneMeasuresSpearman as tb
#tb.do_test()

exit()

import tests.testcomparison as tcomp
tcomp.do_test()

#import tests.test_real_fake as trealfake
#except ImportError:
#    print "Error: the module tests/tests_classifier is missing"

#import tests.test_boxd as tbox
#tmfs.do_test()
#trealfake.do_test()
#tbox.do_test()
