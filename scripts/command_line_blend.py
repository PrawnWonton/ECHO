#Blender command line scripting attempt

import bpy
import numpy as np
import sys
import os

#initialize counter
i=0

for fileName in sys.argv[-1:]:
    bpy.data.scenes["Scene"].frame_start = i
    bpy.data.scenes["Scene"].frame_end = i+5
    bpy.data.textures["dens.002"].voxel_data.filepath = fileName
    print(i)
    print(fileName)
    i=i+1
