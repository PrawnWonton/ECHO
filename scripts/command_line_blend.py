#Blender command line scripting attempt

import bpy
import numpy as np
import sys
import os

#number of frames per timestep
frame_count = 1

#initialize counter
i=0

#color change
color_change = 0.2/(len(sys.argv[5:]))

#print('the last argument is:'+str(sys.argv[5:]))

for fileName in sys.argv[5:len(sys.argv[5:])/3]:
    bpy.data.scenes["Scene"].frame_start = i
    bpy.data.scenes["Scene"].frame_end = i+frame_count
    text_name = '//'+str(fileName)
    bpy.data.textures["Texture"].voxel_data.filepath = str(text_name)


    for j in reversed(range(0,6)):
        bpy.data.textures["Texture"].color_ramp.elements[j].position += color_change
    #print(bpy.data.textures["Texture"].voxel_data.filepath)
    #print(text_name)

    bpy.ops.render.render(animation=True)

    i=i+frame_count+1
