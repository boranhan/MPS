import os
import pickle
import numpy
import scipy.io
import sys

if (len(sys.argv) != 2):
    print "Usage <mosaic_file>"
    exit()

directory = os.path.dirname(sys.argv[1])
if (len(directory) == 0):
    directory = "."

fp = open(sys.argv[1], "r")

file_number = 1
while 1:
    image_name = fp.readline().rstrip()
    if not image_name: break
	
    image_name = image_name[6:]
	
    print "converting:", image_name

    image_dict = pickle.load(open(directory + "/" + image_name))

    mat_dict = {}
    for key in image_dict:
        val_type = type(image_dict[key])
        if (val_type in [type(""), type(0), type(0.0), type(numpy.array([]))]):
            #print key, type(image_dict[key])
            mat_dict[key] = image_dict[key]
        else:
            #print key, str(image_dict[key])
            mat_dict[key] = str(image_dict[key])

    scipy.io.savemat(directory + "/" + image_name[:-4] + ".mat", mat_dict)

    file_number += 1

fp.close()
