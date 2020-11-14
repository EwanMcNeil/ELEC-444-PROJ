import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 
import math
from scipy.spatial import distance

#The classical definition of the NL-means filter considers that each voxel can be linked to all the others, but#
# for practical computational reasons the number of voxels taken into account in the weighted
#average can be limited to the so-called “search volume” Vi of size (2M+1)3, centered at the current voxel xi.


#making the weight function

def weight(tupleI, tupleJ):
    Z = 100
    h = 10
    #h should be 10 times the standard deviation
    #no idea what Z should be 

    dist = distance.euclidean(tupleI,tupleJ)

    div = dist/(h*h)
    exponential = math.exp(div)

    output = (1/Z)*exponential

    return output




def getNewValue(inputTuple, data):

    #need to center around value 
    #with a certain thing but we dont need to go into that rn

    sum = 0;
    for x in range(inputTuple[0]-5, inputTuple[0]+5):
        for y in range(inputTuple[1]-5, inputTuple[1]+5):
                 # for z in range(inputTuple[0]-10, inputTuple[0]-10):
            try:
                w = weight(inputTuple,(x,y))
                #,z))
               
                sum = sum + w*data[x,y] 
                #,z]
            except IndexError:
                break

    
    return sum

img = nib.load('NormalBrains/t1_icbm_normal_1mm_pn3_rf20.mnc')
data = img.get_fdata()
print(data.mean())


#A voxel represents a single sample, or data point, on a regularly spaced, three-dimensional grid.



output = np.zeros(shape=(181,217));
#going through each value and caluclulating their new thing

oneSlice = data[:,:,50]

for x in range(181):
  for y in range(217):
   # for z in range(181):
        output[x,y] = getNewValue((x,y),oneSlice)
        print(x)
      

plt.imshow(data[:,:,50], interpolation = 'nearest')
plt.show()

plt.imshow(output[:,:], interpolation = 'nearest')
plt.show()
