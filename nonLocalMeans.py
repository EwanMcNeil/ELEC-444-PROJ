import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 
import math
from scipy.spatial import distance

from random import seed
from random import randint

#The classical definition of the NL-means filter considers that each voxel can be linked to all the others, but#
# for practical computational reasons the number of voxels taken into account in the weighted
#average can be limited to the so-called “search volume” Vi of size (2M+1)3, centered at the current voxel xi.


#making the weight function

def weight(tupleI, tupleJ):
    Z = 216
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
    for x in range(inputTuple[0]-3, inputTuple[0]+3):
        for y in range(inputTuple[1]-3, inputTuple[1]+3):
                for z in range(inputTuple[2]-3, inputTuple[2]-3):
                    try:
                        w = weight(inputTuple,(x,y))
                        #,z))
                    
                        sum = sum + w*data[x,y] 
                        #,z]
                    except IndexError:
                        break

    
    return sum


def createImage(size):
    seed(1)
    output = np.zeros(shape=(size,size,size));
    for x in range(size):
        for y in range(size):
            for z in range(size):
                output[x,y,z] = randint(0, 10)

    return output



##code for guassian noise
##found online make sure to change later**
def noisy(image):
      row,col,ch= image.shape
      mean = 0
      var = 0.1
      sigma = var**0.5
      gauss = np.random.normal(mean,sigma,(row,col,ch))
      gauss = gauss.reshape(row,col,ch)
      noisy = image + gauss
      return noisy


img = nib.load('NormalBrains/t1_icbm_normal_1mm_pn3_rf20.mnc')
data = img.get_fdata()
print(data.mean())


#A voxel represents a single sample, or data point, on a regularly spaced, three-dimensional grid.



output = np.zeros(shape=(6,6,6));
#going through each value and caluclulating their new thing


groundImage = createImage(6)

noisyImage = noisy(groundImage)


for x in range(6):
 for y in range(6):
    for z in range(6):
        output[x,y,z] = getNewValue((x,y,z),noisyImage)
        print(x)
      




plt.imshow(groundImage[:,:,3], interpolation = 'nearest')
#axs[0].set_title('inputData')
plt.show()




plt.imshow(noisyImage[:,:,3], interpolation = 'nearest')
#axs[0].set_title('noisy')
plt.show()


plt.imshow(output[:,:,3], interpolation = 'nearest')
#axs[0].set_title('noiser')
plt.show()

