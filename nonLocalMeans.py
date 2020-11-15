import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 
import math
from scipy.spatial import distance

import imageio


  

from random import seed
from random import randint

#The classical definition of the NL-means filter considers that each voxel can be linked to all the others, but#
# for practical computational reasons the number of voxels taken into account in the weighted
#average can be limited to the so-called “search volume” Vi of size (2M+1)3, centered at the current voxel xi.


#making the weight function


##needs to be done with the neighboorhood intensities 
def neighboorhoodIntensities(Tuple, data):
    
    
    #M defines the half length of neighboorhood
    output = np.array([])
    count = 0
    for x in range(Tuple[0]-1, Tuple[0]+1):
        for y in range(Tuple[1]-1, Tuple[1]+1):
            for z in range(Tuple[2]-1, Tuple[2]+1):
                try:
                    
                    output = np.append(output, [data[x,y,z]])
                    count = count +1
                except IndexError:
                    output[1,count] = 0
                    count = count +1
                    continue
        
  
   
    return output






def weight(tupleI, tupleJ,data):
    Z = 216
    h = 0.5
    #h should be 10 times the standard deviation
    #no idea what Z should be 
    uNJ = neighboorhoodIntensities(tupleJ, data)
    uNI = neighboorhoodIntensities(tupleI, data)

    dist = np.linalg.norm(uNJ-uNI)
    div = dist/(h*h)
    exponential = math.exp((-1)*div)

    output = (1/Z)*exponential
    return output




def getNewValue(inputTuple, data):

    #need to center around value 
    #with a certain thing but we dont need to go into that rn
    sum = 0;
    global M 
    for x in range(inputTuple[0]-M, inputTuple[0]+M):
        for y in range(inputTuple[1]-M, inputTuple[1]+M):
             for z in range(inputTuple[2]-M, inputTuple[2]+M):
                 try:
                     if(inputTuple[0]-M < -2 or inputTuple[1]-M < -2 or inputTuple[2]-M < -2):
                        sum = sum + 0
                     elif(inputTuple[0]+M > 21 or inputTuple[1]+M > 21 or inputTuple[2]+M > 21):
                        sum = sum + 0
                     else:
                        w = weight(inputTuple,(x,y,z),data)
                        sum = sum + w*data[x,y,z] 
                 except IndexError:
                     sum = sum + 0
                     continue
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
      mean = 128
      var = 50
      sigma = var**0.5
      gauss = np.random.normal(mean,sigma,(row,col,ch))
      gauss = gauss.reshape(row,col,ch)
      noisy = image + gauss
      return noisy


img = nib.load('NormalBrains/t1_icbm_normal_1mm_pn3_rf20.mnc')
data = img.get_fdata()
print(data.mean())


#A voxel represents a single sample, or data point, on a regularly spaced, three-dimensional grid.


##GLobal M to define the neighborhood
M = 3


output = np.zeros(shape=(20,20,20));
#going through each value and caluclulating their new thing

im = imageio.imread('heart.png')
print(im.shape)


imageIN = im[:,:,2]

img3D = np.zeros(shape=(20,20,20))


for x in range(20):
  for y in range(20):
     for z in range(20):
        img3D[x,y,z] = imageIN[x,y]
      
groundImage = createImage(6)

noisyImage = noisy(img3D)


for x in range(20):
  for y in range(20):
     for z in range(20):
        output[x,y,z] = getNewValue((x,y,z),noisyImage)
        print(x)

      



plt.subplot(1,3,1)
plt.imshow(img3D[:,:,3], interpolation = 'nearest')
plt.title("ground")
print(img3D.shape)
print(noisyImage.shape)
print(output.shape)


plt.subplot(1,3,2)
plt.imshow(noisyImage[:,:,3], interpolation = 'nearest')
plt.title("noisy")


plt.subplot(1,3,3)
plt.imshow(output[:,:,3], interpolation = 'nearest')
plt.title("filtered")
plt.show()

