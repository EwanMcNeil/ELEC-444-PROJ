import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 
import math
from scipy.spatial import distance
import concurrent.futures

import imageio

from sigma import sigma

  

from random import seed
from random import randint

#The classical definition of the NL-means filter considers that each voxel can be linked to all the others, but#
# for practical computational reasons the number of voxels taken into account in the weighted
#average can be limited to the so-called “search volume” Vi of size (2M+1)3, centered at the current voxel xi.


#making the weight function


##they use the connected neighborhood not the 26 one
##needs to be done with the neighboorhood intensities 



def neighboorhoodIntensities(tuple, data):
    
    
    #M defines the half length of neighboorhood
    output = np.array([])
    vox1 = (tuple[0]-1,tuple[1],tuple[2])
    vox2 = (tuple[0]+1,tuple[1],tuple[2])
    vox3 = (tuple[0],tuple[1],tuple[2]+1)
    vox4 = (tuple[0],tuple[1],tuple[2]-1)
    vox5 = (tuple[0],tuple[1]+1,tuple[2])
    vox6 = (tuple[0],tuple[1]-1,tuple[2])
    iteration = (vox1,vox2,vox3,vox4,vox5,vox6)

    for vox in iteration:
        try:
            output = np.append(output, [data[vox]])
        except IndexError:
            output = np.append(output, 0)
            continue
    

    return output





def weight(tupleI, tupleJ,data):
    Z = 216
    global h


    #doing part two with the voxel selection 
    #dont know the varience so im just doing the mean
    
    #uNJ = neighboorhoodIntensities(tupleJ, data)
    #uNI = neighboorhoodIntensities(tupleI, data)

    #meanJ = np.sum(uNJ)/len(uNJ)
    #meanI = np.sum(uNI)/len(uNJ)

    #meanDiv = meanJ/meanI
    global blockMeansVal
    uNJ = blockMeansVal[tupleJ]
    uNI = blockMeansVal[tupleI]

    meanDiv = uNI/uNI
    print(uNJ,uNI,meanDiv)

    output = 0

    if(meanDiv < 1):
        dist = np.linalg.norm(uNJ-uNI)
        div = dist/(h*h)
        
        #div = dist/(2*1*h*26)
        exponential = math.exp((-1)*div)
    
        output = (1/Z)*exponential
   
    return output




def getNewValue(inputTuple, data):

    #need to center around value 
    #with a certain thing but we dont need to go into that rn
    global M 
    sum = 0
    for x in range(inputTuple[0]-M, inputTuple[0]+M):
        for y in range(inputTuple[1]-M, inputTuple[1]+M):
             for z in range(inputTuple[2]-M, inputTuple[2]+M):
                 try:
                     w = weight(inputTuple,(x,y,z),data)
                     sum = sum + w*data[x,y,z] 
                 except IndexError:
                     sum = sum + 0
                     continue
    
    #print("sum",sum )
    return sum





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



outputImage = np.zeros(shape=(75,75,75))
#going through each value and caluclulating their new thing

noisyImage = noisy(data)


##GLobal M to define the neighborhood
M = 2
h = 10



##precalculating the neighboorhood values 

            
blockMeansVal = np.zeros(shape=(75,75,75))
blockVectors = np.zeros(shape=(75,75,75))
    

##I think you need to make the blocks first and then knowing the increment you find the 
## array and get the value 
def blockValues(outputImage, noisyImage):
   
    count = 0
    for x in range(0,75):
        for y in range(0,75):
            for z in range(0,75):
                #need to get mean of its self
                sum = np.sum(noisyImage[x:x+5,y:y+5,z:z+5])
                mean = sum/125
                global blockMeansVal
                blockMeansVal[x,y,z] = mean
                
                count = count + 1
    
   


##now using the blockMeans array we do nlm for it

def blockNLM(outputImage, noisyImage):
   
    count = 0
    for x in range(0,75):
        for y in range(0,75):
            for z in range(0,75):
                global blockVectors
                global blockMeansVal
                blockVectors[x,y,z] = getNewValue((x,y,z),blockMeansVal)
              




##now after the block means have been calcaled to find the voxel value we average all the block vectors 
##in its area


""" def loopOperation(outputImage):
   
     for x in range(0,15):
        for y in range(0,15):
            for z in range(0,15):
                for()

                global blockVectors
                global blockMeansVal
                blockVectors[x,y,z] = getNewValue((x,y,z),blockMeansVal) """
               



""" 
items = [((50,52,),outputImage),((53,54),outputImage),((55,56),outputImage),((57,58),outputImage),((59,60),outputImage),((61,62),outputImage),((63,64),outputImage),((65,66),outputImage),((67,68),outputImage),((69,70),outputImage),((71,72),outputImage),((73,75),outputImage)]
executor = concurrent.futures.ProcessPoolExecutor(12)
futures = [executor.submit(loopOperation, item) for item in items]
concurrent.futures.wait(futures) """


#
#loopOperation((0,75),outputImage, noisyImage)

blockValues(outputImage,noisyImage)

blockNLM(outputImage,noisyImage)




""" plt.subplot(1,3,1)
plt.imshow(data[0:75,0:75,60], interpolation = 'nearest')
plt.title("ground")



plt.subplot(1,3,2)
plt.imshow(noisyImage[0:75,0:75,60], interpolation = 'nearest')
plt.title("noisy")



plt.subplot(1,3,3)
plt.imshow(outputImage[:,:,10], interpolation = 'nearest')
plt.title("filtered")
plt.show()
 """
