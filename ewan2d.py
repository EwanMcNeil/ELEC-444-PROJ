import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 
import math
import time
from random import seed
from random import randint
from sklearn.metrics import mean_squared_error

#The classical definition of the NL-means filter considers that each voxel can be linked to all the others, but#
# for practical computational reasons the number of voxels taken into account in the weighted
#average can be limited to the so-called “search volume” Vi of size (2M+1)3, centered at the current voxel xi.
def sigma(data):
    data_small = data
    epsilon = np.zeros(data_small.shape)
    
    for x in range(data_small.shape[0]):
         for row in range(data_small.shape[1]):
             for column in range(data_small.shape[2]):
                 index1=row
                 index2=row+1
                 index3=row-1
                 
                 index4=column
                 index5=column+1
                 index6=column-1 
                 
                 index7=x
                 index8=x+1
                 index9=x-1
                 
                 total=0
                 
                 try:      
                      if index6<0 or index3<0 or index9<0 :
                          total = data_small[index8,index1,index4]  + data_small[index7,index2,index4] + data_small[index7,index1,index5]
                      else:   
                          total = data_small[index8,index1,index4] + data_small[index9,index1,index4] + data_small[index7,index2,index4]+ data_small[index7,index3,index4] + data_small[index7,index1,index5]+ data_small[index7,index1,index6]
                           
                 except IndexError: 
                      total = total +0
        
                 epsilon[x,row,column] = math.pow( math.sqrt(6/7) * ( data_small[x,row,column] - ( (1/6)*total ) ),2 )
                 
    return ( 1/(data_small.size) ) * np.sum(epsilon)

def psnr(truth,denoised):
    
    truth_r = np.reshape(truth,(1,truth.size))
    denoised_r = np.reshape(denoised,(1,denoised.size))
    
    rms = math.sqrt(mean_squared_error(truth_r, denoised_r))
    psnr = 20* (math.log10(255/rms))
    
    return psnr 
    
 

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
    s = sigma(data)
    b = 1 
    N = 27


    output = 0

    uNJ = neighboorhoodIntensities(tupleJ, data)
    uNI = neighboorhoodIntensities(tupleI, data)
    
    dist = np.linalg.norm(uNJ-uNI)
    div = dist/(2*b*s*N)
    exponential = math.exp((-1)*div)
    output = (1/Z)*exponential
   
    return output

def getNewValue(inputTuple, data):

  
    sum = 0;
    global M 
    for x in range(inputTuple[0]-M, inputTuple[0]+M):
        for y in range(inputTuple[1]-M, inputTuple[1]+M):
             for z in range(inputTuple[2]-M, inputTuple[2]+M):
                 try:
                        w = weight(inputTuple,(x,y,z),data)
                        sum = sum + w*data[x,y,z] 
                 except IndexError:
                     sum = sum + 0
                     continue
    return sum

def noisy(image):
      row,col,ch= image.shape
      mean = 273
      var = 3000
      sigma = var**0.5
      gauss = np.random.normal(mean,sigma,(row,col,ch))
      gauss = gauss.reshape(row,col,ch)
      noisy = image + gauss
      return noisy


img = nib.load('NormalBrains/t1_icbm_normal_1mm_pn3_rf20.mnc')


data_big = img.get_fdata() #input data
data = data_big[:6,:6,:6]
print(data.mean())

M = 3


output = np.zeros(data.shape) #output image
noisyImage = noisy(data) # added gaussian noise to original data

start_time = time.time()


for x in range(data.shape[0]):
  for y in range(data.shape[1]):
     for z in range(data.shape[2]):
        output[x,y,z] = getNewValue((x,y,z),noisyImage)
        print(x)

print("--- %s seconds ---" % (time.time() - start_time))     



plt.subplot(1,3,1)
plt.imshow(data[1,:,:], interpolation = 'nearest')
plt.title("ground")



plt.subplot(1,3,2)
plt.imshow(noisyImage[1,:,:], interpolation = 'nearest')
plt.title("noisy")


plt.subplot(1,3,3)
plt.imshow(output[1,:,:], interpolation = 'nearest')
plt.title("filtered")
plt.show()
