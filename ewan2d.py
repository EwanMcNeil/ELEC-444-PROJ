import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 
import math
import time
from sklearn.metrics import mean_squared_error
from multiprocessing import Pool


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
    
    rms = mean_squared_error(truth_r, denoised_r)
    #psnr = 20* (math.log10(255/rms))
    
    return rms 

    


def N_var_mean(N):

    return np.mean(N),np.var(N)
 
    
def neighboorhoodIntensities(a,d, tuplein):
    return [[[ a[x][y][z] if not(x==tuplein[0] and y==tuplein[1] and z==tuplein[2]) and x >= 0 and x < a.shape[0] and y >= 0 and y < a.shape[1] and z >= 0 and z < a.shape[2] else 0          
      for x in range(tuplein[0]-d, tuplein[0]+1+d)] 
          for y in range(tuplein[1]-d, tuplein[1]+1+d)]
              for z in range(tuplein[2]-d, tuplein[2]+1+d)] 
    
 
#weight with voxel selection
def weight_voxel_selection(tupleI, tupleJ, data,Z):
   
    global s
    global M
    b = 1 
    N = 27

    u1 = 0.95
    s1 = 0.5
    
    output = 0
    data_I = data[M:data.shape[0]-M , M:data.shape[1]-M,M:data.shape[2]-M ]
    uNJ = np.asarray ( neighboorhoodIntensities(data,1,tupleJ) )
    uNI = np.asarray ( neighboorhoodIntensities(data_I,1,tupleI) )
    
    uNI_mean_var=uNI
    uNJ_mean_var=uNJ
    
    uNI_mean_var[1,1,1] = data[tupleI]
    uNJ_mean_var[1,1,1] = data[tupleJ]
    
    mean_i,var_i = N_var_mean(uNI_mean_var)
    mean_j,var_j = N_var_mean(uNJ_mean_var)
    
    if( (mean_i==0 and var_i==0) or (mean_j==0 and var_j==0) ):
        return 0
    else:
      div_mean = mean_i/mean_j
      div_var = var_i/var_j
   

    if( ( (div_mean > u1) and (div_mean < 1/u1) ) and ( (div_var > s1) and (div_var < 1/s1) )       ):
        dist = np.linalg.norm(uNJ-uNI)
        div = dist/(2*b*s*N)
        exponential = math.exp((-1)*div)
        output = (1/Z)*exponential
        return output
    else:
        return 0
    
#weight without voxel selection    
def weight(tupleI, tupleJ, data,Z):
    
    
    global s
    b = 1 
    N = 27    
    output = 0
    
    
    uNJ = np.asarray ( neighboorhoodIntensities(data,1,tupleJ) )
    uNI = np.asarray ( neighboorhoodIntensities(data,1,tupleI) )

    dist = np.linalg.norm(uNJ-uNI)
    div = dist/(2*b*s*N)
    exponential = math.exp((-1)*div)
    output = (1/Z)*exponential
    return output

                

# select=1 : use voxel selection selection=0 : do not use voxel selection
def getNewValue(intuple, indata, select,Z):
 
    total = 0;
    global M 
    
    #this defines the search volume
    for x in range(intuple[0]-M,intuple[0]+M+1):
      for y in range(intuple[1]-M,intuple[1]+M+1):
        for z in range(intuple[2]-M,intuple[2]+M+1): 
                if(select==1):
                 w = weight_voxel_selection(intuple,(x,y,z),indata,Z)
                 total = total + w*indata[x,y,z]
                elif(select==0):
                  w = weight(intuple,(x,y,z),indata,Z)
                  total = total + w*indata[x,y,z]
           
                    
    return total

def noisy(image):
      row,col,ch= image.shape
      mean = 273
      var = 3000
      sigma = var**0.5
      np.random.seed(10)
      gauss = np.random.normal(mean,sigma,(row,col,ch))
      gauss = gauss.reshape(row,col,ch)
      noisy = image + gauss
      return noisy

def padding(image,d):
   padded = np.pad(image, ((d, d), (d, d), (d, d)),mode='constant', constant_values=0)
   return padded
    
def Z_calculator(intuple,indata):
    
     total = 0
     for x in range(intuple[0]-M,intuple[0]+M+1):
      for y in range(intuple[1]-M,intuple[1]+M+1):
        for z in range(intuple[2]-M,intuple[2]+M+1): 
               w = weight_z(intuple,(x,y,z),indata)
               total = total + w
      return total
  
def weight_z(tupleI, tupleJ, data):
    
    global s
    b = 1 
    N = 27    
    output = 0
    
    
    uNJ = np.asarray ( neighboorhoodIntensities(data,1,tupleJ) )
    uNI = np.asarray ( neighboorhoodIntensities(data,1,tupleI) )

    dist = np.linalg.norm(uNJ-uNI)
    div = dist/(2*b*s*N)
    exponential = math.exp((-1)*div)
    output = exponential
    return output      
    
    
img = nib.load('NormalBrains/t1_icbm_normal_1mm_pn3_rf20.mnc')


data_big = img.get_fdata() #input data
data = data_big[:3,:217,:181]
print(data.mean())




M = 1


output = np.zeros(data.shape) #output image
noisyImage = noisy(data) # added gaussian noise to original data
s = sigma(noisyImage)
padded = padding(noisyImage,M)

# Without using voxels
start_time_without_voxel = time.time()


for x in range(M,padded.shape[0]-M):
  for y in range(M,padded.shape[1]-M):
     for z in range(M,padded.shape[2]-M):
        Z = Z_calculator((x,y,z),padded)
        output[x-M,y-M,z-M] =  getNewValue((x,y,z),padded,0,Z)
       
        #print(x)

print("--- %s seconds (without voxel sel) ---" % (time.time() - start_time_without_voxel))  
print("PSNR(ground-noise)",psnr(data,noisyImage))
print("PSNR(ground-output)",psnr(data,output)) 

#plot results
plt.subplot(1,3,1)
plt.imshow(data[2,:,:], interpolation = 'nearest')
plt.title("ground")



plt.subplot(1,3,2)
plt.imshow(noisyImage[2,:,:], interpolation = 'nearest')
plt.title("noisy")


plt.subplot(1,3,3)
plt.imshow(output[2,:,:], interpolation = 'nearest')
plt.title("filtered")
plt.show()     


# Using voxels
start_time_with_voxel=time.time()

for x in range(data.shape[0]):
  for y in range(data.shape[1]):
     for z in range(data.shape[2]):
        Z = Z_calculator((x,y,z),padded)
        output[x,y,z] = getNewValue((x,y,z),padded,1,Z)
       
        #print(x)

print("--- %s seconds (with voxel sel) ---" % (time.time() - start_time_with_voxel))  
print("PSNR(ground-noise)",psnr(data,noisyImage))
print("PSNR(ground-output)",psnr(data,output))



#plot results
plt.subplot(1,3,1)
plt.imshow(data[2,:,:], interpolation = 'nearest')
plt.title("ground")



plt.subplot(1,3,2)
plt.imshow(noisyImage[2,:,:], interpolation = 'nearest')
plt.title("noisy")


plt.subplot(1,3,3)
plt.imshow(output[2,:,:], interpolation = 'nearest')
plt.title("filtered")
plt.show()
