from bayes_opt import BayesianOptimization
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 
import math
import time
from sklearn.metrics import mean_squared_error
from multiprocessing import Pool

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

    
    
def weight(tupleI, tupleJ, data,Z,b):
    
    
    global s
   
    N = 27    
    output = 0
    
    
    uNJ = np.asarray ( neighboorhoodIntensities(data,1,tupleJ) )
    uNI = np.asarray ( neighboorhoodIntensities(data,1,tupleI) )

    dist = np.linalg.norm(uNJ-uNI)
    div = dist/(2*b*s*N)
    exponential = math.exp((-1)*div)
    output = (1/Z)*exponential
    return output

def getNewValue(intuple, indata,Z,b,M):
 
    total = 0;
    suma = 0;
     
    
    #this defines the search volume
    for x in range(intuple[0]-M,intuple[0]+M+1):
      for y in range(intuple[1]-M,intuple[1]+M+1):
        for z in range(intuple[2]-M,intuple[2]+M+1):          
                  w = weight(intuple,(x,y,z),indata,Z,b)
                  total = total + w*indata[x,y,z]
                  suma = suma + w
           
                    
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


img = nib.load('NormalBrains/t1_icbm_normal_1mm_pn3_rf20.mnc')


data_big = img.get_fdata() #input data
data = data_big[:3,:25,:25]
noisyImage = noisy(data) # added gaussian noise to original data
s = sigma(noisyImage)



def black_box_function(b, M):
    
  
    M = int(M)
    Z = math.pow(((2*M)+1),3)
    output = np.zeros(data.shape) #output image
    padded = padding(noisyImage,M)
    
    
    for x in range(M,padded.shape[0]-M):
      for y in range(M,padded.shape[1]-M):
         for z in range(M,padded.shape[2]-M):
            output[x-M,y-M,z-M] = getNewValue((x,y,z),padded,Z,b,M)
           
            #print(x)
    
    return -psnr(data,output)

pbounds = {'b': (0, 1), 'M': (1, 3)}

optimizer = BayesianOptimization(
    f=black_box_function,
    pbounds=pbounds,
    random_state=1,
)

optimizer.maximize(
    init_points=10,
    n_iter=50,
)




   