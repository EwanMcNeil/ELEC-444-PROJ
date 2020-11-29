import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 
import math
import time
from sklearn.metrics import mean_squared_error
from multiprocessing import Pool, Process , Manager

# Sigma Function for general estimation of h-parameter
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

#PSNR to calculated difference between results
def psnr(truth,denoised):
    
    truth_r = np.reshape(truth,(1,truth.size))
    denoised_r = np.reshape(denoised,(1,denoised.size))
    
    rms = mean_squared_error(truth_r, denoised_r)
    psnr = 20* (math.log10(255/rms))
    
    return abs(psnr) 

#function that return mean and variance of a block
def N_var_mean(N):

    return np.mean(N),np.var(N)
 
#function that return the neighborhood intensities at touplein  
def neighboorhoodIntensities(a,d, tuplein):
    return [[[ a[x][y][z] if not(x==tuplein[0] and y==tuplein[1] and z==tuplein[2]) and x >= 0 and x < a.shape[0] and y >= 0 and y < a.shape[1] and z >= 0 and z < a.shape[2] else 0          
      for x in range(tuplein[0]-d, tuplein[0]+1+d)] 
          for y in range(tuplein[1]-d, tuplein[1]+1+d)]
              for z in range(tuplein[2]-d, tuplein[2]+1+d)] 
    
 
#weight with voxel selection
def weight_voxel_selection(tupleI, tupleJ, data,Zv):
   
    global s
    global M
    b = 1 
    N = 27

    u1 = 0.95
    s1 = 0.5
    
    output = 0

    uNJ = np.asarray ( neighboorhoodIntensities(data,1,tupleJ) )
    uNI = np.asarray ( neighboorhoodIntensities(data,1,tupleI) )
    
    uNI_mean_var=uNI
    uNJ_mean_var=uNJ
    
    uNI_mean_var[1,1,1] = data[tupleI]
    uNJ_mean_var[1,1,1] = data[tupleJ]
    
    mean_i,var_i = N_var_mean(uNI_mean_var)
    mean_j,var_j = N_var_mean(uNJ_mean_var)
    
    if( (mean_i==0 and var_i==0) or (mean_j==0 and var_j==0) ):
        div_mean = math.inf
        div_var = math.inf 
    else:
     div_mean = mean_i/mean_j
     div_var = var_i/var_j
   

    if( ( (div_mean > u1) and (div_mean < 1/u1) ) and ( (div_var > s1) and (div_var < 1/s1) ) ):
        dist = np.linalg.norm(uNJ-uNI)
        div = dist/(2*b*s*N)
        exponential = math.exp((-1)*div)
        output = (1/Zv)*exponential
        return output
    else:
        return 0


def Z_voxel_selection(tupleI, tupleJ, data):
   
    u1 = 0.95
    s1 = 0.5
    
    z = 0

    uNJ = np.asarray ( neighboorhoodIntensities(data,1,tupleJ) )
    uNI = np.asarray ( neighboorhoodIntensities(data,1,tupleI) )
    
    uNI_mean_var=uNI
    uNJ_mean_var=uNJ
    
    uNI_mean_var[1,1,1] = data[tupleI]
    uNJ_mean_var[1,1,1] = data[tupleJ]
    
    mean_i,var_i = N_var_mean(uNI_mean_var)
    mean_j,var_j = N_var_mean(uNJ_mean_var)
    
    if( (mean_i==0 and var_i==0) or (mean_j==0 and var_j==0) ):
        div_mean = math.inf
        div_var = math.inf 
    else:
     div_mean = mean_i/mean_j
     div_var = var_i/var_j
   

    if( ( (div_mean > u1) and (div_mean < 1/u1) ) and ( (div_var > s1) and (div_var < 1/s1) ) ):
        z = z +1
    else:
        z = z + 0

    return z

    
#weight classical NL-mean   
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

                

#Get new value for classical implementation
def getNewValue(intuple, indata,Z):
 
    total = 0;
    suma = 0;
    global M 
    
    #this defines the search volume
    for x in range(intuple[0]-M,intuple[0]+M+1):
      for y in range(intuple[1]-M,intuple[1]+M+1):
        for z in range(intuple[2]-M,intuple[2]+M+1):          
                  w = weight(intuple,(x,y,z),indata,Z)
                  total = total + w*indata[x,y,z]
                  suma = suma + w
           
                    
    return total

#Get new value for voxel selection implementation
def getNewValue_voxel(intuple, indata):
 
    total = 0;
    suma = 0;
    Z = 0
    global M 
    
    for x in range(intuple[0]-M,intuple[0]+M+1):
      for y in range(intuple[1]-M,intuple[1]+M+1):
        for z in range(intuple[2]-M,intuple[2]+M+1):
            Z = Z + Z_voxel_selection(intuple,(x,y,z),indata)
            
    #this defines the search volume
    for x in range(intuple[0]-M,intuple[0]+M+1):
      for y in range(intuple[1]-M,intuple[1]+M+1):
        for z in range(intuple[2]-M,intuple[2]+M+1):             
                 w = weight_voxel_selection(intuple,(x,y,z),indata,Z)
                 total = total + w*indata[x,y,z]
                 suma = suma + w
    
    return total
           


###BLOCKWISE FUNCTIONS 
def getNewValueBlock(intuple, indata,M,N):

    total = 0;
    suma = 0;
    count = 0
    #blockZ = (N*2)**3
    blockZ = 8
    #this defines the search and only compares on blocks with the step of N
    for x in range(intuple[0]-M,intuple[0]+M+1,N):
        for y in range(intuple[1]-M,intuple[1]+M+1,N):
            for z in range(intuple[2]-M,intuple[2]+M+1,N): 
                    count = count + 1         
                    w = weight(intuple,(x,y,z),indata,blockZ)
                    total = total + w*indata[x,y,z]
                    suma = suma + w         
    return total

def blockwise(M,padded,data):
    N = 2
    A = 1   
    output = np.zeros(data.shape)
    centers = []
    ##this gets all the center values 
    for x in range(M,padded.shape[0]-M,N):
        for y in range(M,padded.shape[1]-M,N):
            for z in range(M,padded.shape[2]-M,N):
                    output[x-M,y-M,z-M] = getNewValueBlock((x,y,z),padded,M,N)
                    centers.append((x-M,y-M,z-M))
                   
    finalOutput = np.zeros(output.shape) #output image

    for x in range(0,output.shape[0]):
        for y in range(0,output.shape[1]):
            for z in range(0,output.shape[2]):
                    if not (x,y,z) in centers:
                        estimator = output[x-A:x+A+1,y-A:y+A+1,z-A:z+A+1]
                        count = np.count_nonzero(estimator)
                        sumOut = np.sum(estimator)
                        div = 0
                        if(count != 0):
                            div = sumOut/count
                            finalOutput[x,y,z] =  div
                        else:
                            finalOutput[x,y,z] =  0


    addition = np.zeros(output.shape) #output image
    addition = np.add(output,finalOutput)
    return addition


def noisy(image):
      row,col,ch= image.shape
      mean = 129
      var = 500
      sigma = var**0.5
      np.random.seed(10)
      gauss = np.random.normal(mean,sigma,(row,col,ch))
      gauss = gauss.reshape(row,col,ch)
      noisy = image + gauss
      return noisy

def padding(image,d):
   padded = np.pad(image, ((d, d), (d, d), (d, d)),mode='constant', constant_values=0)
   return padded




def parallelize_classic(x):
    
    output_thread = np.zeros(data.shape)
    for y in range(M,padded.shape[1]-M):
        for z in range(M,padded.shape[2]-M):
            output_thread[x-M,y-M,z-M] = getNewValue((x,y,z),padded,Z)
    return output_thread  

def parallelize_voxel_selection(x):
    
    output_thread = np.zeros(data.shape)
    for y in range(M,padded.shape[1]-M):
        for z in range(M,padded.shape[2]-M):
            output_thread[x-M,y-M,z-M] = getNewValue_voxel((x,y,z),padded)
    return output_thread 

def parallelize_blockwise(x):
    
            
    finalOutput = np.zeros(output.shape) #output image  
    addition = np.zeros(data.shape)

    A = 1   
    
    print("running parallelize_blockwise")
    for y in range(0,output.shape[1]):
        for z in range(0,output.shape[2]):
            if not (x,y,z) in centers:
                estimator = output[x-A:x+A+1,y-A:y+A+1,z-A:z+A+1]
                count = np.count_nonzero(estimator)
                sumOut = np.sum(estimator)
                div = 0
                if(count != 0):
                    div = sumOut/count
                    finalOutput[x,y,z] =  div
                else:
                    finalOutput[x,y,z] =  0


    addition = np.zeros(output.shape) #output image
    addition = np.add(output,finalOutput)

    return addition

def centers_blockwise(x):
    
    N = 2

    ##this gets all the center values 
    for y in range(M,padded.shape[1]-M,N):
            for z in range(M,padded.shape[2]-M,N):
                    output[x-M,y-M,z-M] = getNewValueBlock((x,y,z),padded,M,N)
                    centers.append((x-M,y-M,z-M))
    return (centers, output)



#######################################MAIN BEGINNING####################################

if __name__ == '__main__':
  
    img = nib.load('NormalBrains/t1_icbm_normal_1mm_pn3_rf20.mnc')

    data_big = img.get_fdata() #input data
    data = data_big[:3,:10,:10]
    noisyImage = noisy(data) # added gaussian noise to original data
    s = sigma(noisyImage)
    
    
    #following 2 values valid only for classic NL-means
    M = 1
    Z = math.pow(((2*M)+1),3)
    
    padded = padding(noisyImage,M)
    
    #create the matrix for outputs
    output_classic = np.zeros(data.shape) 
    output_voxel = np.zeros(data.shape)
    output_blockWise = np.zeros(data.shape)
    
    
    output_classic_thread = np.zeros(data.shape)
    output_voxel_thread = np.zeros(data.shape)
    output_blockwise_thread = np.zeros(data.shape)
    
    print("PSNR(ground-noise)",psnr(data,noisyImage))
    
    #Multiprocessing for classical implementation

    start_time_classic_thread = time.time()
    
    p_classic = Pool(processes=3)
    res_classic = p_classic.map(parallelize_classic,range(M,padded.shape[0]-M))
    p_classic.close()
    p_classic.join()
    
    print("--- %s seconds (classical(thread)) ---" % (time.time() - start_time_classic_thread)) 
    
    for i in range(data.shape[0]):
       output_classic_thread[i]  = np.array(res_classic[i][i])
       
    
    print("PSNR(ground-classic)",psnr(data,output_classic_thread))   
       
       
    #Multiprocessing for voxel selection implementation 
    
    start_time_voxel_selection_thread = time.time()
    
    p_voxel_sel = Pool(processes=3)
    res_voxel_sel = p_voxel_sel.map(parallelize_voxel_selection,range(M,padded.shape[0]-M))
    p_voxel_sel.close()
    p_voxel_sel.join()
    
    print("--- %s seconds (voxel_selection(thread)) ---" % (time.time() - start_time_voxel_selection_thread)) 
    for i in range(data.shape[0]):
       output_voxel_thread[i]  = np.array(res_voxel_sel[i][i])
    
    print("PSNR(ground-classic with voxel selection)",psnr(data,output_voxel_thread)) 
       

    #Multiprocessing for blockwise implementation
    centers = []
    output = np.zeros(data.shape)
    
    start_time_blockwise_thread = time.time()
    
    p_centers = Pool(processes=3)
    test = p_centers.map(centers_blockwise, range(M,padded.shape[0]-M) )
    p_centers.close()
    p_centers.join()
    
    for i in range(len(test)):
       centers = centers + test[i][0]
       
    for i in range(len(test)):
       output = output + test[i][1]
       

    
    p_blockwise = Pool(processes=3)
    test2 = p_blockwise.map(parallelize_blockwise,range(output.shape[0]) )
    p_blockwise.close()
    p_blockwise.join()
    
    blockwiseOutput = np.zeros(output.shape)
    
    for i in range(len(test2)):
        blockwiseOutput[i] =  test2[i][i]
        
    
    print("--- %s seconds (voxel_selection(thread)) ---" % (time.time() - start_time_blockwise_thread))
    
    #for i in range(data.shape[0]):
   #    output_blockwise_thread[i]  = np.array(res_blockwise[i][i])
       
    print("PSNR(ground-blockwise)",psnr(data,blockwiseOutput))

    
# blockWiseOutput = blockwise(M,Z,padded,data)
    


#plot results
fig = plt.figure()

ax1 = fig.add_subplot(5,1,1)
ax1.imshow(data[2,:,:])
#ax1.title("Ground")
ax1.axis('off')

ax2 = fig.add_subplot(5,1,2)
ax2.imshow(noisyImage[2,:,:])
#ax2.title("Noisy")
ax2.axis('off')

ax2 = fig.add_subplot(5,1,3)
ax2.imshow(output_classic_thread[2,:,:])
#ax2.title("Classical")
ax2.axis('off')

ax3 = fig.add_subplot(5,1,4)
ax3.imshow(output_voxel_thread[2,:,:])
#ax3.title("Voxel Selection")
ax3.axis('off')

ax4 = fig.add_subplot(5,1,5)
ax4.imshow(blockwiseOutput[2,:,:])
#ax4.title("Blockwise")
ax4.axis('off')


fig.savefig('result.png', dpi=400)


      
