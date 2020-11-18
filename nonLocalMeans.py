
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 
import math as math

img = nib.load('NormalBrains/t1_icbm_normal_1mm_pn3_rf20.mnc')
data = img.get_fdata()
data_small = data#np.asarray(data[:3,:3,:3])



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
                 total = total = data_small[index8,index1,index4]  + data_small[index7,index2,index4] + data_small[index7,index1,index5]
              else:   
                  total = data_small[index8,index1,index4] + data_small[index9,index1,index4] + data_small[index7,index2,index4]+ data_small[index7,index3,index4] + data_small[index7,index1,index5]+ data_small[index7,index1,index6]
                       
          except IndexError: 
              total = total +0
              
              
                
                   
          epsilon[x,row,column] = math.pow( math.sqrt(6/7) * ( data_small[x,row,column] - ( (1/6)*total ) ),2 )
    
    
sigma = ( 1/(data_small.size) ) * np.sum(epsilon)


