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
                 total=0
                 
                 try:      
                      if index6<0 or index3<0 :
                          total = data_small[x,index1,index5]  + data_small[x,index2,index4]   + data_small[x,index2,index5]
                      else:   
                          total = data_small[x,index1,index5] + data_small[x,index1,index6] + data_small[x,index2,index4] + data_small[x,index3,index4]+ data_small[x,index3,index5] + data_small[x,index2,index6] + data_small[x,index3,index6] + data_small[x,index2,index5]
                           
                 except IndexError: 
                      total = total +0
        
                 epsilon[x,row,column] = math.sqrt(6/7) * ( data_small[x,row,column] - ( (1/6)*total ) )
    
    return ( 1/(data_small.size) ) * math.pow(np.sum(epsilon),2) 