import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib 

img = nib.load('NormalBrains/t1_icbm_normal_1mm_pn3_rf20.mnc')
data = img.get_fdata()
print(data.mean())

plt.imshow(data[:,:,1], interpolation = 'nearest')
plt.show()
