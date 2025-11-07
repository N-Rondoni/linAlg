import numpy as np
import matplotlib.pyplot as plt
import os

errors = np.loadtxt('errors.dat')
x = errors[:,0]
y = errors[:,1]


plt.plot(x, y)
plt.title('Error as a function of number of singular values used in reconstruction')
plt.xlabel('Number of Singular Values')
plt.ylabel('Error')
plt.savefig('error_plot.png')
os.system('cp error_plot.png /mnt/c/Users/nicho/Pictures') #saves picture to your pictures
plt.show()



