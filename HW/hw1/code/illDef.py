import os
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    y = (x-2)**9
    return y

def g(x):
    y = x**9 - 18*x**8 + 144*x**7 - 672*x**6 + 2016*x**5 - 4032*x**4 + 5376*x**3 - 4608*x**2 + 2304*x - 512
    return y

    
if __name__ == '__main__':
    #set up domain
    xx = np.arange(1.920, 2.080, .001)
    #compute outputs, store in arrays
    yy_f = f(xx)
    yy_g = g(xx)
    #plot em
    plt.title('Plot of f(x)=(x-2)^9')
    plt.plot(xx, yy_f, linewidth = 1, label = 'f(x), not expanded')
    plt.plot(xx, yy_g, linewidth = 1, label = 'f(x), expanded')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()
    plt.savefig('ill_conditioned.png')
    #os.system('cp ill_conditioned.png /mnt/c/Users/nicho/Pictures/WSL_Pictures/AM_213a')

