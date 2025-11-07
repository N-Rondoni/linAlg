import numpy as np
import matplotlib.pyplot as plt
import os

errors_Jac = np.loadtxt('iter_error_Jac_D1000.dat')
errors_Sie = np.loadtxt('iter_error_Sie_D1000.dat')

xj = errors_Jac[:, 0]
yj = errors_Jac[:, 1]
yj = np.log(yj)


xs = errors_Sie[:, 0]
ys = errors_Sie[:, 1]
ys = np.log(ys)


if __name__=='__main__':
    plt.plot(xj, yj, '-r', label='Solutions found via Jacobi')
    plt.plot(xs, ys, '--b', label='Solutions found via Seidel')
    plt.legend(loc='upper right')  
    plt.title('Log of Error as a function of iteration, requiring error < 1E-10, D=1000')
    plt.xlabel('Iterations')
    plt.ylabel('Log of Error')
    plt.savefig('error_plot_JS.png')
    os.system('cp error_plot_JS.png /mnt/c/Users/nicho/Pictures') #saves picture to your pictures
    plt.show()



