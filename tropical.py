import matplotlib.pyplot as plt
import numpy as np
import random


c = [random.random() * 3 for i in range(10)]
def ln_H(x,y):
    #return np.max([0 , x , y , 0.01-x-y])
    #data = [0 , x , y , 30-x-y]

    t = 3
    a = 0
    # data = np.array([0 , t+x , y   , t-x-4*y , a + t-2*y ])

    a = np.log(5)
    t = -np.log(2)
    x_sub = y
    y_sub = y + x
    data = np.array([2 * a + 4 * t + 4 * x_sub,
                     2 * a + 4 * t + 3 * x_sub, 
                     a + y_sub,
                     a + x_sub + 2 * t + y_sub,
                     a + 2 * t + 2* x_sub + y_sub,
                     2 * a + 3 * t + 2 * x_sub + y_sub,
                     2 * a + 5 * t + 3 * x_sub + y_sub,
                     3 * a + 6 * t + 4 * x_sub + y_sub,
                     2 * y_sub, 
                     a + 3 * t + x_sub + 2 * y_sub ])
    
    data2 = np.array([c[0] + 4 * x_sub,
                     c[1] + 3 * x_sub, 
                     c[2] + y_sub,
                     c[3] + x_sub + y_sub,
                     c[4] + 2 * x_sub + y_sub,
                     c[5] + 2 * x_sub + y_sub,
                     c[6] + 3 * x_sub + y_sub,
                     c[7] + 4 * x_sub + y_sub,
                     c[8] + 2 * y_sub, 
                     c[9] + x_sub + 2 * y_sub ])
    

    maximum = np.max(data)
    
    for i in range(len(data)):
        #if (abs(data[i] - maximum) < 0.0000001):

        if (abs(data[i] - maximum) < 0.000001):
            data[i] = maximum - 1000;
            break
    return np.min( np.abs(data - maximum))

print(c)
res = 701
array = np.ones((res,res))
bw_array = np.zeros((res,res))

for i in range(-res // 2, -res // 2 + res):
    for j in range(-res // 2, -res // 2 + res):
        val = ln_H(j/res * 40 , i/res * 20)
        array[i+res // 2 , j+res // 2] = val
        if val < 0.1:
            for k in [-1, 0, 1]:
                for l in [-1, 0, 1]:
                    bw_array[i+res // 2 + k, j+res // 2 + l] = 1

plt.imshow(bw_array, cmap='Greys' ,interpolation='nearest')
plt.show()
plt.savefig("graphics/tropical.pdf")

plt.imshow(array)
plt.show()
