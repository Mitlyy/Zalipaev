import numpy as np
from matplotlib.pyplot import *
from matplotlib import pyplot as plt

global c0, c1
c0 = 1
c1 = 0

global a1,a2
a1 = 3
a2 = 2
def func(x):
    return np.cos(x) - c1*a1-c1*x*a2-c0*a2
def K(s,t):
    return a1+a2*(s-t)
def volterra_it(n,x,y):
    p=0
    phi=0*np.array(x)
    for i in range(n):
        k=0
        for j in range(i):
            k+=K(x[i],y[j])*func(y[i])
        phi[i]=-k*(x[-1]-x[1])/len(x)+func(x[i])
    while p<50:
        for i in range(n):
            k=0
            for j in range(i):
                k+=K(x[i],y[j])*phi[j]
            k*= (x[-1]-x[1])/len(x)
            phi[i]=-k+func(x[i])
        p+=1
    return phi

def volterra_test(n):
    a = 0
    b = 10
    x = np.linspace(a,b,n)
    y = np.linspace(a,b,n)
    u = 0*np.array(x)
    phi = volterra_it(n,x,y)
    for i in range(n):
        k=0
        for j in range(i):
            k += phi[j]*(x[i]-y[j])
        u[i] = c0+c1*x[i]+k*(x[-1]-x[1])/len(x)

    tochn = 0*np.array(x)
    for i in range(n):
            tochn[i] = -0.6*np.exp(-2*x[i])+1.5*np.exp(-x[i])+3*np.sin(x[i])/10+np.cos(x[i])/10
    fig, ax = plt.subplots()
    ax.plot(x,u,'--', label = 'Численное решение',color = 'green')
    ax.plot(x,tochn,label = 'Точное решение', color = 'brown') #точное реш
    ax.legend()
    ax.grid()
    plt.show()

volterra_test(100)