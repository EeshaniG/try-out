#Experiment is single stage adsorption
#plotting values of log(x/m) vs logC
import numpy as np
import matplotlib.pyplot as plt
import scipy
data = np.genfromtxt('ads.txt', delimiter=',')#importing data from notepad

x = data[:,][:,0]
y =data[:,][:,1]

a= [[x, np.ones_like(x)]]
b=np.transpose(a)
result= np.lstsq(b, y)
slope= result(0)
intercept= result(1)
n=(1/(result(0)))
k=scipy.exp(result(1))

ynew=slope*x+intercept
plot (x, ynew,"r")
plot (x,y)