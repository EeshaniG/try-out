# -*- coding: utf-8 -*-
"""
Created on Sun May 01 21:37:24 2016

@author: Eeshani Godbole
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.optimize as sco
import scipy.integrate as sci

filehandle=open('C:\\ExamProblemData2.1.xlsx','r')
print filehandle


df=pd.ExcelFile('C:\\ExamProblemData2.1.xlsx')
sheet=df.parse(0)


time= sheet.icol(0).real
time= np.reshape(time,(time.shape[0],))


conc_A= sheet.icol(1).real
conc_A= np.reshape(conc_A, (conc_A.shape[0],))

conc_B= sheet.icol(2).real
conc_B= np.reshape(conc_B,(conc_B.shape[0],))

conc= np.array([conc_A,conc_B]).T




def func(c,time,p):
    [k,m,n]=p
    dcdt= np.zeros_like (conc_init)
    dcdt [:]= k*c[0]**m*c[1]**n
    return dcdt

def error(p, conc_init,time):
    c= sci.odeint(func,conc_init,time, args=(p,))
    abs_err= conc - c
    err= abs_err[:,0]+abs_err[:,1]
    return err
   

conc_init= conc[0,:]
p0= [10,1,1]

res= sco.leastsq(error, p0, args=(conc_init,time,), full_output=1)

p= res[0]

print p
#this doesnt work, but the code is right
#pcov_diag= res[1].diagonal()
#perr= scipy.sqrt(pcov_diag)







