from scipy import signal, optimize, stats
import numpy as np
import sys
from scipy.special import gamma
import csv

def gamma_fit(X,o,p):
    k = p[0]
    theta = p[1]
    a = p[2]
    x_mod = X-o
    res = np.zeros(len(x_mod))
    if k>=1:
        nz = x_mod >= 0
    else:
        nz = x_mod > 0
    res[nz] = a * x_mod[nz]**(k-1) * np.exp(-x_mod[nz]/theta) / (theta **k * gamma(k))
    return res
        
def fit_optimize(X,y,boundaries):
    pranges = ((0.01,10),(0.01,150),(0.01,1))
    res_score = np.ones(boundaries[0]+1)*np.float('inf')
    res_param = [0 for i in range(boundaries[0]+1)]
    for i in range(15,boundaries[0]+1):
        f = lambda p: np.sum((gamma_fit(x,i,p) - y)**2)
        tmpres = optimize.brute(f, pranges,  full_output=True,finish=optimize.fmin)
        res_score[i] = tmpres[1]
        res_param[i] = tmpres[0]
    whichres = np.argmin(res_score)
    res = [res_param[whichres][0],res_param[whichres][1],res_param[whichres][2],whichres]
    return(res)

def read_y():
    x=[]
    with open('y_data.txt','r') as f:
    #next(f) # skip headings
        reader=csv.reader(f,delimiter='\t')
        for i in reader:
            x.extend(map(float, i))
    return(x)


if __name__ == "__main__":
    x= np.arange(35,115)
    test= [1.28774709e-05, 7.08260898e-05   , 3.86324126e-05   , 1.93162063e-05,
     2.57549417e-05  , 1.03019767e-04   , 8.37035606e-05   , 7.72648252e-05,
     9.65810315e-05  , 9.01422961e-05   , 1.28774709e-04   , 1.48090915e-04,
     1.99600798e-04  ,  1.99600798e-04   , 3.15498036e-04   , 2.76865624e-04,
     3.67007920e-04  ,  4.50711480e-04   , 5.34415041e-04   , 6.95383427e-04,
     7.14699633e-04  ,  7.34015839e-04   , 1.04307514e-03   , 1.23623720e-03,
     1.15897238e-03  ,  1.29418582e-03   , 1.72558110e-03   , 1.79640719e-03,
     1.86723328e-03  ,  1.95093684e-03   , 2.35013843e-03   , 2.60768785e-03,
     2.89099221e-03  ,  2.67851394e-03   , 3.23868392e-03   , 3.59281437e-03,
     3.84392505e-03  ,  4.29463653e-03   , 4.18517803e-03   , 4.18517803e-03,
     4.81617410e-03  ,  4.73247054e-03   , 5.29264053e-03   , 5.35702788e-03,
     5.60169983e-03  ,  5.90432039e-03   , 6.38722555e-03   , 6.07172751e-03,
     6.94095680e-03  ,  6.67696864e-03   , 6.63189750e-03   , 7.71360505e-03,
     6.67052991e-03  ,  7.37235207e-03   , 7.23713863e-03   , 7.93896079e-03,
     8.17075526e-03  ,  7.74579873e-03   , 8.42830468e-03   , 8.09992917e-03,
     8.18363273e-03  ,  8.91764857e-03   , 8.66009916e-03   , 8.93696478e-03,
     8.10636791e-03  ,  8.78243513e-03   , 9.26534029e-03   , 8.65366042e-03,
     9.25246282e-03  ,  9.15588178e-03   , 8.32528491e-03   , 8.75024145e-03,
     9.29109523e-03  ,  8.59571180e-03   , 9.56152212e-03   , 8.61502801e-03,
     8.47337583e-03  ,  9.36192132e-03   , 8.70517030e-03   , 6.43873543e-06]
     #out=fit_optimize(x,y,[35,115])
    y=read_y()
    out=fit_optimize(x,y,[35,115])
    print(out)
