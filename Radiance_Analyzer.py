"""
Created on Sun Mar 13 21:12:16 2022
@author: sivakumarvalluri
The txt files in folder expected to be named after flyer velocity with one decimal eg:4.0
The folder name can be anything and all of that will be used as workbook name along with txt file read
For a folder named 'Al-CuO 30 wt.%'
workbook name will be: 'Al-CuO 30 wt.%-4.0km/s.xlsx'
Each sheet in workbook will yield summary of a run: time,radiance, interporlated data, peaks, saddle points, etc
"""

#Libraries used:--------------------------------------------------------------------------------------------

import os
import glob
import io 
import numpy as np
import math
#from scipy.interpolate import interp1d
#from scipy.interpolate import CubicSpline
from scipy.interpolate import PchipInterpolator as pcip
from scipy.signal import find_peaks
import xlwt
from xlwt import Workbook
import matplotlib.pyplot as plt
from scipy import integrate
import pandas as pd
#functions used:--------------------------------------------------------------------------------------------   

def fit_and_derivative(x,y):
    #f = interp1d(x, y, kind='cubic')
    #f = CubicSpline(x,y)
    f=pcip(x, y, axis=0, extrapolate=None)
    x_new= np.logspace(math.log10(x[0]), math.log10(x[-1]), num=500, endpoint=True,base=10.0)
    logx_new=[math.log10(x_new[i]) for i in range(len(x_new))]
    y_new=f(x_new)
    dydlogx=np.gradient(y_new,logx_new)
    dydx=np.gradient(y_new,x_new)       
    return [x_new,y_new,dydlogx,dydx]
    

def peak_and_saddlefinder(x,y):
    pks, _ = find_peaks(y,height=0)
    xx=[x[pks[i]] for i in range(len(pks))]
    yy=[y[pks[i]] for i in range(len(pks))]
    return [xx,yy]
    

def isMonotonic(A):

    return (all(A[i] <= A[i + 1] for i in range(len(A) - 1)) or
            all(A[i] >= A[i + 1] for i in range(len(A) - 1)))

def plotmepls(x1,y1,x2,y2,x3,y3,px,py,x):
            plt.title('Extreme points identified in run'+str((x/2)+1))
            plt.xlabel('time,ns')
            plt.ylabel('Radiance')
            plt.plot(x1,y1,color = 'k')
            plt.plot(x2,y2,color = 'b')
            plt.scatter(x3,y3,marker = 'o',color = 'r')
            plt.scatter(px,py,marker = 's',color = 'g')
            plt.xscale('log')
            #plt.semilogx(x1,y1,x2,y2,x3,y3)
            return plt.show()


#Body------------------------------------------------------------------------------------------------

# Change location as needed and use double back slashes
os.chdir('D:\\Postdoctoral Work\\Experimental Data\\Project III- Study of porous Al-MoO3-KNO3\\Phi-3 (2.83)\\Reference-HMX clusters- Processed 3.5 kms low conc plate\\PMT Data')
path = os.getcwd()
txt_files = glob.glob(os.path.join(path, "*radiance.txt"))

folder_name=path.rpartition('\\')[2] #folder name is assumed sample/condition detail



for t in range(len(txt_files)):
    

    
    #reading txt file: all shots  
    temps=[]
    with io.open(txt_files[t], mode="r") as f:    
        next(f) #label
        next(f) #units
        next(f) #file name
        #copying data
        for line in f:
            temps.append(line.split())
    
    temp=np.array(temps, dtype=np.float32) #actual temporay file 
    temps=[]
    
    #xlsx file named after sample/condition detail and velocity
    #velocity_tested=txt_files[t][-7:].rpartition('.')[0]
    velocity_tested = txt_files[t].split('\\')[-1]
    velocity_tested = velocity_tested.split('.txt')[0]
    #pocessing each run 
    Dmain = pd.DataFrame(columns = ['extreme_x','extreme_y','transition-x','transition-y','Area-Shock rise','Area-decay','Area-Growth']) #Dataframe of radii and respective
    for i in range(0,np.shape(temp)[1],2):
        x=[temp[j,i] for j in range(len(temp[:,i])) if temp[j,i]>0] #accepting non-zero values for time
        y=[temp[j,i+1] for j in range(len(x))] #accepting non-NAN values for radiance 

        
        f_and_d=fit_and_derivative(x,y) #tuple with new fit x and y and dydx
        x_fit=f_and_d[0] #Time in nano seconds
        y_fit=f_and_d[1] #Radiance W/Sr.m^2
        dydlogx=f_and_d[2] #Derivative of Radiance on log of time 
        dydx=f_and_d[3]    #Derivative of Radiance as a fucntion of time
        cd=integrate.cumtrapz(y_fit, x_fit)
        cdf=np.array([0])
        cdf=np.append(cdf,cd)
        
        p_and_s=peak_and_saddlefinder(x_fit,y_fit) #tuple with x and y values of extreme and saddle points
        x_extreme=p_and_s[0]
        y_extreme=p_and_s[1]
        
        
        
        #Peak-1 detail
        pkp1=[x for x in x_extreme if x<10**(-7)]
        index1=[x_extreme.index(x) for x in pkp1]
        pkv1=[y_extreme[x] for x in index1]
        bb=[pkv1[x]/max(pkv1) for x in range(len(pkv1))]
        Peakvalue1=y_extreme[bb.index(next(x for x in bb if x > 0.5))]
        Peakposition1=x_extreme[y_extreme.index(Peakvalue1)]
        A_shockrise=cdf[np.where(x_fit == Peakposition1)]
        #plots
        plotmepls(x,y,x_fit,y_fit,x_extreme,y_extreme,Peakposition1,Peakvalue1, i)
        #saddle point detail
        new=dydlogx-min(dydlogx)+10
        log_new=[math.log10(new[x]) for x in range(len(new))]
        av=np.mean(log_new[0:40])
        d_a1=np.absolute(x_fit-10**(-7))
        i1=d_a1.argmin() #index for where x_fit (time) is closest to 100ns
        d_a2=np.absolute(x_fit-10**(-6))
        i2=d_a2.argmin() #index for where x_fit (time) is closest to 1000ns
        
        difference_array = np.absolute(log_new-av)
        index = difference_array[i1:i2].argmin()+i1
        
        Saddlepoint=x_fit[index]
        Saddlevalue=y_fit[index]
        A_decay=cdf[index]-A_shockrise        
        A_bulk=cdf[-1]-A_decay
        
        
        AllXY = np.column_stack((Peakposition1,Peakvalue1,Saddlepoint,Saddlevalue,A_shockrise,A_decay,A_bulk)) #Matrix of size and location 
        X = pd.DataFrame(AllXY,columns = ['extreme_x','extreme_y','transition-x','transition-y','Area-Shock rise','Area-decay','Area-Growth'])
        Dmain = Dmain.append(X)
        #print("done") 
        
        
        
    
    fname  = "%s.csv" %(velocity_tested)
    Dmain.to_csv(fname)

    print("txt file" +str(t+1)+" complete")              
        
            






