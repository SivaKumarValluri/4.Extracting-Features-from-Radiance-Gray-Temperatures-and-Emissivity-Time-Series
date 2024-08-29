"""
Created on Fri Mar 18 15:57:53 2022

@author: Siva Kumar Valluri
"""

import os
import glob
import io 
import numpy as np
import math
#from scipy.interpolate import interp1d
#from scipy.interpolate import CubicSpline
from scipy.interpolate import PchipInterpolator as pcip
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy import integrate
import pandas as pd
import statistics as stat

#functions used:---------------------------------------------------------------------------------------------------------------------------------------------   

def fit_derivative_and_integral(x,y):
    #f = interp1d(x, y, kind='cubic')
    #f = CubicSpline(x,y)
    f=pcip(x, y, axis=0, extrapolate=None)
    #f = UnivariateSpline(x,y)
    x_new= np.logspace(math.log10(x[0]), math.log10(x[-1]), num=500, endpoint=True,base=10.0)
    logx_new=[math.log10(x_new[i]) for i in range(len(x_new))]
    y_new=f(x_new)
    dydlogx=np.gradient(y_new,logx_new)
    dydx=np.gradient(y_new,x_new)
    cd=integrate.cumtrapz(y_new, x_new)
    cdf=np.array([0])
    cdf=np.append(cdf,cd)
    cd2=integrate.cumtrapz(y_new,logx_new)
    cdf2=np.array([0])
    cdf2=np.append(cdf2,cd2)       
    return [x_new,y_new,dydlogx,dydx,cdf,cdf2]
    

def peak_and_saddlefinder(x,y,dydlogx,x2,y2,z2): #x-fitted time, y-fitted Radiance, dRdlogt, unfitted t_2, unfitted T, unfitted T_error
    pks, _ = find_peaks(y,height=0)
    #extreme points
    xx=[x[pks[i]] for i in range(len(pks))]
    yy=[y[pks[i]] for i in range(len(pks))]
    cutofftime1=100 #use 100 if time scale in ns instead of seconds
    cutofftime2=1000 #use 1000 if time scale in ns instead of seconds
    
    #Peak-1 identification and temperature at peak 1
    pkp1=[x for x in xx if x<cutofftime1] 
    indexp=[xx.index(x) for x in pkp1]
    pkv1=[yy[x] for x in indexp]
    bb=[pkv1[x]/max(pkv1) for x in range(len(pkv1))]
    Peakvalue=yy[bb.index(next(x for x in bb if x > 0.5))] #0.5 value is to ensure that first 'signifcant' peak over 0.5Imax is chosen
    Peakposition=xx[yy.index(Peakvalue)]
    diff=np.absolute(x2-Peakposition)
    index_T=diff.argmin()
    a=y2[index_T-2:index_T+2]
    a_error=z2[index_T-2:index_T+2]
    PeakTemp=np.nanmean(a)
    error1=np.nanmean(a_error)
    
    #Saddle point identification and temperature at saddle point
    new=dydlogx-min(dydlogx)+10
    log_new=[math.log10(new[x]) for x in range(len(new))]
    av=np.mean(log_new[0:40])
    d_a1=np.absolute(x-cutofftime1)
    i1=d_a1.argmin() #index for where x (time) is closest to 100ns
    d_a2=np.absolute(x-cutofftime2)
    i2=d_a2.argmin() #index for where x (time) is closest to 1000ns      
    difference_array = np.absolute(log_new-av)
    indexs = difference_array[i1:i2].argmin()+i1        
    Saddlepoint=x[indexs]
    Saddlevalue=y[indexs]
    diff_2=np.absolute(x2-Saddlepoint)
    index_T2=diff_2.argmin()
    a1=y2[index_T2-2:index_T2+2]
    a1_error=z2[index_T2-2:index_T2+2]
    SaddleTemp=np.nanmean(a1)
    error2=np.nanmean(a1_error)
    BulkTemp=y2[index_T2:-1]
    MaxBulkTemp=max(BulkTemp)
    MedianBulkTemp=stat.median(BulkTemp)
    FinalBulkTemp=np.nanmean(BulkTemp[-5:-1])
          
    return [xx,yy,Peakvalue,Peakposition,PeakTemp,error1,Saddlepoint,Saddlevalue,SaddleTemp,error2,indexs,index_T,index_T2,MaxBulkTemp,MedianBulkTemp,FinalBulkTemp]


def peak_finder(x,dydlogx,cutofftime1,cutofftime2):
    pks, _ = find_peaks(dydlogx,height=0)
    #extreme points
    xx=[x[pks[i]] for i in range(len(pks))]
    yy=[dydlogx[pks[i]] for i in range(len(pks))]
    
    #Peak-1 identification
    pkp1=[x for x in xx if x<cutofftime1] 
    indexp=[xx.index(x) for x in pkp1]
    if len(indexp)==0:
        Peakvalue=math.nan
        Peakposition=math.nan
    else:
        Peakvalue=max([yy[x] for x in indexp])
        Peakposition=xx[yy.index(Peakvalue)]
    #Peak-2 identification
    pkp2=[x for x in xx if x>cutofftime2] 
    indexp2=[xx.index(x) for x in pkp2]
    Peakvalue2=max([yy[x] for x in indexp2])
    Peakposition2=xx[yy.index(Peakvalue2)]    
    return [Peakposition,Peakvalue,Peakposition2,Peakvalue2]

def plotmepls(xr,yr,x2,y2,x_ex,y_ex,px,py,sx,sy,x):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('time,s')
            ax1.title.set_text('Radiance and dR/dlogt plots along with extrema in run '+str((x/2)+1))
            plt.xscale('log')
            ax1.set_ylabel('Radiance, W/Sr.m^2', color = 'black') 
            ax1.plot(xr, yr, color = 'black')
            plt.scatter(x_ex,y_ex,s=15,marker = 'o',color = 'black')
            plt.scatter(px,py,s=40, marker = 'o',color = 'r')
            plt.scatter(sx,sy,s=40, marker = 'o',color = 'r')         
            ax2 = ax1.twinx()
            ax2.set_ylabel('dR/dlogt, arb', color = 'blue') 
            ax2.plot(x2, y2, color = 'blue')
            ax2.tick_params(axis ='y', labelcolor = 'blue')
            return plt.show()

def plotme2pls(x_r,y_r,x_t,y_t,x_ex,y_ex,x_s,y_s,px,py,pt1_x, pt1_y,pt2_x, pt2_y,x):
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('time,s')
            ax1.title.set_text('Radiance and Temperature plots along with extrema in run '+str((x/2)+1))
            plt.xscale('log')
            ax1.set_ylabel('Radiance, W/Sr.m^2', color = 'black') 
            ax1.plot(x_r, y_r, color = 'black') 
            ax1.tick_params(axis ='y', labelcolor = 'black') 
            plt.scatter(x_s,y_s,marker = '^',color = 'blue')
            plt.scatter(x_ex,y_ex,marker = 'o',color = 'black')
            plt.scatter(px,py,marker = '^',color = 'blue')         
            ax2 = ax1.twinx()            
            ax2.set_ylabel('Temperature, K', color = 'red') 
            plt.scatter(x_t, y_t, s=25,color = 'red')
            ax2.plot(x_t, y_t, color = 'red')
            plt.scatter(pt1_x, pt1_y, s=100,marker = '*',color = 'blue')
            plt.scatter(pt2_x, pt2_y, s=100,marker = '*',color = 'blue') 
            ax2.tick_params(axis ='y', labelcolor = 'red')          
            return plt.show()
#Body-------------------------------------------------------------------------------------------------------------------------------------------------------

# Change location as needed and use double back slashes
os.chdir('D:\\Postdoctoral Work\\Experimental Data\\Project III- Study of porous Al-MoO3-KNO3\\Phi-3 (2.83)\\Reference-HMX clusters- Processed 3.5 kms low conc plate\\PMT Data')
path = os.getcwd()
txt_files = glob.glob(os.path.join(path, "*radiance.txt"))
txt_files_2 = glob.glob(os.path.join(path, "*grayTemp.txt"))




folder_name=path.rpartition('\\')[2] #folder name is assumed sample/condition detail


Excelwriter = pd.ExcelWriter(str(folder_name)+'.xlsx', engine='xlsxwriter')
Dlist=[]
namecounter=[]
for t in range(len(txt_files)):
       
    #reading txt file: Radiance data 
    temps=[]
    with io.open(txt_files[t], mode="r") as f:    
        next(f) #label
        next(f) #units
        next(f) #file name
        #copying data
        for line in f:
            temps.append(line.split())
    
    temp=np.array(temps, dtype=np.float32) #actual temporary file 
    temps=[]
    
    #reading txt file: corresponding temperature data  
    temps_2=[]
    with io.open(txt_files_2[t], mode="r") as f:    
        next(f) #label
        next(f) #units
        next(f) #file name
        #copying data
        for line in f:
            temps_2.append(line.split())
    
    temp_2=np.array(temps_2, dtype=np.float32) #actual temporary file 
    temps_2=[]
    
    file_analyzed = txt_files[t].split('\\')[-1]
    file_analyzed = file_analyzed.split('.txt')[0]
    
    #pocessing each sample: has several runs in txt file
    Dmain = pd.DataFrame(columns = ['Peak-I-x','Peak-I-y','Peak T/K','Error','transition-x','transition-y','Transition T','Error2','Area-Shock rise','Area-decay','Area-Growth','NArea-Shock rise','NArea-decay','NArea-Growth','Peak-I-MaxRate-x','Peak-I-Maxrate-y','Peak-II-Maxrate-x','Peak-II-Maxrate-y','MaxBulkT/K','MedianBulkT/K','FinalBulkT/K'])
    for i in range(0,np.shape(temp)[1],2):
        
        #Radiance 
        t_1=[temp[j,i] for j in range(len(temp[:,i])) if temp[j,i]>0] #accepting positive values for time
        R=[temp[j,i+1] for j in range(0,len(t_1),1)] #accepting non-NAN values for radiance 
        
        #Temperature
        corr=int(i*(3/2))       
        t_2=[temp_2[j,corr] for j in range(len(temp_2[:,corr])) if temp_2[j,corr]>0] #accepting positive values for time
        T=[temp_2[j,corr+1] for j in range(len(t_2))] #accepting all Temperature values
        T_error=[temp_2[j,corr+2] for j in range(len(t_2))] #accepting all Temperature error values
        
               
        #Fitting radiance data and arriving at instantaneous differential and cumulative integral values
        f_and_d=fit_derivative_and_integral(t_1,R) #tuple with new fit x and y and dydx
        t_fit=f_and_d[0]   #Time in nano seconds
        R_fit=f_and_d[1]   #Radiance W/Sr.m^2
        dRdlogt=f_and_d[2] #Derivative of Radiance on log of time 
        dRdt=f_and_d[3]    #Derivative of Radiance as a function of time
        Rcdf=f_and_d[4]    #Cumulative Integral on linear timescale
        Rcdflogt=f_and_d[5]   #Cumulative Integral on logscale
        
        
        #Extrema location in Radiance, differnatial radiance and corresponding Temperature at extrema identifed 
        p_and_s=peak_and_saddlefinder(t_fit,R_fit,dRdlogt,t_2,T,T_error)
        t_extreme=p_and_s[0]
        R_extreme=p_and_s[1]
        Peakvalue=p_and_s[2]
        Peakposition=p_and_s[3]
        PeakTemp=p_and_s[4]
        error=p_and_s[5]
        Saddlepoint=p_and_s[6]
        Saddlevalue=p_and_s[7]
        SaddleTemp=p_and_s[8]
        error2=p_and_s[9]
        index=p_and_s[10] #index of saddle point in t_fit data
        index_T=int(p_and_s[11])#index of Peak1 in t_2 data
        index_T2=int(p_and_s[12]) #index of Saddle point in t_2 data
        MaxBulkTemp=p_and_s[13]
        MedianBulkTemp=p_and_s[14]
        FinalBulkTemp=p_and_s[15]

        #Max dRdlogt location and values for nano and micro second peaks
        p2p=peak_finder(t_fit,dRdlogt,Peakposition,Saddlepoint)
        Peak1maxdiffposition=p2p[0]
        Peak1maxdiffvalue=p2p[1]
        Peak2maxdiffposition=p2p[2]
        Peak2maxdiffvalue=p2p[3]

        #Area under curve in three sections identified
        A_shockrise=Rcdf[np.where(t_fit == Peakposition)]
        A_decay=Rcdf[index]-A_shockrise
        A_bulk=Rcdf[-1]-A_decay
        
        #Normalized Area: Integral of Radiance on logtime 
        A_s=Rcdflogt[np.where(t_fit == Peakposition)]
        A_d=Rcdflogt[index]-A_s
        A_b=Rcdflogt[-1]-A_d
               
        #Writing data 
        AllXY = np.column_stack((Peakposition,Peakvalue,PeakTemp,error,Saddlepoint,Saddlevalue,SaddleTemp,error2, A_shockrise,A_decay,A_bulk,A_s,A_d,A_b,Peak1maxdiffposition,Peak1maxdiffvalue,Peak2maxdiffposition,Peak2maxdiffvalue,MaxBulkTemp,MedianBulkTemp,FinalBulkTemp))
        X = pd.DataFrame(AllXY,columns = ['Peak-I-x','Peak-I-y','Peak T/K','Error','transition-x','transition-y','Transition T','Error2','Area-Shock rise','Area-decay','Area-Growth','NArea-Shock rise','NArea-decay','NArea-Growth','Peak-I-MaxRate-x','Peak-I-Maxrate-y','Peak-II-Maxrate-x','Peak-II-Maxrate-y','MaxBulkT/K','MedianBulkT/K','FinalBulkT/K'])
        Dmain = Dmain.append(X)
        

        #print("done") 
        
        #Plots
        #plotmepls(t_fit,R_fit,t_fit,dRdlogt,t_extreme,R_extreme,Peakposition,Peakvalue,Saddlepoint,Saddlevalue, i)
        #plotme2pls(t_fit,R_fit,t_2,T,t_extreme,R_extreme,Saddlepoint,Saddlevalue,Peakposition,Peakvalue,t_2[index_T],PeakTemp,t_2[index_T2],SaddleTemp, i)
    
    #fname  = "%s.csv" %(folder_name)
    #Dmain.to_csv(fname)
    Dlist.append(Dmain)
    namecounter.append(file_analyzed)
    print("txt file "+str(t+1)+" complete")

for i, file in enumerate (Dlist):
    file.to_excel(Excelwriter, sheet_name=str(namecounter[i]),index=False)    

Excelwriter.save()
Excelwriter.close()