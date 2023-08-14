#Libs
import georinex as gr
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from app.bin.allan_variance import allan_variance, params_from_avar
from app.bin.constants import *
import allantools
import base64
import io
import math
import numpy as np


# def

def avAnalysis(file,sats,freq):
    #__init__
    f1=Freq['G'][freq[1][0:2]]
    f2=Freq['G'][freq[3][0:2]]
    landa1=C/f1
    landa2=C/f2
    obsFile=gr.load(file,use=sats,meas=freq)

    #freq 1(C,L)
    meas1=obsFile[freq[0]].data
    meas2=obsFile[freq[1]].data*landa1

    #freq 2(C,L)
    meas3=obsFile[freq[2]].data
    meas4=obsFile[freq[3]].data*landa2

    #geometryFree
    GF=meas4-meas2
    deltaGF=GF[1:,:]-GF[0:-1,:]

    #code-phase combination
    CP1=meas1-meas2
    CP2=meas3-meas4
    deltaCP1=CP1[1:,:]-CP1[0:-1,:]
    deltaCP2=CP2[1:,:]-CP2[0:-1,:]

    #thershold & detect CS
    coef=[0.75,1]
    a0=(coef[0])*(landa1-landa2)
    a1=a0/2
    thershold=np.abs(a0-a1*(math.exp(-1*coef[1])))
    deltaGF[np.abs(deltaGF)>=thershold]=np.nan

    #cal noise 
    ParamsL1=[];ParamsL2=[];avL1=[];avL2=[]
    ParamsC1=[];ParamsC2=[];avC1=[];avC2=[]
    Sats=obsFile.sv.data
    for s in range(deltaGF.shape[1]):
        #phase
        buffer=deltaGF[:,s]
        buffer=buffer[~np.isnan(buffer)]
        #code[f1,f2]
        buffer1=CP1[:,s]
        buffer1=buffer1[~np.isnan(buffer1)]
        buffer2=CP2[:,s]
        buffer2=buffer2[~np.isnan(buffer2)]
        
        if(len(buffer)>10):
            t,v=allan_variance(buffer.copy())
            v_l1=v/(2*(1+(np.power(landa2,2)/np.power(landa1,2))))
            v_l2=v/(2*(1+(np.power(landa1,2)/np.power(landa2,2))))
            avL1.append(v_l1[0])
            avL2.append(v_l2[0])            
            params, av_pred = params_from_avar(t,v_l1)
            ParamsL1.append(params.copy())
            params, av_pred = params_from_avar(t,v_l2)
            ParamsL2.append(params.copy())
        else:
            avL1.append(np.nan);avL2.append(np.nan);ParamsL1.append([np.nan,np.nan,np.nan]);ParamsL2.append([np.nan,np.nan,np.nan])
        if(len(buffer1)>10):
            t1,v1=allan_variance(buffer1.copy())
            v1=v1/2
            avC1.append(v1[0])
            params, av_pred = params_from_avar(t1,v1)
            ParamsC1.append(params.copy())
        else:
            avC1.append(np.nan);ParamsC1.append([np.nan,np.nan,np.nan])
        if(len(buffer2)>10):
            t2,v2=allan_variance(buffer2.copy())
            v2=v2/2
            avC1.append(v2[0])
            params, av_pred = params_from_avar(t2,v2)
            ParamsC2.append(params.copy())
        else:
            avC2.append(np.nan);ParamsC2.append([np.nan,np.nan,np.nan])

    #
    data=io.BytesIO();data2=io.BytesIO();data3=io.BytesIO();data4=io.BytesIO()
    dataEncoded=[]
    # print('ParamsC1'+str(ParamsC1))
    #plots
    fig1,ax1 = plt.subplots(1,1)
    ax1.plot(Sats,[enum[0] for enum in ParamsC1],'r.')
    ax1.plot(Sats,[enum[1] for enum in ParamsC1],'b.')
    ax1.plot(Sats,[enum[2] for enum in ParamsC1],'g.')
    ax1.xaxis.set_tick_params(labelsize=8,rotation=90)
    plt.ylabel('noise parameters');plt.title('noise parameters(%s)'%(freq[0]))
    plt.legend(['quanzition noise','white noise','flicker noise'],fontsize=8,loc="right");fig1.savefig(data, format="png")
    dataEncoded1=base64.b64encode(data.getbuffer())
    dataEncoded.append(dataEncoded1)
    plt.close()

    fig1,ax1 = plt.subplots(1,1)
    ax1.plot(Sats,[enum[0] for enum in ParamsL1],'r.')
    ax1.plot(Sats,[enum[1] for enum in ParamsL1],'b.')
    ax1.plot(Sats,[enum[2] for enum in ParamsL1],'g.')
    ax1.xaxis.set_tick_params(labelsize=8,rotation=90)
    plt.ylabel('noise parameters');plt.title('noise parameters(%s)'%(freq[1]))
    plt.legend(['quanzition noise','white noise','flicker noise'],fontsize=8,loc="right");fig1.savefig(data2, format="png")
    dataEncoded1=base64.b64encode(data2.getbuffer())
    dataEncoded.append(dataEncoded1)
    plt.close()

    fig1,ax1 = plt.subplots(1,1)
    ax1.plot(Sats,[enum[0] for enum in ParamsC2],'r.')
    ax1.plot(Sats,[enum[1] for enum in ParamsC2],'b.')
    ax1.plot(Sats,[enum[2] for enum in ParamsC2],'g.')
    ax1.xaxis.set_tick_params(labelsize=8,rotation=90)
    plt.ylabel('noise parameters');plt.title('noise parameters(%s)'%(freq[2]))
    plt.legend(['quanzition noise','white noise','flicker noise'],fontsize=8,loc="right");fig1.savefig(data3, format="png")
    dataEncoded1=base64.b64encode(data3.getbuffer())
    dataEncoded.append(dataEncoded1)
    plt.close()

    fig1,ax1 = plt.subplots(1,1)
    ax1.plot(Sats,[enum[0] for enum in ParamsL2],'r.')
    ax1.plot(Sats,[enum[1] for enum in ParamsL2],'b.')
    ax1.plot(Sats,[enum[2] for enum in ParamsL2],'g.')
    ax1.xaxis.set_tick_params(labelsize=8,rotation=90)
    plt.ylabel('noise parameters');plt.title('noise parameters(%s)'%(freq[3]))
    plt.legend(['quanzition noise','white noise','flicker noise'],fontsize=8,loc="right");fig1.savefig(data4, format="png")
    dataEncoded1=base64.b64encode(data4.getbuffer())
    dataEncoded.append(dataEncoded1)
    plt.close()
    return dataEncoded


    