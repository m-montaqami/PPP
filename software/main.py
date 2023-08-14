from Atmosphere import atmosphere
from Libs import *
from DataEntry import dataEntry
from Control import control
from Output import *
from PreProcess import preprocess
from PPP import process

OBS,SP3,CTL1,CTL2,OBS_RATE,SP3_RATE,OBS_TIMESEC,SP3_TIMESEC,X,Y,Z,CLK,X_PROC,Y_PROC,Z_PROC,CLK_PROC,Vx_s,Vy_s,Vz_s,EL_PROC,AZ_PROC,OBS_TIMESEC_PROC,SNR,AC_SNR,GIM,VTEC,MF_Iono,MF_Trop,ZHD,L_Matrix,N,X_cap=Start()

def exe():

    dataEntry(OBS,SP3)
    for enum in range(len(OBS)):    
        preproc=preprocess(OBS[enum],SP3[enum],userSats[enum])
        preproc.TimeProc(OBS_TIMESEC,SP3_TIMESEC,OBS_RATE,SP3_RATE)
        preproc.sp3Extract(X,Y,Z,CLK)
        preproc.proc(X[enum],Y[enum],Z[enum],CLK[enum],OBS_TIMESEC[enum],SP3_TIMESEC[enum],OBS_RATE[enum],SP3_RATE[enum],X_PROC,Y_PROC,Z_PROC,CLK_PROC,Vx_s,Vy_s,Vz_s,OBS_TIMESEC_PROC)

        ctl=control(userSats[enum],OBS[enum],SP3[enum],SP3_TIMESEC[enum],OBS_TIMESEC[enum],SP3_RATE[0],OBS_TIMESEC_PROC[enum])        
        ctl.clt1(CTL1,EL_PROC,AZ_PROC,SNR,AC_SNR,X_PROC[enum],Y_PROC[enum],Z_PROC[enum])

        atm1=atmosphere(OBS[enum],OBS_TIMESEC[enum],EL_PROC[enum],AZ_PROC[enum],userSats[enum],preproc,ZHD,MF_Trop,OBS_RATE[enum])
        atm1.iono(GIM,VTEC,MF_Iono)
        atm1.trop(ZHD,MF_Trop)

        process1=process(CTL1[enum],OBS_RATE[enum],OBS[enum],CLK_PROC[enum],MF_Trop[enum],ZHD[0],preproc,userSats[enum])
        process1.CycleSlip(N,CTL2)
        process1.ObsMatrix(L_Matrix,VTEC[enum],MF_Iono[enum])
        
    process1.A_Matrix(L_Matrix,CTL1,X_PROC,Y_PROC,Z_PROC,Vx_s,Vy_s,Vz_s,MF_Trop,N,EL_PROC,SNR,OBS,OBS_TIMESEC,CLK_PROC)
    process1.exe()
    Output(outFile,process1.processed,process1.W_X,process1.EPOCH,process1.sat_num,ZHD,OBS,OBS_RATE[0])

    return process1
        





process1=exe()

