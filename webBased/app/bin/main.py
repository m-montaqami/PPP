from app.bin.Atmosphere3 import atmosphere
from app.bin.Libs import *
from app.bin.DataEntry import dataEntry
from app.bin.Control import control
from app.bin.Output import *
from app.bin.PreProcess4 import preprocess
from app.bin.EqRelated import eqRel


def exe(file_o,file_precise,file_output,user_Sats,user_CutOff,user_SNR_Mask,user_proc_mode,user_CombType,user_Freq,user_SNR_Freq,user_GDOP,plots,user_FileIono=''):
    OBS=[];SP3=[]
    EL=[];AZ=[]
    CTL1=[];CTL2=[]
    DATATIMEGFS=[]
    DATATIMEECMWF=[]
    OBS_TIMESEC=[];SP3_TIMESEC=[];SP3_RATE=[];OBS_RATE=[];OBS_TIMESEC_PROC=[]
    X=[];Y=[];Z=[];CLK=[];Vx=[];Vy=[];Vz=[]
    X_PROC=[];Y_PROC=[];Z_PROC=[];CLK_PROC=[];Vx_s=[];Vy_s=[];Vz_s=[]
    EL_PROC=[];AZ_PROC=[]
    USER_OBS=[]
    SNR=[];AC_SNR=[]
    GIM=[];VTEC=[];MF_Iono=[];MF_Trop=[];ZHD=[]
    L_Matrix=[];N=[]

    dataEntry(OBS,SP3,user_Sats,file_o,file_precise,user_Freq,user_SNR_Freq,'')
    # print('Orginal OBS: '+str(OBS[0].C1C.shape))
    user_TimeInterval=int(OBS[0].interval)
    # print('user_TimeInterval: '+str(user_TimeInterval))
    user_RT=False
    freq=list(user_Freq.keys())
    for enum in range(len(OBS)):
        preproc1=preprocess(OBS[enum],SP3[enum],user_Sats[enum])
        preproc1.TimeProc(OBS_TIMESEC,SP3_TIMESEC,OBS_RATE,SP3_RATE,DATATIMEGFS,DATATIMEECMWF)
        # print('SP3time : '+str(SP3_TIMESEC[enum][0:10]))
        # print('OBStime : '+str(OBS_TIMESEC[enum][0:10]))
        preproc1.sp3Extract(X,Y,Z,Vx,Vy,Vz,CLK)
        preproc1.proc1(X[enum],Y[enum],Z[enum],CLK[enum],OBS_TIMESEC[enum],SP3_TIMESEC[enum],OBS_RATE[enum],SP3_RATE[enum],X_PROC,Y_PROC,Z_PROC,CLK_PROC,Vx_s,Vy_s,Vz_s,OBS_TIMESEC_PROC)
        ctl=control(user_Sats[enum],OBS[enum],SP3[enum],SP3_TIMESEC[enum],OBS_TIMESEC[enum],SP3_RATE[0],OBS_TIMESEC_PROC[enum])
        ctl.clt1(CTL1,OBS_RATE[enum],EL,AZ,EL_PROC,AZ_PROC,preproc1,SNR,AC_SNR,X_PROC[enum],Y_PROC[enum],Z_PROC[enum],user_Freq,user_CombType,user_CutOff,user_SNR_Freq,user_SNR_Mask)
        atm1=atmosphere(OBS[enum],OBS_TIMESEC[enum],EL_PROC[enum],AZ_PROC[enum],user_Sats[enum],preproc1,ZHD,MF_Trop,DATATIMEGFS[enum],DATATIMEECMWF[enum],OBS_RATE[enum],user_Freq,user_CombType)
        atm1.iono(GIM,VTEC,MF_Iono,user_CombType,user_FileIono,user_RT)
        atm1.trop(ZHD,MF_Trop,user_RT)
        eqRel1=eqRel(CTL1[enum],OBS_RATE[enum],OBS[enum],freq[enum],USER_OBS,0,CLK_PROC[enum],MF_Trop[enum],ZHD[0],preproc1,user_Freq,user_CombType)
        eqRel1.CycleSlip(N,CTL2,user_proc_mode,user_Freq)
        eqRel1.ObsMatrix(L_Matrix,user_CombType,user_proc_mode,MF_Iono[enum],VTEC[enum])
        
    eqRel1.A_Matrix(L_Matrix,CTL2,X_PROC,Y_PROC,Z_PROC,Vx_s,Vy_s,Vz_s,MF_Trop,N,EL_PROC,OBS,OBS_TIMESEC,CLK_PROC,user_TimeInterval,user_Sats,user_proc_mode)
    
    eqRel1.exe(user_proc_mode,user_TimeInterval,user_Sats,user_GDOP,SP3_TIMESEC,OBS_TIMESEC,OBS_RATE)
    dataEncoded=Output(file_output,eqRel1.processed,eqRel1.residual,eqRel1.W_X,eqRel1.EPOCH,eqRel1.sat_num,eqRel1.Qx,eqRel1.PHI,ZHD,OBS,OBS_RATE[0],eqRel1.sats_name,eqRel1.out,eqRel1.PREFIT,eqRel1.POSTFIT,user_TimeInterval,user_CombType,file_o,file_precise,user_SNR_Mask,user_CutOff,user_FileIono,plots,SNR,eqRel1.GDOP)
    return dataEncoded
