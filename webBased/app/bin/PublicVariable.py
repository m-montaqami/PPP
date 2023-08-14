
#variables



def Start():
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
    L_Matrix=[];N=[];L_Matrix_AV=[];L_Matrix_AV_params=[]
    X_cap=[]
    return OBS,SP3,EL,AZ,CTL1,CTL2,DATATIMEGFS,DATATIMEECMWF,OBS_RATE,SP3_RATE,OBS_TIMESEC,SP3_TIMESEC,X,Y,Z,CLK,X_PROC,Y_PROC,Z_PROC,CLK_PROC,Vx,Vy,Vz,Vx_s,Vy_s,Vz_s,EL_PROC,AZ_PROC,OBS_TIMESEC_PROC,USER_OBS,SNR,AC_SNR,GIM,VTEC,MF_Iono,MF_Trop,ZHD,L_Matrix,L_Matrix_AV,L_Matrix_AV_params,N,X_cap