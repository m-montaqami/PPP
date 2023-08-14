from Libs import *


class control():
    def __init__(self,sat,obs,sp3,sp3Sec,obsSec,sp3Rate,TimeNewObs):
        self.obs=obs
        self.sp3=sp3
        self.sat=sat
        self.sp3Sec=sp3Sec
        self.obsSec=obsSec
        self.sp3Rate=sp3Rate
        self.TimeNewObs=TimeNewObs

    def El_Az(self,EL_PROC,AZ_PROC,x,y,z):
        #Cal Elevation and Azimuth
        r=self.obs.position_geodetic

        az=pymap3d.ecef2aer(x.T,y.T,z.T,r[0],r[1],r[2], deg=True)[0]        
        el=pymap3d.ecef2aer(x.T,y.T,z.T,r[0],r[1],r[2], deg=True)[1]
        
        EL_PROC.append(el.T)
        AZ_PROC.append(az.T)

    def CutOff(self,el):
        ac_el=el.copy()
        ac_el[ac_el<userCutOff]=np.nan     
        ac_el[ac_el>=userCutOff]=1     #accepted elevation angle(positive and user defined)
        return ac_el

    def SNR(self,SNR,AC_SNR):
        userSNRFreq='S'+userFreq[self.sat][0][1:]
        ac_snr=self.obs[userSNRFreq].data
        SNR.append(ac_snr.copy())
        ac_snr[ac_snr<userSNRMask]=np.nan
        ac_snr[ac_snr>=userSNRMask]=1
        AC_SNR.append(ac_snr.copy())

    def clt1(self,CTL1,EL_PROC,AZ_PROC,SNR,AC_SNR,x,y,z):
        #[P,C]
        control.El_Az(self,EL_PROC,AZ_PROC,x,y,z)
        temp1=np.array([self.obs[j] for j in userFreq[self.sat]]);temp3=temp1.copy();temp3[np.isnan(temp3)]=0;temp3[temp3!=0]=1
        if(ionoMode=='iono-free'):
            temp1F2=np.array([self.obs[j] for j in userFreqF2[self.sat]]);temp3F2=temp1F2.copy();temp3F2[np.isnan(temp3F2)]=0;temp3F2[temp3F2!=0]=1
        temp4_C=np.ones([temp3.shape[1],temp3.shape[2]]);temp4_P=np.ones([temp3.shape[1],temp3.shape[2]])

        #Combination check
        temp4_C*=temp3[0,:,:]
        temp4_P*=temp3[1,:,:]
        if(ionoMode=='iono-free'):
            temp4_C*=temp3F2[0,:,:]
            temp4_P*=temp3F2[1,:,:]
        temp4_C[temp4_C==0.]=np.nan
        temp4_P[temp4_P==0.]=np.nan

        #Elevation Cut-Off
        el=control.CutOff(self,EL_PROC[-1].T)
        temp4_C*=el
        temp4_P*=el

        #SNR Mask
        control.SNR(self,SNR,AC_SNR)
        temp4_P*=AC_SNR[-1]
        temp4_C*=AC_SNR[-1]


        CTL1.append([temp4_P,temp4_C])
 
 