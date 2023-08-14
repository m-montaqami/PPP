from app.bin.Libs import *


class control():
    def __init__(self,sat,obs,sp3,sp3Sec,obsSec,sp3Rate,TimeNewObs):
        self.obs=obs
        self.sp3=sp3
        self.sat=sat
        self.sp3Sec=sp3Sec
        self.obsSec=obsSec
        self.sp3Rate=sp3Rate
        self.TimeNewObs=TimeNewObs

    def El_Az(self,EL,AZ,EL_PROC,AZ_PROC,f1,x,y,z):
        #Cal Elevation and Azimuth
        r=self.obs.position_geodetic
        sat=self.obs.sv.data[0][0];first=list(self.obs.sv.data);second=list(self.sp3.sv.data)
        first_buff=[]
        for k in first:
            if(k[0] in sat and k in second):first_buff.append(second.index(k))
            else:first_buff.append(-1)
        first=first_buff
        a=self.sp3.position.data;b=np.array(np.zeros([a.shape[0],len(first),a.shape[2]]))
        R=np.array([[-math.sin(math.radians(r[1])),-math.sin(math.radians(r[0]))*math.cos(math.radians(r[1])),math.cos(math.radians(r[0]))*math.cos(math.radians(r[1]))],
                        [math.cos(math.radians(r[1])),-math.sin(math.radians(r[0]))*math.sin(math.radians(r[1])),math.cos(math.radians(r[0]))*math.sin(math.radians(r[1]))],
                        [0,math.cos(math.radians(r[0])),math.sin(math.radians(r[0]))]])
        r_ENU=np.array([np.matmul(R.T,((i.T*1000)-np.array([self.obs.position]).T)) for i in a])    
        # az=np.array([np.rad2deg(np.arctan(-i[1,:]/i[0,:])) for i in r_ENU])
        az=pymap3d.ecef2aer(x.T,y.T,z.T,r[0],r[1],r[2], deg=True)[0]
        # el=np.array([np.rad2deg(np.arctan(i[2,:]/np.sqrt(np.power(i[0,:],2)+np.power(i[1,:],2)))) for i in r_ENU])
        el=pymap3d.ecef2aer(x.T,y.T,z.T,r[0],r[1],r[2], deg=True)[1]
        # print(el.shape)
        el2=[];az2=[]
        for i in range(el.shape[0]):
            buff=[];buff2=[]
            for j in range(len(first)):
                try:buff.append(el[i,first[j]])
                except:buff.append(math.nan)
                try:buff2.append(az[i,first[j]])
                except:buff2.append(math.nan)
            el2.append(buff)
            az2.append(buff2)
        EL.append(np.array(el2))
        AZ.append(np.array(az2))
        buff3=np.array(el2).copy()
        buff4=np.array(az2).copy()
        # EL_PROC.append(buff3.T)
        # AZ_PROC.append(buff4.T)
        EL_PROC.append(el.T)
        AZ_PROC.append(az.T)


    def CutOff(self,el,user_CutOff):
        ac_el=el.copy()
        ac_el[ac_el<user_CutOff]=np.nan     
        ac_el[ac_el>=user_CutOff]=1     #accepted elevation angle(positive and user defined)
        return ac_el

    def SNR(self,obsRate,SNR,AC_SNR,user_SNR_Freq,user_SNR_Mask,user_CombType):
        #temp1=np.array([self.obs[j][0:len(self.obs[j]):(user_TimeInterval//obsRate),:] for j in user_SNR_Freq[self.sat]])
        temp1=np.array([self.obs[j] for j in user_SNR_Freq[self.sat]])
        #based on SNR(signal-to-noise raito which is the raito between the expected signal power to received signal power)
        ac_snr=np.array([temp1[i].data for i in range(temp1.shape[0])])
        SNR.append(ac_snr.copy())
        ac_snr[ac_snr<user_SNR_Mask]=np.nan
        ac_snr[ac_snr>=user_SNR_Mask]=1
        if(user_CombType=='iono-free'):
            try:
                AC_SNR.append(ac_snr[0,:,:].copy()*ac_snr[1,:,:].copy())
            except:
                AC_SNR.append(ac_snr[0,:,:].copy())
        elif(user_CombType=='single-freq'):AC_SNR.append(ac_snr[0,:,:].copy())

    def clt1(self,CTL1,obsRate,EL,AZ,EL_PROC,AZ_PROC,preproc1,SNR,AC_SNR,x,y,z,user_Freq,user_CombType,user_CutOff,user_SNR_Freq,user_SNR_Mask):
        #[P,C]
        control.El_Az(self,EL,AZ,EL_PROC,AZ_PROC,preproc1.Interpolation,x,y,z)
        buff1=[user_Freq[self.sat].index(j) for j in user_Freq[self.sat] if not 'L' in j]
        buff2=[user_Freq[self.sat].index(j) for j in user_Freq[self.sat] if 'L' in j]
        #temp1=np.array([self.obs[j][0:len(self.obs[j]):(user_TimeInterval//obsRate),:] for j in user_Freq[self.sat]]);temp3=temp1.copy();temp3[np.isnan(temp3)]=0;temp3[temp3!=0]=1
        temp1=np.array([self.obs[j] for j in user_Freq[self.sat]]);temp3=temp1.copy();temp3[np.isnan(temp3)]=0;temp3[temp3!=0]=1
        temp4_C=np.ones([temp3.shape[1],temp3.shape[2]]);temp4_P=np.ones([temp3.shape[1],temp3.shape[2]])
        #Combination check
        if(user_CombType=='iono-free'):
            try:
                for j in buff1:temp4_C*=temp3[j,:,:]
                for j in buff2:temp4_P*=temp3[j,:,:]
            except:
                temp4_C*=temp3[buff1[0],:,:]
                temp4_P*=temp3[buff2[0],:,:]

        elif(user_CombType=='single-freq'):
            temp4_C*=temp3[buff1[0],:,:]
            temp4_P*=temp3[buff2[0],:,:]
        temp4_C[temp4_C==0.]=np.nan
        temp4_P[temp4_P==0.]=np.nan

        #Elevation Cut-Off
        el=control.CutOff(self,EL_PROC[-1].T,user_CutOff)
        temp4_C*=el
        temp4_P*=el
        #SNR Mask
        control.SNR(self,obsRate,SNR,AC_SNR,user_SNR_Freq,user_SNR_Mask,user_CombType)
        temp4_P*=AC_SNR[-1]
        temp4_C*=AC_SNR[-1]

        CTL1.append([temp4_P,temp4_C])
 
 