from Libs import *

class process():
    def __init__(self,ctl1,obsRate,Obs,CLKs,MF_t,zhd,preproc,sat):
        self.sat=sat
        self.noise=[]
        self.ObsRes=[]
        self.inovation=[]
        self.PREFIT,self.POSTFIT,self.Ns=[],[],[]
        self.GDOP=[]
        self.W_X_Ns=[]
        self.preproc=preproc;self.ZHD=[];self.SP=[];self.residual=[];self.sats_name=[]
        self.out,self.delta_rel,self.delta_rel_e,self.sat_mov,self.earth_rot=[],[],[],[],[]
        self.CLKs=CLKs;self.freq=Freq[self.sat][userFreq[self.sat][1][0:2]];self.ctl1=ctl1;self.obsRate=obsRate
        self.Obs=Obs;self.MF_t=MF_t;self.zhd=zhd
        self.resi=[[],[]];self.resi2=[[],[]];self.resi3=[[],[]];self.resi4=[[],[]] #phase,code
        self.PDOP=[]
        self.satsNum=[]
        self.satsNumProcessed=[]
        self.measNum=[[],[]] #phase - code
        self.measNumAccepted=[[],[]] #phase - code
        self.processed,self.W_X,self.Xcap=[],[],[]
        self.EPOCH,self.sat_num=[],[]
        self.Landa=[C/self.freq]
        # self.corr=C*self.CLKs
        if(ionoMode=='iono-free'):
            self.freq2=Freq[self.sat][userFreqF2[self.sat][1][0:2]]
            self.Landa.append(C/self.freq2)

    def CycleSlip(self,N,CTL2):
        
        if(userProcMode!='code-only'):
            if(ionoMode=='iono-free'):
                Phase1=self.Obs[userFreq[self.sat][1]].data*self.Landa[0]
                Phase2=self.Obs[userFreqF2[self.sat][1]].data*self.Landa[1]

                buff=(Phase1[1:,:]-Phase1[0:-1,:])-(Phase2[1:,:]-Phase2[0:-1,:])
                buff_ctl=buff.copy()
                thershold=self.Landa[0]
                buff[np.abs(buff)>thershold]=np.nan
                buff[np.abs(buff)<=thershold]=0
                buff[np.isnan(buff)]=1
                buff_ctl[np.abs(buff_ctl)>thershold]=np.nan
                buff_ctl[np.abs(buff_ctl)<=thershold]=1

                buff0=self.ctl1[0][0,:].copy();buff0[np.isnan(buff0)]=0
                buff1=self.ctl1[0][1,:].copy();buff1[np.isnan(buff1)]=0
                buff00=np.row_stack([buff0,buff1])
                buff=np.row_stack([buff00,buff[0:-1,:]])
                N.append(buff)
                buff_ctl=np.row_stack([self.ctl1[0][0,:].copy(),buff_ctl])

                CTL2.append([buff_ctl,self.ctl1[1]])
            else:
                userDopplerFreq=userFreq[self.sat][2]
                Phase=self.Obs[userFreq[self.sat][1]].data*self.Landa[0]
                Doppler=self.Obs[userDopplerFreq].data*self.Landa[0]
                midDoppler=(Doppler[0:-1,:]+Doppler[1:,:])/2
                sim_phase=(Phase[1:,:]-Phase[0:-1,:])+midDoppler[0:,:]
                sim_phase[np.abs(sim_phase)>60]=2
                sim_phase[sim_phase<-1*self.Landa[0]]=2
                sim_phase[np.abs(sim_phase)>20]=0

                thershold=self.Landa[0]

                buff=sim_phase.copy()
                buff[np.abs(buff)>=thershold]=np.nan   
                buff[np.abs(buff)<thershold]=0
                buff[np.isnan(buff)]=1
                
                buff_ctl=sim_phase.copy()
                buff_ctl[np.abs(buff_ctl)>=thershold]=np.nan
                buff_ctl[np.abs(buff_ctl)<thershold]=1
            
                buff0=self.ctl1[0][0,:].copy();buff0[np.isnan(buff0)]=0
                buff1=self.ctl1[0][1,:].copy();buff1[np.isnan(buff1)]=0               
                buff00=np.row_stack([buff0,buff1])
                buff=np.row_stack([buff00,buff[0:-1,:]])
                N.append(buff)
                buff_ctl=np.row_stack([self.ctl1[0][0,:].copy(),buff_ctl])
                CTL2.append([buff_ctl,self.ctl1[1]])
               
        else:
            CTL2.append([self.ctl1[0],self.ctl1[1]])

    def ObsMatrix(self,L_Matrix,vTEC,MF_iono):
        if(ionoMode=='iono-free'):
            
            L_P=self.Obs[userFreq[self.sat][1]].data.T*self.Landa[0]
            L_P2=self.Obs[userFreqF2[self.sat][1]].data.T*self.Landa[1]
            L_C=self.Obs[userFreq[self.sat][0]].data.T
            L_C2=self.Obs[userFreqF2[self.sat][0]].data.T
            alpha=np.power(self.freq,2)/(np.power(self.freq,2)-np.power(self.freq2,2))
            beta=np.power(self.freq2,2)/(np.power(self.freq,2)-np.power(self.freq2,2))
            L_P=alpha*L_P-beta*L_P2
            L_C=alpha*L_C-beta*L_C2

        else:
            L_P=self.Obs[userFreq[self.sat][1]].data.T*self.Landa[0]
            L_P-=(MF_iono*np.array([vTEC[:,0]]))
            L_C=self.Obs[userFreq[self.sat][0]].data.T
            L_C-=(MF_iono*np.array([vTEC[:,1]]))
        if(userProcMode=='combined'):
            L_Matrix.append([L_P,L_C])
        elif(userProcMode=='code-only'):
            L_Matrix.append([L_C])
        
    def A_Matrix(self,L_Matrix,ctl,Xs,Ys,Zs,Vx_s,Vy_s,Vz_s,MF_t,N,Elevation,SNR,OBS,OBS_TIME,clk):
        self.L_eq=[];self.X_eq=[];self.Y_eq=[];self.Z_eq=[];self.MF_eq=[]
        self.A_eq=[];self.N_eq=[]
        self.CLK_eq=[];self.time_eq=[]
        self.W_L_eq=[];self.Vx_eq=[];self.Vy_eq=[];self.Vz_eq=[];self.N2=[]
        self.sat_nums=[];self.sat_name=[];self.MF_zhd_eq=[]
        self.N_copy=[];self.L_len=[]
        self.ind_sat_eq=[]

        interv=int(self.obsRate)
        for ep in range(0,L_Matrix[0][0].shape[1]-1,interv):
            ISB=[];isb_range=[0 for i in range(len(userSats))]
            self.L_eq.append([]);self.X_eq.append([]);self.Y_eq.append([]);self.Z_eq.append([]);self.MF_eq.append([])
            self.A_eq.append([]);self.N_eq.append([])
            self.CLK_eq.append([]);self.time_eq.append([])
            self.W_L_eq.append([]);self.Vx_eq.append([])
            self.Vy_eq.append([]);self.Vz_eq.append([]);self.N2.append([]);self.sat_nums.append([]);self.sat_name.append([])
            self.MF_zhd_eq.append([])
            self.N_copy.append([]);self.L_len.append([])
            self.ind_sat_eq.append([])
            for s in range(len(userSats)):
                if(s>0):buff_s=len(self.N_eq[-1])
                else:buff_s=0
                shap=L_Matrix[s][0].shape[0]*interv
                if(userProcMode=='combined' and L_Matrix[s][0][:,ep:ep+interv].shape[0]*L_Matrix[s][0][:,ep:ep+interv].shape[1]==shap):
                    #phase
                    self.sat_nums[-1].append(np.unique(np.where(~np.isnan(ctl[s][0][ep:ep+interv,:]))[1]).shape[0])
                    self.sat_name[-1].append([OBS[s].sv.data[k] for k in np.unique(np.where(~np.isnan(ctl[s][0][ep:ep+interv,:]))[1])])
                    temp_test=N[s][ep:ep+interv,:]*ctl[s][0][ep:ep+interv,:]
                    temp_test2=temp_test.copy();temp_test2[~np.isnan(temp_test2)]=1
                    for i in range(temp_test.shape[1]):
                        if(1 in list(temp_test[:,i])):temp_test[:,i]=1
                    temp_test=temp_test*temp_test2
                    self.ind_sat_eq[-1].append(np.reshape(ctl[s][0][ep:ep+interv,:],[shap,1]))
                    if(s==0):
                        self.N_copy[-1].append(temp_test)
                    else:
                        self.N_copy[-1].append(temp_test.copy())
                    temp=N[s][ep:ep+interv,:]*ctl[s][0][ep:ep+interv,:];
                    temp2=temp.copy();temp2[~np.isnan(temp2)]=1
                    for i in range(temp.shape[1]):
                        if(1 in temp[:,i]):temp[:,i]=1
                    temp=temp*temp2
                    temp[temp==-1]=1
                    [self.N2[-1].append(enum) for enum in temp]
                    [self.N_eq[-1].append(enum+buff_s) for enum in temp[~np.isnan(temp)]]
                    temp=np.reshape(L_Matrix[s][0][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    isb_range[s]+=len(temp)
                    [self.L_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(clk[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.CLK_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Xs[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.X_eq[-1].append(enum) for enum in temp]
                    temp=OBS_TIME[s][ep:ep+interv]*len(temp)
                    [self.time_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Ys[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Y_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Zs[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Z_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Vx_s[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Vx_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Vy_s[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Vy_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Vz_s[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Vz_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(MF_t[s][0][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.MF_zhd_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(MF_t[s][1][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.MF_eq[-1].append(enum) for enum in temp]

                    if(userMeasNoise=='CN0'):
                        SNR=SNR.T
                        temp=np.reshape(SNR[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power((customNoise[1])/np.sin(np.radians(enum)),2)) for enum in temp]
                    elif(userMeasNoise=='elev'):
                        temp=np.reshape(Elevation[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power((customNoise[1])/np.sin(np.radians(enum)),2)) for enum in temp]
                    elif(userMeasNoise=='constant'):
                        temp=np.reshape(Elevation[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power(customNoise[1],2)) for enum in temp]
                    elif(userMeasNoise=='custom'):
                        temp=np.array([customNoise[1][s]])
                        temp*=ctl[s][0][ep,:]
                        temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power(enum,2)) for enum in temp]


                    #code
                    self.ind_sat_eq[-1].append(np.reshape(ctl[s][1][ep:ep+interv,:],[shap,1]))
                    self.sat_nums[-1].append(np.unique(np.where(~np.isnan(ctl[s][1][ep:ep+interv,:]))[1]).shape[0])
                    self.sat_name[-1].append([OBS[s].sv.data[k] for k in np.unique(np.where(~np.isnan(ctl[s][1][ep:ep+interv,:]))[1])])
                    temp=np.reshape(L_Matrix[s][1][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    isb_range[s]+=len(temp)
                    [self.L_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(clk[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.CLK_eq[-1].append(enum) for enum in temp]
                    [self.N_eq[-1].append(0) for enum in range(len(temp))]
                    temp=np.reshape(Xs[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.X_eq[-1].append(enum) for enum in temp]
                    self.L_len[-1].append(temp)
                    temp=OBS_TIME[s][ep:ep+interv]*len(temp)
                    [self.time_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Ys[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Y_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Zs[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Z_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Vx_s[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Vx_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Vy_s[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Vy_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Vz_s[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Vz_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(MF_t[s][0][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.MF_zhd_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(MF_t[s][1][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.MF_eq[-1].append(enum) for enum in temp]

                    if(userMeasNoise=='CN0'):
                        SNR=SNR.T
                        temp=np.reshape(SNR[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power((customNoise[0])/np.sin(np.radians(enum)),2)) for enum in temp]
                    elif(userMeasNoise=='elev'):
                        temp=np.reshape(Elevation[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power((customNoise[0])/np.sin(np.radians(enum)),2)) for enum in temp]
                    elif(userMeasNoise=='constant'):
                        temp=np.reshape(Elevation[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power(customNoise[0],2)) for enum in temp]
                    elif(userMeasNoise=='custom'):
                        temp=np.array([customNoise[0][s]])
                        temp*=ctl[s][1][ep,:]
                        temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power(enum,2)) for enum in temp]

                if(userProcMode=='code-only' and L_Matrix[s][0][:,ep:ep+interv].shape[0]*L_Matrix[s][0][:,ep:ep+interv].shape[1]==shap):
                    self.ind_sat_eq[-1].append(np.reshape(ctl[s][1][ep:ep+interv,:],[shap,1]))
                    self.sat_nums[-1].append(np.unique(np.where(~np.isnan(ctl[s][1][ep:ep+interv,:]))[1]).shape[0])
                    self.sat_name[-1].append([OBS[s].sv.data[k] for k in np.unique(np.where(~np.isnan(ctl[s][1][ep:ep+interv,:]))[1])])
                    temp=np.reshape(L_Matrix[s][0][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    isb_range[s]+=len(temp)
                    [self.L_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(clk[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.CLK_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Xs[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.X_eq[-1].append(enum) for enum in temp]
                    temp=OBS_TIME[s][ep:ep+interv]*len(temp)
                    [self.time_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Ys[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Y_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Zs[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Z_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Vx_s[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Vx_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Vy_s[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Vy_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(Vz_s[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.Vz_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(MF_t[s][0][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.MF_zhd_eq[-1].append(enum) for enum in temp]
                    temp=np.reshape(MF_t[s][1][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                    [self.MF_eq[-1].append(enum) for enum in temp]

                    if(userMeasNoise=='CN0'):
                        SNR=SNR.T
                        temp=np.reshape(SNR[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power((customNoise[0])/np.sin(np.radians(enum)),2)) for enum in temp]
                    elif(userMeasNoise=='elev'):
                        temp=np.reshape(Elevation[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power((customNoise[0])/np.sin(np.radians(enum)),2)) for enum in temp]
                    elif(userMeasNoise=='constant'):
                        temp=np.reshape(Elevation[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power(customNoise[0],2)) for enum in temp]
                    elif(userMeasNoise=='custom'):
                        temp=np.array([customNoise[0][s]])
                        temp*=ctl[s][1][ep,:]
                        temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power(enum,2)) for enum in temp]

                if(len(userSats)>1):
                    if(userSats[s]=='G'):
                        if(len(temp>0)):
                            ISB.append([0 for ii in range(isb_range[s])])
                    else:
                        if(len(temp>0)):
                            ISB.append([1 for ii in range(isb_range[s])])
                        

            self.A_eq[-1]=np.array([self.MF_eq[-1]]).T
            self.A_eq[-1]=np.column_stack([self.A_eq[-1],np.ones([self.A_eq[-1].shape[0],1])])
            if(len(ISB)>1):

                ISB=ISB[0]+ISB[1]
                ISB=np.array([ISB]).T
                self.A_eq[-1]=np.column_stack([self.A_eq[-1],ISB])


            if(userProcMode!='code-only'):
                temp0=self.N_copy[-1][0].copy()
                temp0[~np.isnan(temp0)]=1
                buff1=temp0[:,np.any(~np.isnan(temp0), axis=0)]
                buff1[buff1==0]=1
                temp0=np.diag(list(buff1[0,:]))
                for enum0 in range(1,buff1.shape[0]):temp0=np.row_stack([temp0,np.diag(list(buff1[enum0,:]))])
                temp0=np.row_stack([temp0,np.zeros([np.shape(self.L_len[-1][0])[0],temp0.shape[1]])])

                for enum in range(1,len(userSats)):
                    temp0_size=temp0.shape[1]
                    buffer=self.N_copy[-1][enum].copy()
                    buffer[~np.isnan(buffer)]=1
                    buff1=buffer[:,np.any(~np.isnan(buffer), axis=0)]
                    temp0=np.column_stack([temp0,np.zeros([temp0.shape[0],buff1.shape[1]])])
                    temp00=np.diag(list(buff1[0,:]))
                    for enum0 in range(1,buff1.shape[0]):temp00=np.row_stack([temp00,np.diag(list(buff1[enum0,:]))])
                    temp00=np.column_stack([np.zeros([temp00.shape[0],temp0_size]),temp00])
                    temp0=np.row_stack([temp0,temp00])
                    temp0=np.row_stack([temp0,np.zeros([np.shape(self.L_len[-1][enum])[0],temp0.shape[1]])])

                temp0=temp0[np.all(~np.isnan(temp0), axis=1)]
                buff2=temp0.copy()
                if(buff2.shape[1]!=0):
                    self.A_eq[-1]=np.column_stack([self.A_eq[-1],buff2])
                     
    def PostResControl(self,A,L,X0,W_L,N_copy,Lres,ind_sat,init,c_v,res):

        W_L=list(np.diag(W_L))
        X0=X0[0:(5+len(userSats)-1)]
        if(userProcMode=='combined'):N_copy=np.column_stack([i for i in N_copy])
        L2,A2,W_L2,N_buff2=[],[],[],[]
        phi=self.PHI.copy();qx=self.Qx.copy()
        enum3,enum4=[],[]
        enum5=0

        if(init=='normal'):W_X=self.W_X[-1][0:(5+len(userSats)-1),0:(5+len(userSats)-1)].copy()
        
        elif(init=='init'):
            buffer=kfParameters['P0'].copy()[0:-2]
            #isb
            [buffer.append(kfParameters['P0'].copy()[-2]) for i in range(len(userSats)-1)]

            W_X=np.diag(buffer)
        k=-1;k2=-1
        buffer_measNumAccepted0=0
        buffer_measNumAccepted1=0
        for i in range(len(ind_sat)):
            for j in range(len(ind_sat[i])):

                if(np.isnan(ind_sat[i][j])==False):
                    k+=1
                    if(userProcMode=='combined'):
                        
                        if(i==0 or i==2):
                            k2+=1

                            thresholdPhase=phaseThreshold*np.sqrt(W_L[k])
                            print('threPhase:'+str(thresholdPhase))

                            validation=np.abs(Lres[k,0]/(1/np.sqrt(W_L[k])))
                            print('Phase-validation:'+str(validation))

                            if(validation<=thresholdPhase):
                                buffer_measNumAccepted0+=1
                                L2.append(L[k,0])
                                A2.append(A[k,:])
                                W_L2.append(W_L[k])
                                N_buff2.append(1)
                            else:
                                ind_sat[i][j]=np.nan
                                self.W_X_N[0,k2]=0
                                N_copy[0,k2]=np.nan
                        if(i==1 or i==3):
                            thresholdCode=codeThreshold*np.sqrt(W_L[k])
                            print('threCode:'+str(thresholdCode))

                            validation=np.abs(Lres[k]/(np.sqrt(W_L[k])))
                            print('Code-validation:'+str(validation))

                            if(validation<=thresholdCode):
                                buffer_measNumAccepted1+=1
                                L2.append(L[k,0])
                                A2.append(A[k,:])
                                W_L2.append(W_L[k])
                            else:
                                ind_sat[i][j]=np.nan
                    elif(userProcMode=='code-only'):

                        thresholdCode=codeThreshold
                        validation=np.abs(Lres[k]/(np.sqrt(W_L[k])))
                        print('validation:'+str(validation))
                        if(validation<=thresholdCode):
                            buffer_measNumAccepted1+=1
                            L2.append(L[k,0])
                            A2.append(A[k,:])
                            W_L2.append(W_L[k])
                        else:
                            ind_sat[i][j]=np.nan
                else:
                    if(userProcMode=='combined'):
                        if(i==0 or i==2):
                            k2+=1
                            N_buff2.append(np.nan)
                            
        if(userProcMode!='code-only'):
            for enum1 in range(len(N_copy[0,:])):
                if(np.isnan(N_copy[0,enum1])==False):
                    if(N_copy[0,enum1]==1):
                        enum4.append(0)
                        enum3.append(kfParameters['P0'][-1])
                        self.W_X_N[0,enum1]!=0
                    elif(N_copy[0,enum1]==0 and self.W_X_N[0,enum1]!=0):
                        enum4.append(self.N[0,enum1])
                        enum3.append(self.W_X_N[0,enum1])
                        enum5+=1
                    elif(N_copy[0,enum1]==0 and self.W_X_N[0,enum1]==0):
                        enum4.append(0)
                        enum3.append(kfParameters['P0'][-1])
                        self.W_X_N[0,enum1]!=0
            N_copy=N_copy.T
            W_X1=np.zeros([W_X.shape[0]+len(enum3),W_X.shape[1]+len(enum3)])
            
            X0=np.array([np.append(X0,enum4)]).T
            phi=list(np.append(phi,[1 for i in range(len(enum3))]))
            qx=list(np.append(qx,[0 for i in range(len(enum3))]))
            W_X1[0:(5+len(userSats)-1),0:(5+len(userSats)-1)]=W_X.copy()
            if(enum3!=[]):np.fill_diagonal(W_X1[(5+len(userSats)-1):,(5+len(userSats)-1):],enum3)
            W_X=W_X1.copy()
        W_L2=np.diag(W_L2)
        A2=np.array(A2)
        A2=A2[:,~np.all(A2==0,axis=0)]
        print(A2.shape)
        if(A2.shape[0]<A2.shape[1] or len(A2)==0):
            if(res==False):
                self.satsNumProcessed.append(np.nan)
                self.measNumAccepted[1].append(np.nan)
                self.measNumAccepted[0].append(np.nan)
            return [],[],[],[],[],[],[],[]
        phi=np.diag(phi)

        qx=np.diag(qx)
        if(A2.shape[1]!=phi.shape[0]):
            if(res==False):
                self.satsNumProcessed.append(np.nan)
                self.measNumAccepted[1].append(np.nan)
                self.measNumAccepted[0].append(np.nan)
            return [],[],[],[],[],[],[],[]
        if(res==False):
            self.satsNumProcessed.append(buffer_measNumAccepted1)
            self.measNumAccepted[1].append(buffer_measNumAccepted1)
            self.measNumAccepted[0].append(buffer_measNumAccepted0)

        pre_x=np.matmul(phi,X0)
        pre_px=np.matmul(np.matmul(phi,W_X),phi.T)+qx
        L2=np.array([L2]).T
        L002=L2.copy()
        L2=L2.copy()-np.matmul(A2,pre_x)
        K=np.matmul(np.matmul(pre_px,A2.T),np.linalg.pinv(np.matmul(np.matmul(A2,pre_px),A2.T)+W_L2))
        Xcap=pre_x+np.matmul(K,L2.copy())
        post_px=np.matmul((np.eye(A2.shape[1])-np.matmul(K,A2)),pre_px)
        post_px=np.diag(list(np.diag(post_px)))

        if(res==True):
            self.residual.append(np.matmul(K,L2.copy()))
        return N_buff2,Xcap,post_px,A2,L002,N_copy,W_L2,ind_sat

    def LS(self,A,L,W_L,ind_sat,buffer):
        inoVector=L-buffer
        A2,L2_p,L2_L,W_L2,W_P2=[],[],[],[],[]
        W_L=list(np.diag(W_L))
        k=-1
        for i in range(len(ind_sat)):
            for j in range(len(ind_sat[i])):

                if(np.isnan(ind_sat[i][j])==False):
                    k+=1
                    if(userProcMode=='combined'):
                        if(i==1 or i==3):
                            L2_p.append(inoVector[k,0])
                            A2.append(A[k,:])
                            W_P2.append(W_L[k])
                        if(i==0 or i==2):
                            L2_L.append(inoVector[k,0])
                            W_L2.append(W_L[k])
                    elif(userProcMode=='code-only'):
                        L2_p.append(inoVector[k,0])
                        A2.append(A[k,:])
                        W_L2.append(W_L[k])
                        W_L2.append(W_L[k])
        
        A2=np.array(A2)
        if(userProcMode=='code-only'):
            L2_L=L2_p.copy()
            W_L2=W_P2.copy()
        
        
        L2_p=np.array([L2_p]).T
        L2_L=np.array([L2_L]).T
        self.measNum[1].append(L2_p.shape[0])
        self.measNum[0].append(L2_L.shape[0])
        W_P2=np.diag(W_P2)
        W_L2=np.diag(W_L2)
        if(A2.shape[0]<A2.shape[1] or len(A2)==0):
            self.PDOP.append(np.nan)
            return np.nan,np.nan,np.nan,np.nan,np.nan
        else:
            buffPDOP=np.diag(np.linalg.inv(np.matmul(A.T,A)))
            PDOP=np.sqrt(buffPDOP[0]+buffPDOP[1]+buffPDOP[2])
            self.PDOP.append(PDOP)
            X_cap=np.matmul(np.linalg.inv(np.matmul(A2.T,A2)),np.matmul(A2.T,L2_p))
            return X_cap[3],L2_p-X_cap[3],L2_L-X_cap[3],W_P2,W_L2

    def KalmanFilter(self,X0,L,A,Xs,Ys,Zs,Vx,Vy,Vz,W_L,init,MF_zhd,clk,N_copy,ind_sat):
        prefit,postfit=[],[]        
        phi=self.PHI.copy();qx=self.Qx.copy()
        
        enum3,enum4=[],[]
        enum5=0
        if(init=='normal'):W_X=self.W_X[-1][0:(5+len(userSats)-1),0:(5+len(userSats)-1)].copy()
        
        elif(init=='init'):
            buffer=kfParameters['P0'].copy()[0:-2]
            #isb
            [buffer.append(kfParameters['P0'].copy()[-2]) for i in range(len(userSats)-1)]
            W_X=np.diag(buffer)
        N_buff2=[]
        if(userProcMode!='code-only'):
            for enum0 in range(len(userSats)):
                N_buff=N_copy[enum0][0].copy()
                for enum1 in range(len(N_buff)):
                    if(np.isnan(N_buff[enum1])==False):
                        if(N_buff[enum1]==1):
                            enum4.append(0)
                            enum3.append(kfParameters['P0'][-1])
                            #!
                        elif(N_buff[enum1]==0 and self.W_X_N[0,enum1+(len(N_buff)*enum0)]!=0):
                            enum4.append(self.N[0,enum1+(len(N_buff)*enum0)])
                            enum3.append(self.W_X_N[0,enum1+(len(N_buff)*enum0)])
                            enum5+=1
                        elif(N_buff[enum1]==0 and self.W_X_N[0,enum1+(len(N_buff)*enum0)]==0):
                            enum4.append(0)
                            enum3.append(kfParameters['P0'][-1])


            W_X1=np.zeros([W_X.shape[0]+len(enum3),W_X.shape[1]+len(enum3)])
            X0=X0[0:(5+len(userSats)-1)]
            X0=np.array([np.append(X0,enum4)]).T
            phi=list(np.append(phi,[1 for i in range(len(enum3))]))
            qx=list(np.append(qx,[0 for i in range(len(enum3))]))
            W_X1[0:(5+len(userSats)-1),0:(5+len(userSats)-1)]=W_X.copy()
            if(enum3!=[]):np.fill_diagonal(W_X1[(5+len(userSats)-1):,(5+len(userSats)-1):],enum3)
            W_X=W_X1.copy()

        phi=np.diag(phi)
        qx=np.diag(qx)
        x0=X0.copy()
        A1=A.copy()

        ru=np.sqrt(np.power(x0[0]-Xs,2)+np.power(x0[1]-Ys,2)+np.power(x0[2]-Zs,2))
        prefit.append(L.copy())
        prefit.append(ru.copy())
        cor_rot=((x0[0]-Xs)*(-OMEGA_Earth*x0[1])+(x0[1]-Ys)*(OMEGA_Earth*x0[0]))/C
        self.earth_rot.append(cor_rot)
        self.sat_mov.append([-1*(ru/C)*Vx,-1*(ru/C)*Vy,-1*(ru/C)*Vz])
        Xs+=self.sat_mov[-1][0]
        Ys+=self.sat_mov[-1][1]
        Zs+=self.sat_mov[-1][2]
        A1=np.column_stack([np.column_stack([(x0[0]-Xs)/ru,
                        (x0[1]-Ys)/ru,
                        (x0[2]-Zs)/ru]),A1])
        ru_s=np.sqrt(np.power(Xs,2)+np.power(Ys,2)+np.power(Zs,2))
        ru_r=np.sqrt(np.power(x0[0],2)+np.power(x0[1],2)+np.power(x0[2],2))
        self.delta_rel_e.append((2/C)*(Xs*Vx+Ys*Vy+Zs*Vz))
        self.delta_rel.append(((2*MU)/np.power(C,2))*np.log((ru_r+ru_s+ru)/(ru_r+ru_s-ru)))
        ru=np.sqrt(np.power(x0[0]-Xs,2)+np.power(x0[1]-Ys,2)+np.power(x0[2]-Zs,2))-clk+self.delta_rel[-1]+self.delta_rel_e[-1]+self.earth_rot[-1]+(self.zhd*MF_zhd)
        prefit.append(ru.copy())
        L1=L.copy()
        L1=L1.copy()-ru
        L00=L1.copy()

        if(init=='init'):
            self.W_X.append(W_X)
            
        x0[0:3]-=self.processed[0][0:3]
        pre_x=np.matmul(phi,x0)
    
        bufferPre_px=np.zeros([W_X.shape[0],W_X.shape[1]])
        if(init!='init'):bufferPre_px[0,0]=W_X[0,0];bufferPre_px[1,1]=W_X[1,1];bufferPre_px[2,2]=W_X[2,2]

        pre_px=np.matmul(np.matmul(phi,W_X),phi.T)
        c_v=np.matmul(np.matmul(A1,pre_px),A1.T)+W_L.copy()
        c_v=np.diag(c_v)
        # print('c_v: '+str(c_v))

        Lres=L00-np.matmul(A1,pre_x)
        LSdef=self.LS(np.column_stack([A1[:,0:3],A1[:,4]]),Lres,W_L,ind_sat,np.matmul(A1,pre_x))

        if(np.isnan(LSdef[0][0])==False):
            Lres-=LSdef[0][0]
            
        N_buff2,Xcap,post_px,A2,L2,N_copy2,W_L2,ind_sat2=self.PostResControl(A1,L00,x0.copy(),W_L,N_copy,Lres,ind_sat,init,c_v,res=True)
        if(Xcap==[]):
            return 'FDE_not_accepted' 

        Xcap[0:3]+=self.processed[0][0:3]
        if(A2.shape[0]>A2.shape[1]):
            if(userProcMode=='code-only'):
                self.POSTFIT.append(postfit)
            self.PREFIT.append(prefit)
            if(userProcMode!='code-only'):
                self.W_X_Ns.append(self.W_X_N.copy())
                self.POSTFIT.append(postfit)
                buff_Nindex=list(Xcap[(5+len(userSats)-1):][:,0])
                buff_W_N_index=list(np.diag(post_px[(5+len(userSats)-1):,(5+len(userSats)-1):]))
                jj=0   
                for enum1 in range(len(N_buff2)):
                    if(np.isnan(N_buff2[enum1])==False):
                        self.N[0,enum1]=buff_Nindex[jj]
                        self.W_X_N[0,enum1]=buff_W_N_index[jj]
                        jj+=1
            print('Xcap: '+str(Xcap))
            self.processed.append(Xcap)
            self.W_X.append(post_px)
            self.Ns.append(self.N.copy())
            return 'FDE_accepted'
        else:
            return 'FDE_not_accepted'    

    def exe(self):
        self.Qx=kfParameters['Qx'][0:-2]
        self.PHI=kfParameters['PHI'][0:-2]
        #isb
        for i in range(len(userSats)-1):
            self.Qx.append(kfParameters['Qx'][-2])
            self.PHI.append(kfParameters['PHI'][-2])

        if(initPos==[]):
            buffer=[self.Obs.position[0],self.Obs.position[1],self.Obs.position[2],0,0]
        else:
            buffer=initPos+[0,0]
        for i in range(len(userSats)-1):
            buffer.append(0)
        self.processed.append(np.array([buffer]).T)
        if(userProcMode!='code-only'):
            self.N,self.W_X_N=np.zeros([1,np.sum([i.shape[1] for i in self.N_copy[0]])]),np.zeros([1,np.sum([i.shape[1] for i in self.N_copy[0]])])
        else:
            self.N,self.W_X_N=[],[]
        enum2=0

        for enum in range(0,len(self.X_eq)-1):
            x0=self.processed[enum2]
            Xs=np.array([self.X_eq[enum]]).T
            Ys=np.array([self.Y_eq[enum]]).T
            Zs=np.array([self.Z_eq[enum]]).T
            Vx=np.array([self.Vx_eq[enum]]).T
            Vy=np.array([self.Vy_eq[enum]]).T
            Vz=np.array([self.Vz_eq[enum]]).T
            MF_zhd=np.array([self.MF_zhd_eq[enum]]).T
            clk=np.array([self.CLK_eq[enum]]).T*C
            if(userProcMode=='combined'):
                N_copy=self.N_copy[enum].copy()
            else:
                N_copy=[]

            L=np.array([self.L_eq[enum]]).T
            
            A=self.A_eq[enum]
            W_L=np.diag(self.W_L_eq[enum])
            if(L.shape[0]>=(3+A.shape[1])):
                print('ep: '+str(enum))
                init='normal'
                if(enum2==0):
                    init='init'            
                res=process.KalmanFilter(self,x0,L,A,Xs,Ys,Zs,Vx,Vy,Vz,W_L,init,MF_zhd,clk,N_copy,self.ind_sat_eq[enum])
                if(res=='FDE_accepted'):
                    self.EPOCH.append(enum);self.sat_num.append(self.sat_nums[enum])
                    # temp=[[self.sats_name[-1].append(j) for j in i] for i in self.sat_name[enum]]

                    enum2+=1
                elif(res=='FDE_not_accepted' and self.W_X_N!=[]):
                    self.N[self.N!=0]=0
                    self.W_X_N[self.W_X_N!=0]=0

            else:
                if(self.W_X_N!=[]):
                    self.N[self.N!=0]=0
                    self.W_X_N[self.W_X_N!=0]=0