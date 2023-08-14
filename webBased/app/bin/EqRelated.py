from app.bin.Libs import *



class eqRel():
    def __init__(self,ctl1,obsRate,Obs,freq,USER_OBS,Corr,CLKs,MF_t,zhd,preproc,user_Freq,user_CombType):
        user_CycleSlipType='GF'
        self.tau,self.av,self.av_params,self.AV_Var=[],[],[],[]
        self.processed_AV=[]
        self.ObsRes=[]
        self.inovation=[]
        self.PREFIT,self.POSTFIT,self.Ns=[],[],[]
        self.SolidEarthTide=[]
        self.user_CombType=user_CombType
        self.clk_init=0
        
        
        self.GDOP=[]
        self.preproc=preproc;self.ZHD=[];self.SP=[];self.residual=[];self.sats_name=[]
        self.out,self.delta_rel,self.delta_rel_e,self.sat_mov,self.earth_rot=[],[],[],[],[]
        self.CLKs=CLKs;self.freq=freq;self.ctl1=ctl1;self.obsRate=obsRate
        self.Obs=Obs;self.corr=Corr;self.MF_t=MF_t;self.zhd=zhd;
        self.Qx=[0,0,0,1.00e-004,9.00e+10]
        self.PHI=[1,1,1,1,0]
        self.processed,self.W_X,self.Xcap=[],[],[]
        self.freq_val=user_Freq[self.freq]
        self.freq_val=[enum[0:2] for enum in self.freq_val]
        self.L=[user_Freq[self.freq].index(i) for i in user_Freq[self.freq] if 'L' in i]   #carrier-phase index
        self.P=[user_Freq[self.freq].index(i) for i in user_Freq[self.freq] if not 'L' in i]    #code index
        self.obs=np.array([Obs[j] for j in user_Freq[self.freq]])
        self.EPOCH,self.sat_num=[],[]
        USER_OBS.append(self.obs)
        self.Landa=[C/Freq[self.freq][enum[0:2]] for enum in user_Freq[self.freq] if 'L' in enum]
        self.corr=(C*self.CLKs)
        if(user_CombType=='iono-free'):
            self.alpha=(math.pow(Freq[self.freq][self.freq_val[self.L[0]]],2))/(math.pow(Freq[self.freq][self.freq_val[self.L[0]]],2)-math.pow(Freq[self.freq][self.freq_val[self.L[1]]],2))
            self.beta=(math.pow(Freq[self.freq][self.freq_val[self.L[1]]],2))/(math.pow(Freq[self.freq][self.freq_val[self.L[0]]],2)-math.pow(Freq[self.freq][self.freq_val[self.L[1]]],2))
            #freq Ionospheric-free
            self.freqComb=self.alpha*Freq[self.freq][self.freq_val[self.L[0]]]-self.beta*Freq[self.freq][self.freq_val[self.L[1]]]
            self.cs=user_CycleSlipType
        else:
            self.freqComb=Freq[self.freq][self.freq_val[self.L[0]]]
            self.cs='phase-code'

    def CycleSlip(self,N,CTL2,user_proc_mode,user_Freq):
        buff3=[]

        if(user_proc_mode!='code-only' and self.user_CombType=='iono-free'):
            buff2_ctl=np.array([i[:,:] for i in self.obs])
            buff2_ctl=C/Freq[self.freq][user_Freq[self.freq][self.L[0]][0:2]]*buff2_ctl[self.L[0],:,:]-C/Freq[self.freq][user_Freq[self.freq][self.L[1]][0:2]]*buff2_ctl[self.L[1],:,:]

            buff=np.array([i[:-1,:]-i[1:,:] for i in self.obs])
            buff_ctl_code=buff[self.P[0],:,:]-buff[self.P[1],:,:]
            buff=C/Freq[self.freq][user_Freq[self.freq][self.L[0]][0:2]]*buff[self.L[0],:,:]-C/Freq[self.freq][user_Freq[self.freq][self.L[1]][0:2]]*buff[self.L[1],:,:]
            buff_ctl=buff.copy()
            coef=[0.25,1]
            #thershold
            a0=(coef[0])*(C/Freq[self.freq][user_Freq[self.freq][self.L[1]][0:2]]-C/Freq[self.freq][user_Freq[self.freq][self.L[0]][0:2]])
            a1=a0/2
            thershold=np.abs(a0-a1*(math.exp(-1*coef[1])))
            thershold_code=10
            buff[np.abs(buff)>thershold]=np.nan
            buff[np.abs(buff)<=thershold]=0
            buff[np.isnan(buff)]=1
            buff2_ctl[np.abs(buff2_ctl)>10]=np.nan
            buff2_ctl[np.abs(buff2_ctl)<=10]=1
            buff_ctl[np.abs(buff_ctl)>thershold]=np.nan
            buff_ctl[np.abs(buff_ctl)<=thershold]=1
            buff_ctl_code[np.abs(buff_ctl_code)>thershold_code]=np.nan
            buff_ctl_code[np.abs(buff_ctl_code)<=thershold_code]=1                
            buff0=self.ctl1[0][0,:].copy();buff0[np.isnan(buff0)]=0
            buff1=self.ctl1[0][1,:].copy();buff1[np.isnan(buff1)]=0
            buff00=np.row_stack([buff0,buff1])
            buff=np.row_stack([buff00,buff[0:-1,:]])
            N.append(buff)
            buff_ctl=np.row_stack([self.ctl1[0][0,:].copy(),buff_ctl])
            buff_ctl*=buff2_ctl.copy()
            buff_ctl_code=np.row_stack([self.ctl1[1][0,:].copy(),buff_ctl_code])
            CTL2.append([buff_ctl,buff_ctl_code])
        elif(user_proc_mode!='code-only' and self.user_CombType=='single-freq'):
                buff_ctl=np.array([i[:,:] for i in self.obs])
                buff_ctl_phase=C/Freq[self.freq][user_Freq[self.freq][self.L[0]][0:2]]*buff_ctl[self.L[0],:,:]-buff_ctl[self.P[0],:,:]
                buff2_ctl=np.array([i[:,:] for i in self.obs])
                buff2_ctl=C/Freq[self.freq][user_Freq[self.freq][self.L[0]][0:2]]*buff2_ctl[self.L[0],:,:]-buff2_ctl[self.P[0],:,:]

                buff=np.array([i[:-1,:]-i[1:,:] for i in self.obs])

                buff=C/Freq[self.freq][user_Freq[self.freq][self.L[0]][0:2]]*buff[self.L[0],:,:]-buff[self.P[0],:,:]
                buff_ctl=buff.copy()
                thershold=4
                buff[np.abs(buff)>thershold]=np.nan
                buff[np.abs(buff)<=thershold]=0
                buff[np.isnan(buff)]=1
                buff2_ctl[np.abs(buff2_ctl)>10]=np.nan
                buff2_ctl[np.abs(buff2_ctl)<=10]=1
                buff_ctl[np.abs(buff_ctl)>thershold]=np.nan
                buff_ctl[np.abs(buff_ctl)<=thershold]=1
                buff0=self.ctl1[0][0,:].copy();buff0[np.isnan(buff0)]=0
                buff1=self.ctl1[0][1,:].copy();buff1[np.isnan(buff1)]=0
                buff00=np.row_stack([buff0,buff1])
                buff=np.row_stack([buff00,buff[0:-1,:]])
                N.append(buff)
                buff_ctl=np.row_stack([self.ctl1[0][0,:].copy(),buff_ctl])
                buff_ctl*=buff2_ctl.copy()   
                CTL2.append([self.ctl1[0],self.ctl1[1]])
        else:
            CTL2.append(self.ctl1)

    def ObsMatrix(self,L_Matrix,user_CombType,user_proc_mode,MF_iono,vTEC):
        
        if(user_CombType=='single-freq'):
            # print(vTEC)
            L_P=np.array([self.Landa[0]*self.obs[self.L[0]]]).T
            L_P[:,:,0]-=(MF_iono*np.array([vTEC[:,0]]))
            L_C=np.array([self.obs[self.P[0]]]).T
            L_C[:,:,0]-=(MF_iono*np.array([vTEC[:,1]]))
        elif(user_CombType=='iono-free'):
            #L=[P,C]
            L_P=np.array([((self.alpha*self.Landa[0]*self.obs[self.L[0]])-(self.beta*self.Landa[1]*self.obs[self.L[1]]))]).T
            L_C=np.array([((self.alpha*self.obs[self.P[0]])-(self.beta*self.obs[self.P[1]]))]).T


        if(user_proc_mode=='combined'):
            L_Matrix.append([L_P[:,:,0],L_C[:,:,0]])
        elif(user_proc_mode=='code-only'):
            L_Matrix.append([L_C[:,:,0]])
        


        
        self.L_eq,self.X_eq,self.Y_eq,self.Z_eq,self.MF_eq,self.MF_gradE_eq,self.MF_gradN_eq=[],[],[],[],[],[],[]
        self.A_eq,self.N_eq,self.CLK_eq,self.time_eq=[],[],[],[]
        self.W_L_eq,self.Vx_eq,self.Vy_eq,self.Vz_eq,self.N2=[],[],[],[],[]
        self.sat_nums,self.sat_name,self.MF_zhd_eq=[],[],[]
        self.test0,self.test1,self.ind_sat_eq=[],[],[]

    def A_Matrix(self,L_Matrix,ctl,Xs,Ys,Zs,Vx_s,Vy_s,Vz_s,MF_t,N,Elevation,OBS,OBS_TIME,clk,user_TimeInterval,user_Sats,user_proc_mode):
        user_move='static'
        self.L_eq=[];self.X_eq=[];self.Y_eq=[];self.Z_eq=[];self.MF_eq=[];self.MF_gradE_eq=[];self.MF_gradN_eq=[]
        self.A_eq=[];self.N_eq=[];
        self.CLK_eq=[];self.time_eq=[]
        self.W_L_eq=[];self.Vx_eq=[];self.Vy_eq=[];self.Vz_eq=[];self.N2=[]
        self.sat_nums=[];self.sat_name=[];self.MF_zhd_eq=[]
        self.test0=[];self.test1=[]
        self.ind_sat_eq=[]
        interv=int(user_TimeInterval/self.obsRate)
        # print('interv: '+str(interv))
        if(user_move=='static'):
            for ep in range(0,L_Matrix[0][0].shape[1]-1,interv):
                self.L_eq.append([]);self.X_eq.append([]);self.Y_eq.append([]);self.Z_eq.append([]);self.MF_eq.append([]);self.MF_gradE_eq.append([]);self.MF_gradN_eq.append([])
                self.A_eq.append([]);self.N_eq.append([]);
                self.CLK_eq.append([]);self.time_eq.append([])
                self.W_L_eq.append([]);self.Vx_eq.append([])
                self.Vy_eq.append([]);self.Vz_eq.append([]);self.N2.append([]);self.sat_nums.append([]);self.sat_name.append([])
                self.MF_zhd_eq.append([])
                self.test0.append([]);self.test1.append([])
                self.ind_sat_eq.append([])
                for s in range(len(user_Sats)):
                    if(s>0):buff_s=len(self.N_eq[-1])
                    else:buff_s=0
                    shap=L_Matrix[s][0].shape[0]*interv
                    if(user_proc_mode=='combined' and L_Matrix[s][0][:,ep:ep+interv].shape[0]*L_Matrix[s][0][:,ep:ep+interv].shape[1]==shap):
                        #phase
                        self.sat_nums[-1].append(np.unique(np.where(~np.isnan(ctl[s][0][ep:ep+interv,:]))[1]).shape[0])
                        self.sat_name[-1].append([OBS[s].sv.data[k] for k in np.unique(np.where(~np.isnan(ctl[s][0][ep:ep+interv,:]))[1])])
                        temp_test=N[s][ep:ep+interv,:]*ctl[s][0][ep:ep+interv,:]
                        temp_test2=temp_test.copy();temp_test2[~np.isnan(temp_test2)]=1
                        for i in range(temp_test.shape[1]):
                            if(1 in list(temp_test[:,i])):temp_test[:,i]=1
                        temp_test=temp_test*temp_test2
                        if(s==0):
                            self.test0[-1].append(temp_test)
                            self.ind_sat_eq[-1].append(np.reshape(ctl[s][0][ep:ep+interv,:],[shap,1]))
                        else:
                            self.test0[-1].append(temp_test.copy())
                            self.ind_sat_eq[-1][0]=np.row_stack([self.ind_sat_eq[-1][0],np.reshape(ctl[s][0][ep:ep+interv,:],[shap,1])])
                        temp=N[s][ep:ep+interv,:]*ctl[s][0][ep:ep+interv,:];
                        temp2=temp.copy();temp2[~np.isnan(temp2)]=1
                        for i in range(temp.shape[1]):
                            if(1 in temp[:,i]):temp[:,i]=1
                        temp=temp*temp2
                        temp[temp==-1]=1
                        [self.N2[-1].append(enum) for enum in temp]
                        [self.N_eq[-1].append(enum+buff_s) for enum in temp[~np.isnan(temp)]]
                        temp=np.reshape(L_Matrix[s][0][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
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
                        temp=np.reshape(MF_t[s][2][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.MF_gradE_eq[-1].append(enum) for enum in temp]
                        temp=np.reshape(MF_t[s][3][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.MF_gradN_eq[-1].append(enum) for enum in temp]
                        temp=np.reshape(Elevation[s][:,ep:ep+interv].T*ctl[s][0][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power(0.03,2)) for enum in temp]
                        

                        #code
                        if(s==0):
                            self.ind_sat_eq[-1].append(np.reshape(ctl[s][1][ep:ep+interv,:],[shap,1]))
                        else:
                            self.ind_sat_eq[-1][1]=np.row_stack([self.ind_sat_eq[-1][1],np.reshape(ctl[s][1][ep:ep+interv,:],[shap,1])])
                        self.sat_nums[-1].append(np.unique(np.where(~np.isnan(ctl[s][1][ep:ep+interv,:]))[1]).shape[0])
                        self.sat_name[-1].append([OBS[s].sv.data[k] for k in np.unique(np.where(~np.isnan(ctl[s][1][ep:ep+interv,:]))[1])])
                        temp=np.reshape(L_Matrix[s][1][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.L_eq[-1].append(enum) for enum in temp]
                        temp=np.reshape(clk[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.CLK_eq[-1].append(enum) for enum in temp]
                        [self.N_eq[-1].append(0) for enum in range(len(temp))]
                        temp=np.reshape(Xs[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.X_eq[-1].append(enum) for enum in temp]
                        self.test1[-1].append(temp)
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
                        temp=np.reshape(MF_t[s][2][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.MF_gradE_eq[-1].append(enum) for enum in temp]
                        temp=np.reshape(MF_t[s][3][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.MF_gradN_eq[-1].append(enum) for enum in temp]
                        temp=np.reshape(Elevation[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.W_L_eq[-1].append(np.power(3,2)) for enum in temp]


                    if(user_proc_mode=='code-only' and L_Matrix[s][0][:,ep:ep+interv].shape[0]*L_Matrix[s][0][:,ep:ep+interv].shape[1]==shap):
                        if(s==0):
                            self.ind_sat_eq[-1].append(np.reshape(ctl[s][1][ep:ep+interv,:],[shap,1]))
                        else:
                            self.ind_sat_eq[-1][0]=np.row_stack([self.ind_sat_eq[-1][0],np.reshape(ctl[s][1][ep:ep+interv,:],[shap,1])])

                        self.sat_nums[-1].append(np.unique(np.where(~np.isnan(ctl[s][1][ep:ep+interv,:]))[1]).shape[0])
                        self.sat_name[-1].append([OBS[s].sv.data[k] for k in np.unique(np.where(~np.isnan(ctl[s][1][ep:ep+interv,:]))[1])])
                        temp=np.reshape(L_Matrix[s][0][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
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
                        temp=np.reshape(MF_t[s][2][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.MF_gradE_eq[-1].append(enum) for enum in temp]
                        temp=np.reshape(MF_t[s][3][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]
                        [self.MF_gradN_eq[-1].append(enum) for enum in temp]
                        temp=np.reshape(Elevation[s][:,ep:ep+interv].T*ctl[s][1][ep:ep+interv,:],[shap,1]);temp=temp[~np.isnan(temp)]

                        [self.W_L_eq[-1].append(9) for enum in temp]

                self.A_eq[-1]=np.array([self.MF_eq[-1]]).T
                self.A_eq[-1]=np.column_stack([self.A_eq[-1],np.ones([self.A_eq[-1].shape[0],1])])

                if(user_proc_mode!='code-only'):
                    temp0=self.test0[-1][0].copy()
                    buff1=temp0[:,np.any(~np.isnan(temp0), axis=0)]
                    buff1[buff1==0]=1
                    temp0=np.diag(list(buff1[0,:]))
                    for enum0 in range(1,buff1.shape[0]):temp0=np.row_stack([temp0,np.diag(list(buff1[enum0,:]))])
                    for enum in range(1,len(user_Sats)):
                        temp0_size=temp0.shape[1]
                        buffer=self.test0[-1][enum].copy()
                        buff1=buffer[:,np.any(~np.isnan(buffer), axis=0)]
                        temp0=np.column_stack([temp0,np.zeros([temp0.shape[0],buff1.shape[1]])])
                        temp00=np.diag(list(buff1[0,:]))
                        for enum0 in range(1,buff1.shape[0]):temp00=np.row_stack([temp00,np.diag(list(buff1[enum0,:]))])
                        temp00=np.column_stack([np.zeros([temp00.shape[0],temp0_size]),temp00])
                        temp0=np.row_stack([temp0,temp00])
                    
                    temp0=temp0[np.all(~np.isnan(temp0), axis=1)]
                    buff2=np.row_stack([temp0,np.zeros([self.A_eq[-1].shape[0]-temp0.shape[0],temp0.shape[1]])])
                    if(buff2.shape[1]!=0):
                        self.A_eq[-1]=np.column_stack([self.A_eq[-1],buff2])

    def Noise(self,mode,ind_sat,L,A,W_L,N_buff,x0,L00,W_X,user_TimeInterval,user_GDOP):
        user_OutlierThreshold=2+np.abs(self.clk_init)
        # user_OutlierThreshold=10

        # global_test=np.sqrt(np.matmul(L.T,np.matmul(np.linalg.inv(W_L),L)))
        # global_test+=np.abs(self.clk_init)
        # print(global_test)
        # print('*****')
        
        # print(L)
        # print('###############')
        if(mode=='init' and self.W_X_N==[]):
            thershold=1e7
        elif(np.sum(self.W_X_N)==0 and self.W_X_N!=[]):
            thershold=1e6
        else:
            thershold=1e6

        W_L=list(np.diag(W_L))
        N_buff2=[];[[N_buff2.append(j) for j in i] for i in N_buff]
        ech=user_TimeInterval//self.obsRate
        temp3=[]
        j=0;A_buff=[];L_buff=[];W_L_buff=[];L00_buff=[]
        phi=self.PHI;qx=self.Qx
        for kk in range(len(ind_sat)):
            temp2=[]
            for jj in range(0,ech):
                temp=[]
                for i in range(0,len(ind_sat[kk])//ech):
                    if(np.isnan(ind_sat[kk][i+(jj*(len(ind_sat[kk])//ech))][0])==False):
                        if(np.abs(L[j])<thershold):
                        # if(global_test[0][0]<thershold):
                            temp.append(L[j][0])
                            L_buff.append(L[j][0])
                            L00_buff.append(L00[j][0])
                            A_buff.append(list(A[j,:]))
                            W_L_buff.append(W_L[j])
                            if(len(ind_sat)>1 and kk==0):
                                if(N_buff2[i]==3 or N_buff2[i]==0):N_buff2[i]=1000
                                elif(N_buff2[i]==4 or N_buff2[i]==1):N_buff2[i]=1001
                        else:
                            temp.append(np.nan)
                            if(len(ind_sat)>1 and kk==0):
                                if(N_buff2[i]==0):N_buff2[i]=3
                                elif(N_buff2[i]==1):N_buff2[i]=4
                        j+=1
                    else:
                        temp.append(np.nan)
                 
                temp2.append(temp)
            temp3.append(temp2)
        for i in range(len(N_buff2)):
            if(N_buff2[i]==3 or N_buff2[i]==4):N_buff2[i]=np.nan
            if(N_buff2[i]==1000):N_buff2[i]=0
            if(N_buff2[i]==1001):N_buff2[i]=1

        X0=x0[0:5]
        if(len(ind_sat)>1):
            enum3,enum4=[],[]
            enum5=0
            for enum0 in range(len(N_buff2)):
                if(np.isnan(N_buff2[enum0])==False):
                    if(N_buff2[enum0]==1):
                        enum4.append(0)
                        enum3.append(400)
                    elif(N_buff2[enum0]==0 and self.W_X_N[0,enum0]!=0):
                        enum4.append(self.N[0,enum0])
                        enum3.append(self.W_X_N[0,enum0])
                        enum5+=1
                    elif(N_buff2[enum0]==0 and self.W_X_N[0,enum0]==0):
                        enum4.append(0)
                        enum3.append(400)
                    else:
                        enum4.append(0)
                        enum3.append(400)
                else:
                    self.W_X_N[0,enum0]=0
                    self.N[0,enum0]=0
            W_X1=np.zeros([W_X.shape[0]+len(enum3),W_X.shape[1]+len(enum3)])
            X0=x0[0:5]
            X0=np.array([np.append(X0,enum4)]).T
            phi=list(np.append(phi,[1 for i in range(len(enum3))]))
            qx=list(np.append(qx,[0 for i in range(len(enum3))]))
            W_X1[0:5,0:5]=W_X.copy()
            if(enum3!=[]):np.fill_diagonal(W_X1[5:,5:],enum3)
            W_X=W_X1.copy()

        A_buff=np.array(A_buff)
        A_buff=A_buff[:,~np.all(A_buff==0,axis=0)]
        if(A_buff.shape[0]<A_buff.shape[1] or len(A_buff)==0):
            return [],[],[],[],[],[],[],[],[]

        phi=np.diag(phi)
        qx=np.diag(qx)
        x0=X0.copy()
        L_buff=np.array([L_buff]).T

        L00_buff=np.array([L00_buff]).T
        pre_x=np.matmul(phi,x0)
        pre_px=np.matmul(np.matmul(phi,W_X),phi.T)+qx
        L_buff=L00_buff-np.matmul(A_buff,pre_x)
        A0=np.column_stack([A_buff[:,0:3],A_buff[:,4]])
        buffer=np.diag(np.linalg.inv(np.matmul(A0.T,A0)))
        GDOP=np.sqrt(buffer[0]+buffer[1]+buffer[2]+buffer[3])
        if(GDOP>user_GDOP):
            return [],[],[],[],[],[],[],[],[]
        self.GDOP.append(GDOP)
        if(len(A_buff)!=0):
            if(A_buff.shape[0]>=A_buff.shape[1]):
                self.inovation.append(temp3)
        return L_buff,L00_buff,A_buff,np.diag(W_L_buff),W_X,x0,pre_x,pre_px,N_buff2

    def LeastSquare(self,L,A,W_L):
        err=1e12
        df=A.shape[0]-A.shape[1]
        Xcap=np.zeros([A.shape[1],1])
        while(err>1e-5):
            buffer=np.matmul(np.linalg.inv(np.matmul(np.matmul(A.T,W_L),A)),(np.matmul(np.matmul(A.T,W_L),L)))
            err=np.abs(Xcap[0,0]-buffer[0,0])
            Xcap+=buffer.copy()
            sigma=np.matmul(np.matmul(L.T,W_L),L)/df
            W_X=np.linalg.inv(np.matmul(np.matmul(A.T,W_L),A))*sigma
            # W_X=np.matmul(np.matmul(A.T,W_L),A)
        return Xcap,W_X

    def KalmanFilter(self,X0,L,A,Xs,Ys,Zs,Vx,Vy,Vz,N2,W_L,mode,init,MF_zhd,MF_wet,clk,ind_sat,test0,user_TimeInterval,user_Sats,user_GDOP):
        user_move='static'
        prefit,postfit=[],[]        
        phi=self.PHI;qx=self.Qx
        enum3,enum4=[],[]
        enum5=0
        if(init=='normal'):W_X=self.W_X[-1][0:5,0:5].copy()
        elif(init=='init'):W_X=np.diag([1.00e+008,1.00e+008,1.00e+008,0.25,9.00e+010])
        N_buff2=[]
        if(mode!='code-only'):
            for enum0 in range(len(user_Sats)):
                N_buff=test0[enum0][0].copy()
                N_buff[np.isnan(N_buff)]=2
                for enum_N in range(0,len(test0[enum0])):
                    buff=test0[enum0][enum_N].copy()
                    buff[np.isnan(buff)]=1000
                    N_buff*=buff;N_buff[N_buff==2]=1;N_buff[N_buff==1000]=1;N_buff[N_buff==2000]=2
                N_buff[N_buff==2]=np.nan
                N_buff2.append(N_buff)
                for enum1 in range(len(N_buff)):
                    if(np.isnan(N_buff[enum1])==False):
                        if(N_buff[enum1]==1):
                            enum4.append(0)
                            enum3.append(400)
                        elif(N_buff[enum1]==0 and self.W_X_N[0,enum1+(len(N_buff)*enum0)]!=0):
                            enum4.append(self.N[0,enum1+(len(N_buff)*enum0)])
                            enum3.append(self.W_X_N[0,enum1+(len(N_buff)*enum0)])
                            enum5+=1
                        elif(N_buff[enum1]==0 and self.W_X_N[0,enum1+(len(N_buff)*enum0)]==0):
                            enum4.append(0)
                            enum3.append(400)

            W_X1=np.zeros([W_X.shape[0]+len(enum3),W_X.shape[1]+len(enum3)])
            X0=X0[0:5]
            X0=np.array([np.append(X0,enum4)]).T
            phi=list(np.append(phi,[1 for i in range(len(enum3))]))
            qx=list(np.append(qx,[0 for i in range(len(enum3))]))
            W_X1[0:5,0:5]=W_X.copy()
            if(enum3!=[]):np.fill_diagonal(W_X1[5:,5:],enum3)
            W_X=W_X1.copy()
        phi=np.diag(phi)
        qx=np.diag(qx)
        x0=X0.copy()
        A1=A.copy()
        geodetic=pymap3d.ecef2geodetic(x0[0],x0[1],x0[2])
        self.ZHD.append((2.3)*np.exp(-0.116*math.pow(10,-3)*(geodetic[2])))
        ru=np.sqrt(np.power(x0[0]-Xs,2)+np.power(x0[1]-Ys,2)+np.power(x0[2]-Zs,2))
        prefit.append(L.copy())
        prefit.append(ru.copy())
        if(user_move=='static'):
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
        # print('sat: '+str(Xs))
        # print('pre_ru: '+str(L.copy()-ru.copy()))
        # print('L: '+str(L.copy()))
        ru=np.sqrt(np.power(x0[0]-Xs,2)+np.power(x0[1]-Ys,2)+np.power(x0[2]-Zs,2))-clk+self.delta_rel[-1]+self.delta_rel_e[-1]+self.earth_rot[-1]+(self.ZHD[-1]*MF_zhd)
        # print('next_ru: '+str(L.copy()-ru.copy()))
        prefit.append(ru.copy())
        L1=L.copy()
        L1=L1.copy()-ru
        L00=L1.copy()

        if(init=='init'):
            self.W_X.append(W_X)
        x0[0:3]-=self.processed[0][0:3]
        pre_x=np.matmul(phi,x0)
        pre_px=np.matmul(np.matmul(phi,W_X),phi.T)+qx
        L1=L1.copy()-np.matmul(A1,pre_x)
        L2,L002,A2,W_L2,W_X,x0,pre_x,pre_px,N_buff2=self.Noise(init,ind_sat,L1,A1,W_L,N_buff2,x0[0:5].copy(),L00,W_X[0:5,0:5],user_TimeInterval,user_GDOP)
        # print(A1.shape)
        # print(A2.shape)
        if(L2==[]):
            return 'FDE_not_accepted'
        if(len(A2)==0):
            return 'FDE_not_accepted'        
        if(A2.shape[0]<A2.shape[1]):
            return 'FDE_not_accepted'
        K=np.matmul(np.matmul(pre_px,A2.T),np.linalg.pinv(np.matmul(np.matmul(A2,pre_px),A2.T)+W_L2))
        Xcap=pre_x+np.matmul(K,L2.copy())
        self.clk_init=Xcap[4].copy()
        self.residual.append(np.matmul(K,L2.copy()))
        prefit.append(L2.copy())
        postfit.append(L2.copy()+np.matmul(A2,pre_x)-np.matmul(A2,Xcap))
        temp=[];temp2=[];j=0;SS=L002.copy()-np.matmul(A2,Xcap)
        for enum0 in range(len(self.inovation[-1])):
            temp2=[]
            for enum1 in range(len(self.inovation[-1][enum0][0])):
                if(np.isnan(self.inovation[-1][enum0][0][enum1])==False):
                    temp2.append(SS[j][0])
                    j+=1
                else:
                    temp2.append(np.nan)
            temp.append(temp2)

        self.POSTFIT.append(temp)

        post_px=np.matmul((np.eye(A2.shape[1])-np.matmul(K,A2)),pre_px)
        post_px=np.diag(list(np.diag(post_px)))
        # print('xvap: '+str(Xcap.copy()))
        # print('#'*10)
        Xcap[0:3]+=self.processed[0][0:3]

        if(A2.shape[0]>A2.shape[1]):
            if(mode=='code-only'):
                self.POSTFIT.append(postfit)
            self.PREFIT.append(prefit)
            if(mode!='code-only'):
                buff_Nindex=list(Xcap[5:][:,0])
                postfit.append(buff_Nindex.copy())
                buff_W_N_index=list(np.diag(post_px[5:,5:]))
                jj=0   
                for enum1 in range(len(N_buff2)):
                    if(np.isnan(N_buff2[enum1])==False):
                        self.N[0,enum1]=buff_Nindex[jj]
                        self.W_X_N[0,enum1]=buff_W_N_index[jj]
                        jj+=1

            self.processed.append(Xcap)
            self.W_X.append(post_px)
            self.Ns.append(self.N.copy())
            return 'FDE_accepted'
        else:
            return 'FDE_not_accepted'    

    def exe(self,user_proc_mode,user_TimeInterval,user_Sats,user_GDOP,SP3_TIMESEC,OBS_TIMESEC,OBS_RATE):
        self.processed.append(np.array([[self.Obs.position[0],self.Obs.position[1],self.Obs.position[2],0,0]]).T)
        if(user_proc_mode!='code-only'):
            self.N,self.W_X_N=np.zeros([1,np.sum([i.shape[1] for i in self.test0[0]])]),np.zeros([1,np.sum([i.shape[1] for i in self.test0[0]])])
        else:
            self.N,self.W_X_N=[],[]
        enum2=0
        if(SP3_TIMESEC[0][-1]<OBS_TIMESEC[0][-1]):
            lens=len(self.X_eq)-int((OBS_TIMESEC[0][-1]-SP3_TIMESEC[0][-1])/OBS_RATE[0])
        else:
            lens=len(self.X_eq)        
        for enum in range(0,lens):
            x0=self.processed[enum2]
            Xs=np.array([self.X_eq[enum]]).T
            Ys=np.array([self.Y_eq[enum]]).T
            Zs=np.array([self.Z_eq[enum]]).T
            Vx=np.array([self.Vx_eq[enum]]).T
            Vy=np.array([self.Vy_eq[enum]]).T
            Vz=np.array([self.Vz_eq[enum]]).T
            MF_zhd=np.array([self.MF_zhd_eq[enum]]).T
            MF_wet=np.array([self.MF_eq[enum]]).T
            clk=np.array([self.CLK_eq[enum]]).T*C
            test0=self.test0[enum].copy()
            if(user_proc_mode=='code-only'):ind_sat=[self.ind_sat_eq[enum][0]]
            else:ind_sat=[self.ind_sat_eq[enum][0],self.ind_sat_eq[enum][1]]

            L=np.array([self.L_eq[enum]]).T
            A=self.A_eq[enum]
            N2=[]
            if(user_proc_mode=='combined'):N2=self.N2[enum]
            W_L=np.diag(self.W_L_eq[enum])
            if(L.shape[0]>=(3+A.shape[1])):
                # print('ep: '+str(enum))
                init='normal'
                if(enum2==0):
                    init='init'
                
                res=eqRel.KalmanFilter(self,x0,L,A,Xs,Ys,Zs,Vx,Vy,Vz,N2,W_L,user_proc_mode,init,MF_zhd,MF_wet,clk,ind_sat,test0,user_TimeInterval,user_Sats,user_GDOP)
                if(res=='FDE_accepted'):
                    self.EPOCH.append(enum);self.sat_num.append(self.sat_nums[enum]);self.sats_name.append([])
                    temp=[[self.sats_name[-1].append(j) for j in i] for i in self.sat_name[enum]]

                    self.out.append([self.X_eq[-1],self.Y_eq[-1],self.Z_eq[-1],
                                     self.sat_mov[-1][0],self.sat_mov[-1][1],self.sat_mov[-1][2],
                                     self.Vx_eq[-1],self.Vy_eq[-1],self.Vz_eq[-1],self.earth_rot[-1],
                                     self.delta_rel[-1],self.delta_rel_e[-1]])

                    enum2+=1
                elif(res=='FDE_not_accepted' and self.W_X_N!=[]):
                    self.N[self.N!=0]=0
                    self.W_X_N[self.W_X_N!=0]=0

            else:
                if(self.W_X_N!=[]):
                    self.N[self.N!=0]=0
                    self.W_X_N[self.W_X_N!=0]=0
