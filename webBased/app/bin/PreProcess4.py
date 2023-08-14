from app.bin.Libs import *




class preprocess():
    def __init__(self,obs,sp3,sat):
        self.obs=obs
        # print('self.obs: '+str(self.obs.C1C.data.shape))
        self.sp3=sp3
        self.sat=sat
        self.sats_obs=list(obs.sv.data)

    def TimeProc(self,OBS_TIMESEC,SP3_TIMESEC,OBS_RATE,SP3_RATE,DATATIMEGFS,DATATIMEECMWF):
        

        buff1=self.obs.time.data[0:self.obs.time.data.shape[0]]
        # print('self.obs.time.data.shape[0]: '+str(self.obs.time.data.shape))
        # print('buff1: '+str(len(buff1)))

        buff2=self.sp3.time.data[0:self.sp3.time.data.shape[0]]
        if(str(self.obs.version)[0]=='2'):
            OBS_TIMESEC.append([int(enum0.timestamp()) for enum0 in self.obs.time.data])
            
            # OBS_TIMESEC.append([int(enum0.strftime("%s")) for enum0 in self.obs.time.data])
            # SP3_TIMESEC.append([int(enum0.strftime("%s")) for enum0 in self.sp3.time.data])
            # SP3_TIMESEC.append([int(datetime.timestamp(datetime.strptime(str(enum0).split('.')[0],"%Y-%m-%dT%H:%M:%S"))) for enum0 in buff2])
            SP3_TIMESEC.append([int(str(enum0.item())[0:-9]) for enum0 in self.sp3.time.data]) 
        elif(str(self.obs.version)[0]=='3'):
            # OBS_TIMESEC.append([np.round(datetime.timestamp(datetime.strptime(str(enum0).split('.')[0],"%Y-%m-%dT%H:%M:%S"))) for enum0 in buff1])
            # SP3_TIMESEC.append([np.round(datetime.timestamp(datetime.strptime(str(enum0).split('.')[0],"%Y-%m-%dT%H:%M:%S"))) for enum0 in buff2])
            OBS_TIMESEC.append([int(np.round(datetime.timestamp(datetime.strptime(str(enum0)[:-3],"%Y-%m-%dT%H:%M:%S.%f")))) for enum0 in buff1])
            SP3_TIMESEC.append([int(np.round(datetime.timestamp(datetime.strptime(str(enum0)[:-3],"%Y-%m-%dT%H:%M:%S.%f")))) for enum0 in buff2])

        # SP3_TIMESEC.append([int(str(enum0.item())[0:-9]) for enum0 in self.sp3.time.data])
        # OBS_RATE.append(np.round(OBS_TIMESEC[-1][1]-OBS_TIMESEC[-1][0]))
        OBS_RATE.append(int(self.obs.interval))
        SP3_RATE.append(np.round(SP3_TIMESEC[-1][1]-SP3_TIMESEC[-1][0]))
        ShiftedNow=datetime.now()-timedelta(days=0, seconds=0, microseconds=0, milliseconds=0, minutes=0, hours=9, weeks=0)
        GFSTime=ShiftedNow.strftime("%Y%m%d")
        series=[0,6,12,18];temp=[abs(ShiftedNow.hour-enum) for enum in series if ShiftedNow.hour>=enum];temp=series[temp.index(min(temp))]
        if(len(str(temp))==1):temp='0'+str(temp)
        else:temp=str(temp)
        DATATIMEGFS.append([GFSTime,temp,OBS_TIMESEC[-1][0]])
        if(str(self.obs.version)[0]=='2'):
            DATATIMEECMWF.append([buff1[0].strftime("%Y"),buff1[0].strftime("%m"),buff1[0].strftime("%d"),buff1[0]>=(datetime.now()-timedelta(days=7, seconds=0, microseconds=0, milliseconds=0, minutes=0, hours=0, weeks=0))])
        elif(str(self.obs.version)[0]=='3'):
            buff=datetime.strptime(str(buff1[0]).split('.')[0],"%Y-%m-%dT%H:%M:%S")
            DATATIMEECMWF.append([buff.strftime("%Y"),buff.strftime("%m"),buff.strftime("%d"),buff>=(datetime.now()-timedelta(days=7, seconds=0, microseconds=0, milliseconds=0, minutes=0, hours=0, weeks=0))])

        
    def sp3Extract(self,X,Y,Z,Vx,Vy,Vz,CLK):
        buff1=[];buff2=[];buff3=[];buff4=[]
        clk_ed=self.sp3.clock.data
        clk_ed[clk_ed==999999.999999]=np.nan
        for i in self.sats_obs:
            if(self.sat in i):
                try:buff1.append(self.sp3.position[:,list(self.sp3.sv.data).index(i),0].data*math.pow(10,3))
                except:buff1.append(np.nan*np.ones([self.sp3.position.data.shape[0]]))
                try:buff2.append(self.sp3.position[:,list(self.sp3.sv.data).index(i),1].data*math.pow(10,3))
                except:buff2.append(np.nan*np.ones([self.sp3.position.data.shape[0]]))
                try:buff3.append(self.sp3.position[:,list(self.sp3.sv.data).index(i),2].data*math.pow(10,3))
                except:buff3.append(np.nan*np.ones([self.sp3.position.data.shape[0]]))
                try:buff4.append(clk_ed[:,list(self.sp3.sv.data).index(i)]*math.pow(10,-6))
                except:buff4.append(np.nan*np.ones([self.sp3.clock.data.shape[0]]))
        X.append(np.array(buff1));Y.append(np.array(buff2));Z.append(np.array(buff3));CLK.append(np.array(buff4))

    def Interpolation(self,items,vals,res,dim,Kind='linear'):
        if(dim=='1'):
            f=interpolate.interp1d(items,vals,kind=Kind)
            return f(res)
        elif(dim=='2'):
            f=interpolate.interp2d(items[0],items[1],vals,kind='cubic')
            return f(res[0],res[1])

        elif(dim=='3'):
            f=interpolate.splrep(items,vals,k=5)
            return interpolate.splev(res,f)

        elif(dim=='4'):
            f=interpolate.lagrange(items,vals)
            return f(res)


    def proc1(self,x,y,z,clk,obsSec,sp3Sec,obsRate,sp3Rate,X_PROC,Y_PROC,Z_PROC,CLK_PROC,Vx_s,Vy_s,Vz_s,OBS_TIMESEC_PROC): 

        user_Sp3InterpolationDegree=10
        user_ClkInterpolationDegree=2

        scale=(obsSec[0]-sp3Sec[0])
        scale2=int(np.floor(scale/sp3Rate))
        # scale=scale2*sp3Rate
        TimeNewObs=list(range(0,len(obsSec)*obsRate,obsRate))
        # print('obsSec: '+str(len(obsSec)))
        # print('TimeNewObs: '+str(len(TimeNewObs)))
        TimeNewSp3=list(range(0,len(sp3Sec)*sp3Rate,sp3Rate))
        # print('TimeNewSp3'+str(TimeNewSp3))
        OBS_TIMESEC_PROC.append(TimeNewObs)

        # print(scale)
        # print('scale: '+str(scale))
        # print('scale2: '+str(scale2))
        temp0=[[] for enum in range(x.shape[0])]
        temp1=[[] for enum in range(x.shape[0])]
        temp2=[[] for enum in range(x.shape[0])]
        temp3=[[] for enum in range(x.shape[0])]
        for i in range(scale2,len(TimeNewSp3)):
            # print('epoch:' +str(i))
            if(i<int(user_Sp3InterpolationDegree/2)+1):
                for j in range(x.shape[0]):

                    f=np.poly1d(np.polyfit(TimeNewSp3[0:user_Sp3InterpolationDegree+1],x[j,0:user_Sp3InterpolationDegree+1],user_Sp3InterpolationDegree))
                    
                    temp=[temp0[j].append(f(TimeNewObs[enum]+scale)) for enum in range(len(TimeNewObs)) if(obsSec[enum]>=sp3Sec[i] and obsSec[enum]<sp3Sec[i+1])]

                    f=np.poly1d(np.polyfit(TimeNewSp3[0:user_Sp3InterpolationDegree+1],y[j,0:user_Sp3InterpolationDegree+1],user_Sp3InterpolationDegree))
                    temp=[temp1[j].append(f(TimeNewObs[enum]+scale)) for enum in range(len(TimeNewObs)) if(obsSec[enum]>=sp3Sec[i] and obsSec[enum]<sp3Sec[i+1])]

                    f=np.poly1d(np.polyfit(TimeNewSp3[0:user_Sp3InterpolationDegree+1],z[j,0:user_Sp3InterpolationDegree+1],user_Sp3InterpolationDegree))
                    temp=[temp2[j].append(f(TimeNewObs[enum]+scale)) for enum in range(len(TimeNewObs)) if(obsSec[enum]>=sp3Sec[i] and obsSec[enum]<sp3Sec[i+1])]

            else:
                

                for j in range(x.shape[0]):
                    try:
                        timee=TimeNewSp3[i-int((user_Sp3InterpolationDegree)/2):i+int(user_Sp3InterpolationDegree/2)+1]
                        # print(timee)
                        timee=[k-TimeNewSp3[i] for k in timee]
                        # print(timee)
                        # print([(TimeNewObs[enum]-TimeNewSp3[i]+scale) for enum in range(len(TimeNewObs)) if((obsSec[enum])>=sp3Sec[i] and (obsSec[enum])<sp3Sec[i+1])])
                        f=np.poly1d(np.polyfit(timee,x[j,i-int((user_Sp3InterpolationDegree)/2):i+int(user_Sp3InterpolationDegree/2)+1],len(timee)-1))  
                        temp=[temp0[j].append(f((TimeNewObs[enum]-TimeNewSp3[i]+scale))) for enum in range(len(TimeNewObs)) if((obsSec[enum])>=sp3Sec[i] and (obsSec[enum])<sp3Sec[i+1])]

                        f=np.poly1d(np.polyfit(timee,y[j,i-int((user_Sp3InterpolationDegree)/2):i+int(user_Sp3InterpolationDegree/2)+1],len(timee)-1))  
                        temp=[temp1[j].append(f((TimeNewObs[enum]-TimeNewSp3[i]+scale))) for enum in range(len(TimeNewObs)) if(obsSec[enum]>=sp3Sec[i] and obsSec[enum]<sp3Sec[i+1])]

                        f=np.poly1d(np.polyfit(timee,z[j,i-int((user_Sp3InterpolationDegree)/2):i+int(user_Sp3InterpolationDegree/2)+1],len(timee)-1))   
                        temp=[temp2[j].append(f((TimeNewObs[enum]-TimeNewSp3[i]+scale))) for enum in range(len(TimeNewObs)) if((obsSec[enum])>=sp3Sec[i] and (obsSec[enum])<sp3Sec[i+1])]
                    except:
                        # print('no')
                        timee=TimeNewSp3[-6:]
                        timee=[k-TimeNewSp3[-1] for k in timee]
                        i2=-1
                        f=np.poly1d(np.polyfit(timee,x[j,i-int((user_Sp3InterpolationDegree)/2):i+int(user_Sp3InterpolationDegree/2)+1],len(timee)-1))  
                        temp=[temp0[j].append(f((TimeNewObs[enum]-TimeNewSp3[i2]+scale))) for enum in range(len(TimeNewObs)) if((obsSec[enum])>=sp3Sec[i2])]
                        f=np.poly1d(np.polyfit(timee,y[j,i-int((user_Sp3InterpolationDegree)/2):i+int(user_Sp3InterpolationDegree/2)+1],len(timee)-1))  
                        temp=[temp1[j].append(f((TimeNewObs[enum]-TimeNewSp3[i2]+scale))) for enum in range(len(TimeNewObs)) if(obsSec[enum]>=sp3Sec[i2])]
                        f=np.poly1d(np.polyfit(timee,z[j,i-int((user_Sp3InterpolationDegree)/2):i+int(user_Sp3InterpolationDegree/2)+1],len(timee)-1))   
                        temp=[temp2[j].append(f((TimeNewObs[enum]-TimeNewSp3[i2]+scale))) for enum in range(len(TimeNewObs)) if((obsSec[enum])>=sp3Sec[i2])]
            

                    # f=np.poly1d(np.polyfit(timee,y[j,i-int((user_Sp3InterpolationDegree)/2):i+int(user_Sp3InterpolationDegree/2)+1],len(timee)-1))  
                    # temp=[temp1[j].append(f((TimeNewObs[enum]-TimeNewSp3[i]+scale))) for enum in range(len(TimeNewObs)) if(obsSec[enum]>=sp3Sec[i] and obsSec[enum]<sp3Sec[i+1])]

                    # f=np.poly1d(np.polyfit(timee,z[j,i-int((user_Sp3InterpolationDegree)/2):i+int(user_Sp3InterpolationDegree/2)+1],len(timee)-1))   
                    # temp=[temp2[j].append(f((TimeNewObs[enum]-TimeNewSp3[i]+scale))) for enum in range(len(TimeNewObs)) if((obsSec[enum])>=sp3Sec[i] and (obsSec[enum])<sp3Sec[i+1])]
            
            if(i<int(user_ClkInterpolationDegree/2)+1):
                for j in range(x.shape[0]):

                    f=np.poly1d(np.polyfit(TimeNewSp3[0:user_ClkInterpolationDegree+1],clk[j,0:user_ClkInterpolationDegree+1],user_ClkInterpolationDegree))
                    temp=[temp3[j].append(f(TimeNewObs[enum]+scale)) for enum in range(len(TimeNewObs)) if(obsSec[enum]>=sp3Sec[i] and obsSec[enum]<sp3Sec[i+1])]
            
            else:
                    for j in range(x.shape[0]):
                        try:
                            timee=TimeNewSp3[i-int((user_ClkInterpolationDegree)/2):i+int(user_ClkInterpolationDegree/2)+1]
                            timee=[k-TimeNewSp3[i] for k in timee]
                            f=np.poly1d(np.polyfit(timee,clk[j,i-int((user_ClkInterpolationDegree)/2):i+int(user_ClkInterpolationDegree/2)+1],len(timee)-1))   
                            temp=[temp3[j].append(f(TimeNewObs[enum]-TimeNewSp3[i]+scale)) for enum in range(len(TimeNewObs)) if(obsSec[enum]>=sp3Sec[i] and obsSec[enum]<sp3Sec[i+1])]
                        except:
                            timee=TimeNewSp3[-2:]
                            timee=[k-TimeNewSp3[i] for k in timee]
                            i2=-1
                            f=interpolate.interp1d(timee,clk[j,i-int((user_ClkInterpolationDegree)/2):i+int(user_ClkInterpolationDegree/2)+1],fill_value='extrapolate')
                            temp=[temp3[j].append(f(TimeNewObs[enum]-TimeNewSp3[i2]+scale)) for enum in range(len(TimeNewObs)) if(obsSec[enum]>=sp3Sec[i2])]
            

        X_PROC.append(np.array(temp0))
        Y_PROC.append(np.array(temp1))
        Z_PROC.append(np.array(temp2))
        CLK_PROC.append(np.array(temp3))

 
        Vx=(X_PROC[-1][:,2:]-X_PROC[-1][:,0:-2])/(2*obsRate);Vx=np.column_stack([Vx[:,0],Vx]);Vx=np.column_stack([Vx,Vx[:,-1]]);Vx_s.append(Vx)
        Vy=(Y_PROC[-1][:,2:]-Y_PROC[-1][:,0:-2])/(2*obsRate);Vy=np.column_stack([Vy[:,0],Vy]);Vy=np.column_stack([Vy,Vy[:,-1]]);Vy_s.append(Vy)
        Vz=(Z_PROC[-1][:,2:]-Z_PROC[-1][:,0:-2])/(2*obsRate);Vz=np.column_stack([Vz[:,0],Vz]);Vz=np.column_stack([Vz,Vz[:,-1]]);Vz_s.append(Vz)
