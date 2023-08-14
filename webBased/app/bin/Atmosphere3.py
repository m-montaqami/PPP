from app.bin.Libs import *
from app.bin.IonRead import *
from app.bin.config_admin import Files

class atmosphere():
    def __init__(self,obs,obsSec,elevation,azimuth,freq,preproc,ZHD,MF_Trop,DATATIMEGFS,DATATIMEECMWF,obs_rate,user_Freq,mode):
        self.obs=obs
        self.obs_rate=obs_rate
        self.obsSec=obsSec
        # print('atmosphere3 obsSec: '+str(len(obsSec)))
        self.TimeSec=[self.obsSec[0],self.obsSec[-1]]
        self.preproc=preproc
        self.Geodetic=obs.position_geodetic
        if(str(self.obs.version)[0]=='2'):
            self.DOY=self.obs.time.data[0].strftime('%j')
        elif(str(self.obs.version)[0]=='3'):
            buff=datetime.strptime(str(self.obs.time.data[0]).split('.')[0],"%Y-%m-%dT%H:%M:%S")
            self.DOY=buff.strftime('%j')
        self.lat=self.Geodetic[0]
        self.lon=self.Geodetic[1]
        self.h=self.Geodetic[2]  #geodetic(meters)
        self.elev=elevation
        self.azimuth=azimuth
        self.DATATIMEGFS=DATATIMEGFS
        self.DATATIMEECMWF=DATATIMEECMWF
        ginterpolator = pygeodesy.GeoidKarney(Files['EGM'])
        geoidHeight=ginterpolator(pygeodesy.ellipsoidalKarney.LatLon(self.lat,self.lon))
        self.h_ortho=self.h-geoidHeight
        self.freq=freq
        self.freq_val=user_Freq[self.freq]
        self.freq_val=[enum[0:2] for enum in self.freq_val]
        self.L=[user_Freq[self.freq].index(i) for i in user_Freq[self.freq] if 'L' in i]   #code index
        self.P=[user_Freq[self.freq].index(i) for i in user_Freq[self.freq] if not 'L' in i]    #carrier-phase index
        
        if(mode=='iono-free'):
            self.alpha=(math.pow(Freq[self.freq][self.freq_val[self.L[0]]],2))/(math.pow(Freq[self.freq][self.freq_val[self.L[0]]],2)-math.pow(Freq[self.freq][self.freq_val[self.L[1]]],2))
            self.beta=(math.pow(Freq[self.freq][self.freq_val[self.L[1]]],2))/(math.pow(Freq[self.freq][self.freq_val[self.L[0]]],2)-math.pow(Freq[self.freq][self.freq_val[self.L[1]]],2))
            #freq Ionospheric-free
            self.freqComb=self.alpha*Freq[self.freq][self.freq_val[self.L[0]]]-self.beta*Freq[self.freq][self.freq_val[self.L[1]]]
        else:
            #freq Single-frequency
            self.freqComb=Freq[self.freq][self.freq_val[self.L[0]]]


    def iono(self,GIM,VTEC,MF_Iono,user_CombType,user_FileIono,user_RT):
        def ReadFile():
            if(user_CombType=='iono-free'):
                VTEC.append([])

            elif(user_RT==True and user_CombType=='single-freq'):
                res=requests.get(Ex_Links['itrg'],allow_redirects=True,timeout=Ex_Links['RequestTimeOut'])
                FileName='latest-%s'%(str(int(time.time())))+'.Z'
                file=open(FileName,'wb');file.write(res.content);file.close()
                if(software['OS']=='linux'):os.system("uncompress %s"%(FileName));os.system("mv %s archive/iono/"%(FileName[:-2]))
                elif(software['OS']=='windows'):pass
                
            elif(user_RT==False and user_CombType=='single-freq'):
                FileName=user_FileIono
            if(user_CombType=='single-freq'):
                if(user_RT==True):
                    RawData=list(open('archive/iono/'+FileName[:-2]))
                elif(user_RT==False):
                    RawData=list(open(FileName))

                inds=[]
                for enum in range(len(RawData)):
                    if("EPOCH OF CURRENT MAP" in RawData[enum]):
                        inds.append(enum+RawData[enum:].index(RawData[enum]))
                if(user_RT==False):
                    TECDate=[RawData[enum0] for enum0 in inds]
                    TECDate=[enum.split(' ') for enum in TECDate]
                    TECDate=[[enum0 for enum0 in enum if enum0!=''][0:6] for enum in TECDate]  #format: yr,month,day,hour,min,sec

                    # print(TECDate[0:3])
                    TECDate=[int(np.round(datetime.timestamp(datetime.strptime('-'.join([enum for enum in enum2]),"%Y-%m-%d-%H-%M-%S")))) for enum2 in TECDate]
                    # print(TECDate)
                    TECDate_rate=int(TECDate[1]-TECDate[0])

                    starter=[np.abs(enum-self.TimeSec[0]) for enum in TECDate];starter=starter.index(min(starter))
                    end=[np.abs(enum-self.TimeSec[-1]) for enum in TECDate];end=end.index(min(end))
                    if(starter!=0):
                        starter-=1
                    if((end-starter)<2):
                        end+=2
                    # print('starter'+str(starter))
                    # print('end: '+str(end))
                    # interpol_time=list(range(0,(int(self.TimeSec[-1]-self.TimeSec[0])*self.obs_rate)+1,self.obs_rate))
                    # print('obs_rate: '+str(self.obs_rate))
                    # print('TimeSec: '+str(self.TimeSec))
                    # interpol_time=list(range(0,(int(self.TimeSec[-1]-self.TimeSec[0]))*self.obs_rate,self.obs_rate))
                    interpol_time=list(range(0,(int(len(self.obsSec)))*self.obs_rate,self.obs_rate))
                    TECDate=TECDate[starter:end]
                    Tscale=TECDate[0]
                    # print(interpol_time[-1])
                    if(interpol_time[-1]>TECDate[-1]):
                        interpol_time=[enum2-Tscale for enum2 in interpol_time]
                    # print(interpol_time[-1])
                    TECDate=[enum2-Tscale for enum2 in TECDate]
                    # print(TECDate)



                inds.append(len(RawData)-1)
                VTEC_buff=[]
                for indx in range(len(inds)-1):
                    RawData2=RawData[inds[indx]+1:inds[indx+1]].copy()
                    geomap={str(float(RawData2[enum][3:8])):''.join(enum2[:-1] for enum2 in RawData2[enum+1:enum+6]) for enum in range(0,len(RawData2)-6,6)}
                    Lons=[str(float(enum)) for enum in range(-180,181,5)]
                    Lats=list(geomap.keys())
                    GIM.append([np.array([[float(j) for j in i.split(' ') if j!=''] for i in list(geomap.values())])])
                    buff0=[abs(self.lat-float(enum0)) for enum0 in Lats];buff0=buff0.index(min(buff0))
                    buff02=[abs(self.lon-float(enum0)) for enum0 in Lons];buff02=buff02.index(min(buff02))
                    if(buff0<5):buff0=5
                    elif(buff0>65):buff0=65
                    if(buff02<5):buff02=5
                    elif(buff02>67):buff02=67            
                    Grid=[enum0.split(' ') for enum0 in list(geomap.values())];GIM.append([[int(enum02) for enum02 in enum0 if enum02!=''] for enum0 in Grid])
                    DesiredGrid=[geomap[Lats[enum0]] for enum0 in range(buff0-5,buff0+5)]
                    DesiredGrid=[enum0.split(' ') for enum0 in DesiredGrid]
                    DesiredGrid=np.array([[int(enum02) for enum02 in enum0 if enum02!=''] for enum0 in DesiredGrid]) #TECU
                    DesiredGrid=DesiredGrid[:,buff02-5:buff02+5]
                    GridX=[float(Lats[enum0]) for enum0 in range(buff0-5,buff0+5)]
                    GridY=[float(Lons[enum0]) for enum0 in range(buff02-5,buff02+5)]
                    DesiredTEC=self.preproc.Interpolation([GridX,GridY],DesiredGrid,[self.lat,self.lon],'2')


                    # VTEC_buff.append([((-40.3)/math.pow(self.freqComb,2))*DesiredTEC[0],((40.3)/math.pow(self.freqComb,2))*DesiredTEC[0]])  #[C,P])
                    # VTEC_buff.append(DesiredTEC[0]*1e16)
                    VTEC_buff.append(DesiredTEC[0])
                if(user_RT==False):
                    VTEC_buff=VTEC_buff[starter:end]
                    f=interpolate.interp1d(TECDate,VTEC_buff,fill_value='extrapolate')
                    VTEC_buff=f(interpol_time)
                    VTEC.append(np.array([[((-40.3*1e15)/math.pow(self.freqComb,2))*enum2,((40.3*1e15)/math.pow(self.freqComb,2))*enum2] for enum2 in VTEC_buff]))
                else:
                    VTEC.append(np.array([[((-40.3*1e15)/math.pow(self.freqComb,2))*VTEC_buff[0],((40.3*1e15)/math.pow(self.freqComb,2))*VTEC_buff[0]]]))
                    
                # print(VTEC_buff)
                
            
        def MF_I():
            user_MFIono='SLM'
            if(user_MFIono=='SLM'):
                # print(self.azimuth)
                z_prim=np.arcsin((R_Earth*np.sin(np.radians(self.azimuth)))/(R_Earth+450))
                MF_Iono.append(1/np.cos(z_prim))

        ReadFile()
        MF_I()

    def trop(self,ZHD,MF_Trop,user_RT):
        user_TropMF='neill'
        user_SurfacePressure='saastamoinen'

        def ReadFile():
            pass

        def ZHD1():            
            
            #Sastamainen
            #h in meters(h from ellipsoid(wgs84))
            #P0 in hPa
            #lat in degree
            if(ZHD==[]):
                if(user_RT==True):
                    res=requests.get(Ex_Links['GFS']['url'][0]+self.DATATIMEGFS[-1][0]+'/'+self.DATATIMEGFS[-1][1]+Ex_Links['GFS']['url'][1]+Ex_Links['GFS']['item'][0]+self.DATATIMEGFS[-1][1]+Ex_Links['GFS']['item'][1])
                    file=open(Ex_Links['GFS']['downloadPath'],'wb');file.write(res.content);file.close()


                # elif(user_SurfacePressure=='ECMWF' and user_RT==False):
                #     c=cdsapi.Client()
                #     req=Ex_Links['ECMWF']['request']
                #     req['year']=self.DATATIMEECMWF[0]
                #     req['month']=self.DATATIMEECMWF[1]
                #     req['day']=self.DATATIMEECMWF[2]

                #     req['area']=[self.lon+5,self.lat-5,self.lon-5,self.lat+5]
                #     c.retrieve(Ex_Links['ECMWF']['product'],req,Ex_Links['ECMWF']['downloadPath'])
                #     b=cfgrib.open_datasets(Ex_Links['ECMWF']['downloadPath'])
                #     lat=b[0].latitude;lon=b[0].longitude;z=b[0].z
                #     p=[float(enum0) for enum0 in Ex_Links['ECMWF']['request']['pressure_level']]
                #     orth=[self.preproc.Interpolation([lat,lon],z[enum0,:,:],[self.lat,self.lon],'2')[0]/9.8 for enum0 in range(z.shape[0])]
                #     orth.reverse()
                #     P0=self.preproc.Interpolation(orth,p,self.h_ortho,'1')


                elif(user_SurfacePressure=='MOP'):
                    if(abs(self.lat)>75):
                        P0=[CONST['Trop']['ZHD']['P0']['MOP'][75]]
                    elif(abs(self.lat)<15):
                        P0=[CONST['Trop']['ZHD']['P0']['MOP'][15]]
                    else:
                        ids=[abs(float(enum)-float(self.lat)) for enum in list(CONST['Trop']['ZHD']['P0']['MOP'].keys())]
                        ids=ids.index(min(ids))
                        vals=list(CONST['Trop']['ZHD']['P0']['MOP'].keys())
                        P0=[CONST['Trop']['ZHD']['P0']['MOP'][vals[enum]] for enum in range(ids-1,ids+2)]
                        #interpolate
                        P0=float(self.preproc.Interpolation([float(enum) for enum in vals[ids-1:ids+2]],P0,float(self.lat),'1'))

                elif(user_SurfacePressure=='saastamoinen'):
                    P0=1013.25*np.power((1-(2.2557*math.pow(10,-5)*self.h_ortho)),5.2568)

                # ZHD.append((0.0022768*P0)/((1-((0.00266)*math.cos(2*math.radians(self.lat))))-(0.00028*(self.h*math.pow(10,-3)))))  #  h should be geodetic

                ZHD.append((2.3)*np.exp(-0.116*math.pow(10,-3)*(self.h)))

        def MF_T():
            elev=np.radians(self.elev)
            az=np.radians(self.azimuth)
            if(user_TropMF=='neill'):
                if(self.lat<0):DOY0=211
                else:DOY0=28
                if(abs(self.lat)>75):
                    params_dry=CONST['Trop']['MF']['neill']['dry'][75]
                    params_wet=CONST['Trop']['MF']['neill']['wet'][75]
                elif(abs(self.lat)<15):
                    params_dry=CONST['Trop']['MF']['neill']['dry'][15]
                    params_wet=CONST['Trop']['MF']['neill']['wet'][15]
                else:
                    ids=[abs(float(enum)-float(np.abs(self.lat))) for enum in list(CONST['Trop']['MF']['neill']['dry'].keys())[:-3]]
                    ids=ids.index(min(ids))
                    if(ids==0):
                        ids=1
                    elif(ids==4):
                        ids=3
                    vals=list(CONST['Trop']['MF']['neill']['dry'].keys())[:-3]
                    params_dry=[np.array(CONST['Trop']['MF']['neill']['dry'][vals[enum]]) for enum in range(len(vals))]
                    params_wet=[np.array(CONST['Trop']['MF']['neill']['wet'][vals[enum]]) for enum in range(len(vals))]

                    params_dry=[float(self.preproc.Interpolation(vals,[enum2[enum] for enum2 in params_dry],float(self.lat),'1')) for enum in range(len(params_dry[0]))]
                    params_wet=[float(self.preproc.Interpolation(vals,[enum2[enum] for enum2 in params_wet],float(self.lat),'1')) for enum in range(len(params_wet[0]))]
                  
                a_h=np.array(CONST['Trop']['MF']['neill']['dry']['a_h'])*math.pow(10,-5)
                b_h=np.array(CONST['Trop']['MF']['neill']['dry']['b_h'])*math.pow(10,-3)
                c_h=np.array(CONST['Trop']['MF']['neill']['dry']['c_h'])*math.pow(10,-3)
                #interpolation

                a_dry=(params_dry[0]*math.pow(10,-3))-(params_dry[3]*math.pow(10,-5)*math.cos(2*math.pi*((float(self.DOY)-DOY0)/365.25)))
                b_dry=(params_dry[1]*math.pow(10,-3))-(params_dry[4]*math.pow(10,-5)*math.cos(2*math.pi*((float(self.DOY)-DOY0)/365.25)))
                c_dry=(params_dry[2]*math.pow(10,-3))-(params_dry[5]*math.pow(10,-5)*math.cos(2*math.pi*((float(self.DOY)-DOY0)/365.25)))
                a_wet=params_wet[0]*math.pow(10,-4);b_wet=params_wet[1]*math.pow(10,-3);c_wet=params_wet[2]*math.pow(10,-2)
                
                delta_m_d=((1/np.sin(elev))-((1+(a_h/(1+(b_h/(1+c_h)))))/(np.sin(elev)+(a_h/(np.sin(elev)+(b_h/(np.sin(elev)+c_h)))))))*self.h_ortho*math.pow(10,-3)
                mf_d=((1+(a_dry/(1+(b_dry/(1+c_dry)))))/(np.sin(elev)+(a_dry/(np.sin(elev)+(b_dry/(np.sin(elev)+c_dry))))))+delta_m_d
                self.n=mf_d
                mf_w=((1+(a_wet/(1+(b_wet/(1+c_wet)))))/(np.sin(elev)+(a_wet/(np.sin(elev)+(b_wet/(np.sin(elev)+c_wet))))))
                MF_Trop.append([mf_d,mf_w,mf_w*(1/(np.tan(elev)))*np.cos(az),mf_w*(1/(np.tan(elev)))*np.sin(az)])

        ZHD1();MF_T()


