from Libs import *
from IonRead import *

class atmosphere():
    def __init__(self,obs,obsSec,elevation,azimuth,freq,preproc,ZHD,MF_Trop,obs_rate):
        self.obs=obs
        self.obs_rate=obs_rate
        self.TimeSec=[obsSec[0],obsSec[-1]]
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
        ginterpolator = pygeodesy.GeoidKarney(Files['EGM'])
        geoidHeight=ginterpolator(pygeodesy.ellipsoidalKarney.LatLon(self.lat,self.lon))
        self.h_ortho=self.h-geoidHeight
        self.freq=freq
        self.freq_val=userFreq[self.freq]
        self.freq_val=[enum[0:2] for enum in self.freq_val]
        self.L=[userFreq[self.freq].index(i) for i in userFreq[self.freq] if 'L' in i]   #code index
        self.P=[userFreq[self.freq].index(i) for i in userFreq[self.freq] if not 'L' in i]    #carrier-phase index
        
        self.freqComb=Freq[self.freq][self.freq_val[self.L[0]]]


    def iono(self,GIM,VTEC,MF_Iono):
        def ReadFile():

            if ionoMode=='iono-free':
                VTEC.append([])
            else:
                RawData=list(open(ionoFile))
                inds=[]
                for enum in range(len(RawData)):
                    if("EPOCH OF CURRENT MAP" in RawData[enum]):
                        inds.append(enum+RawData[enum:].index(RawData[enum]))
                TECDate=[RawData[enum0] for enum0 in inds]
                TECDate=[enum.split(' ') for enum in TECDate]
                TECDate=[[enum0 for enum0 in enum if enum0!=''][0:6] for enum in TECDate]  #format: yr,month,day,hour,min,sec

                TECDate=[int(np.round(datetime.timestamp(datetime.strptime('-'.join([enum for enum in enum2]),"%Y-%m-%d-%H-%M-%S")))) for enum2 in TECDate]
                TECDate_rate=int(TECDate[1]-TECDate[0])

                starter=[np.abs(enum-self.TimeSec[0]) for enum in TECDate];starter=starter.index(min(starter))
                end=[np.abs(enum-self.TimeSec[1]) for enum in TECDate];end=end.index(min(end))
                if(starter!=0):
                    starter-=1
                if(end<len(TECDate)-1):
                    end+=2

                interpol_time=list(range(0,(int(self.TimeSec[-1]-self.TimeSec[0])*self.obs_rate)+1,self.obs_rate))
                TECDate=TECDate[starter:end]
                Tscale=TECDate[0]
                TECDate=[enum2-Tscale for enum2 in TECDate]

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
                    VTEC_buff.append(DesiredTEC[0])
                VTEC_buff=VTEC_buff[starter:end]
                f=interpolate.interp1d(TECDate,VTEC_buff)
                VTEC_buff=f(interpol_time)
                VTEC.append(np.array([[((40.3*1e15)/math.pow(self.freqComb,2))*enum2,((-40.3*1e15)/math.pow(self.freqComb,2))*enum2] for enum2 in VTEC_buff]))
                
            
        def MF_I():
            z_prim=np.arcsin((R_Earth*np.sin(np.radians(self.azimuth)))/(R_Earth+450))
            MF_Iono.append(1/np.cos(z_prim))

        ReadFile()
        MF_I()

    def trop(self,ZHD,MF_Trop):

        def ZHD1():            
            
            #Sastamainen
            #h in meters(h from ellipsoid(wgs84))
            #P0 in hPa
            #lat in degree
            if(ZHD==[]):
                P0=1013.25*np.power((1-(2.2557*math.pow(10,-5)*self.h_ortho)),5.2568)
                ZHD.append((0.0022768*P0)/((1-((0.00266)*math.cos(2*math.radians(self.lat))))-(0.00028*(self.h*math.pow(10,-3)))))  #  h should be geodetic

        def MF_T():
            elev=np.radians(self.elev)
            az=np.radians(self.azimuth)
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
            #if(len(params_dry)==1):params_dry=params_dry[0];params_wet=params_wet[0]
            #else:pass   #!   interpolate
            
            delta_m_d=((1/np.sin(elev))-((1+(a_h/(1+(b_h/(1+c_h)))))/(np.sin(elev)+(a_h/(np.sin(elev)+(b_h/(np.sin(elev)+c_h)))))))*self.h_ortho*math.pow(10,-3)
            mf_d=((1+(a_dry/(1+(b_dry/(1+c_dry)))))/(np.sin(elev)+(a_dry/(np.sin(elev)+(b_dry/(np.sin(elev)+c_dry))))))+delta_m_d
            self.n=mf_d
            mf_w=((1+(a_wet/(1+(b_wet/(1+c_wet)))))/(np.sin(elev)+(a_wet/(np.sin(elev)+(b_wet/(np.sin(elev)+c_wet))))))
            MF_Trop.append([mf_d,mf_w,mf_w*(1/(np.tan(elev)))*np.cos(az),mf_w*(1/(np.tan(elev)))*np.sin(az)])
                    

        ZHD1();MF_T()


