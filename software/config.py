'''
     user configuration file

'''


#INPUTS

obsFile='./testFile/GEOP140I_A.19o'

sp3File="./testFile/COD0MGXFIN_20191400000_01D_05M_ORB.SP3"

outFile='./output/'

ionoFile='./testFile/data140'

userTime=['2019-05-20T08:11','2019-05-20T08:20']

userSats=['G','E']
userCutOff=15
userSNRMask=30

#if initPos==[] => the initialization position for receiver will be chosen based on "".obs file"
#otherwise please select this parameters: initPos=[x(ECEF),y(ECEF),z(ECEF)]
initPos=[]

#ionospheric model (items:[single-freq,iono-free])
ionoMode='single-freq'

#process mode (items:[code-only,combined])
userProcMode='combined'


userFreq={'G':['C1C','L1C','D1C','S1C'],'E':['C1C','L1C','D1C','S1C']}
#if you're using single-freq you can set userFreqF2={}
userFreqF2={'G':['C5Q','L5Q','D5Q','S5Q'],'E':['C5Q','L5Q','D5Q','S5Q']}



#Measurement Noise [ CN0 (if so; customNoise=[[STD of codeNoise],[STD of phaseNoise]]), 
#                    elev (if so; customNoise=[[STD of codeNoise],[STD of phaseNoise]]),
#                    constant (if so; customNoise=[[STD of codeNoise],[STD of phaseNoise]]),
#                    custom (if so; customNoise=[[STD of codeNoise for each sat],[STD of phaseNoise for each sat]])
#                   ]
userMeasNoise='constant'

customNoise=[3,0.003]
# customNoise=[[[3.01, 2.49, 3, 2.08, 1.76, 2.16, 1.59, 1.94, 1.76, 2.04, 3, 2.04], [3,1.72, 2.19, 3, 1.86, 2.23, 1.15]],
#              [[0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.004, 0.004, 0.007, 0.005],[0.007, 0.005, 0.003, 0.007, 0.005, 0.002, 0.004]]] 


#Threshold to remove outliers from code Measurement residual ( codeThreshold= n*sigmaNoise (meters))
codeThreshold=4
#Threshold to remove outliers from phase Measurement residual( phaseThreshold= n*sigmaNoise (meters))
phaseThreshold=5

# kalman filter common parameters [xReceiver,yReceiver,zReceiver,zwd,clk,isb,N]
kfParameters={'Qx':[0,0,0,0.25,9e5,9e3,100],
              'PHI':[1,1,1,0,0,0,1],
              'P0':[1.00e+2,1.00e+2,1.00e+2,0.25,9e5,9e3,100]
              }


Files={
    'EGM':"./geoids/egm2008-5.pgm",
    'ephemeris':'./eph/de405.bsp'
    }