from Libs import *
import requests
import os

def SendRequest(url,timeout):
    res=requests.get(url,allow_redirects=True,timeout=timeout)
    if(res.status_code==200):
        print("url: %s fetched"%(url))
        return res
    else:
        print("can not fetch the URL(please Check internet Connection or Input URL)")
        return -1


def Write2File(data):
    if(data!=-1):
        FileName='latest-%s'%(str(int(time.time())))+'.Z'
        file=open(FileName,'wb')
        file.write(data.content)
        file.close()
        print("file with name:%s saved to local"%(FileName))
        return FileName
    elif(data==-1):
        return data


def UnZip(FileName):
    try:
        os.system("uncompress %s"%(FileName))
        print("Uncompressed the file")
        return 0
    except:
        print("can not uncompress file")
        -2  #can not uncompress file

def TextProcess(RawData,lat,lon):
    #read IONEX 
    date=RawData[[RawData.index(enum) for enum in RawData if "EPOCH OF CURRENT MAP" in enum][0]]
    RawData=RawData[[RawData.index(enum) for enum in RawData if "EPOCH OF CURRENT MAP" in enum][0]+1:]
    geomap={str(float(RawData[enum][3:8])):''.join(enum2[:-1] for enum2 in RawData[enum+1:enum+6]) for enum in range(0,len(RawData)-6,6)}
    Lons=[str(float(enum)) for enum in range(-180,181,5)]
    buffer=geomap[lat].split(' ')
    buffer=[enum for enum in buffer if enum!='']
    TECU=buffer[Lons.index(lon)]
    return TECU,date

def exe(url="",lat="",lon="",timeout=5):
    if(url==""):
        print("enter RawData url First")
        return -1
    data=SendRequest(url,timeout)
    FileName=Write2File(data)
    res=UnZip(FileName)
    if(res==-2):return res
    try:RawData=list(open(FileName[:-2]))
    except:
        return -2
    TECU,date=TextProcess(RawData,lat,lon)
    #os.system('rm %s'%(FileName[:-2]))
    return TECU,date
    
    