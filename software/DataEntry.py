from Libs import *



def dataEntry(OBS,SP3):
    for sat in userSats:
        OBS.append(gr.load(obsFile,use=[sat],meas=userFreq[sat]+userFreqF2[sat],tlim=userTime))
        SP3.append(gr.load(sp3File,use=[sat]))