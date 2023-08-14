from app.bin.Libs import *



def dataEntry(OBS,SP3,user_Sats,file_o,file_precise,user_Freq,user_SNR_Freq,user_Time):
    
    if(software['remote-mode']=='online'):
        pass
    elif(software['remote-mode']=='offline'):
        for sat in user_Sats:
            OBS.append(gr.load(file_o,use=[sat],meas=user_Freq[sat]+user_SNR_Freq[sat]))
            SP3.append(gr.load(file_precise,use=[sat]))
        return 