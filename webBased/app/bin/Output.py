from app.bin.Libs import *
from matplotlib.figure import Figure
import io
import base64


def Output(path1,processed,UNresidual,P,epoch,sat_num,q,phi,zhd,obs,OBSinterv,sat_name,epoch_item,pre_fit,post_fit,user_TimeInterval,user_CombType,file_o,file_precise,user_SNR_Mask,user_CutOff,user_FileIono,plots,SNR,gdop):
    user_move='static'
    if(path1!=''):
        
        
        time_interv=int(user_TimeInterval/OBSinterv);interv=0

        indentation=' '*6
        file=open('app/static/results/%s.pos'%(path1),'w')
        file_o=file_o.split('/')[-1]
        file_precise=file_precise.split('/')[-1]
        user_FileIono=user_FileIono.split('/')[-1]
        file.write('kntuPPP'+'\n\n\n')
        # file.write('observation file : '+file_o+'\n')
        # file.write('sp3 file : '+file_precise+'\n')
        if(user_CombType=='iono-free'):file.write('ionospheric correction : '+'first order iono-free combination\n')
        elif(user_CombType=='single-freq'):file.write('ionospheric correction : '+user_FileIono+'\n')
        file.write('kalman filter parameters: \nQ'+indentation*3+str(q)+'\nP0'+indentation*3+str(list(np.diag(P[0])[0:7]))+'\nPHI'+indentation*3+str(phi)+'\n')
        file.write('Receiver a priori position: ['+str(processed[0][0][0])+' '+str(processed[0][1][0])+' '+str(processed[0][2][0])+' ]\n')
        file.write('Cut-Off angle(degree): '+str(user_CutOff)+'\n')
        file.write('SNR mask: '+str(user_SNR_Mask)+'\n')
        file.write('calculated ZHD value: '+str(round(zhd[0],3))+'\n')
        file.write('\n'*3)
        file.write('year'+indentation+'mon'+indentation+'day'+indentation+'hour'+indentation+'min'+indentation+'sec'+indentation+'       x       '+
                           indentation+'  y  '+indentation+'   z   '+
                           indentation+'sat_num'+indentation+'  sigma_x  '+indentation+'  sigma_y  '+indentation+' sigma_z  '+'\n\n')
        for i in range(len(epoch)):
            if(str(obs[0].version)[0]=='2'):
                if(i!=0):interv+=time_interv
                time=obs[0].time.data[epoch[i]*time_interv]
            elif(str(obs[0].version)[0]=='3'):
                time=datetime.strptime(str(obs[0].time.data[epoch[i]*time_interv]).split('.')[0],"%Y-%m-%dT%H:%M:%S")
            p=np.diag(P[i+1])
            file.write(str(time.strftime('%Y'))+indentation+str(time.strftime('%m'))+indentation+                       
                       str(time.strftime('%d'))+indentation+str(time.strftime('%H'))+indentation+
                       str(time.strftime('%M'))+indentation+str(time.strftime('%S'))+indentation+
                       str(round(processed[i+1][0][0],4))+indentation+str(round(processed[i+1][1][0],4))+
                       indentation+str(round(processed[i+1][2][0],4))+
                       indentation+str(sum(sat_num[i]))+indentation+str(round(np.sqrt(p[0]),6))+indentation+str(round(np.sqrt(p[1]),6))+indentation+str(round(np.sqrt(p[2]),6))+'\n')

        file.close()
        geo_pos=[pymap3d.ecef2geodetic(enum[0],enum[1],enum[2]) for enum in processed]

        data=io.BytesIO();data2=io.BytesIO();data3=io.BytesIO();data4=io.BytesIO();data5=io.BytesIO()
        dataEncoded=[]
        #

        if('plt1' in plots):
            # fig1 = Figure()
            fig1,ax1 = plt.subplots(1,1)
            [ax1.plot(SNR[sat][0,:,:]) for sat in range(len(SNR))]
            sats=[];[[sats.append(enum2) for enum2 in enum.sv.data] for enum in obs]
            if('Geo' in sats):sats.remove('Geo')
            plt.xlabel('Time(epochs)');plt.ylabel('value(meters)');plt.title('SNR(dbHz)');plt.legend(sats,fontsize=8,loc="right");fig1.savefig(data, format="png")
            dataEncoded1=base64.b64encode(data.getbuffer())
            dataEncoded.append(dataEncoded1)
            plt.close()

        if('plt2' in plots):
            # fig2 = Figure()
            fig2,ax2 = plt.subplots()
            ax2.plot([enum[0] for enum in UNresidual[20:]]);ax2.plot([enum[1] for enum in UNresidual[20:]]);ax2.plot([enum[2] for enum in UNresidual[20:]]);plt.legend(['X','Y','Z'])
            plt.xlabel('Time(epochs)');plt.ylabel('value(meters)');plt.title('residual');fig2.savefig(data2, format="png")
            dataEncoded2=base64.b64encode(data2.getbuffer())
            dataEncoded.append(dataEncoded2)
            plt.close()


        if('plt3' in plots):
            # fig3 = Figure()
            fig3,ax3 = plt.subplots()
            ax3.plot([i[0][0] for i in geo_pos[20:]],[i[1][0] for i in geo_pos[20:]],'.')
            plt.xlabel('East(ECEF coordinate)');plt.ylabel('North(ECEF coordinate)');plt.title('Horizontal coordiante');fig3.savefig(data3, format="png")
            dataEncoded3=base64.b64encode(data3.getbuffer())
            dataEncoded.append(dataEncoded3)
            plt.close()


        if('plt4' in plots):
            # fig4 = Figure()
            fig4,ax4 = plt.subplots()
            ax4.plot([(enum[4][0]/C)*np.power(10,6) for enum in processed[10:]],'.')
            plt.xlabel('Time(epochs)');plt.ylabel('value(micro-seconds)');plt.title('reciever-clock');fig4.savefig(data4, format="png")
            dataEncoded4=base64.b64encode(data4.getbuffer())
            dataEncoded.append(dataEncoded4)
            plt.close()

        if('plt5' in plots):
            # fig4 = Figure()
            fig5,ax5 = plt.subplots()
            ax5.plot(gdop[20:])
            plt.xlabel('Time(epochs)');plt.ylabel('value(gdop)');plt.title('GDOP');fig5.savefig(data5, format="png")
            dataEncoded5=base64.b64encode(data5.getbuffer())
            dataEncoded.append(dataEncoded5)
            plt.close()

        return dataEncoded
    
        
