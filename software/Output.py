from Libs import *



def Output(path1,processed,P,epoch,sat_num,zhd,obs,obsRate):
    if(path1!=''):
        user_move='static'
        time_interv=int(obsRate);interv=0

        file=open(path1+'out.pos','w')
        indentation=' '*6
        file.write('observation file path: '+obsFile+'\n')
        file.write('observation file path: '+sp3File+'\n')
        file.write('ionospheric correction: '+'model\n')
        file.write('kalman filter parameters: \nQ'+indentation*3+str(kfParameters['Qx'])+'\nP0'+indentation*3+str(kfParameters['P0'])+'\nPHI'+indentation*3+str(kfParameters['PHI'])+'\n')
        file.write('reciever mode: '+user_move+'\n')
        file.write('Receiver a priori position: ['+str(processed[0][0][0])+' '+str(processed[0][1][0])+' '+str(processed[0][2][0])+' ]\n')
        file.write('Cut-Off angle(degree): '+str(userCutOff)+'\n')
        file.write('SNR mask: '+str(userSNRMask)+'\n')
        file.write('calculated ZHD value: '+str(zhd[0])+'\n')
        file.write('\n'*3)
        file.write('year'+indentation+'mon'+indentation+'day'+indentation+'hour'+indentation+'min'+indentation+'sec'+indentation+'       x       '+
                           indentation+'  y  '+indentation+'   z   '+indentation+'  ztd  '+
                           indentation+'sat_num'+indentation+'  var_x  '+indentation+'  var_y  '+indentation+' var_z  '+'\n\n')
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
                       indentation+str(round(processed[i+1][2][0],4))+indentation+str(round(processed[i+1][3][0]+zhd[0],4))+

                       indentation+str(sum(sat_num[i]))+indentation+str(round(p[0],6))+indentation+str(round(p[1],6))+indentation+str(round(p[2],6))+'\n')

        file.close()
