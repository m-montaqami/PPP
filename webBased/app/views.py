from app import app 
from flask import render_template, request, redirect
from app.bin.main import exe
from app.bin.AV import avAnalysis
import os

def proc(form):
    #user parameters
    plots=[]
    freq={'1':'C1','2':'P1','3':'C1C','4':'C1-C2','5':'P1-P2','6':'C1-C5','7':'C1C-C5X','8':'C1C-C5Q'}
    freq=freq[form['freq']]
    if('-' in freq):
        iono='iono-free'
        freq=freq.split('-')
        L_freq=['L'+enum[1:] for enum in freq]
        S_freq=['S'+enum[1:] for enum in freq]
        Freq={'G':[freq[0],L_freq[0],freq[1],L_freq[1]]}
        SNR_Freq={'G':[S_freq[0],S_freq[1]]}
    else:
        iono='single-freq'
        L_freq='L'+freq[1:]
        S_freq='S'+freq[1:]
        Freq={'G':[freq,L_freq,freq,L_freq]}
        SNR_Freq={'G':[S_freq,S_freq]}
    gdop=int(form['gdop'])
    snr_mask=int(form['snr'])
    elv_mask=int(form['elev'])
    if('use_phase' in list(form.keys())):
        comb_type='combined'
    else:
        comb_type='code-only'
    for enum in ['plt1','plt2','plt3','plt4','plt5']:
        if(enum in list(form.keys())):plots.append(enum)
    # print('plots: '+str(plots))
    #save uploaded file by user
    name=form['email'].split('@')[0].split('.')[0]
    obs_file = request.files['obs_file']
    obs_file.save('app/bin/downloads/downs_%s.o'%(name))
    obs_file.close()
    sp3_file = request.files['sp3_file']
    sp3_file.save('app/bin/downloads/downs_%s.sp3'%(name))
    sp3_file.close()
    ionex_file,dcb_file='',''
    if(int(form['freq'])<=3):
        ionex_file=request.files['ionex_file']
        ionex_file.save('app/bin/downloads/downs_%s.ionex'%(name))
        ionex_file="app/bin/downloads/downs_%s.ionex"%(name)
    
    if(request.files['dcb_file'].filename!=''):
        dcb_file=request.files['dcb_file']
        dcb_file.save('app/bin/downloads/downs_%s.dcb'%(name))
    dataEncoded=exe("app/bin/downloads/downs_%s.o"%(name),"app/bin/downloads/downs_%s.sp3"%(name),name+'_out',['G'],elv_mask,snr_mask,comb_type,iono,Freq,SNR_Freq,gdop,plots,ionex_file)

    os.system('rm %s'%('app/bin/downloads/downs_%s.o'%(name)))
    os.system('rm %s'%('app/bin/downloads/downs_%s.sp3'%(name)))
    if not('-' in freq):os.system('rm %s'%('app/bin/downloads/downs_%s.ionex'%(name)))
    return name,dataEncoded
    
@app.route('/',methods=['GET','POST'])
def index():
    return render_template('public/index.html')

@app.route('/error')
def error():

    return render_template('public/error.html')


@app.route('/contact')
def contact():
    return render_template('public/contact.html')


@app.route('/kntuapp',methods=['GET','POST'])
def kntuapp():
    if(request.method=='POST'):
        try:
            name,dataEncoded=proc(request.form)
            if(len(dataEncoded)==0):
                return render_template('public/result.html',
                            res1="../../static/results/%s_out.pos "%(name))
            elif(len(dataEncoded)==1):
                return render_template('public/result.html',
                            res1="../../static/results/%s_out.pos "%(name),
                            fig1=dataEncoded[0].decode("ascii"))
            elif(len(dataEncoded)==2):
                return render_template('public/result.html',
                            res1="../../static/results/%s_out.pos "%(name),
                            fig1=dataEncoded[0].decode("ascii"),
                            fig2=dataEncoded[1].decode("ascii"))
            elif(len(dataEncoded)==3):
                return render_template('public/result.html',
                            res1="../../static/results/%s_out.pos "%(name),
                            fig1=dataEncoded[0].decode("ascii"),
                            fig2=dataEncoded[1].decode("ascii"),
                            fig3=dataEncoded[2].decode("ascii"))
            elif(len(dataEncoded)==4):
                return render_template('public/result.html',
                            res1="../../static/results/%s_out.pos "%(name),
                            fig1=dataEncoded[0].decode("ascii"),
                            fig2=dataEncoded[1].decode("ascii"),
                            fig3=dataEncoded[2].decode("ascii"),
                            fig4=dataEncoded[3].decode("ascii"))
            elif(len(dataEncoded)==5):
                return render_template('public/result.html',
                            res1="../../static/results/%s_out.pos "%(name),
                            fig1=dataEncoded[0].decode("ascii"),
                            fig2=dataEncoded[1].decode("ascii"),
                            fig3=dataEncoded[2].decode("ascii"),
                            fig4=dataEncoded[3].decode("ascii"),
                            fig5=dataEncoded[4].decode("ascii"))

        except:
            return redirect('/error')
    return render_template('public/app.html')


@app.route('/av',methods=['GET','POST'])
def av():
    if(request.method=='POST'):
        try:
            form_val=list(request.form.values())
            name=form_val[0].split('@')[0].split('.')[0]
            file=request.files['obs_file']
            file.save('app/bin/downloads/downs_%s.o'%(name))
            file.close()
            freq={'4':'C1-C2','5':'P1-P2','6':'C1-C5','7':'C1C-C5X','8':'C1C-C5Q'}
            freq=freq[form_val[1]]
            freq=list(freq.split('-'))
            freq=[freq[0],'L'+freq[0][1:],freq[1],'L'+freq[1][1:]]
            res=avAnalysis('app/bin/downloads/downs_%s.o'%(name),['G'],freq)
            os.system('rm %s'%('app/bin/downloads/downs_%s.o'%(name)))
            return render_template('public/result2.html',fig1=res[0].decode("ascii"),
                                                        fig2=res[1].decode("ascii"),
                                                        fig3=res[2].decode("ascii"),
                                                        fig4=res[3].decode("ascii"))
            
        except:
            return redirect('/error')
    if(request.method=='GET'):
        return render_template('public/av.html')

