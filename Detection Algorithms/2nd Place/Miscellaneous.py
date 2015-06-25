
# coding: utf-8

# In[ ]:

###################################################################
## Generate spectrogram features using short-time Fourier transform
###################################################################
import os
from scipy.io import loadmat
from scipy.signal import resample
import pickle

import scipy, numpy as np

def stft(x, fftsize=1024, overlap=4):   
    hop = fftsize / overlap
    w = scipy.hanning(fftsize+1)[:-1]      
    return np.array([np.fft.rfft(w*x[i:i+fftsize]) for i in range(0, len(x)-fftsize, hop)])

for i in range(4):
    print i
    sys.stdout.flush()
    dn = './data/clips/Dog_%d/'%(i+1)
    fns = [fn for fn in os.listdir(dn) if '.mat' in fn]
    
    feat = []
    labels = []
    
    for fn in fns:
        m = loadmat(dn+fn)
        d = m['data']
        d = resample(d, int(d.shape[1]*500/m['freq']), axis=1)
        s = array([stft(x,64,8) for x in d])
        s = sqrt(real(s)**2+imag(s)**2)    
        if 'inter' in fn:
            l = 0
        elif '_ictal' in fn:
            l = 1
        else:
            l = -1
        labels.append(l)
        feat.append(s.flatten())    
        
    feat = array(feat)
    labels = array(labels)
    fns = array(fns)
    
    d = {}
    d['feat'] = feat
    d['labels'] = labels
    d['feat_norm'] = (feat-feat.mean(0))/feat.std(0)
    d['fns'] = fns
    
    savez('./data/data_Dog_%d'%(i+1), **d)

for i in range(8):
    print i
    sys.stdout.flush()
    dn = './data/clips/Patient_%d/'%(i+1)
    fns = [fn for fn in os.listdir(dn) if '.mat' in fn]
    
    feat = []
    labels = []
    
    for fn in fns:
        m = loadmat(dn+fn)
        d = m['data']
        d = resample(d, int(d.shape[1]*500/m['freq']), axis=1)
        s = array([stft(x,64,8) for x in d])
        s = sqrt(real(s)**2+imag(s)**2)    
        if 'inter' in fn:
            l = 0
        elif '_ictal' in fn:
            l = 1
        else:
            l = -1
        labels.append(l)
        feat.append(s.flatten())    
        
    feat = array(feat)
    labels = array(labels)
    fns = array(fns)
    
    d = {}
    d['feat'] = feat
    d['labels'] = labels
    d['feat_norm'] = (feat-feat.mean(0))/feat.std(0)
    d['fns'] = fns
    
    savez('./data/data_Patient_%d'%(i+1), **d)


# In[ ]:

####################################################################
## Make predictions with logistic regression on spectrogram features
####################################################################
preds = []

for i in range(4):
    print 'Dog %d'%(i+1)
    sys.stdout.flush()    
    d = load('./data/data_Dog_%d.npz'%(i+1))
    
    labels = d['labels'][d['labels'] != -1]
    data_normed = d['feat_norm'][d['labels'] != -1]
    
    data_test = d['feat_norm'][d['labels']==-1]
    fn_test = d['fns'][d['labels']==-1]    

    clf = sklearn.linear_model.LogisticRegression()
    clf.fit(data_normed, labels)
    print clf.score(data_normed, labels)
    
    pred_test = clf.predict_proba(data_test)
    for fn, p in zip(fn_test, pred_test):
        preds.append('%s,%f,0'%(fn,p[1]))
        
for i in range(8):
    print 'Patient %d'%(i+1)
    sys.stdout.flush()    
    d = load('./data/data_Patient_%d.npz'%(i+1))
    
    labels = d['labels'][d['labels'] != -1]
    data_normed = d['feat_norm'][d['labels'] != -1]
    
    data_test = d['feat_norm'][d['labels']==-1]
    fn_test = d['fns'][d['labels']==-1]    

    clf = sklearn.linear_model.LogisticRegression()
    clf.fit(data_normed, labels)
    print clf.score(data_normed, labels)
    
    pred_test = clf.predict_proba(data_test)
    for fn, p in zip(fn_test, pred_test):
        preds.append('%s,%f,0'%(fn,p[1]))
        
of = open('submission.csv','w')
of.write('clip,seizure,early\n')
for p in preds:
    of.write(p+'\n')
of.close()


# In[ ]:

###################################################################
## Generate maximal cross-correlation features
###################################################################

import os, sys
from scipy.io import loadmat
import numpy as np
from scipy.signal import resample

from joblib import Parallel, delayed

def max_ccors(X):
    N = np.correlate(np.ones(X.shape[1]),np.ones(X.shape[1]),'full')
    ac = [np.correlate(X[i],X[i])[0] for i in range(X.shape[0])]
    mcc = []
    for i in range(X.shape[0]):
        for j in range(i):
            cc = np.correlate(X[i],X[j],'full')/N*1./np.sqrt(ac[i]*ac[j]+1e-6)
            mcc.append(cc.max())
    return np.array(mcc)

def process(fn):
    m = loadmat(fn)
    d = m['data']
    d = resample(d, int(d.shape[1]*500/m['freq']), axis=1)
    s = max_ccors(d)
    if 'inter' in fn:
        l = 0
    elif '_ictal' in fn:
        l = 1
    else:
        l = -1
    return l, s

for i in range(4):
    print i
    sys.stdout.flush()
    dn = './data/clips/Dog_%d/'%(i+1)
    fns = [fn for fn in os.listdir(dn) if '.mat' in fn]
    
    results = Parallel(n_jobs=-1)(delayed(process)(dn+fns[i]) for i in xrange(len(fns)))
    labels = np.array([x[0] for x in results])
    feat = np.array([x[1] for x in results])
    fns = np.array(fns)
    
    d = {}
    d['labels'] = labels
    d['feat_norm'] = (feat-feat.mean(0))/feat.std(0)
    d['fns'] = fns
    d['feat'] = feat
    
    np.savez('./data/data_ccor_Dog_%d'%(i+1), **d)
    
for i in range(8):
    print i
    sys.stdout.flush()
    dn = './data/clips/Patient_%d/'%(i+1)
    fns = [fn for fn in os.listdir(dn) if '.mat' in fn]
    
    results = Parallel(n_jobs=-1)(delayed(process)(dn+fns[i]) for i in xrange(len(fns)))
    labels = np.array([x[0] for x in results])
    feat = np.array([x[1] for x in results])
    fns = np.array(fns)
    
    d = {}
    d['labels'] = labels
    d['feat_norm'] = (feat-feat.mean(0))/feat.std(0)
    d['fns'] = fns
    d['feat'] = feat
    
    np.savez('./data/data_ccor_Patient_%d'%(i+1), **d)


# In[ ]:

#############################################################################
## Prepare files for scattering coefficient extraction (using Matlab/Octave)
#############################################################################

import sys, os
import scipy.io
import scipy.signal

for i in range(1,5):
    dn = 'Dog_%d'%i
    fns = os.listdir('./data/clips/'+dn)
    ifns = ['./data/clips/'+dn+'/'+fn for fn in fns]
    ofns = ['./data/scattering/'+dn+'/'+fn for fn in fns]

    scipy.io.savemat('./data/scattering/'+dn+'_fns.mat',{'ofns':ofns,'ifns':ifns})

for i in range(1,9):
    dn = 'Patient_%d'%i
    try:
        os.makedirs('./data/resample/'+dn)
    except:
        continue
    print dn
    sys.stdout.flush()
    fns = os.listdir('./data/clips/'+dn)
    ifns = ['./data/resample/'+dn+'/'+fn for fn in fns]
    ofns = ['./data/scattering/'+dn+'/'+fn for fn in fns]
    
    for fn in fns:
        data = scipy.io.loadmat('./data/clips/'+dn+'/'+fn)['data']
        data = scipy.signal.resample(data,400,axis=1)
        scipy.io.savemat('./data/resample/'+dn+'/'+fn, {'data':data})
    scipy.io.savemat('./data/scattering/'+dn+'_fns.mat',{'ofns':ofns,'ifns':ifns})


# Matlab code to calculate scattering coefficients from the clips
# ```
# % Filter length
# N = 400;
# 
# % Calculate coefficients with averaging scale of 8192 samples (~370 ms @
# % 22050 Hz sampling rate).
# T = 40;
# 
# % First-order filter bank with 8 wavelets per octave. Second-order filter bank
# % with 1 wavelet per octave.
# filt_opt.Q = [8 1];
# % Calculate maximal wavelet scale so that largest wavelet will be of bandwidth 
# % T.
# filt_opt.J = T_to_J(T, filt_opt);
# 
# % Only calculate scattering coefficients up to second order.
# scat_opt.M = 2;
# 
# % Prepare wavelet transforms to be used for scattering.
# Wop = wavelet_factory_1d(N, filt_opt, scat_opt);
# 
# feature_fun = @(x)(format_scat(log_scat(renorm_scat(scat(x, Wop)))));
# 
# 
# mats = {'Dog_3_fns.mat', 'Dog_4_fns.mat', 'Patient_1_fns.mat','Patient_2_fns.mat','Patient_3_fns.mat','Patient_4_fns.mat','Patient_5_fns.mat','Patient_6_fns.mat','Patient_7_fns.mat','Patient_8_fns.mat'}
# 
# for j=2:2
#     load(mats{j})
# 
#     for i=1:length(ifns)
#         load(strtrim(ifns(i,:)));
#         d = reshape(transpose(data),400,1,16);
#         f = feature_fun(d);
#         save(strtrim(ofns(i,:)), 'f'); 
#     end
# end
# ```

# In[ ]:



