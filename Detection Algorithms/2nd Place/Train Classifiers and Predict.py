
# coding: utf-8

# In[1]:

import os, sys, itertools

import pickle

import joblib

import sklearn.cross_validation
import sklearn.metrics


# In[2]:

from simplednn.nets import *


# In[3]:

def oprint(s):
    os.system('echo "%s"'%s)

def train_net(label, w, n):
    """
        Given a sample label and filter set label w, train a neural network and make predictions. 
        Output is saved to a pickle indexed by parameter n.
    """
    ff = open('./data/cov_opt_%s_%s.pickle'%(label,w),'rb')    
    d = pickle.load(ff)
    fff = np.load(ff)

    X = d['covs'].astype('float32')
    y = d['y'].astype('int32')
   
    N_CHANNELS = 12

    # Select N_CHANNELS at random from the full covariance matrices
    channelidx = range(X.shape[2])
    numpy.random.shuffle(channelidx)
    channelidx = channelidx[:N_CHANNELS]
    X = X[:,:,channelidx,:][:,:,:,channelidx]
    X = X.reshape((X.shape[0],-1))
    
    # Partition the labeled data and split into train/test
    X_traintest = X[y!=-1]
    y_traintest = y[y!=-1]

    test_size = 0.1
    rseed = randint(2**32)
    tri, tei = next(iter(sklearn.cross_validation.ShuffleSplit(
                         X_traintest.shape[0], n_iter=1, test_size=test_size, random_state=rseed)))
    X_train = X_traintest[tri]
    X_test = X_traintest[tei]
    y_train = y_traintest[tri]
    y_test = y_traintest[tei]
    
    oprint("%s, %s"%(label, w))
    oprint("Selecting EEG channels: %s"%channelidx)
    oprint("Total dataset size:")
    oprint("n labels: %d" % y_traintest.shape[0])
    oprint("n samples: %d" % X_traintest.shape[0])
    oprint("n features: %d" % X_traintest.shape[1])
    oprint("n classes: %d" % len(set(y_traintest)))
    
    # Initialize and train network
    dnn1 = DropoutNet(numpy_rng=numpy.random.RandomState(rseed), n_ins=X_traintest.shape[1],
                    layers_types=[ReLU, ReLU, LogisticRegression],
                    layers_sizes=[200, 100],
                    trainables=[1, 1, 1, 1, 1, 1],
                    dropout_rates=[0., 0.5, 0.5],
                    # TODO if you have a big enough GPU, use these:
                    #layers_types=[ReLU, ReLU, ReLU, ReLU, LogisticRegression],
                    #layers_sizes=[2000, 2000, 2000, 2000],
                    #dropout_rates=[0., 0.5, 0.5, 0.5, 0.5],
                    n_outs=2,
                    max_norm=4.,
                    fast_drop=True,
                    debugprint=0)

    dnn1.fit(X_train, y_train, max_epochs=100, method='adadelta', verbose=False, test=(X_test is not None), plot=False, save=None)
    
    # Test on CV data. Fails if no seizure clips in split.
    try:
        cv = sklearn.metrics.roc_auc_score(y_test, dnn1.predict_proba(X_test)[:,1])
        os.system('echo "%s"'%('%s %s CV AUC: %f'%(label, w, cv)))
    except:
        pass
    
    # Make predictions on full dataset and save to pickle
    pred = dnn1.predict_proba(X)[:,1]
    d = {'fns':d['fns'], 'y':d['y'], 'pred':pred, 'rseed':rseed}
    pickle.dump(d, open('./output/pred_%d_%s_%s'%(n, label, w),'w'))


# In[ ]:

labels = ['patient_0']

ws = 'ABC'

# Train multiple networks on each processed subject dataset (50 x 12 x 3)
# This will take several hours
r = joblib.Parallel(n_jobs=-1)(joblib.delayed(train_net)(label, w, n) for n, label, w in 
                              itertools.product(range(50),array(labels),ws))



# In[ ]:



