
# coding: utf-8

# In[ ]:

import joblib, os, sys, re, pickle


# In[ ]:

def warp(x):
    """
        Logarithmically rescale values with p < 0.1 and (1-p) < 0.1 to avoid loss of precision
    """
    x = copy(x)
    low = x<0.1
    high = x>0.9

    xl = x[low]
    pmin = log10(xl[xl!=0].min())-1
    lx = log10(xl)
    lx[isinf(lx)] = pmin
    lx = (lx+1)*(3./(lx.max()-lx.min()))-1
    xl = 10**lx
    x[low] = xl
    
    xh = x[high]
    pmin = log10(1-xh[xh!=1].max())-1
    lx = log10(1-xh)
    lx[isinf(lx)] = pmin
    lx = (lx+1)*(1./(lx.max()-lx.min()))-1
    xh = 1-10**lx
    x[high] = xh

    return x


# In[ ]:

labels = ['patient_7']

of = open('./submission.csv','w')
of.write('clip,seizure,early\n')

lines = 0
for label in labels:
    fns = [fn for fn in os.listdir('./output/') if label in fn]
    print "Averaging %d predictions for subject %s" % (len(fns), label)
    if len(fns) == 0:
        print "Error: missing predictions for subject %s" % label
        continue
        
    preds = vstack([warp(pickle.load(open('./output/'+fn))['pred']) for fn in fns])
    d = pickle.load(open('./output/'+fn))
    fns = array(d['fns'])
    predfns = fns[d['y']==-1]
    p = preds.mean(0)
    a = dict(zip(fns,p))

    for fn in predfns:
        of.write('%s,%f,%f\n'%(fn,a[fn],a[fn]))
        lines += 1
of.close()

if lines != 32915:
    print "Warning: incorrect number of predictions in submission"

