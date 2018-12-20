
#make nicer histograms:
#https://stackoverflow.com/questions/26218704/matplotlib-histogram-with-collection-bin-for-high-values

#normalize histogram:
#https://stackoverflow.com/questions/5498008/pylab-histdata-normed-1-normalization-seems-to-work-incorrect


#import itertools as it
#import multiprocessing
#from matplotlib import gridspec


from parametersgs import *
from load_fano_from_sim import *

import numpy as np

#fanoLrpB=np.loadtxt('/Users/sdebuyl/Dropbox/fanoLrpB.txt')
#fano3DS=np.loadtxt('/Users/sdebuyl/Dropbox/fano3DS.txt')
#fano2DS=np.loadtxt('/Users/sdebuyl/Dropbox/fano2DS.txt')
#fanoMDS=np.loadtxt('/Users/sdebuyl/Dropbox/fanoMDS.txt')
#np.savetxt('/Users/sdebuyl/fanoLrpB.txt',fanoLrpB)



fano1_name='LrpB'
print(fano1_name)

fano1=load_ff(mymodel='3',oscillators=DDDS,file='/Users/sdebuyl/Dropbox/stoch_LrpB/')

fano2_name='3DS'
print(fano2_name)
#fanoMDS = load_ff(mymodel='M', oscillators=MDS, file='/Users/sdebuyl/stoch_MDS/')
fano3DS = load_ff(mymodel='3', oscillators=DDDS, file='/Users/sdebuyl/stoch_3DS/')
fano2=fano3DS


bins=np.linspace(0,20,100)
bins=np.append(bins,np.max([np.max(fano1),np.max(fano2)]))
if 1==1:
    fig = plt.figure()
    #plt.hist(fanoMDS, bins=[-1.0,1,2,10,max(fanoMDS)],alpha=.3,label='MDS')
    # #plt.hist(fano2DS, bins='auto',alpha=.3,label='2DS')
    plt.hist(fano1, bins=bins,alpha=.3,label=fano1_name,normed=True)
    plt.xlim([0,20])
    plt.xlabel('fano factor')
    plt.ylabel('freq. of solutions')
    plt.hist(fano2, bins=bins,alpha=.3,label=fano2_name,normed=True)
    plt.legend(loc='upper right')
    plt.savefig('hist.pdf')

#print('lrpb',fanoLrpB.shape)
print(fano1_name,fano1.shape)
#print('2ds',fano2DS.shape)
print(fano2_name,fano2.shape)
print(fano1_name+' bigger than 2',len(np.where(fano1>1.2)[0])/len(fano1))
print(fano2_name+'bigger than 2',len(np.where(fano2>1.2)[0])/len(fano2))
#
#
# fano1rest=fano1[0:800]
# fano2rest=fano2[0:800]
# print(fano1_name+' bigger than 1 restricted',len(np.where(fano1rest>1.5)[0])/len(fano1rest))
# print(fano2_name+' bigger than 1 restricted',len(np.where(fano2rest>1.5)[0])/len(fano2rest))