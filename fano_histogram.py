
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



fano1_name='MDSn'
print(fano1_name)
fanoMDSn = load_ff(mymodel='M',oscillators=MDS,file='/Users/sdebuyl/stoch_MDSn/')
fano1=fanoMDSn

fano2_name='3DS'
print(fano2_name)
#fanoMDS = load_ff(mymodel='M', oscillators=MDS, file='/Users/sdebuyl/stoch_MDS/')
fano3DS = load_ff(mymodel='3', oscillators=DDDS, file='/Users/sdebuyl/sgo/stoch_3DS/')
fano2=fano3DS

if 1==9:
    fanoLrpB=load_ff(mymodel='3',oscillators=DDDS,file='/Users/sdebuyl/Dropbox/stoch_LrpB/')
    fano2DS=load_ff('2',DDS,'/Users/sdebuyl/stoch_2DS/')
    fanoMDS=load_ff('M',MDS,'/Users/sdebuyl/stoch_MDS/')
    fanoLrpB = load_ff(mymodel='3', oscillators=DDDS, file='/Users/sdebuyl/Dropbox/stoch_LrpB/')

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
print(fano1_name+' bigger than 2',len(np.where(fano1>2.0)[0])/len(fano1))
print(fano2_name+'bigger than 2',len(np.where(fano2>2.0)[0])/len(fano2))

#print('lrpB bigger between 1 and 10',len(np.where(np.logical_and(2<fanoLrpB, fanoLrpB<10.0))[0])/len(fanoLrpB))
#print('3DS bigger ',len(np.where(np.logical_and(2<fano3DS, fano3DS<10.0))[0])/len(fano3DS))

# bins= [-1.0,1,2,10,max(fano3DS)]
# data = fano3DS
# hist, bin_edges = np.histogram(data,bins) # make the histogram
#
# fig,ax = plt.subplots()
#
# # Plot the histogram heights against integers on the x axis
# ax.bar(range(len(hist)),hist,width=1)
#
# # Set the ticks to the middle of the bars
# ax.set_xticks([0.5+i for i,j in enumerate(hist)])
#
# # Set the xticklabels to a string that tells us what the bin edges were
# ax.set_xticklabels(['{} - {}'.format(bins[i],bins[i+1]) for i,j in enumerate(hist)])
#
# plt.savefig('hist2DS.pdf')


