
import matplotlib.pylab as plt
import numpy as np


fanoLrpB=np.loadtxt('/Users/sdebuyl/Dropbox/fanoLrpB.txt')
fano3DS=np.loadtxt('/Users/sdebuyl/Dropbox/fano3DS.txt')
fano2DS=np.loadtxt('/Users/sdebuyl/Dropbox/fano2DS.txt')
fanoMDS=np.loadtxt('/Users/sdebuyl/Dropbox/fanoMDS.txt')



fig = plt.figure()
 #plt.hist(fanoMDS, bins=[-1.0,1,2,10,max(fanoMDS)],alpha=.3,label='MDS')
# #plt.hist(fano2DS, bins='auto',alpha=.3,label='2DS')
plt.hist(fano3DS, bins='auto',alpha=.3,label='3DS',normed=True)
plt.xlim([0,20])
plt.hist(fanoLrpB, bins='auto',alpha=.3,label='LrpB',normed=True)
plt.legend(loc='upper right')
plt.savefig('hist.pdf')

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

print('lrpb',fanoLrpB.shape)
print('3ds',fano3DS.shape)