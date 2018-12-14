
def load_ff(mymodel,oscillators,file):
    #### extract data 2DS


    vartoplot=-1



    #load LrpB 3DS data
    #file = '/Users/sdebuyl/stoch_2DS/'  # '/Users/sdebuyl/sgo/stoch_' + mymodel + 'DS_LrpB/'

    LrpB_shift = 0  # this variable was introduced to be able to change label the position of numbers keeping track of different simulations
    namefiletosavedata = 'stoch_' + mymodel + 'DS'
    # actually import the data:
    mylist = ospack.listdir(file)
    # create list for each type of non-monotonicity
    list0 = list()
    list1 = list()
    list2 = list()
    list3 = list()
    # create a counter for spiky solution
    counter_spikes_LrpB=0
    for kk in mylist:

        if kk.find('stoch_' + mymodel + 'DS_parms_') != -1:
            end = kk.find('.txt')
            #print('looking', kk[16 + LrpB_shift:17 + LrpB_shift], int(kk[18 + LrpB_shift:end]))
            if kk[16 + LrpB_shift:17 + LrpB_shift] == '0':
                list0.append(int(kk[18 + LrpB_shift:end]))
            if kk[16 + LrpB_shift:17 + LrpB_shift] == '1':
                list1.append(int(kk[18 + LrpB_shift:end]))
            if kk[16 + LrpB_shift:17 + LrpB_shift] == '2':
                list2.append(int(kk[18 + LrpB_shift:end]))
            if kk[16 + LrpB_shift:17 + LrpB_shift] == '3':
                list3.append(int(kk[18 + LrpB_shift:end]))

    lists = [list0, list1, list2, list3]

    print('len list 0', len(list0))
    print('len list 1', len(list1))
    print('len list 2', len(list2))
    print('len list 3', len(list3))
    tot_simulations_LrpB=len(list0)+len(list1)+len(list2)+len(list3)

    # we will plot the fano factor as a function of the mean value of the dimer concentration
    # we first create empty arrays to be filled with means and variances of each simulation
    mymeans0 = np.zeros(len(list0))
    mymeans1 = np.zeros(len(list1))
    mymeans2 = np.zeros(len(list2))
    mymeans3 = np.zeros(len(list3))

    mymeans_wt = [mymeans0, mymeans1, mymeans2, mymeans3]

    myvars0 = np.zeros(len(list0))
    myvars1 = np.zeros(len(list1))
    myvars2 = np.zeros(len(list2))
    myvars3 = np.zeros(len(list3))

    myvars_wt = [myvars0, myvars1, myvars2, myvars3]



    myfanofactor3=np.array([])

    jj=1

    for k in range(len(lists[jj])):  #np.random.choice(len(lists[jj]),160)
        mymeans_wt = np.loadtxt(file + namefiletosavedata + '_means_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
        myvars_wt = np.loadtxt(file + namefiletosavedata + '_vars_' + str(jj) + '_' + str(lists[jj][k]) + '.txt')
        temp = mymeans_wt[vartoplot]
        if temp!=0:
            myfanofactor3=np.append(myfanofactor3,myvars_wt[vartoplot] / mymeans_wt[vartoplot])



    ##############
    ##############
    ##############


    return myfanofactor3