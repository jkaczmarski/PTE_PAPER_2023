import sys
import matplotlib.pyplot as plt
import numpy as np
import os 
import subprocess
import matplotlib as mtl
import glob
sys.path.append('/Users/u4855540/Desktop/PTE_Data/ER/')
import PyEnsembleRefinement as er 


def split_chains(list_of_tuples):
    resids = [int(list_of_tuples[i][0]) for i in range(len(list_of_tuples))]
    breaks = []
    for i in range(len(resids)):
        if (len(str(resids[i-1]))!=len(str(resids[i]))):
            if (resids[i]!=100) and i!=0 :
                breaks.append(i)
    breaks.append(-1)
    
    chains = []
    for i in range(len(breaks)):
        if i==0:
            chain = list_of_tuples[:breaks[i]]

        else:
            chain = list_of_tuples[breaks[i-1]:breaks[i]]
        chains.append(chain)
    
    return chains

def processrmsf_1(xvg_file):
    """Takes an input xvg file (string) with rmsfs; splits chains using split_chains(); returns RMSF data by chain
    infile(xvgfile) --> list of lists([chain[[residue, rmsf]],])"""

    data = open(xvg_file)
    data = data.readlines()
    data = [i.strip() for i in data]
    data = [i for i in data if ('@' not in i) and ('#' not in i)]
    proc_data = []

    for i in data: 
        r = i.split()
        proc_data.append(r)
    
    chain_data = split_chains(proc_data)
    
    nan_chains = []
    
    
    ####THE FOLLOWING CODE IS JUNK######
    for i in chain_data: 
        nan_chain = []
        chain = np.array(i)
        resids = chain[:,0]
        rmsfs  = chain[:,1]
        
        for i in range(len(resids)): 
            if i != 0:
                current =  (int(resids[i]))
                previous =  (int(resids[i-1]))

                previous_c = previous+1
                if current == previous_c:
                    nan_chain.append((resids[i], rmsfs[i]))
                elif current != previous_c:
                    print (current)
                    difference = current-previous
                    print (difference)
                    for k in range(difference-1):
                        nan_chain.append(('nan', 'nan'))
                    nan_chain.append((resids[i], rmsfs[i]))
        nan_chains.append(np.array(nan_chain))
 ########################################################
 
    return chain_data


def processrmsf_2(inlist):
    """Takes the raw data from processmsf_1 and transforms it into numerical array-like data.
    list of lists([chain[[residue, rmsf]],]) --> [[[resid],[rmsfs]]chain]"""
    data_series=[]
    for i in inlist:
        t_array = np.array(i)
        x_data=(list(t_array[:,0]))
        y_data=(list(t_array[:,1]))
        x_data=[float(i) for i in x_data]
        y_data=[float(i) for i in y_data]

        series = [x_data,y_data]

        data_series.append(series)
        
    return data_series

def processrmsf_3(inlist, eliminate=False):

    """If a residue position exists in both chains, this returns a 'tuple' like [resid, rmsfchain1, rmsfchain2].
    If residue position doesn't exist in both chains, this is omitted."""
    chains_resids = []
    for i in inlist:
        chains_resids.append(i[0])

    if eliminate == True:
        average_resids = list(set.intersection(*map(set,chains_resids)))
    else: 
        average_resids = list(set.union(*map(set,chains_resids)))


    print (average_resids)
    print(inlist)
    to_average_data = []
    
    if eliminate:
    
        for resid in average_resids:
            res_data_to_average = [resid]
            for chain in inlist:
                data_index = chain[0].index(resid)
                rmsf_dat = chain[1][data_index]
                res_data_to_average.append(rmsf_dat)
            to_average_data.append(res_data_to_average)

        return to_average_data 
    
    else: 
        for resid in average_resids: 
            res_data_to_average = [resid]
            for chain in inlist: 
                if resid in chain[0]:
                    data_index = chain[0].index(resid)
                    rmsf_dat = chain[1][data_index]
                    res_data_to_average.append(rmsf_dat)
                    
                else:
                    res_data_to_average.append(float('nan'))
                    
            to_average_data.append(res_data_to_average)
            print (res_data_to_average)
        return to_average_data
                    
           
        
        
        

def processrmsf_4(inlist):
    averaged_data = []
    for i in inlist:
        res_id = i[0]
        rmsf_data = i[1:]
        av_rmsf = np.nanmean(rmsf_data)
        std_rmsf = np.nanstd(rmsf_data)
        av_list = [res_id,av_rmsf, std_rmsf]
        #print (av_list)
        averaged_data.append(av_list)
    return averaged_data

def analyseXVG(xvgfile, elim=False):
    return np.array(processrmsf_4(processrmsf_3(processrmsf_2(processrmsf_1(xvgfile)), eliminate=elim)))


def analyseXVG_notav(xvgfile, elim=False):
    return np.array(processrmsf_3(processrmsf_2(processrmsf_1(xvgfile)), eliminate=elim))




def pooledSTD(array):
    return np.sqrt(np.sum(np.square(array))/len(array)) 

def find_min_Rfrees(parentdir, best_nums):
    import glob
    os.chdir(parentdir)
    print(parentdir)
    data_list = []
    rep_dirs = [i for i in os.listdir() if os.path.isdir(i) and ('replicate' in i) ]
    for i in rep_dirs:
        print (i)
        os.chdir(str(i))
        logfile = glob.glob("*.log")
        logfile = open(str(logfile[0]), 'r')
        lines = logfile.readlines()
        for j in lines: 
            if ('FINAL' in j) and ('Rfree' in j):
                r_free = j.split(' ')[-4]
                print(r_free)
                r_work = j.split(' ')[3]
                print (r_work)
                
            elif 'Ensemble size : ' in j:
                print (j)
                ens_size = int(j.split(':')[-1])
            
            

        rep_name = i
        out_tuple = (i,r_free,r_work, ens_size)
        data_list.append(out_tuple)
        os.chdir(parentdir)
        r_frees = [data_list[i][1] for i in range(len(data_list))]
        repl_name = [data_list[i][0] for i in range(len(data_list))]

    #print (data_list)
    
    data_array = np.array(data_list)
    rpls = data_array[:,0]
    rfrs = data_array[:,1]
    models = data_array[:,3]
    rwrk = data_array[:,2]
    
    idx = np.argsort(rfrs)
    
    best_rpls = [rpls[i] for i in idx]
    print(parentdir)
    print(i)
    print(best_rpls)
    best_rfrs = [rfrs[i] for i in idx]
    best_models = [models[i] for i in idx]
    best_rwrk = [rwrk[i] for i in idx]
    #print (best_rpls)
    #print (best_rfrs)
    return (best_rpls[:best_nums], best_rfrs[:best_nums],best_rwrk[:best_nums],
            best_models[:best_nums])

def grandRMSF(dirlist):
    """All subdirs in dirlist, computes the mean and standard devation in RMSF for each residue from both chains
    of each replicate; then computes the grand mean and grand standard deviation, and returns this figure, 
    by residue."""
    replicate_resids = []
    replicate_rmsfs  = []
    replicate_stds   = []
    for i in dirlist:
        xvg_data = analyseXVG(i+'/rmsf.xvg')
        
        replicate_resids.append(xvg_data[:,0])
        replicate_rmsfs.append(xvg_data[:,1])
        replicate_stds.append(xvg_data[:,2])

    usable_resids = list(set.union(*map(set,replicate_resids)))
    usable_data = []
    for i in range(len(replicate_resids)):
        data_list = []
        for j in range(len(replicate_resids[i])):
            if replicate_resids[i][j] in usable_resids:
                data_list.append([replicate_resids[i][j],replicate_rmsfs[i][j], replicate_stds[i][j]])
        usable_data.append(data_list)
    usable_data = np.array(usable_data)            
    restructured_means = []
    restructured_stds  = []
    
    for replicate in usable_data:
        means  = replicate[:,1]
        stdevs = replicate[:,2]
        restructured_means.append(means)
        restructured_stds.append(stdevs)
        
    
    grand_means   = np.nanmean(np.array(restructured_means), axis=0)
    pooled_stdevs = np.apply_along_axis(pooledSTD, 0,np.array(restructured_stds))
            
    
    
    return (usable_data, grand_means, pooled_stdevs)

def doAverageAnalysis(parentdir):
    os.chdir(parentdir)
    best_dirs = find_min_Rfrees(parentdir,5)
    
    print (best_dirs)

    return grandRMSF(best_dirs[0])
    

def getRefineParams(ensemble_dir):

    os.chdir(ensemble_dir)
    os.chdir('replicate_1/')
    data_list = []
    logfile = glob.glob("*.log")
    logfile = open(str(logfile[0]), 'r')
    lines = logfile.readlines()
    for j in lines: 
        if ('FINAL' in j) and ('Rfree' in j):
            print (j)
            r_free = float(j[28:36])
        elif 'tx =' in j and '_tx' not in j:
            print (j)
            tx = float(j.split('=')[-1])
        elif 'wxray_coupled_tbath_offset' in j:
            print(j)
            wxray = float(j.split('=')[-1])
        elif 'ptls =' in j and '_ptls' not in j:
            print(j)
            ptls = float(j.split('=')[-1])
    

  
    out_tuple = [ ensemble_dir.split('/')[-2],ptls, wxray,tx]
  
    os.chdir(ensemble_dir)
    return out_tuple


def get_mean_Rfrees(dirlist):
	mean_rFrees = []
	for d in dirlist:
	    data =  find_min_Rfrees(d,5)
	    rfree_str = data[1]
	    rfrees = [float(i) for i in rfree_str]
	    mean_rfree = np.mean(rfrees)
	    stdev_refree = np.std(rfrees)
	   
	    rwork_str = data[2]
	    rwrks = [float(i) for i in rwork_str]
	    mean_rwrk = np.mean(rwrks)
	    std_rwrk = np.std(rwrks)
	    
	    
	    ens_size = data[3]
	    ens_size = [int(i) for i in ens_size]
	    mean_size = np.mean(ens_size)
	    std_size  = np.std(ens_size)
	    mean_rFrees.append((d.split('/')[-2],mean_rfree,stdev_refree,mean_rwrk, std_rwrk,mean_size,std_size ))
	    
	mean_rFrees=np.array(mean_rFrees)
	return mean_rFrees

	
