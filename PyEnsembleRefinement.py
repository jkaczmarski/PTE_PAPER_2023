#Mahakaran Sandhu, RSC, 2019#

import matplotlib.pyplot as plt
import os
import numpy as np
import subprocess
import random 
import multiprocessing
import glob
import shutil
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mtl

"""
This is a Python wrapper to the Ensemble Refinement (ER) protocol in PHENIX. It is indended to be used with Jupyter Notebooks. 
The motivation behind this project was to streamline the ER protocol. Doing everything manually is very tedious and prone to mistakes. 
This project is designed to cut setup and analysis times and to improve consistency. The package calls on a number of dependencies,
phenix.enseble_refine being the obvious one. Other dependencies are: GROMACS (to calculate RMSF), and grep. The package is built around
command-line tools so it often communicates with shells via the Python package subprocess. The multiprocessing functionality allows the 
user to run multiple simulations across multiple cores from the same Jupyter Notebook, which is much better than manually running simulations
in many different shells.

"""




 
def runEnsembleRefinement(dirname):
    print ('Now running:'+str(dirname))
    os.system('cd '+ str(dirname) + '&& sh *sh')
    print ('Finished running:'+str(dirname))

def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out



###########################################################GRID SEARCH SETUP AND ANALYSIS MODULES#############################################


def setupGridSearch(parentdir,modelpath, datapath, cifpath,ligandNames):
    """ Sets up a 45-subdirectory grid-search for the optimal values of the parameters pTLS, tx and w-xray. Values for these are
    hard-coded, but can be changed by changing this file. 
    parentdirectory(type:str)--> writes directories with appropriate bash (.sh) run scripts
    For modelpath(pdb), datapath(.mtz) and cifpath(ligand restraint file), its best to write the absolute path rather than the relative 
    path. It's recommended that the absolute path for the pdb and the mtz file is in the same directory, meaning you'll need to copy these
    files into the generated subdirectories (this function will do that, so ensure that you have these files in the 
    parentdir). ligandNames is a list (type list) with strings that denote names of
    the ligands that you need to restrain harmonically."""

    #initialise pTLS
    p = [0.6, 0.7, 0.8, 0.9, 1.0]
    #initialise w_xray_weight
    w = [2.5, 5, 10]
    #initialise t_x
    t = [0.3, 0.6, 1]


    #generate possible combinations of above parameters (45 different combinations)
    paramCombinations = []
    for i in p: 
        for j in w:
            for k in t: 
                out_tuple = (i,j,k)
                paramCombinations.append(out_tuple)

    
    os.chdir('paramOptimisation')


    restr_sel = ['resname '+i +' or' if ligandNames.index(i)!=ligandNames.index(ligandNames[-1]) else 'resname '+i for i in ligandNames]
    restr_str = ' '.join(restr_sel)
    
    



    for i in paramCombinations:
        dirname = str(i[0]) + '_'+ str(i[1]) + '_' + str(i[2])
        outfilename = 'refine_'+str(i[0]) + '_'+ str(i[1]) + '_' + str(i[2]) + '.sh'
        os.mkdir(dirname)
        filepath = os.path.join(dirname,outfilename)
        subprocess.call('cp ../*.pdb . && cp ../*.mtz .', shell=True) 
        outfile = open(filepath, 'w+')
        outfile.write('#!/bin/bash \n \n')
        outfile.write('p=' + str(i[0]) + '\n')
        outfile.write('w=' + str(i[1]) + '\n')
        outfile.write('t=' + str(i[2]) + '\n \n')
        outfile.write('model='+modelpath+' \n')
        outfile.write('data='+datapath+' \n')
        outfile.write('cif='+cifpath+' \n \n')
        outfile.write("phenix.ensemble_refinement $model $data $cif harmonic_restraints.selections='"+str(restr_str)+"' ptls=$p wxray_coupled_tbath_offset=$w tx=$t  nproc=1")
        outfile.close()



def runGridSearch(parentdir, cores):
    """ Runs the grid search protocol in a parentdir that has been setup with setupGridSearch. Use cores number of logical cores. """


    replicate_dirs = [i for i in os.listdir(parentdir) if os.path.isdir(i) and '0' in i]
    completedER = []
    
    
    x = len(replicate_dirs)/cores
    
    
    
    processedDirList = chunkIt(replicate_dirs, x)
    print (processedDirList)

    for i in processedDirList:
        procs = []
        for j in i:
            print(j)
            proc = multiprocessing.Process(target=runEnsembleRefinement, args=(j,))
            procs.append(proc)
            proc.start()
        for k in procs:
            k.join()



























def analyseGridSearch(parentdir):
    """Analyse the ensemble refinement grid search."""
    os.chdir(parentdir)
    data_list = []
    dirs = [i for  i in os.listdir(parentdir) if os.path.isdir(i) and '0' in i]
    for i in dirs:
        print (i)
        os.chdir(str(i))
        logfile = glob.glob("*.log")
        logfile = open(str(logfile[0]), 'r')
        lines = logfile.readlines()
        for j in lines: 
            if ('FINAL' in j) and ('Rfree' in j):
                print (j)
                r_free = float(j[28:36])
        parameter_data = i.split('_')
        ptls = float(parameter_data[0])
        w = float(parameter_data[1])
        tx = float(parameter_data[2])
        out_tuple = (r_free, ptls, w,tx, i)
        data_list.append(out_tuple)
        os.chdir(parentdir)

    r_free = [i[0] for i in data_list ]
    rf_min_ind = r_free.index(min(r_free))
    ptls = [i[1] for i in data_list]
    w = [i[2] for i in data_list]
    tx = [i[3] for i in data_list]
    min_params = (ptls[rf_min_ind],tx[rf_min_ind], w[rf_min_ind])
    print('Best pTLS: '+str(min_params[0])+'\n Best tx:'+str(min_params[1])+'\n Best w-xray:' + str(min_params[2]))
    return min_params













############################################################################################REPLICATE RUNNING MODULES##############################################################################################

def replicates_setup(parentdir, p_opt, w_opt, t_opt, modelpath, datapath, cifpath, ligandNames, nreps=10):
    """Ligand names is a list of names of ligands to restrain."""

    replicate_dirs = list(range(1,nreps+1))

    restr_sel = ['resname '+i +' or' if ligandNames.index(i)!=ligandNames.index(ligandNames[-1]) else 'resname '+i for i in ligandNames]
    restr_str = ' '.join(restr_sel)
    
    



    for i in replicate_dirs:
        random_seed = random.randint(1000000, 9000000)
        os.chdir(parentdir)
        dirname = 'replicate_' + str(i)
        outfilename = 'refine_' + str(i) +  '.sh'
        os.mkdir(dirname)
        filepath = os.path.join(dirname,outfilename)
        outfile = open(filepath, 'w+')
        outfile.write('#!/bin/bash \n \n')
        outfile.write('p=' + str(p_opt) + '\n')
        outfile.write('w=' + str(w_opt) + '\n')
        outfile.write('t=' + str(t_opt) + '\n \n')
        outfile.write('model='+modelpath+' \n')
        outfile.write('data='+datapath+' \n')
        outfile.write('cif='+cifpath+' \n \n')
        outfile.write("phenix.ensemble_refinement $model $data $cif harmonic_restraints.selections='"+str(restr_str)+"' ptls=$p wxray_coupled_tbath_offset=$w tx=$t  nproc=1 seed=" + str(random_seed))
        outfile.close()





def run_replicates(parentdir, cores):
    """For a list of directories, use cores number of cores to run the ER protocol"""
    replicate_dirs = [i for i in os.listdir(parentdir) if os.path.isdir(i) and 'replicate' in i]
    completedER = []
    
    
    x = len(replicate_dirs)/cores
    
    
    
    processedDirList = chunkIt(replicate_dirs, x)
    print (processedDirList)

    for i in processedDirList:
        procs = []
        for j in i:
            print(j)
            proc = multiprocessing.Process(target=runEnsembleRefinement, args=(j,))
            procs.append(proc)
            proc.start()
        for k in procs:
            k.join()









#################################################################################################################ANALYSIS MODULES##############################################################################################################


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
    """Takes an input xvg file (string) with rmsfs"""

    data = open(xvg_file)
    data = data.readlines()
    data = [i.strip() for i in data]
    data = [i for i in data if ('@' not in i) and ('#' not in i)]
    proc_data = []

    for i in data: 
        r = i.split()
        proc_data.append(r)
    
    chain_data = split_chains(proc_data)
 
    return chain_data




def processrmsf_2(inlist):

    data_series=[]
    for i in inlist:
        t_array = np.array(i)
        x_data=(list(t_array[:,0]))
        y_data=(list(t_array[:,1]))
        x_data=[int(i) for i in x_data]
        y_data=[float(i) for i in y_data]

        series = [x_data,y_data]

        data_series.append(series)
        
    return data_series


def processrmsf_3(inlist):
    chains_resids = []
    for i in inlist:
        chains_resids.append(i[0])

    average_resids = list(set.intersection(*map(set,chains_resids)))

   
    to_average_data = []
    for resid in average_resids:
        res_data_to_average = [resid]
        for chain in inlist:
            data_index = chain[0].index(resid)
            rmsf_dat = chain[1][data_index]
            res_data_to_average.append(rmsf_dat)
        to_average_data.append(res_data_to_average)

    return to_average_data            



def processrmsf_4(inlist):
    averaged_data = []
    for i in inlist:
        res_id = i[0]
        rmsf_data = i[1:]
        av_rmsf = np.mean(rmsf_data)
        av_list = [res_id,av_rmsf]
        #print (av_list)
        averaged_data.append(av_list)
    return averaged_data


def plot_xvg(xvg_file):
    """Takes an input xvg file (string), label for the x-axis (string) and for the y-axis (string);
    plots a plot"""

    data = open(xvg_file)
    data = data.readlines()
    data = [i.strip() for i in data]
    data = [i for i in data if ('@' not in i) and ('#' not in i)]
    x_axis=[]
    y_axis = []
    for i in data: 
        r = i.split()
        x_axis.append(float(r[0]))
        y_axis.append(float(r[1]))
        
    return (x_axis, y_axis)
    




def plot_replicate_RMSFs(parentdir):

    ER_directory = parentdir
    subplt_dims = (5,2)
    figure,axs = plt.subplots(subplt_dims[0], subplt_dims[1], figsize=(15,15))

    rows = [i for i in range(subplt_dims[0])]
    cols = [i for i in range(subplt_dims[1])]

    sbplt_size = len(rows)*len(cols)

    j=0
    k=0
    for i in range(1,sbplt_size+1):
        try: 
            row_num = rows[j]
            col_num = cols[k]
        except:
            j=0
            k+=1
            row_num = rows[j]
            col_num = cols[k]
        j+=1    

        try:
            rmsf_path = ER_directory+'/replicate_'+str(i)+'/rmsf.xvg'
            rmsf_data = plot_xvg(rmsf_path)
            axs[row_num,col_num].plot(rmsf_data[0], rmsf_data[1])
            axs[row_num,col_num].set_title('replicate_'+str(i))
        except:
            pass


        



def ensemble_RMSF(replicate_dir):
    return np.array(processrmsf_4(processrmsf_3(processrmsf_2(processrmsf_1(replicate_dir+'rmsf.xvg')))))




def optimisationAnalyser(parentdir):
    os.chdir(parentdir)
    data_list = []
    rep_dirs = [i for i in os.listdir() if os.path.isdir(i) and ('0' in i) ]
    for i in rep_dirs:
        os.chdir(str(i))
        logfile = glob.glob("*.log")
        logfile = open(str(logfile[0]), 'r')
        lines = logfile.readlines()
        for j in lines: 
            if ('FINAL' in j) and ('Rfree' in j):
                r_free = float(j[28:36])
        rep_name = i
        out_tuple = (i,r_free)
        data_list.append(out_tuple)
        os.chdir(parentdir)

    
    return data_list





def find_min_Rfree(parentdir):
    import glob
    os.chdir(parentdir)
    data_list = []
    rep_dirs = [i for i in os.listdir() if os.path.isdir(i) and ('replicate' in i) ]
    for i in rep_dirs:
        os.chdir(str(i))
        logfile = glob.glob("*.log")
        logfile = open(str(logfile[0]), 'r')
        lines = logfile.readlines()
        for j in lines: 
            if ('FINAL' in j) and ('Rfree' in j):
                r_free = float(j[28:36])
        rep_name = i
        out_tuple = (i,r_free)
        data_list.append(out_tuple)
        os.chdir(parentdir)
        r_frees = [data_list[i][1] for i in range(len(data_list))]
        repl_name = [data_list[i][0] for i in range(len(data_list))]
    print (repl_name)
    print (r_frees)
    return (data_list[r_frees.index(min(r_frees))][0], data_list[r_frees.index(min(r_frees))][1])







def get_rmsfs(parent_dir, ligRemove): 
    os.chdir(parent_dir)
    rep_dirs = [i for i in os.listdir() if os.path.isdir(i) and ('replicate' in i) ]

    remove_sel = [i +'|' if ligRemove.index(i)!=ligRemove.index(ligRemove[-1]) else i for i in ligRemove]
    remove_str = ''.join(remove_sel)
    
    for i in rep_dirs:
        os.chdir(i)
        subprocess.call('rm clean_ensemble.pdb', shell=True)
        print('Unzipping in directory '+str(i)+'...')
        subprocess.call('gunzip *.gz', shell=True)
        print('grepping: '+ 'grep -vwE "('+remove_str+')" *ensemble.pdb > clean_ensemble.pdb')
        subprocess.call('grep -vwE "('+remove_str+')" *ensemble.pdb > clean_ensemble.pdb', shell=True)
        subprocess.call('sed -i.bak \'s/KCX/ALA/g\' clean_ensemble.pdb > clean_ensemble2.pdb' , shell=True)
        print('Running GMX RMSF...')
        subprocess.call("echo '3'|gmx rmsf -f clean_ensemble.pdb -s clean_ensemble2.pdb -o rmsf.xvg -res", shell=True)
        print('Done')
        os.chdir(parent_dir)







def analyse(data_list):
    out_list = []
    for i in range(len(data_list[:,0])):
        ptls  = float(data_list[:,0][i].split('_')[0])
        wxray = float(data_list[:,0][i].split('_')[1])
        tx    = float(data_list[:,0][i].split('_')[2])
        rfree = float(data_list[:,1][i])
        out_list.append((ptls, wxray, tx, rfree))
    return np.array(out_list)


def plotOptimisation(parentdir, title):
    param_data = analyse(np.array(optimisationAnalyser(parentdir)))
    fig = plt.figure(figsize=(10,8))
    ax = fig.gca(projection='3d')
    p=ax.scatter(list(param_data[:,0]), list(param_data[:,1]), list(param_data[:,2]),c=param_data[:,3],       s=100, alpha=1,cmap=mtl.cm.hot,
                vmin=0.15, vmax=0.25)
    ax.set_xlabel(r'$p_{TLS}$', fontsize=20)
    ax.set_ylabel(r'$w_{xray}$', fontsize=20)
    ax.set_zlabel(r'$\tau_x$', fontsize=20)
    ax.set_title('Grid Search Results for '+ title, fontsize=20)
    cbar = fig.colorbar(p)
    cbar.ax.set_ylabel(r'$R_{free}$', fontsize=20)
