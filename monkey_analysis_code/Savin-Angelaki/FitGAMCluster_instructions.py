# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 15:29:07 2021

@author: ges6
"""
# google unix/linux to find commands

# login to cluster: ssh ges6@greene.hpc.nyu.edu
# cd /scratch/ges6
# ls

# cd GAM_fits
# ls GAM_fits
# vim PGAM_fit.sh

# cp: copy
# delete folder and contents: rm -r FolderName


# we need analysis code script, batch script, fit list, e.g. ls PGAM_fit.sh ; ls >> fit_with_cluster.py  PGAM_fit.sh  sim_list.npy

# copy .npz files to cluster (must not be connected to cluster):
scp C:\Data\firefly_analysis\LFP_band\concatenation_with_accel\m71*.npz ges6@greene.hpc.nyu.edu:/scratch/ges6/dataset_firefly/

# Create file in Python with instructions of what to fit:
dtype_dict = {'names':('session', 'unit_num') ,'formats':('U20',int)}    
table = np.zeros(0,dtype=dtype_dict)
tmp = np.zeros(1,dtype=dtype_dict)
tmp['session'] = 'm53s12'
table = np.hstack((table, tmp))

# copy batch list to cluster (generate it with create_batch_list.py)
scp "G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\preprocessing_pipeline\util_preproc\batch_list.npy"  ges6@greene.hpc.nyu.edu:/scratch/ges6/GAM_fits/
#the batch list array is called within the fitting function fit_with_cluster.py

# copy fitting script to cluster
scp "G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\fitting\fit_with_cluster.py"  ges6@greene.hpc.nyu.edu:/scratch/ges6/GAM_fits/

# run fit with cluster
sbatch --array=1-458  PGAM_fit.sh # batchscript: specs about cluster job
# where 1-400 is the length of the table I'm fitting


# status of job: squeue -u ges6
# edit batch script: vim PGAM_fit.sh
# escape/insertion
# quit:hitEsc,colon (shift+colon), q, enter (to discard changes do q! instead)


# to make a test run before submitting the job to the cluster:
module purge
#. /etc/profile.d/modules.sh
module load anaconda3/2020.02
module load r/intel/4.0.3
#source activate /home/jpn5/.local/lib/pycuda3.6
source /scratch/eb162/venv/bin/activate

srun --mem=20GB --pty --cpus-per-task=1 python -m pdb fit_with_cluster.py 1
# then hit 'c' to continue because it is in debug mode (pdb)




# get fits from cluster (must not be connected to cluster)
scp ges6@greene.hpc.nyu.edu:/scratch/ges6/GAM_fits/*_m71*/* "G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\fits_from_cluster"
scp ges6@greene.hpc.nyu.edu:/scratch/ges6/GAM_fits/gam_m71*/* "G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\fits_from_cluster"
scp ges6@greene.hpc.nyu.edu:/scratch/ges6/GAM_fits/*_m71s*_fit_info*tte_mag* "G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\fits_from_cluster"

# to list all files with a certain pattern in naming
grep m73s26

# to check the time a file was created
ls -lrt | grep m73s26
 




# send matlab extraction code
scp $(find "G:\My Drive\MATLAB\Code\Savin-Angelaki\" -type f ! -name "PGAM") ges6@greene.hpc.nyu.edu:/scratch/ges6/extract_data/extraction_code/
scp -r "G:\My Drive\MATLAB\Code\Savin-Angelaki\!(*PGAM*)" ges6@greene.hpc.nyu.edu:/scratch/ges6/extract_data/extraction_code/
rsync -av -e ssh --exclude='PGAM' "G:\My Drive\MATLAB\Code\Savin-Angelaki" ges6@greene.hpc.nyu.edu:/scratch/ges6/extract_data/extraction_code/

scp -r "G:\My Drive\MATLAB\Code\Savin-Angelaki\firefly-monkey" ges6@greene.hpc.nyu.edu:/scratch/ges6/extract_data/extraction_code/
scp -r "G:\My Drive\MATLAB\Code\Savin-Angelaki\nonparametric-regression" ges6@greene.hpc.nyu.edu:/scratch/ges6/extract_data/extraction_code/
scp -r "G:\My Drive\MATLAB\Code\Savin-Angelaki\neuroGAM" ges6@greene.hpc.nyu.edu:/scratch/ges6/extract_data/extraction_code/
scp -r "G:\My Drive\MATLAB\Code\Savin-Angelaki\matlab-util-tools" ges6@greene.hpc.nyu.edu:/scratch/ges6/extract_data/extraction_code/
scp -r "G:\My Drive\MATLAB\Code\Savin-Angelaki\OfflineSDK" ges6@greene.hpc.nyu.edu:/scratch/ges6/extract_data/extraction_code/

scp -r "Z:\Data\Monkey2_newzdrive\Viktor\Utah Array\Dec*2021"  ges6@greene.hpc.nyu.edu:"/scratch/ges6/extract_data/raw_data/Viktor/Utah Array"

scp -r "G:\My Drive\MATLAB\Code\Savin-Angelaki\firefly-monkey\Analysis\monkeyInfoFile_joysticktask_cluster.m" ges6@greene.hpc.nyu.edu:/scratch/ges6/extract_data/extraction_code/firefly-monkey\Analysis\


# vim commands
# :q! : exit without saving
# :x : exit and save 
# :w : save
