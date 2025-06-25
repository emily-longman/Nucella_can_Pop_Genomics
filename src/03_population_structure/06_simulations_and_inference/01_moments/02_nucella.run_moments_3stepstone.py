#module load python3.12-anaconda/2024.06-1;conda activate moments_jcbn; python

# import packages that'll be used
import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
import dadi
from dadi import Misc
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from datetime import datetime

###
i=int(sys.argv[1])
print("now processing", i)
###

indat = pd.read_csv("trios_guide.txt",header=0,delimiter="\t")
Ancestral_id1=indat.Parent1[i]
Ancestral_id2=indat.Parent2[i]
Derived_id=indat.Derived[i]
root_dadi="dadi_objects_trios"
namsam = ["probs",Ancestral_id1, Ancestral_id2, Derived_id,"delim"]
paths = [root_dadi,".".join(namsam)]
"/".join(paths)
fs_file = "/".join(paths)
###
Pair_name = Derived_id
###
### collect the metadata
root_meta="./L_meta_objects_trios"
namsam2 = ["probs",Ancestral_id1, Ancestral_id2, Derived_id,"meta"]
paths2 = [root_meta,".".join(namsam2)]
"/".join(paths2)
meta_file = "/".join(paths2)
inmet = pd.read_csv(meta_file,header=0,delimiter="\t")
Ancestral_n1=inmet.an1_ne[0]
Ancestral_n2=inmet.an2_ne[0]
Derived_n=inmet.Der_ne[0]
L = inmet.L[0]

print("fs file is", fs_file )
print("Pair_name is", Pair_name )
print("names are", Ancestral_id1,Ancestral_id2,Derived_id )

#constants
mu = 2.8e-9 #from Keightley et al. 2014
#g = 0.06666667 #equals 15 gen/year --- See Pool 2015
iterations=10
#####

####
print('loading data')
dd = dadi.Misc.make_data_dict(fs_file) #reads in genomalicious SNP file
data = pd.read_csv(fs_file, sep="\t", nrows=1)

#####
#setting pop id's and projections from if/else
pop_id=["A_parent1",  "B_parent2",  "C_derived"]
pools=[Ancestral_n1, Ancestral_n2, Derived_n ] 

print('folding sfs')
fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=pools, polarized=False) #takes data dict and folds
ns = fs_folded.sample_sizes #gets sample sizes of dataset
S = fs_folded.S()

####
#opening output file to give column names
PMmod=open('o3step/%s_output.3step.txt' % Pair_name,'w')
PMmod.write(
            str("Pair_name")+'\t'+ #print pair name
            str("L")+'\t'+ #double checking L is working as I want it to
            str("theta")+'\t'+ 
			str("nu1")+'\t'+ 
			str("nu2")+'\t'+ 
			str("nu3")+'\t'+ 
			str("T1")+'\t'+ 
			str("T2")+'\t'+ 
			str("m12")+'\t'+ 
			str("m23")+'\t'+ 
            str("fs_sanitycheck")+'\t'+
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()

print('defining functions')
####

def twosplits_step(params, ns, pop_ids=None):
    nu1, nu_intermediate, nu2, nu3, T1, T2, m12, m23 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = fs.split(0, ns[0], ns[1] + ns[2])
    # Migration matrix: only neighbors connected
    mig_mat2 = [
        [0,   m12],
        [m12, 0  ]
    ]
    fs.integrate([nu1, nu_intermediate], T1, m=mig_mat2)
    fs = fs.split(1, ns[1], ns[2])
     # Migration matrix: only neighbors connected
    mig_mat3 = [
        [0,   m12, 0],
        [m12, 0,   m23],
        [0,   m23, 0]
    ]
    fs.integrate([nu1, nu2, nu3], T2, m=mig_mat3)
    fs.pop_ids = pop_ids
    return fs

func_moments = twosplits_step
###   idx0, idx1, num_lineages, proportion, new_id=None
#params_guess = [1.0, 1.0, 1.0, 1.0, 100, 50, 0.01, 0.01]
lower_bounds = [0.001, 0.001, 0.001,  0.001, 10, 10, 0.001, 0.001]
upper_bounds = [100, 100, 100, 100, 1000, 1000, 0.1, 0.1]

#upper_bound = [100, 100, 100, 100, 100, 1]
#lower_bound = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 0]


print('optimization loop')
for i in range(int(iterations)): #iterations is imported from sys. argument #1
	print("starting optimization "+str(i))
#This number is 5. i.e., count parameters. 
	#Start the run by picking random parameters from a uniform distribution.
	#popt=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(params)]
	params = len(["nu1", "nu_intermediate", "nu2", "nu3", "T1", "T2", "m12", "m23"])    
	#This is the optimization step for moments.
	#popt is the prior.
	#fs folded is a tranform SFS by folding it. The original SFS loaded in sys.arg #1 is quasi-folded. i.e. polarized to reference genome.
	#Folding it is a must, because reference is not 100% ancestral
	params_guess=[np.random.uniform(lower_bounds[x],upper_bounds[x]) for x in range(params)]
	popt=moments.Inference.optimize_log(params_guess, fs_folded, 
	func_moments,
	lower_bound=lower_bounds, 
	upper_bound=upper_bounds,
	verbose=True, 
	maxiter=100)
	
	model = func_moments(popt, ns)
	
	#Calculate log likelihood of the model fit
	ll_model=moments.Inference.ll_multinom(model, fs_folded)
	#Now calculate AIC of model fit
	aic = 2*params - 2*ll_model
	print('Maximum log composite likelihood: {0}'.format(ll_model))	
	#Now estimate theta from model fit
	theta = moments.Inference.optimal_sfs_scaling(model, fs_folded)
	#reducing complexity of calculations to follow by adding variables in lieu of expressions/esoteric df calls
	Nref= theta/(4*mu*L)
	nu1=popt[0]
	nu2=popt[2]
	nu3=popt[3] 
	T1=popt[4] 
	T2=popt[5]
	m12=popt[6]
	m23=popt[7]
	
	#Open the output file
	PMmod=open('o3step/%s_output.3step.txt' % Pair_name,'a')
	    #Dumping output ot outfile
	PMmod.write(
            str(Pair_name)+'\t'+ #print pair name
            str(L)+'\t'+ #double checking L is working as I want it to
            str(theta)+'\t'+ 
			str(nu1)+'\t'+ 
			str(nu2)+'\t'+ 
			str(nu3)+'\t'+ 
			str(T1)+'\t'+ 
			str(T2)+'\t'+ 
			str(m12)+'\t'+ 
			str(m23)+'\t'+ 
            str(fs_file)+'\t'+
            str(ll_model)+'\t'+
            str(aic)+'\n')
	PMmod.close()
	
	print("Moments finished running")
