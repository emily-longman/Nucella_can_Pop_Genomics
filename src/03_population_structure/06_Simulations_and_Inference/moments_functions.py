def get_mig_mat(n_pops: int, m: float) -> np.ndarray:
    """
    :param n_pops int: Number of populations
    :param m float: Migration rate, in coalescent units of `1 / 2N_anc`
    :return: Square `n_pops` x `n_pops` matrix representing symmetric 
        migration, equals `m` everywhere except for zeros along the diagonal
    :rtype: np.ndarray
    """
    return m * np.ones([n_pops] * 2) - np.diag([m] * n_pops)

def twosplits(params, ns, pop_ids=None):
    nu1, nu_intermediate, nu2, nu3, T1, T2, m2, m3 = params
    mig_mat2 = get_mig_mat(2, m2)
    mig_mat3 = get_mig_mat(3, m3) 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = fs.split(0, ns[0], ns[1] + ns[2])
    fs.integrate([nu1, nu_intermediate], T1, m=mig_mat2)
    fs = fs.split(1, ns[1], ns[2])
    fs.integrate([nu1, nu2, nu3], T2, m=mig_mat3)
    fs.pop_ids = pop_ids
    return fs

    
#Just make sure T_split2 < T_split1 when you evaluate parameters.  
lower_bounds = [0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001]
upper_bounds = [10, 10, 10, 10, 1000, 1000, 0.5, 0.5]
params_guess = [0.5, 0.5, 1.0, 1.0, 100, 50, 0.05, 0.05]

popt=moments.Inference.optimize_log(params_guess, fs_folded, 
twosplits,
lower_bound=lower_bounds, 
upper_bound=upper_bounds,
verbose=True, 
maxiter=100)


#####
def admixture(params, ns, pop_ids=None):
    nu1, nu2, nu_admix, T_split, T_admix, m2, m3, admix_prop = params
    mig_mat2 = get_mig_mat(2, m2)
    mig_mat3 = get_mig_mat(3, m3)
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + 2 * ns[2])
    fs = moments.Spectrum(sts)
    fs = fs.split(0, ns[0] + ns[2], ns[1] + ns[2])
    fs.integrate([nu1, nu2], T_split, m=mig_mat2)
    fs = fs.admix(0, 1, ns[2], admix_prop)
    fs.integrate([nu1, nu2, nu_admix], T_admix, m=mig_mat3)
    fs.pop_ids = pop_ids
    return fs.fold()

params_guess = [1.0, 1.0, 1.0, 100, 50, 0.01, 0.01, 0.2]
lower_bounds = [0.001, 0.001, 0.001, 10, 10, 0.001, 0.001, 0.001]
upper_bounds = [100, 100, 100, 1000, 1000, 0.1, 0.1, 0.99]

popt=moments.Inference.optimize_log(params_guess, fs_folded, 
admixture,
lower_bound=lower_bounds, 
upper_bound=upper_bounds,
verbose=True, 
maxiter=100) 

### serial splits
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

params_guess = [1.0, 1.0, 1.0, 1.0, 100, 50, 0.01, 0.01]
lower_bounds = [0.001, 0.001, 0.001,  0.001, 10, 10, 0.001, 0.001]
upper_bounds = [100, 100, 100, 100, 1000, 1000, 0.1, 0.1]


params = len(["nu1", "nu_intermediate", "nu2", "nu3", "T1", "T2", "m12", "m23"])
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
nu2=popt[1]
nu3=popt[2] 
T1=popt[4] 
T2=popt[5]  
