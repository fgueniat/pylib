import numpy as np
import matplotlib.pyplot as plt
import plot_tools as pt
from cantera_tools import chemreac
import obs_tools as ot
import sys
import div_tools as dt
norm = np.linalg.norm

def params(argv=None,inputs = None):
 	if inputs is None: 
			inputs = {
			'input_spec':'H2',
			'input_n':4,
			'input_n_square':0,
			'input_n_cubic':0,
			'input_tp':(1200,101325),
			'input_ci':np.array([1.e-03,1.0e-09,1.e-03,1.e-09,1.e-9,1.e-09,1.e-9,1.e-9]),
			'input_nt':150,
			'input_tmin':None,
			'input_tmax':None,
			'input_rho':None,
			'input_eps':None,
			'input_compute_random':False,
			'input_mechanism':'boivin.cti',
			'input_compute_sm':True,
                        'input_method':'Nelder-Mead',
			'input_ci_random':None,
			'input_tp_random':None,
			'input_minimization':None,
			'input_verbose':False,
                        'input_index_sorted':None,
			'input_plots':[False,False,False,False,False,False],
			'ploc':{},
			'outputs':{},
			}
			inputs['ploc']['debug']={}
			#plots: concentration, mass fraction,rates,
			#ploc: local parameters
	try:
		inputs['input_spec' ] = argv[1]
		inputs['input_n'] = int(argv[2])
		inputs['input_n_square'] = int(argv[3])
		inputs['input_n_cubic'] = int(argv[4])
		inputs['input_ci'] = dt.split_list(argv[5])
		inputs['input_tp'] = dt.split_tuple(argv[6])
		inputs['input_ci_random'] = dt.split_list(argv[7])
		inputs['input_tp_random'] = dt.split_tuple(argv[8])
	except:
		pass
	return inputs

def extended_sm(inputs,specs=None,verbose=False):
	#compute algebraic constraints for all species
	if specs is None:
		r = reaction_init(inputs)
		specs = [r.list[i] for i in xrange(r.n_species)] 
	n_s = len(specs)
	dict_param = {}
	# verif param are 0:
	ploc = inputs['ploc']
	verbose = inputs['input_verbose']
	ploc={}
	inputs['input_compute_sm'] = True
	inputs['input_compute_random'] = False
	inputs['input_minimization'] = None
	#
	order_observability(inputs)
	#
	for ispec,spec in enumerate(specs):
		if verbose is True:print 'spec # ' + str(ispec) + ' over ' + str(n_s) 
		inputs['input_spec']=spec
		id_algebraic_constrains(inputs)
		dict_param[spec] = inputs['outputs'][spec].copy()
	return dict_param

def observable(z,alpha,inputs):
	ploc = inputs[ 'ploc']
	y =  z[ploc['ind2rec']]
	y += np.sum([alpha[i]*z[ploc['ind2use'][i]] for i in xrange(inputs['input_n'])])
	y += np.sum([alpha[i+inputs['input_n']]*(z[ploc['ind2use'][i]]**2) for i in xrange(inputs['input_n_square'])])
	y += np.sum([alpha[i+inputs['input_n']+inputs['input_n_square']]*(z[ploc['ind2use'][i]]**3) for i in xrange(inputs['input_n_cubic'])])
	return y

def specs2inds(inputs,spec2use=None):

	ploc = inputs['ploc']
	spec2rec = inputs['input_spec']
	ind2rec = inputs['outputs']['main_reaction'].list[spec2rec]	
	ploc['ind2rec'] = ind2rec
	ploc['spec2rec'] = spec2rec
	# more efficient to use directly the species in the order of observability
	obs_prop = inputs['outputs']['observable_sorted']
	if spec2use is None:
	#spec to reconstruct not being present
		spec2use = [spec for spec in obs_prop[0:inputs['input_n']+1] if spec not in [spec2rec]]
		if len(spec2use)>inputs['input_n']:spec2use.pop() 		
	#
	ind2use = np.array([inputs['outputs']['main_reaction'].list[spec] for spec in spec2use if spec not in [spec2rec]])
	ploc['ind2use']=ind2use
	spec2use = [inputs['outputs']['main_reaction'].list[ind] for ind in ind2use]
	ploc['spec2use']=spec2use

	
def reaction_init(inputs,israndom=False):
	ploc=inputs['ploc']
	verbose = inputs['input_verbose']
	#
	#mechanism = 'boivin.cti'# mechanism
	mechanism = inputs['input_mechanism']# mechanism
	prob = ('0D','pressure') #type of problem
	#
#	n_s = len(inputs['input_ci']) # number of species
	########################################
	# initial conditions
	V = 1.e3 * (1.e-2)**3 # initial volume
	if israndom is False: 
		c_init = inputs['input_ci']	
		#ci_random = inputs['input_ci_random']
		tpinit = inputs['input_tp']
	else:
		c_init = inputs['input_ci_random']	
		#ci_random = inputs['input_ci_random']
		tpinit = inputs['input_tp_random']
	#######################################
	#define the time for integration
	nt = inputs['input_nt']
	inputs['ploc']['nt'] = nt
	if inputs['input_tmin'] is None:tmin = 1.e-7
	else:tmin =inputs['input_tmin']
	if inputs['input_tmax'] is None:tmax = 1.e-2
	else:tmax =inputs['input_tmax']
	
	inputs['ploc']['tmin'] = tmin
	inputs['ploc']['tmax'] = tmax
	offset_ts = .4
	ploc['offset_ts'] = offset_ts
	timelog = np.logspace(np.log10(tmin),np.log10(tmax),nt)
	ploc['time'] = timelog

	#######################################
	# define the reaction
	reaction = chemreac(cti = mechanism, problem = prob,V=V)
	reaction.TP(tp=tpinit,reset = False)
	reaction.concentrations(c=c_init,reset = True)
	# integrate the system
	reaction.integrate(time = timelog)
	return reaction

def scale_(zh,z,offset='end'):
	if offset == 'begin':
		off_zh = zh[0]
		off_z  = z[0]
	if offset == 'end':
            off_zh = np.mean(zh[-5:])
            off_z  = np.mean(z[-5:])
	zh = zh - off_zh
	zh=zh+off_z
 	return zh

def reconstruction(a,inputs,offset='end',r=None):
	time = inputs['ploc']['time']
	
	ind2rec = inputs['ploc']['ind2rec']
	if r is None: 
		try: r=inputs['ploc']['reaction']
		except:r=reaction_init(inputs)		
	# observable
	obs = np.array([observable(r.h_z[i+1],a,inputs) for i in xrange(len(time) )])
	# reconstruction
	z1 = r.time_series('z',ind2rec)[1:]
	# obs is (supposedly) approx a constant.
	# instead of using the constrains, we use obs for the reconstruction
	z1hat = -obs + z1
	z1hat_offset = scale_(z1hat,z1)	
	z1_offset=z1
	###########################
	return obs,z1hat,z1,z1hat_offset

def order_observability(inputs):
	ploc = inputs['ploc']
	##################################################
	# This section sets the parameters of the mechanism
	reaction = reaction_init(inputs)
	inputs['outputs']['main_reaction']=reaction
	n_s = reaction.n_species
	################################################################
	# This section sets the parameters of the observability analysis
	#eps_,rho=5.e-2,1.e-9
	if inputs['input_rho'] is None:rho = 1.e-9
	else:rho = inputs['input_rho']
	if inputs['input_eps'] is None:eps_ = 5.e-2
	else:eps_ = inputs['input_eps']

	ploc['eps_']=eps_
	ploc['rho']=rho
	#eps_ = 0.01*c_init
	#
	mm = 'trace' # measure
	ploc['measure'] = mm
	##############################
	# exploration:
	c_init = inputs['input_ci']
	time = ploc['time']
	exp_dx,eps_eff = ot.f_explore(reaction,c_init,time,eps_ = eps_)
	ploc['exp_dx']=exp_dx
	ploc['eps_eff']=eps_eff
	# order of observability
	fobs = lambda x: np.diag(x)
	w=ot.EG(fobs,exp_dx,eps_eff,time)
	w0=np.mean(w,axis=0)
	obs_prop = np.argsort(w0)[::-1]
	obs_sorted = [reaction.list[obs_prop[i]] for i in xrange(n_s)]
	if inputs['input_verbose'] is True:
		print 'from most observable to least observable:'
		print ', '. join(obs_sorted)
	#ploc['observable_sorted'] = obs_sorted
	inputs['outputs']['observable_sorted'] = obs_sorted

def id_algebraic_constrains(inputs,reaction=None,spec2use=None):
	ploc = inputs['ploc']
	outputs=inputs['outputs']
	#if 'reaction' not in ploc:order_observability(inputs)
	if reaction is None:reaction = outputs['main_reaction']
	measure = ploc['measure']
	eps_eff,rho = ploc['eps_eff'],ploc['rho']
	exp_dx,offset_ts=ploc['exp_dx'],ploc['offset_ts']
	time = ploc['time']
	#if input_n is not None: n_obs,n_obs_square,n_obs_cubic = (input_n,input_n_square,input_n_cubic)
	n_obs,n_obs_square,n_obs_cubic = inputs['input_n'],inputs['input_n_square'],inputs['input_n_cubic']

	specs2inds(inputs)
	
	spec2rec = ploc['spec2rec']
	ind2rec = ploc['ind2rec']
	ind2use = ploc['ind2use']
	spec2use = ploc['spec2use']
        method = inputs['input_method']
	#########################################
	# This section computes the slow manifold
	compute_sm = inputs['input_compute_sm']
	if compute_sm is True:
		cinit_min = np.zeros(n_obs+n_obs_square+n_obs_cubic)
		if inputs['input_minimization'] is None:
			fobs = lambda x,a: observable(x,a,inputs)
			ploc['debug']['fobs']=fobs
			a = ot.slow_manifold(cinit_min,dx = exp_dx, f_def_obs = fobs, eps_ = eps_eff,time = time,rho=rho,method = 'BFGS',measure = measure,verbose = False,offset = offset_ts,maxiter=[1000,100],n_recuit = 10, isrecuit=False,r=reaction)
		if inputs['input_verbose'] is True:
			print a
		gram=ot.EG(fobs = lambda x:observable(x,a.x,inputs),dx=exp_dx,eps_=eps_eff,time=time,offset=offset_ts)
	#
	else:
		a = inputs['input_minimization']
	inputs['outputs'][spec2rec] = {'output_minimization':a,'spec2use':spec2use,}
	

def plot_mech(inputs,reaction=None,xlim=None):
	ploc = inputs['ploc']
	outputs = inputs['outputs']
	#
	if reaction is None:
		if 'main_reaction' not in outputs.keys():order_observability(inputs)
		reaction = outputs['main_reaction']
	#time
	tmin,tmax = ploc['tmin'],ploc['tmax']
	n_s = reaction.n_species
	tmini = int(np.floor(np.log10(tmin)))
	tmaxi = int(np.ceil(np.log10(tmax)))
	#main parameters
	pd = pt.Paradraw()
	if n_s<10:pd.legend = [reaction.list[i] for i in xrange(n_s)]
	pd.x_scale = 'log'
	if xlim is None:
		pd.xlim = [tmin,tmax]
		pd.x_tick_label = [1.*10**i for i in xrange(tmini+1,tmaxi+1,2)]
	else: pd.xlim = xlim
	pd.x_label = 'time'
	import matplotlib.colors as colors
	pd.colors = []
	for ic,c in enumerate(colors.cnames):
		if ic<n_s:pd.colors.append(str(c))		
#	pd.colors = ['tan','g','r','b','k','y','magenta','aqua']
		
	cc = [reaction.time_series('concentrations',i) for i in xrange(n_s)]
	tt = (reaction.time,)*n_s
	pd.y_label = 'n'
	pd.ylim = [np.min(cc),np.max(cc)]
	pt.multiplot2(tt,cc,pd)
#
	zz = [reaction.time_series('z',i) for i in xrange(n_s)]
	pd.y_scale = 'log'
	pd.y_label = 'z'
	pd.ylim = [np.min(zz),np.max(zz)]
	pt.multiplot2(tt,zz,pd)
#
	pd.y_scale = 'symlog'
	pd.y_label = 'p rate'
	pp = [reaction.time_series('rates',i) for i in xrange(n_s)]
	pd.ylim = [np.min(pp),np.max(pp)]
	pt.multiplot2(tt,pp,pd)


def rec(inputs=params()):
	
	ploc = inputs['ploc']
	specs2inds(inputs)
	###########################################################
	# This section just plots the reconstruction of one species
	offs='end'
	if inputs['input_minimization'] is None:
		if inputs['input_spec'] in inputs['outputs'].keys():
			a = inputs['outputs'][inputs['input_spec']]['output_minimization']
		else:
			a = id_algebraic_constrains(inputs)
	else:
		a = inputs['input_minimization']
	
	try:
		reaction = ploc['reaction']
	except:
		reaction = reaction_init(inputs)
	ret = reconstruction(a.x,inputs,offset=offs,r=reaction)

	###########################################################
	# This section just plots the reconstruction of one species
	#, with (strongly) different initial conditions
	rand_test,rand_plot = inputs['input_compute_random'],inputs['input_plots'][4]
	rand_test = False
	if rand_test is True:
		nr=1
		if inputs['input_tp_random'] is None:tpinit_random = tpinit
		else:tpinit_random = inputs['input_tp_random']
		reaction2 = reaction_init(inputs,israndom=True)
		ret2 = reconstruction(a.x,inputs,r=reaction2)
		if inputs['input_verbose'] is True:
			print 'spec mole: ' + str(reaction2.h_z[0])
	else:
		ret2 = None
		reaction2=None

#	gram =lambda a: ot.EG(fobs = lambda x:observable(x,a),dx=exp_dx,eps_=eps_eff,time=timelog,offset=offset_ts)
#	rr = lambda a: reconstruction(a,r=reaction,isplot=False)

#
#
#	ret_ = {
#		'reaction':reaction,'retour':ret,
#		'reaction_pert':reaction2,'retour_pert':ret2,
#		'result_min':a,
#		'indexes':ind2use,
#		'time':timelog,
#		'obs_index':obs_sorted,
#		'gramian':w,
#		'parameters':inputs
#		}

	ret_ = {
		'reaction':reaction,'retour':ret,
		'reaction_pert':reaction2,'retour_pert':ret2,
		'result_min':a,
		}


	return ret_
	#return reaction,ret,reaction2,ret2,a,ind2use[0:n_obs],timelog,obs_sorted

def rec_Temp(inputs):
    if 'main_reaction' not in inputs['outputs'].keys():
        order_observability(inputs)
    reaction = inputs['outputs']['main_reaction']
    specs = inputs['outputs']['observable_sorted']
    d={}
   # inputs['ploc'].pop('reaction',None)
    for spec in specs:
        inputs['input_spec'] = spec
#        inputs['input_minimization'] = inputs['outputs'][spec]['output_minimization']
        d[spec] = rec(inputs)
    z = [d[spec]['retour'][2] for spec in specs]
    reaction = d['H']['reaction']
    z_hat = [d[spec]['retour'][3] for spec in specs]
    time = inputs['ploc']['time']
    #reaction = inputs['ploc']['reaction']
    t_real = reaction.h_T[1:]
    w = reaction.gas.molecular_weights
    t_hat = []
    v_real = reaction.h_volume
    #v_pert = reaction_pert.h_volume
    for it,tt in enumerate(time):
        zh = [z_hat[i][it] for i in xrange(reaction.n_species)]
        yh = zh*w
        t_hat.append(reaction.y2t(y=yh,v=v_real[it+1]))
    t = reaction.h_T[1:]
    return t,t_hat,time

if __name__ == "__main__":
	inputs = params(sys.argv)
	main(inputs=inputs)
