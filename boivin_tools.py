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
			'input_compute_random':False,
			'input_compute_sm':True,
			'input_ci_random':None,
			'input_tp_random':None,
			'input_minimization':None,
			'input_verbose':False,
			'input_plots':[False,False,False,False,False,False],
			}
			#plots: concentration, mass fraction,rates,
	try:
		inputs['input_spec'] = argv[1]
		inputs['input_n'] = int(argv[2])
		inputs['input_n_square'] = int(argv[3])
		inputs['input_n_cubic'] = int(argv[4])
		inputs['input_ci'] = dt.split_list(argv[5])
		inputs['input_tp'] = dt.split_tuple(argv[6])
		inputs['input_ci_random'] = dt.split_list(argv[7])
		inputs['input_tp_random'] = dt.split_tuple(argv[8])
	except:
		pass
	inputs['spec'] = inputs['input_spec']
	return inputs


def main(inputs=params()):
#	
	##################################################
	# This section sets the parameters of the mechanism
	if True:
		verbose = inputs['input_verbose']

		mechanism = 'boivin.cti'# mechanism
		prob = ('0D','pressure') #type of problem
		#
		n_s = 8 # number of species
		########################################
		# initial conditions
		V = 1.e3 * (1.e-2)**3 # initial volume
		#n_h2,n_o2,n_h2o,n_h,n_o,n_oh,n_ho2,n_h2o2,n_n = 1.e-3,1.e-3,1.e-9,1.e-9,1.e-9,1.e-9,1.e-9,1.e-9,1.e-5	
		#c_init = np.array([n_h2,n_h,n_o2,n_o,n_oh,n_ho2,n_h2o,n_h2o2 ])
		c_init = inputs['input_ci']		
		ispert = False
		if ispert is True:pert = 0.10*(np.random.rand(n_s)-0.5)*2.
		else:pert = 0.
		c_init = c_init*(np.ones(n_s)+pert)
		#ci_random =  np.array([ 9.6e-04, 1.0e-09, 6.5e-04, 1.3e-09, 1.1e-9, 1.2e-09, 4.2e-10,   2.3e-10])
		ci_random = inputs['input_ci_random']
		tpinit = inputs['input_tp']
		#######################################
		#define the time for integration
		nt = 150
		tmin = 1.e-7
		tmax = 1.e-2
		offset_ts = .4
		timelog = np.logspace(np.log10(tmin),np.log10(tmax),nt)
		#######################################
		# define the reaction
		reaction = chemreac(cti = mechanism, problem = prob,V=V)
		reaction.TP(tp=tpinit,reset = False)
		reaction.concentrations(c=c_init,reset = True)
		# integrate the system
		reaction.integrate(time = timelog)

	def param_plot():# params for plots
		#time
		tmini = int(np.floor(np.log10(tmin)))
		tmaxi = int(np.ceil(np.log10(tmax)))
		#main parameters
		pd = pt.Paradraw()
		pd.legend = [reaction.list[i] for i in xrange(n_s)]
		pd.x_scale = 'log'
		pd.xlim = [tmin,tmax]
		pd.x_tick_label = [1.*10**i for i in xrange(tmini+1,tmaxi+1,2)]
		pd.x_label = 'time'
		pd.colors = ['tan','g','r','b','k','y','magenta','aqua']
		return pd

	def reconstruction(a,title=None,std=False,offset='end',r=None,isobs=True,isplot=False):#plot reconstruction
			if r  is None: return -1
			def scale_(zh,z):
				if offset == 'begin':
					off_zh = zh[0]
					off_z  = z[0]
				if offset == 'end':
					off_zh = zh[-1]
					off_z  = z[-1]
				zh = zh - off_zh
				zh=zh+off_z
				return zh
			# observable
			obs = np.array([observable(r.h_z[i+1],a) for i in xrange(len(timelog) )  ])
			p = param_plot()
			p.y_label = 'observable'
			p.colors = ['k']
			p.thickness = [2]
			p.legend=[False]
			if isobs is True and isplot is True:pt.plot2((timelog,obs),p)
			# reconstruction
			z1 = r.time_series('z',ind2rec)[1:]
			# obs is (supposedly) approx a constant.
			# instead of using the constrains, we use obs for the reconstruction
			z1hat = -obs + z1
			z1hat_offset = scale_(z1hat,z1)	
			z1_offset=z1
			###########################
			# plot
			p=param_plot()
			p.title = title
			p.marks = ['','-']
			p.markers = ['o','']
			p.colors = ['r','k']
			p.thickness = [2,2]
			p.legend=[False]
			p.y_label = 'z'
			p.y_scale = 'symlog'
			p.legend=['$z_{'+spec2rec + '}$','$\hat{z}_{'+ spec2rec + '}$']
			tt = (timelog,)*2
			if isplot is True:pt.multiplot2(tt,(z1_offset,z1hat_offset),p)
			#
			return obs,z1hat,z1,z1hat_offset
	#
	################################################################
	# This section sets the parameters of the observability analysis
	if True:
		eps_,rho=5.e-2,1.e-9
	#	eps_ = 0.01*c_init
		eps_eff = eps_

		mm = 'trace' # measure
		##############################
		# exploration:
		exp_dx,eps_eff = ot.f_explore(reaction,c_init,timelog,eps_ = eps_)
		# order of observability
		fobs = lambda x: np.diag(x)
		w=ot.EG(fobs,exp_dx,eps_eff,timelog)
		w0=np.mean(w,axis=0)
		obs_prop = np.argsort(w0)[::-1]
		obs_sorted = [reaction.list[obs_prop[i]] for i in xrange(n_s)]
		if verbose is True:
			print 'from most observable to least observable:'
			print ', '.join(obs_sorted)
	#	
	##################################################
	# This section sets the parameters of the manifold
	if True:
		spec2use = None 
		
	#	if input_n is not None: n_obs,n_obs_square,n_obs_cubic = (input_n,input_n_square,input_n_cubic)
		n_obs,n_obs_square,n_obs_cubic = inputs['input_n'],inputs['input_n_square'],inputs['input_n_cubic']
		spec2rec = inputs['input_spec']
		ind2rec = reaction.list[spec2rec]
		# more efficient to use directly the species in the order of observability
		if spec2use is None:
		#spec to reconstruct not being present
			ind2use = np.array([ind for ind in obs_prop if ind not in [ind2rec]])
		else:
		# or specify:
			ind2use = np.array([reaction.list[spec] for spec in spec2use if spec not in [spec2rec]])
		spec2use = [reaction.list[ind] for ind in ind2use]	
		# Observable function
		def observable(z,alpha):
			y =  z[ind2rec]
			y += np.sum([alpha[i]*z[ind2use[i]] for i in xrange(n_obs)])
			y += np.sum([alpha[i+n_obs]*(z[ind2use[i]]**2) for i in xrange(n_obs_square)])
			y += np.sum([alpha[i+n_obs+n_obs_square]*(z[ind2use[i]]**3) for i in xrange(n_obs_cubic)])
			return y
		def dobservable(z,alpha):
			#dy = np.zeros(n_obs+n_obs_square+n_obs_cubic)
			#for i in xrange(n_obs+n_obs_square+n_obs_cubic):
			dy = [z[ind2use[i]] for i in xrange(n_obs)]
			dy += [z[ind2use[i]]**2 for i in xrange(n_obs_square)]
			dy += [z[ind2use[i]]**3 for i in xrange(n_obs_cubic)]
			return np.array(dy)
	#	 
		def J(alpha): 
			''' cost function '''
			return ot.cost_gramian(alpha,f_def_obs=observable,rho=rho,dx=exp_dx,eps_=eps_eff,time=timelog, verbose = False, measure = mm)
	#
	########################################
	# This section just plots the  mechanism
	if False:
		plot_concentration = inputs['input_plots'][0]
		if plot_concentration is True:
			cc = [reaction.time_series('concentrations',i) for i in xrange(reaction.n_species)]
			tt = (reaction.time,)*n_s
			pd =param_plot()
			pd.y_label = 'n'
			pt.multiplot2(tt,cc,pd)
	#
		plot_mass_fraction = inputs['input_plots'][1]
		if plot_mass_fraction is True:
			zz = [reaction.time_series('z',i) for i in xrange(reaction.n_species)]
			pd = param_plot()
			pd.y_scale = 'log'
			pd.y_label = 'z'
			pt.multiplot2(tt,zz,pd)
	#
		plot_rates= inputs['input_plots'][2]
		if plot_rates is True:
			pd = param_plot()
			pd.y_scale = 'symlog'
			pd.y_label = 'p rate'
			pp = [reaction.time_series('rates',i) for i in xrange(n_s)]
			pt.multiplot2(tt,pp,pd)
	#
	#########################################
	# This section computes the slow manifold
	compute_sm = inputs['input_compute_sm']
	if compute_sm is True:
		cinit_min = np.zeros(n_obs+n_obs_square+n_obs_cubic)
		#cinit_min = np.r_[np.random.rand(n_obs)-0.5,10*(np.random.rand(n_obs_square)-0.5),100*(np.random.rand(n_obs_cubic)-0.5)]
		if inputs['input_minimization'] is None:a = ot.slow_manifold(cinit_min,dx = exp_dx, f_def_obs = observable, eps_ = eps_eff,time = timelog,rho=rho,method = 'BFGS',measure = mm,verbose = False,offset = offset_ts,maxiter=[1000,100],n_recuit = 10, isrecuit=False,r=reaction)
		else: a = inputs['input_minimization']
		if verbose is True:
			print a
			print cinit_min
		#ot.cost_gramian(a.x,f_def_obs = observable,rho=rho,dx=dx,eps_=eps_,time=timelog,offset=offset_ts,measure=mm,verbose=True)
		gram=ot.EG(fobs = lambda x:observable(x,a.x),dx=exp_dx,eps_=eps_eff,time=timelog,offset=offset_ts)
	else:
		a = inputs['input_minimization']


	###########################################################
	# This section just plots the reconstruction of one species
	isrec,isplot = True,inputs['input_plots'][3]
	
	if isrec is True:#reconstruction of the species
		offs='end' 
		ret = reconstruction(a.x,'',offset=offs,r=reaction,isplot=False)
		pd = param_plot()
		pd.marks = ['-']
		pd.colors = ['k']
		pd.y_scale = 'symlog'
		pd.y_label = 'err in reconstruction'
		pd.legend = [False]
		if isplot is True: pt.plot2((timelog[1:-1],np.abs(ret[2][1:-1]-ret[3][1:-1])/np.abs(ret[2][1:-1])),pd)
		if verbose is True:print J(a.x)
		is_a_random = False
		if is_a_random is True:
			abis = 1-2*np.random.rand(a.x.size)
			#print(J(abis))
			reconstruction(abis,'alpha random',r=reaction)

	############################################################
	# This section compares the results with a pure minimization
	min_test = False
	if min_test is True:#compare with minimization
		def jj(a):
			x=reaction.time_series('z',ind2rec) + a[0]*reaction.time_series('z',ind2use[0]) + a[1]*reaction.time_series('z',ind2use[1])
		#	x=np.dot(x[1:-1],np.diff(timelog))
			j=np.linalg.norm(x)
			return j

		from scipy.optimize import minimize as solve

		aa = solve(jj,np.zeros(2),method ='Powell')
		tt = (timelog,timelog)
		yy = (reaction.time_series('z',ind2rec)[1:], -aa.x[0]*reaction.time_series('z',ind2use[0])[1:] - aa.x[1]*reaction.time_series('z',ind2use[1])[1:]) 
		pd = pt.Paradraw()
		p=pt.Paradraw()
		pd.title = 'test'
		pd.x_scale = 'log'
		pd.x_label = 'time'
		pd.marks = ['-']
		pd.colors = ['k']
		pd.y_label = 'z'
		pd.y_scale = 'log'
		pd.ylim = False
		pd.x_tick_label = [1.*10**i for i in xrange(tmini,tmaxi,2)]
		#pd.x_tick_label = [1.*10**i for i in xrange(-9,-1,2)]
		pd.xlim = [1e-9,1e-3]
		pd.thickness = [2]
		pd.ylim = False
		pd.legend=[False]
		pd.y_scale = 'symlog'

		pt.multiplot2(tt,yy,pd)

	###########################################################
	# This section just plots the reconstruction of one species
	#, with (strongly) different initial conditions
	rand_test,rand_plot = inputs['input_compute_random'],inputs['input_plots'][4]
	if rand_test is True:
		if inputs['input_ci_random'] is not None:nr=1
		else:nr = 5
		ci_save=[]
		for irand in xrange(nr):
			#tpinit_random = (1467,101325)
			if inputs['input_tp_random'] is None:tpinit_random = tpinit
			else:tpinit_random = inputs['input_tp_random']
			reaction2 = chemreac(cti = mechanism, problem = prob,V=V)
			reaction2.TP(tp=tpinit_random,reset = False)
			if ci_random is None:
		 		ci2 = c_init*(np.ones(n_s) +  .90*(np.random.rand(n_s)-0.5)*2. )
				ci_save.append(ci2)
			else:
				ci2 = ci_random
			#ci2 = np.array([0.0151,0.0033,0.0177,0.0147,0.0051,0.0023,0.0049,0.0131])
			reaction2.concentrations(c=ci2,reset = True)
			reaction2.integrate(time = timelog)
			#ret2 = reconstruction(a.x,'ci random',offset=offs,r=reaction2)
			ret2 = reconstruction(a.x,'',r=reaction2,isobs=False,isplot=rand_plot)
			if verbose is True:
				print 'spec mole: ' + str(reaction2.h_z[0])
			if ci_random is not False:
				obj_ = (timelog,ret,ret2,a,(n_obs,n_obs_square,n_obs_cubic))
 				#dt.save(filename=spec2rec + '.rec',obj=obj_)
	else:
		ret2 = None
		reaction2=None

	name = 'plop'
	path = '/home/fgueniat/Documents/productions/UIUC/XPACC_chat/janvier/text'

	gram =lambda a: ot.EG(fobs = lambda x:observable(x,a),dx=exp_dx,eps_=eps_eff,time=timelog,offset=offset_ts)
	rr = lambda a: reconstruction(a,r=reaction,isplot=False)

	isderiv = False
	if isderiv is True:
		def testderiv(err):
			dd = np.logspace(-10,-1,100)
			jj = [(J(cinit_min + ddd)-J(cinit_min))/ddd for ddd in dd]
			pd = pt.Paradraw()
			pd.x_scale = 'log'
			pd.y_scale = 'symlog'
			pd.y_label = 'dj/dalpha'
			pd.x_label = 'dalpha'
			pt.plot2((dd,jj),pd)
			fobs = lambda z: observable(z,a.x)
			y1 = np.array([fobs(z) for z in reaction.h_z[1:-1]])
			fobs = lambda z: observable(z,a.x+err)
			y2 = np.array([fobs(z) for z in reaction.h_z[1:-1]])
			pd.marks = ['','-']
			pd.markers = ['o','']
			pd.colors = ['r','k']
			pd.y_label = 'obs'
			pd.x_label = 'time'
			pd.title = str(np.linalg.norm((y1-y2)*np.diff(timelog))/err/np.diff(timelog).mean())
			pt.multiplot2((timelog[1:],timelog[1:]),(y1,y2),pd)
			return jj,y1,y2
		dj,y1,y2=testderiv(1.e-2)
#
	def getW(a):
		fobs = lambda z: observable(z,a)
		w = ot.EG(fobs,exp_dx,eps_eff,timelog,offset = offset_ts,r=reaction)
		return w
#
	issurface = False
	if issurface is True:
		nx,ny = 50,50
		dx,dy = 1.e-3,1.e-3
		
		indx,indy=np.argsort(a.x)[0:2]
		indx,indy = 0,1
		jj = np.zeros((nx,ny))
		x = np.linspace(a.x[indx]-dx,a.x[indx]+dx,nx)
		y = np.linspace(a.x[indy]-dy,a.x[indy]+dy,ny)
		for ix,xx in enumerate(x):
			print 'pourcentage: ' + '%.3f' %(100.*ix/nx)
			for iy,yy in enumerate(y):
				ax = np.copy(a.x)
				ax[indx] = xx
				ax[indy] = yy
				jj[ix,iy] = J(ax)
		pd=pt.Paradraw()
		pd.x_label = spec2use[indx]
		pd.y_label = spec2use[indy]
		pt.pcolor((y,x,np.log10(np.abs(jj))),pd)
#
	ret_ = {
		'reaction':reaction,'retour':ret,
		'reaction_pert':reaction2,'retour_pert':ret2,
		'result_min':a,
		'indexes':ind2use,
		'time':timelog,
		'obs_index':obs_sorted,
		'gramian':w,
		'parameters':inputs
		}
	return ret_
	#return reaction,ret,reaction2,ret2,a,ind2use[0:n_obs],timelog,obs_sorted

if __name__ == "__main__":
	inputs = params(sys.argv)
	main(inputs=inputs)
