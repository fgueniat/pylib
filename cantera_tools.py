import cantera as ct
import numpy as np
'''
API for cantera
'''

class TwoWayDict(dict):
        ''' Dictionary that works both way, as a bijection
        '''
        def __init__(self,dic={}):
                self.merge(dic)
        def __setitem__(self, key, value):
                # Remove any previous connections with these values
                if key in self:
                        del self[key]
                if value in self:
                        del self[value]
                dict.__setitem__(self, key, value)
                dict.__setitem__(self, value, key)

        def __delitem__(self, key):
                dict.__delitem__(self, self[key])
                dict.__delitem__(self, key)

        def __len__(self):
                """Returns the number of connections"""
                return dict.__len__(self) // 2

        def merge(self,dic):
                """merge the dic with another one"""
                for key in dic:
                        self[key] = dic[key]
                

class chemreac():
    '''
    Class for managing a reaction object from cantera
    '''
    def __init__(self,cti='gri30.cti',TP=False,V=1.e-3, concentrations = False, X=False, Y=False, problem=('0D','volume')):
        ''' 
        initialization.
            cti string, cti file
            TP tuple, temperature and pressure
            V real, initial volume
            concentrations array of reals, initial concentrations
            X 
        '''
        #define the problem             
        self.problem = problem
        #define gaz
        self.cti = cti
        self.gas = ct.Solution(cti)
        self.V = V
        self.n_species = self.gas.n_species
        self.n_elements = self.gas.n_elements
        self.construct_E()
        #initialization
        if TP is not False:self.TP(TP)
        if X is not False: self.X(x=X)
        if Y is not False: self.Y(y=Y)
        if concentrations is not False: self.concentrations(c=concentrations)
        if (X is False) and (Y is False) and (concentrations is False):
            #print 'Initial conditions: iso molar fraction'
            self.Y(y=np.ones(self.n_species)) 
            # dictionaries for the species

        dd = {}
        for k,species in enumerate(self.gas.species()):
            dd[k] = species.name
        self.list = TwoWayDict(dd)

        self.reset()
        
    def construct_E(self):
        '''
        Construct matrix of mass
        '''
        self.E = np.zeros((self.n_species,self.n_elements))
        for i_e in range(self.n_elements):
            for i_s in range(self.n_species):
                n = self.gas.n_atoms(i_s,i_e)
                #if n>0:n=1./n
                self.E[i_s,i_e] = n
            #self.E[i_e,:] = self.E[i_e,:]/sum(self.E[i_e,:])
                    
    def reset(self):
        self.time = [0]
        self.h_concentrations = [self.concentrations()]
        self.h_T = [self.temperature()]
        self.h_P = [self.pressure()]
        self.h_x = [self.X()]
        self.h_y = [self.Y()]
        self.h_z = [self.Z()]
        try:
            self.h_volume = [self.reactor.volume]
        except:
            self.h_volume = [self.V]
        self.h_rates = [self.rates()]
        self.tpinit = self.gas.TP
        self.set_case(self.problem)
        self.conv_y2c= self.h_concentrations[0]/(self.h_y[0]/(self.h_volume[0]*self.gas.molecular_weights))
    def y2c(self,y=None,v=None):
        '''
        concert mol frac to concentration
        '''
        if y is None: y = self.Y()
        if v is None:v = self.volume()
        c = self.conv_y2c * y / (v*self.gas.molecular_weights)
        return c
    def y2t(self,y=None,v=None,p=None):
        '''
        compute the temperature, using the law of perfect gas
        '''
        if y is None: y = self.Y()
        if v is None:v = self.volume()
        if p is None: p = self.pressure()
        c = self.y2c(y=y,v=v)
        T = p/(np.sum(c)*ct.gas_constant)
        return T
    def set_case(self,problem):
        ''' function that set up the scenario'''
        if isinstance(problem[0],str) is False: problem[0] = '0D'
        if isinstance(problem[1],str) is False: problem[1] = 'volume'
        if problem[0].lower() == '0d':
            if problem[1].lower() == 'pressure':
                self.reactor = ct.IdealGasConstPressureReactor(self.gas)
                self.problemis = '0D, iso bar'
            if problem[1].lower() == 'volume':
                self.reactor = ct.IdealGasReactor(self.gas)
                self.reactor.volume = self.V
                self.problemis = '0D, iso volume'
            self.system = ct.ReactorNet((self.reactor,))
                
    def step(self,t=1,save=True):
        '''
        save the state of the system
        '''
        #       dt = self.system.step(t)
        #   self.time.append(self.time[-1]+dt)
        tp = self.system.step(t)
        if save is True:
            self.time.append(tp)
            self.h_concentrations.append(self.concentrations())
            self.h_x.append(self.X())
            self.h_y.append(self.Y())
            self.h_z.append(self.Z())
            self.h_T.append(self.temperature())
            self.h_P.append(self.pressure())
            self.h_rates.append(self.rates())
                
    def integrate(self,time = None,save=True):
        '''
            advance to time.
            if time is an array {t_1,...,t_N}, the method will take advance the system to each t_i
        '''

        if time is None: return -1
#               dt = self.system.step(t)
#               self.time.append(self.time[-1]+dt)
        for t in time:
            try:
                self.system.advance(t)
            except:
                self.system.set_max_time_step(1)
                self.system.advance(t)
            if save is True:
                self.time.append(t)
                self.h_concentrations.append(self.concentrations())
                self.h_x.append(self.X())
                self.h_y.append(self.Y())
                self.h_z.append(self.Z())
                self.h_T.append(self.temperature())
                self.h_P.append(self.pressure())
                self.h_volume.append(self.volume())
                self.h_rates.append(self.rates())
                        
                            
    def which(self,i=None):
        '''
            convert number to species
        '''
        if i is None:
            for specie in range(self.n_species):
                self.which(i=specie)
        else:
            if isinstance(i,int):
                s = 'the ' + str(i) + 'th species is ' + self.list[i]
                print(s)
            else:
                s = 'the species ' + i + ' is associated with number ' + str(self.list[i])
                print(s)

    def concentrations(self,c=None,reset = True):
        '''
        set or get the concentrations
        '''
        if c is None: 
                return self.gas.concentrations
        else: 
                if len(c) == self.n_species:
                        self.Y(y=c/np.sum(c),reset = reset)
                else:
                        print('some species are missing')

    def TP(self,tp=None,reset = True):
        '''
        set or get the temperature and pressure as a tuple.
        '''
        if tp is None: 
                return self.gas.TP
        else:
                if len(tp) == 2:
                        self.gas.TP = tp
                        if reset is True:self.reset()
                else: 
                        print('tuplet (T,P) is needed')
                    
    def X(self,x=None,reset=True):
        ''' set or get the specific moles
        '''
        return self.fraction_mol(x=x,reset=reset)
    def fraction_mol(self,x=None,reset=True):
        if x is None: 
                return self.gas.X
        else:
                if len(x) == self.n_species:
                        self.gas.X = x
                        if reset is True:self.reset()
                else: 
                        print('some species are missing')

    def fraction_mass(self,y=None,reset=True):
        if y is None: 
                return self.gas.Y
        else:
                if len(y) == self.n_species:
                        self.gas.Y = y
                        if reset is True:self.reset()
                else: 
                        print('some species are missing')
    def Y(self,y=None,reset=True):
        return self.fraction_mass(y=y,reset=reset)

    def specific_mole(self,z=None,reset=True):
        if z is None: 
                return self.gas.Y/self.gas.molecular_weights
        else:
                if len(z) == self.n_species:
                        self.gas.Y = z*self.gas.molecular_weights
                        if reset is True:self.reset()
                else: 
                        print('some species are missing')
    def Z(self,z=None,reset=True):return self.specific_mole(z=z,reset=reset)


    def temperature(self): return self.gas.T
    def pressure(self): return self.gas.P
    def volume(self): return self.reactor.volume
    def rates(self): return self.gas.net_production_rates
    def time_series(self,quantity=False,i=0):
        # check for errors in inputs
        if isinstance(i,int) is False:i=0
        if (quantity is False) or (isinstance(quantity,str) is False):
                print('temperature')
                return self.ts_temperature(i=i)

        if quantity.lower() == 'concentrations': return self.ts_concentrations(i=i)             
        elif quantity.lower() == 'x': return self.ts_X(i=i)             
        elif quantity.lower() == 'y': return self.ts_Y(i=i)             
        elif quantity.lower() == 'z': return self.ts_Z(i=i)             
        elif quantity.lower() == 'rates': return self.ts_rates(i=i)
        elif quantity.lower() == 't': return self.ts_temperature()              
        elif quantity.lower() == 'temperature': return self.ts_temperature()            
        elif quantity.lower() == 'p': return self.ts_pressure()         
        elif quantity == 'pressure': return self.ts_pressure()          
        elif quantity == 'volume': return self.ts_volume()              
        else: return self.ts_temperature()

    def ts_pressure(self):return np.array(self.h_P).flatten()
    def ts_volume(self):return np.array(self.h_volume).flatten()
    def ts_temperature(self):return np.array(self.h_T).flatten()
    def ts_concentrations(self,i=0):return np.array(self.h_concentrations)[:,i]
    def ts_X(self,i=0):return np.array(self.h_x)[:,i]
    def ts_Y(self,i=0):return np.array(self.h_y)[:,i]
    def ts_Z(self,i=0):return np.array(self.h_z)[:,i]
    def ts_rates(self,i=0):return np.array(self.h_rates)[:,i]
    def plot(self,style='z'):
        import plot_tools as pt
        pd = pt.Paradraw()
        pd.colors = pt.get_colors(self.n_species)
        pd.x_scale = 'symlog'
        pd.y_scale = 'symlog'
        pd.x_label = 't'
        if style == 'rates':
                pd.y_label = 'rates'
                zs = [self.ts_rates(i)[1:] for i in range(self.n_species)]
        else:
                pd.y_label = 'z'
                zs = [self.ts_Z(i)[1:] for i in range(self.n_species)]
        time = (self.time[1:],)*self.n_species
        figs = pt.multiplot2(time,zs,pd)
        return figs

    def stoich(self,fuel='H2',oxidizer='O2'):
        ind_ox = self.list[oxidizer]
        try:
            ind_fuel = self.list[fuel]
        except:
            print('fuel does not exist')
            return -1
        ratio_oxfuel = self.h_y[0][ind_ox]/self.h_y[0][ind_fuel]
        print('Ox/Fuel ratio by mass: ' + str(ratio_oxfuel))
        if fuel == 'H2': st = 34.3
        elif fuel == 'CH4': st = 17.19
        else:
            return -1
        print('Fuel/ox equivalent ratio' + str(ratio_oxfuel/st))

    
    def __getitem__(self,key):
        return self.list[key]
        



def compute_T(P,V,zs,r):
    R = 8.31
    T = np.sum([P*V/(n*R) for n in ns])
    return T

