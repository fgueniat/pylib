import pickle
import numpy as np
import scipy as sp

def identify_inflection(time_series,grid = False,verbose=True):
    x = np.arange(len(time_series))
    spline = sp.interpolate.UnivariateSpline(x,time_series)
    if grid is not False:
        g = sp.interpolate.UnivariateSpline(x,grid)
    dspline = spline.derivative()
    def j(x):
        return -dspline(x)
    x0 = (x[-1]+x[0])/2.
    r = sp.optimize.minimize(j,x0)
    alpha = r.x - int(r.x)
    tm = (1-alpha)*time_series[int(r.x)] + alpha*time_series[int(r.x)+1]
    if grid is not False:
        if verbose is True:
            print('grid point: ' + str(r.x))
        return g(r.x),tm
    else:
        return r.x,tm

def split_list(str_):
        '''
        sys.argv to list
        '''
        temp = str_.split('[')
        temp = temp[1].split(']')
        return [float(s) for s in temp[0].split(',')]
def split_tuple(str_):
        '''
        sys.argv to tuple
        '''
        return tuple(split_list(str_))

def str_none(str_):
        pass
def print_(prec=3):
#       print('np.set_printoptions(precision=prec,suppress=True)')
        np.set_printoptions(precision=prec,suppress=True)

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
                

def Cleandic(twodic,list_):
        ''' clean a dic of useless keys
        '''
        td = TwoWayDict()
        for key in twodic:
                if key in list_:
                        try:
                                td[key] = twodic[key]
                        except:
                                pass
        return td


def save(obj = False, filename = False,verbose=True):
        if obj is False:
                print('nothing to save')
                return
        if filename is False:
                from datetime import date
                filename = 'data_' + date.today().isoformat() + '.dat'
        output = open(filename, 'wb')
        pickle.dump(obj,output)
        output.close()
        if verbose is True:print('file saved in:' + filename)
        return filename


def label2csv(y=False,x=False,filename =False,verbose = True,fmt = '%15.5e',delim = '\t'):
    '''
    dump labels in a csv type file.
    used to plot with tikz
    '''
    if y is False:
        print('nothing to save')
        return
    #    # Check if matrix
    # check dimensions
    ny=len(y)
    if x is False:x = np.arange(ny)
    # Check filename
    if filename is False:
        from datetime import date
        filename = 'data_' + date.today().isoformat() + '.csv'
    # start writing
    array_file = open(filename, "w")
    s = ''
    s += '{:>10}'.format('xlabel')
    s += '{:>10}'.format('ylabel')
    s += '\n'
    array_file.write(s)    
    for iy,pos_y in enumerate(y):
        s = ''
        s += '{:>10}'.format(x[iy])
        s += '{:>50}'.format(pos_y)
        s += '\n'
        array_file.write(s)
        #
    #
    array_file.close() 
    # verbose
    if verbose is True:print('file saved in:' + filename)
    #
    return filename


def arrays2csv(y=False,x=False,filename =False,verbose = True,fmt = '%15.5e',delim = '\t'):
    '''
    dump arrays in a csv type file.
    used to plot with tikz
    '''
    if y is False:
        print('nothing to save')
        return
    #    # Check if matrix
    try:
        ny = y.shape[0]
        if ny != y.size:
            print('y has to be an array, not a matrix')
            return -1
    except:
        print('y has to be a numpy array')
        return -1
    # check dimensions
    if x is False:x = np.arange(ny)
    # Check filename
    if filename is False:
        from datetime import date
        filename = 'data_' + date.today().isoformat() + '.csv'
    # check format    
    if fmt[0] == '%':fmt_ = fmt[1:]
    else: fmt_ = fmt
    #
    try:
        _ = '{a:{f}}'.format(a=np.sqrt(3),f=fmt_)
    except:
        print('invalid format')
        return -1
    # start writing
    array_file = open(filename, "w")
    for iy,pos_y in enumerate(y):
        s = ''
        s += '{x:{f}}'.format(x=x[iy],f=fmt_)
        s += '{y:{f}}'.format(y=pos_y,f=fmt_)
        s += '\n'
        array_file.write(s)
        #
    #
    array_file.close() 
    # verbose
    if verbose is True:print('file saved in:' + filename)
    #
    return filename


def matrix2csv(matrix = False,x=False,y=False, filename = False,verbose=True,fmt = '%15.5e',delim = '\t'):
    '''
    dump a matrix in a csv type file.
    used to plot with tikz
    '''
    if matrix is False:
        print('nothing to save')
        return
    #    # Check if matrix
    try:        
        nx,ny = matrix.shape
    except:
        print('matrix has to be a numpy matrix')
        return -1
    # check dimensions
    if x is False:x = np.arange(nx)
    if y is False:y = np.arange(ny)
    # Check filename
    if filename is False:
        from datetime import date
        filename = 'data_' + date.today().isoformat() + '.csv'
    # check format    
    if fmt[0] == '%':fmt_ = fmt[1:]
    else: fmt_ = fmt
    #
    try:
        _ = '{a:{f}}'.format(a=np.sqrt(3),f=fmt_)
    except:
        print('invalid format')
        return -1
    # start writing
    matrix_file = open(filename, "w")
    for iy,pos_y in enumerate(y):
        for ix,pos_x in enumerate(x):
            s = ''
            s += '{x:{f}}'.format(x=pos_x,f=fmt_)
            s += '{y:{f}}'.format(y=pos_y,f=fmt_)
            s += '{m:{f}}'.format(m=matrix[ix,iy],f=fmt_)
            s += '\n'
            matrix_file.write(s)
        #
        matrix_file.write('\n')
    #
    matrix_file.close() 
    # verbose
    if verbose is True:print('file saved in:' + filename)
    #
    return filename

def save2csv(obj = False, filename = False,verbose=True,fmt = '%10.5f',delim = '\t'):
        if obj is False:
                print('nothing to save')
                return
        if filename is False:
                from datetime import date
                filename = 'data_' + date.today().isoformat() + '.csv'
        np.savetxt(filename,obj,delimiter = delim,fmt=fmt)
        if verbose is True:print('file saved in:' + filename)
        return filename


def load(filename = False):
        if filename is False:
                print('nothing to load')
                return
        input_ = open(filename, 'rb')
        obj = pickle.load(input_)
        input_.close()
        return obj


def log_std(y,std,sign='+'):
    n = len(y)
    if sign == '-':
        sgn = -1.
    else:
        sgn = 1.
    #
    def f(yy,ss):
        r = ss/yy 
        x = np.log(yy) + sgn * 0.434 * r
        fx = np.exp(np.log(yy) + sgn * 0.434 * ss/yy)
        return fx
    logy = np.array(
            [f(y[i],std[i]) for i in xrange(n)]
            )
    return logy

def num2str_prec(n,prec=2):
    param = '.' +str(int(prec)) + 'f'
    s = format(n,param)
    return s

def n2s(n,prec=2):
    return num2str_prec(n,prec=prec)




def M_swap(i,j,n):
        M = np.eye(n)
        M[i,i] = 0
        M[j,j] = 0
        M[i,j] = 1
        M[j,i] = 1
        return M


def reorder(T):
    #
    N=T.shape[0]
    T_ordo = np.copy(T)
    perm = np.arange(N,dtype = int)
    for i in xrange(N-1):
        perm_temp = i+1 + np.argsort(T_ordo[i][i+1:])[::-1]
        #
        if perm_temp[0]>i+1:
            ind_temp = perm[i+1]
            perm[i+1] = perm_temp[0]
            perm[perm_temp[0]] = ind_temp
            #
            T_ordo = np.dot(T_ordo,M_swap(i+1,perm[i+1],N))
        #
    #
    return T_ordo,perm





