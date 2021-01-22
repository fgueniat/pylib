import numpy as np
import matplotlib.pylab as plt
from scipy import stats
from collections import *
from itertools import *

#
# Convenience Class for Time Series
# 
# Series Objects can be added and subtracted using infix operators.
# > X = Series([1,2,3])
# > X + 10 -> [11,12,13]
#

import collections

def IterMap(f, I):
    for y in I: yield f(y)

def IterAdd(x, I):
    return IterMap(lambda y: y+x, I)

def IterMul(x, I):
    return IterMap(lambda y: x*y, I)

def IterIterAdd(I, J):
    while True:
        yield I.next() + J.next()

def IterIterMul(I, J):
    while True:
        yield I.next() * J.next()

class Series(collections.Iterator):
    I = None

    def __init__(self, I):
        # convert list to iterator
        self.I = iter(I)
        
    def next(self): return self.I.next()
    
    def __add__(self, x):
        if isinstance(x, (int, long, float, complex)):
            return Series(IterAdd(x, self))
        elif isinstance(x, Series):
            return Series(IterIterAdd(x, self))
        else: raise TypeError()

    def __mul__(self, x):
        if isinstance(x, (int, long, float, complex)):
            return Series(IterMul(x, self))
        elif isinstance(x, Series):
            return Series(IterIterMul(x, self))
        else: raise TypeError()
    
    def __sub__(self, x): return __add__(self, -x)

    def __div__(self, x): return __mul__(self, 1./x)
    
    __radd__ = __add__
    __rsub__ = __sub__
    __rmul__ = __mul__

    
def ToSeries(func):
    def SeriesWrapper(*args,**kwargs):
        return Series(func(*args,**kwargs))
    return SeriesWrapper
    
    
def Sample(I,N=1000):
    return np.array([ y for y in islice(I,N) ])

def Plot(I, *args, **kwargs):
    N = kwargs.pop("N", 1000)

    plt.plot(Sample(I,N),*args, **kwargs)
    
def Hist(I, *args, **kwargs):
    N = kwargs.pop("N", 1000)
    kwargs['bins'] = kwargs.get('bins', np.sqrt(N))

    H = plt.hist(Sample(I,N), *args, **kwargs)

def FancyPlot(I, *args, **kwargs):
    N = kwargs.pop("N", 1000)
    y = Sample(I,N)

    # rects [left, bottom, width, height]
    rect_time = [0,  0, 0.9, 1]
    rect_histy= [0.905,0, 0.095, 1]

    # start with a rectangular Figure
    plt.figure(1, figsize=(20,4))

    axHisty = plt.axes(rect_histy)
    axTime  = plt.axes(rect_time)

    # create plots
    Plot(Series(y),N=N, *args, **kwargs)
    axHisty.hist(y, orientation='horizontal', bins=np.sqrt(N))

    # align y axis
    if (max(y) == axTime.get_ylim()[1]):
        axTime.set_ylim((axTime.get_ylim()[0], max(y) * 1.1))
    axHisty.set_ylim( axTime.get_ylim() )
