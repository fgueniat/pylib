# biblio plot3
#import pylab
from matplotlib import pyplot
import pylab as pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mtick
from os.path import expanduser
from matplotlib.colors import LogNorm
from scipy.stats import gaussian_kde
import graphviz
from itertools import cycle
#import la_tools as la
##################################################################################
#########################  Parameters  ###########################################
##################################################################################
path = '/home/fgueniat/Documents/productions/'
cycol = cycle('bgrcmk')

class Paradraw():
    def __init__(self,marks=['-'], colors = ['k'], markers = [''], thickness=[1],x_label = 'axis 1',y_label = 'axis 2',z_label = 'axis 3',c_label = 'colors', colmap = 'hot_r',xlim = False,ylim = False,zlim = False,clim = False, x_scale = 'linear',y_scale = 'linear',z_scale = 'linear',c_scale = 'linear',xlintresh=1.e-10,ylintresh=1.e-10,zlintresh=1.e-10,clintresh=1.e-10,title='',iscolorbar = False,fontsize = 28,ticksize = 22, markerssize = [3], fillstyle=['none'], x_tick_label = False,y_tick_label = False,z_tick_label = False,cbar_tick_label = False,ax_x_format = '%.2f',ax_y_format = '%.2f',ax_z_format = '%.2f',cbformat = '%.2f',figure=False, bar_col = [0.0,0.0,1.0],transparancy = False,ncont = 15,pcol_type = 'contourf',tight_layout=True, stem = False,legend=[False],legloc = 'best',axis=None,antialiasing=False,levels = False):

#plot style
        self.marks = marks 
        self.colors = colors
        self.markers = markers
        self.fillstyle = fillstyle
        self.markerssize = markerssize
        self.thickness = thickness
        self.markeredge = False
        self.colmap = colmap
        self.ncont = ncont
        self.pcol_type = pcol_type
        self.bar_col = bar_col
        self.transparancy = transparancy
        self.tight_layout = tight_layout
        self.stem  = stem
        self.axis = axis
        self.fig,self.ax,self.cax,self.cbar = None,None,None,None
        self.antialiasing=antialiasing
        self.levels = levels

#legend
        self.legend=legend
        self.legloc = legloc
#axis labels
        self.x_label = x_label
        self.y_label = y_label
        self.z_label = z_label
        self.c_label = c_label
        self.x_tick_label = x_tick_label
        self.y_tick_label = y_tick_label
        self.z_tick_label = z_tick_label
        self.cbar_tick_label = cbar_tick_label
        self.title = title
#axis style
        self.x_scale = x_scale
        self.y_scale = y_scale
        self.z_scale = z_scale
        self.c_scale = c_scale
        self.xlintresh = xlintresh
        self.ylintresh = ylintresh
        self.zlintresh = zlintresh
        self.clintresh = clintresh
#axis limits
        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim
        self.clim = clim

        self.iscolorbar = iscolorbar
#fonts
        self.fontsize = fontsize
        #ax_x_format='%.2f',
        #ax_y_format = '%.0f',
        self.ticksize = ticksize
        self.cbformat = cbformat
        self.ax_x_format = ax_x_format 
        self.ax_y_format = ax_y_format
        self.ax_z_format = ax_z_format
        #other exemples : '%.2f'
        self.figure = figure

        self.error_collected = []
    def help(self):
        return self.__dict__    
    
    def __len__(self):
        return 1
    def __getitem__(self, key):
        return self
    
    @property
    def help(self):
        return self.__dict__
    @property
    def r(self): self.figure=False
    
    @property
    def closeall(self): closeall()


def remove_inf(x,dim=None):
    x = np.array(x)
    if dim is None:
        try:
            dim = x.shape[1]
        except:
            dim = None
    if dim is None:
        m = np.isfinite(x)
        return x[m]
    else:
        M = np.isfinite(x[:,0])
        for i in range(1,dim):
            m = np.isfinite(x[:,i])
            M = (M) & (m)
        return x[M,:]
    
    
##################################################################################
############################  PLOT 2D  ###########################################
##################################################################################

def rawplot(data,param,fig,ax):
    x = data[0]
    y = data[1]
    if param is False:
        param = Paradraw()
    param.iscolorbar = False
    if param.markeredge is False:
        ax.plot( x, y, linestyle = param.marks[0], marker = param.markers[0], ms = param.markerssize[0], color=param.colors[0],markeredgecolor='none',linewidth = param.thickness[0])
    else:
        ax.plot( x, y, linestyle = param.marks[0], marker = param.markers[0], ms = param.markerssize[0], color=param.colors[0],linewidth = param.thickness[0])

    # set the limits of the plot to the limits of the data
    Options(ax,(x,y),param)

    return fig,ax

def plot1(data,param=False):
    if param is False:
        param = Paradraw()
    param.iscolorbar = False
    x = np.arange(len(data))
    y = data
    data2 = (x,y)
    plot2(data2,param)
    return param

def plot2(data,param=False):
    x = data[0]
    y = data[1]
    if param is False:
        param = Paradraw()
    param.iscolorbar = False
    pylab.ion()
    if param.figure is False:
        fig = pylab.figure()
        ax = fig.add_subplot(111)
    else:
        fig = param.figure[0]
        ax = param.figure[1]
    if param.markeredge is False:
        if param.stem is False:ax.plot( x, y, linestyle = param.marks[0], marker = param.markers[0], fillstyle=param.fillstyle[0], ms = param.markerssize[0], color = param.colors[0],markeredgecolor='none',linewidth = param.thickness[0],label=param.legend[0])
        else:ax.stem( x, y, linestyle = param.marks[0], marker = param.markers[0], fillstyle=param.fillstyle[0], ms = param.markerssize[0], color = param.colors[0], markeredgecolor='none',linewidth = param.thickness[0],label=param.legend[0])
    else:
        if param.stem is False:ax.plot( x, y, linestyle = param.marks[0], marker = param.markers[0], fillstyle=param.fillstyle[0], ms = param.markerssize[0], color = param.colors[0], linewidth = param.thickness[0],label=param.legend[0])
        else:ax.stem( x, y, linestyle = param.marks[0], marker = param.markers[0], fillstyle=param.fillstyle[0], ms = param.markerssize[0], color = param.colors[0], linewidth = param.thickness[0],label=param.legend[0])
    # set the limits of the plot to the limits of the data
    param.fig=fig
    param.ax = ax
    Options((x,y),param)
    for err in param.error_collected:
        import matplotlib.pyplot as plt
        if err[0] == 'xticks':plt.xticks(err[1],err[2])
        if err[0] == 'yticks':plt.yticks(err[1],err[2])
        
    return param


def multiplot1(y,param=False):
    pylab.ion()
    if param is False:
        param = Paradraw()
    param.iscolorbar = False
    figure = param.figure

    x = (np.array(range(1,len(y[0])+1)),)
    for i in range(1,len(y)):
        x = x + (np.array(range(1,len(y[i])+1)),)
    if figure is False:
        fig, ax = pyplot.subplots()
    else:
        fig = figure[0]
        ax = figure[1]
    param.figure = (fig,ax)
    param = multiplot2(x,y,param)
    return param

def multiplot2(x,y,param=False):

    if param is False:
        param = Paradraw()
    param.iscolorbar = False
    if param.xlim is False:
        param.xlim = [np.nanmin(x[0]),np.nanmax(x[0])]
        for i in range(1,len(x)):
            m=np.nanmin(x [i])
            if m < param.xlim[0]: param.xlim[0] = m
            m=np.nanmax(x[i])
            if m > param.xlim[1]: param.xlim[1] = m 
    if param.ylim is False:
        param.ylim = [np.nanmin (y[0]),np.nanmax(y[0])]
        for i in range(1,len(y)):
            m=np.nanmin(y[i]) 
            if m < param.ylim[0]: param.ylim[0] = m
            m=np.nanmax(y[i])
            if m > param.ylim[1]: param.ylim[1] = m 

    figure = param.figure

    pylab.ion()
    if figure is False:
        fig, ax = pyplot.subplots()
    else:
        fig = figure[0] 
        ax = figure[1]
    if param.transparancy is False:
        transparancy = 1.00
    else:
        transparancy = 0.3
    for i in range(0,len(y)):

        j = len(y)-i-1 
        if j==0:
            transparancy = 1.0
        try:
            sty = param.marks[j]
        except:
            sty = '-'
        try:
            mi = param.markers[j]
        except:
            mi = ''
        try:
            ci = param.colors[j]
        except:
            ci = next(cycol)
        try:
            linewidth = param.thickness[j]
        except:
            try:
                linewidth = param.thickness[0]
            except:
                linewidth = 3
        try:
            leg=param.legend[j]
        except:
            leg=False
        try:
            ms = param.markerssize[j]
        except:
            ms = ms = param.markerssize[0]
        if np.sum(np.isnan(y[j]))<len(y[j]):
            if param.markeredge is False:
                if param.stem is False:
                    ax.plot(x[j],y[j],linestyle = sty, marker = mi, ms=ms, color = ci, markeredgecolor='none',linewidth = linewidth,alpha=transparancy,label=leg)
                else: 
                    ax.stem(x[j],y[j],linestyle = sty, marker = mi, ms=ms, color = ci, markeredgecolor='none',linewidth = linewidth,alpha=transparancy,label=leg)
            else:
                if param.stem is False:
                    ax.plot( x[j],y[j],linestyle = sty, marker = mi, ms=ms, color = ci, linewidth = linewidth,alpha=transparancy,label=leg)
                else:
                    ax.stem(x[j],y[j],linestyle = sty, marker = mi,ms=ms, color = ci,linewidth = linewidth,alpha=transparancy,label=leg)
    
    ax.set_xlabel(param.x_label)
    ax.set_ylabel(param.y_label)
    param.fig = fig
    param.ax = ax
    Options((nanminmax(x),nanminmax(y)),param)

    return param
def nanminmax(xx):
    xm,xM = 1.e30,-1.e30
    for x in xx:
        if np.nanmin(x)<xm:xm = np.nanmin(x)
        if np.nanmax(x)>xM:xM = np.nanmax(x)
    return [xm,xM]


##################################################################################
############################  PLOT 3D  ###########################################
##################################################################################
            
def plot3(x,y,z, param = False,figure=False):
    if param is False:
        param = Paradraw()
    param.iscolorbar = False
    if param.xlim is False:
        param.xlim = [np.nanmin(x[0]),np.nanmax(x[0])]
        for i in range(1,len(x)):
            m=np.nanmin(x [i])
            if m < param.xlim[0]: param.xlim[0] = m
            m=np.nanmax(x[i])
            if m > param.xlim[1]: param.xlim[1] = m 
    if param.ylim is False:
        param.ylim = [np.nanmin (y[0]),np.nanmax(y[0])]
        for i in range(1,len(y)):
            m=np.nanmin(y[i]) 
            if m < param.ylim[0]: param.ylim[0] = m
            m=np.nanmax(y[i])
            if m > param.ylim[1]: param.ylim[1] = m 

    figure = param.figure

    pylab.ion()
    if figure is False:
        fig, ax = pyplot.subplots()
    else:
        fig = figure[0] 
        ax = figure[1]
    if param.transparancy is False:
        transparancy = 1.00
    else:
        transparancy = 0.3
    for i in range(0,len(y)):

        j = len(y)-i-1 
        if j==0:
            transparancy = 1.0
        try:
            sty = param.marks[j]
        except:
            sty = '-'
        try:
            mi = param.markers[j]
        except:
            mi = ''
        try:
            fi = param.fillstyle[j]
        except:
            fi = 'none'
        try:
            ci = param.colors[j]
        except:
            ci = next(cycol)
        try:
            linewidth = param.thickness[j]
        except:
            try:
                linewidth = param.thickness[0]
            except:
                linewidth = 3
        try:
            leg=param.legend[j]
        except:
            leg=False
        try:
            ms = param.markerssize[j]
        except:
            ms = ms = param.markerssize[0]
        if np.sum(np.isnan(y[j]))<len(y[j]):
            if param.markeredge is False:
                if param.stem is False:
                    ax.plot(x[j],y[j],linestyle = sty, marker = mi, fillstyle=fi, ms=ms, color = ci, markeredgecolor='none',linewidth = linewidth,alpha=transparancy,label=leg)
                else: 
                    ax.stem(x[j],y[j],linestyle = sty, marker = mi, fillstyle=fi, ms=ms, color = ci, markeredgecolor='none',linewidth = linewidth,alpha=transparancy,label=leg)
            else:
                if param.stem is False:
                    ax.plot( x[j],y[j],linestyle = sty, marker = mi, fillstyle=fi, ms=ms, color = ci, linewidth = linewidth,alpha=transparancy,label=leg)
                else:
                    ax.stem(x[j],y[j],linestyle = sty, marker = mi, fillstyle=fi, ms=ms, color = ci,linewidth = linewidth,alpha=transparancy,label=leg)
    
    ax.set_xlabel(param.x_label)
    ax.set_ylabel(param.y_label)
    param.fig = fig
    param.ax = ax
    Options((nanminmax(x),nanminmax(y)),param)

    return param
def nanminmax(xx):
    xm,xM = 1.e30,-1.e30
    for x in xx:
        if np.nanmin(x)<xm:xm = np.nanmin(x)
        if np.nanmax(x)>xM:xM = np.nanmax(x)
    return [xm,xM]


##################################################################################
############################  PLOT 3D  ###########################################
##################################################################################
            
def plot3(x,y,z, param = False,figure=False):

    if param is False:
        param = Paradraw()
    param.iscolorbar = False
    pylab.ion()
    if figure is False:
        #fig = pylab.figure()
        #ax = Axes3D(fig)
        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = figure[0]
        ax = figure[1]
    mi = param.markers[0]
    fi = param.fillstyle[0]
    sty = param.marks[0]
    ci = param.colors[0]
    ms = param.markerssize[0]
#        if param.markeredge is False:
#        ax.plot(a[j], b[j], c[j],mci, markeredgecolor='none')
#        else:
    if param.markeredge is False:
        ax.plot(x, y, z, linestyle = sty, color = ci, marker = mi, fillstyle=fi, ms = ms, markeredgecolor='none',label=param.legend[0])
    else:
        ax.plot(x, y, z, linestyle = sty, color = ci, marker = mi, fillstyle=fi, ms = ms, label=param.legend[0])
    param.fig = fig
    param.ax = ax
    Options((x,y,z),param)
    param.ax.xaxis.labelpad = 12
    param.ax.yaxis.labelpad = 14
    param.ax.zaxis.labelpad = 10
    pylab.show()
    return param


def multiplot3(x,y,z,figure=False, param = False):
    if param is False:
        param = Paradraw()
    pylab.ion()
    param.iscolorbar = False
    if figure is False:
        #fig = pylab.figure()
        #ax = Axes3D(fig)        
        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = figure[0]
        ax = figure[1]
    for i in range(0,len(x)):
        j = len(x)-i-1
        try:
            mi = param.markers[j]
        except:
            mi = ''
        try:
            fi = param.fillstyle[j]
        except:
            fi = 'none'
        try:
            sty = param.marks[j]
        except:
            sty = '-'
        try: 
            ci = param.colors[j]
        except:
            ci = next(cycol)
        try:
            ms = param.markerssize[j]
        except:
            ms = param.markerssize[0]
        try:
            leg=param.legend[j]
        except:
            leg=False
#        if param.markeredge is False:
#        ax.plot(a[j], b[j], c[j],mci, markeredgecolor='none')
#        else:
        if param.markeredge is False:
            ax.plot(x[j], y[j ], z[j],linestyle = sty, marker = mi, fillstyle = fi, ms = ms, color = ci, markeredgecolor='none',label=leg)
        else:
            ax.plot(x[j], y[j], z[j], linestyle = sty, marker = mi, fillstyle = fi, ms = ms, color = ci,label=leg)

#        ax.plot(x[j], y[j], z[j],mci,label=leg)

#        ax.plot(x[j], y[j], z[j],mci)
    param.ax = ax
    param.fig = fig
    Options((x[0],y[0],z[0]),param)    
    param.ax.xaxis.labelpad = 12
    param.ax.yaxis.labelpad = 14
    param.ax.zaxis.labelpad = 10
    pylab.show()
    return param

##################################################################################
##########################  PLOT COLOR  ##########################################
##################################################################################

def fill(x,y,y_plus,y_minus,param=False,color="#006699"):
    if param is False:
        param = Paradraw()
    param = plot2((x,y),param)
    #param.figure = (fig,ax)
    param.ax.fill_between(x, y_plus,y_minus, color=color)
    #Options((x,y),param)
    return param

def multifill(xs,ys,y_plus,y_minus,param=False,colors=["#006699",'#FF4136']):
    if param is False:
        param = Paradraw()
    param = multiplot2(xs,ys,param)
    #param.figure = (fig,ax)
    ym = nanminmax( [ nanminmax(ys), nanminmax(y_plus), nanminmax(y_minus) ] )
    Options((nanminmax(xs),ym),param)
    n=len(ys)
    for i in range(n):
        try:
            color = colors[i]
        except:
            color= colors[0]
        #
        param.ax.fill_between(xs[i], y_plus[i],y_minus[i], color=color,alpha=0.4)
    return param

def pcolor(data,param=False):
    if param is False:
        param = Paradraw()

    if len(data)==3:
        x = data[0] 
        y = data[1]
        z = data[2]
    else:
        z = data+0.0 
        z = np.flipud(z)
        y = np.arange(z.shape[0])
        x = np.arange(z.shape[1])
    if param.zlim is False:
        if param.clim is False:  
            zmin = np.nanmin(z) 
            zmax = np.nanmax(z)
            param.clim = [zmin,zmax]
        else:
            zmin = param.clim[0] 
            zmax = param.clim[1]
            param.clim = [zmin,zmax]
    else:
        zmin = param.zlim[0] 
        zmax = param.zlim[1]

    pylab.ion()
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    if param.axis is not None:
        ax.axis(param.axis)
    if param.clim is not False:
        cmin = param.clim[0] 
        cmax = param.clim[1]
    else:
        cmin = zmin 
        cmax = zmax

    if param.pcol_type == 'pcolor':
        if param.c_scale == 'log': 
            cax = ax.pcolor(x,y,z, cmap=param.colmap, norm=LogNorm(vmin=cmin, vmax=cmax), vmin=cmin, vmax=cmax)
        else:
            cax = ax.pcolor(x,y,z, cmap=param.colmap, vmin=zmin, vmax=cmax)
    else:
        if param.c_scale == 'log': 
            if param.levels is False:levels = np.logpace(cmin,cmax,param.ncont)
            else:levels=param.levels
            cax = ax.contourf(x,y,z, param.ncont, antialiased = param.antialiasing, cmap=param.colmap, norm=LogNorm(vmin=cmin, vmax=cmax), vmin=cmin, vmax=cmax,levels=levels,extend="both")
        else:
            if param.levels is False:levels = np.linspace(cmin,cmax,param.ncont)
            else:levels=param.levels
            cax = ax.contourf(x,y,z, param.ncont, antialiased = param.antialiasing, cmap=param.colmap, vmin=cmin, vmax=cmax,levels=levels,extend="both")
    param.cax = cax
    # set the limits of the plot to the limits of the data
    param.fig = fig
    param.ax = ax
#        fig.delaxes(fig.axes[1]) #cbar IS axes[1] while ax is axes[0]
    Options((x,y),param)
    pylab.show()
    return param

def plot3c(figure,a,b,c,z,a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):
    p = np.argsort(z)
    print(p.shape)
    print(a.shape)
    print(b.shape)
    print(c.shape)
    print(p.max())
    macol = pyplot.cm.bwr(np.arange(len(z)))
    print(macol.shape)
    pylab.ion()
    if figure is False:
        fig = pylab.figure() 
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = figure[0] 
        ax = figure[1]
    ax.scatter(a[p], b[p], c[p], c=macol,edgecolor='')
    ax.set_xlabel(a_label)
    ax.set_ylabel(b_label)
    ax.set_zlabel(c_label)

    return (fig,ax)

def plot3cluster(a,b,c,keys,codomain,order,a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):
    pylab.ion()
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    colbase = pyplot.cm.jet(np.arange(256))
    colbase = colbase[np.int64(np.linspace(0,255,len(codomain))),:]
    for i in range(len(codomain)):
        p = np.where(keys==codomain[order[i]-1])[0]
        macol = np.array([np.random.random(),np.random.random(),np.random.random(),1.])
        macol = colbase[i,:]
        ax.scatter(a[p], b[p], c[p], c=macol,edgecolor='')
    ax.set_xlabel(a_label)
    ax.set_ylabel(b_label)
    ax.set_zlabel(c_label)

    return (fig,ax)

def plot3cluster2(X,clusters,keys,colors,param=False):
    
    if param is False:
        param = Paradraw()
    pylab.ion()
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')


    for p in range(len(colors)):
        try:
            mtype=param.markers[p]
        except:
            mtype=param.markers[0]
        #      
        try:
            s=param.markerssize[p]
        except:
            s=param.markerssize[0]
        #
        ax.scatter(X[p,0], X[p,1], X[p,2], c=colors[p],edgecolor='none',s=s,marker=mtype)
    for key in keys:
        k= [key[kk] for kk in range(key.size)] 
        xx = []
        norm_xx = []
        for i in range(len(clusters)):
            clust=clusters[i]
            c= [clust[kk] for kk in range(len(clust))]
            if c == k:
                col = colors[i]
                xx.append(X[i])
            #
        #
        xx = np.array(xx)
        try:
            axx.plot_trisurf(xx[:,0], xx[:,1], xx[:,2], linewidth=0.0, antialiased=False,edgecolor='none',color = col,shade=False)
        except:
            pass
    param.fig = fig
    param.ax = ax
    param.iscolorbar = False
    Options((X[:,0],X[:,1],X[:,2]),param)
    return (fig,ax)


def density2D(x,y,param=False):
    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    try:x, y, z = x[idx], y[idx], z[idx]
    except:x, y, z = np.array(x)[idx], np.array(y)[idx], np.array(z)[idx]


    if param is False:
        param = Paradraw()
    pylab.ion()
    if param.figure is False:
        fig, ax = pyplot.subplots()
        #fig = pylab.figure()
        #ax = Axes3D(fig)
    else:
        fig = figure[0]
        ax = figure[1]
    ax.scatter(x, y, c=z, s=50, edgecolor='')
    param.ax = ax
    param.fig = fig
    Options((x,y),param)
    
    return param

##################################################################################
#########################  PLOT STATS  ###########################################
##################################################################################

def pwelch(data=False,dt=False,param=False,normalize = None,normalize_around = None,dnorm=None,strouhal = 1.,isplot=True):
    '''
    normalize:
        False: nothing
        'absolute': normalize
        float: by the max aroud float
    '''
    if data is False:
        print('no data')
        return -1
    import scipy.fftpack
    from scipy import signal
    # Number of samplepoints
    if param is False:
        param=Paradraw()
        param.ax_x_format = '%.3f'
        param.ax_y_format = '%.3f'
        param.y_scale = 'log'
        param.ylabel = 'PSD'
        param.xlabel = 'f'
    N = len(data)
    # sample spacing
    if dt is False:dt = 1.
    df = 1./dt

    xf, yf = signal.welch(data, df, nperseg=1024)
    xf = xf*strouhal

    if normalize == True:
        yf = yf/np.nanmax(yf)
    if normalize_around is not None:
        from numbers import Number
        if isinstance(normalize_around,Number):
            ind = np.argmin(np.abs(xf-normalize_around))
            if isinstance(dnorm,Number):
                df = np.argmin(np.abs(xf -  dnorm  ) )
                if df>ind:df = ind
                print(df,ind)
            else:
                df = int(0.1*normalize_around*(N*dt)/strouhal) # N*dt = 1/df
            indm = np.nanmax([0,ind-df])
            indM = np.nanmin([N,ind+df])
            rn = np.nanmax(yf[indm:indM])
            yf = yf/rn
    if isplot is True:
        fig,ax = plot2((xf, yf),param)
        param.ax = ax
        param.fig = fig
        return param
    else:
        return (xf,yf)

def fft(data=False,time=False,param=False,normalize = False,strouhal=1.,isplot=True):
    '''
    normalize:
        False: nothing
        True: normalize at 1
        float: by the max aroud float
    '''
    if data is False:
        print('no data')
        return -1
    import scipy.fftpack
    # Number of samplepoints
    if param is False:
        param=Paradraw()
        param.ax_x_format = '%.3f'
        param.ax_y_format = '%.3f'
        param.ylabel = 'PSD'
        param.xlabel = 'f'
    
        
    N = len(data)
    # sample spacing
    if time is False:time = np.linspace(0,1,N)
    dt = time[1]-time[0]
    yf = scipy.fftpack.fft(data)
    xf = np.linspace(0.0, 1.0/(2.0*dt), N/2)*strouhal
#    print len(xf)
#    print len(yf)
    yf = 2.0/N * np.abs(yf[:N//2])
    if normalize is True:
        yf = yf/np.nanmax(yf)
    elif normalize is False:
        pass
    else:
        from numbers import Number
        if isinstance(normalize,Number):
            ind = np.argmin(np.abs(xf-normalize))
            df = int(0.1*normalize*(N*dt)/strouhal) # N*dt = 1/df
            indm = np.nanmax([0,ind-df])
            indM = np.nanmin([N,ind+df])
            rn = np.nanmax(yf[indm:indM])
            yf = yf/rn
            
    if isplot is True:
        fig,ax = plot2((xf, yf),param)
        param.fig = fig
        param.ax = ax
        return param
    else:
        return (xf,yf)

def plot3c(figure,a,b,c,z,a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):
    p = np.argsort(z)
    print(p.shape)
    print(a.shape)
    print(b.shape)
    print(c.shape)
    print(p.max())
    macol = pyplot.cm.bwr(np.arange(len(z)))
    print(macol.shape)
    pylab.ion()
    if figure is False:
        fig = pylab.figure() 
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig = figure[0] 
        ax = figure[1]
    ax.scatter(a[p], b[p], c[p], c=macol,edgecolor='')
    ax.set_xlabel(a_label)
    ax.set_ylabel(b_label)
    ax.set_zlabel(c_label)

    return (fig,ax)

def plot3cluster(a,b,c,keys,codomain,order,a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):
    pylab.ion()
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    colbase = pyplot.cm.jet(np.arange(256))
    colbase = colbase[np.int64(np.linspace(0,255,len(codomain))),:]
    for i in range(len(codomain)):
        p = np.where(keys==codomain[order[i]-1])[0]
        macol = np.array([np.random.random(),np.random.random(),np.random.random(),1.])
        macol = colbase[i,:]
        ax.scatter(a[p], b[p], c[p], c=macol,edgecolor='')
    ax.set_xlabel(a_label)
    ax.set_ylabel(b_label)
    ax.set_zlabel(c_label)

    return (fig,ax)

def plot3cluster2(X,clusters,keys,colors,param=False):
    
    if param is False:
        param = Paradraw()
    pylab.ion()
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')


    for p in range(len(colors)):
        try:
            mtype=param.markers[p]
        except:
            mtype=param.markers[0]
        #      
        try:
            s=param.markerssize[p]
        except:
            s=param.markerssize[0]
        #
        ax.scatter(X[p,0], X[p,1], X[p,2], c=colors[p],edgecolor='none',s=s,marker=mtype)
    for key in keys:
        k= [key[kk] for kk in range(key.size)] 
        xx = []
        norm_xx = []
        for i in range(len(clusters)):
            clust=clusters[i]
            c= [clust[kk] for kk in range(len(clust))]
            if c == k:
                col = colors[i]
                xx.append(X[i])
            #
        #
        xx = np.array(xx)
        try:
            axx.plot_trisurf(xx[:,0], xx[:,1], xx[:,2], linewidth=0.0, antialiased=False,edgecolor='none',color = col,shade=False)
        except:
            pass
    param.fig = fig
    param.ax = ax
    param.iscolorbar = False
    Options((X[:,0],X[:,1],X[:,2]),param)
    return (fig,ax)


def density2D(x,y,param=False):
    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    try:x, y, z = x[idx], y[idx], z[idx]
    except:x, y, z = np.array(x)[idx], np.array(y)[idx], np.array(z)[idx]


    if param is False:
        param = Paradraw()
    pylab.ion()
    if param.figure is False:
        fig, ax = pyplot.subplots()
        #fig = pylab.figure()
        #ax = Axes3D(fig)
    else:
        fig = figure[0]
        ax = figure[1]
    ax.scatter(x, y, c=z, s=50, edgecolor='')
    param.ax = ax
    param.fig = fig
    Options((x,y),param)
    
    return param

##################################################################################
#########################  PLOT STATS  ###########################################
##################################################################################

def pwelch(data=False,dt=False,param=False,normalize = None,normalize_around = None,dnorm=None,strouhal = 1.,isplot=True,nperseg=1024):
    '''
    normalize:
        False: nothing
        'absolute': normalize
        float: by the max aroud float
    '''
    if data is False:
        print('no data')
        return -1
    import scipy.fftpack
    from scipy import signal
    # Number of samplepoints
    if param is False:
        param=Paradraw()
        param.ax_x_format = '%.3f'
        param.ax_y_format = '%.3f'
        param.y_scale = 'log'
        param.ylabel = 'PSD'
        param.xlabel = 'f'
    N = len(data)
    # sample spacing
    if dt is False:dt = 1.
    df = 1./dt

    xf, yf = signal.welch(data, df, nperseg=nperseg)
    xf = xf*strouhal

    if normalize == True:
        yf = yf/np.nanmax(yf)
    if normalize_around is not None:
        from numbers import Number
        if isinstance(normalize_around,Number):
            ind = np.argmin(np.abs(xf-normalize_around))
            if isinstance(dnorm,Number):
                df = np.argmin(np.abs(xf -  dnorm  ) )
                if df>ind:df = ind
                print(df,ind)
            else:
                df = int(0.1*normalize_around*(N*dt)/strouhal) # N*dt = 1/df
            indm = np.nanmax([0,ind-df])
            indM = np.nanmin([N,ind+df])
            rn = np.nanmax(yf[indm:indM])
            yf = yf/rn
    if isplot is True:
        param = plot2((xf, yf),param)
        return param
    else:
        return (xf,yf)

def fft(data=False,time=False,param=False,normalize = False,strouhal=1.,isplot=True):
    '''
    normalize:
        False: nothing
        True: normalize at 1
        float: by the max aroud float
    '''
    if data is False:
        print('no data')
        return -1
    import scipy.fftpack
    # Number of samplepoints
    if param is False:
        param=Paradraw()
        param.ax_x_format = '%.3f'
        param.ax_y_format = '%.3f'
        param.ylabel = 'PSD'
        param.xlabel = 'f'
    
        
    N = len(data)
    # sample spacing
    if time is False:time = np.linspace(0,1,N)
    dt = time[1]-time[0]
    yf = scipy.fftpack.fft(data)
    xf = np.linspace(0.0, 1.0/(2.0*dt), N/2)*strouhal
#    print len(xf)
#    print len(yf)
    yf = 2.0/N * np.abs(yf[:N//2])
    if normalize is True:
        yf = yf/np.nanmax(yf)
    elif normalize is False:
        pass
    else:
        from numbers import Number
        if isinstance(normalize,Number):
            ind = np.argmin(np.abs(xf-normalize))
            df = int(0.1*normalize*(N*dt)/strouhal) # N*dt = 1/df
            indm = np.nanmax([0,ind-df])
            indM = np.nanmin([N,ind+df])
            rn = np.nanmax(yf[indm:indM])
            yf = yf/rn
            
    if isplot is True:
        fig,ax = plot2((xf, yf),param)
        param.fig = fig
        param.ax = ax
        return param
    else:
        return (xf,yf)

def hist(data,n=20,tresh = False,param = False):
    data2 = data[~np.isnan(data)]
    if tresh is not False:
        try:
            tresh_min = tresh[0]
            tresh_max = tresh[1]
        except:
            tresh_max = tresh
            tresh_min = -tresh
        data2[data2>tresh_max] = tresh_max
        data2[data2<tresh_min] = tresh_min
    if param is False:
        print('default parameters') 
        param = Paradraw()
        param.marks = ['k','or']
        param.x_label = 'bins'
        param.y_label = 'Frequency'
        param.ax_x_format = '%.2f'
        param.ax_y_format = '%.0f'
        param.iscolorbar = False
#        if is3d is True:
#            param.ax_z_format = '%.0f'
    hist, bins = np.histogram(data2, bins=n)
#    hist = hist
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    pylab.ion()
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.bar(center, hist, align='center', width=width, color = param.bar_col)
    param.ax = ax
    param.fig = fig
    Options((center,hist),param)
    pylab.show()
    return param,(hist,bins)

def plot_Mahalanobis(d,r,param=False):
    if param is False:
        print('default parameters') 
        param = Paradraw()
        param.marks = ['k','or']
        param.x_label = 'Normal theoritical quantiles'
        param.y_label = 'Ordered Mahalanobis distances'
        param.ax_x_format = '%.0f'
        param.ax_y_format = '%.0f'

    x2=np.linspace(  np.nanmin( d[0] ),np.nanmax( d[0] ), len(d[0])  )
    y2 = r[1] + r[0]*x2

    x = ( x2, d[0] )
    y = ( y2, d[1] )
    multiplot2(x,y,param)
    return param
##################################################################################
############################  OPTIONS  ###########################################
##################################################################################

def Options(X,param=False):
    ax   = param.ax
    fig  = param.fig
    cbar = param.cbar
    cax  = param.cax
    x=X[0]
    y=X[1]
    try:
        z=X[2]
        if cbar is None:
            is3d=True
        else:
            is3d = False
    except:
        is3d = False
#might have a bug with colorbar
    if param is False:
        param = Paradraw()
    if param.xlim is False:
        ax.set_xlim(  [ np.nanmin(x),np.nanmax(x) ]  )
        xm = np.nanmin(x)
        xM = np.nanmax(x)
    else:
        ax.set_xlim(  [ param.xlim[0],param.xlim[1] ]  )
        xm = param.xlim[0]
        xM = param.xlim[1]
    if param.ylim is False:
        ym = np.nanmin(y)
        yM = np.nanmax(y)
        ax.set_ylim(  [ ym,yM ]  )
    else:
        ax.set_ylim(  [ param.ylim[0],param.ylim[1] ]  )
        ym = param.ylim[0]
        yM = param.ylim[1]
    if is3d is True:
        try:
            if param.zlim is False:
                zm = np.nanmin(z)
                zM = np.nanmax(z)
                ax.set_zlim(  [ zm,zM ]  )
            else:
                ax.set_zlim(  [ param.zlim[0],param.zlim[1] ]  )
                zm = param.zlim[0]
                zM = param.zlim[1]
        except:
            pass
    #set scales
    ax.set_xscale(param.x_scale)
    ax.set_yscale(param.y_scale)
    #ax.set_xscale(param.x_scale,linthreshx = param.xlintresh)
    #ax.set_yscale(param.y_scale,linthreshy = param.ylintresh)
    if is3d is True:
        try:
            ax.set_zscale(param.z_scale)
        except:
            pass
    ax.set_xlabel(param.x_label,None,None,fontsize=param.fontsize)
    ax.set_ylabel(param.y_label,None,None,fontsize=param.fontsize)
    if is3d is True:
        try:
            ax.set_zlabel(param.z_label,None,None,fontsize=param.fontsize)
        except:
            pass

    if param.x_tick_label is None:
        ax.xaxis.set_visible(False)
    else:
        if param.x_tick_label is False:
            if param.x_scale == 'log':
                xlabel = np.exp(np.linspace(np.log(xm),np.log(xM),3))
            elif param.x_scale == 'symlog':
                if xm<0:xlabel = np.array([xm,0,xM])
                else:xlabel = np.exp(np.linspace(np.log(xm),np.log(xM),3))
            else:
                xlabel = np.linspace(xm,xM,3)
        else:
            xlabel = param.x_tick_label
        try:
            ax.set_xticks(xlabel)
        except:
            param.error_collected.append(['xticks',x,xlabel])
        ax.set_xticklabels(ax.get_xticks(),None,None,fontsize=param.ticksize)
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter(param.ax_x_format))

    if param.y_tick_label is None:
        ax.yaxis.set_visible(False)
    else:
        if param.y_tick_label is False:
            if param.y_scale == 'log':
                ylabel = np.exp(np.linspace(np.log(ym),np.log(yM),3))
            elif param.y_scale == 'symlog':
                if ym<0:
                    ylabel = np.array([ym,0,yM])
                else:
                    ylabel = np.exp(np.linspace(np.log(ym),np.log(yM),3))
            else:
                ylabel = np.linspace(ym,yM,3)
        else:
            ylabel = param.y_tick_label
        try:
            ax.set_yticks(ylabel)
        except:
            param.error_collected.append(['yticks',y,ylabel]) 
        ax.set_yticklabels(ax.get_yticks(),None,None,fontsize=param.ticksize)
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter(param.ax_y_format))

    if is3d is True:
        try:
            if param.z_tick_label is None:
                ax.zaxis.set_visible(False)
            else:
                if param.z_tick_label is False:
                    if param.z_scale == 'log':
                        zlabel = np.exp(np.linspace(np.log(np.nanmin(z)),np.log(np.nanmax(z)),3))
                    else:
                        zlabel = np.linspace(zm,zM,3)
                else:
                    zlabel = param.z_tick_label
                ax.set_zticks(zlabel)
                ax.set_zticklabels(ax.get_zticks(),fontsize=param.ticksize)
                ax.zaxis.set_major_formatter(mtick.FormatStrFormatter(param.ax_z_format))
        except:
            pass

    if is3d is True:
        try:
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False
        except:
            pass

    ax.set_title(param.title)
    if param.iscolorbar is True:
        if is3d is True:
            c=X[2]
        else:
            c=X[1]        
        cm,cM = np.nanmin(c),np.nanmax(c)
        if param.cbar_tick_label is False:
            if param.clim is False:
                clabel = np.linspace(np.nanmin(c),np.nanmax(c),3)
            else:
                clabel = np.linspace(param.clim[0],param.clim[1],3)
                cm,cM = param.clim[0],param.clim[1]
        else:
            clabel = param.cbar_tick_label
        #cbar = f.colorbar(cf, ticks=[3000,4000,5000,6000])
        #cbar.ax.set_ylim([cbar.norm(3000), cbar.norm(6000)])
        #cbar.outline.set_ydata([cbar.norm(3000)] * 2 + [cbar.norm(6000)] * 4 + [cbar.norm(3000)] * 3)
        #cbar.ax.set_aspect(60)
#        cbar.set_ticklabels(clabel)
        cbar = fig.colorbar(cax,format=param.cbformat,ax = ax, ticks = clabel)
        cbar.set_label(param.c_label,fontsize=param.fontsize)
        for t in cbar.ax.get_yticklabels():
            t.set_fontsize(param.ticksize)
#        cbar.set_clim(cm, cM)
#        cbar.set_ticks(clabel)
#        cbar.set_ticklabels(clabel)
#        cbar.ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))


    if is3d is True:
        if param.z_tick_label is None: ax.zaxis.set_visible(True)

    try:
        if all(legend is False for legend in param.legend) is False:
            leg = ax.legend(loc=param.legloc,fontsize = 'large')
            leg.draggable(state=True)
    except:
        pass
                
    if is3d is False:    
        if param.tight_layout is True:
            try:
                pyplot.tight_layout()
            except:
                print('*********************************************************')
                print('not tight')
                print('most probably, figure parameter in Paradraw has to be set to False')
                ax.tight_layout()
#    ax.tick_params(axis='both', which='minor', labelsize=param.ticksize)




def change_options():
    ax = pyplot.gca()
#    ax.tick_params(axis='both', which='minor', labelsize=param.ticksize)
    ax.xaxis.set_minor_locator(pyplot.FixedLocator([2,4]))

##################################################################################
###########################  TOOLS  ##############################################
##################################################################################

def closeall(n=50):
    for i in range(n):
        pyplot.close()
def close():
    pyplot.close()

def subplot(y,x=False,params=False):
    import copy
    if params is False:
        params = Paradraw()
    n = len(y)
    f, a = pyplot.subplots(nrows=n, ncols=1)
    for i in range(n):
        params[i].tight_layout = False
        ax = a[i]
        params[i].figure = (f,ax)
        if x is False:xx = np.arange(len(y[i]))
        else: xx = x[i]
        plot2((xx,y[i]),param = params[i])

def fig_list():
    figlist = pyplot.get_fignums()
    if len(figlist) == 0:
        return 0
    else:
        return figlist
        
def save(path = None,name = None, fig = None):
    ''' Save a fig'''
    if path is None:
        path = expanduser("~") + '/Desktop/'
    if path[-1] != '/':
        path = path + '/'
    if name is None:
        import datetime
        name = 'fig_' + datetime.datetime.now().strftime('%Y%m%d-%H%M%S')

    if fig is not None:
        try:
            pyplot.savefig(path + name + fig,transparent=True)
        except:
            pyplot.savefig(path + name + '.eps',transparent=True)
            pyplot.savefig(path + name + '.pdf',transparent=True)
            pyplot.savefig(path + name + '.png',transparent=True)
    else:
        pyplot.savefig(path + name + '.eps',transparent=True)
        pyplot.savefig(path + name + '.pdf',transparent=True)
        pyplot.savefig(path + name + '.png',transparent=True)

##################################################################################
###########################  DIVERS  #############################################
##################################################################################
def markov(A,filename = None,param = None,states = None):
    '''bind to plot_markov'''
    return plot_markov2(A,filename=filename,param=param,states=states)



def plot_markov2(A,filename = None,param = None,states = None):
    '''
    Plot a markov process.
    A is the transition matrix.
    filename to save the file
    param: dictionnary accepting:
        'treshold':0,
        'fontsize':8,
        'max_edges':2,
        'color_edge':'red',
        'format':'eps',
        'save':True
    '''
    param_init = {'treshold':0,
        'fontsize':8,
        'max_edges':2,
        'color_edge':'red',
        'format':'eps',
        'save':True,
        'ratio':.7}    
    if param is None: param = {}
    for key in param_init.keys():
        if key not in param.keys():param[key]=param_init[key]   
    #
    if filename is None:
        filename = expanduser("~") + '/Desktop/graph.' + param['format']
    #
    f = graphviz.Graph('finite_state_machine', filename=filename)
    f.attr(rankdir='LR',size='7,5')
    f.attr('node', shape='doublecircle',fontsize = str(param['fontsize']))
    f.attr('edge',fontsize = str(param['fontsize']),color=param['color_edge'])
    nstate = A.shape[0]
    if param['max_edges'] is None or param['max_edges'] < 1:
        m_edge = nstate
    else:
        m_edge = param['max_edges']
    #
    if states is None:
        states = [str(i) for i in range(nstate)]
    else:
        if len(states)!=nstate:
            print('not enough states')
            return -1
    for state,name_state in enumerate(states):
        f.node(name_state)
        to_draw = np.argsort(A[state])[::-1][:m_edge]
        for to_state in range(state,nstate):
            to_draw_reciprocal = np.argsort(A[to_state])[::-1][:m_edge] 
            #  
            is_draw = False        
            is_left = False
            is_right = False
            if A[state,to_state]>param['treshold']:
                if to_state in to_draw:
                    is_right = True
                    is_draw = True
                #
            #
            if A[to_state,state]>param['treshold']:
                if state in to_draw_reciprocal:
                    is_left = True
                    is_draw = True
                #
            #
            if is_left is True and is_right is True:
                is_double = True
            else:
                is_double = False
            #
            if is_draw is True:
                if is_double is True:
                    f.attr('edge',dir='both')
                    if state == to_state:label_str = '%.2f' % A[state,to_state]
                    else:label_str = '%.2f' % A[to_state,state] + ', ' + '%.2f' % A[state,to_state]
                else:
                    if is_right is True:
                        f.attr('edge',dir='forward')
                        label_str = '%.2f' % A[state,to_state]
                    else:
                        f.attr('edge',dir='back')
                        label_str = '%.2f' % A[to_state,state]
                    #
                #
                #f.edge(str(state), str(to_state), label = label_str)
                f.edge(states[state], states[to_state], label = label_str)


    f.view()
    if param['save'] is True:
        file_ = open(filename, "w")
        file_.write(f.pipe(param['format']))
        file_.close()
    return f



def plot_markov(A,filename = None,param = None,states = None):
    '''
    Plot a markov process.
    A is the transition matrix.
    filename to save the file
    param: dictionnary accepting:
        'treshold':0,
        'fontsize':8,
        'max_edges':2,
        'color_edge':'red',
        'format':'eps',
        'save':True
    '''
    param_init = {'treshold':0,
        'fontsize':8,
        'max_edges':2,
        'color_edge':'red',
        'format':'eps',
        'save':True}    
    if param is None: param = {}
    for key in param_init.keys():
        if key not in param.keys():param[key]=param_init[key]   
    #
    if filename is None:
        filename = expanduser("~") + '/Desktop/graph.' + param['format']
    #
    f =graphviz.Graph('finite_state_machine', filename=filename)
    f.attr(rankdir='LR', size='8,5')
    f.attr('node', shape='doublecircle')
    f.attr('edge',fontsize = str(param['fontsize']),color=param['color_edge'])
    nstate = A.shape[0]
    if param['max_edges'] is None or param['max_edges'] < 1:
        m_edge = nstate
    else:
        m_edge = param['max_edges']
    #
    if states is None:
        states = [str(i) for i in range(nstate)]
    else:
        if len(states)!=nstate:
            print('not enough states')
            return -1
    for state,name_state in enumerate(states):
        f.node(name_state)
        perm = np.argsort(A[state,:])[::-1]
        for to_state in range(state,nstate):
            perm2 = np.argsort(A[to_state,:])[::-1] 
            #  
            is_draw = False        
            is_left = False
            is_right = False
            if A[state,to_state]>param['treshold']:
                if perm[to_state]<m_edge:
                    is_right = True
                    is_draw = True
                #
            #
            if A[to_state,state]>param['treshold']:
                if perm2[state]<m_edge:
                    is_left = True
                    is_draw = True
                #
            #
            if is_left is True and is_right is True:
                is_double = True
            else:
                is_double = False
            #
            if is_draw is True:
                if is_double is True:
                    f.attr('edge',dir='both')
                    if state == to_state:label_str = '%.2f' % A[state,to_state]
                    else:label_str = '%.2f' % A[to_state,state] + ', ' + '%.2f' % A[state,to_state]
                else:
                    if is_right is True:
                        f.attr('edge',dir='forward')
                        label_str = '%.2f' % A[state,to_state]
                    else:
                        f.attr('edge',dir='back')
                        label_str = '%.2f' % A[to_state,state]
                    #
                #
                #f.edge(str(state), str(to_state), label = label_str)
                f.edge(states[state], states[to_state], label = label_str)


    f.view()
    if param['save'] is True:
        file_ = open(filename, "w")
        file_.write(f.pipe(param['format']))
        file_.close()
    return f

def get_colors(n,cname = 'hot'):
    import matplotlib
    cmap = matplotlib.pyplot.get_cmap(cname,n)
    colors = [cmap(i) for i in range(n)]
    return colors

def get_cmap(n=1):
    import matplotlib.colors as colors
    cols = []
    for ic,c in enumerate(colors.CSS4_COLORS):
        if ic<n:cols.append(str(c))
    return cols

def flim1(figure,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):
    if figure is False:
        fig = pyplot.figure( 1 )
        ax = fig.add_subplot( 111 )
        im = ax.imshow( np.zeros( ( 256, 256, 3 ) ) ) # Blank starting image
        figure = (fig,ax,im)
        fig.draw()
        im.axes.figure.canvas.draw()
    else:
        fig = figure[0]
        ax = figure[1]
        im = figure[2]
    im.set_data( b )
    im.axes.figure.canvas.draw()
    return (fig,ax,im)


#ani = matplotlib.animation.FuncAnimation(fig, animate, init_func=init_animation, frames=50)
#ani.save('test_anim.gif', writer='imagemagick', fps=30)





