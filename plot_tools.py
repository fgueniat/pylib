# biblio plot3
import pylab
from matplotlib import pyplot
import pylab as pylab
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mtick
from os.path import expanduser
from matplotlib.colors import LogNorm



class Paradraw():
	def __init__(self,marks="-k",thickness=[1],x_label = 'axis 1',y_label = 'axis 2',z_label = 'axis 3',c_label = 'colors', colmap = 'hot_r',xlim = False,ylim = False,zlim = False,clim = False, x_scale = 'linear',y_scale = 'linear',z_scale = 'linear',c_scale = 'linear',xlintresh=0.,ylintresh=0.,zlintresh=0.,clintresh=0.,title='',iscolorbar = True,fontsize = 20,ticksize = 16, x_tick_label = False,y_tick_label = False,z_tick_label = False,cbar_tick_label = False,ax_x_format = '%.2e',ax_y_format = '%.2e',ax_z_format = '%.2f',cbformat = '%.2f',figure=False, bar_col = [0.0,0.0,1.0],transparancy = False,ncont = 15,pcol_type = 'contourf',tight_layout=True, stem = False,legend=[False],legloc = 'best'):


#plot style
		self.marks = marks
		self.mark = self.marks[0]
		self.thickness = thickness
		self.markeredge = True
		self.colmap = colmap
		self.ncont = ncont
		self.pcol_type = pcol_type
		self.bar_col = bar_col
		self.transparancy = transparancy
		self.tight_layout = tight_layout
		self.stem  = stem
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
		self.ticksize = ticksize
		self.cbformat = cbformat
		self.ax_x_format = ax_x_format 
		self.ax_y_format = ax_y_format
		self.ax_z_format = ax_z_format
		#other exemples : '%.2f'
		self.figure = figure
	def __len__(self):
		return 1
	def __getitem__(self, key):
		return self
	@property
	def help(self):
		return self.__dict__
	@property
	def r(self): self.figure=False

	def help(self):
		return self.__dict__
def help():
	s = "Paradraw: \n marks='-k',x_label = 'axis 1',y_label = 'axis 2',z_label = 'axis 3',c_label = 'colors', \n colmap = 'hot_r', \n xlim = False, ylim = False,zlim = False,clim = False, \n x_scale = 'linear',y_scale = 'linear',z_scale = 'linear', \n title='',iscolorbar = True,fontsize = 20,ticksize = 16, x_tick_label = False,y_tick_label = False,z_tick_label = False,cbar_tick_label = False,ax_x_format = '%.1e',ax_y_format = '%.1e',ax_z_format = '%.1e',cbformat = '%.1e', \n figure=False, \n bar_col = [0.0,0.0,1.0], \n transparancy = False, \n ncont = 15,pcol_type = 'contourf' "
	print(s)
				
def plot3c(figure,a,b,c,z,a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):
	p = np.argsort(z)
	print(p.shape)
	print(a.shape)
	print(b.shape)
	print(c.shape)
	print(p.max())
	macol = pyplot.cm.bwr(np.arange(z.size))
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

def plot3cluster(a,b,c,keys,codomain,order,figure=False,a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):
	pylab.ion()
	if figure is False:
		fig = pylab.figure()
		ax = fig.add_subplot(111, projection='3d')
	else:
		fig = figure[0]
		ax = figure[1]
	colbase = pyplot.cm.jet(np.arange(256))
	colbase = colbase[np.int64(np.linspace(0,255,codomain.size)),:]
	for i in range(codomain.size):
		p = np.where(keys==codomain[order[i]-1])[0]
		macol = np.array([np.random.random(),np.random.random(),np.random.random(),1.])
		macol = colbase[i,:]
		ax.scatter(a[p], b[p], c[p], c=macol,edgecolor='')
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)
	ax.set_zlabel(c_label)

	return (fig,ax)


def plot3old(a,b,c,figure=False,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2',c_label = 'axis 3'):

	pylab.ion()

	if figure is False:
		fig = pylab.figure()
		ax = Axes3D(fig)
	else:
		fig = figure[0]
		ax = figure[1]
	ax.plot(a, b, c,mark_col)
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)
	ax.set_zlabel(c_label)
	return (fig,ax)

def plot3(x,y,z, param = False,figure=False):

	if param is False:
		param = Paradraw()
	pylab.ion()
	if figure is False:
		fig = pylab.figure()
		ax = Axes3D(fig)
	else:
		fig = figure[0]
		ax = figure[1]
	mci = param.marks[0]
#		if param.markeredge is False:
#		ax.plot(a[j], b[j], c[j],mci, markeredgecolor='none')
#		else:
	if param.markeredge is False:
		ax.plot(x, y, z, mci, markeredgecolor='none',label=param.legend[0])
	else:
		ax.plot(x, y, z, mci, label=param.legend[0])
	Options(ax,(x,y,z),param)
	return (fig,ax)


#def plot2(figure,a,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

#	pylab.ion()
#	if figure is False:
#		fig, ax = pyplot.subplots()
#	else:
#		fig = figure[0]
#		ax = figure[1]
##	fig = pylab.figure()
#	pylab.plot(a, b,mark_col)
#	return (fig,ax)

#def plot1(figure,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

#	a = range(0,np.size(b))
#	if figure is False:
#		fig, ax = pyplot.subplots()
#	else:
#		fig = figure[0]
#		ax = figure[1]
#	figure = (fig,ax)
#	plot2(figure,a,b,mark_col,a_label,b_label)
#	return  (fig,ax)


def multiplot1bis(figure,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):
	a = (np.array(range(0,b[0].size)),)
	for i in range(1,len(b)):
		a = a + (np.array(range(0,b[i].size)),)
	if figure is False:
		fig, ax = pyplot.subplots()
	else:
		fig = figure[0]
		ax = figure[1]
	figure = (fig,ax)
	figure = multiplot2(figure,a,b,mark_col,a_label,b_label)
	return figure

def multiplot1(y,param=False):
	pylab.ion()
	if param is False:
		param = Paradraw()
	figure = param.figure

	x = (np.array(range(0,y[0].size)),)
	for i in range(1,len(y)):
		x = x + (np.array(range(0,y[i].size)),)
	if figure is False:
		fig, ax = pyplot.subplots()
	else:
		fig = figure[0]
		ax = figure[1]
	param.figure = (fig,ax)
	figure = multiplot2(x,y,param)
	return figure

def multiplot2bis(figure,a,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):


	pylab.ion()
	if figure is False:
		fig, ax = pyplot.subplots()
	else:
		fig = figure[0]
		ax = figure[1]

	for i in range(0,len(b)):
		j = len(b)-i-1
		try:
			mci = mark_col[j]
		except:
			mci = '-k'
		ax.plot(a[j],b[j],mci)
	ax.set_xlabel(a_label)
	ax.set_ylabel(b_label)

	return (fig,ax)


def multiplot2(x,y,param=False):

	if param is False:
		param = Paradraw()
	if param.xlim is False:
		param.xlim = [np.min(x[0]),np.max(x[0])]
		for i in xrange(1,len(x)):
			m=np.min(x[i])
			if m < param.xlim[0]: param.xlim[0] = m
			m=np.max(x[i])
			if m > param.xlim[1]: param.xlim[1] = m 
	if param.ylim is False:
		param.ylim = [np.min(y[0]),np.max(y[0])]
		for i in xrange(1,len(y)):
			m=np.min(y[i])
			if m < param.ylim[0]: param.ylim[0] = m
			m=np.max(y[i])
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
			mci = param.marks[j]
		except:
			mci = '-k'
		try:
			linewidth = param.thickness[j]
		except:
			linewidth = 2
		try:
			leg=param.legend[j]
		except:
			leg=False
		if param.markeredge is False:
			if param.stem is False:ax.plot(x[j],y[j],mci, markeredgecolor='none',linewidth = linewidth,alpha=transparancy,label=leg)
			else: ax.stem(x[j],y[j],mci, markeredgecolor='none',linewidth = linewidth,alpha=transparancy,label=leg)
		else:
			if param.stem is False:ax.plot(x[j],y[j],mci,linewidth = linewidth,alpha=transparancy,label=leg)
			else:ax.stem(x[j],y[j],mci,linewidth = linewidth,alpha=transparancy,label=leg)
	
	ax.set_xlabel(param.x_label)
	ax.set_ylabel(param.y_label)
	Options(ax,(x[0],y[0]),param)

	return (fig,ax)

def multiplot3(x,y,z,figure=False, param = False):
	if param is False:
		param = Paradraw()
	pylab.ion()
	if figure is False:
		fig = pylab.figure()
		ax = Axes3D(fig)
	else:
		fig = figure[0]
		ax = figure[1]
	for i in range(0,len(x)):
		j = len(x)-i-1
		try:
			mci = param.marks[j]
		except:
			mci = '-k'
		try:
			leg=param.legend[j]
		except:
			leg=False
#		if param.markeredge is False:
#		ax.plot(a[j], b[j], c[j],mci, markeredgecolor='none')
#		else:
		if param.markeredge is False:
			ax.plot(x[j], y[j], z[j], mci, markeredgecolor='none',label=leg)
		else:
			ax.plot(x[j], y[j], z[j], mci,label=leg)

#		ax.plot(x[j], y[j], z[j],mci,label=leg)

		ax.plot(x[j], y[j], z[j],mci)
	Options(ax,(x[0],y[0],z[0]),param)
	return (fig,ax)

def closeall(n=50):
	for i in range(n):
		pyplot.close()
		

def flim1(figure,b,mark_col="-k",a_label = 'axis 1',b_label = 'axis 2'):

	if figure is False:
		fig = pyplot.figure( 1 )
		ax = fig.add_subplot( 111 )
		im = ax.imshow( np.zeros( ( 256, 256, 3 ) ) ) # Blank starting image
		figure = (fig,ax,im)
		fig.show()
		im.axes.figure.canvas.draw()
	else:
		fig = figure[0]
		ax = figure[1]
		im = figure[2]
	im.set_data( b )
	im.axes.figure.canvas.draw()
	return (fig,ax,im)

def save(path,i=-1):
	if i == -1:
		s = path
	else:
		s = path + 'temp_0000'
		s2 = str(i)
		s = s[:-len(s2)] + s2
		s = s + '.png'
	pyplot.savefig(s)


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
			zmin = np.min(z)
			zmax = np.max(z)
		else:
			zmin = param.clim[0]
			zmax = param.clim[1]
	else:
		zmin = param.zlim[0]
		zmax = param.zlim[1]

	pylab.ion()
	fig = pylab.figure()
	ax = fig.add_subplot(111)
	if param.clim is not False:
		cmin = param.clim[0]
		cmax = param.clim[1]
	else:
		cmin = zmin
		cmax = zmax

	if param.pcol_type == 'pcolor':
		if param.c_scale == 'log':
			cax = ax.pcolor(x,y,z, cmap=param.colmap, norm=LogNorm(vmin=cmin, vmax=cmax), vmin=zmin, vmax=zmax)
		else:
			cax = ax.pcolor(x,y,z, cmap=param.colmap, vmin=zmin, vmax=zmax)
	else:
		if param.c_scale == 'log':
			cax = ax.contourf(x,y,z, param.ncont, antialiased = True, cmap=param.colmap, norm=LogNorm(vmin=cmin, vmax=cmax), vmin=zmin, vmax=zmax)
		else:
			cax = ax.contourf(x,y,z, param.ncont, antialiased = True, cmap=param.colmap, vmin=zmin, vmax=zmax)
	# set the limits of the plot to the limits of the data

	if param.iscolorbar is False:
#		fig.delaxes(fig.axes[1]) #cbar IS axes[1] while ax is axes[0]
		Options(ax,(x,y),param)
		pylab.show()
		return fig,ax
	else:
		cbar = fig.colorbar(cax,format=param.cbformat)
		Options(ax,(x,y,z),param,cbar)
		pylab.show()
		return fig,ax,cbar

def plot1(data,param=False):
	if param is False:
		param = Paradraw()
	x = np.arange(len(data))
	y = data
	data2 = (x,y)
	fig,ax = plot2(data2,param)
	return fig,ax


def multiplot2bis(data,param=False):
	if param is False:
		param = Paradraw()

	x = data[0]
	y = data[1]

	pylab.ion()
	fig = pylab.figure()
	ax = fig.add_subplot(111)
	ax.plot( x, y, param.mark)
	# set the limits of the plot to the limits of the data
	Options(ax,(x,y),param)

	return fig,ax

def rawplot(data,param,fig,ax):
	x = data[0]
	y = data[1]
	if param is False:
		param = Paradraw()
	if param.markeredge is False:
		ax.plot( x, y, param.mark, markeredgecolor='none',linewidth = param.thickness[0])
	else:
		ax.plot( x, y, param.mark,linewidth = param.thickness[0])

	# set the limits of the plot to the limits of the data
	Options(ax,(x,y),param)

	return fig,ax

def plot2(data,param=False):
	x = data[0]
	y = data[1]
	if param is False:
		param = Paradraw()
	pylab.ion()
	if param.figure is False:
		fig = pylab.figure()
		ax = fig.add_subplot(111)
	else:
		fig = param.figure[0]
		ax = param.figure[1]
	if param.markeredge is False:
		if param.stem is False:ax.plot( x, y, param.marks[0], markeredgecolor='none',linewidth = param.thickness[0],label=param.legend[0])
		else:ax.stem( x, y, param.marks[0], markeredgecolor='none',linewidth = param.thickness[0],label=param.legend[0])
	else:
		if param.stem is False:ax.plot( x, y, param.marks[0],linewidth = param.thickness[0],label=param.legend[0])
		else:ax.stem( x, y, param.marks[0],linewidth = param.thickness[0],label=param.legend[0])
	# set the limits of the plot to the limits of the data
	Options(ax,(x,y),param)

	return fig,ax

def Options(ax,X,param=False, cbar=False):
	x=X[0]
	y=X[1]
	try:
		z=X[2]
		if cbar is False:
			is3d=True
		else:
			is3d = False
	except:
		is3d = False
#might have a bug with colorbar
	if param is False:
		param = Paradraw()
	if param.xlim is False:
		ax.set_xlim(  [ np.min(x),np.max(x) ]  )
		xm = np.min(x)
		xM = np.max(x)
	else:
		ax.set_xlim(  [ param.xlim[0],param.xlim[1] ]  )
		xm = param.xlim[0]
		xM = param.xlim[1]
	if param.ylim is False:
		ym = np.min(y)
		yM = np.max(y)
		ax.set_ylim(  [ ym,yM ]  )
	else:
		ax.set_ylim(  [ param.ylim[0],param.ylim[1] ]  )
		ym = param.ylim[0]
		yM = param.ylim[1]
	if is3d is True:
		if param.zlim is False:
			zm = np.min(z)
			zM = np.max(z)
			ax.set_zlim(  [ zm,zM ]  )
		else:
			ax.set_zlim(  [ param.zlim[0],param.zlim[1] ]  )
			zm = param.zlim[0]
			zM = param.zlim[1]


	if is3d is True:
		ax.set_zscale(param.z_scale)
	ax.set_xlabel(param.x_label,None,None,fontsize=param.fontsize)
	ax.set_ylabel(param.y_label,None,None,fontsize=param.fontsize)
	if is3d is True:
		ax.set_zlabel(param.z_label,None,None,fontsize=param.fontsize)

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
		ax.set_xticks(xlabel)
		ax.set_xticklabels(ax.get_xticks(),None,None,fontsize=param.ticksize)
		ax.xaxis.set_major_formatter(mtick.FormatStrFormatter(param.ax_x_format))

	if param.y_tick_label is None:
		ax.yaxis.set_visible(False)
	else:
		if param.y_tick_label is False:
			if param.y_scale == 'log':
				ylabel = np.exp(np.linspace(np.log(ym),np.log(yM),3))
			elif param.y_scale == 'symlog':
				print ym,yM
				if ym<0:
					print '$$$$$$$$$$$$$$$44'
					ylabel = np.array([ym,0,yM])
					print ylabel
				else:
					ylabel = np.exp(np.linspace(np.log(ym),np.log(yM),3))
			else:
				ylabel = np.linspace(ym,yM,3)
		else:
			ylabel = param.y_tick_label
		ax.set_yticks(ylabel)
		ax.set_yticklabels(ax.get_yticks(),None,None,fontsize=param.ticksize)
		ax.yaxis.set_major_formatter(mtick.FormatStrFormatter(param.ax_y_format))
	ax.axes.set_xscale(param.x_scale,lintreshx = param.xlintresh)
	ax.axes.set_yscale(param.y_scale,lintreshy = param.ylintresh)
	if is3d is True:
		if param.z_tick_label is None:
			ax.zaxis.set_visible(False)
		else:
			if param.z_tick_label is False:
				if param.z_scale == 'log':
					zlabel = np.exp(np.linspace(np.log(np.min(z)),np.log(np.max(z)),3))
				else:
					zlabel = np.linspace(zm,zM,3)
			else:
				zlabel = param.z_tick_label
			ax.set_zticks(zlabel)
			ax.set_zticklabels(ax.get_zticks(),None,None,fontsize=param.ticksize)
			ax.zaxis.set_major_formatter(mtick.FormatStrFormatter(param.ax_z_format))

	ax.set_title(param.title)
	if is3d is True:
		ax.xaxis.pane.fill = False
		ax.yaxis.pane.fill = False
		ax.zaxis.pane.fill = False
	if cbar is not False:
		if is3d is True:
			c=X[3]
		else:
			c=X[2]		
		cbar.set_label(param.c_label,fontsize=param.fontsize)
		if param.cbar_tick_label is False:
			if param.clim is False:
				clabel = np.linspace(np.min(c),np.max(c),3)
			else:
				clabel = np.linspace(param.clim[0],param.clim[1],3)
				
		else:
			clabel = param.cbar_tick_label


#		cbar.set_ticklabels(clabel)
		for t in cbar.ax.get_yticklabels():
			t.set_fontsize(param.ticksize)
		cbar.set_ticks(clabel)
#		cbar.set_ticklabels(clabel)
#		cbar.ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
		cbar.update_ticks()

	if param.y_tick_label is None: ax.yaxis.set_visible(True)
	if is3d is True:
		if param.z_tick_label is None: ax.zaxis.set_visible(True)

	if all(legend is False for legend in param.legend) is False:
		ax.legend(loc=param.legloc)
	if is3d is False:	
		if param.tight_layout is True:
			try:
					pyplot.tight_layout()
			except:
				print('*********************************************************')
				print('not tight')
				print('most probably, figure parameter in Paradraw has to be set to False')
				ax.tight_layout()
#	ax.tick_params(axis='both', which='minor', labelsize=param.ticksize)


def change_options():
	ax = pyplot.gca()
#	ax.tick_params(axis='both', which='minor', labelsize=param.ticksize)
	ax.xaxis.set_minor_locator(pyplot.FixedLocator([2,4]))





def subplot(y,x=False,params=False):
	import copy
	if params is False:
		params = Paradraw()
	n = len(y)
	f, a = pyplot.subplots(nrows=n, ncols=1)
	for i in xrange(n):
		params[i].tight_layout = False
		ax = a[i]
		params[i].figure = (f,ax)
		if x is False:xx = np.arange(len(y[i]))
		else: xx = x[i]
		plot2((xx,y[i]),param = params[i])


def hist(data,n=20,tresh = False,param = False):
	data2 = data+0.0
	if tresh is not False:
		data2[data2>tresh] = tresh
		data2[data2<-tresh] = -tresh
	if param is False:
		print('default parameters') 
		param = Paradraw()
		param.marks = ['k','or']
		param.x_label = 'bins'
		param.y_label = 'Frequency'
		param.ax_x_format = '%.2f'
		param.ax_y_format = '%.0f'
#		if is3d is True:
#			param.ax_z_format = '%.0f'
	hist, bins = np.histogram(data2, bins=n)
	hist = hist
	width = 0.7 * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	pylab.ion()
	fig = pylab.figure()
	ax = fig.add_subplot(111)
	ax.bar(center, hist, align='center', width=width, color = param.bar_col)
	Options(ax,(center,hist),param)
	pylab.show()
	return fig,ax

def plot_Mahalanobis(d,r,param=False):
	if param is False:
		print('default parameters') 
		param = Paradraw()
		param.marks = ['k','or']
		param.x_label = 'Normal theoritical quantiles'
		param.y_label = 'Ordered Mahalanobis distances'
		param.ax_x_format = '%.0f'
		param.ax_y_format = '%.0f'

	x2=np.linspace(  np.min( d[0] ),np.max( d[0] ), len(d[0])  )
	y2 = r[1] + r[0]*x2

	x = ( x2, d[0] )
	y = ( y2, d[1] )
	multiplot2(x,y,param)


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
		name = 'fig_' + datetime.datetime.now().strftime('%Y%m%d-%H%M')

	if fig is not None:
		try:
			pyplot.savefig(path + name + fig)
		except:
			pyplot.savefig(path + name + '.eps')
			pyplot.savefig(path + name + '.pdf')
			pyplot.savefig(path + name + '.png')
	else:
		pyplot.savefig(path + name + '.eps')
		pyplot.savefig(path + name + '.pdf')
		pyplot.savefig(path + name + '.png')

