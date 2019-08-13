import pickle
import numpy as np
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
#	print('np.set_printoptions(precision=prec,suppress=True)')
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
		print 'nothing to save'
		return
	if filename is False:
		from datetime import date
		filename = 'data_' + date.today().isoformat() + '.dat'
	output = open(filename, 'wb')
	pickle.dump(obj,output)
	output.close()
	if verbose is True:print('file saved in:' + filename)
	return filename

def load(filename = False):
	if filename is False:
		print 'nothing to load'
		return
	input_ = open(filename, 'rb')
	obj = pickle.load(input_)
	input_.close()
	return obj
 
