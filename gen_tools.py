import numpy as np # for computations
import uuid # for random names
from random import randint, random, choice, gauss
from treelib import Tree, Node #pip treelib
from copy import deepcopy

from div_tools import TwoWayDict
from joblib import Parallel, delayed #pip joblib

# lists
# modifying the lists is not enough: RPN and the Rand_xxx functions have to be updated as well.
list_function = ['cos','sin','tanh','exp']
list_operator = [ '+', '*', '/']
list_reductor = ['sum','x','y','z']
list_possibilities = ['number','operator','function','reductor']

def Rand_type(forbiden = [], chances = None):
        '''
        Function that select a random type for a node: a function, a number, an operator or a reductor
        '''
        if chances is None:
                chances = [20,50,80]
        r = random()*100
        if r<chances[0]:
                type_ = list_possibilities[0]
        elif r<chances[1]:
                type_ = list_possibilities[1]
        elif r<chances[2]:
                type_ = list_possibilities[2]
        else: 
                type_ = list_possibilities[3]
        while type_ in forbiden:
                type_ = Rand_type(forbiden)
        return type_

def Rand_function(forbiden = [], chances = None):
        '''
        Function that select a random function for a node
        '''
        if chances is None:
                chances = [20,70,98]
        r = random()*100
        if r<chances[0]:
                func = list_function[0]
        elif r<chances[1]:
                func = list_function[1]
        elif r<chances[2]:
                func = list_function[2]
        else: 
                func = list_function[3]
        while func in forbiden:
                func = Rand_function(forbiden)
        return func

def Rand_depth(ind):
        keys = [node.identifier for node in ind.tree.all_nodes()]
        d=ind.tree.DEPTH
        s = 0.85*d/2
        l = int(gauss((d-1)/2.,s))-1
        if l>d: l=d
        if l<1: l=1
        id_ = choice(keys)
        while np.abs(ind.tree.level(id_) - l)>2:
                id_ = choice(keys)
        return id_


def Rand_reductor(forbiden = [], chances = None):
        '''
        Function that select a random reductor for a node
        '''
        if chances is None:
                chances = [0,100,100,100]
        r = random()*100
        if r<chances[0]:
                red = list_reductor[0]
        elif r<chances[1]:
                red = list_reductor[1]
        elif r<chances[2]:
                red = list_reductor[2]
        else: 
                red = list_reductor[3]
        while red in forbiden:
                red = Rand_reductor(forbiden)
        return red

def Rand_operator(forbiden = [], chances = None):
        '''
        Function that select a random operator for a node
        '''
        if chances is None:     
                chances = [50,100]
        r = random()*100
        if r<chances[0]:
                ope = list_operator[0]
        elif r<chances[1]:
                ope = list_operator[1]
        else: 
                ope = list_operator[2]
        while ope in forbiden:
                ope = Rand_operator(forbiden)
        return ope

                
def Population(number=100,depth = 8):
        '''
        Create a population of Individual, i.e., random functions.
        '''
        pop = []
        for i in range(number):
                pop.append(Individual(depth = depth))
        return pop


class Node_value:
        '''
        Node of the tree describing the function. It can be a scalar (a component of x or a random number), an operator (i.e., + ) or a function (i.e., cos)
        *Mutate() changes the value of the node.
        '''
        def __init__(self,_type = None,_value = None):
                if _type is not None:
                        self._type = _type
                else:
                        _type = randint(0,3)
                        self._type = list_possibilities[_type]
                
                if self._type == 'number': self.n_child = 0
                if self._type == 'operator': self.n_child = 2
                if self._type == 'function': self.n_child = 1
                if self._type == 'reductor': self.n_child = 0

                if _value is not None:
                        self._value = _value
                else:
                        self.Mutate()

                
        def __str__(self):
                return str(self._value)
                
        def Mutate(self):
                if self._type == 'number': self._value = gauss(0,3.5)
                if self._type == 'operator': self._value = Rand_operator()
                if self._type == 'function': self._value = Rand_function()
                if self._type == 'reductor': self._value = Rand_reductor()


class Individual():
        '''
        It represents a function f. It is basically a tree. Each node represents a step for describing the function in the Reverse Polish Notations.
        *Copy(indiv) change the tree or return a copy of the tree
        *Mutate(id_) mutates a random or given node.
        *CrossOver(father) does the random breeding of the individual with an other one.
        *Evaluate(x) gives back f(x). If x=None, evaluate f(0).
        *show gives a visual representation of the tree
        *depth gives back the depth of the function
        *rpn gives back the reverse polish expression of the function
        
        '''
        def __init__(self,depth=8):
                self.max_depth = depth
                self.tree = Small_tree()
                for i_depth in range(1,depth):
                        if i_depth == depth-1:
                                isleaves = True
                        else:
                                isleaves = False
                                
                        leaves = [n.identifier for n in self.tree.leaves()]
                        for leave in leaves:
                                if self.tree.level(leave) == i_depth:
                                        st = Small_tree(id_ = leave, isleaves = isleaves, node = self.tree[leave].tag)
                                        parent_id = self.tree[leave].bpointer
                                        self.tree.remove_node(leave)
                                        self.tree.paste(nid = parent_id,new_tree = st)
                self.rpn = T2E(self.tree)
        def Copy(self,indiv=None):
                if indiv is None:       return deepcopy(self)
                else:   
                        self = deepcopy(indiv)
                        self.rpn = T2E(self.tree)


        def Mutate(self, id_ = None):
                ''' Mutate a node. if id_ is None, the node is randomly chosen ''' 
                if id_ is None:
                        id_ = choice([n.identifier for n in self.tree.all_nodes()])
                else:
                        if self.tree.contains(id_) is False:
                                        print('id invalid: ' + str(id_)+ ' is not a valid node')
                                        return
                self.tree[id_].tag.Mutate()
                self.rpn = T2E(self.tree)

                
        def CrossOver(self,father = None):
                ''' Crossover between the present indivudual - i.e., the mother- and a father.
                A random branch from the father is pasted to a random branch of the mother
                 ''' 
                #self.show
                if father is None:
                        print('need a father for the breeding')
                        return

                #random branch of the mother
                mother_id = choice([n.identifier for n in self.tree.all_nodes()])
                depth_mother = self.tree.level(mother_id)

                #random branch of the father
                father_id = choice([n.identifier for n in father.tree.all_nodes()])
                depth_father = father.tree.subtree(father_id).depth()

                #The total depth max_depth should not be exceeded
                nit = 0
                nitmax = 10
                while (self.max_depth < depth_mother + depth_father) and (nit < nitmax) :
                        nit +=1
                        father_id = choice([n.identifier for n in father.tree.all_nodes()])
                        depth_father = father.tree.subtree(father_id).depth()
                if nit >= nitmax: #breeding failed
                        return

                #The node before has to be identified
                parentid = self.tree[mother_id].bpointer
                if parentid is None:
                        # the mother_id is actually the root of the tree
                        self.tree = ReId(father.tree.subtree(father_id))
                else:           
                        # the branch is removed and the father branch is pasted. It is renamed for avoiding collusions
                        self.tree.remove_node(mother_id)
                        self.tree.paste(nid = parentid, new_tree = ReId(father.tree.subtree(father_id)))

                self.rpn = T2E(self.tree)
        
        def Evaluate(self,input_=0):
                return RPN(self.rpn,input_)
        @property
        def show(self):
                self.tree.show()
        @property
        def depth(self):
                return self.tree.depth()
        @property
        def expression(self):
                return self.rpn
                                
def ReId(tree):
        ''' create a copy of the tree with new unique ID for each nodes, in order to avoid collisions of node IDs.
        '''     
        new_tree = Tree()
        dic = TwoWayDict()
        for node in tree.all_nodes():
                dic[node.identifier] = str(uuid.uuid4())
                
        for depth in range(tree.depth()+1):    
                if depth == 0:
                        root_id = dic[tree[tree.root].identifier]
                        root_node = tree[tree.root].tag
                        new_tree.create_node(root_node, root_id, parent = None)
                else:
                        nodes = []
                        for node in tree.all_nodes():
                                if tree.level(node.identifier) == depth: nodes.append(node)
                        for node in nodes:
                                node_id = dic[node.identifier]
                                leaf_node = node.tag
                                new_tree.create_node(leaf_node, node_id, parent = dic[node.bpointer])

        return new_tree

def Small_tree(id_ = None, isleaves = False, node = None):
        ''' create a small tree with maximum two leaves, that will be added on a node, in order to create a random function.
        -> "isleaves" states this is the last possible depth. only reduction and number can be accepted
        -> "node" is the value of the root. None means a random one will be created
        '''
        #initialization
        st = Tree()
        if id_ is None:
                root_id = str(uuid.uuid4())
        else:
                root_id = id_
        if node is None:
                root_node = Node_value(Rand_type())     
        else:
                root_node = node                
        st.create_node(root_node, root_id, parent = None)
        for i_child in range(root_node.n_child):
                leave_id = str(uuid.uuid4())
                if isleaves is True:
                        leave_type = Rand_type(forbiden = ['function','operator'])
                else:
                        leave_type = Rand_type()
                st.create_node(Node_value(leave_type), leave_id, parent = root_id)
        return st


        
def T2E(tree):
        return Tree2expr(tree)

def Tree2expr(tree):
        """ Express the tree as a reverse polish expression. """
        exp =  [tree[node].tag._value for node in tree.expand_tree()]
        return exp[::-1]

def Minimax(n,tresh = 100):
        return np.min([tresh,np.max([-tresh,n])])
                
def RPN(expression,input_):
        """ Evaluate a reverse polish notation """
        stack = []
        try:
                for val in expression:
                        #basic operation
                        if val in list_operator:
                                op1 = stack.pop()
                                op2 = stack.pop()
                                if val=='-': result = op2 - op1
                                if val=='+': result = op2 + op1
                                if val=='*': result = op2 * op1
                                if val=='/': result = op2 / op1
                                stack.append(Minimax(result))
                        #basic functions
                        elif val in list_function:
                                op1 = stack.pop()
        #                       print(op1)
                                if val=='cos': result = np.cos(op1)
                                if val=='sin': result = np.sin(op1)
                                if val=='tanh': result = np.tanh(op1)
                                if val=='exp': result = np.exp(Minimax(op1)) # avoiding overflow
                                stack.append(result)
                        #basic reductions
                        elif val in list_reductor:
#                               op1 = stack.pop()
                                op1 = input_
                                if val=='sum': result = np.sum(op1)
                                if val=='x': 
                                        try:
                                                result = op1[0]
                                        except:
                                                result = op1
                                if val=='y': 
                                        try:
                                                result = op1[1]
                                        except:
                                                result = op1
                                if val=='z':
                                        try:
                                                result = op1[2]
                                        except:
                                                result = op1
                                stack.append(result)
                        else:
                                stack.append(val)
        except :
                print('error, expression not valid')
                print("Unexpected error:", sys.exc_info()[0])
                stack = [np.nan]
                
        return stack.pop()

def Evolve(pop, f_fit=None, f_mut=None, f_cross=None, f_new = False, retain=0.35, random_select=0.15, mutate=0.1, new = 0.01, verbose = False,  npar = 0,parameters = False,is_final = False):
        """ Evolution of a population. 
        * f_fit: fitness function, ie the function that computes the score -cost- associated to an individu.
        * f_mut: mutation function, ie the function that describes how mutation of an individu is done.
        * f_cross: crossover function, ie the function that describes how breeding betwin two individus is done.
        * retain: fraction of best individuals of the population which is kept.
        * random_select: chances for an individual to be kept, to promote diversity.
        * mutate: chances for an individual to mutate.
        """

        if (f_fit == None) or(f_cross == None) or(f_mut == None):
                print('functions needed')
                print('f_fit: fitness function, ie the function that computes the score -reward- associated to an individu. It should take an individual as argument and return a score as a scalar - the lower, the better.')
                print('f_mut: mutation function, ie the function that describes how mutation of an individu is done. It should take an individual as argument, and return a mutated individual.')
                print('f_cross: crossover function, ie the function that describes how breeding betwin two individus is done. It should take two individuals - parents - as arguments, and returned a individual -- the child -  as the crossover of the two parents.')
                return pop 
        #compute scores 
        if npar == 0:
                scores = [ (f_fit(individual), individual) for individual in pop]
        else:
                scores_temp = Parallel(n_jobs=npar)(delayed(f_fit)(individual) for individual in pop)
                scores = [ (scores_temp[i],pop[i]) for i in range(len(pop))]
#
        scores = [ individual for individual in sorted(scores, key=lambda indiv: indiv[0])]
        if verbose is True:
                s = scores[0][0]
                score_median = np.median([score[0] for score in scores]) 
                score_std = np.std([score[0] for score in scores]) 
                print('score of the best individual :' + str(s))
                print('median :' + str(score_median) + ' and std:' + str(score_std))

            
        # keep the best
        retain_length = int(len(scores)*retain)
        if retain_length<3:retain_length=3
        parents = [scores[i_keep][1] for i_keep in range(retain_length)]

        # randomly add other individuals to
        # promote genetic diversity
        for individual in scores[retain_length:]:
                if random_select > random():
                        parents.append(individual[1])
                if f_new is not False:
                        if new > random():
                                parents.append(f_new())

        # mutate some individuals
        for individual in parents:
                if mutate > random():
                        individual = f_mut(individual)


        # crossover parents to create children
        parents_length = len(parents)
        desired_length = len(pop) - parents_length
        children = []
        while len(children) < desired_length:
                male = np.random.randint(0, parents_length-1)
                female = np.random.randint(0, parents_length-2)
                if male != female:
                        child = f_cross(parents,male,female)
                        children.append(child)
                else:
                        female = male+1
                        child = f_cross(parents,male,female)
                        children.append(child)

        parents.extend(children)
        if is_final is True:
            return scores[0][1],parents,[score[0] for score in scores] 
        else:
            return parents,[score[0] for score in scores] 
