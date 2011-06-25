from __future__ import division, with_statement
from time import time
import itertools as it,operator as op, numpy as np 
import re,cPickle, os
from os.path import join

##########################
####### Decorators #######
##########################



def simple_decorator(decorator):
    """This decorator can be used to turn simple functions
    into well-behaved decorators, so long as the decorators
    are fairly simple. If a decorator expects a function and
    returns a function (no descriptors), and if it doesn't
    modify function attributes or docstring, then it is
    eligible to use this. Simply apply @simple_decorator to
    your decorator and it will automatically preserve the
    docstring and function attributes of functions to which
    it is applied."""
    def new_decorator(f):
        g = decorator(f)
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        g.__dict__.update(f.__dict__)
        return g
    # Now a few lines needed to make simple_decorator itself
    # be a well-behaved decorator.
    new_decorator.__name__ = decorator.__name__
    new_decorator.__doc__ = decorator.__doc__
    new_decorator.__dict__.update(decorator.__dict__)
    return new_decorator



@simple_decorator
def time_fun(func):
    def _timed_ver(*args,**kw):
        __name__ = func.__name__
        t_start = time()
        print('%s: '%(func.__name__))
        res = func(*args, **kw)        
        print('done: %f.1 seconds'%(t_start-time()))
        return res 
    return _timed_ver
        
def first_time(code,globals_,locals_):
    """use like this:
    @first_time('foo=3',globals(),locals()
    def bar():
        pass
    The first time bar is called, the code will be executed in the current scope."""
    
    class decorator(object):
        def __init__(self,func):
            self.func = func
            self.dispatch = self._first_call
        def _first_call(self,*args,**kw):
            exec code in globals_,locals_
            self.dispatch = self._rest_call
            return self.func(*args,**kw)
        def _rest_call(self,*args,**kw):
            return self.func(*args,**kw)
        def __call__(self,*args,**kw):
            return self.dispatch(*args,**kw)
            
    return decorator

class memoized(object):
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        try:
            return self.cache[args]
        except KeyError:
            self.cache[args] = value = self.func(*args)
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)
    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__


#####################################
######### functional programming ####
#####################################

def List(n,fill=None):
    return [fill for _ in xrange(n)]

def mymap(func,seqs,start=0,end=None,func_args={}):
    n = len(seqs[0])
    Out = List(n)
    SeqArgs = zip(*seqs)    
    for i in xrange(start,n or end):
        Out[i] = func(*SeqArgs[i],**func_args)
    return Out

class _fakeiterator():
    def __iter__(self):
        return self
    def next(self):
        return self.val
    def send(self,val):
        self.val = val
        
class consumerize(object):
    def __init__(self,gen_func,*IterTargets,**kw):
        self.its = []
        for targ_str in IterTargets:
            kw[targ_str] = _fakeiterator()
            self.its.append(kw[targ_str])
        self.gen = gen_func(**kw)
        
    def send(self,*vals):
        for val,iterator in zip(vals,self.its):
            iterator.send(val)
        return self.gen.next() 
    def close(self):
        self.gen.close()
            
def helicase(Iterator,n_values):
    ValList = list(it.islice(Iterator,n_values))
    return ValList,it.chain(ValList,Iterator)


def chain_from_iterable(iterables):
    for seq in iterables:
        for element in seq:
            yield element             

def identity(x):
    return x
            
def mapdict(Dict,keymapper=None,valmapper=None):
    keymapper = keymapper or identity
    valmapper = valmapper or identity
    return dict([(keymapper(k),valmapper(v)) for k,v in Dict.items()])
    
def append(lol):
    return reduce(op.__add__,lol)

########################
####### Arrays/math ####
########################

def flatdist(N):
    return np.ones(N,dtype=np.float32)/N

def max2min(li):
    return max(li)-min(li)

def tofloat32(Arr):
    return Arr.astype(np.float32)

def to2d(Arr):
    return Arr.reshape(Arr.shape[0],-1)

def inrange(x,lo,hi):
    return x >= lo and x < hi

def unraveled_argmax(X_hm):
    return np.unravel_index(X_hm.argmax(),X_hm.shape)

def unraveled_argmin(X_hm):
    return np.unravel_index(X_hm.argmin(),X_hm.shape)

def get_padded(Arr,Start,End):
    if Start < 0:
        StartZeros = np.zeros((-Start,Arr.shape[1]),dtype=Arr.dtype)
        return np.vstack((StartZeros,Arr[:End]))
    elif End > Arr.shape[0]:
        EndZeros = np.zeros((End-Arr.shape[0],Arr.shape[1]),dtype=Arr.dtype)
        return np.vstack((Arr[Start:],EndZeros))
    else:
        return Arr[Start:End]

def naive_maximize(func,domain,func_args={},MinimizeInstead=False):
    values = mymap(func,seqs=(domain,),func_args=func_args)    
    return domain[np.argmin(values) if MinimizeInstead else np.argmax(values)]

####################
#### Li and Stri ###
####################


def first(li): return li[0]
def second(li): return li[1]
def is_numerical(x):
    try:
        float(x)
        return True
    except Exception:
        return False

def is_bool(x):
    return x is True or x is False

def dict_append(old_dict,new_dict):
    old_dict.update(new_dict)
    return old_dict

def avg(seq):
    iterator = iter(seq)
    count = 1
    total = iterator.next()
    for el in iterator:
        count += 1
        total += el
    return total/count

def index_safe(el,li):
    try: return li.index(el)
    except ValueError: return None
        
########################
###### Files ###########
########################

def dump(filename,obj):
    with open(filename,"w") as fd: cPickle.dump(obj,fd)
    
    
def load(filename):
    with open(filename,"r") as fd: return cPickle.load(fd)

class indir(object):
    def __init__(self,new_dir):
        self.orig_dir = os.getcwd()
        self.new_dir = new_dir
    def __enter__(self):
        print("entering %s"%self.new_dir)
        mkdir_and_enter(self.new_dir)
    def __exit__(self,*exc_info):
        print("exiting %s"%self.new_dir)        
        os.chdir(self.orig_dir)

def mkdir_and_enter(DirName):
    if not os.path.exists(DirName):
        os.mkdir(DirName)
        os.chdir(DirName)
    else:
        os.chdir(DirName)
        # is_repeat = re.match("(.+_)(\d+)$",DirName)
        # if is_repeat:
        #     DirName = is_repeat.group(1)+str(int(is_repeat.group(2))+1)
        # else:
        #     DirName = DirName+"_1"
        # DirName = mkdir_and_enter(DirName)
    return DirName

def mkdir_maybe_rm(DirName):
    if os.path.exists(DirName):
        os.system("rm -rf %s"%DirName)  
    os.mkdir(DirName)



def find_file_with_ext(directory,ext,ex_if_not_found = False):
    matches = [fname for fname in os.listdir(directory) if fname.endswith(ext)]
    if len(matches) > 1:
        print("multiple files found in %s ending with %s:"%(directory,ext))
        for (i,filename) in enumerate(matches):
            print("\t%i. %s"%(i,filename))
        choice=input("which one do you want?")
        return join(directory,matches[choice])
    elif len(matches) == 1:
        return join(directory,matches[0])
    else:
        if ex_if_not_found:
            raise Exception("Tried to find file with extension %s in directory %s: not found."%(ext,directory))
        else:
            return None

def splitext(filename):
    m = re.match(r"(.+)\.(\w+\.\d+)",filename)
    return m.groups() if m is not None else os.path.splitext(filename)
    
    
        
def switch_ext(filepath,new_ext):
    return splitext(filepath)[0]+"."+new_ext

def basename_noext(filepath):
    return splitext(os.path.basename(filepath))[0]

def parent_dir(filepath):
    return os.path.dirname(filepath)

OUTPUT_FILE = None
def oprint(string):
    if OUTPUT_FILE is None:
        print(string)
    else:
        print(string)
        OUTPUT_FILE.write(string+"\n")               
        
def test_mymap():
    def add(x,y,negate):
        return -(x+y) if negate else x+y
    print mymap(add,(range(0,5),range(5,10)),start=1,func_args=dict(negate=True))
    
def test_maximize():
    fn = lambda x: 1-x**2
    print naive_maximize(fn,range(-5,5))
    
if __name__ == '__main__':
    print os.listdir(".")
    with indir("/home/joschu/Documents"):
        print os.listdir(".")
    print os.listdir(".")
   
    
    
    
