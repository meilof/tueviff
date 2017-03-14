from viff.runtime import Share, gather_shares
from twisted.python.failure import Failure
from twisted.internet.defer import Deferred, _inlineCallbacks, returnValue

def all_shares(obj): # could be replaced by standard Python library call?
    if type(obj) is list:
        return reduce(lambda x,y: y+x, [all_shares(x) for x in obj])
    return [obj]
    
def mkstruct(ret,x):
    if type(ret) is list:
        return [mkstruct(reti, x) for reti in ret]
    return x.pop()


class pc_wrapper(object):
    """ Decorator inside which the program counter is forked. """
    def __init__(self, runtime, f):
        self.runtime = runtime
        self.runtime.increment_pc()
        self.runtime.fork_pc()
        self.saved_pc = self.runtime.program_counter[:]
        self.runtime.unfork_pc()
        self.f = f
        
    def send(self, value):
        current_pc = self.runtime.program_counter[:]
        self.runtime.program_counter[:] = self.saved_pc
        try:
            reto = self.f.send(value)
            
            # Recursively collect all shares in the yielded value. This enables
            # to yield on (nested) lists without explicitly having to call
            # gather_shares. Moreover, note that directly returning "ret" does
            # _not_ work since Twisted's inlinecallback mechanism consumes the
            # returned default which may have side-effects in the calling code,
            # see https://github.com/twisted/twisted/pull/707.
            ret = gather_shares(all_shares(reto)) 

            # Because of gather_shares, ret becomes a list of shares: here we
            # transform it back to the original shape (single Share, list, or
            # nested list).
            ret.addCallback(lambda x: mkstruct(reto, x))
            return ret
        finally:
            self.saved_pc = self.runtime.program_counter[:]
            self.runtime.program_counter[:] = current_pc
        
def reconcile(decl,givn):
    if isinstance(givn,Failure): # detect failures from the inline callback
       givn.raiseException()
    elif decl is None:
        return
    elif type(decl) is list:
        for (d,g) in zip(decl,givn): reconcile(d,g)
    elif isinstance(givn, Deferred):
        givn.chainDeferred(decl)
    else:
        decl.callback(givn)
        
def deepcopy(val):
    return [deepcopy(x) for x in val] if isinstance(val,list) else val
    
def declareReturn(rt, Zp, *args):
    def make(ix):
        if ix>=len(args): return Share(rt, Zp)
        return [make(ix+1) for _ in xrange(args[ix])]
    return rt, make(0)
    
def declareReturnNop(rt, Zp):
    return rt, None
            
def viffinlinecb(f):
    """
    VIFF inline callback mechanism.
    
    The decorated function should be a generator that:
    
     - yields as first value (rt,shs), where rt is the used runtime and shs
       is a Share/(nested) list of Shares that will be the deferred result
       of this function. For convenience, (rt,shs) can be generated using
       declareReturn or declareReturnNop: use declareReturn(rt,Zp) for
       functions that return a single value; declareReturn(rt,Zp,n1,n2,...) for
       functions returning n1xn2x...-shaped nested list, and declareReturnNop
       for functions without a return value.
       
     - next, yields any Shares or (nested) lists of shares that it needs to
       have the values of. The return of the yield is the value of that share.
       
     - finally, returns values or Shares with returnValue.    
    """
    
    def unwindGenerator(*args, **kwargs):
        gen = f(*args, **kwargs)
        rt, ret = gen.send(None)
        defr = _inlineCallbacks(None, pc_wrapper(rt, gen), Deferred())
        defr.addCallback(lambda x: reconcile(ret, x))
        return deepcopy(ret)
        
    return unwindGenerator
