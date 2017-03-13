
from twisted.internet import reactor

from viff.field import FieldElement
from viff.runtime import Runtime, Share, gather_shares

import random


class LocalRuntime(Runtime):
    def __init__(self, player, threshold, options=None):
        Runtime.__init__(self, player, threshold, options)
        self.res = options.res
        
    def field_and_shares(self, share_a, share_b):
        field = getattr(share_a, "field", getattr(share_b, "field", None))
        if not isinstance(share_a, Share):
            if isinstance(share_a, int): share_a = field(share_a)
            share_a = Share(self, field, share_a)
        if not isinstance(share_b, Share):
            if isinstance(share_b, int): share_b = field(share_b)
            share_b = Share(self, field, share_b)    
                
        return field, share_a, share_b
        
    def prss_share_random(self, Zp):
        return Share(self, Zp, Zp(random.randrange(0, Zp.modulus)))
    
    def truncate(self, val):
        if isinstance(val, FieldElement) and val.fp:
            if (val.r >= self.res):
                return val.field(int(round(float(val.signed()) / (2 ** (val.r - self.res)))), self.res)
            else:
                print "*** Possible loss of precision"
                return val.field(int(round(float(val.signed()) * (2 ** (self.res - val.r)))), self.res)
        else:
            return val
        
    def clone_share_nofp(self, share):
            """ clone a Share with fp-flag disabled."""

            clone = share.clone_nofp()
            clone.addCallback(lambda x: share.field(x.value))
            return clone                
        
    def clone_share_fp(self, share):
            """ clone a Share with fp-flag enabled."""

            clone = share.clone_fp()

            clone.addCallback(lambda x: share.field(x.value, True))
            return clone                
        

    def truncPR(self, sh, f, in_fp=False):
        def doTrunc(val, f, fp):
            sv = float(val.signed()) / float(2 ** f)
            return val.field(int(round(sv)), fp)

        field = sh.field
        l = self.options.bit_length
        k = self.options.security_parameter

        # ensure normal addition
        fp = in_fp or sh.fp
        if isinstance(sh, Share):
                share = self.clone_share_nofp(sh)
        else:
                share = Share(self, field, sh.field(sh.value))
                 
        share.addCallback(doTrunc, f, fp)

        if fp: share = self.clone_share_fp(share)
        return share     
    

    def calc2(self, share_a, share_b, fn):
        field, share_a, share_b = self.field_and_shares(share_a, share_b)
        
        result = gather_shares([share_a, share_b])
        result.addCallback(fn)
            
        return result
    
                
    def greater_than_equal(self, share_a, share_b):
        return self.calc2(share_a, share_b, lambda (x, y): x - x + (1 if x >= y else 0))
 
    def sub(self, share_a, share_b):
        return self.calc2(share_a, share_b, lambda (x, y): x - y)

    def add(self, share_a, share_b):
        return self.calc2(share_a, share_b, lambda (x, y): x + y)
    
    def mul(self, share_a, share_b, dotrunc=True):        
        trunc = False
                
        if not isinstance(share_b, Share):
            # Local multiplication. share_a always is a Share by
            # operator overloading in Share. We clone share_a first
            # to avoid changing it.
            result = share_a.clone()
            result.addCallback(lambda a: share_b * a)
            return result
        elif share_a.fp and share_b.fp and dotrunc:
            trunc = True
                
        result = gather_shares([share_a, share_b])
        result.addCallback(lambda (a, b): a * b)
        if trunc: self.schedule_callback(result, self.truncPR, self.res, True)

        return result                        
    
    def div(self, share_a, share_b):
        def divfn((x, y)):
            # retv=int(round(float(self.get_val(x))/float(self.get_val(y))/(2**(self.res))))
            if isinstance(x, FieldElement):
#                print float(self.get_val(x))
#                print float(self.get_val(y))
                return x.field(int(round(float(self.get_val(x)) / float(self.get_val(y)) * (2 ** (self.res)))), self.res)
            else:
                return int(round(float(self.get_val(x)) / float(self.get_val(y))))
            
#            if isinstance(x,FieldElement) and x.fp and isinstance(y,FieldElement) and y.fp:
#                iv=int(round(float(x.signed())/float(y.signed())*(2**y.r)))
#                return x.field(iv,x.r)
#            else:
#                iv=int(round(float(x.signed())/float(y.signed())*(2**y.r)))
#                return x.field(iv,self.res)
            
        return self.calc2(share_a, share_b, divfn)
    
    def pow(self, share_a, exp):
        result = share_a.clone()
        result.addCallback(lambda a: a ** exp)
        return result     
    
    def gtefun(self, x, y):
        if x > x.field(2 ** self.options.bit_length):
            if y > x.field(2 ** self.options.bit_length):
                ret = 1 if x >= y else 0    # both are negative
            else:
                ret = 0    # y is positive
        else:
            if y > x.field(2 ** self.options.bit_length):
                ret = 1    # y is negative
            else:
                ret = 1 if x >= y else 0    # both are positive
        return ret
    
    def greater_than_equal(self, share_a, share_b):
        return self.calc2(share_a, share_b, lambda (x, y): self.gtefun(x, y))
    

    def equal(self, share_a, share_b):
        return self.calc2(share_a, share_b, lambda (x, y): x.field(1) if x == y else x.field(0))    
    
    # from passive.py
    def lin_comb(self, coefficients, shares):
        """Linear combination of shares.

        Communication cost: none. Saves the construction of unnecessary shares
        compared to using add() and mul()."""

        assert len(coefficients) == len(shares), \
                "Number of coefficients and shares should be equal."

        assert all(map(lambda s: not isinstance(s, Share), coefficients)), \
                        "Coefficients should not be shares."

        assert all(map(lambda s: isinstance(s, Share), shares)), \
                        "Shares should be shares."

        def computation(shares, coefficients):
            return sum(map(operator.mul, coefficients, shares))

        result = gather_shares(shares)
        result.addCallback(computation, coefficients)
        return result

    # from passive.py    
    def bin_comb(self, shares):
            """Binary combination of shares."""
            def computation(shares):
                    sum = 0
                    for i in xrange(len(shares) - 1, -1, -1):
                            sum = 2 * sum + shares[i]        
                    return sum
            result = gather_shares(shares)
            result.addCallback(computation)
            return result    
    
    def in_prod(self, shares_a, shares_b):
            def computation(share_list):
                    n = len(share_list) / 2
                    s = 0
                    for i in xrange(n):
                            s += share_list[i] * share_list[n + i] 
                    return s

            result = gather_shares(shares_a + shares_b)
            result.addCallback(computation)
            result.addCallback(self.truncate)
            return result    
        
    def neg(self, share_a):
        result = share_a.clone()
        result.addCallback(lambda a:-a)
        return result        
    
    def invert(self, share_a):
        result = share_a.clone()
        result.addCallback(lambda a:~a)
        return result        
    
    
    def gauss(self, A, d, b, c):
        def computation(x, m):
                n = len(x) / (m + 1) - 1
                mn = m * n
                for i in xrange(m):
                        for j in xrange(n):
                                x[n * i + j] = x[n * i + j] * x[mn] - x[mn + 1 + i] * x[mn + 1 + m + j]
                del x[mn:]
                return x
        
        def peek(s, k):
                return s[k]
                
        def wrap(s, m, n):
                return [[s.clone().addCallback(peek, n * i + j) for j in xrange(n)] for i in xrange(m)]
                 
        result = gather_shares(reduce(operator.add, A) + [d] + b + c)
        self.schedule_callback(result, computation, len(A))
                
        return wrap(result, len(A), len(A[0]))
            
    def scalar_mul(self, a, x):
        def computation(x):
            a = x.pop()
            return [self.truncate(xi * a) for xi in x]
        
        result = gather_shares(x + [a])
        result.addCallback(computation)
        
        def peek(s, k): return s[k]        
        def wrap(s, n): return [s.clone().addCallback(peek, i) for i in xrange(n)]        
        return wrap(result, len(x))

    def open(self, share):
        return share    
        
    def bit_decompose(self, val, l):
        def peek(s, k):
            return s.field((s.value >> k) & 1)
                        
        def wrap(s):
                return [s.clone().addCallback(peek, i) for i in xrange(l)]         
            
        return wrap(val)
    
    def prefix_orC(self, bits):
        ret = [bits[0]]
        for i in range(len(bits)):
            ret = ret + [ret[-1] + bits[i] - ret[-1] * bits[i]]
        return ret

    def get_val(self, val): 
#                print val.fp 
#                return val.signed()
#        print "val", type(val)
#        print "FP?", val.fp
        if not isinstance(val, FieldElement): return val
            
        if val.fp:
            return float(val.signed()) / (2 ** val.r)
        else:
            return val.signed()    
        
    def get_value(self, val): 
#                print val.fp        
#                return val.signed()
#        print "val", type(val)
#        print "FP?", val.fp
        if val.fp:
            return float(val.signed()) / (2 ** val.r)
        else:
            return val.signed()                
        
    def shutdown(self):
        reactor.stop()

    def shamir_share(self, inputters, field, number=None, threshold=None, defaults=None):
        def val(i):
            if i == self.id:
                assert number != None, "No value given for own share"
                return field(number, fp)
            else:
                assert defaults != None and defaults[i] != None, "No default value given for share of other player " + str(i) + " " + str(self.id)
                return defaults[i]    # should be field element
            
        fp = isinstance(number, float)
        if fp:
            number = int(2 ** self.res * number)
        
        results = [Share(self, field, val(i), fp) for i in inputters]

        if len(results) == 1:
            return results[0]
        else:
            return results        


def localversion(localf):
    """ Decorator for a function in Shamir shares that has a LocalRuntime
        alternative localf. """
    def localversion_impl(realf):
        def runner(*args, **kwargs):
            return localf(*args, **kwargs) if isinstance(args[0], LocalRuntime) else realf(*args, **kwargs)
        return runner
    return localversion_impl 