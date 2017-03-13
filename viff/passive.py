# -*- coding: utf-8 -*-
#
# Copyright 2008 VIFF Development Team.
#
# This file is part of VIFF, the Virtual Ideal Functionality Framework.
#
# VIFF is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License (LGPL) as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# VIFF is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General
# Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with VIFF. If not, see <http://www.gnu.org/licenses/>.

"""FP Passively secure VIFF runtime."""

import operator, math, copy

from viff import shamir
from viff.runtime import Runtime, Share, FPShare, ShareList, gather_shares, preprocess
from viff.prss import prss, prss_double, prss_lsb, prss_zero, prss_multi
from viff.field import GF, GF256, FieldElement
from viff.util import rand, profile

import viff.inlinecb

from twisted.internet.defer import gatherResults

from gmpy import mpz


class PassiveRuntime(Runtime):
    """The VIFF runtime.

    The runtime is used for sharing values (:meth:`shamir_share` or
    :meth:`prss_share`) into :class:`~viff.runtime.Share` object and
    opening such shares (:meth:`open`) again. Calculations on shares
    is normally done through overloaded arithmetic operations, but it
    is also possible to call :meth:`add`, :meth:`mul`, etc. directly
    if one prefers.

    Each player in the protocol uses a :class:`~viff.runtime.Runtime`
    object. To create an instance and connect it correctly with the
    other players, please use the :func:`~viff.runtime.create_runtime`
    function instead of instantiating a :class:`~viff.runtime.Runtime`
    directly. The :func:`~viff.runtime.create_runtime` function will
    take care of setting up network connections and return a
    :class:`Deferred` which triggers with the
    :class:`~viff.runtime.Runtime` object when it is ready.
    """

    def __init__(self, player, threshold, options=None):
        """Initialize runtime."""
 
        assert hasattr(options, "res"), "no resolution for FP runtime defined"

        Runtime.__init__(self, player, threshold, options)
        self.res = options.res
		
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
		
    def output(self, share, receivers=None, threshold=None):
        return self.open(share, receivers, threshold)

    def open(self, share, receivers=None, threshold=None):
        """Open a secret sharing.

        The *receivers* are the players that will eventually obtain
        the opened result. The default is to let everybody know the
        result. By default the :attr:`threshold` + 1 shares are
        reconstructed, but *threshold* can be used to override this.

        Communication cost: every player sends one share to each
        receiving player.
        """
        assert isinstance(share, Share)
        # all players receive result by default
        if receivers is None:
            receivers = self.players.keys()
        if threshold is None:
            threshold = self.threshold

        def filter_good_shares(results):
            # Filter results, which is a list of (success, share)
            # pairs.
            return [result[1] for result in results
                    if result is not None and result[0]][:threshold+1]

        def recombine(shares):
            assert len(shares) > threshold
            result = ShareList(shares, threshold+1)
            result.addCallback(filter_good_shares)
            result.addCallback(shamir.recombine)
            return result

        def exchange(share):
#            print "share fp: " + str(share.fp)
            # Send share to all receivers.
            for peer_id in receivers:
                if peer_id != self.id:
                    pc = tuple(self.program_counter)
                    self.protocols[peer_id].sendShare(pc, share)
            # Receive and recombine shares if this player is a receiver.
            if self.id in receivers:
                deferreds = []
                for peer_id in self.players:
                    if peer_id == self.id:
                        d = Share(self, share.field, (share.field(peer_id), share), share.fp)
#                        print "own d is fp? " +str(d.fp)
                    else:
                        d = self._expect_share(peer_id, share.field, share.r)
#                        print "d is fp? "+str(d.fp)
                        d.addCallback(lambda s, peer_id: (s.field(peer_id), s), peer_id)
                    deferreds.append(d)
                return recombine(deferreds)

        result = share.clone()
#        print "result fp: " + str(result.fp)
        self.schedule_callback(result, exchange)

        # do actual communication
        self.activate_reactor()

        if self.id in receivers:
            return result
			
    def get_value(self, val): 
#        print val.fp	
#        return val.signed()
        if val.fp:
            return float(val.signed())/(2**val.r)
        else:
            return val.signed()

    def trunc_value(self, val, res):        
        if val.fp:
            return float(val.signed())/(2**(self.res+res))
        else:
            return float(val.signed())/(2**res)


    def truncPR(self, sh, f, in_fp=False):
        """ Random truncation: precision f bits.
        Protocol 
        """
        field = sh.field
        l = self.options.bit_length
        k = self.options.security_parameter
        ln = (int)(math.log(self.num_players, 2)+1)
        # ensure normal addition
        fp = in_fp or sh.fp
        if isinstance(sh, Share):
            share = self.clone_share_nofp(sh)
        else:
            share = sh.field(sh.value)

        assert field.modulus > 2**(l+f+1)+2**(l+f+k+ln+1), "field is too small"
			
        r_bits = [self.prss_share_random(field, True) for _ in xrange(f)]        
        r_modf = self.bin_comb(r_bits)
        r_divf = self.prss_share_random_max(field, 2**(k+l+1))
        c = share + 2**(l+f) + r_modf + 2**f*r_divf
 		
        def finish(c):
            return reduce(operator.add, [2**i*c.bit(i) for i in xrange(f)])
		
        c = self.open(c)
        self.schedule_callback(c, finish)

        result = (share - c + r_modf)*~field(2**f, fp)*field(1, fp)
        if fp:
            result = self.clone_share_fp(result)
        return result
				
    def mul_public(self, share_a, share_b):
        assert isinstance(share_a, Share), "share_a should be a share"
        assert isinstance(share_b, Share), "share_b should be a share"

        comp_result = gather_shares([share_a, share_b])
        comp_result.addCallback(lambda (a, b): a * b)
        zero_share = self.prss_share_zero(share_a.field, 1)		

        if (not share_a.fp) or (not share_b.fp):
            return self.open(comp_result+zero_share[0], None, 2*self.threshold) 
		
        def truncate(value):
            return value.field(value.value/(2**self.res), True)

        result = self.open(comp_result+zero_share[0], None, 2*self.threshold) 
        self.schedule_callback(result,truncate)
        return result		

    def in_prod_public(self, shares_a, shares_b):
        """Computing inner product of shares_a with shares_b
        using only one round of communication."""
        
        assert len(shares_a) == len(shares_b), \
            "Number of shares_a and shares_b should be equal."

        trunc = reduce(operator.or_, [shares_a[i].fp and shares_b[i].fp for i in xrange(len(shares_a))])

        def computation(share_list):
            n = len(share_list)/2
            s = 0
            for i in xrange(n):
                s += share_list[i] * share_list[n+i]
            return s
		
        comp_result = gather_shares(shares_a+shares_b)
        comp_result.addCallback(computation)
        
        # security fix: add random zero sharing
        zero_share = self.prss_share_zero(shares_a[0].field, 1)		
        if not trunc:
            return self.open(comp_result+zero_share[0], None, 2*self.threshold) 
		
        def truncate(value):
            return value.field(value.value/(2**self.res), True)

        result = self.open(comp_result+zero_share[0], None, 2*self.threshold) 
        self.schedule_callback(result,truncate)
        return result		

    def equal_public(self, share_a, share_b):
        """Public equality test.

        Communication cost: 1 opening (ie. each player sends a share)."""

        return self.equal_zero_public(share_a-share_b)

    def equal_zero_public(self, share_a):
        r = self.prss_share_random(share_a.field)
        ra = self.open(r*share_a)
        return ra.addCallback(lambda x: x == 0)

    @profile
    def neg(self, share_a):
        """Negation of shares.

        Communication cost: none.
        """
        result = share_a.clone()
        result.addCallback(lambda a: -a)
        return result

    @profile
    def add(self, share_a, share_b):
        """FP Addition of shares.

        Communication cost: none.
        """ 

        if not isinstance(share_b, Share):
            if not share_a.fp:
                # Addition with constant. share_a always is a Share by
                # operator overloading in Share. Clone share_a to avoid
                # changing it.
                result = share_a.clone()
                result.addCallback(lambda a, b: b + a, share_b)
                return result
            else:                
                result = share_a.clone()
                result.addCallback(lambda a, b: b + a, share_b)
                return result
			
        result = gather_shares([share_a, share_b])
        result.addCallback(lambda (a, b): a + b)
        return result

    @profile
    def sub(self, share_a, share_b):
        """Subtraction of shares.

        Communication cost: none.
        """
        if not isinstance(share_a, Share):
            share_a = Share(self, share_b.field, share_a)
        if not isinstance(share_b, Share):
            share_b = Share(self, share_a.field, share_b)
        result = gather_shares([share_a, share_b])
        result.addCallback(lambda (a, b): a - b)
        return result

    @profile
    def lin_comb(self, coefficients, shares):
        """Linear combination of shares.

        Communication cost: none. Saves the construction of unnecessary shares
        compared to using add() and mul()."""

        assert len(coefficients) == len(shares), \
            "Number of coefficients and shares should be equal."

        assert all(map(lambda s: not isinstance(s, Share),coefficients)), \
                "Coefficients should not be shares."

        assert all(map(lambda s: isinstance(s, Share),shares)), \
                "Shares should be shares."

        def computation(shares, coefficients):
            return sum(map(operator.mul, coefficients, shares))

        result = gather_shares(shares)
        result.addCallback(computation, coefficients)
        return result

    @profile
    def bin_comb(self, shares):
        """Binary combination of shares."""

        def computation(shares):
            sum = 0
            for i in xrange(len(shares)-1, -1, -1):
                sum = 2*sum + shares[i]    
            return sum

        result = gather_shares(shares)
        result.addCallback(computation)
        return result

    @profile
    def mul(self, share_a, share_b, dotrunc=True):
        """Multiplication of FP shares.

        Communication cost: 1 Shamir sharing plus plus.
        """
        
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
    	
		# At this point both share_a and share_b must be Share
        # objects. So we wait on them, multiply and reshare.

        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        result = gather_shares([share_a, share_b])
        result.addCallback(lambda (a, b): a*b)
        self.schedule_callback(result, share_recombine)
        if trunc:
            self.schedule_callback(result, self.truncPR, self.res, True)

        # do actual communication
        self.activate_reactor()

        return result

    def invert(self, share_a):
        r = self.prss_share_random(share_a.field)
        ra = self.open(r*self.clone_share_nofp(share_a))
        return ra.addCallback(lambda ra, r: ~ra*r*2**(self.res*share_a.fp), r)

    @profile
    def gauss(self, A, d, b, c):
        """Gaussian elimination A:= A d - b c

        Communication cost: m * n Shamir sharings.
        """
        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)

            return result

        def computation(x, m):
            n = len(x)/(m+1) - 1
            mn = m*n
            for i in xrange(m):
                for j in xrange(n):
                    trunc = (x[n*i+j].fp and x[mn].fp) or (x[mn+1+i].fp and x[mn+1+m+j].fp)
                    x[n*i+j] = share_recombine(x[n*i+j] * x[mn] - x[mn+1+i] * x[mn+1+m+j])
                    if trunc:
                        x[n*i+j] = self.truncPR(x[n*i+j], self.res, True)
 
            del x[mn:]
            return x
        
        def peek(s, k):
            return s[k]
            
        def wrap(s, m, n):
            return [[s.clone().addCallback(peek, n*i+j) for j in xrange(n)] for i in xrange(m)]
            
        result = gather_shares(reduce(operator.add, A) + [d] + b + c)
        self.schedule_callback(result, computation, len(A))

        # do actual communication
        self.activate_reactor()
        
        return wrap(result, len(A), len(A[0]))
        #return result

    @profile
    def scalar_mul(self, a, x):
        """Scalar multiplication x:= a x

        Communication cost: n Shamir sharings.
        """
        
        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        def computation(x):
            a = x.pop()
            for i in xrange(len(x)):
                trunc = a.fp and x[i].fp
                x[i] = share_recombine(x[i]*a)
                if trunc:
                    x[i]= self.truncPR(x[i], self.res, True)
            return x    
         
        def peek(s, k):
            return s[k]
            
        def wrap(s, n):
            return [s.clone().addCallback(peek, i) for i in xrange(n)]
            
        result = gather_shares(x+[a])
        self.schedule_callback(result, computation)
        
        # do actual communication
        self.activate_reactor()
        
        return wrap(result, len(x))

    @profile
    def scalar_prod(self, x, y):
        """Scalar multiplication of vectors x and y

        Communication cost: n Shamir sharings.
        """
        
        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        def computation(x):
            n = len(x)/2
            for i in xrange(n):
                trunc = x[i].fp and x[n+i].fp
                x[i] = share_recombine(x[i]*x[n+i])
                if trunc:
                    x[i]= self.truncPR(x[i], self.res, True)
            del x[n:]
            return x    
        
        def peek(s, k):
            return s[k]
            
        def wrap(s, n):
            return [s.clone().addCallback(peek, i) for i in xrange(n)]
            
        result = gather_shares(x+y)
        self.schedule_callback(result, computation)
        
        # do actual communication
        self.activate_reactor()
        
        return wrap(result, len(x))
  
    def matrix_prod(self, A, B):
        """Computing matrix product of A with transposed B
        using only one round of communication."""
		
        if (isinstance(A[0], list) and isinstance(B[0], list)):	
            trunc = [[reduce(operator.or_, [A[i][j].fp and B[k][j].fp for j in xrange(len(A[0]))]) for k in xrange(len(B))] for i in xrange(len(A))]      
        else:
            if(isinstance(A[0], list)):
                trunc = [[reduce(operator.or_, [A[i][j].fp and B[j].fp for j in xrange(len(A[i]))])] for i in xrange(len(A))]      
            else:
                if(isinstance(B[0], list)):
                    trunc = [[reduce(operator.or_, [A[j].fp and B[k][j].fp for j in xrange(len(A[i]))]) for k in xrange(len(B[i]))]]      
                else:
                    trunc = [[reduce(operator.or_, [A[j].fp and B[j].fp for j in xrange(len(A))])]]      


				
        def computation(x, ma, mb, trunc):
            n = len(x)/(ma+mb)
            C = []
            for ia in xrange(ma):
                for ib in xrange(mb):
                    s = 0
                    for i in xrange(n):
                        s += x[ia*n+i]*x[(ma+ib)*n+i]
                    z = share_recombine(s)
                    #z.fp = s.fp
                    if trunc[ia][ib]:
                        z = self.truncPR(z, self.res, True)
                    C.append(z)
            return C

        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        def peek(s, k):
            return s[k]
            
        def wrap(s, ma, mb):
            return [[s.clone().addCallback(peek, mb*ia+ib) for ib in xrange(mb)] for ia in xrange(ma)]

        result = gather_shares(reduce(operator.add, A)+reduce(operator.add, B))
        self.schedule_callback(result, computation, len(A), len(B), trunc)

        # do actual communication
        self.activate_reactor()
        
        return wrap(result, len(A), len(B))
                            
    def in_prod(self, shares_a, shares_b):
        """Computing inner product of shares_a with shares_b
        using only one round of communication."""
        
        assert len(shares_a) == len(shares_b), \
            "Number of shares_a and shares_b should be equal."

        trunc = reduce(operator.or_, [shares_a[i].fp and shares_b[i].fp for i in xrange(len(shares_a))])

        # avoid errors due to inconsistent sharings

        def computation(share_list):
            n = len(share_list)/2
            s = 0
            for i in xrange(n):
                s += share_list[i] * share_list[n+i] 
            return s

        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        result = gather_shares(shares_a+shares_b)
        result.addCallback(computation)
        self.schedule_callback(result, share_recombine)
        if trunc:
            self.schedule_callback(result, self.truncPR, self.res, True)
        
        # do actual communication
        self.activate_reactor()

        return result
        
    def in_prod2(self, shares_a):
        """Computing inner product of shares_a with itself
        using only one round of communication."""

        trunc = reduce(operator.or_, [shares_a[i].fp for i in xrange(len(shares_a))])

        # avoid errors due to inconsistent sharings
        
        def computation(share_list):
            n = len(share_list)
            s = 0
            for i in xrange(n):
                s += share_list[i]*share_list[i]
            return s

        def share_recombine(number):
            shares = shamir.share(number, self.threshold, self.num_players)

            exchanged_shares = []
            for peer_id, share in shares:
                d = self._exchange_shares(peer_id.value, share)
                d.addCallback(lambda share, peer_id: (peer_id, share), peer_id)
                exchanged_shares.append(d)

            # Recombine the first 2t+1 shares.
            result = gather_shares(exchanged_shares[:2*self.threshold+1])
            result.addCallback(shamir.recombine)
            return result

        result = gather_shares(shares_a)
        result.addCallback(computation)
        self.schedule_callback(result, share_recombine)
        if trunc:
            self.schedule_callback(result, self.truncPR, self.res, True)
      
        # do actual communication
        self.activate_reactor()

        return result        

    def fan_in_mulL(self, shares, cp=True):
        """
		Logarithmic rounds multiplication
        """

        assert isinstance(shares, list), "fan-in-mulL: require list as input"         
        assert len(shares) > 0, "fan-in-MulL: Cannot multiply an empty list"

        if cp:
            tmp = copy.copy(shares)
        else:
            tmp = shares
        
        while len(tmp) > 1:
            h = []
            while len(tmp) > 1: h.append(tmp.pop() * tmp.pop())
            h.extend(tmp)
            tmp = h
        return tmp[0]

    def fan_in_mulC(self, shares):
        """
        Unbounded Fan-in Multiplication. Protocol 4.16 in Thesis [H12]
        """
        assert isinstance(shares, list), "prefix-mulC: require list as input"         
        assert len(shares) > 0, "prefix-MulC: Cannot multiply an empty list"
        		
        field = shares[0].field
        n_shares = len(shares)
        
        if n_shares==1:
            return shares[0]
		# TODO: check for zeroness
		
        def u_list_ready(ulist):

            wlist = [rlist[0]]+map(operator.mul, vlist, map(operator.invert, ulist[:-1]))
            m_list = map(self.mul_public, shares, wlist)

            def m_list_ready(mlist):

                rslt = reduce(operator.mul, mlist)
                rslt = rslt * slist[-1]*~ulist[-1]
			
                trunc = reduce(operator.add, [share.fp for share in shares])
                print trunc

                def truncate(val, tr):
                    if tr > 1:
                       val = self.truncPR(val, (tr-1)*self.res, True)
                    return val

                rslt.addCallback(truncate, trunc)
                return rslt

            tmp = gather_shares(m_list)
            self.schedule_callback(tmp, m_list_ready)		
            return tmp

        rlist = [self.prss_share_random(field, False) for i in xrange(n_shares)]
        slist = [self.prss_share_random(field, False) for i in xrange(n_shares)]
        u_list = map(self.mul_public,rlist,slist)
        vlist = map(self.mul, rlist[1:], slist[:-1])

        result = gather_shares(u_list)        
        self.schedule_callback(result,u_list_ready)		
        return result

    def prefix_mulC(self, shares):
        """
        Unbounded Fan-in Prefix Multiplication. Protocol 4.17 in Thesis [H12]
        """

        assert isinstance(shares, list), "prefix-mulC: require list as input"         
        assert len(shares) > 0, "prefix-MulC: Cannot multiply an empty list"
        		
        field = shares[0].field
        n_shares = len(shares)
        
        if n_shares==1:
            return shares[0]
		# TODO: check for zeroness
		
        def u_list_ready(ulist):

            wlist = [rlist[0]]+map(operator.mul, vlist, map(operator.invert, ulist[:-1]))
            m_list = map(self.mul_public, shares, wlist)

            def m_list_ready(mlist):

                for i in xrange(1, n_shares):
	                mlist[i] = mlist[i-1]*mlist[i]
        
                results = [shares[0]]+[slist[i]*~ulist[i]*mlist[i] for i in xrange(1,n_shares)]
                trunc = [reduce(operator.add, [bool(shares[i].fp) for i in xrange(j+1)]) for j in xrange(1,n_shares)]

                def truncate(list, tr):
                    for i in xrange(1,len(list)):
                        if tr[i-1] > 1:
                            list[i] = self.truncPR(list[i], (tr[i-1]-1)*self.res, True)
                    return list

                results = gather_shares(results)
                results.addCallback(truncate, trunc)

                return results
            tmp = gather_shares(m_list)
            self.schedule_callback(tmp, m_list_ready)		
            return tmp

        def peek(s, k):
            return s[k]
            
        def wrap(s):
            return [s.clone().addCallback(peek, i) for i in xrange(n_shares)]		

        rlist = [self.prss_share_random(field, False) for i in xrange(n_shares)]
        slist = [self.prss_share_random(field, False) for i in xrange(n_shares)]
        u_list = map(self.mul_public,rlist,slist)
        vlist = map(self.mul, rlist[1:], slist[:-1])

        result = gather_shares(u_list)        
        self.schedule_callback(result,u_list_ready)		
        return wrap(result)

    def lsb(self, sh):
        """
        Least Significant Bit Gate [ST06]
        """
        share = self.clone_share_nofp(sh)
		
        field = share.field
        l = self.options.bit_length
        k = self.options.security_parameter
        ln = (int)(math.log(self.num_players,2)+1)
		
        assert share.field.modulus > 2**(l+1)+2**(l+k+ln+1), "field is too small"
		
        b = self.prss_share_random(field, True)
        r = self.prss_share_random_max(field, 2**(l+k))

        def finish_computation(c, b):
            return c.bit(0)+b-2*c.bit(0)*b

        c = self.open(share+b+2*r)
        self.schedule_callback(c, finish_computation, b)
        return c

    def pow(self, share, exponent):
        """Exponentation of a share to an integer by square-and-multiply."""

        assert isinstance(exponent, (int, long)), "Exponent must be an integer"
        assert exponent >= 0, "Exponent must be non-negative"

        if exponent == 0:
            return 1
        elif exponent % 2 == 0:
            tmp = share ** (exponent / 2)
            return tmp * tmp
        else:
            return share * (share ** (exponent-1))

    def xor(self, share_a, share_b):
        field = share_a.field
        if not isinstance(share_b, Share):
            if not isinstance(share_b, FieldElement):
                share_b = field(share_b)
            share_b = Share(self, field, share_b)

        if field is GF256:
            return share_a + share_b
        else:
            return share_a + share_b - 2 * share_a * share_b

    def prss_key(self):
        """Create unique key for PRSS.

        This increments the program counter and returns it as a tuple.
        Each straight-line program (typically a callback attached to
        some :class:`Deferred`) is executed in a context with unique
        starting program counter. This ensures that consequetive calls
        to PRSS-related methods will use unique program counters.
        """

        # This is called by every function using PRSS, so do it here.
        # If the assertion is not met, things go wrong, i.e. the PRSS
        # functions generate shares with higher degrees than what
        # open() and mul() expect.
        assert self.threshold >= \
               len(self.players) - len(self.players[self.id].keys.keys()[0]), \
               "PRSS functions have higher threshold than the runtime."

        self.increment_pc()
        return tuple(self.program_counter)

    def prss_share(self, inputters, field, element=None):
        """Creates pseudo-random secret sharings.

        This protocol creates a secret sharing for each player in the
        subset of players specified in *inputters*. Each inputter
        provides an integer. The result is a list of shares, one for
        each inputter.

        The protocol uses the pseudo-random secret sharing technique
        described in the paper "Share Conversion, Pseudorandom
        Secret-Sharing and Applications to Secure Computation" by
        Ronald Cramer, Ivan Damgård, and Yuval Ishai in Proc. of TCC
        2005, LNCS 3378. `Download
        <http://www.cs.technion.ac.il/~yuvali/pubs/CDI05.ps>`__

        Communication cost: Each inputter does one broadcast.
        """
        # Verifying parameters.
        if element is None:
            assert self.id not in inputters, "No element given."
        else:
            assert self.id in inputters, \
                "Element given, but we are not sharing?"

        n = self.num_players

        # Key used for PRSS.
        key = self.prss_key()

        # The shares for which we have all the keys.
        all_shares = []

        # Shares we calculate from doing PRSS with the other players.
        tmp_shares = {}

        prfs = self.players[self.id].dealer_prfs(field.modulus)

        # Compute and broadcast correction value.
        if self.id in inputters:
            for player in self.players:
                share = prss(n, player, field, prfs[self.id], key)
                all_shares.append((field(player), share))
            shared = shamir.recombine(all_shares[:self.threshold+1])
            correction = element - shared
            # if this player is inputter then broadcast correction value
            # TODO: more efficient broadcast?
            pc = tuple(self.program_counter)
            for peer_id in self.players:
                if self.id != peer_id:
                    self.protocols[peer_id].sendShare(pc, correction)

        # Receive correction value from inputters and compute share.
        result = []
        for player in inputters:
            tmp_shares[player] = prss(n, self.id, field, prfs[player], key)
            if player == self.id:
                d = Share(self, field, correction)
            else:
                d = self._expect_share(player, field)
            d.addCallback(lambda c, s: s + c, tmp_shares[player])
            result.append(d)

        # Unpack a singleton list.
        if len(result) == 1:
            return result[0]
        else:
            return result

    def prss_share_random_bit(self, field):
        return self.prss_share_random(field, True)
            
    def prss_share_random(self, field, binary=False):
        """Generate shares of a uniformly random element from the field given.

        If binary is True, a 0/1 element is generated. No player
        learns the value of the element.

        Communication cost: none if binary=False, 1 open otherwise.
        """
        if field is GF256 and binary:
            modulus = 2
        else:
            modulus = field.modulus

        # Key used for PRSS.
        prss_key = self.prss_key()
        prfs = self.players[self.id].prfs(modulus)
        share = prss(self.num_players, self.id, field, prfs, prss_key)

        if field is GF256 or not binary:
            return Share(self, field, share)

        # Open the square and compute a square-root
        result = self.open(Share(self, field, share*share),
                           threshold=2*self.threshold)
        
        def finish(square, share, binary):
            if square == 0:
                # We were unlucky, try again...
                return self.prss_share_random(field, binary)
            else:
                # We can finish the calculation
                root = square.sqrt()
                # When the root is computed, we divide the share and
                # convert the resulting -1/1 share into a 0/1 share.
                return Share(self, field, (share/root + 1) / 2)

        self.schedule_callback(result, finish, share, binary)
        return result

    def prss_share_random_max(self, field, max):
        # Key used for PRSS.
        prss_key = self.prss_key()
        prfs = self.players[self.id].prfs(max)
        share = prss(self.num_players, self.id, field, prfs, prss_key)
        return Share(self, field, share)

    def prss_share_random_double_max(self, field1, field2, max):
        # Key used for PRSS.
        prss_key = self.prss_key()
        prfs = self.players[self.id].prfs(max)
        share1, share2 = prss_double(self.num_players, self.id, field1, field2, prfs, prss_key)
        return Share(self, field1, share1), Share(self, field2, share2)

    def prss_share_random_multi(self, field, quantity, binary=False):
        """Does the same as calling *quantity* times :meth:`prss_share_random`,
        but with less calls to the PRF. Sampling of a binary element is only
        possible if the field is :class:`GF256`.

        Communication cost: none.
        """
        assert not binary or field == GF256, "Binary sampling not possible " \
            "for this field, use prss_share_random()."

        if field is GF256 and binary:
            modulus = 2
        else:
            modulus = field.modulus

        # Key used for PRSS.
        prss_key = self.prss_key()
        prfs = self.players[self.id].prfs(modulus ** quantity)
        shares = prss_multi(self.num_players, self.id, field, prfs, prss_key,
                            modulus, quantity)
        return [Share(self, field, share) for share in shares]

    def prss_share_zero(self, field, quantity):
        """Generate *quantity* shares of the zero element from the
        field given.

        Communication cost: none.
        """
        # Key used for PRSS.
        prss_key = self.prss_key()
        prfs = self.players[self.id].prfs(field.modulus)
        zero_share = prss_zero(self.num_players, self.threshold, self.id,
                               field, prfs, prss_key, quantity)
        return [Share(self, field, zero_share[i]) for i in range(quantity)]

    def prss_double_share(self, field, quantity):
        """Make *quantity* double-sharings using PRSS.

        The pair of shares will have degree t and 2t where t is the
        default threshold for the runtime.
        """
        r_t = self.prss_share_random_multi(field, quantity)
        z_2t = self.prss_share_zero(field, quantity)
        return (r_t, [r_t[i] + z_2t[i] for i in range(quantity)])

    def prss_share_bit_double(self, field):
        """Share a random bit over *field* and GF256.

        The protocol is described in "Efficient Conversion of
        Secret-shared Values Between Different Fields" by Ivan Damgård
        and Rune Thorbek available as `Cryptology ePrint Archive,
        Report 2008/221 <http://eprint.iacr.org/2008/221>`__.
        """
        n = self.num_players
        k = self.options.security_parameter
        prfs = self.players[self.id].prfs(2**k)
        prss_key = self.prss_key()

        b_p = self.prss_share_random(field, binary=True)
        r_p, r_lsb = prss_lsb(n, self.id, field, prfs, prss_key)

        b = self.open(b_p + r_p)
        # Extract least significant bit and change field to GF256.
        b.addCallback(lambda i: GF256(i.value & 1))
        b.field = GF256

        # Use r_lsb to flip b as needed.
        return (b_p, b ^ r_lsb)

    def prss_shamir_share_bit_double(self, field):
        """Shamir share a random bit over *field* and GF256."""
        n = self.num_players
        k = self.options.security_parameter
        prfs = self.players[self.id].prfs(2**k)
        prss_key = self.prss_key()
        inputters = range(1, self.num_players + 1)

        ri = rand.randint(0, 2**k - 1)
        ri_p = self.shamir_share(inputters, field, ri)
        ri_lsb = self.shamir_share(inputters, GF256, ri & 1)

        r_p = reduce(self.add, ri_p)
        r_lsb = reduce(self.add, ri_lsb)

        b_p = self.prss_share_random(field, binary=True)
        b = self.open(b_p + r_p)
        # Extract least significant bit and change field to GF256.
        b.addCallback(lambda i: GF256(i.value & 1))
        b.field = GF256

        # Use r_lsb to flip b as needed.
        return (b_p, b ^ r_lsb)

    def powerchain(self, share, max):
        """Returns the list ``[share, share**2, share**4, ..., share**(2**max)]``."""
        result = [share]
        for i in range(max):
            share = share * share
            result.append(share)
        return result

    @preprocess("prss_powerchains")
    def prss_powerchain(self, max=7):
        """Generate a random secret share in GF256 and returns
        ``[share, share**2, share**4, ..., share**(2**max)]``."""
        share = self.prss_share_random(GF256)
        return self.powerchain(share, max)

    def prss_powerchains(self, max=7, quantity=20):
        """Does *quantity* times the same as :meth:`prss_powerchain`.
        Used for preprocessing."""
        quantity = min(quantity, 20)
        shares = self.prss_share_random_multi(GF256, quantity)
        return [gatherResults(self.powerchain(share, max)) for share in shares]

    def input(self, inputters, field, number=None, threshold=None):
        """Input *number* to the computation.

        The input is shared using the :meth:`shamir_share` method.
        """
        return self.shamir_share(inputters, field, number, threshold)

    # TODO: can we remove fp functionality here or does that break stuff?
    def shamir_share(self, inputters, field, number=None, threshold=None):
        """Secret share *number* over *field* using Shamir's method.

        The number is shared using polynomial of degree *threshold*
        (defaults to :attr:`threshold`). Returns a list of shares
        unless there is only one inputter in which case the
        share is returned directly.

        In code it is used like this::

            a, b, c = runtime.shamir_share([1, 2, 3], Zp, x)

        where ``Zp`` is a field and ``x`` is a Python integer holding
        the input of each player (three inputs in total).

        If only a subset of the players provide input it looks like
        this::

            if runtime.myid == 1:
                a = runtime.shamir_share([1], Zp, x)
            else:
                a = runtime.shamir_share([1], Zp)

        Instead of branching when calling :meth:`shamir_share`, one
        can give ``None`` as input::

            if runtime.myid == 1:
                x = int(raw_input("Input x: "))
            else:
                x = None
            a = runtime.shamir_share([1], Zp, x)

        which might be practical in some cases.

        Communication cost: n elements transmitted.
        """
        assert number is None or self.id in inputters
        if threshold is None:
            threshold = self.threshold

        fp = isinstance(number, float)
        if fp:
            number = int(2**self.res*number)
            
        results = []
        for peer_id in inputters:
            # Unique program counter per input.
            self.increment_pc()

            if peer_id == self.id:
                pc = tuple(self.program_counter)
                shares = shamir.share(field(number), threshold,
                                      self.num_players)
                shares = [[myid, val.value << 1 ^ int(fp)] for myid, val in shares]
                for other_id, share in shares:
                    if other_id.value == self.id:
                        results.append(Share(self, field, field(share >> 1, fp), fp))
                    else:
                        self.protocols[other_id.value].sendShare_shamir(pc, share)
            else:
                results.append(self._expect_share_shamir(peer_id, field))

        # do actual communication
        self.activate_reactor()

        # Unpack a singleton list.
        if len(results) == 1:
            return results[0]
        else:
            return results

    def shamir_share_fp(self, inputters, field, number=None, threshold=None):
        assert number is None or self.id in inputters
        if threshold is None:
            threshold = self.threshold

        if number!=None: number = int(2**self.res*number)

        results = []
        for peer_id in inputters:
            # Unique program counter per input.
            self.increment_pc()

            if peer_id == self.id:
                pc = tuple(self.program_counter)
                shares = shamir.share(field(number), threshold,
                                      self.num_players)
                shares = [[myid, val.value << 1 ^ True] for myid, val in shares]
                for other_id, share in shares:
                    if other_id.value == self.id:
                        results.append(Share(self, field, field(share >> 1, True), True))
                    else:
                        self.protocols[other_id.value].sendShare_shamir(pc, share)
            else:
                results.append(self._expect_share_shamir(peer_id, field,True))

        # do actual communication
        self.activate_reactor()

        # Unpack a singleton list.
        if len(results) == 1:
            return results[0]
        else:
            return results

    def fpval(self, Zp, val):
        return Zp(int(round(val*(1<<self.options.res))), self.options.res)
        
    def lsb(self, sh):
        """
        Least Significant Bit Gate [ST06]
        """
        share = sh.clone()
                
        field = share.field
        l = self.options.bit_length-10 # TODO: fix
        k = self.options.security_parameter
        ln = (int)(math.log(self.num_players,2)+1)
                
        assert share.field.modulus > 2**(l+1)+2**(l+k+ln+1), "field is too small"
                
        b = self.prss_share_random(field, True)
        r = self.prss_share_random_max(field, 2**(l+k))

        def finish_computation(c, b):
            return c.bit(0)+b-2*c.bit(0)*b

        c = self.open(share+b+2*r)
        self.schedule_callback(c, finish_computation, b)
        return c
      
      
    def prefix_carry(self, bits_x, bits_y):
        """Compute the prefix carries when adding the binary 
        representation of share_x + share_y.
        Protocol 4.23 of [SH12].
        """
        
        assert len(bits_x) == len(bits_y), "prefix carry: lists should have equal length"

        k = len(bits_x)
        carry_list = [bits_x[i]*bits_y[i] for i in xrange(k)]
        xor_list = [bits_x[i]+bits_y[i] - 2*carry_list[i] for i in xrange(k)]
                
        for i in xrange(1,int(math.ceil(math.log(k,2)))+1):
            for j in xrange(1,int(math.ceil(float(k)/(2**i)))+1):
                u = 2**(i-1)+(j-1)*2**i-1
                v = min([2**(i-1), k-u-1])
                if v > 0:
                    for z in xrange(1,v+1):
                        carry_list[u+z] = carry_list[u+z]+xor_list[u+z]*carry_list[u]
                                                
                        if j != 1:
                            xor_list[u+z] = xor_list[u+z]*xor_list[u]

        return carry_list
    
    def addBitwise(self, bits_x, bits_y):
        """given the binary representation of share_x and share_y 
        compute the binary representation of share_x + share_y. 
        Protocol 4.24 of [SH12].
        """

        assert len(bits_x) == len(bits_y), "prefix carry: lists should have equal length"
        k = len(bits_x)
                                
        carry_list = self.prefix_carry(bits_x, bits_y)
        
        s_init = bits_x[0] + bits_y[0] - 2*carry_list[0]

        return [s_init]+[bits_x[i]+bits_y[i] + carry_list[i-1] - 2*carry_list[i] for i in xrange(1,k)]+[carry_list[k-1]]
        
    def zp_secret_index(self, N):
        self.N=N
        
        if N%4 == 0: N += 1
        p = 1 + N*(1 + (N**2)%4 + 4*(2**(self.options.bit_length+self.options.security_parameter+1)/N))
        while not mpz(p).is_prime(): p += 4 * N
        
        self.Zp=GF(p)

        w = self.Zp(2)**((p-1)/N)  # lucky choice 2
        self.w_powers = [self.Zp(1)]
        for _ in xrange(N-1): self.w_powers.append(self.w_powers[-1] * w)
        
    def vc_secret_index(self, N):
        self.N=2  # 2^k|p-1
        while self.N < N: self.N=self.N*2
        w=self.Zp(7)**((self.Zp.modulus-1)/self.N)
        self.w_powers = [self.Zp(1)]
        for _ in xrange(self.N-1): self.w_powers.append(self.w_powers[-1] * w)
        
    def secret_index(self, ix):
        return Share(self, self.Zp, self.w_powers[-ix])
        
    def shareZp(self, val):
        return Share(self, self.Zp, self.Zp(val))
        
    def secret_index_get(self, values, indices, start, end):
        def pow_list(a, x, n):
            if n==1:
                return [a]
            xs = []
            x2s = pow_list(a, x**2, n/2)
            for x2 in x2s:
                xs.append(x2)
                xs.append(x2*x)
            if n%2 == 1:
                xs.append(xs[-1]*x)
            return xs
        solution = [0]*(end-start)
        for i in xrange(len(indices)):
            x_powers = pow_list(values[i]*~self.Zp(self.N), indices[i], self.N)
            for j in xrange(start, end):
                coefs = [self.w_powers[(j*k)%self.N] for k in xrange(self.N)]
                solution[j-start] += self.lin_comb(coefs,x_powers)
        return solution
    
    def bit_decompose(self, share, l=-1):
        """Extract the bits of share given that the modulus 
        much larger than the bit-length of share. This is 
        the bit-decomposition protocol of [ST06].
        """
        field = share.field

        a = self.clone_share_nofp(share)
                
        if l==-1: l = self.options.bit_length
        k = self.options.security_parameter
        ln = (int)(math.log(self.num_players,2)+1)
        
#        print "minimal modulus:", 2**(l+1)+2**(l+k+ln+1)

#        assert field.modulus > 2**(l+1)+2**(l+k+ln+1), "field is too small" TODO
                
        r_bits = [self.prss_share_random(field, True) for _ in xrange(l)]   
        r_modl = self.bin_comb(r_bits)
                                                        
        r_divl = self.prss_share_random_max(field, 2**k)

        c = self.open(a + 2**(l) - r_modl + 2**(l)*r_divl)

        def extract_bits(value, size):
            return [field(value.bit(i)) for i in xrange(size)]

        self.schedule_callback(c, extract_bits, l)
        self.schedule_callback(c, self.addBitwise, r_bits)

                                
        def peek(s, k):
            return s[k]
            
        def wrap(s):
            return [s.clone().addCallback(peek, i) for i in xrange(l)]          

        return wrap(c)      
        
    def argmin(self, xs, gte_min):
        n = len(xs)
        if n == 1:
            return (xs[0], [1])
        if n == 2:
            c, min = gte_min(xs[0], xs[1])
            return (min, [1-c, c])
        x2s = []
        c2s = []
        while len(xs) > 1:       
            a1 = xs.pop()
            a0 = xs.pop()
            c, min = gte_min(a0, a1)
            x2s.append(min)
            c2s.append(c)
        if n%2 == 1:
            x2s.append(xs.pop())
        min, index2 = self.argmin(x2s, gte_min)
        index = []
        if n%2 == 1:
            index.append(index2[n/2])
        for i in xrange(n/2-1, -1, -1):
            b = index2[i]*c2s[i]
            index.append(index2[i] - b)
            index.append(b)
        return (min, index)

    def int_minimal(self, xs):
        # def printer(res):
            # print "comparison result",
            # for i in res:
                # if i.value > 21888242871839275222246405745257275088548364400416034343698204186575808495617/2:
                    # print i.value - 21888242871839275222246405745257275088548364400416034343698204186575808495617,
                # else:
                    # print i.value,
            # print
        def int_gte_min(x0, x1):
            c = x0 >= x1
            # res = gather_shares(map(self.open, [x0, x1, c]))
            # res.addCallback(printer)
            return (c, c * (x1 - x0) + x0)
        return self.argmin(xs, int_gte_min)

    def frac_minimal(self, xs):
        def frac_gte_min((n0, d0), (n1, d1)):
            c = self.in_prod([n0,d0],[d1,-n1]) >= 0
            h = self.scalar_mul(c, [n1-n0, d1-d0])
            return (c, (h[0]+n0, h[1]+d0))
        return self.argmin(xs, frac_gte_min)                
        
    # Protocol 4.26 from [dH12]
    # TODO: small optimization: if know that bitlength is at most a certain value, then do it more easily
    # https://eprint.iacr.org/2014/137.pdf
    @viff.inlinecb.viffinlinecb
    def eqz(self, x, k=-1):
        yield self, Share(self, x.field)
        
        if k==1: viff.inlinecb.returnValue(x)
        if k==2: viff.inlinecb.returnValue(1-(x-1)*(x-2)*(x-3)*(~x.field(-6)))
        if k==-1: k = self.options.bit_length
        
        r_bits = [self.prss_share_random(x.field, True) for _ in xrange(k)]
        r_divl = self.prss_share_random_max(x.field, 2**self.options.security_parameter)

        c = yield self.open(x + sum([(2**i)*r_bits[i] for i in xrange(k)]) + (2**(k))*r_divl)
        c_bits = [x.field(c.bit(i)) for i in xrange(k)]
        
        d = sum([c_bits[i]+r_bits[i]-2*c_bits[i]*r_bits[i] for i in xrange(k)])
        
        rec = yield self.eqz(d, k.bit_length())    
        viff.inlinecb.returnValue(rec)        
        