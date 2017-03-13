#SH new mixin for secure fixed point division 
# Copyright 2013 TUE.
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

"""Division (related) protocols. The mixin classes defined here provide a
:meth:`division`,
:meth:`reciprocal`,
:meth:`bit-decompose`,
:meth:`addBitwise`,
:meth:`prefix_carry`
methods to the
:class:`~viff.runtime.Runtime` they are mixed with.
"""

import math, operator

from viff.util import rand, profile, find_prime
from viff.runtime import Share, gather_shares
from viff.passive import PassiveRuntime
from viff.active import ActiveRuntime
from viff.field import GF256, GF, FieldElement

from twisted.internet.defer import Deferred, DeferredList, gatherResults

class DivisionSH12Mixin:

    """Efficient comparison by Catrina and de Hoogh. This mixin provides a
    :meth:`division`,
    :meth:`reciprocal`,
    :meth:`bit-decompose`,
    :meth:`addBitwise`,
    :meth:`prefix_carry`, method which can compare Zp field
    elements and gives a secret result shared over Zp.
    """

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

        assert field.modulus > 2**(l+1)+2**(l+k+ln+1), "field is too small"
 		
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

    def prefix_orC(self, bits):
        """Compute prefix Or in Constant Rounds
        Protocol 4.18 of [SH12]"""

        z_list = [bit+1 for bit in bits]
        z_list = self.prefix_mulC(z_list)
        return [bits[0]]+[1-self.lsb(entry) for entry in z_list[1:]]

    def norm(self, sh):
        """ Normalization of share. Preprocessing stage of division
        Protocol 4.36 of [SH12]."""
		
        share = self.clone_share_nofp(sh)
        bits = self.bit_decompose(share)

        d_list = self.prefix_orC([bits[-i] for i in xrange(1,len(bits)+1)])
        b_list = [d_list[0]]+[d_list[i] - d_list[i-1] for i in xrange(1,len(bits))]

        v = self.bin_comb(b_list)
				
        return [share*v, v] 

    def reciprocal(self, share):
        """compute 1/share with runtime.res precision. Newton Rhapson Method.
        Protocol 4.35 of [SH12]
        """

        l = self.options.bit_length
		
        assert share.field.modulus > 2**(3*l+self.res+3)
		        
        x, v = self.norm(share)
	
        theta = int(math.ceil(math.log(float(l+1)/3.5,2)))
        y = int(2.9142*2**l)-2*x
        for i in xrange(theta):
            y = y*(2**(2*l+1)-y*x)
            y = self.truncPR(y, 2*l)
        result_tmp = y*v
        result = self.clone_share_fp(result_tmp)

        if share.fp:
            return self.truncPR(result, 2*l-2*self.res)
        return self.truncPR(result, 2*l-self.res)

    def rec(self, share):
        """compute 1/share. with options.bit_length precision. 
        Newton Rhapson Method. Protocol 4.35 of [SH12]. Used as a submodule.
        """

        l = self.options.bit_length
		
        assert share.field.modulus > 2**(3*l+self.res+3)
		        
        x, v = self.norm(share)
	
        theta = int(math.ceil(math.log(float(l+1)/3.5,2)))
        y = int(2.9142*2**l)-2*x
        for i in xrange(theta):
            y = y*(2**(2*l+1)-y*x)
            y = self.truncPR(y, 2*l)
        return y*v

    def div(self, share_a, share_b):
        """Secure division using Newton Rhapson Method. 
        Protocol 4.34 of [SH12].
        """

        l = self.options.bit_length
        field = getattr(share_a, "field", getattr(share_b, "field", None))
		
        if isinstance(share_a, int) and share_a == 1:
            return self.reciprocal(share_b)
		
        if not isinstance(share_b, Share):
            value = share_b
            if hasattr(share_b, "value"):
                value = share_b.signed()
            value = field(int(float(1)/value*2**self.res), True)
            result = share_a * value
            if share_a.fp:
                return self.truncPR(result, self.res)
            return result

        result_tmp = self.rec(share_b)*share_a
        result = self.clone_share_fp(result_tmp)
				
        if not getattr(share_a, "fp", False):
            if share_b.fp:
                return self.truncPR(result, 2*l-2*self.res, True)
            return self.truncPR(result, 2*l-self.res, True)				
        if share_b.fp:
            return self.truncPR(result, 2*l-self.res, True)        
        return self.truncPR(result, 2*l, True)

    def floordiv(self, share_a, share_b):
        """Secure division using Newton Rhapson Method. 
        Protocol 4.34 of [SH12].
        """

        l = self.options.bit_length
        field = getattr(share_a, "field", getattr(share_b, "field", None))
		
        if isinstance(share_a, int) and share_a == 1:
            return self.reciprocal(share_b)
		
        if not isinstance(share_b, Share):
            value = share_b
            if hasattr(share_b, "value"):
                value = share_b.signed()
            value = field(int(float(1)/value*2**self.res), True)
            result = share_a * value
            if share_a.fp:
                return self.truncPR(result, 2*self.res)
            return self.truncPR(result, self.res)

        result = self.rec(share_b)*share_a
        result.fp = False
        if not getattr(share_a, "fp", False):
            if share_b.fp:
                return self.truncPR(result, 2*l-self.res)
            return self.truncPR(result, 2*l)				
        if share_b.fp:
            return self.truncPR(result, 2*l)        
        return self.truncPR(result, 2*l+self.res)