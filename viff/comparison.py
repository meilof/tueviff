#BS Improved and cleaned up version of Toft07 
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

"""Comparison protocols. The mixin classes defined here provide a
:meth:`greater_than_equal` method to the
:class:`~viff.runtime.Runtime` they are mixed with.
"""

import math, operator

from viff.util import rand, profile, find_prime
from viff.runtime import Share, gather_shares
from viff.passive import PassiveRuntime
from viff.active import ActiveRuntime
from viff.field import GF256, GF, FieldElement

class ComparisonToft07Mixin:

    """Efficient comparison by Tomas Toft 2007. This mixin provides a
    :meth:`greater_than_equal` method which can compare Zp field
    elements and gives a secret result shared over Zp.
    """

    def convert_bit_share(self, bit, dst_field):
        """Convert a 0/1 share into *dst_field*."""
        #k = self.options.security_parameter
        k = 29
        dst_share, src_share = self.prss_share_random_double_max(dst_field, bit.field, 2**k)
        #BS bit.field larger than log n * 2**k
        tmp = self.open(bit + src_share)
        tmp.addCallback(lambda i: dst_field(i.value))
        tmp.field = dst_field
        return tmp - dst_share

    def greater_than_equal(self, share_a, share_b):
        """Compute ``share_a >= share_b``.

        Both arguments must be shares from the same field. The result
        is a new 0/1 share from the field.
        """

        field = getattr(share_a, "field", getattr(share_b, "field", None))
#        if not isinstance(share_a, Share):
#            if not isinstance(share_a, FieldElement):
#                share_a = field(share_a)
#            share_a = Share(self, field, share_a)
#        if not isinstance(share_b, Share):
#            if not isinstance(share_b, FieldElement):
#                share_b = field(share_b)
#            share_b = Share(self, field, share_b)

        b = share_a - share_b # We will check if a == 0
        a = self.clone_share_nofp(b)

        smallField = GF(2147483647)  #Mersenne prime 2^31-1, bit length > 30.9759 old: GF(2703455651027L)
        
        l = self.options.bit_length + (self.options.res if a.fp or b.fp else 0)
        #k = self.options.security_parameter
        k = 29
        
        # TODO: verify asserts are correct...
        #assert field.modulus > 2**(l+2) + 2**(l+k), "Field too small"
        #assert smallField.modulus > 3 + 3*l, "smallField too small"
        #BS assert ook smallField.modulus  > k + log(n) .... but PRSS has many values added.
        
        r_bits = [self.prss_share_random(smallField, True) for _ in xrange(l)]
        r_bitsField = [self.convert_bit_share(bit, field) for bit in r_bits]
        r_modl = self.bin_comb(r_bitsField)
        r_divl = self.prss_share_random_max(field, 2**k)
        
        z_rmodl = a + 2**l + r_modl
        c = self.open(z_rmodl + 2**l*r_divl)
        self.schedule_callback(c, self._finish_greater_than_equal,
                               field, smallField, r_bits, z_rmodl)
        return c

    def _finish_greater_than_equal(self, c, field, smallField, r_bits, z_rmodl):
        l = len(r_bits)
        s_bit = self.prss_share_random(smallField, True)
        s_sign = 1 - 2 * s_bit
        # mask: uniformly random -- should be non-zero, failure prob. 1/2^k
        mask = self.prss_share_random(smallField, False) 
        #BS mask = mask * mask + 1, assuming -1 is in NQR. 
        
        E = [mask]
        sumXORs = 0
        for i in xrange(l-1, -1, -1):
            c_bit = c.bit(i)
            E.append(s_sign + r_bits[i] - c_bit + 3*sumXORs)
            sumXORs += r_bits[i] ^ c_bit
        E.append(s_sign - 1 + 3*sumXORs)
        
        E = self.fan_in_mulL(E, False)
            
        E = self.open(E)
        E.addCallback(lambda bit: field(bit.value != 0))
        UF = E ^ self.convert_bit_share(s_bit, field)
        z_rmodl = self.open(z_rmodl)
        # return  (z - (c%2**l - r%2**l + UF*2**l)) * 2^(-l)
        return (z_rmodl - (c.value%2**l + UF*2**l)) * ~field(2**l)

class Toft07Runtime(ComparisonToft07Mixin, PassiveRuntime):
    """Default mix of :class:`~viff.comparison.ComparisonToft07Mixin`
    and :class:`~viff.passive.PassiveRuntime`.
    """
    pass

class ComparisonConstantRoundMixin:

    """Efficient Constant Rounds Protocols. This mixin provides a
	
    - :meth:`greater_than_equal` method which can compare Zp field
    elements and gives a secret result shared over Zp.
	

    """


    def greater_than_equal(self, share_a, share_b):
        """Compute ``share_a >= share_b``.

        Both arguments must be shares from the same field. The result
        is a new 0/1 share from the field.
        """

        field = getattr(share_a, "field", getattr(share_b, "field", None))

        b = share_a - share_b
        a = self.clone_share_nofp(b)
		
        l = self.options.bit_length + self.options.res
        k = self.options.security_parameter
        ln = (int)(math.log(self.num_players,2)+1)
		
        assert field.modulus > 2**(l+1)+2**(l+k+ln+1), "field is too small"
		
        r_bits = [self.prss_share_random(field, True) for _ in xrange(l)]
        
        r_modl = self.bin_comb(r_bits)
        r_divl = self.prss_share_random_max(field, 2**k)
        
        z_rmodl = a + 2**l + r_modl
        c = self.open(z_rmodl + 2**l*r_divl)
        self.schedule_callback(c, self._finish_greater_than_equal, field, r_bits, z_rmodl)

        return c

    def _finish_greater_than_equal(self, c, field, r_bits, z_rmodl):
        l = self.options.bit_length

        plist = self.prefix_mulC([c.bit(l-1-i)+r_bits[l-1-i] -2*c.bit(l-1-i)*r_bits[l-1-i] + 1 for i in xrange(l)])
        slist = [plist[l-1-i]-plist[l-i-2] for i in xrange(l-1)]+[plist[0]-1]

        E = reduce(operator.add, map(operator.mul, map(lambda z: 1-z, [c.bit(i) for i in xrange(l)]), slist))
        UF = self.lsb(E)        

        return (z_rmodl - (c.value%2**l + UF*2**l)) * ~field(2**l)

class ConstantRoundRuntime(ComparisonConstantRoundMixin, PassiveRuntime):
    """Default mix of :class:`~viff.comparison.ConstantRoundMixin`
    and :class:`~viff.passive.PassiveRuntime`.
    """
    pass

