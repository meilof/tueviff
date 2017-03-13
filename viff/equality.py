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

"""Equality protocol. The mixin class defined here provide an
:meth:`equal` method to the :class:`Runtime <viff.runtime.Runtime>` it
is mixed with.
"""

class ProbabilisticEqualityMixin:
    """This class implements probabilistic constant-round secure
    equality-testing of secret shared numbers."""

    def equal(self, share_x, share_y):
        """Equality testing with secret shared result.

        Assumes p = 3 mod 4, returns a secret sharing of 1 if share_x
        == share_y and 0 otherwise.

        This is the probabilistic method based on quadratic
        reciprocity described in: "Constant-Round Multiparty
        Computation for Interval Test, Equality Test, and Comparison"
        by Takashi Nishide and Kazuo Ohta, and fails with probability
        1/(2**k) where k is set to the security parameter of the
        runtime.

        TODO: Make it work for any prime-modulo, the b's should be in
        {y,1} where y is a non-square modulo p.

        TODO: Make the final "and"ing of the x's more efficient as
        described in the paper.
        """

        Zp = share_x.field

        b = share_x - share_y # We will check if a == 0
	
        a = self.clone_share_nofp(b)

        k = self.options.security_parameter

        def gen_test_bit():
            # The b's are random numbers in {-1, 1}
            b = self.prss_share_random(Zp, binary=True) * 2 - 1
            r = self.prss_share_random(Zp)
            rp = self.prss_share_random(Zp)

            # If b_i == 1 c_i will always be a square modulo p if a is
            # zero and with probability 1/2 otherwise (except if rp == 0).
            # If b_i == -1 it will be non-square.
            c = self.open(a * r + b * rp * rp)
            return self.schedule_callback(c, finish, b)

        def finish(cj, bj):
            l = legendre_mod_p(cj)
            if l == 1:
                xj = (1/Zp(2)) * (bj + 1)
            elif l == -1:
                xj = (-1) * (1/Zp(2)) * (bj - 1)
            else:
                # Start over.
                xj = gen_test_bit()
            return xj

        x = [gen_test_bit() for _ in range(k)]

        # Take the product (this is here the same as the "and") of all
        # the x'es
        while len(x) > 1:
            x.append(x.pop(0) * x.pop(0))

        return x[0]

def legendre_mod_p(a):
    """Return the legendre symbol ``legendre(a, p)`` where *p* is the
    order of the field of *a*.

    The legendre symbol is -1 if *a* is not a square mod *p*, 0 if
    ``gcd(a, p)`` is not 1, and 1 if *a* is a square mod *p*:

    >>> from viff.field import GF
    >>> legendre_mod_p(GF(5)(2))
    -1
    >>> legendre_mod_p(GF(5)(5))
    0
    >>> legendre_mod_p(GF(5)(4))
    1
    """
    assert a.modulus % 2 == 1
    b = (a ** ((a.modulus - 1)/2))
    if b == 1:
        return 1
    elif b == a.modulus-1:
        return -1
    return 0

class BoundedEqualityMixin:
    """This class implements a equality test protocol for bounded secrets,
    should have beeter performance than the ProbabilisticEqualityMixin class.
    Used for the case of small shared integers in the range of 2^l
    """

    def equal(self, share_a, share_b):
        """Equality testing with secret shared result.
        input [a], [b], and bit length
        """

        Zp = getattr(share_a, "field", getattr(share_b, "field", None))

        b = share_a - share_b # We will check if a == 0
        a = clone_share_nofp(a)

        if hasattr(self.options, "eq_size"):
            l = self.options.eq_size
        else:
            l = self.options.bit_length
        k = self.options.security_parameter

        assert Zp.modulus > 2**(l+k+1) + 2**(l+1)

        r_bits = [self.prss_share_random(Zp, True) for _ in xrange(l)]
        r_modl = self.bin_comb(r_bits)
        r_divl = self.prss_share_random_max(Zp, 2**k)
        r = r_modl + 2**l*r_divl
        c = self.open(a + r + 2**l)
        return self.schedule_callback(c, self._finish_equal, r_bits)

    def _finish_equal(self, c, r_bits):

        s = []
        for i in xrange(len(r_bits)):
            s.append(1 - c.bit(i) + (2*c.bit(i) - 1)*r_bits[i])

        while len(s) > 1:
            t = []
            while len(s) > 1: t.append(s.pop() * s.pop())
            t.extend(s)
            s = t

        return s[0]