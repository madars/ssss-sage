#!/usr/bin/env sage

#
# Auditable reimplementation of ssss.c in Sage.
# Author: Madars Virza <madars@mit.edu>
#

# constants copied from ssss.c
max_degree = 1024
irred_coeff = [
  4,3,1,5,3,1,4,3,1,7,3,2,5,4,3,5,3,2,7,4,2,4,3,1,10,9,3,9,4,2,7,6,2,10,9,
  6,4,3,1,5,4,3,4,3,1,7,2,1,5,3,2,7,4,2,6,3,2,5,3,2,15,3,2,11,3,2,9,8,7,7,
  2,1,5,3,2,9,3,1,7,3,1,9,8,3,9,4,2,8,5,3,15,14,10,10,5,2,9,6,2,9,3,2,9,5,
  2,11,10,1,7,3,2,11,2,1,9,7,4,4,3,1,8,3,1,7,4,1,7,2,1,13,11,6,5,3,2,7,3,2,
  8,7,5,12,3,2,13,10,6,5,3,2,5,3,2,9,5,2,9,7,2,13,4,3,4,3,1,11,6,4,18,9,6,
  19,18,13,11,3,2,15,9,6,4,3,1,16,5,2,15,14,6,8,5,2,15,11,2,11,6,2,7,5,3,8,
  3,1,19,16,9,11,9,6,15,7,6,13,4,3,14,13,3,13,6,3,9,5,2,19,13,6,19,10,3,11,
  6,5,9,2,1,14,3,2,13,3,1,7,5,4,11,9,8,11,6,5,23,16,9,19,14,6,23,10,2,8,3,
  2,5,4,3,9,6,4,4,3,2,13,8,6,13,11,1,13,10,3,11,6,5,19,17,4,15,14,7,13,9,6,
  9,7,3,9,7,1,14,3,2,11,8,2,11,6,4,13,5,2,11,5,1,11,4,1,19,10,3,21,10,6,13,
  3,1,15,7,5,19,18,10,7,5,3,12,7,2,7,5,1,14,9,6,10,3,2,15,13,12,12,11,9,16,
  9,7,12,9,3,9,5,2,17,10,6,24,9,3,17,15,13,5,4,3,19,17,8,15,6,3,19,6,1 ]

# fresh code below

def get_field(deg):
    """Return a binary field GF(2^deg) from a hard-coded list of polynomials.

This roughly corresponds to field_init function of ssss.c
"""
    if not (deg >= 8 and deg <= max_degree and deg % 8 == 0):
        raise "Invalid extension degree (see field_size_valid in ssss.c)"
    GF2.<x> = GF(2)[]

    c1 = irred_coeff[3 * (deg / 8 - 1) + 0]
    c2 = irred_coeff[3 * (deg / 8 - 1) + 1]
    c3 = irred_coeff[3 * (deg / 8 - 1) + 2]

    poly = x^deg + x^c1 + x^c2 + x^c3 + 1
    # done anyway inside constructor below
    assert poly.is_irreducible()
    GFext.<a> = GF(2^deg, modulus=poly)

    return GFext

def test_fields():
    for deg in xrange(8, max_degree+8, 8):
        F = get_field(deg)
        if F.degree() != deg:
            raise "Did not get a finite field of %d elements." % deg

if __name__ == '__main__':
    test_fields()
