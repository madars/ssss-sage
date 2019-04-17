#!/usr/bin/env sage

#
# Auditable reimplementation of ssss.c in Sage.
# Author: Madars Virza <madars@mit.edu>
#

# Notes:

# - ssss.c uses binary extension fields and stores their elements in
# GMP's mpz_t data type, using the natural encoding.

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
        raise ValueError("Invalid extension degree (see field_size_valid in ssss.c)")
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
            raise Exception("Did not get a finite field of %d elements." % deg)

def cprng_init(fn):
    """Opens a file to read the randomness from and returns its handle. In
normal operation you'll want to set fn to /dev/random; but to a fixed
file for equivalence testing between shamir.sage and ssss.c.


This function corresponds to cprng_init of ssss.c."""
    return open(fn, "rb")

def mpz_import_impl(buf, F):
    """Implements mpz_import(x, degree / 8, 1, 1, 0, 0, buf); see
cprng_read() for details."""

    deg = F.degree()
    if len(buf) != deg/8:
        raise ValueError("This function requires buffer of the same size as a field element representation.")

    a = F.gen()
    el = F(0)
    for i in range(deg/8):
        # shift in the current byte
        el *= a^8
        for j in range(8):
            if ord(buf[i]) & (1<<j):
                el += a^j
    return el

def test_mpz_import_impl():
    F = get_field(16)

    # 0x4083 = 0100 0000 1000 0011
    #        = a^14 + a^7 + a + 1

    el = mpz_import_impl([chr(0x40), chr(0x83)], F)

    a = F.gen()
    el_expected = a^14 + a^7 + a + 1
    if el != el_expected:
        raise Exception("Simple manual mpz_import_impl test failed.")

def cprng_read(handle, F):
    """Reads a field element from a randomness source previously opened by
cprng_init. The principal call in ssss.c is

   mpz_import(x, degree / 8, 1, 1, 0, 0, buf)

Where arguments mean interpreting the degree / 8 byte buffer as an
array of 1 byte "words", stored with the most significant word first,
and not using GMP's nail bits.

(See https://gmplib.org/manual/Integer-Import-and-Export.html for details.)

This function corresponds to cprng_deinit of ssss.c"""
    deg = F.degree()

    # read deg/8 bytes
    buf = []
    for i in range(deg/8):
        char = handle.read(1)
        if char == "":
            raise Exception("Unexpected end of file.")
        buf.append(char)

    # convert them into a field element
    el = mpz_import_impl(buf, F)
    return el

        
if __name__ == '__main__':
    test_fields()
    test_mpz_import_impl()
