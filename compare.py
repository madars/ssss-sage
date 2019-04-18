#!/usr/bin/env python3

import random
import string
import os
from subprocess import Popen, PIPE, STDOUT

NUM_SHARES=1000
THRESHOLD=400
MIN_BITS=16
MAX_BITS=1024

def pick_test_secret(n):
    """Generate a test secret consisting of n ASCII characters."""
    return ''.join(random.choices(string.ascii_letters +
                                  string.digits, k=n))

def get_ssss_split_shares(secret, threshold, num_shares):
    """Invoke a ssss-spit -t THRESHOLD -n NUM_SHARES -D, inputting the
specified secret, and return the list of shares output by it."""
    p = Popen(['./ssss-split',
               '-t', str(THRESHOLD),
               '-n', str(NUM_SHARES),
               '-D'],
              stdout=PIPE, stdin=PIPE, stderr=PIPE)
    stdout = p.communicate(input=bytes(secret, "ascii")+b"\n")[0]
    # skip the first line which is purely diagnostic
    shares = stdout.split(b"\n")[1:]
    return shares

def get_sage_split_shares(secret, threshold, num_shares):
    """Invoke a shamir.sage -t THRESHOLD -n NUM_SHARES, inputting the
specified secret, and return the list of shares output by it."""
    p = Popen(['./shamir.sage',
               '-t', str(THRESHOLD),
               '-n', str(NUM_SHARES),
               'split'],
              stdout=PIPE, stdin=PIPE, stderr=PIPE)
    stdout = p.communicate(input=bytes(secret, "ascii")+b"\n")[0]
    shares = stdout.split(b"\n")
    return shares

if __name__ == '__main__':
    num_random_bytes = MAX_BITS * (THRESHOLD - 1) // 8
    with open("fixed-random", "wb") as f:
        f.write(os.urandom(num_random_bytes))

    all_pass = True

    for num_secret_bytes in range(MIN_BITS//8, (MAX_BITS//8) + 1):
        print("Testing secret length {}...".format(num_secret_bytes))
        secret = pick_test_secret(num_secret_bytes)
        ssss_shares = get_ssss_split_shares(secret, THRESHOLD, NUM_SHARES)
        sage_shares = get_sage_split_shares(secret, THRESHOLD, NUM_SHARES)
        if ssss_shares != sage_shares:
            all_pass = False
            print(("Did not get the same shares for secret='{secret:s}', " +
                   "t={threshold:d}, n={number:d}.").format(
                       secret=secret,
                       threshold=THRESHOLD,
                       number=NUM_SHARES))
    if all_pass:
        print("All comparisons passed!")
    else:
        print("Implementations did not return the same results.")
