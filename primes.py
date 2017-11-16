#!/usr/bin/env python
#import psyco; psyco.full()
from math import sqrt, ceil
#import numpy as np

def rwh_primes(n):
    # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    """ Returns  a list of primes < n """
    sieve = [True] * n
    for i in xrange(3,int(n**0.5)+1,2):
        if sieve[i]:
            sieve[i*i::2*i]=[False]*((n-i*i-1)/(2*i)+1)
    return [2] + [i for i in xrange(3,n,2) if sieve[i]]

def rwh_primes1(n):
    # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    """ Returns  a list of primes < n """
    sieve = [True] * (n/2)
    for i in xrange(3,int(n**0.5)+1,2):
        if sieve[i/2]:
            sieve[i*i/2::i] = [False] * ((n-i*i-1)/(2*i)+1)
    return [2] + [2*i+1 for i in xrange(1,n/2) if sieve[i]]

def rwh_primes3(n):
    # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    """ Returns  a list of primes < n """
    sieve = [True] * (n/2)
    for i in xrange(3,int(n**0.5)+1,2):
        if sieve[i/2]:
            sieve[i*i/2::i] = [False] * ((n-i*i-1)/(2*i)+1)
    all = [2]

    for i in xrange(1, n / 2):
        if sieve[i]:
            all.append(2*i+1)
            print 'p3:{}\n'.format(2*i+1)

    print 'count:{}\n'.format(len(all))
    return all

def rwh_primes2(n):
    # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    """ Input n>=6, Returns a list of primes, 2 <= p < n """
    correction = (n%6>1)
    n = {0:n,1:n-1,2:n+4,3:n+3,4:n+2,5:n+1}[n%6]
    sieve = [True] * (n/3)
    sieve[0] = False
    for i in xrange(int(n**0.5)/3+1):
      if sieve[i]:
        k=3*i+1|1
        sieve[      ((k*k)/3)      ::2*k]=[False]*((n/6-(k*k)/6-1)/k+1)
        sieve[(k*k+4*k-2*k*(i&1))/3::2*k]=[False]*((n/6-(k*k+4*k-2*k*(i&1))/6-1)/k+1)
    return [2,3] + [3*i+1|1 for i in xrange(1,n/3-correction) if sieve[i]]

def sieve_wheel_30(N):
    # http://zerovolt.com/?p=88
    ''' Returns a list of primes <= N using wheel criterion 2*3*5 = 30

Copyright 2009 by zerovolt.com
This code is free for non-commercial purposes, in which case you can just leave this comment as a credit for my work.
If you need this code for commercial purposes, please contact me by sending an email to: info [at] zerovolt [dot] com.'''
    __smallp = ( 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
    61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
    149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
    313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
    409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
    499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
    601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683,
    691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
    809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887,
    907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997)

    wheel = (2, 3, 5)
    const = 30
    if N < 2:
        return []
    if N <= const:
        pos = 0
        while __smallp[pos] <= N:
            pos += 1
        return list(__smallp[:pos])
    # make the offsets list
    offsets = (7, 11, 13, 17, 19, 23, 29, 1)
    # prepare the list
    p = [2, 3, 5]
    dim = 2 + N // const
    tk1  = [True] * dim
    tk7  = [True] * dim
    tk11 = [True] * dim
    tk13 = [True] * dim
    tk17 = [True] * dim
    tk19 = [True] * dim
    tk23 = [True] * dim
    tk29 = [True] * dim
    tk1[0] = False
    # help dictionary d
    # d[a , b] = c  ==> if I want to find the smallest useful multiple of (30*pos)+a
    # on tkc, then I need the index given by the product of [(30*pos)+a][(30*pos)+b]
    # in general. If b < a, I need [(30*pos)+a][(30*(pos+1))+b]
    d = {}
    for x in offsets:
        for y in offsets:
            res = (x*y) % const
            if res in offsets:
                d[(x, res)] = y
    # another help dictionary: gives tkx calling tmptk[x]
    tmptk = {1:tk1, 7:tk7, 11:tk11, 13:tk13, 17:tk17, 19:tk19, 23:tk23, 29:tk29}
    pos, prime, lastadded, stop = 0, 0, 0, int(ceil(sqrt(N)))
    # inner functions definition
    def del_mult(tk, start, step):
        for k in xrange(start, len(tk), step):
            tk[k] = False
    # end of inner functions definition
    cpos = const * pos
    while prime < stop:
        # 30k + 7
        if tk7[pos]:
            prime = cpos + 7
            p.append(prime)
            lastadded = 7
            for off in offsets:
                tmp = d[(7, off)]
                start = (pos + prime) if off == 7 else (prime * (const * (pos + 1 if tmp < 7 else 0) + tmp) )//const
                del_mult(tmptk[off], start, prime)
        # 30k + 11
        if tk11[pos]:
            prime = cpos + 11
            p.append(prime)
            lastadded = 11
            for off in offsets:
                tmp = d[(11, off)]
                start = (pos + prime) if off == 11 else (prime * (const * (pos + 1 if tmp < 11 else 0) + tmp) )//const
                del_mult(tmptk[off], start, prime)
        # 30k + 13
        if tk13[pos]:
            prime = cpos + 13
            p.append(prime)
            lastadded = 13
            for off in offsets:
                tmp = d[(13, off)]
                start = (pos + prime) if off == 13 else (prime * (const * (pos + 1 if tmp < 13 else 0) + tmp) )//const
                del_mult(tmptk[off], start, prime)
        # 30k + 17
        if tk17[pos]:
            prime = cpos + 17
            p.append(prime)
            lastadded = 17
            for off in offsets:
                tmp = d[(17, off)]
                start = (pos + prime) if off == 17 else (prime * (const * (pos + 1 if tmp < 17 else 0) + tmp) )//const
                del_mult(tmptk[off], start, prime)
        # 30k + 19
        if tk19[pos]:
            prime = cpos + 19
            p.append(prime)
            lastadded = 19
            for off in offsets:
                tmp = d[(19, off)]
                start = (pos + prime) if off == 19 else (prime * (const * (pos + 1 if tmp < 19 else 0) + tmp) )//const
                del_mult(tmptk[off], start, prime)
        # 30k + 23
        if tk23[pos]:
            prime = cpos + 23
            p.append(prime)
            lastadded = 23
            for off in offsets:
                tmp = d[(23, off)]
                start = (pos + prime) if off == 23 else (prime * (const * (pos + 1 if tmp < 23 else 0) + tmp) )//const
                del_mult(tmptk[off], start, prime)
        # 30k + 29
        if tk29[pos]:
            prime = cpos + 29
            p.append(prime)
            lastadded = 29
            for off in offsets:
                tmp = d[(29, off)]
                start = (pos + prime) if off == 29 else (prime * (const * (pos + 1 if tmp < 29 else 0) + tmp) )//const
                del_mult(tmptk[off], start, prime)
        # now we go back to top tk1, so we need to increase pos by 1
        pos += 1
        cpos = const * pos
        # 30k + 1
        if tk1[pos]:
            prime = cpos + 1
            p.append(prime)
            lastadded = 1
            for off in offsets:
                tmp = d[(1, off)]
                start = (pos + prime) if off == 1 else (prime * (const * pos + tmp) )//const
                del_mult(tmptk[off], start, prime)
    # time to add remaining primes
    # if lastadded == 1, remove last element and start adding them from tk1
    # this way we don't need an "if" within the last while
    if lastadded == 1:
        p.pop()
    # now complete for every other possible prime
    while pos < len(tk1):
        cpos = const * pos
        if tk1[pos]: p.append(cpos + 1)
        if tk7[pos]: p.append(cpos + 7)
        if tk11[pos]: p.append(cpos + 11)
        if tk13[pos]: p.append(cpos + 13)
        if tk17[pos]: p.append(cpos + 17)
        if tk19[pos]: p.append(cpos + 19)
        if tk23[pos]: p.append(cpos + 23)
        if tk29[pos]: p.append(cpos + 29)
        pos += 1
    # remove exceeding if present
    pos = len(p) - 1
    while p[pos] > N:
        pos -= 1
    if pos < len(p) - 1:
        del p[pos+1:]
    # return p list
    return p

def sieveOfEratosthenes(n):
    """sieveOfEratosthenes(n): return the list of the primes < n."""
    # Code from: <dickinsm@gmail.com>, Nov 30 2006
    # http://groups.google.com/group/comp.lang.python/msg/f1f10ced88c68c2d
    if n <= 2:
        return []
    sieve = range(3, n, 2)
    top = len(sieve)
    for si in sieve:
        if si:
            bottom = (si*si - 3) // 2
            if bottom >= top:
                break
            sieve[bottom::si] = [0] * -((bottom - top) // si)
    return [2] + [el for el in sieve if el]

def sieveOfAtkin(end):
    """sieveOfAtkin(end): return a list of all the prime numbers <end
    using the Sieve of Atkin."""
    # Code by Steve Krenzel, <Sgk284@gmail.com>, improved
    # Code: https://web.archive.org/web/20080324064651/http://krenzel.info/?p=83
    # Info: http://en.wikipedia.org/wiki/Sieve_of_Atkin
    assert end > 0
    lng = ((end-1) // 2)
    sieve = [False] * (lng + 1)

    x_max, x2, xd = int(sqrt((end-1)/4.0)), 0, 4
    for xd in xrange(4, 8*x_max + 2, 8):
        x2 += xd
        y_max = int(sqrt(end-x2))
        n, n_diff = x2 + y_max*y_max, (y_max << 1) - 1
        if not (n & 1):
            n -= n_diff
            n_diff -= 2
        for d in xrange((n_diff - 1) << 1, -1, -8):
            m = n % 12
            if m == 1 or m == 5:
                m = n >> 1
                sieve[m] = not sieve[m]
            n -= d

    x_max, x2, xd = int(sqrt((end-1) / 3.0)), 0, 3
    for xd in xrange(3, 6 * x_max + 2, 6):
        x2 += xd
        y_max = int(sqrt(end-x2))
        n, n_diff = x2 + y_max*y_max, (y_max << 1) - 1
        if not(n & 1):
            n -= n_diff
            n_diff -= 2
        for d in xrange((n_diff - 1) << 1, -1, -8):
            if n % 12 == 7:
                m = n >> 1
                sieve[m] = not sieve[m]
            n -= d

    x_max, y_min, x2, xd = int((2 + sqrt(4-8*(1-end)))/4), -1, 0, 3
    for x in xrange(1, x_max + 1):
        x2 += xd
        xd += 6
        if x2 >= end: y_min = (((int(ceil(sqrt(x2 - end))) - 1) << 1) - 2) << 1
        n, n_diff = ((x*x + x) << 1) - 1, (((x-1) << 1) - 2) << 1
        for d in xrange(n_diff, y_min, -8):
            if n % 12 == 11:
                m = n >> 1
                sieve[m] = not sieve[m]
            n += d

    primes = [2, 3]
    if end <= 3:
        return primes[:max(0,end-2)]

    for n in xrange(5 >> 1, (int(sqrt(end))+1) >> 1):
        if sieve[n]:
            primes.append((n << 1) + 1)
            aux = (n << 1) + 1
            aux *= aux
            for k in xrange(aux, end, 2 * aux):
                sieve[k >> 1] = False

    s  = int(sqrt(end)) + 1
    if s  % 2 == 0:
        s += 1
    primes.extend([i for i in xrange(s, end, 2) if sieve[i >> 1]])

    return primes

def ambi_sieve_plain(n):
    s = range(3, n, 2)
    for m in xrange(3, int(n**0.5)+1, 2):
        if s[(m-3)/2]:
            for t in xrange((m*m-3)/2,(n>>1)-1,m):
                s[t]=0
    return [2]+[t for t in s if t>0]

def sundaram3(max_n):
    # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/2073279#2073279
    numbers = range(3, max_n+1, 2)
    half = (max_n)//2
    initial = 4

    for step in xrange(3, max_n+1, 2):
        for i in xrange(initial, half, step):
            numbers[i-1] = 0
        initial += 2*(step+1)

        if initial > half:
            return [2] + filter(None, numbers)


'''
from primes import primecc
from primes import rwh_primes3
primecc(5105100)
rwh_primes3(5105100)
'''
def primecc(n):
    #lprimes = [2]
    #cprimes = [2,3]
    #nprimes = [2, 3, 5]
    #np_c_primes = [5]
    all = [2, 3, 5]

    #lp = 2
    #cp = 3
    #np = 5

    lcc = 2
    cc = 6 # 2 * 3
    ncc = 30 # 2 * 3 * 5
    c_index = 2 # lprimes = all[0:c_index-1], cprimes = all[0:c_index], nprimes = all[0:c_index + 1], np_c_primes = all[c_index:]
    done = False
    while cc < n:
        for i in xrange(cc + 1, ncc, 2):
            if i > n:
                done = True
                break
            found = False
            for j in all[0:c_index]:
                if i % j == 0:
                    found = True
                    break

            for j in all[c_index:]:
                if i % j == 0:
                    found = True
                    break
                #if j > lcc:
                 #   break
                if i / j < j:
                    break
            if not found:
                print 'cc:{}\n'.format(i)
                all.append(i)
        if done:
            break
        lcc = cc
        cc = ncc
        #cp = np
        c_index = c_index + 1
        np = all[c_index]
        ncc = cc * np

    return all

def get_multiple_set(array, min, max ):
    d = {}
    for i in array:
        d[i] = 0
    #for i in array:

'''
from primes import cr
list(cr([11, 13, 17, 19, 23, 29],210))
from primes import rwh_primes3
'''
import math
def cr(iterable, max):
    # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
    pool = tuple(iterable)

    for i in pool:
        yield i
    i = 0
    s = pool[0]
    while s < max:
        s = s * pool[0]
        i = i + 1

    r = i

    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    sum = 1
    for i in indices:
        sum = sum * pool[i]
    if sum < max:
        yield sum
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        sum = 1
        for i in indices:
            sum = sum * pool[i]
        if sum < max:
            yield sum
        #yield tuple(pool[i] for i in indices)

'''
from primes import primepcn
primepcn(1)
from primes import rwh_primes3
rwh_primes3(510510)
from primes import sieveOfAtkin
sieveOfAtkin(510510)
'''
def primepcn(x): # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
    c = 1
    c_n = 2
    n_c_l = [1]
    p_l = []
    p_n = 2
    for n in xrange(1,x+1):
        l_n_c_l = len(n_c_l)

        i = 1
        j = 1
        s = []
        for i in n_c_l[2:]:
            for j in n_c_l[1:i-1]:
                t = i * j
                if t > c and t < c_n:
                    s.append(t)
                if t > c_n:
                    break


        s.sort()
        #print s
        si = 0
        for k in xrange(1, p_n):
            for i in n_c_l:
                if i > c:
                    break
                t = c * k + i
                if t % p_n != 0:
                    n_c_l.append(t)
                elif t != p_n:
                    continue

                while si < len(s) and s[si] < t:
                    si = si + 1
                if len(s) > si:
                    if t != s[si]:
                        p_l.append(t)
                    else:
                        si = si + 1
                else:
                    p_l.append(t)
                    #print 'pcn:{}\n'.format(t)

        for i in n_c_l:
             if i % p_n == 0:
                 n_c_l.remove(i)

        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]


        c = c_n
        c_n = c * p_n

    print 'count:{}\n'.format(len(p_l))
    return p_l


def primepcn9(x): # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
    c = 1
    c_n = 2
    n_c_l = [1]
    p_l = []
    p_n = 2
    for n in xrange(1,x+1):
        l_n_c_l = len(n_c_l)

        i = 1
        j = 1
        s = []
        for i in n_c_l[2:]:
            for j in n_c_l[1:i-1]:
                t = i * j
                if t > c and t < c_n:
                    s.append(t)

        s.sort()
        #print s
        si = 0
        for k in xrange(1, p_n):
            for i in n_c_l:
                if i > c:
                    break
                t = c * k + i
                if t % p_n != 0:
                    n_c_l.append(t)
                elif t != p_n:
                    continue

                while si < len(s) and s[si] < t:
                    si = si + 1
                if len(s) > si:
                    if t != s[si]:
                        p_l.append(t)
                    else:
                        si = si + 1
                else:
                    p_l.append(t)
                    #print 'pcn:{}\n'.format(t)

        for i in n_c_l:
             if i % p_n == 0:
                 n_c_l.remove(i)

        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]


        c = c_n
        c_n = c * p_n

    print 'count:{}\n'.format(len(p_l))
    return p_l


def primepcn8(x): # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
    c = 1
    c_n = 2
    n_c_l = [1]
    p_l = []
    p_n = 2
    for n in xrange(1,x+1):
        l_n_c_l = len(n_c_l)

        i = 1
        j = 1
        s = []
        for i in n_c_l[1:]:
            for j in n_c_l[1:i-1]:
                t = i * j
                if t > c and t < c_n:
                    s.append(t)

        s.sort()
        print s
        si = 0
        for k in xrange(1, p_n):
            for i in n_c_l:
                if i > c:
                    break
                t = c * k + i
                if t % p_n != 0:
                    n_c_l.append(t)
                elif t != p_n:
                    continue

                #found = False

                # for j in n_c_l[1:l_n_c_l]:
                #     if t % j == 0:
                #         found = True
                #         break
                #if not found:
                while si < len(s) and s[si] < t:
                    si = si + 1
                if len(s) > si:
                    if t != s[si]:
                        p_l.append(t)
                    else:
                        si = si + 1
                else:
                    p_l.append(t)
                #if t not in s:
                    #p_l.append(t)
                    #print 'pcn:{}\n'.format(t)

        for i in n_c_l:
             if i % p_n == 0:
                 n_c_l.remove(i)

        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]


        c = c_n
        c_n = c * p_n

    print 'count:{}\n'.format(len(p_l))
    return p_l


def primepcn7(x): # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
    c = 1
    c_n = 2
    n_c_l = [1]
    p_l = []
    p_n = 2
    for n in xrange(1,x+1):
        l_n_c_l = len(n_c_l)
        for k in xrange(1, p_n):
            for i in n_c_l:
                if i > c:
                    break
                t = c * k + i
                if t % p_n != 0:
                    n_c_l.append(t)
                elif t != p_n:
                    continue

                found = False

                for j in n_c_l[1:l_n_c_l]:
                    if t % j == 0:
                        found = True
                        break
                if not found:
                    p_l.append(t)
                    #print 'pcn:{}\n'.format(t)

        for i in n_c_l:
             if i % p_n == 0:
                 n_c_l.remove(i)



        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]


        c = c_n
        c_n = c * p_n

    print 'count:{}\n'.format(len(p_l))
    return p_l


def primepcn6(x): # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
    c = 1
    c_n = 2
    n_c_l = [1]
    p_l = []
    p_n = 2
    p_j_mark = 1
    for n in xrange(1,x+1):
        l_n_c_l = len(n_c_l)
        for k in xrange(1, p_n):
            for i in n_c_l:
                if i > c:
                    break
                t = c * k + i
                if t % p_n != 0:
                    n_c_l.append(t)
                elif t != p_n:
                    continue

                found = False

                for j in n_c_l[p_j_mark:l_n_c_l]:
                    if t % j == 0:
                        found = True
                        break
                if not found:
                    p_l.append(t)
                    #print 'pcn:{}\n'.format(t)

        for i in n_c_l:
             if i % p_n == 0:
                 n_c_l.remove(i)

        for i in xrange(1,len(n_c_l) - 1):
            if n_c_l[i]  > p_n:
                p_j_mark = i
                break


        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]


        c = c_n
        c_n = c * p_n

    print 'count:{}\n'.format(len(p_l))
    return p_l

def primepcn5(x): # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
    c = 1
    c_n = 2
    n_c_l = [1]
    p_l = []
    #p = 1
    p_n = 2
    p_n_mark = 0
    p_j_mark = 1
    p_j_n_mark = 0
    for n in xrange(1,x+1):
        for k in xrange(1, p_n):
            for i in n_c_l:
                if i > c:
                    break
                t = c * k + i
                if t % p_n != 0:
                    n_c_l.append(t)
                elif t != p_n:
                    continue

                found = False
                t_mark = p_j_n_mark
                if t_mark == 0:
                    if n_c_l[-1] == t:
                        t_mark = len(n_c_l) - 1
                    else:
                        t_mark = len(n_c_l)
                for j in n_c_l[p_j_mark:t_mark]:
                    if t % j == 0:
                        found = True
                        break
                if not found:
                    p_l.append(t)
                    print 'pcn:{}\n'.format(t)

        # if p_n in n_c_l:
        #    n_c_l.remove(p_n)
        # for i in n_c_l:
        #     if i % p_n == 0:
        #         n_c_l.remove(i)

        # for i in n_c_l[:p_n_mark]:
        #     if i * p_n in n_c_l:
        #         n_c_l.remove(i * p_n)
        if p_n in n_c_l:
            n_c_l.remove(p_n)

        found_p_j_mark = False

        for i in xrange(1,len(n_c_l) - 1):
            if not found_p_j_mark and n_c_l[i]  > p_n:
                p_j_mark = i
                found_p_j_mark = True
            if n_c_l[i] >= c_n:
                p_j_n_mark = i
                break


        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]


        c = c_n
        c_n = c * p_n

        # for i in xrange(len(n_c_l) - 1,0,-1):
        #     if n_c_l[i] * p_n < c_n:
        #          p_n_mark = i
        #          break

        #print 'p_n_mark:{}\n'.format(p_n_mark)
    print 'count:{}\n'.format(len(p_l))
    return p_l

def primepcn4(x): # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
    c = 1
    c_n = 2
    n_c_l = [1]
    p_l = []
    #p = 1
    p_n = 2
    p_n_mark = 0
    for n in xrange(1,x+1):
        for k in xrange(1, p_n):
            for i in n_c_l:
                if i > c:
                    break;
                t = c * k + i
                if t % p_n != 0:
                    n_c_l.append(t)

                found = False
                for j in n_c_l[1:]:
                    if j > c:
                        break
                    if t % j == 0:
                        found = True
                        break
                if not found:
                    p_l.append(t)
                    print 'pcn:{}\n'.format(t)


        for i in n_c_l[:p_n_mark]:
                n_c_l.remove(i * p_n)

        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]

        c = c_n
        c_n = c * p_n

        for i in xrange(len(n_c_l) - 1,0,-1):
            if n_c_l[i] * p_n < c_n:
                p_n_mark = i

    return p_l

def primepcn3(x): # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
    c = 1
    c_n = 2
    n_c_l = [1]
    p_l = []
    #p = 1
    p_n = 2
    for n in xrange(1,x+1):
        for k in xrange(1, p_n):
            for i in n_c_l:
                if i > c:
                    break;
                t = c * k + i
                n_c_l.append(t)
                found = False
                for j in n_c_l[1:]:
                    if j > c:
                        break
                    if t % j == 0:
                        found = True
                        break
                if not found:
                    p_l.append(t)
                    print 'pcn:{}\n'.format(t)


        for i in n_c_l:
            if i % p_n == 0:
                n_c_l.remove(i)

        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]
        c = c_n
        c_n = c * p_n

    return p_l

def primepcn2(x):  # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
    all = [2, 3, 5]
    s_p = [1]
    c_p = 2
    c = 6
    p = 3
    p_n = 5

    for n in xrange(3, x + 1):
        s = []
        # len_p = len(all)
        # s(p(n)) = {C(n-1)*j + s(p(n-1))| 0=<j<p(n) exp {p(n) * k|k = {N}}
        for j in xrange(0, p):
            for i in s_p:
                t = c_p * j + i
                if (t % p) > 0:
                    s.append(t)

        # f(p(n)) = {C(n) * k + s(p(n)) | 0 <= k < p(n + 1)}
        # SP(n) = f(p(n)) exp {*p(k)| p(n)<k<C(n) and  C(n) <*p(k) < C(n+1)} = f(p(n)) exp e(n)
        for k in xrange(1, p_n):
            for i in s:
                t = c * k + i
                found = False
                for j in xrange(n - 1, len(all)):
                    if all[j] > c:
                        break;
                    if t % all[j] == 0:
                        found = True
                        break
                if not found:
                    all.append(t)
                    print 'pcn:{}\n'.format(t)

        s_p = s
        p = all[n - 1]
        c_p = c
        c = c * p
        p_n = all[n]

    return all

    '''
    from primes import primepcn
    primepcn(2)
    '''

    def primepcn_old(x):  # x is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30}
        all = [1, 2, 3, 5]
        de = {5: [25]}
        s_p = [1]
        c_p = 2
        c = 6
        p = 3
        p_n = 5
        for n in xrange(3, x):
            s = []
            # len_p = len(all)
            # s(p(n)) = {C(n-1)*j + s(p(n-1))| 0=<j<p(n) exp {p(n) * k|k = {N}}
            for j in xrange(0, p):
                for i in s_p:
                    t = c_p * j + i
                    if (t % p) > 0:
                        s.append(t)

            # f(p(n)) = {C(n) * k + s(p(n)) | 0 <= k < p(n + 1)}
            # SP(n) = f(p(n)) exp {*p(k)| p(n)<k<C(n) and  C(n) <*p(k) < C(n+1)} = f(p(n)) exp e(n)
            for k in xrange(1, p_n):
                for i in s:
                    t = c * k + i

                    found = False
                    for u in de:
                        if u > p_n:
                            break;
                        for v in de[u]:
                            if v == t:
                                found = True
                                break
                        if found:
                            break

                    if not found:
                        all.append(t)

                    for j in xrange(n + 1, len(all)):
                        aj = all[j]
                        cj = 1
                        for m in all[1:j]:
                            cj = cj * m
                        if t < cj:
                            if aj not in de:
                                de[aj] = []
                            de[aj].append(aj * t)
                        else:
                            break
                        if aj in de:
                            for v in de[aj]:
                                if v * t < cj * aj:
                                    de[aj].append(v * t)

            for i in de.keys():
                if i % p_n == 0:
                    de.pop(i, None)
            s_p = s
            p = all[n]
            c_p = c
            c = c * p
            p_n = all[n + 1]

        return all








        ################################################################################
# Using Numpy:
# def ambi_sieve(n):
#     # http://tommih.blogspot.com/2009/04/fast-prime-number-generator.html
#     s = np.arange(3, n, 2)
#     for m in xrange(3, int(n ** 0.5)+1, 2):
#         if s[(m-3)/2]:
#             s[(m*m-3)/2::m]=0
#     return np.r_[2, s[s>0]]
#
# def primesfrom3to(n):
#     # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
#     """ Returns a array of primes, p < n """
#     assert n>=2
#     sieve = np.ones(n/2, dtype=np.bool)
#     for i in xrange(3,int(n**0.5)+1,2):
#         if sieve[i/2]:
#             sieve[i*i/2::i] = False
#     return np.r_[2, 2*np.nonzero(sieve)[0][1::]+1]
#
# def primesfrom2to(n):
#     # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
#     """ Input n>=6, Returns a array of primes, 2 <= p < n """
#     sieve = np.ones(n/3 + (n%6==2), dtype=np.bool)
#     sieve[0] = False
#     for i in xrange(int(n**0.5)/3+1):
#         if sieve[i]:
#             k=3*i+1|1
#             sieve[      ((k*k)/3)      ::2*k] = False
#             sieve[(k*k+4*k-2*k*(i&1))/3::2*k] = False
#     return np.r_[2,3,((3*np.nonzero(sieve)[0]+1)|1)]



if __name__=='__main__':
    import itertools
    import sys

    def test(f1,f2,num):
        a =[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
         109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
         233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,
         367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
         499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641,
         643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787,
         797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
         947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
         1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213,
         1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321,
         1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481,
         1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601,
         1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
         1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877,
         1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017,
         2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143,
         2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297,
         2309]
        s = 1
        for i in xrange(0,num):
            s = s * a[i]

        print('Testing {f1} and {f2} return same results'.format(
            f1=f1.func_name,
            f2=f2.func_name))
        if not all([a==b for a,b in itertools.izip_longest(f1(s),f2(num))]):
            sys.exit("Error: %s(%s) != %s(%s)"%(f1.func_name,num,f2.func_name,s))

    n=6
    a = primepcn(n)
    print a
    print len(a)
    #test(sieveOfAtkin,primepcn,n)
    #test(sieveOfAtkin,ambi_sieve,n)
    #test(sieveOfAtkin,ambi_sieve_plain,n)
    #test(sieveOfAtkin,sundaram3,n)
    #test(sieveOfAtkin,sieve_wheel_30,n)
    #test(sieveOfAtkin,primesfrom3to,n)
    #test(sieveOfAtkin,primesfrom2to,n)
    #test(sieveOfAtkin,rwh_primes,n)
    #test(sieveOfAtkin,rwh_primes1,n)
    #test(sieveOfAtkin,rwh_primes2,n)
    #test(primecc,rwh_primes3, n)