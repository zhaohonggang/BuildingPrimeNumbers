'''
Theory:

P - Prime Numbers: P = {p(n)|n>=1,p(n) is the n-th prime number}
C - Prefect Composite Number: C = {c(p(n))| c(p(n)) = *p(i)|1<=i<=p(n)}
NC - Numbers co-prime with C: NC(c(p(n))) = {nc(c(p(n)))| nc(c(p(n))) is co-prime with c(p(n))}

NC{ c(p(n))..c(p(n+1)) } - Co-prime List. They are all possible numbers between c(p(n))and c(p(n+1)) which is co-prime with c(p(n)):

    NC{ c(p(n))..c(p(n+1)) } = {nc(c(p(n)))| c(p(n)) < nc(c(p(n))) < c(p(n+1)), nc(c(p(n))) is co-prime with c(p(n))}
                             = {c(p(n)) * k + s(c(p(n)))| 0<k<p(n+1), 0<s(c(p(n)))<c(p(n)), s(c(p(n))) is co-prime with c(p(n)) }
                             = {c(p(n)) * k + s(c(p(n)))| 0<k<p(n+1), 0<s(c(p(n)))<c(p(n)), s(c(p(n))) = *p(i)|p(n+1)<=p(i)<c(p(n)) }

NCC{ c(p(n))..c(p(n+1)) } - Except List. They are all possible composite numbers between c(p(n))and c(p(n+1)) which is co-prime with c(p(n)):

    NCC{ c(p(n))..c(p(n+1)) } = {ncc(c(p(n)))| c(p(n)) < ncc(c(p(n))) < c(p(n+1)), ncc(c(p(n))) is co-prime with c(p(n)), ncc(c(p(n))) is composite number }
                              = {ncc(c(p(n)))| c(p(n)) < ncc(c(p(n))) < c(p(n+1)), ncc(c(p(n))) = *p(i)|p(n+1)<=p(i)<c(p(n))}

Finally,
P{ c(p(n))..c(p(n+1)) } - Primes List. They are all possible prime numbers between c(p(n))and c(p(n+1)) which is co-prime with c(p(n)):

    P{ c(p(n))..c(p(n+1)) } = NC{ c(p(n))..c(p(n+1)) } - NCC{ c(p(n))..c(p(n+1)) }

That's all
'''
from math import log
'''
from primepcn import getPrimes
a = getPrimes(6)
'''
def getPrimes(level):
    """primepcn(level): return a list of all the prime numbers using Prefect Composite Number Theory founded by Honggang(Grant), Zhao from 2017-11-11"""
    # parameter is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30} 5|{210}, 6|2310,7:30030,8:510510,9:9699690,10:223092870
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>, from 2017-11-16
    # Code: https://github.com/zhaohonggang/BuildingPrimeNumbers/primepcn.py
    # Info: https://github.com/zhaohonggang/BuildingPrimeNumbers/blob/master/building%20prime%20number.txt

    c = 1 # current Prefect Composite Number
    c_n = 2 # next Prefect Composite Number
    p_l = [] # Prime number list
    p_n = 2 # Next Prime number

    for n in xrange(1,level):
        print "min {} max {}".format(c, c_n)

        # Prepare Composite Number Except List will appear in the generation
        e = [] if n-1 >= len(p_l) else list(set(getArrayProducts(p_l[n-1:], c, c_n)))
        print "Except List Prepared:{}"#.format(v)

        # Sort the Except list to be used in the next stage
        e.sort()
        print "Except List Sorted"

        s = [1] + ([] if n-1 >= len(p_l) else list(set(getArrayProducts(p_l[n-1:], 2, c))))
        print "Seed List Prepared:{}"#.format(s)

        # Sort the validation list to be used in the next stage
        s.sort()
        print "Seed List Sorted"

        e_i = 0 # Except List Index
        # generate prime numbers based on the Prefect Composite Number and the Co-prime list
        for k in xrange(1, p_n):
            for i in s:
                t = c * k + i

                # prime validation. If t is not in Except List, t is prime.
                if len(e) > e_i:
                    if t != e[e_i]:
                        p_l.append(t)
                    else:
                        e_i = e_i + 1
                else:
                    p_l.append(t)

        print "primes generated"

        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]

        c = c_n
        c_n = c * p_n

        del e[:]
        del s[:]
        print "ready for the next round"

    print 'count:{}\n'.format(len(p_l))
    return p_l

'''
a = [7, 11, 13]
b = list(getArrayProducts([5], 6, 30))
from primepcn import getArrayProducts
'''
def getArrayProducts(array, min, max):
    # multiply sorted array numbers any times including 0,1,2..n.., yield results between min and max
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>
    if len(array) == 0:
        return

    array_0 = array[0]

    max_1 = float(max) / array_0

    max_level = int(log(max, array[0])) - 1

    min_level = int(log(min, array[-1]))

    t_max = max / (array_0 ** min_level)

    max_i = 0

    while max_i < len(array) and array[max_i] < t_max:
        max_i = max_i + 1

    #print "array:{},min:{},max:{},max_i:{},t_max:{}".format(array, min, max, max_i, t_max)

    l = 0
    stack = []

    if min_level == 0:# push or yield
        for i in xrange(0, max_i):
            t = array[i]
            if t > max:
                break
            stack.append(t)
            if min < t:
                yield t
    else:# push only
        for i in xrange(0, max_i):
            t = array[i]
            if t > max:
                break
            stack.append(t)

    while l <= max_level - 1:
        #print stack
        stack_n = []
        if l < min_level - 1: # push only
            for j in xrange(0, len(stack)):
                a = stack[j]
                for i in xrange(0, max_i):
                    if array[i] > a:
                        break
                    t = a * array[i]
                    if t > max:
                        break
                    if t < max_1:
                        stack_n.append(t)
        elif l == max_level - 1: # yield only
            for j in xrange(0, len(stack)):
                a = stack[j]
                for i in xrange(0, max_i):
                    if array[i] > a:
                        break
                    t = a * array[i]
                    if t > max:
                        break
                    if min < t:
                        yield t
        else: # push or yield
            for j in xrange(0, len(stack)):
                a = stack[j]
                for i in xrange(0, max_i):
                    if array[i] > a:
                        break
                    t = a * array[i]
                    if t > max:
                        break
                    if min < t:
                        yield t
                    if t < max_1:
                        stack_n.append(t)

        l = l + 1
        del stack[:]
        stack = stack_n

'''
from primepcn import primepcn_v1
primepcn_v1(5)
a = primepcn_v1(10)
'''
def primepcn_v1(level):
    # it is the first version of primepcn. n_c_l is accumulating in this version.
    """primepcn_v1(level): return a list of all the prime numbers using Prefect Composite Number Theory founded by Honggang(Grant), Zhao from 2017-11-11"""
    # parameter is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30} 5|{210}, 6|2310,7:30030,8:510510,9:9699690,10:223092870
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>, from 2017-11-16
    # Code: https://github.com/zhaohonggang/BuildingPrimeNumbers/primepcn.py
    # Info: https://github.com/zhaohonggang/BuildingPrimeNumbers/blob/master/building%20prime%20number.txt

    c = 1 # current Prefect Composite Number
    c_n = 2 # next Prefect Composite Number
    n_c_l = [1] #lisk of numbers Co-prime with c
    p_l = [] #Prime number lisk
    p_n = 2 #Next Prime number

    for n in xrange(1,level):
        print "min {} max {}".format(c, c_n)

        # Prepare Composite Number validation List will apprear in the generation
        v = list(set(get_array_product_once(n_c_l[2:], c, c_n)))
        print "validation List Prepared"

        # Sort the validation list to be used in the next stage
        v.sort()

        print "validation List Sorted"

        v_i = 0
        len_n_c_l = len(n_c_l)

        # generate prime numbers based on the Prefect Composite Number and the Co-prime list
        for k in xrange(1, p_n):
            for i in n_c_l[:len_n_c_l]:

                t = c * k + i
                if t % p_n != 0:
                    n_c_l.append(t)
                elif t != p_n:
                    continue

                # prime validation
                while v_i < len(v) and v[v_i] < t:
                    v_i = v_i + 1
                if len(v) > v_i:
                    if t != v[v_i]:
                        p_l.append(t)
                    else:
                        v_i = v_i + 1
                else:
                    p_l.append(t)

        print "primes generated"

        # exit door
        if n == level - 1:
            break;

        # remove Next Prime number from Co-prime list since in the next round it will the factor of Prefect Composite Number
        for i in n_c_l[:len_n_c_l]:
            if i * p_n < c and i * p_n in n_c_l:
                n_c_l.remove(i * p_n)

        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]

        c = c_n
        c_n = c * p_n

        print "ready for the next round"

    print 'count:{}\n'.format(len(p_l))
    return p_l

def get_array_product_once(array, min, max):
    # multiply sorted array numbers one time, yield results between min and max
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>
    l = len(array)
    if l == 0:
        return

    start = 0

    while [start] * array[start] < min:
        start = start + 1

    end = start

    while end < l - 1 and array[end] * array[end] < max :
        end = end + 1

    if array[end] * array[end] > max:
        end = end - 1

    for i in xrange(start,l):
       for j in xrange(0, (i if i < end else end) + 1):
            t = array[i] * array[j]
            if min < t < max:
                yield t
            if t > max:
                break



        
if __name__ == '__main__':
    b = list(getArrayProducts([5], 6, 30))
