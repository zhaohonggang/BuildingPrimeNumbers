from math import log
'''
from primepcn import primepcn
primepcn(5)
a = primepcn(10)
'''
def primepcn(level):
    # compare to v1, it is simpler and clearer describe the theory
    """primepcn(level): return a list of all the prime numbers using Prefect Composite Number Theory founded by Honggang(Grant), Zhao from 2017-11-11"""
    # parameter is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30} 5|{210}, 6|2310,7:30030,8:510510,9:9699690,10:223092870
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>, from 2017-11-16
    # Code: https://github.com/zhaohonggang/BuildingPrimeNumbers/
    # Info: https://github.com/zhaohonggang/BuildingPrimeNumbers/blob/master/building%20prime%20number.txt
    # Math: https://github.com/zhaohonggang/BuildingPrimeNumbers/blob/master/prime%20algorithm.txt

    c = 1 # current Prefect Composite Number
    c_n = 2 # next Prefect Composite Number
    p_l = [] #Prime number lisk
    p_n = 2 #Next Prime number

    for n in xrange(1,level):
        print "min {} max {}".format(c, c_n)

        # Prepare Composite Number validation List will apprear in the generation
        v = [] if n-1 >= len(p_l) else list(set(array_product_all_nr(p_l[n-1:], c, c_n)))

        print "validation List Prepared:{}"#.format(v)

        # Sort the validation list to be used in the next stage
        v.sort()
        print "validation List Sorted"

        s = [1] if n-1 >= len(p_l) else list(set(array_product_all_nr_and_one(p_l[n-1:], 2, c)))

        print "seed List Prepared:{}"#.format(s)

        # Sort the validation list to be used in the next stage
        s.sort()
        print "seed List Sorted"

        v_i = 0
        # generate prime numbers based on the Prefect Composite Number and the Co-prime list
        for k in xrange(1, p_n):
            for i in s:
                t = c * k + i

                # prime validation
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

        if n == 1:
            p_n = 3
        else:
            p_n = p_l[n]

        c = c_n
        c_n = c * p_n

        print "ready for the next round"

    print 'count:{}\n'.format(len(p_l))
    return p_l

'''
from primepcn import primepcn_v1
primepcn_v1(5)
a = primepcn_v1(10)
'''
def primepcn_v1(level):
    """primepcn(level): return a list of all the prime numbers using Prefect Composite Number Theory founded by Honggang(Grant), Zhao from 2017-11-11"""
    # parameter is the level of profect composite number 1|{1}, 2|{2}, 3|{6}, 4|{30} 5|{210}, 6|2310,7:30030,8:510510,9:9699690,10:223092870
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>, from 2017-11-16
    # Code: https://github.com/zhaohonggang/BuildingPrimeNumbers/
    # Info: https://github.com/zhaohonggang/BuildingPrimeNumbers/blob/master/building%20prime%20number.txt
    # Math: https://github.com/zhaohonggang/BuildingPrimeNumbers/blob/master/prime%20algorithm.txt

    c = 1 # current Prefect Composite Number
    c_n = 2 # next Prefect Composite Number
    n_c_l = [1] #lisk of numbers Co-prime with c
    p_l = [] #Prime number lisk
    p_n = 2 #Next Prime number

    for n in xrange(1,level):
        print "min {} max {}".format(c, c_n)

        # Prepare Composite Number validation List will apprear in the generation
        v = list(set(array_product_2(n_c_l[2:], c, c_n)))
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

def array_product_2(array, min, max):
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

def array_product_all_n(array, min, max, n, mn):
    # multiply sorted array numbers for any possible ways, yield results between min and max
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>
    if n == 0:
        for a in array:
            if a < min:
                yield 1
                return
            if a > max:
                yield 1
                return
            yield a
    else:
        max_t = max **(1.0/(mn+1))
        for a in array:
            if a > max_t:
                yield 1
                return
            for b in array_product_all_n(array, min/a , max/a + 1, n - 1, mn):
                yield a * b

'''
a = [7, 11, 13]
b = list(array_product_all([7, 11, 13], 0, 200))
from primepcn import array_product_all
'''
def array_product_all(array, min, max):
    # multiply sorted array numbers for any possible ways, yield results between min and max
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>
    if len(array) == 0:
        return

    n = int(log(max, array[0])) - 1

    mn = int(log(min, array[-1])) - 1

    if mn < 0:
        mn = 0

    #yield (item for item in array_product_all_n(array, min, max, n) if min <item <max)
    for a in array_product_all_n(array, min, max, n, mn):
        if min <a <max:
            yield a

def array_product_all_and_one(array, min, max):
    # multiply sorted array numbers for any possible ways, yield results between min and max
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>
    yield 1
    for a in array_product_all(array, min, max):
        yield a

def array_product_all_nr_and_one(array, min, max):
    # multiply sorted array numbers for any possible ways, yield results between min and max
    # Code by Honggang(Grant), Zhao, <honggangz@gmail.com>
    yield 1
    for a in array_product_all_nr(array, min, max):
        yield a

'''
a = [7, 11, 13]
b = list(array_product_all_nr([5], 6, 30))
from primepcn import array_product_all_nr
'''
def array_product_all_nr(array, min, max):
    if len(array) == 0:
        return

    array_0 = array[0]

    max_level = int(log(max, array[0])) - 1

    min_level = int(log(min, array[-1])) - 1

    if min_level < 0:
        min_level = 0

    t_max = max ** (1.0 / (min_level + 1))

    max_i = 0

    while max_i < len(array) and array[max_i] < t_max:
        max_i = max_i + 1

    #print "array:{},min:{},max:{},max_i:{}".format(array, min, max, max_i)

    l = 0
    stack = []
    for i in xrange(0, max_i):
        t = array[i]
        if t > max:
            break
        if min < t < max:
            yield t
        if t * array_0 < max:
            stack.append(t)
    len_p = 0

    while l <= max_level - 1:
        #print stack
        len_n = len(stack)
        for j in xrange(len_p, len_n):
            a = stack[j]
            for i in xrange(0, max_i):
                if array[i] > a:
                    break
                t = a * array[i]
                if t > max:
                    break
                if min < t < max:
                    yield t
                if l < max_level - 1 and t * array_0 < max:
                    stack.append(t)

        l = l + 1
        len_p = len_n



'''
a = [7, 11, 13]
b = list(array_product_all_nr_no_index([5], 6, 30))
from primepcn import array_product_all_nr_no_index
'''
def array_product_all_nr_no_index(array, min, max):
    if len(array) == 0:
        return

    array_0 = array[0]

    max_level = int(log(max, array[0])) - 1

    min_level = int(log(min, array[-1])) - 1

    if min_level < 0:
        min_level = 0

    t_max = max ** (1.0 / (min_level + 1))

    #print "array:{},min:{},max:{},max_i:{}".format(array, min, max, max_i)
    max_i = 0

    while max_i < len(array) and array[max_i] < t_max:
        max_i = max_i + 1

    l = 0
    stack = []
    for i in array[0:max_i]:
        if min < i:
            yield i
        if i * array_0 < max:
            stack.append(i)
    len_p = 0

    while l <= max_level - 1:
        #print stack
        len_n = len(stack)
        for a in stack[len_p:len_n]:
            for i in array[0:max_i]:
                if i > a:
                    break
                t = a * i
                if t > max:
                    break
                if min < t:
                    yield t
                if l < max_level - 1 and t * array_0 < max:
                    stack.append(t)
        l = l + 1
        len_p = len_n


        
if __name__ == '__main__':
    b = list(array_product_all([7, 11, 13], 0, 200))

    primepcn(5)
