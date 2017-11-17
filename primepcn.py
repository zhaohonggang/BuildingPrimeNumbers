def multiply_array(array, min, max):
    # multiply array numbers, yield results between min and max
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



'''
from primes import primepcn
primepcn(1)
'''
def primepcn(level):
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
        v = list(multiply_array(n_c_l[2:], c, c_n))

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