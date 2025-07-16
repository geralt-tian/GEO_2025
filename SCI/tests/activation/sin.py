
import numpy as np
import math
import random
import matplotlib.pyplot as plt


def unsigned_to_signed(x, N):
    if x >= N/2:
        return x - N
    return x



def compute_MW_plain(x0, x1, N):
    if 0 <= x0 + x1 < N/2:
        MW = 0
    if N/2 <= x0 + x1 < 3*N/2:
        MW = 1
    if 3*N/2 <= x0 + x1 < 2*N:
        MW = 2
    return MW



N = 2**21
f = 11



def sin_test():
    print("*********** sin test ***********")

    x_real = 1.5
    x_ring = x_real * 2**f

    x0 = 1048676
    x1 = 1051546

    # x0 = 3000
    # x1 = 70

    MW = compute_MW_plain(x0, x1, N)
    print("MW: ", MW)


    print("math.sin(x0/2**f): ", math.sin(x0/2**f))
    print("math.cos(x0/2**f): ", math.cos(x0/2**f))
    print("math.sin(x1/2**f): ", math.sin(x1/2**f))
    print("math.cos(x1/2**f): ", math.cos(x1/2**f))

    print("math.sin(-MW * N/2**f): ", math.sin(-MW * N/2**f))
    print("math.cos(-MW * N/2**f): ", math.cos(-MW * N/2**f))


    sin_SS = math.sin(x0/2**f) * math.cos(x1/2**f) * math.cos(-MW * N/2**f) + \
            math.cos(x0/2**f) * math.sin(x1/2**f) * math.cos(-MW * N/2**f) + \
            math.cos(x0/2**f) * math.cos(x1/2**f) * math.sin(-MW * N/2**f) - \
            math.sin(x0/2**f) * math.sin(x1/2**f) * math.sin(-MW * N/2**f) 

    print("sin(x_real): ", math.sin(x_real))
    print("sin_SS: ", sin_SS)
    print()


    k = 4
    print("*********** round ***********")

    print("math.sin(x0/2**f): ", round(math.sin(x0/2**f),k))
    print("math.cos(x0/2**f): ", round(math.cos(x0/2**f),k))
    print("math.sin(x1/2**f): ", round(math.sin(x1/2**f),k))
    print("math.cos(x1/2**f): ", round(math.cos(x1/2**f),k))

    print("math.sin(-MW * N/2**f): ", round(math.sin(-MW * N/2**f),k))
    print("math.cos(-MW * N/2**f): ", round(math.cos(-MW * N/2**f),k))

    sin_SS = round(math.sin(x0/2**f),k) * round(math.cos(x1/2**f),k) * round(math.cos(-MW * N/2**f),k) + \
            round(math.cos(x0/2**f),k) * round(math.sin(x1/2**f),k) * round(math.cos(-MW * N/2**f),k) + \
            round(math.cos(x0/2**f),k) * round(math.cos(x1/2**f),k) * round(math.sin(-MW * N/2**f),k) - \
            round(math.sin(x0/2**f),k) * round(math.sin(x1/2**f),k) * round(math.sin(-MW * N/2**f),k) 

    print("sin(x_real): ", math.sin(x_real))
    print("sin_SS: ", sin_SS)
    print()


# sin_test()





def exp_test():
    print("*********** exp test ***********")

    N = 2**21
    f = 11

    x_real = 1.5
    x_ring = x_real * 2**f

    # x0 = 1048676
    # x1 = 1051546

    # x0 = 3000
    # x1 = 70

    # MW = compute_MW_plain(x0, x1, N)
    # print("MW: ", MW)


    N = 2**21
    f = 19

    x_real = -1.5
    x_ring = x_real * 2**f

    r = random.randint(0,N)
    r = 1
    x0 = r % N
    x1 = (x_ring - r) % N
    print("x0, x1: ", x0, x1)

    MW = compute_MW_plain(x0, x1, N)
    print("MW: ", MW)


    print("2**(x0/2**f): ", 2**(x0/2**f))
    print("2**(x1/2**f): ", 2**(x1/2**f))
    print("2**(-MW*N/2**f): ", 2**(-MW*N/2**f))


    # MW 
    # 一次乘法，一次trun-reduce，
    # 又一次乘法，一次trun-reduce，

    # 优势：精度高点？开销应该会高
    exp_SS = 2**(x0/2**f) * 2**(x1/2**f) * 2**(-MW*N/2**f)

    print("exp(x_real): ", 2**x_real)
    print("exp_SS: ", exp_SS)
    print()


    k = 4
    print("*********** round ***********")

    print("2**(x0/2**f): ", round(2**(x0/2**f),k))
    print("2**(x1/2**f): ", round(2**(x1/2**f),k))
    print("2**(-MW*N/2**f): ", round(2**(-MW*N/2**f),k))

    exp_SS = round(2**(x0/2**f),k) * round(2**(x1/2**f),k) * round(2**(-MW*N/2**f),k)

    print("exp(x_real): ", 2**x_real)
    print("exp_SS: ", exp_SS)



exp_test()





