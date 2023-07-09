import math
from time import time

from numpy import mat
# k=[]
# b=[]
# strip=1
# def approximation_log(x):
#     z=0
#     for i in range(len(b)):
#         if x>=b[i] and x<b[i+1]:
#             z=i
#             break
#     return math.log(b[i])+k[i]*(z-b[i])+1
# def approximation2_log(x):
#     z=0
#     for i in range(1,30):
#         if x>=math.pow(2,i)+math.pow(2,i-1) and x<math.pow(2,i+1)+math.pow(2,i):
#             z=i+1
#             break
#     res=x/math.pow(2,z)-1
#     r=res-res*res/2+res**3/3

#     return z*math.log(2)+r
# strip=1
# for i in range(1,10,strip):
#       k.append((math.log(i+strip)-math.log(i))/strip)
#       b.append(i)

# strip=1
# for i in range(10,100,strip):
#       k.append((math.log(i+strip)-math.log(i))/strip)
#       b.append(i)
# strip=10
# for i in range(100,1000,strip):
#       k.append((math.log(i+strip)-math.log(i))/strip)
#       b.append(i)
# strip=100
# for i in range(1000,10000,strip):
#       k.append((math.log(i+strip)-math.log(i))/strip)
#       b.append(i)
# strip=1000
# for i in range(10000,100000,strip):
#       k.append((math.log(i+strip)-math.log(i))/strip)
#       b.append(i)

# strip=10000
# for i in range(100000,1000000,strip):
#       k.append((math.log(i+strip)-math.log(i))/strip)
#       b.append(i)

# strip=100000
# for i in range(1000000,10000000,strip):
#       k.append((math.log(i+strip)-math.log(i))/strip)
#       b.append(i)

# strip=1000000
# for i in range(10000000,100000000,strip):
#       k.append((math.log(i+strip)-math.log(i))/strip)
#       b.append(i)



# max=0
# for i in range(1,100000000):
#     t=time.clock()
#     re=approximation_log(i)
#     # if re>max:
#     #     max=re
#     #     print(max)
#     #     print(i)
#     if i%100000==0:

#         print('temp'+str(re))


# print(max)
# x=2**(-30)
# while (1):
#     x=x*(2-(2**25)*x)
#     print(x)
#     if(math.fabs(x-1/(2**25)))<0.0001:
#         print('successful')
#         break
# print(1055530156032/math.pow(2,22))
# print(math.log(math.pow(2,24)-1))
# print(math.log(0.5))

def apexp(x):
    k=(2**8)
    return (1+(x/k))**256
def aplogari(x):
    y=x/120+apexp(-2*x-1)+3
    for z in range(2):
        h=1-x*apexp(-y)
        y-=polia(h)
    return y
def polia(h):
    k=0
    for i in range(1,9):
        k+=math.pow(h,i)/i
    return k

print(aplogari(1000))
print(math.log(1000))