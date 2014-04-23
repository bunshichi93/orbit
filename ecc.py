#draft 
from math import sqrt, atan2 
from cmath import polar
G = 6.67E-11
ms = 1.9891E+30 # massa Sole
tau = 3600
r_o = [46001272000, 107476002000, 147098074000, 206644545000]
v_o = [58980, 35259, 30287, 26499]

def calc_ecc(rin, vin):

    r_x = [0]
    r_y = [0]
    phi = [0]
    r = [0]
    # velocit #
    # in m/s   #
    v_x = [0]
    v_y = [0]
    # accelerazioni #
    # in N          #
    a_x = [0]
    a_y = [0]
    
    # condizioni iniziali istante 0
    r_x[0] = rin
    r_y[0] = 0

    r[0] = sqrt(r_x[0]**2 + r_y[0]**2)

    v_x[0] = 0
    v_y[0] = vin 

    a_x[0] = -G*ms*r_x[0]/r[0]**3
    a_y[0] = -G*ms*r_y[0]/r[0]**3
    time = 0
    i = 0
    coo_p = [0, 1]
    while coo_p[1] >= 0:
        r_x.append(r_x[i] + v_x[i]*tau + 0.5*a_x[i]*tau**2)
        r_y.append(r_y[i] + v_y[i]*tau + 0.5*a_y[i]*tau**2)
        coo_p = polar(r_x[i+1] + r_y[i+1]*1j)
        r.append(coo_p[0])
        a_x.append(-G*ms*r_x[i+1]/r[i+1]**3)
        a_y.append(-G*ms*r_y[i+1]/r[i+1]**3)
        v_x.append(v_x[i] + 0.5*tau*(a_x[i] + a_x[i+1]))
        v_y.append(v_y[i] + 0.5*tau*(a_y[i] + a_y[i+1]))
        time += tau
        i += 1


    a = min(r)
    b = max(r)
    return a, b, time


per = [0]*len(r_o)
afe = [0]*len(r_o)
T = [0]*len(r_o)

for i in range(len(r_o)):
     c = calc_ecc(r_o[i], v_o[i])
     per[i] = c[0]
     afe[i] = c[1]
     T[i] = c[2]*2/3600./24.
print per
print afe
print T
print "\n"
