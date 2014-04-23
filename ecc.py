#draft 
from math import sqrt, atan2 
G = 6.67E-11
ms = 1.9891E+30 # massa Sole
Nstep = 1000000
tau = 3600
r_o = [46001272000, 107476002000, 147098074000, 206644545000]
v_o = [58980, 35259, 30287, 26499]

def calc_ecc(rin, vin):

    r_x = [0]*Nstep
    r_y = [0]*Nstep
    phi = [0]*Nstep
    r = [0]*Nstep
    # velocit #
    # in m/s   #
    v_x = [0]*Nstep
    v_y = [0]*Nstep
    # accelerazioni #
    # in N          #
    a_x = [0]*Nstep
    a_y = [0]*Nstep
    
    # condizioni iniziali istante 0
    r_x[0] = rin
    r_y[0] = 0

    r[0] = sqrt(r_x[0]**2 + r_y[0]**2)

    v_x[0] = 0
    v_y[0] = vin 

    a_x[0] = -G*ms*r_x[0]/r[0]**3
    a_y[0] = -G*ms*r_y[0]/r[0]**3


    for i in range(Nstep-1):
     
        r_x[i+1] = r_x[i] + v_x[i]*tau + 0.5*a_x[i]*tau**2
        r_y[i+1] = r_y[i] + v_y[i]*tau + 0.5*a_y[i]*tau**2
        r[i+1] = sqrt(r_x[i+1]**2 + r_y[i+1]**2)

        a_x[i+1] = -G*ms*r_x[i+1]/r[i+1]**3
        a_y[i+1] = -G*ms*r_y[i+1]/r[i+1]**3
        v_x[i+1] = v_x[i] + 0.5*tau*(a_x[i] + a_x[i+1])
        v_y[i+1] = v_y[i] + 0.5*tau*(a_y[i] + a_y[i+1])
     
        
    a = min(r)
    b = max(r)
    return a, b

calc_ecc(1.47E11, 3.0287E4)
per = [0]*len(r_o)
afe = [0]*len(r_o)
for i in range(len(r_o)):
     c = calc_ecc(r_o[i], v_o[i])
     per[i] = c[0]
     afe[i] = c[1]
    
print per
print afe
