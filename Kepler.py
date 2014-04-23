from cmath import polar
Nstep = 4500000   
tau = 600


M =[0]*8
R =[0]*8
v =[0]*8
A =[0]*8
pi=3.14159265358979323846

R[0]=0.30749951
R[1]=0.71843270
R[2]=0.9832898912
R[3]=1.38133346
R[4]=4.95155843
R[5]=9.02063224
R[6]=18.28605596
R[7]=29.81079527

v[0]=58980
v[1]=35259
v[2]=30287
v[3]=26499
v[4]=13712
v[5]=10183
v[6]=7110
v[7]=5479



F_x = [0]*Nstep
F_y = [0]*Nstep

ms = 1.9891E+30
G = 6.67384*1.E-11

for k in range(8):     
    
    R[k] =R[k]*149597870700
    r_x = R[k]
    r_y = 0

    r = r_x
    phi = 0

    v_x = 0
    v_y = v[k]

    
    F_x[0] = -G*ms*r_x/r**3
    F_y[0] = -G*ms*r_y/r**3
    
    for i in range(Nstep-1):
        
        if phi >= 0:
            
            time = tau*i
        
            r_x = r_x + v_x*tau + 0.5*F_x[i]*tau**2
            r_y = r_y + v_y*tau + 0.5*F_y[i]*tau**2
        
            coo_p = polar(r_x + r_y*1j)
            r = coo_p[0]
            phi = coo_p[1]
        
        
            F_x[i+1] = -G*ms*r_x/r**3
            F_y[i+1] = -G*ms*r_y/r**3
            v_x += 0.5*tau*(F_x[i] + F_x[i+1])
            v_y += 0.5*tau*(F_y[i] + F_y[i+1])
        
     
            semiassecubo = (r+R[k])**3/8
            periodoquadro = 4*time**2
            A[k] = periodoquadro/semiassecubo
                           
print 'Mercurio', A[0]
print 'Venere', A[1]
print 'Terra', A[2]
print 'Marte' , A[3]
print 'Giove', A[4]
print 'Saturno' , A[5]
print 'Urano', A[6]
print 'Nettuno', A[7]

print 'In teoria' , 4*pi**2/(G*ms)
