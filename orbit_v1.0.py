# -*- coding: utf-8 -*-
from cmath import polar
from math import sqrt, pi
import sys

####################
#####   DATI   #####
####################
G = 6.67E-11
# masse #
# in Kg #
ms = 1.9891E+30 # massa Sole
mt = 5.9743E+24 # massa Terra
mp = 1.8986E+27 # massa Giove
# !!!perielio!!!
# velocità è forza gravitazionale sono perpendicolari
# raggi iniziali #
# in m           #
r0_xt = 1.47098074E+11
r0_yt = 0
r0_xp = 7.40742598E+11
r0_yp = 0
# moduli #
rt = sqrt(r0_xt**2 + r0_yt**2)
rp = sqrt(r0_xp**2 + r0_yp**2)
# velocità iniziali #
# in m/s            #
v0_xt = 0
v0_yt = 3.0287E+4
v0_xp = 0
v0_yp = 1.3712E+4


def ask_params():
    global Nstep, tau, mt, mp, r0_xt, r0_yt, r0_xp, r0_yp, v0_xt, v0_yt, v0_xp, v0_yp

    print "Inserire passo temporale tau (in ore):"
    tau = float(raw_input(">>>"))*3600
    print 

    if operation == "3":
        return None

    best_Nstep = 365.*24.*3600/tau
    print "Una singola orbita Terrestre verrà eseguita in ~%.3f passi con tau di %s ore):" % (best_Nstep, tau/3600)
    print "Inserire numero di passi:"
    Nstep = int(raw_input(">>>"))
    print 

    if Nstep < best_Nstep:

        print "Il programma non sarà in grado di completare una singola orbita terrestre."
        print 

    print "Usare parametri predefiniti? [Y/n/q]"
    while True:
        answer = raw_input(">>>")
        if answer == "Y" or answer == "y":
            if operation == "2":
                mp = 0
                r0_xp = 1
                r0_yp = 1
                v0_xp = 0
                v0_yp = 0
            break
        elif answer == "N" or answer == "n":
            print
            print "Inserire i parametri della simulazione."
            print
            print "Massa pianeta #1 (kg)"
            mp = float(raw_input(">>>"))
            print "Raggio iniziale orbita #1 (km)"
            print "[componente x]"
            r0_xp = float(raw_input(">>>"))*1000
            print "[componente y]"
            r0_yp = float(raw_input(">>>"))*1000
            print "Velocità iniziale pianeta #1 (km/s)"
            print "[componente x]"
            v0_xp = float(raw_input(">>>"))*1000
            print "[componente y]"
            v0_yp = float(raw_input(">>>"))*1000
            if operation == "2":
                mt = mp
                r0_xt = r0_xp
                r0_yt = r0_yp
                v0_xt = v0_xp
                v0_yt = v0_yp
                mp = 0
                r0_xp = 1
                r0_yp = 1
                v0_xp = 0
                v0_yp = 0
            break
        elif answer == "q":
            print
            print "Bye!"
            break
        else:
            print
            print "Opzione non valida."


def drawProgressBar(percent, barLen):
    """Funzione che disegnerà una prograss bar sul terminale"""

    sys.stdout.write("\r")
    progress = "="*int(barLen * percent) + ">" + " "*(barLen - int(barLen * percent))
    sys.stdout.write("[%s] %.2f%%" % (progress, percent * 100))
    sys.stdout.flush()


def orbit(cond):
    """Funzione che simula l'orbita di Giove e della Terra, tenendo conto 
    dell'interazione Giove-Terra"""
    global r, v_x, v_y
    ##############
    # Crea matrici 2*Nstep (righe*colonne) con zeri. 
    # Prima riga dati Terra, seconda riga dati Giove.
    ###############
    # coordinate #
    # in m e rad #
    # Coordinate cartesiane (r_x, r_y) e polari (phi, r) dei pianeti.
    r_x = [[0]*2, [0]*2]
    r_y = [[0]*2, [0]*2]
    phi = [[0]*Nstep, [0]*Nstep] 
    r = [[0]*Nstep, [0]*Nstep]   
    # velocità #
    # in m/s   #
    v_x = [[0]*Nstep, [0]*Nstep]
    v_y = [[0]*Nstep, [0]*Nstep]
    # accelerazioni #
    # in N          #
    a_x = [[0]*2, [0]*2]
    a_y = [[0]*2, [0]*2]
    
    # condizioni iniziali istante 0
    r_x[0][0] = r0_xt
    r_y[0][0] = r0_yt
    r_x[1][0] = r0_xp
    r_y[1][0] = r0_yp 

    for i in range(2):
        coo_p = polar(r_x[i][0] + r_y[i][0]*1j)
        r[i][0] = coo_p[0]
        phi[i][0] = coo_p[1]

    r_rel = sqrt((r_x[1][0] - r_x[0][0])**2 + (r_y[1][0] - r_y[0][0])**2)    

    v_x[0][0] = v0_xt 
    v_y[0][0] = v0_yt 
    v_x[1][0] = v0_xp 
    v_y[1][0] = v0_yp 

    a_x[0][0] = -G*mp*(r_x[1][0] - r_x[0][0])/r_rel**3 -G*ms*r_x[0][0]/r[0][0]**3
    a_y[0][0] = -G*mp*(r_y[1][0] - r_y[0][0])/r_rel**3 -G*ms*r_y[0][0]/r[0][0]**3
    a_x[1][0] = cond*(-G*mt*(r_x[0][0] - r_x[1][0])/r_rel**3 -G*ms*r_x[1][0]/r[1][0]**3)
    a_y[1][0] = cond*(-G*mt*(r_y[0][0] - r_y[1][0])/r_rel**3 -G*ms*r_y[1][0]/r[1][0]**3)

    print
    print "Calcolo dati:"
    for i in range(Nstep-1):
        ##############################
        ##############################
        ### Terra
        # calcolo dati sfruttando Velocity-Verlet
        # parte 1
        r_x[0][1] = r_x[0][0] + v_x[0][i]*tau + 0.5*a_x[0][0]*tau**2
        r_y[0][1] = r_y[0][0] + v_y[0][i]*tau + 0.5*a_y[0][0]*tau**2
        # conversione coordinate cartesiane in polari
        coo_p = polar(r_x[0][1] + r_y[0][1]*1j)
        r[0][i+1] = coo_p[0]
        phi[0][i+1] = coo_p[1]
        ### Giove
        # calcolo dati sfruttando Velocity-Verlet
        # parte 1
        r_x[1][1] = r_x[1][0] + v_x[1][i]*tau + 0.5*a_x[1][0]*tau**2
        r_y[1][1] = r_y[1][0] + v_y[1][i]*tau + 0.5*a_y[1][0]*tau**2
        # conversione coordinate cartesiane in polari
        coo_p = polar(r_x[1][1] + r_y[1][1]*1j)
        r[1][i+1] = coo_p[0]
        phi[1][i+1] = coo_p[1]

        r_rel = sqrt((r_x[0][1] - r_x[1][1])**2 + (r_y[0][1] - r_y[1][1])**2)

        # calcolo dati sfruttando Velocity-Verlet
        # parte 2
        ### Terra
        a_x[0][1] = -G*ms*r_x[0][1]/r[0][i+1]**3 -G*mp*(r_x[1][1] - r_x[0][1])/r_rel**3
        a_y[0][1] = -G*ms*r_y[0][1]/r[0][i+1]**3 -G*mp*(r_y[1][1] - r_y[0][1])/r_rel**3
        v_x[0][i+1] = v_x[0][i] + 0.5*tau*(a_x[0][0] + a_x[0][1])
        v_y[0][i+1] = v_y[0][i] + 0.5*tau*(a_y[0][0] + a_y[0][1])
        ### Giove
        a_x[1][1] = cond*(-G*ms*r_x[1][1]/r[1][i+1]**3 - G*mt*(r_x[0][1] - r_x[1][1])/r_rel**3)
        a_y[1][1] = cond*(-G*ms*r_y[1][1]/r[1][i+1]**3 - G*mt*(r_y[0][1] - r_y[1][1])/r_rel**3)
        v_x[1][i+1] = v_x[1][i] + 0.5*tau*(a_x[1][0] + a_x[1][1])
        v_y[1][i+1] = v_y[1][i] + 0.5*tau*(a_y[1][0] + a_y[1][1])
        ##############################
        ##############################

        r_x[0][0] = r_x[0][1]
        r_y[0][0] = r_y[0][1]
        r_x[1][0] = r_x[1][1]
        r_y[1][0] = r_y[1][1]
        a_x[0][0] = a_x[0][1]
        a_y[0][0] = a_y[0][1]
        a_x[1][0] = a_x[1][1]
        a_y[1][0] = a_y[1][1]

        drawProgressBar(float(i)/Nstep, 60)

    print
    print "Fine calcolo dati."
    print
    

    # apre i file su cui scrivere i dati con la proprietà 'w'
    orbit_g = open('orbit_p.dat', 'w') # Pianeta
    orbit_t = open('orbit_t.dat', 'w') # Terra
    
    for j in range(1+1):#cond +1):
        if j == 0:
            planet = "Terra"
            filename = orbit_t                
        elif j == 1:
            planet = ""
            filename = orbit_g
        print "Scrittura dati %s su file:" % planet

        for i in range(Nstep):
            # scrive i dati su file
            # j = 0 > Terra
            # j = 1 > Giove
            filename.write(str(phi[j][i]))
            filename.write("\t")
            filename.write(str(r[j][i]))
            filename.write("\n")
            
            drawProgressBar(float(i)/Nstep, 60)

        print
        print "Fine scrittura dati %s su file." % planet
        print
            
    # chiude i file su cui sono stati scritti i dati
    orbit_g.close()
    orbit_t.close()
    
    ape_t = max(r[0])
    per_t = min(r[0])
    ecc_t = (ape_t - per_t)/(ape_t + per_t)
    ape_g = max(r[1])
    per_g = min(r[1])
    ecc_g = (ape_g - per_g)/(ape_g + per_g)        
    print ape_t, per_t, ecc_t
    print ape_g, per_g, ecc_g


def calc_data():
    vel_ang = open('omega.dat', 'w')
    kin = open('e_cinetica.dat', 'w')
    pot = open('e_potenziale', 'w')
    tot = open('e_totale', 'w')
    time = [0]*Nstep
    v = [0]*Nstep
    omega = [0]*Nstep
    K = [0]*Nstep
    U = [0]*Nstep
    Tot = [0]*Nstep
    for k in range(Nstep):
        time[k] = k*tau
        v = sqrt(v_x[0][k]**2 + v_y[0][k]**2)
        omega[k] = v/r[0][k]
        K[k] = 0.5*omega[k]**2*mt*r[0][k]**2
        U[k] = -G*ms*mt/r[0][k]
        Tot[k] = K[k] + U[k]
    
        vel_ang.write(str(r[0][k]))
        vel_ang.write("\t")
        vel_ang.write(str(omega[k]))
        vel_ang.write("\n")
        
        kin.write(str(time[k]))
        kin.write("\t")
        kin.write(str(K[k]))
        kin.write("\n")
        
        pot.write(str(time[k]))
        pot.write("\t")
        pot.write(str(U[k]))
        pot.write("\n")
        
        tot.write(str(time[k]))
        tot.write("\t")
        tot.write(str(Tot[k]))
        tot.write("\n")


def kepler(rin, vin, plan):

    tau = 600
    r_x = rin
    r_y = 0
    v_x = 0
    v_y = vin
    coo_p = polar(r_x + r_y*1j)
    r = coo_p[0]
    a_x = [-G*ms*r_x/r**3,[] ]
    a_y = [-G*ms*r_y/r**3,[] ]
    time = 0
    while coo_p[1] >= 0:
        r_x += v_x*tau + 0.5*a_x[0]*tau**2
        r_y += v_y*tau + 0.5*a_y[0]*tau**2

        coo_p = polar(r_x + r_y*1j)
        r = coo_p[0]

        a_x[1] = -G*ms*r_x/r**3
        a_y[1] = -G*ms*r_y/r**3
        v_x += 0.5*tau*(a_x[0] + a_x[1])
        v_y += 0.5*tau*(a_y[0] + a_y[1])

        a_x[0] = a_x[1]
        a_y[0] = a_y[1]
        
        time += tau

    print "Dati relativi a ", plan
    T = (time)*2
    print "Periodo orbitale:\t", T/(24.*60*60), " giorni"
    semiax = (abs(rin) + abs(r_x))/2
    print "Semiasse maggiore:\t", semiax/1000, " km"

    cost_K = T**2/semiax**3
    print "Costante Keplero:\t", cost_K

print """Programme bla bla bla
1) Simulare l'orbita di un sistema Sole - Terra - pianeta
2) Simulare l'orbita e calcolare i parametri di un sistema Sole - pianeta
3) Verificare la legge di Keplero per i pianeti del Sistema Solare
q) Per uscire"""

while True:
    operation = raw_input("Scegliere l'operazione da eseguire: ")
    print
    if operation == "1":
        print "Simulazione orbita sistema Sole - Terra - pianeta"
        ask_params()
        orbit(1)
        break
    elif operation == "2":
        print "Simulazione orbita e calcolo parametri di un sistema Sole - pianeta"
        ask_params()
        orbit(0)
        calc_data()
        break
    elif operation == "3":
        ask_params()
        plan = ["Mercurio", "Venere", "Terra", "Marte", "Giove", "Saturno", "Urano", "Nettuno"]
        rad = [0.46001272E+11, 1.07476002E+11, 1.47098074E+11, 2.06664545e+11, 7.40742598E+11, 13.49467375E+11, 27.35555035E+11, 44.59631496E+11]
        vel = [5.898E+4, 3.5020E+4, 3.0287E+4, 2.499E+4, 1.3712E+4, 1.0183E+4, 0.711E+4, 0.5479E+4]
        for i in range(len(rad)):
            kepler(rad[i], vel[i], plan[i])
        print "Teorica =", 4*pi**2/(G*ms)
        break
    elif operation =="q":
        break
    else:
        print "Opzione non valida.\n"
