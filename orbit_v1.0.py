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
    print "Inserire numero di passi (singola orbita Terrestre - 9000 passi con tau == 1):"
    Nstep = int(raw_input(">>>"))
    
    print "Usare parametri predefiniti? [Y/n/q]"
    while True:
        answer = raw_input(">>>")
        if answer == "Y" or answer == "y":
            if operation == "2":
                mp = 0
                v0_xp = 0
                v0_yp = 0
            break
        elif answer == "N" or answer == "n":
            print "Inserire i parametri della simulazione."
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
            print "Bye!"
            break
        else:
            print "Opzione non valida."


def ask_data_on_file():
    global data_on_file
    print "Scrivere dati su file? [Y/n/q]"
    while True:
        answer = raw_input(">>>")
        if answer == "Y" or answer == "y":
            data_on_file = True
            break
        elif answer == "N" or answer == "n":
            data_on_file = False
            break
        elif answer == "q":
            print "Bye!"
            break
        else:
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

    ##############
    # Crea matrici 2*Nstep (righe*colonne) con zeri. 
    # Prima riga dati Terra, seconda riga dati Giove.
    ###############
    # coordinate #
    # in m e rad #
    # Coordinate cartesiane (r_x, r_y) e polari (phi, r) dei pianeti.
    r_x = [[0]*Nstep, [0]*Nstep]
    r_y = [[0]*Nstep, [0]*Nstep]
    phi = [[0]*Nstep, [0]*Nstep] 
    r = [[0]*Nstep, [0]*Nstep]   
    # velocità #
    # in m/s   #
    v_x = [[0]*Nstep, [0]*Nstep]
    v_y = [[0]*Nstep, [0]*Nstep]
    # accelerazioni #
    # in N          #
    a_x = [[0]*Nstep, [0]*Nstep]
    a_y = [[0]*Nstep, [0]*Nstep]
    
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
    a_x[1][0] = -G*ms*r_x[1][0]/r[1][0]**3*cond
    a_y[1][0] = -G*ms*r_y[1][0]/r[1][0]**3*cond

    print "Calcolo dati:"
    for i in range(Nstep-1):
        ##############################
        ##############################
        ### Giove
        # calcolo dati sfruttando Velocity-Verlet
        # parte 1
        r_x[1][i+1] = r_x[1][i] + v_x[1][i]*tau + 0.5*a_x[1][i]*tau**2
        r_y[1][i+1] = r_y[1][i] + v_y[1][i]*tau + 0.5*a_y[1][i]*tau**2
        # conversione coordinate cartesiane in polari
        coo_p = polar(r_x[1][i+1] + r_y[1][i+1]*1j)
        r[1][i+1] = coo_p[0]
        phi[1][i+1] = coo_p[1]
        # calcolo dati sfruttando Velocity-Verlet
        # parte 2
        a_x[1][i+1] = -G*ms*r_x[1][i+1]/r[1][i+1]**3*cond
        a_y[1][i+1] = -G*ms*r_y[1][i+1]/r[1][i+1]**3*cond
        v_x[1][i+1] = v_x[1][i] + 0.5*tau*(a_x[1][i] + a_x[1][i+1])
        v_y[1][i+1] = v_y[1][i] + 0.5*tau*(a_y[1][i] + a_y[1][i+1])
        ##############################
        ##############################
        
        ##############################
        ##############################
        ### Terra
        # calcolo dati sfruttando Velocity-Verlet
        # parte 1
        r_x[0][i+1] = r_x[0][i] + v_x[0][i]*tau + 0.5*a_x[0][i]*tau**2
        r_y[0][i+1] = r_y[0][i] + v_y[0][i]*tau + 0.5*a_y[0][i]*tau**2
        # conversione coordinate cartesiane in polari
        coo_p = polar(r_x[0][i+1] + r_y[0][i+1]*1j)
        r[0][i+1] = coo_p[0]
        phi[0][i+1] = coo_p[1]
        # calcolo dati sfruttando Velocity-Verlet
        # parte 2
        r_rel = sqrt((r_x[0][i+1] - r_x[1][i+1])**2 + (r_y[0][i+1] - r_y[1][i+1])**2)
        a_x[0][i+1] = -G*ms*r_x[0][i+1]/r[0][i+1]**3 - G*mp*(r_x[1][i+1] - r_x[0][i+1])/r_rel**3
        a_y[0][i+1] = -G*ms*r_y[0][i+1]/r[0][i+1]**3 - G*mp*(r_y[1][i+1] - r_y[0][i+1])/r_rel**3
        v_x[0][i+1] = v_x[0][i] + 0.5*tau*(a_x[0][i] + a_x[0][i+1])
        v_y[0][i+1] = v_y[0][i] + 0.5*tau*(a_y[0][i] + a_y[0][i+1])
        ##############################
        ##############################
        
        drawProgressBar(float(i)/Nstep, 60)
    print "\nFine calcolo dati."
            
    if data_on_file:

        # apre i file su cui scrivere i dati con la proprietà 'w'
        orbit_g = open('orbit_g.dat', 'w') # Giove
        orbit_t = open('orbit_t.dat', 'w') # Terra

        for j in range(1+1):#cond +1):
            if j == 0:
                planet = "Terra"
                filename = orbit_t                
            elif j == 1:
                planet = "Giove"
                filename = orbit_g
            print "\nScrittura dati %s su file:" % planet
            for i in range(Nstep):
                # scrive i dati su file
                # j = 0 > Terra
                # j = 1 > Giove
                filename.write(str(phi[j][i]))
                filename.write("\t")
                filename.write(str(r[j][i]))
                filename.write("\n")
                
                drawProgressBar(float(i)/Nstep, 60)
            print "\nFine scrittura dati %s su file." % planet
            
        # chiude i file su cui sono stati scritti i dati
        orbit_g.close()
        orbit_t.close()
    



print """Programme bla bla bla
1) Simulare orbita sistema Sole - Terra - pianeta
2) Simulare e calcolare orbita e parametri di un sistema Sole - pianeta
q) Per uscire
"""

while True:
    operation = raw_input("Scegliere l'operazione da eseguire: ")
    if operation == "1":
        print "Simulazione orbita sistema Sole - Terra - pianeta"
        ask_params()
        ask_data_on_file()
        orbit(1)
        break
    elif operation == "2":
        print "Simulazione orbita e calcolo parametri di un sistema Sole - pianeta"
        ask_params()
        ask_data_on_file()
        orbit(0)
        break
    else:
        print "Opzione non valida.\n"
        ans = True
