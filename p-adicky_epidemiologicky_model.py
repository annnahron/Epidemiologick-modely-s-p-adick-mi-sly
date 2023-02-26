import numpy as np
import matplotlib.pyplot as plt

#základní parametry modelu
p = 2
pocetDeleni = 5
pocetClusteru = p**pocetDeleni
pocetPrvkuClusteru = 16

#převrácená hodnota infekční doby
beta = 1/10

maticePrechodu = np.zeros((pocetClusteru, pocetClusteru))

#definování p-adické vzdálenosti
def pVzdalenost(a, b):
    if (a==b):
        vzdalenost = 1/64
    else:
        rozdil = abs(a-b)
        k = pocetDeleni
        while k > -1:
            if (rozdil % (p**k) == 0):
                pvaluace = k
                k = -1
            else:
                k = k-1
        vzdalenost = p**(-pvaluace)
    return vzdalenost

#koeficienty p-adického modelu
alpha = 0.5
soucet = 0
for promenna in range(pocetClusteru):
    soucet = soucet + 1/pow((pVzdalenost(2*promenna, 2*0)), 1+alpha)
Calpha = 1/soucet

#vepsání hodnot do matice přechodu
for l in range(pocetClusteru):
    for m in range(pocetClusteru):
        maticePrechodu[l][m] = Calpha/(pow(pVzdalenost(2*l,2*m),(1+alpha)))
        

#počáteční podmínky
vPravdepodobnosti0 = np.zeros(pocetClusteru)
vPravdepodobnosti0[0] = 5/16
vPravdepodobnosti0[2] = 5/16

#diferenční rovnice
def deriv(vPravdepodobnosti, S, I, R, t, mPrechodu, beta):
    dvPravdepodobnosti = np.zeros(pocetClusteru)
    pravdepodobnostnakazycloveka = np.zeros(pocetClusteru)

    for q in range(pocetClusteru):
        kladnysoucet = 0
        zapornysoucet = 0
        for r in range(pocetClusteru):
            kladnysoucet = kladnysoucet + mPrechodu[r][q]*vPravdepodobnosti[r]
            zapornysoucet = zapornysoucet + mPrechodu[q][r]*vPravdepodobnosti[q]
        dvPravdepodobnosti[q] = kladnysoucet - zapornysoucet

    dSdt = np.zeros(pocetClusteru)
    dIdt = np.zeros(pocetClusteru)
    dRdt = np.zeros(pocetClusteru)

    for q in range(pocetClusteru):
        for r in range(pocetClusteru):
            pravdepodobnostnakazycloveka[q] = pravdepodobnostnakazycloveka[q] + maticePrechodu[q][r]*vPravdepodobnosti[r]
    
    for q in range(pocetClusteru):
        dSdt[q] = -pravdepodobnostnakazycloveka[q]* S[q]
        dIdt[q] = pravdepodobnostnakazycloveka[q]* S[q] - beta * I[q]
        dRdt[q] = beta * I[q]

    return dvPravdepodobnosti, dSdt, dIdt, dRdt


#počet kroků
tSteps = 400
tMax = 200
dt = tMax/tSteps

#matice na ukládání vypočítaných hodnot
evolMatrixVect = np.zeros((tSteps, pocetClusteru))
evolMatrixVect[0] = vPravdepodobnosti0
evolS = np.zeros((tSteps, pocetClusteru))
evolI = np.zeros((tSteps, pocetClusteru))
evolR = np.zeros((tSteps, pocetClusteru ))

for q in range(pocetClusteru):
    evolS[0][q] = pocetPrvkuClusteru

#počáteční podmínky
evolS[0][0] = 11
evolI[0][0] = pocetPrvkuClusteru - evolS[0][0]
evolS[0][2] = 11
evolI[0][2] = pocetPrvkuClusteru - evolS[0][2]

#vepisování hodnot do evol-matic
for i in range(tSteps-1):
    t = tMax/tSteps * i
    dvPravdepodobnosti, dSdt, dIdt, dRdt = deriv(evolMatrixVect[i, :], evolS[i, :], evolI[i, :], evolR[i, :], t, maticePrechodu, beta)
    evolMatrixVect[i+1, :] = evolMatrixVect[i, :] + dvPravdepodobnosti * dt
    evolS[i+1] = evolS[i] + dSdt * dt
    evolI[i+1] = evolI[i] + dIdt * dt
    evolR[i+1] = evolR[i] + dRdt * dt

#vykreslení grafů
plt.plot(np.linspace(0, tMax, tSteps), evolS[:,0]/pocetPrvkuClusteru*100, linestyle="-", label = 'S')
plt.plot(np.linspace(0, tMax, tSteps), evolI[:,0]/pocetPrvkuClusteru*100, linestyle="-", label = 'I')
plt.plot(np.linspace(0, tMax, tSteps), evolR[:,0]/pocetPrvkuClusteru*100, linestyle="-", label = 'R')
plt.title('Šíření epidemie v clusteru C_0')
plt.xlabel('Čas')
plt.ylabel('Počet studentů [%]')
plt.legend()
plt.show()

plt.plot(np.linspace(0, tMax, tSteps), evolS[:,16]/pocetPrvkuClusteru*100, linestyle="-", label = 'S')
plt.plot(np.linspace(0, tMax, tSteps), evolI[:,16]/pocetPrvkuClusteru*100, linestyle="-", label = 'I')
plt.plot(np.linspace(0, tMax, tSteps), evolR[:,16]/pocetPrvkuClusteru*100, linestyle="-", label = 'R')
plt.title('Šíření epidemie v clusteru C_16')
plt.xlabel('Čas')
plt.ylabel('Počet studentů [%]')
plt.legend()
plt.show()

plt.plot(np.linspace(0, tMax, tSteps), evolS[:,3]/pocetPrvkuClusteru*100, linestyle="-", label = 'S')
plt.plot(np.linspace(0, tMax, tSteps), evolI[:,3]/pocetPrvkuClusteru*100, linestyle="-", label = 'I')
plt.plot(np.linspace(0, tMax, tSteps), evolR[:,3]/pocetPrvkuClusteru*100, linestyle="-", label = 'R')
plt.title('Šíření epidemie v clusteru C_3')
plt.xlabel('Čas')
plt.ylabel('Počet studentů [%]')
plt.legend()
plt.show()


