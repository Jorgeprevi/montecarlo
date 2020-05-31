import numpy as np
import matplotlib.pyplot as plt
import random

pasos=100                      # número de pasos montecarlo
N=150                            # número de POSICIONES????????
nr=50                            # número de pasos rechazados entre cada medida
pterm=200                        # número de pasos de termalización
x=np.zeros(N, dtype=float)       # vector de ???????
xpre=np.zeros(N, dtype=float)
a=0.1                            # parámetro de muestreo
lam=0.0
m=1.0
w=1.0
muestra=np.zeros(0, dtype=float)

def pbc(num, L):
    '''
    Función que implementa condiciones de contorno periódicas
    '''
    if num>L:
        return num-L-1 
    elif num<0:
        return L+num+1
    else:
        return num

def accion(x):
    '''
    Función que devuelve la acción a partir del vector x, m y omega
    '''
    S=0.0
    for i in range(N):
        #S+=0.5*m*(x[pbc(i+1,N-1)]-x[i])**2.0/a+(0.5*m*w*w*x[i]**2.0+lam*x[i]**4.0)/a
        S+=0.5*m*((x[pbc(i+1,N-1)]-x[i])**2.0/a+w*w*x[i]**2.0)*a
    return S

def acepto(x0, x1):
    '''
    Función booleana que decide si se acepta un paso a partir del factor de Boltzmann
    '''
    S0=accion(x0)
    S1=accion(x1)
    if random.uniform(0,1)<(min(1.0,np.exp(S0-S1))):
        return True 
    else: 
        return False

# se inicia el proceso de termalización
for i in range(pterm*N):
    x[random.randint(0, N-1)]+=random.uniform(-a,a)
    if acepto(xpre, x)==True:
        xpre=np.copy(x)
    else:
        x=np.copy(xpre)
    # if i%(N*nr)==0:
    #     plt.plot(x)
    #     plt.show()
    

# se inicia el proceso principal de medidas
xp2=0.0
count=0
for i in range(pasos*N):
    x[random.randint(0, N-1)]+=random.uniform(-a,a)
    if acepto(xpre, x)==True:
        xpre=np.copy(x)
    else:
        x=np.copy(xpre)
    if i%(N*nr)==0:
        for j in range(N):
            xp2+=x[j]*x[j]/float(N)
        muestra=np.append(muestra, x)
        count+=1
        # plt.plot(x)
        # plt.show()
xp2/=float(count)
energia=m*w*w*xp2

# ploteo resultados
pesos=np.ones_like(muestra)/float(len(muestra))
plt.hist(muestra, weights=pesos, bins=50) 
plt.title("funcion de onda") 
plt.show()
print(count)
print(np.amax(x))
print("Energía: %s" % (energia))