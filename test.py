#Test if in the limit of high A_z we get Ising again

import random
import matplotlib.pyplot as plt
import math
import numpy as np
import time

start=time.time()

L= 50        #Size of lattice
N= 1000    #Number of MC steps
J= 1         #Exchage coupling
A_z= 0.4       #Anisotropic coupling
N_eq=100   #Steps for MC to equilibrate

fh=open('output.txt','w')

fh.write("High Anisotrpoy Limit of Heisenberg Hamiltonian")
fh.write('\n')
fh.write("Length= " + str(L) )
fh.write('\n')
fh.write('\n')
fh.flush()
  
TEMP=[]
MAG_avg=[]

def func(x):
    return(x%L)

def add(x,y):              #Adds two arrays
    x_temp=x[0]+y[0]
    y_temp=x[1]+y[1]
    z_temp=x[2]+y[2]
    return(x_temp,y_temp,z_temp)
    
def sub(x,y):
    x_temp=x[0]-y[0]
    y_temp=x[1]-y[1]
    z_temp=x[2]-y[2]
    return(x_temp,y_temp,z_temp)    

def div(x,y):              #Divides array x by number y
    x_temp=x[0]/y
    y_temp=x[1]/y
    z_temp=x[2]/y
    return(x_temp,y_temp,z_temp)
    
def absolute(x):
    x_temp=abs(x[0])
    y_temp=abs(x[1])
    z_temp=abs(x[2])
    return(x_temp,y_temp,z_temp)
    

#Now we look over temperatures as well

for t in range(0,11):
    kT=0.4 + t/10
    fh.write("Temperature=" + str(kT))
    fh.write('\n')
    spin=[]                         #At each T we initalise the spins
    for i in range(0,L):
        temp=[]
        j=0
        while(j<L):
            x=random.uniform(-1,1)        #Now instead of up or down we have a vector 
            y=random.uniform(-1,1)
            z=random.uniform(-1,1)
            s=(x**2+y**2+z**2)**0.5
            if(s<1):
                spin_site=[(x/s),(y/s),(z/s)]
                j=j+1
                temp.append(spin_site)
        spin.append(temp)
    fh.write("Spin orientation initialised")
    fh.write('\n')

    mag=[0,0,0]   
    energy=0                        #And compute avergae magnetisation and energy per site
    for i in range (0,L):              
        for j in range(0,L):
            mag=add(mag,spin[i][j])
            energy= energy - 0.5*J*((np.dot(spin[i][j],spin[func(i+1)][j]))+(np.dot(spin[i][j],spin[func(i-1)][j]))+(np.dot(spin[i][j],spin[i][func(j+1)]))+(np.dot(spin[i][j],spin[i][func(j-1)]))) - A_z*(spin[i][j][2])**2 
    fh.write('Ini_mag='+ str(div(mag,L**2)))
    fh.write('\n')
    fh.write('Ini_en=' + str(energy/L**2))
    fh.write('\n')
    fh.write('\n')
    fh.write('Monte Carlo steps begin')
    fh.write('\n')

    MAG=[]                          #MAG stores the magnetisation over a MC run at a fixed temperature
    for k in range(0,N):
        if(k%100==0):
            fh.write('Steps: ' + str(k))
            fh.write('\n')
            fh.flush()
        for a in range(0,L):
            for b in range(0,L):
                i=random.randint(0,L-1)
                j=random.randint(0,L-1)
                E_ini= -J*((np.dot(spin[i][j],spin[func(i+1)][j]))+(np.dot(spin[i][j],spin[func(i-1)][j]))+(np.dot(spin[i][j],spin[i][func(j+1)]))+(np.dot(spin[i][j],spin[i][func(j-1)]))) - A_z*(spin[i][j][2])**2
                spin_temp=spin[i][j]
                o=0
                while(o<1):
                    x=random.uniform(-1,1)
                    y=random.uniform(-1,1)
                    z=random.uniform(-1,1)
                    s=(x**2+y**2+z**2)**0.5
                    if(s<1):
                        spin[i][j]=[(x/s),(y/s),(z/s)]
                        o=o+1
                E_fin= -J*((np.dot(spin[i][j],spin[func(i+1)][j]))+(np.dot(spin[i][j],spin[func(i-1)][j]))+(np.dot(spin[i][j],spin[i][func(j+1)]))+(np.dot(spin[i][j],spin[i][func(j-1)]))) - A_z*(spin[i][j][2])**2
                dE=E_fin-E_ini
                dmag=sub(spin[i][j],spin_temp)
                if(dE<=0):
                    energy=energy+(dE)
                    mag=add(mag,dmag)
                else:
                    u=math.exp(-dE/kT)
                    h=np.random.uniform(0,1)
                    if(h<u):
                        energy=energy+(dE)
                        mag=add(mag,dmag)
                    else:
                        spin[i][j]=spin_temp
        if(k>N_eq):
            MAG.append(absolute(mag))
        
    MAG_temp=[0,0,0]
    for i in range(0,len(MAG)):             #Calculate and write average magnetisation (absolute) over N-N_eq steps of the MC run
        MAG_temp=add(MAG_temp,MAG[i])
    m_avg=div(MAG_temp,len(MAG))
    fh.write('\n')
    fh.write("Average Magnetisation" + str(m_avg) )
    fh.write('\n')
    fh.write('\n')
    fh.write('\n')
    MAG_avg.append(m_avg)
    TEMP.append(kT)

MAG_avg_tot=[]
for i in range(0,len(MAG_avg)):
    MAG_avg_tot.append(((MAG_avg[i][0]**2)+(MAG_avg[i][1]**2)+(MAG_avg[i][2]**2))**0.5)
    
#plot abs(M) vs kT

plt.plot(TEMP,MAG_avg_tot,'ro')
plt.xlabel("Temperature (kT)")
plt.ylabel("Absolute Mean Magetisation")
plt.savefig("Anisotropic Heisenberg.png")

end=time.time()
time=end-start

fh.write('Run Time=' + str(time) )

fh.close()
