
# coding: utf-8

# In[1]:

class vortexParticle:
    def __init__(self, z=0.0+0.0j, gamma=0.0, r=0.1):
        self.z=z
        self.gamma=gamma
        self.v=0.0+0.0j
        self.r=r
    def Set_z(self,z):
        self.z=z
    def Set_gamma(self,gamma):
        self.gamma=gamma
    def vel(self,z):
        if(abs(z-self.z)<self.r):
            return (self.gamma)*(z-self.z)*(1j)/(2*3.14159*(self.r)*(self.r))
        else:
            return (self.gamma)*(1j)/(2*3.14159*(z-self.z).conjugate())


# In[3]:

import numpy as np


# In[3]:

with open('VortexMerging_Wall+Doublet.txt') as f:
    Data = []
    for line in f:
        line = line.split()
        if line:
            line = [float(x) for x in line]
            Data.append(line)
n=Data[0][0]
dt=Data[0][1]
nSteps=Data[0][2]
i=1
vp=[]
while i<n+1:
    vp.append(Data[i])
    i+=1


# In[4]:

n=int(n)
nSteps=int(nSteps)


# In[2]:

Patch=np.empty((nSteps,2*n),dtype=object)


# In[1]:

Patch[1][1]


# In[6]:

i=0
a=0.1
while i<n:
    Patch[0,i]=vortexParticle(vp[i][1]+vp[i][2]*1j,vp[i][0],a)
    Patch[0,i+n]=vortexParticle(Patch[0,i].z.conjugate(),vp[i][0],a)
    i+=1


# In[7]:

print "Computaion begins"
t=0
while t<nSteps-1:
    i=0
    while i<n:
        j=0
        while j<2*n:
            Patch[t,i].v+=Patch[t,j].vel(Patch[t,i].z)
            j+=1
        Patch[t+1,i]=vortexParticle(Patch[t,i].z+Patch[t,i].v*dt,Patch[t,i].gamma,Patch[t,i].r)
        if(Patch[t+1,i].z.imag<0):
            Patch[t+1,i].z=Patch[t+1,i].z.conjugate()
            Patch[t+1,i].v=Patch[t+1,i].v.conjugate()
        i+=1
    while i<2*n:
        Patch[t+1,i]=vortexParticle((Patch[t+1,i-n].z).conjugate(),-Patch[t+1,i-n].gamma,Patch[t+1,i-n].r)
        i+=1
    t+=1
print "Computation ends"


# In[ ]:

print "Entering data to a text file"
f=open('./data.txt','w')
print >> f, "Number of time steps is:",nSteps,"and number of particles is:\n",n
t=0
while t<nSteps:
    i=0
    while i<n:
        print >> f, Patch[t,i].z,"\t"
        i+=1
    print >> f, "\n"
    t+=1


# In[ ]:

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt


# In[ ]:

import math


# In[ ]:

print "Variances, CoMs and Corelation"
t=0
z1_bar=np.empty(nSteps,dtype=complex)
z2_bar=np.empty(nSteps,dtype=complex)
var1=np.zeros(nSteps)
var2=np.zeros(nSteps)
var=np.zeros(nSteps)
while t<nSteps:
    z1_bar[t]=0
    z2_bar[t]=0
    var1[t]=0
    var2[t]=0
    var[t]=0
    i=0
    while i<n/2:
        z1_bar[t]+=Patch[t,i].z
        i+=1
    z1_bar[t]/=n/2
    while i<n:
        z2_bar[t]+=Patch[t,i].z
        i+=1
    z2_bar[t]/=n/2
    i=0
    while i<n/2:
        var1[t]+=math.pow(abs(Patch[t,i].z-z1_bar[t]),2)
        var[t]+=math.pow(abs(Patch[t,i].z-z1_bar[t]),2)
        i+=1
    var1[t]/=n/2
    while i<n:
        var2[t]+=math.pow(abs(Patch[t,i].z-z2_bar[t]),2)
        var[t]+=math.pow(abs(Patch[t,i].z-z2_bar[t]),2)
        i+=1
    var2[t]/=n/2
    var[t]/=n
    t+=1

fig=plt.figure()
plt.scatter(z1_bar.real,z1_bar.imag)
fig.savefig("center_traj_1.png")

fig=plt.figure()
plt.scatter(z2_bar.real,z2_bar.imag)
fig.savefig("center_traj_2.png")

fig=plt.figure()
plt.plot(var1)
fig.savefig("Variance1.png")

fig=plt.figure()
plt.plot(var2)
fig.savefig("Variance2.png")


# In[ ]:

print "Vortex Momemtum"
I=np.zeros(nSteps,dtype=complex)
t=0
while t<nSteps:
    I[t]=0
    i=0
    while i<n:
        I[t]+=-1j*Patch[t,i].z*Patch[t,i].gamma
        i+=1
    t+=1

fig=plt.figure()    
plt.plot(abs(I))
fig.savefig("VortexMomentum.png")


# In[ ]:

print "Scatter Plot"
t=0
while t<nSteps:
    i=0
    fig=plt.figure()
    plt.xlim(-5,10)
    plt.ylim(-5,10)
    while i<n/2:
        plt.scatter(Patch[t,i].z.real,Patch[t,i].z.imag,color='r')
        i+=1
    while i<n:
        plt.scatter(Patch[t,i].z.real,Patch[t,i].z.imag,color='b')
        i+=1
    s="VortexMerge_Wall+Doublet_%d.png"%t 
    fig.savefig(s)
    plt.close(fig)
    t+=nSteps/100

