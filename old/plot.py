import numpy as np
from mayavi import mlab
phafile=open('pha.dat','r')
pha=phafile.readlines()
a=[]
b=[]
for item in pha:
    if not item == '\n':
        list=item.split(' ')
        a.append(float(list[0]))
        b.append(float(list[1]))
phafile.close()

n=1

Nx = 100
Ny = 100
Nz = 1

a=np.array(a).reshape(int(Nx),int(Ny),int(Nz))
b=np.array(b).reshape(int(Nx),int(Ny),int(Nz))

amatrix=np.zeros((int(Nx)*n,int(Ny)*n,int(Nz)*n))
bmatrix=np.zeros((int(Nx)*n,int(Ny)*n,int(Nz)*n))
for i in range(int(Nx)*n):
    for j in range(int(Ny)*n):
        for k in range(int(Nz)*n):
            amatrix[i][j][k]=a[i%int(Nx)][j%int(Ny)][k%int(Nz)]
            bmatrix[i][j][k]=b[i%int(Nx)][j%int(Ny)][k%int(Nz)]
            

src = mlab.pipeline.scalar_field(bmatrix)
mlab.pipeline.iso_surface(src, opacity=1,contours=[0.5],color=(1,0,0))

src = mlab.pipeline.scalar_field(amatrix)
mlab.pipeline.iso_surface(src, opacity=1,contours=[0.5],color=(0,1,0))