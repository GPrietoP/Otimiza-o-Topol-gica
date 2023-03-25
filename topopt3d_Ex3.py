## Otimização topológica para estruturas em 3D
from __future__ import division
import numpy as np
from numpy.core.fromnumeric import size
from scipy.sparse import coo_matrix
import scipy.sparse.linalg as ssl
from matplotlib import colors
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

# Programa principal
def main(nelx,nely,nelz,volfrac,penal,rmin,ft):	
	print("Problema de Otimização Topológica")
	print("Malha: " + str(nelx) + " x " + str(nely) + " x " + str(nelz))
	print("volfrac: " + str(volfrac) + ", rmin: " + str(rmin) + ", penal: " + str(penal))
	print("Método de filtragem: " + ["Sensibilidade","Densidade"][ft])
	ini1 = time.time()
	# Parâmetros pré-definidos para o looping
	maxloop = 200		# Número máximo de iterações
	tolx = 0.01			# Critério de parada
	# Propriedade do material
	E0 = 1.0         # Módulo de Young do material sólido
	Emin = 1e-9      # Módulo de Young para as regiões em vazio
	nu = 0.3         # Coeficiente do Poisson
	# Definindo o local de aplicação da carga
	il = np.array([nelx, nely]); jl = np.array([0, nely]); kl = np.array([nelz/2, nelz/2])          # Coordenadas
	loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl) 						# IDs dos nós
	loaddof = 3*loadnid[:] - 2                             						# Graus de liberdade
	# Definindo os graus de liberdade fixos
	iif, jf, kf = np.meshgrid(0, np.linspace(0,nely,num=nely+1), np.linspace(0,nelz,num=nelz+1))        # Coordenadas
	fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf) 											# IDs dos nós
	fixeddof = np.array([[3*fixednid[:]-1], [3*fixednid[:]-2], [3*fixednid[:]-3]])						# Graus de liberdade
	# PREPARANDO A ANÁLISE DE ELEMENTOS FINITOS
	nele = nelx*nely*nelz 					# Número total de elementos
	ndof = 3*(nelx+1)*(nely+1)*(nelz+1) 	# Número total de graus de liberdade
	F = np.zeros((ndof,2))
	n1 = int(loaddof[0]); n2 = int(loaddof[1])
	F[n1,0] = -1; F[n2,1] = 1      # Definindo a intensidade da força                         
	U = np.zeros((ndof,2))	
	VetDof = np.arange(ndof)	
	freedofs = np.setdiff1d(VetDof,fixeddof)	
	KE = lk_H8(nu)  
	edofMat=np.zeros((nelx*nely*nelz,24),dtype=int)
	for elz in range(nelz):
		for elx in range(nelx):
			for ely in range(nely):
				el = ely+elx*nely+elz*nely*nelx # De cima para baixo, de trás para frente				
				n1 = el+2+elx+(nelx+nely+1)*elz
				n2 = n1+(nely+1)
				n3 = n1+nely
				n4 = n1-1
				n5 = n1+(nelx+1)*(nely+1)
				n6 = n2+(nelx+1)*(nely+1)
				n7 = n3+(nelx+1)*(nely+1)
				n8 = n4+(nelx+1)*(nely+1)
				edofMat[el,:]=np.array([3*n1-2, 3*n1-1, 3*n1, 3*n2-2, 3*n2-1, 3*n2, 3*n3-2, 3*n3-1, 3*n3, 3*n4-2, 3*n4-1, 3*n4,
				 3*n5-2, 3*n5-1, 3*n5, 3*n6-2, 3*n6-1, 3*n6, 3*n7-2, 3*n7-1, 3*n7, 3*n8-2, 3*n8-1, 3*n8])-1                      # Matriz que armazena os graus de liberdade para cada elemento
	iK = np.kron(edofMat,np.ones((24,1))).flatten()
	jK = np.kron(edofMat,np.ones((1,24))).flatten()  
	#  PREPARAÇÃO DO FILTRO
	nfilter=int(nele*((2*(np.ceil(rmin)-1)+1)**3))
	iH = np.zeros(nfilter)	
	jH = np.zeros(nfilter)
	sH = np.zeros(nfilter)	
	cc=0
	for k1 in range(nelz):
		for i1 in range(nelx):
			for j1 in range(nely):
				row=k1*nelx*nely+i1*nely+j1
				kk1=int(np.maximum(k1-(np.ceil(rmin)-1),0))
				kk2=int(np.minimum(k1+np.ceil(rmin),nelz))
				ii1=int(np.maximum(i1-(np.ceil(rmin)-1),0))
				ii2=int(np.minimum(i1+np.ceil(rmin),nelx))
				jj1=int(np.maximum(j1-(np.ceil(rmin)-1),0))
				jj2=int(np.minimum(j1+np.ceil(rmin),nely))
				for k2 in range(kk1,kk2):
					for i2 in range(ii1,ii2):
						for j2 in range(jj1,jj2):
							col=k2*nelx*nely+i2*nely+j2
							fac=rmin-np.sqrt((i1-i2)*(i1-i2)+(j1-j2)*(j1-j2)+(k1-k2)*(k1-k2)) # raio minimo - a distância de centro a centro até o elemento "e"
							iH[cc]=row
							jH[cc]=col
							sH[cc]=np.maximum(0.0,fac)
							cc=cc+1
	# Finalização da montagem e a conversão para o "Compressed Sparse Column Format"
	H=coo_matrix((sH,(iH,jH)),shape=(nele,nele)).tocsc()	
	Hs=H.sum(1)

	# Inicializar e alocar variáveis
	x=volfrac * np.ones(nele,dtype=float)
	xold=x.copy()
	xPhys=x.copy()
	g=0
	dc=np.zeros((nely,nelx,nelz), dtype=float)

	# inicialização da iteração
	loop=0
	change=1
	dv = np.ones(nele)
	dc = np.ones(nele)
	ce = np.ones(nele)   
	while change>tolx and loop<maxloop:
		loop=loop+1		
		# Configurar o problema de elementos finitos
		sK=((KE.flatten()[np.newaxis]).T*(Emin+(xPhys)**penal*(E0-Emin))).flatten(order='F')
		K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()		
		# Remover deslocamentos restritos da matriz
		K = K[freedofs,:][:,freedofs]										
		# Solução do sistema
		ini2 = time.time()										
		U[freedofs,0], iter = ssl.cg(K,F[freedofs,0], x0=None, tol=1e-8, maxiter=8000, M=None)
		U[freedofs,1], iter = ssl.cg(K,F[freedofs,1], x0=None, tol=1e-8, maxiter=8000, M=None)	
		fim2 = time.time()
		print(iter)
		# Função Objetivo e Sensibilidade
		dc = np.zeros(nele)
		obj = 0
		for i in range(np.size(F,1)):			
			Ue = U[:,i]
			A = Ue[edofMat].reshape(nele,24)			
			ce[:] = (np.dot(A,KE)*A).sum(1)		
			obj = obj + ((Emin+xPhys**penal*(E0-Emin))*ce).sum().sum().sum()
			dc[:] = dc[:] - (-penal*xPhys**(penal-1)*(E0-Emin))*ce			
		dv[:] = np.ones(nele)   
		# Filtro Sensibilidade
		if ft==0:
			dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
		# Filtro Densidade
		elif ft==1:
			dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
			dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]
		# Critério de Optimalidade
		xold[:]=x
		(x[:],g)=oc(nelx,nely,nelz,x,volfrac,dc,dv,g)
		# Filtrando as variáveis de projeto
		if ft==0:   xPhys[:]=x
		elif ft==1:	xPhys[:]=np.asarray(H*x[np.newaxis].T/Hs)[:,0]
		# Calcular a mudança através de (valor de densidade antigo - o valor atual)
		change=np.linalg.norm(x.reshape(nele,1)-xold.reshape(nele,1),np.inf)
		# Printar o histórico de iteração na tela
		print("it.: {0} , obj.: {1:.3f} Vol.: {2:.3f}, ch.: {3:.3f}".format(\
					loop,obj,(g+volfrac*nele)/(nele),change))
		print(fim2-ini2)		
	fim1 = time.time()
	print(fim1-ini1)

	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	max_scale=max(nelx, nelz, nely)
	scale_x=nelx/max_scale
	scale_y=nelz/max_scale
	scale_z=nely/max_scale
	ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 0.5]))

	x,y,z = np.indices((nelx,nelz,nely))
	total = (x==nelx) & (y==nelz) & (z==nely)
	colors = np.zeros(total.shape + (3,))
	cont = 0
	for i in range(nelz):
		for j in range(nelx):
			for k in range(nely):
				if xPhys[cont]>0.5:
					cube = (x==j) & (y==i) & (z==k)
					colors[cube] = 0.2+0.8*(1-xPhys[cont])
					total = total | cube										
				cont = cont+1
	ax.voxels(total, facecolors=colors)
	plt.axis('off')
	plt.show()
	
#  Função que gera a matriz de rigidez de um elemento
def lk_H8(nu):
	A = np.matrix([[32, 6, -8, 6, -6, 4, 3, -6, -10, 3, -3, -3, -4, -8],
	[-48, 0, 0, -24, 24, 0, 0, 0, 12, -12, 0, 12, 12, 12]])

	k = 1/144*A.T*np.matrix([[1], [nu]])
	
	K1 = np.array([ [k[0,0], k[1,0], k[1,0], k[2,0], k[4,0], k[4,0]],
	[k[1,0], k[0,0], k[1,0], k[3,0], k[5,0], k[6,0]],
	[k[1,0], k[1,0], k[0,0], k[3,0], k[6,0], k[5,0]],
	[k[2,0], k[3,0], k[3,0], k[0,0], k[7,0], k[7,0]],
	[k[4,0], k[5,0], k[6,0], k[7,0], k[0,0], k[1,0]],
	[k[4,0], k[6,0], k[5,0], k[7,0], k[1,0], k[0,0]] ])   

	K2 = np.array([ [k[8,0],  k[7,0],  k[11,0], k[5,0],  k[3,0],  k[6,0]],
	[k[7,0],  k[8,0],  k[11,0], k[4,0],  k[2,0],  k[4,0]],
	[k[9,0], k[9,0], k[12,0], k[6,0],  k[3,0],  k[5,0]],
	[k[5,0],  k[4,0],  k[10,0], k[8,0],  k[1,0],  k[9,0]],
	[k[3,0],  k[2,0],  k[4,0],  k[1,0],  k[8,0],  k[11,0]],
	[k[10,0], k[3,0],  k[5,0], k[11,0],  k[9,0], k[12,0]] ])   

	K3 = np.array([ [k[5,0],  k[6,0],  k[3,0],  k[8,0],  k[11,0], k[7,0]],
	[k[6,0],  k[5,0],  k[3,0],  k[9,0], k[12,0], k[9,0]],
	[k[4,0],  k[4,0],  k[2,0],  k[7,0],  k[11,0], k[8,0]],
	[k[8,0],  k[9,0], k[1,0],  k[5,0],  k[10,0], k[4,0]],
	[k[11,0], k[12,0], k[9,0], k[10,0], k[5,0],  k[3,0]],
	[k[1,0],  k[11,0], k[8,0],  k[3,0],  k[4,0],  k[2,0]] ])  

	K4 = np.array([ [k[13,0], k[10,0], k[10,0], k[12,0], k[9,0], k[9,0]],
	[k[10,0], k[13,0], k[10,0], k[11,0], k[8,0],  k[7,0]],
	[k[10,0], k[10,0], k[13,0], k[11,0], k[7,0],  k[8,0]],
	[k[12,0], k[11,0], k[11,0], k[13,0], k[6,0],  k[6,0]],
	[k[9,0], k[8,0],  k[7,0],  k[6,0],  k[13,0], k[10,0]],
	[k[9,0], k[7,0],  k[8,0],  k[6,0],  k[10,0], k[13,0]] ]) 

	K5 = np.array([ [k[0,0], k[1,0],  k[7,0],  k[2,0], k[4,0],  k[3,0]],
	[k[1,0], k[0,0],  k[7,0],  k[3,0], k[5,0],  k[10,0]],
	[k[7,0], k[7,0],  k[0,0],  k[4,0], k[10,0], k[5,0]],
	[k[2,0], k[3,0],  k[4,0],  k[0,0], k[7,0],  k[1,0]],
	[k[4,0], k[5,0],  k[10,0], k[7,0], k[0,0],  k[7,0]],
	[k[3,0], k[10,0], k[5,0],  k[1,0], k[7,0],  k[0,0]] ])

	K6 = np.array([ [k[13,0], k[10,0], k[6,0],  k[12,0], k[9,0], k[11,0]],
	[k[10,0], k[13,0], k[6,0],  k[11,0], k[8,0],  k[1,0]],
	[k[6,0],  k[6,0],  k[13,0], k[9,0], k[1,0],  k[8,0]],
	[k[12,0], k[11,0], k[9,0], k[13,0], k[6,0],  k[10,0]],
	[k[9,0], k[8,0],  k[1,0],  k[6,0],  k[13,0], k[6,0]],
	[k[11,0], k[1,0],  k[8,0],  k[10,0], k[6,0],  k[13,0]] ])

	KE = 1/((nu+1)*(1-2*nu))*np.concatenate((np.concatenate((K1,K2,K3,K4),axis=1), np.concatenate((K2.T,K5,K6,K3.T),axis=1), np.concatenate((K3.T,K6,K5.T,K2.T),axis=1), np.concatenate((K4,K3,K2,K1.T),axis=1)), axis=0)

	return (KE)

# Função critério de Optimalidade
def oc(nelx,nely,nelz,x,volfrac,dc,dv,g):
	l1=0
	l2=1e9
	move=0.2
	# Remodelar para realizar operções vetoriais
	xnew=np.zeros(nelx*nely*nelz)
	while (l2-l1)/(l1+l2)>1e-3:
		lmid=0.5*(l2+l1)
		xnew[:]= np.maximum(0.0,np.maximum(x-move,np.minimum(1.0,np.minimum(x+move,x*np.sqrt(-dc/dv/lmid)))))
		gt=g+np.sum((dv*(xnew-x)))
		if gt>0 :
			l1=lmid
		else:
			l2=lmid
	return (xnew,gt)

	# The real main driver    
if __name__ == "__main__":
	# Default input parameters
	nelx = 30 #60
	nely = 30 #20
	nelz = 2 #4
	volfrac = 0.4
	rmin = 1.5
	penal = 3.0
	ft = 1 # ft==0 -> sens, ft==1 -> dens
	import sys
	if len(sys.argv)>1: nelx   =int(sys.argv[1])
	if len(sys.argv)>2: nely   =int(sys.argv[2])
	if len(sys.argv)>3: nelz=float(sys.argv[3])
	if len(sys.argv)>4: volfrac=float(sys.argv[4])
	if len(sys.argv)>5: rmin   =float(sys.argv[5])
	if len(sys.argv)>6: penal  =float(sys.argv[6])
	if len(sys.argv)>7: ft     =int(sys.argv[7])
	main(nelx,nely,nelz,volfrac,penal,rmin,ft)