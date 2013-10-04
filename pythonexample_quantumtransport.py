from numpy import *
from numpy import dot, outer
from numpy.linalg import eig
from random import uniform
from pylab import *
#from scipy import *
import scipy.linalg as linalg

lol=open('stability/TEST6dots100MSele09DEPH005stabilitaA.txt','w')
#coordinates generation
a=1
n_dots=6
#gamma=0.1*1*2.0/6.52967771124
gamma=.05
delta=0.05
totti=100

##loading process:from file or random?

if a==0:			#generate random positions
	coord=[]
	coord+=[-1.0,-1.0,-1.0,1.0,1.0,1.0]

	for ii in range(3*(n_dots-2)):
		coord.append( float(uniform(-1.0,1.0)) )

if a==1:			#load input position file
	data=[map(float,iddu.split()) for iddu in open('coordinates/6dots100MSele09.txt').readlines()]

	
if a==3:
	data=[map(float,iddu.split()) for iddu in open('../matchingstructures/clusterTRY/cl19str9.txt').readlines()]
	foe=[1,2,8,9,10]
	coord=[]
	for ii in range(len(foe)):
		for ii2 in range(3):
			coord.append( data[0][(foe[ii]-1)*3+ii2] )

if a==2:
	coord=[-1.0,0.0,0.0,1.0,0.0,0.0,-0.5,0.0,0.0,+0.5,0.0,0.0,0.,0.35,0.0,0.0,-0.35,0.0]





##matrices calculation

s=zeros((n_dots**2,n_dots,n_dots),complex)

k=0
for i in range(n_dots):
	for j in range(n_dots-1-i): 
		s[k][i][j+i+1]=1.0/sqrt(2)
		s[k][j+i+1][i]=1.0/sqrt(2)
		k+=1

for i in range(n_dots):
	for jj in range(n_dots-1-i):
		s[k][i][jj+i+1]=1.j/sqrt(2)
		s[k][jj+i+1][i]=-1.j/sqrt(2)
		k+=1

for i in range(n_dots):
	s[k][i][i]=1.0
	k+=1


sigmaz=zeros((n_dots,n_dots,n_dots),complex)
for i in range(n_dots):
	for j in range(n_dots):
		sigmaz[i][j][j]=1
	sigmaz[i][i][i]=-1


#print s[-1][:][:],s[-2][:][:],s[-3][:][:],s[-4][:][:],s[-5][:][:],s[-6][:][:]

######################## beginning the real calculation

for ii00 in range(len(data)):

	coord=( [float(data[ii00][ii]) for ii in range(3*n_dots)] )
	coord0=( [float(data[ii00][ii]) for ii in range(3*n_dots)] )

	effii=[]
	for iii2 in range(totti):

		for iii3 in range((n_dots)*3):
			coord[iii3]=coord0[iii3]+delta*float(uniform(-1.0,1.0))


		dist=zeros((n_dots,n_dots),complex)
		for ii in range(n_dots):
			for ii2 in range(int(n_dots-ii-1)):
				dist[ii][ii+ii2+1]=sqrt((coord[ii*3]-coord[(ii+ii2+1)*3])**2+(coord[ii*3+1]-coord[(ii+ii2+1)*3+1])**2+(coord[ii*3+2]-coord[(ii+ii2+1)*3+2])**2)
	
#lambda matrix
		lambd=zeros((n_dots,n_dots),complex)
		alpha=1.0
		for ii in range(n_dots):
			for ii2 in range(int(n_dots-ii-1)):
				lambd[ii][ii+ii2+1]=alpha/float(dist[ii][ii+ii2+1]**3)
				lambd[ii+ii2+1][ii]=alpha/float(dist[ii][ii+ii2+1]**3)
#---------------------

		for ii0 in range(1):
			basis=eye(n_dots)		#projections on states
			HH=zeros((n_dots,n_dots),complex)	#hamiltonian

			HH+=sum(lambd[ii][ii2]*outer(basis[ii:ii+1],basis[ii2:ii2+1]) for ii in range(n_dots) for ii2 in range(n_dots))		#sinthetic way to do the previous 3 lines

			evoluto=zeros((n_dots**2,n_dots,n_dots),complex)
			for jj in range(n_dots**2):
				plop=zeros((n_dots,n_dots),complex)
				for ii in range(n_dots):
					plop+=dot(dot(sigmaz[ii],s[jj]),sigmaz[ii])-s[jj]
				evoluto[jj]=-1.j*(dot(HH,s[jj])-dot(s[jj],HH))+gamma*(plop)

			LL=zeros((n_dots**2,n_dots**2),complex)
			for ii in range(n_dots**2):
				for jj in range(n_dots**2):
					LL[jj][ii]=trace(dot(evoluto[jj],s[ii]))

#---------------------- diagonalization
			eigval2=eig(LL)[0]
			eigvecaa2=eig(LL)[1]	#the eigenvectors in columns
			eigvec2=[]		#the eigenvectors in lines

			for ii in range(n_dots**2):
				eigvec2.append( eigvecaa2[:,ii] )

#		print eigvec2

#---------------------- interaction matrix

			coeffi=zeros((n_dots**2,n_dots**2),complex)
		
			for ii in range(n_dots**2):
				for ii2 in range(n_dots**2):
					coeffi[ii][ii2]=dot(eigvec2[ii],eigvec2[ii2])
			coeffi2=linalg.inv(coeffi)

#		print coeffi

			basis2=eye(n_dots**2)
			inter2=zeros((n_dots**2,n_dots**2),complex)		
			for ii in range(n_dots**2):
				for ii2 in range(n_dots**2):
#				print (n_dots**2-n_dots),(n_dots**2-n_dots+1)
					inter2[ii][:]+=coeffi2[ii][ii2]*eigvec2[ii][:]*dot(eigvec2[ii2],basis2[(n_dots**2-n_dots):(n_dots**2-n_dots+1)][0])
	
#		print inter2
#		print inter2

#----------------------dynamics
			taum=0.05*pi*dist[0][1]**3	#tau

#		print taum

			length=500
			efft2=[[] for ii in range(n_dots)]

			polipo=zeros((n_dots,n_dots,n_dots),complex)	#6matrices, 6x6 with 1 in the nth diagonal place
			for ii in range(n_dots):		#projectors for the various sites
				polipo[ii][ii][ii]=1.0
		
			for ii in range(length):
				tt=float(taum/length*ii)
#			for ii2 in range(n_dots):
				for ii2 in range(1):
					aux=zeros(n_dots**2,complex)
					for ii3 in range(n_dots**2):
						aux[:]+=exp(1*eigval2[ii3]*tt)*inter2[ii3][:]
					rhoout=zeros((n_dots,n_dots),complex)
					for ii3 in range(n_dots**2):
						rhoout+=aux[ii3]*s[ii3]
#				efft2[ii2].append( float(trace(dot(rhoout,polipo[ii2]))) )
					efft2[1].append( float(trace(dot(rhoout,polipo[1]))) )

#			print aux[:]
#			print rhoout-rhoout.conjugate().transpose()
			eff2=max(efft2[1])
			effii.append( eff2 )
#		print eff2

#		lol.write('%s\t%s\n'%(ii00+1,eff2))
	lol.write('%s\t%f\t%f\n'%(ii00+1,mean(effii),std(effii)))
#		for ii in range(n_dots):
#		print rhoout-rhoout.conjugate().transpose()#rhoout.conjugate().transpose()[ii][ii]#, rhoout-rhoout.conjugate().transpose()
#		print trace(rhoout)

#		summo=0
#		for ii in range(n_dots):
#			summo+=efft2[ii][200]
#		print summo, trace(rhoout)

#		for ii in range(n_dots):
#			plot(efft2[ii], '-',label='%s'%(ii+1))
#	xlabel('Time')
#	ylabel('Excit')
#	ylim([0,1])
#		legend()
#		grid(1)
#		show()
