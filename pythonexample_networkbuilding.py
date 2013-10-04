#from pylab import *
from numpy import *


a=2

if a==0:

#	data=[map(str,iddu) for iddu in open('tryrun1ALLco10.rnk').readlines()]

	deit=open('../RMSDclust-allruns-unbok-fitrealprot-2002-2006-co08.rnk')
	lollo=open('RMSDclust-allruns-unbok-fitrealprot-2002-2006-co08.mm','w')
	###################### insert 2 digits to increase space for RESID number

#	print len(data)


#	ntot=max(data)
#	print ntot
#	ntot=6497
	ntot=57834
#	matrix=zeros((ntot+1,ntot+1))

	aa=[]
	for ii in deit:
		oo=ii.split()
		aa.append( int(oo[0]) )


	hh={}
	for ii in range(len(aa)-1):
		try:
			hh['%d\t%d'%(aa[ii],aa[ii+1])]+=1
		except:
			hh['%d\t%d'%(aa[ii],aa[ii+1])]=1
	

	for ii in hh:
		lollo.write('%s\t%d\n'%(ii,hh[ii]))


if a==1:		#connectivity

	data=[map(int,iddu.split()) for iddu in open('RMSDclust-allruns-unbok-fitrealprot-2002-2006-co08.mm').readlines()]
#	print len(data), len(data[0])
	ntot=57834
	connect=zeros(ntot)
	
	for ii in range(len(data)):
		connect[data[ii][0]-1]+=1
		connect[data[ii][1]-1]+=1

	yy1,xx1=histogram(connect, bins=50, normed=True)	#range=(0.,0.6),
	
	for ii in range(len(yy1)):
		print (xx1[ii]+xx1[ii+1])*0.5,yy1[ii]	

if a==2:	#build the network only with RMSD <4.5

#	deit=open('../RMSDclust-allruns-unbok-fitrealprot-co08.rnk')
	deit=[map(int,iddu.split()) for iddu in open('RMSDclust-allruns-unbok-fitligreducted-2002-2006-co04.rnk').readlines()]
	deitR=[map(float,iddu.split()) for iddu in open('../rmsdTOT/rmsdfitprotTOTcut.txt').readlines()]
#	deitR=open('../rmsdTOT/rmsdfitprotTOTcut.txt')
	lollo=open('RMSDclust-allruns-unbok-fitligreducted-2002-2006-co04co45A.mm','w')
	lollo2=open('RMSDclust-allruns-unbok-fitligreducted-2002-2006-co04co45A.rank','w')
		
	tot=0
	hh={}
	pp={}

	iniz=0		#tells us if it's the beginning of a series of OK
	for ii in range(len(deitR)-1):
		if float(deitR[ii+1][1])<4.5 and float(deitR[ii][1])<4.5:
			tot+=1
			try:
				hh['%d\t%d'%(deit[ii][0],deit[ii+1][0])]+=1
			except:
				hh['%d\t%d'%(deit[ii][0],deit[ii+1][0])]=1

###########################################################################
###			population
###			if it's the first of a series
			if iniz==0:
				try:
					pp['%s'%(deit[ii][0])]+=1
				except:
					pp['%s'%(deit[ii][0])]=1
				try:
					pp['%s'%(deit[ii+1][0])]+=1
				except:
					pp['%s'%(deit[ii+1][0])]=1

###			if not
			else:
				try:
					pp['%s'%(deit[ii+1][0])]+=1
				except:
					pp['%s'%(deit[ii+1][0])]=1
###########################################################################
			iniz=1
		else:
			iniz=0

	globi=0
	for ii2 in pp:
		globi+=pp[ii2]

	print tot,len(pp), globi

	for ii in hh:
		lollo.write('%s\t%d\n'%(ii,hh[ii]))

	for ii2 in pp:
		lollo2.write('%s\t%d\n'%(ii2,pp[ii2]))

if a==3:	#controprova
	deit=[map(int,iddu.split()) for iddu in open('../RMSDclust-allruns-unbok-fitrealprot-2002-2006-co08.rnk').readlines()]
	deitR=[map(float,iddu.split()) for iddu in open('../rmsdTOT/rmsdfitprotTOTcut.txt').readlines()]
		
	middle=0
	globa=0
	compl=0
	dx=0
	sx=0

	for ii in range(len(deit)-2):
		if float(deitR[ii+1][1])<4.5 and float(deitR[ii][1])>4.5 and float(deitR[ii+2][1])>4.5:
			middle+=1
		if float(deitR[ii+1][1])<4.5:
			globa+=1
		if float(deitR[ii+1][1])<4.5 and float(deitR[ii][1])<4.5 and float(deitR[ii+2][1])<4.5:
			compl+=1
		if float(deitR[ii+1][1])<4.5 and float(deitR[ii][1])>4.5 and float(deitR[ii+2][1])<4.5:
			dx+=1
		if float(deitR[ii+1][1])<4.5 and float(deitR[ii][1])<4.5 and float(deitR[ii+2][1])>4.5:
			sx+=1

	print dx,sx,compl,middle,globa
