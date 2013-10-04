program q_dots

implicit none

TYPE array_pointer
INTEGER,DIMENSION(:),POINTER::p1
END TYPE array_pointer

integer::i,j,k,g,h,ii,jj,kk,gg,unit,frame,total_frames,total_dots,bestlabel,lo,lo2,angle,kkk,contttin
integer::iii1,iii2,iii3,iii4,iii5,iii6,iii7,iii8
integer::ttt,qqq

double precision,dimension(:),ALLOCATABLE:: conformation,effic
real,dimension(:,:,:),ALLOCATABLE:: backup_conf
real:: prov_limite,limite,distancia,norm,cutoff,th,rmsd,rmsdgl,alpha,totaal,uttolo
integer,dimension(:),ALLOCATABLE::conf_a,conf_b,vect_aux
integer,dimension(:,:),ALLOCATABLE::listado
real,dimension(:),ALLOCATABLE::aux1
real,dimension(:,:),ALLOCATABLE::matriz_a,matriz_b,frame_a,frame_b,rotato,rmsdprint
real,dimension(:,:,:),ALLOCATABLE::matrizes,bestfit
real,dimension(:,:,:),ALLOCATABLE::coordin
logical,dimension(:),ALLOCATABLE::filtro
integer,dimension(:),ALLOCATABLE::order,contador,auxiliar
logical::interruptor
TYPE(array_pointer),DIMENSION(:),POINTER::links
cutoff=0.5

!OPEN(101,FILE="../coordinates/6dotsSele09M10N4active.txt",STATUS="OLD",ACTION="READ")
OPEN(101,FILE="coordinates/6dots100MSele09.txt",STATUS="OLD",ACTION="READ")

total_frames=14280
total_dots=6

ALLOCATE(filtro(total_dots),conformation(total_dots*3),contador(total_frames))
ALLOCATE(vect_aux(8),aux1(3),rmsdprint(total_frames,total_frames))
ALLOCATE(effic(1),backup_conf(total_frames,total_dots,3))
ALLOCATE(frame_a(total_dots,3),frame_b(total_dots,3),rotato(total_dots,3))
ALLOCATE(matriz_a(total_dots,total_dots),matriz_b(total_dots,total_dots))
ALLOCATE(conf_a(total_dots),conf_b(total_dots),matrizes(total_frames,total_dots,total_dots),bestfit(total_frames,total_dots,3))

conformation=0.0d0
effic=0.0d0
uttolo=0.05

frame=0
effic=0.0d0
backup_conf=0.0d0

frame_a=0.0d0
frame_b=0.0d0
matriz_a=0.0d0
matriz_b=0.0d0
matrizes=0.0d0
rmsdprint=0.0d0

DO i=1,total_frames
   frame=frame+1
   READ(101,*) conformation(:), effic(:)	!loading process
   
   g=0
   DO k=1,total_dots
      DO h=1,3
         g=g+1
         backup_conf(i,k,h)=conformation(g)
      END DO
   END DO

   DO ii=1,total_dots
      DO jj=1,total_dots
         aux1=backup_conf(i,ii,:)-backup_conf(i,jj,:)
         norm=sqrt(dot_product(aux1,aux1))
         matrizes(i,ii,jj)=norm
      END DO
   END DO

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE(listado(24,4))
listado=0
ALLOCATE(links(total_frames))

OPEN(102,FILE="listado6.inp",STATUS="OLD",ACTION="READ")	!here we take a file with the permutation list
DO i=1,24
   READ (102,'(8I3)') listado(i,:)
END DO
CLOSE(102)

conf_b(1)=1
conf_b(2)=2

!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!


contador=0 

DO i=1,total_frames				!you compare any with any

!	bestfit(i,:,:)=backup_conf(1,:,:)!record the first structure
	matriz_a(:,:)=matrizes(i,:,:)		!distances
	DO h=i+1,total_frames

		interruptor=.false.
		matriz_b(:,:)=matrizes(h,:,:)   !distances
      		rmsdgl=1000.0d0
		DO j=1,24  
			conf_b(3:total_dots)=listado(j,:)
			DO angle=1,180		!we rotate the second structure around the axis (1,1,1) and we keep the first structure still
				th=2.0*3.1415/360.0*angle*2
				DO kkk=0,1
					DO k=1,total_dots
						rotato(k,1)=(1.0/3.0+2.0/3.0*cos(th))*(backup_conf(h,conf_b(k),1)*kkk+backup_conf(h,conf_b(k),2)*(1-kkk))+(1.0/3.0-1.0/3.0*cos(th)-sin(th)/(sqrt(3.0)))*(backup_conf(h,conf_b(k),2)*kkk+backup_conf(h,conf_b(k),1)*(1-kkk))+(1.0/3.0-1.0/3.0*cos(th)+sin(th)/(sqrt(3.0)))*backup_conf(h,conf_b(k),3)*(1)
						rotato(k,2)=(1.0/3.0-1.0/3.0*cos(th)+sin(th)/(sqrt(3.0)))*(backup_conf(h,conf_b(k),1)*kkk+backup_conf(h,conf_b(k),2)*(1-kkk))+(1.0/3.0+2.0/3.0*cos(th))*(backup_conf(h,conf_b(k),2)*kkk+backup_conf(h,conf_b(k),1)*(1-kkk))+(1.0/3.0-1.0/3.0*cos(th)-sin(th)/(sqrt(3.0)))*backup_conf(h,conf_b(k),3)*(1)
						rotato(k,3)=(1.0/3.0-1.0/3.0*cos(th)-sin(th)/(sqrt(3.0)))*(backup_conf(h,conf_b(k),1)*kkk+backup_conf(h,conf_b(k),2)*(1-kkk))+(1.0/3.0-1.0/3.0*cos(th)+sin(th)/(sqrt(3.0)))*(backup_conf(h,conf_b(k),2)*kkk+backup_conf(h,conf_b(k),1)*(1-kkk))+(1.0/3.0+2.0/3.0*cos(th))*backup_conf(h,conf_b(k),3)*(1)
					END DO	!here i have the first structure and the other rotated

					rmsd=0
					DO k=1,total_dots			
						rmsd=rmsd+(dot_product(backup_conf(i,k,:)-rotato(k,:),backup_conf(i,k,:)-rotato(k,:)))/total_dots
					END DO

					IF (rmsd<rmsdgl) THEN
						rmsdgl=rmsd
!						bestfit(h,:,:)=rotato(:,:)
					END IF
				END DO
			END DO   		
      		END DO
!		print*,rmsdgl				
!		rmsdprint(i,h)=rmsdgl
		print*,i,h,rmsdgl
	END DO
END DO

END program q_dots
