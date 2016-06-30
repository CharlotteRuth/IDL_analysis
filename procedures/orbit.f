C program to put two n-body systems in snapshot format on a given orbit
C 5/1993 version by F.Governato (governato@astmib.astro.it)

	parameter(idim=34000)
        real xp(3),xdotp(3),xdot(3,idim),m1,m2,
     &	x(3,idim)
	real xc(3),xcdot(3),a(8),body(idim),eps(idim)
	character*1 orbit
        character *70, ifile

	pi2=1.57079632

111	format(70a)
        write (*,*) 'Enter output file name'
        read (*,111) ifile

	write(6,*) 'input m1,m2,p,r,ecc '
	read(5,*)m1,m2,p,r,ecc
c m1: mass of first system; m2: mass of 2nd system
c p: pericenter distance; r: point at which collision starts;
c ecc :eccentricity of orbit


	write(6,*)' input type of orbit: e,p,h '
	read(5,1)orbit
1       format(a1)
	if(orbit.eq.'e'.or.orbit.eq.'E')then
	write(25,*)' elliptic orbit'

c elliptic orbit

	c=1.-ecc
	semi=p/c
	safe=semi*(1.+ecc)
44      if(r.gt.safe)then
	write(6,*)'you are starting from beyond apocentre!
     >  try again'
	write(6,*)' input r'
	read(5,*)r
	goto 44
	endif


	c1=semi*(1-ecc**2)
	v=sqrt((m1+m2)*(2./r-1./semi))
	f=acos((1./ecc)*(c1/r-1))
	phi=pi2+acos(sqrt(semi*c1/(r*(2.*semi-r))))
	theta=phi-f
	xdotp(1)=v*cos(theta)
	xdotp(2)=v*sin(theta)
	xdotp(3)=0.

c  position of perturber

	if(f.lt.pi2)then
	xp(1)=-r*cos(f)
	else
	xp(1)=r*cos(f)
	endif
	xp(2)=-r*sin(f)
	xp(3)=0.
	goto 999

	else

c parabolic orbit

	if(orbit.eq.'p'.or.orbit.eq.'P')then
	write(25,*)' parabolic orbit'
	v=sqrt(2.*(m1+m2)/r)
	p=p*2
	f=acos(p/r-1)
	phi=pi2+acos(p/(2.*r))
	theta=phi-f
	xdotp(1)=v*cos(theta)
	xdotp(2)=v*sin(theta)
	xdotp(3)=0.

c xy position

	xp(1)=r*cos(f)
	xp(2)=-r*sin(f)
	xp(3)=0.
	goto 999

	else

c hyperbolic orbit

	write(25,*) ' hyperbolic orbit '

	d=p/(ecc-1)
	semi=d*(ecc**2-1)
	v=sqrt((m1+m2)*(2./r+1./d))
	f=acos((semi/r-1)/ecc)
	phi=pi2+acos(sqrt((d**2*(ecc**2-1))/(r*(2.*d+r))))
	theta=phi-f
	xdotp(1)=v*cos(theta)
	xdotp(2)=v*sin(theta)
	xdotp(3)=0.


c xy position

	xp(1)=r*cos(f)
	xp(2)=-r*sin(f)
	xp(3)=0.
	endif
	endif

999     continue
	

c unit 33 contains data for first system, unit 34 for second
	open(33,file='sys1.dat',status='old',form='formatted')
        open(34,file='sys2.dat',status='old',form='formatted')


	read(33,*)n
	read(33,*)ajunk1
	read(33,*)ajunk2
	do j=1,n
	read(33,*)body(j)
	enddo
	do j=1,n
	read(33,*)x(1,j),x(2,j),x(3,j)
	enddo
	do j=1,n
	read(33,*)xdot(1,j),xdot(2,j),xdot(3,j)
	enddo
	do j=1,n
	read(33,*)eps(j)
	enddo

	read(34,*)n1
	read(34,*)ajunk1
	read(34,*)ajunk2
	do j=n+1,n+n1
	read(34,*)body(j)
	enddo
	do j=n+1,n+n1
	read(34,*)x(1,j),x(2,j),x(3,j)
	enddo
	do j=n+1,n+n1
	read(34,*)xdot(1,j),xdot(2,j),xdot(3,j)
	enddo
	do j=n+1,n+n1
	read(34,*)eps(j)
	enddo


	write(25,22)m1,m2,p,r,ecc
22	format(/' m1= ',f8.3,' m2= ',f8.3,' p= ',f8.3,' r= ',f8.3,
     >  ' ecc= ',f8.3)
	write(25,11) (xp(k),k=1,3),(xdotp(k),k=1,3)
11      format(1x,//' coordinates of perturber',6f8.3//)


c   for 2nd system
	nt=n+n1
	do 34  k=1,3
	do 34 i=n+1,nt
	x(k,i)=x(k,i)+xp(k)
34      xdot(k,i)=xdot(k,i)+xdotp(k)

c for whole system

	do 30 k = 1,3
	xc(k)=0.
30      xcdot(k)=0.



        bd=0.
	do  i=1,nt
	bd=bd+body(i)
	enddo	
	write(25,23)bd
23	format(' mass of sys1 + sy2 (from sum of body(i)) ',f7.1)

	do 369 i = 1,nt
	do 369 k = 1,3
	xc(k)=xc(k)+body(i)*x(k,i)
369     xcdot(k)=xcdot(k)+body(i)*xdot(k,i)



	do 365 k = 1, 3
	xc(k)=xc(k)/bd
365     xcdot(k)=xcdot(k)/bd

	do 909 k=1,3
	do 909 i=1,nt
	x(k,i)=x(k,i)-xc(k)
909     xdot(k,i)=xdot(k,i)-xcdot(k)



	write(25,112)xc,xcdot
112     format(1x/,' coordinates shift for cm of whole system'
     >  ,6f10.5/)

        open (26,file=ifile,status='new',form='formatted')
	ndim=3
	tnow=0
	write(26,*)nt
	write(26,*)ndim
	writE(26,*)tnow          
	do j=1,nt
	write(26,*)body(j)
	enddo
	do j=1,nt
	write(26,*)x(1,j),x(2,j),x(3,j)
	enddo
	do j=1,nt
	write(26,*)xdot(1,j),xdot(2,j),xdot(3,j)
	enddo
	do j=1,nt
	write(26,*)eps(j)    
	enddo
   


36      continue
	end
