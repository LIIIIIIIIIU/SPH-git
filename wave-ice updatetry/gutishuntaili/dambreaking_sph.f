      program dambreaking_2d

	!Improved by Zheng Xing in HEU 2013_3_18
      implicit none
      include 'menu.inc'
      include 'omp_lib.h'  

      integer tsn,yesorno

	double precision s1,s2,fq,tq,cq
	open(145,file='force.dat')
	open(42,file='animate.dat')
	write(42,*)'TITLE="DISTRIBUTION"'
	write(42,*)'VARIABLES="X" "Y" "VX" "VY" "dVX" "p" "qd" "cs"'
      open(45,file='ice.dat')
	write(45,*)'TITLE="DISTRIBUTION"'
	write(45,*)'VARIABLES="X" "Y" "VX" "VY" "dVX" "p" "qd" "cs"'
      open(401,file='freesurface.dat')
      write(401,*)'TITLE="DISTRIBUTION"'
	write(401,*)'VARIABLES="X" "Y" "VX" "VY" "dVX" "p" "qd" "cs"'
      open(10,file='wave1.dat')
      open(11,file='wave2.dat')
      open(12,file='wave3.dat')
      open(13,file='eng.dat')
	open(43,file='fre_s_animate.dat')
!     open(55,file='es.dat')
	write(43,*)'TITLE="DISTRIBUTION"'
	write(43,*)'VARIABLES="X" "Y" "VX" "VY"'
!      open(66,file='ylyb.dat')
      wq=1
      px_0=0.0   !original position on x axis
	py_0=0.0   !original position on y axis
	
	xl=3.0 !width of dam
	yl=0.02  !heigth of dam
      xsl=15.9  !width of dam
	ysl=0.4   !heigth of dam
      
      
      
c	dl=0.1

	x0=16.  !Range of x direction
	y0=0.6    !Range of y direction

      px_1=x0    !max range of x direction
	py_1=y0    !max range of y direction
          
      x_n=600 !particle number on x direction
	y_n=4  !particle number on y direction 
      xs_n=3180  !particle number on x direction
	ys_n=80  !particle number on y direction 
      !     y_dn=20
      !     x_dn=800

	dx=xl/x_n  !particle scale size of x direction
	dy=yl/y_n  !particle scale size of y direction
      dx1=xsl/xs_n  !particle scale size of x direction
	dy1=ysl/ys_n  !particle scale size of y direction
      !!!!!!!!!!!!!!!!!!!!!!
         fq=0.396e+5
      cq=0.637e+5!2.0*fq
      tq=(1.-sin(jiao))/(1.+sin(jiao))*cq
      q0=cq*tq/(cq-tq)*tan(jiao)
      write(*,*)"h",q0
      nc=int(pi*bj/dx)+1
       arf=4*sin(jiao)/3.**0.5/(3.-sin(jiao))! !tan(jiao)/(9.+12*tan(jiao)**2)**0.5
c      a_rf=6*sin(jiao)/3.**0.5/(3.-sin(jiao))
      a_rf=4*sin(jzj)/3.**0.5/(3.-sin(jzj))! !sin(jiao)/3.**0.5/(3.-sin(jiao))
      sgm=6*cos(jiao)/3.**0.5/(3.-sin(jiao))! !sin(jiao)/3.**0.5/(3.-sin(jiao))
      !q0=0.9e+6/3.**0.5 !*3./(9.+12.*tan(jiao)**2)**0.5 !*(-log(gx))**(1./2.5)  !((x(i,1)-0.1)**2+x(i,2)**2)**0.5)
      JG=E/(2.*(1+miu))
      KG=E/(3.*(1-2.*miu))
      
      !!!!!!!!!!!!!!!!!!!!!!
      k=0

	do i=1,x_n
	  do j=1,y_n	
      k=k+1
	
	x(1,k)=(i-1)*dx+dx/2.+6.0!+dx/3.
      x(2,k)=(j)*dy+0.4-4.0*dx!
	xc(1,k)=x(1,k)
	xc(2,k)=x(2,k)
	
	vx(1,k)=0. !a*omigaw/sinh(k0*h)*(cosh(k0*(x(2,k)+h))*cos(k0*x(1,k)))
	vx(2,k)=0. !a*omigaw/sinh(k0*h)*(sinh(k0*(x(2,k)+h))*sin(k0*x(1,k)))
     
	rho(k)=rho0  !*(1+omiga*9.81*(yl-x(2,k))/(co*co))**(1./omiga)
	 !p(k)=E/(3.*(1.-2.*miu))*((rho(k)/rho0-1.))	
	sx(k)=1
	 mass(k)=rho0*dx*dy
        ip_f(k)=1
       jc(k)=1
        enddo
      enddo
      hsml=s_rad*dx/2.
     
      
      ntotalb=k
      do i=1,ntotalb 
          qfqd(i)=q0!*3./(9.+12.*tan(jiao)**2)**0.5  !*(-log(gx))**(1./2.5)  !((x(i,1)-0.1)**2+x(i,2)**2)**0.5)
c     &*6*cos(jiao)/3.**0.5/(3.-sin(jiao)) 
          ljsxyb(i)=0.0
c               ys(i)=qfqd(i)*sgm
      enddo

       do i=1,ntotalb 
          ip_f(i)=1
          
          if(i.ge.1.and.i.le.5.or.i.ge.596.and.i.le.600.or.mod(i,5).eq.1
     &     .or.mod(i,5).eq.0)then
          ip_f(i)=0
          n_f=n_f+1
         fre_num(n_f)=i
!          qfqd(i)=q0!*3./(9.+12.*tan(jiao)**2)**0.5  !*(-log(gx))**(1./2.5)  !((x(i,1)-0.1)**2+x(i,2)**2)**0.5)
!c     &*6*cos(jiao)/3.**0.5/(3.-sin(jiao)) 
!          ljsxyb(i)=0.0
!c               ys(i)=qfqd(i)*sgm
          endif
      enddo
      
      
      
       co=(4.*E/(2.*(1+miu))/3/rho0)**0.5 !(E*(1-miu)/(1+miu)/(1-2.*miu)/rho0)**0.5 !(E*(1-0.3)/(1+0.3)/(1.-2.*0.3)/rho0)**0.5
      write(*,*) co
      
      
      c01=5.
      k1=0
      do i=1,xs_n
	  do j=1,ys_n	
      k1=k1+1
	
	x(1,ntotalb+k1)=(i-1)*dx1+dx1/2.+0.1
      x(2,ntotalb+k1)=(j)*dy1 !+dy/2.
	
	vx(1,ntotalb+k1)=0. !a*omigaw/sinh(k0*h)*(cosh(k0*(x(2,k)+h))*cos(k0*x(1,k)))
	vx(2,ntotalb+k1)=0. !a*omigaw/sinh(k0*h)*(sinh(k0*(x(2,k)+h))*sin(k0*x(1,k)))
      
	rho(ntotalb+k1)=rho01*(1+omiga*9.81*(ysl-x(2,ntotalb+k1))/(c01*c01))
     & **(1./omiga)
		sx(ntotalb+k1)=2
	
      mass(ntotalb+k1)=rho(ntotalb+k1)*dx1*dy1
	
        enddo
      enddo
      hsmls=s_rad1*dx/2.
      !k1=0
      
      ntotals=k1
       do i=ntotalb+1,ntotalb+ntotals
      if(x(1,i).ge.(6.0).and.x(1,i).le.(9.0)!+dx/3.
     &.and.x(2,i).gt.(0.4-4.*dx))then !.and.x(2,i)& .le.0.5
      wq(i)=0
      endif
      enddo
      k1=0
       do i=ntotalb+1,ntotalb+ntotals
           if(wq(i)==1)then
       k1=k1+1
	
	x(1,ntotalb+k1)=x(1,i)
      x(2,ntotalb+k1)=x(2,i)
	
	vx(1,ntotalb+k1)=0. !a*omigaw/sinh(k0*h)*(cosh(k0*(x(2,k)+h))*cos(k0*x(1,k)))
	vx(2,ntotalb+k1)=0. !a*omigaw/sinh(k0*h)*(sinh(k0*(x(2,k)+h))*sin(k0*x(1,k)))
      
	rho(ntotalb+k1)=rho01*(1+omiga*9.81*(ysl-x(2,i))/(c01*c01))
     & **(1./omiga)
		sx(ntotalb+k1)=2
	
      mass(ntotalb+k1)=rho(i)*dx1*dy1
	endif
        enddo
       ntotals=k1
      
      
	
     
c	do i=1,x_dn-x_n
c	  do j=1,y_dn	
c      k=k+1
c	
c	x(1,k)=xl+(i-1)*dx+dx/2.
c      x(2,k)=(j-1)*dy+dy/2.
	
c	vx(1,k)=0. !a*omigaw/sinh(k0*h)*(cosh(k0*(x(2,k)+h))*cos(k0*x(1,k)))
c	vx(2,k)=0. !a*omigaw/sinh(k0*h)*(sinh(k0*(x(2,k)+h))*sin(k0*x(1,k)))
      
c	rho(k)=rho0*(1+omiga*9.81*(dl-x(2,k))/(c0*c0))**(1./omiga)
		
c	hsml(k)=dx
c      mass(k)=rho(k)*dx*dy
	
c        enddo
c      enddo

      ntotal=k+k1
      
         
c	open(4,file='int.dat')
c      do i=1,fpn
c	j=fp(i)
c      write(3,'(e14.8,'','',e14.8)')x(1,j),x(2,j)
c	enddo
c      write(4,*)'TITLE="INITITAL DISTRIBUTION"'
c	write(4,*)'VARIABLES="X" "Y" "vx" "vy"'
c	write(4,*)'ZONE I=800 F=POINT'
c      do i=1,ntotal
      
c	write(4,'(4e16.6)')x(1,i),x(2,i),vx(1,i),vx(2,i)
c	enddo
c	close(4)
      call s_grid

  10  write(*,*)'input the the total time step =tsn'
      read(*,*)tsn
	call time_print
      
	call time_elapsed(s1)
      call omp_set_num_threads(31)
	call timeintegration(tsn)
      
	call time_elapsed(s2)
      write(*,*)'expend time',(s2-s1),'seconds'
      write(*,*)'Do you want to do more calculation '
	write(*,*)'input 1: yes; input 0: no'
      read(*,*)yesorno
      if(yesorno.eq.1)go to 10

	
	
	end



      


