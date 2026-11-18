       subroutine timeintegration(tsn)
      implicit none
      include 'menu.inc'
c	tsn: time step number

    
      integer it,tsn,thist

      double precision xs,dis,v_fraction,v_c,dis1,dis2

      integer ppi,i_s_f,zp1,zp2,zp3,j0
      double precision time,dt1,time1 ,wh,wn(11),period,err,wbs,umax,
     &vmax,wa,wp,wl,wti,dis_min1,zd1,dis_min2,zd2,zd3,eng,eng1,dis_min3
     & ,dwdz,wmax,dymax,ps,x_i,y_i,z_i,x_j,y_j,z_j,pstx,psty,pstz,wij
     & ,rzi,rci,r0_h,w,dwdx,dwdy,fr1,fr2,fr3,rk1,rk2,rk3 
      
      
	thist=0
    
      do it=current+1,current+tsn

      thist=thist+1
	
	time=it*dt
      !bs=12
!-----------------------------------------------------------------
!	if(mod(it,screan_p).eq.0)then
!c      write(*,*)'************************************************'
!	write(*,*)'time step number=   ',it,'  current time=',real(time)
!c      write(*,*)'************************************************'
!      endif
!-----------------------------------------------------------------
	!wh=0.15  !wave height
 !       period=1.5
 !       wp=2.*pi/period !!tuibanyundongpinlv
 !       wn(1)=2.*pi/period/(g*ysl)**0.5    !!wave number
 !       
 !       do i=1,10
 !         wn(i+1)=wn(i)-(wn(i)-wp**2/(g*tanh(wn(i)*ysl)))/
 !    &   (1.+wp**2*ysl/g/sinh(wn(i)*ysl)**2)
 !          err=abs(wn(i+1)-wn(i))
 !         if(err.lt.1.0e-8) exit
 !         wbs=wn(i+1)
 !       enddo
 !       !write(*,*)wbs
 !       wa=wh*(2.*wbs*ysl+sinh(2.*wbs*ysl))/(4.*sinh(wbs*ysl)**2) !!tuibanfuzhi
 !       
 !     if(it==1)then
 !         wti=dt
 !         wmv=wa*wp/2.*cos(wp*wti)
	!wma=-wa*(wp**2)/2.*sin(wp*wti)
	!wmx=wa/2.*sin(wp*time)
 !         else if (it.ne.1)then
 !     if(mod(it,bs).eq.1)then
 !       dt1=dt*bs
 !          wti=(int(it/bs))*dt1+dt
 !          wmv=wa*wp/2.*cos(wp*wti)
	!wma=-wa*(wp**2)/2.*sin(wp*wti)
	!wmx=wa/2.*sin(wp*wti)
 !     endif
 !     endif
      
          !write(*,*) wa*wp/2.
        
       
          
c       wmv=wa*wp*sin(wp*t)
c	wma=wa*(wp**2)*cos(wp*t)
c	wmx=wa*(1-cos(wp*t))
      
      
      if(it.eq.1)then
	 
	 do i=1,ntotalb
      rho_0(i)=rho(i)
      vx_0(1,i)=vx(1,i)
	vx_0(2,i)=vx(2,i)
	x_0(1,i)=x(1,i)
	x_0(2,i)=x(2,i)
	 enddo
	

      call singlestep(it,time)
c k1       
	 do i=1,ntotalb
      

	rho(i)=rho_0(i)+drho(i)*dt
	vx(1,i)=vx_0(1,i)+dvx(1,i)*dt  +va(1,i)
	vx(2,i)=vx_0(2,i)+dvx(2,i)*dt  +va(2,i)
	x(1,i)=x_0(1,i)+vx(1,i)*dt
	x(2,i)=x_0(2,i)+vx(2,i)*dt
	k1vx(1,i)=vx(1,i)
      k1vx(2,i)=vx(2,i)

	k1drho(i)=drho(i)
	k1dvx(1,i)=dvx(1,i)
	k1dvx(2,i)=dvx(2,i)
	 enddo
      

      call singlestep(it,time)
c k2
       do i=1,ntotalb
	rho(i)=rho_0(i)+0.5*dt*(k1drho(i)+drho(i))
	vx(1,i)=vx_0(1,i)+dt*(0.5*k1dvx(1,i)+0.5*dvx(1,i))+va(1,i)
	vx(2,i)=vx_0(2,i)+dt*(0.5*k1dvx(2,i)+0.5*dvx(2,i))+va(2,i)
	x(1,i)=x_0(1,i)+dt*(0.5*k1vx(1,i)+0.5*vx(1,i))+0.5*dt*dt*dvx(1,i)
	x(2,i)=x_0(2,i)+dt*(0.5*k1vx(2,i)+0.5*vx(2,i))+0.5*dt*dt*dvx(2,i)
       enddo
       
!       dt1=dt*bs
!        time1=it*dt1
!        
       do i=ntotalb+1,ntotal
      rho_0(i)=rho(i)
      vx_0(1,i)=vx(1,i)
	vx_0(2,i)=vx(2,i)
	x_0(1,i)=x(1,i)
	x_0(2,i)=x(2,i)
	 enddo
	

      call singlestep1(it,time)
c k1       
	 do i=ntotalb+1,ntotal
	rho(i)=rho_0(i)+drho(i)*dt
	vx(1,i)=vx_0(1,i)+dvx(1,i)*dt  +va(1,i)
	vx(2,i)=vx_0(2,i)+dvx(2,i)*dt  +va(2,i)
 !      if(x(1,i)<=lf.and.x(1,i)>=(lf-ld)) then
	!   v(1,i)=0.05*(1-cos(pi*(x(1,i)-lf+ld)/ld))
	!   vx(1,i)=vx(1,i)-v(1,i)*vx(1,i)
	!   vx(2,i)=vx(2,i)-v(1,i)*vx(2,i)
	!end if
      
      
	x(1,i)=x_0(1,i)+vx(1,i)*dt
	x(2,i)=x_0(2,i)+vx(2,i)*dt
	k1vx(1,i)=vx(1,i)
      k1vx(2,i)=vx(2,i)

	k1drho(i)=drho(i)
	k1dvx(1,i)=dvx(1,i)
	k1dvx(2,i)=dvx(2,i)
	 enddo
      

      call singlestep1(it,time)
c k2
       do i=ntotalb+1,ntotal
	rho(i)=rho_0(i)+0.5*dt*(k1drho(i)+drho(i))
	vx(1,i)=vx_0(1,i)+0.5*dt*(k1dvx(1,i)+dvx(1,i))+va(1,i)
	vx(2,i)=vx_0(2,i)+0.5*dt*(k1dvx(2,i)+dvx(2,i))+va(2,i)
 !      if(x(1,i)<=lf.and.x(1,i)>=(lf-ld)) then
	!   v(i)=0.1*(1-cos(pi*(x(1,i)-lf+ld)/ld))
	!   vx(1,i)=vx(1,i)-v(i)*vx(1,i)
	!   vx(2,i)=vx(2,i)-v(i)*vx(2,i)
	!end if
	x(1,i)=x_0(1,i)+0.5*dt*(k1vx(1,i)+vx(1,i))+0.5*dt*dt*dvx(1,i)
	x(2,i)=x_0(2,i)+0.5*dt*(k1vx(2,i)+vx(2,i))+0.5*dt*dt*dvx(2,i)
       enddo

        do i=ntotalb+1,ntotal
        if(x(1,i)<=lf.and.x(1,i)>=(lf-ld)) then
	   !v(1,i)=0.1*(1-cos(pi*(x(1,i)-lf+ld)/ld))
      v(i)=1.0*((exp(((x(1,i)-lf+ld)/ld)**np)-1)/(exp(1.)-1))
	   vx(1,i)=vx(1,i)-v(i)*vx(1,i)
	   vx(2,i)=vx(2,i)-v(i)*vx(2,i)
	end if
        enddo
        
      elseif(it.ne.1)then
         !$omp parallel private(i) 
      
      
      !$omp do schedule(dynamic,50)   
          
        do i=1,ntotalb
        rho_0(i)=rho(i)
        vx_0(1,i)=vx(1,i)
	  vx_0(2,i)=vx(2,i)
	  x_0(1,i)=x(1,i)
	  x_0(2,i)=x(2,i)
c k1	 
        

	  rho(i)=rho_0(i)+drho(i)*dt
	  vx(1,i)=vx_0(1,i)+dvx(1,i)*dt +va(1,i)
	  vx(2,i)=vx_0(2,i)+dvx(2,i)*dt +va(2,i)
	  x(1,i)=x_0(1,i)+vx(1,i)*dt
	  x(2,i)=x_0(2,i)+vx(2,i)*dt
	  
	  k1drho(i)=drho(i)
	  k1dvx(1,i)=dvx(1,i)
	  k1dvx(2,i)=dvx(2,i)
	  k1vx(1,i)=vx(1,i)
	  k1vx(2,i)=vx(2,i)
	 enddo
	 !$omp end do 
    ! 
        !$omp end parallel

      call singlestep(it,time)
c k2
      !$omp parallel private(i) 
      
      
      !$omp do schedule(dynamic,50)   
        do i=1,ntotalb
	rho(i)=rho_0(i)+0.5*dt*(k1drho(i)+drho(i))
	vx(1,i)=vx_0(1,i)+dt*(0.5*k1dvx(1,i)+0.5*dvx(1,i))+va(1,i)
	vx(2,i)=vx_0(2,i)+dt*(0.5*k1dvx(2,i)+0.5*dvx(2,i))+va(2,i)
	x(1,i)=x_0(1,i)+dt*(0.5*k1vx(1,i)+0.5*vx(1,i))+0.5*dt*dt*dvx(1,i)
	x(2,i)=x_0(2,i)+dt*(0.5*k1vx(2,i)+0.5*vx(2,i))+0.5*dt*dt*dvx(2,i)
        enddo
       !$omp end do 
    ! 
        !$omp end parallel
       
 !        if(mod(it,28).eq.1)then
 !       do i=ntotalb+1,ntotal
 !
 !      if(x(1,i)<=lf.and.x(1,i)>=(lf-ld)) then
	!   v(1,i)=0.1*(1-cos(pi*(x(1,i)-lf+ld)/ld))
	!   vx(1,i)=vx(1,i)-v(1,i)*vx(1,i)
	!   vx(2,i)=vx(2,i)-v(1,i)*vx(2,i)
	!end if
	!
 !      enddo
 !       endif
        
       
       if(mod(it,bs).eq.1)then
           
  !     !$omp parallel private(i) 
  !      
  !    !$omp do  
  !        do i=ntotalb+1,ntotal
  !          !Fluid  particle near solid wall should add damping zone
	 ! if((x(1,i)-wmx).le.0.7*dx.and.vx(1,i).lt.0.0)then
  !    vx(1,i)=0.2*vx(1,i)
  !      endif
  !   
	 !if((x(2,i)-py_0).le.0.7*dx.and.vx(2,i).lt.0.0)then
  !    vx(2,i)=0.2*vx(2,i)
  !     endif
	 !if((py_1-x(2,i)).le.0.7*dx.and.vx(2,i).gt.0.0)then
  !    vx(2,i)=0.2*vx(2,i)
  !     endif
  !    if((px_1-x(1,i)).le.0.7*dx.and.vx(1,i).gt.0.0)then
  !    vx(1,i)=0.2*vx(1,i)
  !    endif
  !    
	 !    enddo 
  !        !$omp end do 
  !      !$omp end parallel  
 !          !$omp parallel private(i) 
 !     
 !     
 !     !$omp do
 !         do i=ntotalb+1,ntotal
 !
 !      if(x(1,i)<=lf.and.x(1,i)>=(lf-ld)) then
	!   v(1,i)=0.1*(1-cos(pi*(x(1,i)-lf+ld)/ld))
	!   vx(1,i)=vx(1,i)-v(1,i)*vx(1,i)
	!   vx(2,i)=vx(2,i)-v(1,i)*vx(2,i)
	!end if
	!
 !         enddo
 !        !$omp end do 
 !   ! 
 !       !$omp end parallel  
          
          
           dt1=dt*bs
           time1=(int(it/bs))*dt1+dt
           
         !$omp parallel private(i) 
      
      
      !$omp do schedule(dynamic,50)  
        do i=ntotalb+1,ntotal
        rho_0(i)=rho(i)
        vx_0(1,i)=vx(1,i)
	  vx_0(2,i)=vx(2,i)
	  x_0(1,i)=x(1,i)
	  x_0(2,i)=x(2,i)
c k1	 
        

	  rho(i)=rho_0(i)+drho(i)*dt1
	  vx(1,i)=vx_0(1,i)+dvx(1,i)*dt1 +va(1,i)
	  vx(2,i)=vx_0(2,i)+dvx(2,i)*dt1 +va(2,i)
 !        if(x(1,i)<=lf.and.x(1,i)>=(lf-ld)) then
	!   v(1,i)=0.05*(1-cos(pi*(x(1,i)-lf+ld)/ld))
	!   vx(1,i)=vx(1,i)-v(1,i)*vx(1,i)
	!   vx(2,i)=vx(2,i)-v(1,i)*vx(2,i)
	!end if
	  x(1,i)=x_0(1,i)+vx(1,i)*dt1
	  x(2,i)=x_0(2,i)+vx(2,i)*dt1
	  
	  k1drho(i)=drho(i)
	  k1dvx(1,i)=dvx(1,i)
	  k1dvx(2,i)=dvx(2,i)
	  k1vx(1,i)=vx(1,i)
	  k1vx(2,i)=vx(2,i)
	 enddo
	 !$omp end do 
    ! 
        !$omp end parallel  
      !call singlestep(it,time)
      call singlestep1(int(it/bs),time1)
c k2
        !$omp parallel private(i) 
      
      
      !$omp do schedule(dynamic,50)  
       do i=ntotalb+1,ntotal
	rho(i)=rho_0(i)+0.5*dt1*(k1drho(i)+drho(i))
	vx(1,i)=vx_0(1,i)+0.5*dt1*(k1dvx(1,i)+dvx(1,i))+va(1,i)
	vx(2,i)=vx_0(2,i)+0.5*dt1*(k1dvx(2,i)+dvx(2,i))+va(2,i)
 !      if(x(1,i)<=lf.and.x(1,i)>=(lf-ld)) then
	!   v(1,i)=0.05*(1-cos(pi*(x(1,i)-lf+ld)/ld))
	!   vx(1,i)=vx(1,i)-v(1,i)*vx(1,i)
	!   vx(2,i)=vx(2,i)-v(1,i)*vx(2,i)
	!end if
	x(1,i)=x_0(1,i)+0.5*dt1*(k1vx(1,i)+vx(1,i))+0.5*dt1*dt1*dvx(1,i)
	x(2,i)=x_0(2,i)+0.5*dt1*(k1vx(2,i)+vx(2,i))+0.5*dt1*dt1*dvx(2,i)
       enddo
       !$omp end do 
    ! 
        !$omp end parallel 
         do i=ntotalb+1,ntotal
        if(x(1,i)<=lf.and.x(1,i)>=(lf-ld)) then
	   !v(1,i)=0.1*(1-cos(pi*(x(1,i)-lf+ld)/ld))
      v(i)=1.0*((exp(((x(1,i)-lf+ld)/ld)**np)-1)/(exp(1.)-1))
	   vx(1,i)=vx(1,i)-v(i)*vx(1,i)
	   vx(2,i)=vx(2,i)-v(i)*vx(2,i)
	end if
         enddo
       endif
       
       
       endif
       
       
       
       
       
       if(mod(it,bs).eq.1)then
       dt1=dt*bs   
       
        umax=0.0
	vmax=0.0
      !!wmax=0.0
       dymax=0.0
       ps=0.0
       
      !! !$omp parallel private(i)   
      !!
      !!!$omp do schedule(dynamic,50)  
        do i=ntotalb+1,ntotal
     !!!!  if(x(1,i).lt.px_1.and.x(2,i).lt.py_1.and.x(1,i).gt.wmx.and.x(2,i)
     !!!!& .gt.py_0)then         
	if(umax.le.abs(vx(1,i))) then
	umax=abs(vx(1,i))
	end if
      if(vmax.le.abs(vx(2,i))) then
	vmax=abs(vx(2,i))
      end if
 !     if(wmax.le.abs(vx(3,i))) then
	!wmax=abs(vx(3,i))
 !     end if
 !!!!     if(dymax.le.dyl(i)) then
	!!!!dymax=dyl(i)
 !!!!     end if
 !!!!     endif
 !!!!     ps=ps+sqrt(vx(1,i)**2+vx(2,i)**2)
      end do
 !!!!     if(mod(it,1000).eq.0)then
 !!!!     write(*,*)  umax,vmax
 !!!!     endif
 !!!!     
 !!!!     write(555,'(3e14.6)') time,dymax/(ps/ntotal)**2,!(umax**2+vmax**2),!
 !!!!    & ps/ntotal/sqrt(9.81*sqrt(yl*xl))!ps/ntotalsqrt(umax**2+vmax**2)
 !!!!     
      !!!mass0=rho0*dx*dx*dx
      
      
 !!!!       umax=0.07
	!!!!vmax=0.14
 !!!!     
      
      x_i=0.0
       y_i=0.0
        !!z_i=0.0
        x_j=0.0
       y_j=0.0
       !!z_j=0.0
       r0_h=hsmls  
      call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,r0_h) 
      wij=w 
      
      
      !$omp parallel  private(i,j,j0,x_i,y_i)
     &private(x_j,y_j,r0_h,fr1,fr2,rci)
     & private(rk1,rk2,w,dwdx,dwdy,dis,rzi,pstx,psty)  
      !$omp do schedule(dynamic) 
 	 do i=ntotalb+1,ntotal
      x_i=x(1,i)
	 y_i=x(2,i)
       !!z_i=x(3,i)
      r0_h=hsmls    ! 0.5*s_rad*dx  
 !!     dudx=0
	!!dvdy=0
 !!     dvdz=0.0
      rzi=0.0
	fr1=0.0
	fr2=0.0
      !!fr3=0.0
      pstx=0.0
      psty=0.0
      !!pstz=0.0
 !     do j=1,rn(i)
	! j0=rpn(j,i)
	! if(j0.ne.i)then     
 !      x_j=x(1,j0)
	!  y_j=x(2,j0)
 !     dis=sqrt( (x_i-x_j)**2+(y_i-y_j)**2 )
 !     call kernel(x_i,y_i,x_j,y_j,w,dwdx,dwdy,r0_h)
 !     dudx=dudx+vx(1,j0)*mass0/rho0*dwdx !
 !     dvdy=dvdy+vx(2,j0)*mass0/rho0*dwdy ! 
	!end if   
	!end do

!!      do j=1,rn(i)
!!	 j0=rpn(j,i)
!!c	 if(ip_f(j0).ne.1)then 
!!      x_j=x(1,j0)
!!	  y_j=x(2,j0)
!!        z_j=x(3,j0)
!!      dis=sqrt( (x_i-x_j)**2+(y_i-y_j)**2+(z_i-z_j)**2)  
!!      rzi=rzi+dis
!!c	endif
!!      end do
      
       do j=1,rn1(i)
      j0=rpn1(j,i)
	 if(j0.ne.i)then
      x_j=x(1,j0)
	  y_j=x(2,j0)
        !!z_j=x(3,j0)
         dis=sqrt( (x_i-x_j)**2+(y_i-y_j)**2)  
      rzi=rzi+dis
      call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,r0_h)  
        
       pstx=pstx+(1.+0.2*(w/wij)**4)*dx1*dx1/2.*dwdx 
       psty=psty+(1.+0.2*(w/wij)**4)*dx1*dx1/2.*dwdy 
       !!pstz=pstz+(1.+0.2*(w/wij)**4)*dx1*dx1*dx1/2.*dwdz 
       
      rk1=x(1,i)-x(1,j0)
      rk2=x(2,i)-x(2,j0)
      !!rk3=x(3,i)-x(3,j0)
	dis=sqrt( (x_i-x_j)**2+(y_i-y_j)**2)
      if(dis.ne.0.0) then 
      fr1=fr1+rk1/dis**3
      fr2=fr2+rk2/dis**3
      !!fr3=fr3+rk3/dis**3
	end if
	end if
      end do
      if(rn1(i).ne.0)then
	rci=rzi/rn1(i)
      endif
       x(1,i)=x(1,i)+0.00005*fr1*rci**2*sqrt(umax**2+vmax**2)
     &   *dt1
	 x(2,i)=x(2,i)+0.00005*fr2*rci**2*sqrt(umax**2+vmax**2)
     &  *dt1
     !!  x(3,i)=x(3,i)+0.00001*fr3*rci**2*sqrt(umax**2+vmax**2+wmax**2)
     !!&  *dt1
  !!     x(1,i)=x(1,i)-sqrt(umax**2+vmax**2+wmax**2)*dt*4.*r0_h*pstx
	 !!x(2,i)=x(2,i)-sqrt(umax**2+vmax**2+wmax**2)*dt*4.*r0_h*psty
  !!     x(3,i)=x(3,i)-sqrt(umax**2+vmax**2+wmax**2)*dt*4.*r0_h*pstz


  !!     x(1,i)=x(1,i)+0.01*fr1*rci**2*umax*dt
	 !!x(2,i)=x(2,i)+0.01*fr2*rci**2*vmax*dt
  !!     x(3,i)=x(3,i)+0.01*fr3*rci**2*wmax*dt
      
C        x(1,i)=x(1,i)+0.05*fr1*rc(i)**2*umax*dt  
C	 x(2,i)=x(2,i)+0.05*fr2*rc(i)**2*vmax*dt
c       vx(1,i)=vx(1,i)+dudx*0.1*fr1*rc(i)**2*umax*dt   
c        vx(2,i)=vx(2,i)+dvdy*0.1*fr2*rc(i)**2*vmax*dt
c        vx(1,i)=vx(1,i)+vg(1,i)   
c        vx(2,i)=vx(2,i)+vg(2,i)
	     
      end do
  
      !$omp end do
      !$omp end parallel            
          
      endif
c      write(*,*)it,rho(6400),p(6400),vx(6400,2),dvx(6400,2)
 
	! calculating the point pressure
!----------------------------------------------------      
c	 ppx1=cx3
c	ppy1=cy3+0.15

c	pp1=0.0
c      pp2=0.0
c      zyb=0.0
c      zyl=0.0
c	ppi=0
c	do i=1,ntotal

c	dis=sqrt((x(1,i)-ppx1)**2+(x(2,i)-ppy1)**2)

c      if(dis.le.(2.1*dx))then
c	ppi=ppi+1
c	pp1=pp1+jsd(i)*mass
c      pp2=pp2-p(i)+sts(2,i)
c	endif
c	enddo 
c      if(ppi.ge.1)then
c	pp1=pp1/ppi
c      pp2=pp2/ppi
c	else
c	pp1=0.0
c      pp2=0.0
c      endif
c      pp3=-p(364)+sts(2,364)+sts(3,364)
c      pp4=-p(364)+sts(2,364)
c          li=0.0
c      do i=1,ntotal
c          zyb=zyb+sts(2,i)
c          zyl=zyl+mass*9.81*x(2,i)
c      enddo
       if(mod(it,100).eq.0)then
      write(145,'(7e14.6)')it*dt,x(1,600),x(2,600),x(1,604),
     &x(2,604),ybl(1,161)+ybl(2,161),ybl(1,164)+ybl(2,164)!,jsd(298)*mass(298)
         endif 
     
       

!-------------------------------------------------
         
      if(it.eq.1.or.mod(it,100).eq.1)then 
          dis_min1=1.3*dx
	zp1=0
      zd1=0.0
      
       dis_min2=1.3*dx
	zp2=0
      zd2=0.0
      dis_min3=1.3*dx
      zp3=0
      zd3=0.0
      do i=1,ntotal
         if(x(1,i).gt.8.5.and.x(1,i).lt.23.5)then
	dis=abs(x(1,i)-9.0)
      if(dis.le.1.3*dx)then
	 if(zd1.le.x(2,i))then
	zd1=x(2,i)
      zp1=i      
      end if
      endif
      
      dis1=abs(x(1,i)-14.0)
      if(dis1.le.1.3*dx)then
	 if(zd2.le.x(2,i))then
	zd2=x(2,i)
      zp2=i      
      end if
      endif
      
      dis2=abs(x(1,i)-22.)
      if(dis2.le.dis_min3)then
          if(zd3.le.x(2,i))then
           !dis_min3=dis2
	 
	zd3=x(2,i)
      zp3=i      
      end if
      endif
      
      endif
	enddo
	if(zp1.ne.0)then
	write(10,'(2e16.6)')time,x(2,zp1)
      endif 
      
      if(zp2.ne.0)then
	write(11,'(2e16.6)')time,x(2,zp2)
      endif 
      
      if(zp3.ne.0)then
	write(12,'(2e16.6)')time,x(2,zp3)
      endif 
      
      eng=0.0
       do i=1,ntotalb
      eng=eng+rho0*dx*dx*g*x(2,i)
      eng1=eng1+0.5*rho0*dx*dx*(vx(1,i)**2+vx(2,i)**2)
      enddo
      
      write(13,'(4e16.6)')time,eng,eng1,eng+eng1
      
        endif  
 !        
 !     if(it.eq.1.or.mod(it,50).eq.0)then
	!!rho_t=0.0
	!i_s_f=0
	!!p_n=0.0
	!do i=1,ntotal
	!  !p_n=p_n+mass*9.81*x(2,i)
 !  !     rho_t=rho_t+rho(i)
	!  v_c=0
	!  do j=1,rn1(i)
 !       v_c=v_c+dx*dy
	!  enddo
 !       !p_n=p_n+0.5*mass*(vx(1,i)**2+vx(2,i)**2)
	!v_fraction=v_c/(pi*(s_rad1*dx)**2)
 !
	!if(v_fraction.le.0.73)then
 !
 !
	!i_s_f=i_s_f+1
	!free_surface(i_s_f)=i
	!endif
 !     enddo
 !     
 !    
 !     
 !     write(43,*)'ZONE I=',i_s_f,' F=POINT'
	!do i=1,i_s_f
	!j=free_surface(i)
 !     write(43,'(4e16.6)')x(1,j),x(2,j),vx(1,j),vx(2,j)
	!enddo
 !     endif
 !     
      

       if(it.eq.1.or.it.eq.3188571.or.it.eq.4200000.or.it.eq.5580000
     &      .or.mod(it,100000).eq.0)then
	!rho_t=0.0
	!i_s_f=0
	!p_n=0.0
	!do i=1,ntotal
	!  p_n=p_n+mass*9.81*x(2,i)
 !       rho_t=rho_t+rho(i)
	!  v_c=0
	!  do j=1,rn(i)
 !       v_c=v_c+dx*dy
	!  enddo
 !       p_n=p_n+0.5*mass*(vx(1,i)**2+vx(2,i)**2)
	!v_fraction=v_c/(pi*(s_rad*dx)**2)
 !
	!if(v_fraction.le.0.73)then
 !
 !
 !
 !
 !
 !
 !
 !
 !
 !
	!i_s_f=i_s_f+1
	!free_surface(i_s_f)=i
	!endif
	!enddo
c      rho_t=rho_t/ntotal 
c	write(43,*)'ZONE I=',i_s_f,' F=POINT'
c	do i=1,i_s_f
c	j=free_surface(i)
c      write(43,'(4e16.6)')x(1,j),x(2,j)+0.5*dx,vx(1,j),vx(2,j)
c	enddo
	write(42,*)'ZONE I=',ntotals,' F=POINT'
	do i=ntotalb+1,ntotal
      write(42,'(7e16.6,i2)')x(1,i),x(2,i),vx(1,i),rho(i),
     &p(i),tsts(1,i),ljsxyb(i),fs(i)
	enddo
       write(45,*)'ZONE I=',ntotalb,' F=POINT'
	do i=1,ntotalb
      write(45,'(6e16.6,2i2)')x(1,i),x(2,i),rho(i),
     &dvx(2,i),p(i),ljsxyb(i),sx(i),fs(i)
	enddo 
        write(401,*)'ZONE I=',n_f,' F=POINT'
	do i=1,n_f
          j=fre_num(i)
      write(401,'(7e16.6,i2)')x(1,j),x(2,j),vx(1,j),vx(2,j),
     &tsts(1,j),p(j),ljsxyb(j),jc(j)
	enddo

c	write(71,*)'ZONE I=',y_n+y_dn,' J=',x_n,' F=POINT'
c	do i=1,ntotal1
c      write(71,'(4e16.6)')x(1,i),x(2,i),p(i)
c	enddo
c	do i=ntotal1+1,ntotal
c      write(71,'(4e16.6)')x(1,i),x(2,i),p(i)
c	enddo
	endif

c	write(40,'(i5,e14.6)')it,dvx(2,200)
c	write(45,'(3e14.6)')time,rho_t,p_n
c      if(it.eq.1.or.mod(it,50).eq.0)then
c      write(66,'(3e16.6)')it*dt,x(2,6380),x(2,6381)
c      endif
	enddo
      
      current=current+thist
	
	
	end



