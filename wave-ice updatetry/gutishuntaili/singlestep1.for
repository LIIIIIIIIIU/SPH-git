      subroutine singlestep1(t,ti)
      implicit none
	include 'menu.inc'
  
      double precision ti
c      real a0,b0,c0,ss,ccs,et(maxn),ph(maxn),tao,dis
	integer t,fpn,nj
      integer coor_1,coor_2,coor_3,coor_4,coor_x_1,
     &coor_x_2,coor_x_3,coor_x_4,ip0,ip1,js,js1
      double precision cx1,cy1,cx2,cy2,cx3,cy3,cx4,cy4,jiao0,juli
     &	,dis
     & ,pp,tsts1,tsts2,tsts3,tsts4,x_i,y_i
     &,x_j,y_j,xx_j,yy_j,n0,dwdx,dwdy,hsml0,w
     &,xi(dim),xj(dim) ,wh,wn(11),period,err,wbs,wti,
     & wa,wp,wl,time,dt1,dis_i,sin_sita,cos_sita

      nvirt1=0
	nvirt2=0
	nvirt3=0
      cx3=0.35  !1.0+1.6625 !1.+66.5*dx
	cy3=yl+bj+ysl+vs*ti+0.5*dx
!      if(t==1)then
!          time=dt
!          else if (t.ne.1)then
!      if(mod(t,10).eq.1)then
!        dt1=dt*10
!           time=(int(t/10))*dt1
!      endif
!      endif
 !     bs=12
         wh=0.025 !wave height
        period=1.2
        wp=2.*pi/period !!tuibanyundongpinlv
        wn(1)=2.*pi/period/(g*ysl)**0.5    !!wave number
        
        do i=1,10
          wn(i+1)=wn(i)-(wn(i)-wp**2/(g*tanh(wn(i)*ysl)))/
     &   (1.+wp**2*ysl/g/sinh(wn(i)*ysl)**2)
           err=abs(wn(i+1)-wn(i))
          if(err.lt.1.0e-8) exit
          wbs=wn(i+1)
        enddo
        !write(*,*)wbs,i
        !wh=0.3/wbs
        
        wa=wh*(2.*wbs*ysl+sinh(2.*wbs*ysl))/(4.*sinh(wbs*ysl)**2) !!tuibanfuzhi
        
 !     if(t==1)then
 !         wti=dt
 !         wmv=wa*wp/2.*cos(wp*wti)
	!wma=-wa*(wp**2)/2.*sin(wp*wti)
	!wmx=wa/2.*sin(wp*time)
 !         else if (t.ne.1)then
 !     if(mod(t,bs).eq.1)then
 !       dt1=dt*bs
           wti=ti!(int(t/bs))*dt1+dt
           wmv=wa*wp/2.*cos(wp*wti)
	wma=-wa*(wp**2)/2.*sin(wp*wti)
	wmx=wa/2.*sin(wp*wti)+0.1
      !call wavemaker(ti)
    
      
 !     U=aa*bb*sin(bb*ti)
	!dU=aa*(bb**2)*cos(bb*ti)
	!S=aa*(1-cos(bb*ti))
 !     
!             do i=1,ntotalb
!c       if(i_p(i).eq.0)then
!c	p(i)=0.0
!c      rho(i)=rho0
!c	 else
!      !p(i)=-(sts(1,i)+sts(2,i))/2.
!      p(i)=E/(3*(1-2.*miu))*((rho(i)/rho0-1.))  !-1./3.*(tsts(i,2)+tsts(i,2))c0*c0*(rho(i)-rho0) !(1-0.5*tao*et(i))*ph(i)+tao*rho(i)*E  !c0*c0*rho0/omiga*((rho(i)/rho0)**omiga-1)
!c	 endif
!	!mass(i)=dx*dy*rho(i)
!      !write(*,*) p(6400),mass(6400)
!      enddo
!      
      
      
       do i=ntotalb+1,ntotal
c       if(i_p(i).eq.0)then
c	p(i)=0.0
c      rho(i)=rho0
c	 else
      p(i)=c01*c01*rho01/omiga*((rho(i)/rho01)**omiga-1)
c	 endif
c	mass(i)=dx*dy*rho(i)
       enddo
        
      do i=ntotalb+1,ntotal
       if(x(2,i).gt.dy1/2..and.(x(2,i)-dy1/2.).le.s_rad1*dx1)then
      nvirt1=nvirt1+1
      x(1,ntotal+nvirt1)=x(1,i) 
	x(2,ntotal+nvirt1)=2.*dy1/2.-x(2,i)  
	!hsml(ntotal+nvirt1)=hsml(i)
	rho(ntotal+nvirt1)=rho(i)
	mass(ntotal+nvirt1)=dx1*dy1*rho(ntotal+nvirt1)
	p(ntotal+nvirt1)=p(i)+rho01*g*(x(2,i)-x(2,ntotal+nvirt1))
      !sts(1,ntotal+nvirt1)=sts(1,i)
      !sts(2,ntotal+nvirt1)=sts(2,i)
      !sts(3,ntotal+nvirt1)=sts(3,i)
      !sts(4,ntotal+nvirt1)=sts(4,i)
	vx(1,ntotal+nvirt1)=vx(1,i) 
	vx(2,ntotal+nvirt1)=-vx(2,i) 
      sx(ntotal+nvirt1)=3
         elseif((y0-x(2,i)).le.s_rad1*dx1)then
      nvirt1=nvirt1+1
      x(1,ntotal+nvirt1)=x(1,i) 
	x(2,ntotal+nvirt1)=2.*y0-x(2,i)  
	!hsml(ntotal+nvirt1)=hsml(i)
	rho(ntotal+nvirt1)=rho(i)
	mass(ntotal+nvirt1)=dx1*dy1*rho(ntotal+nvirt1)
	p(ntotal+nvirt1)=p(i)
      ! sts(1,ntotal+nvirt1)=sts(1,i)
      !sts(2,ntotal+nvirt1)=sts(2,i)
      !sts(3,ntotal+nvirt1)=sts(3,i)
      !sts(4,ntotal+nvirt1)=sts(4,i)
	vx(1,ntotal+nvirt1)=vx(1,i) 
	vx(2,ntotal+nvirt1)=-vx(2,i) 
	sx(ntotal+nvirt1)=3
	   
	   endif
 !      if(x(1,i)<=lf.and.x(1,i)>=(lf-ld)) then
	!   v(i)=0.2*(1-cos(pi*(x(1,i)-lf+ld)/ld))
	!   vx(1,i)=vx(1,i)-v(i)*vx(1,i)
	!   vx(2,i)=vx(2,i)-v(i)*vx(2,i)
	!end if

	 enddo

c	yd=ao-ao*cos(omigao*ti)
c      vd=ao*omigao*sin(omigao*ti)

	do i=ntotalb+1, ntotal+nvirt1

      if(x(1,i).gt.wmx.and.(x(1,i)-wmx).le.s_rad1*dx1)then

      nvirt2=nvirt2+1
      x(1,ntotal+nvirt1+nvirt2)=2.*wmx-x(1,i) 
	x(2,ntotal+nvirt1+nvirt2)=x(2,i)  
	!hsml(ntotal+nvirt1+nvirt2)=hsml(i)
	rho(ntotal+nvirt1+nvirt2)=rho(i)
	mass(ntotal+nvirt1+nvirt2)=dx1*dy1*rho(ntotal+nvirt1+nvirt2)
	p(ntotal+nvirt1+nvirt2)=p(i)-rho01*(x(1,ntotal+nvirt1+nvirt2)-x(1,i))
     &*wma
      ! sts(1,ntotal+nvirt1+nvirt2)=sts(1,i)
      !sts(2,ntotal+nvirt1+nvirt2)=sts(2,i)
      !sts(3,ntotal+nvirt1+nvirt2)=sts(3,i)
      !sts(4,ntotal+nvirt1+nvirt2)=sts(4,i)
	vx(1,ntotal+nvirt1+nvirt2)=2.*wmv-vx(1,i) 
	vx(2,ntotal+nvirt1+nvirt2)=vx(2,i) 
        sx(ntotal+nvirt1+nvirt2)=3
	  endif
      
	  if(x(1,i).lt.x0.and.(x0-x(1,i)).le.s_rad1*dx1)then
      nvirt2=nvirt2+1
      x(1,ntotal+nvirt1+nvirt2)=2.*x0-x(1,i) 
	x(2,ntotal+nvirt1+nvirt2)=x(2,i)  
	!hsml(ntotal+nvirt1+nvirt2)=hsml(i)
	rho(ntotal+nvirt1+nvirt2)=rho(i)
	mass(ntotal+nvirt1+nvirt2)=dx1*dy1*rho(ntotal+nvirt1+nvirt2)
	p(ntotal+nvirt1+nvirt2)=p(i)
      ! sts(1,ntotal+nvirt1+nvirt2)=sts(1,i)
      !sts(2,ntotal+nvirt1+nvirt2)=sts(2,i)
      !sts(3,ntotal+nvirt1+nvirt2)=sts(3,i)
      !sts(4,ntotal+nvirt1+nvirt2)=sts(4,i)
	vx(1,ntotal+nvirt1+nvirt2)=-2.*0.0-vx(1,i) 
	vx(2,ntotal+nvirt1+nvirt2)=vx(2,i) 
      sx(ntotal+nvirt1+nvirt2)=3
	  endif


      enddo
!      
!       do i=1,nc
!      nvirt3=nvirt3+1
!      x(1,ntotal+nvirt1+nvirt2+nvirt3)=cx3-bj*cos((i-1)*pi/(nc-1))
!	x(2,ntotal+nvirt1+nvirt2+nvirt3)=cy3-bj*sin((i-1)*pi/(nc-1))
!c	hsml(ntotal+nvirt1+nvirt2)=hsml(i)
!	rho(ntotal+nvirt1+nvirt2+nvirt3)=rho0
!c	mass(ntotal+nvirt1+nvirt2)=dx*dy*rho(ntotal+nvirt1+nvirt2)
!       k=0
!      pp=0.0
!      tsts1=0.0
!      tsts2=0.0
!      tsts3=0.0
!      tsts4=0.0
!      n0=0.0
!      x_i=x(1,ntotal+nvirt1+nvirt2+nvirt3)
!            y_i=x(2,ntotal+nvirt1+nvirt2+nvirt3)
!            
!       do j=1,rn(ntotal+nvirt1+nvirt2+nvirt3)
!c           k=k+1
! 	      nj=rpn(j,ntotal+nvirt1+nvirt2+nvirt3)
!            if(nj.le.ntotalb)then
!            x_j=x(1,nj)
!            y_j=x(2,nj)
!             hsml0=hsml
!            call kernel(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
!c            pp=pp+p(nj)*w
!             tsts1=tsts1+tsts(1,nj)*w
!            tsts2=tsts2+tsts(2,nj)*w
!            tsts3=tsts3+tsts(3,nj)*w
!            tsts4=tsts4+tsts(4,nj)*w
!            n0=n0+w
!            endif
!            enddo
!	  if(n0.ne.0)then
!c	p(ntotal+nvirt1+nvirt2)=pp/n0
!      tsts(1,ntotal+nvirt1+nvirt2+nvirt3)=tsts1/n0
!      tsts(2,ntotal+nvirt1+nvirt2+nvirt3)=tsts2/n0
!      tsts(3,ntotal+nvirt1+nvirt2+nvirt3)=tsts3/n0
!      tsts(4,ntotal+nvirt1+nvirt2+nvirt3)=tsts4/n0
!      endif
!	vx(1,ntotal+nvirt1+nvirt2+nvirt3)=0.0
!	vx(2,ntotal+nvirt1+nvirt2+nvirt3)=vs
!      
!      
!       sx(ntotal+nvirt1+nvirt2+nvirt3)=4
!      
!      enddo
!      
      
c      if(mod(t,3).eq.0.or.t.eq.1)then
c      call sp(t,ntotalb) 
      call sp1(t,ntotal+nvirt1+nvirt2) 	!+nvirt3
      
      if(mod(t,10000).eq.0.or.t.eq.1)then
      
      
	write(*,*)t,"sp"

	endif
c       do i=1,ntotal+nvirt1+nvirt2+nvirt3
c	  epsino_i(i)=max(epsino*sqrt(vx(1,i)**2+vx(2,i)**2)/
c     &	          sqrt(2*9.81*yl), epsino)
c	  enddo
	

     
c      if(mod(t,sa_n).eq.0.or.t.eq.1)then
c	open(5,file='freeparticle.dat')
c	write(5,*)'TITLE="INITITAL DISTRIBUTION"'
c	write(5,*)'VARIABLES="X" "Y" "VX" "VY"'
c	write(5,*)'ZONE I=',ntotal,' F=POINT'
c	do i=1,ntotal
c      write(5,'(4e16.6)')x(1,i),x(2,i),vx(1,i),vx(2,i) 
c	enddo
c	close(5)
      
c	open(6,file='bottomparticle.dat')
c	write(6,*)'TITLE="INITITAL DISTRIBUTION"'
c	write(6,*)'VARIABLES="X" "Y" "VX" "VY"'
c	write(6,*)'ZONE I=',nvirt1,' F=POINT'
c	do i=1+ntotal,ntotal+nvirt1
c      write(6,'(4e16.6)')x(1,i),x(2,i),vx(1,i),vx(2,i) 
c	!(5,'(e14.8,'','',e14.8)')
c	enddo
c	close(6)

c	open(7,file='wallparticle.dat')
c	write(7,*)'TITLE="INITITAL DISTRIBUTION"'
c	write(7,*)'VARIABLES="X" "Y" "VX" "VY"'
c	write(7,*)'ZONE I=',nvirt2,' F=POINT'
c	do i=1+ntotal+nvirt1,ntotal+nvirt1+nvirt2
c      write(7,'(4e16.6)')x(1,i),x(2,i),vx(1,i),vx(2,i) 
	!(5,'(e14.8,'','',e14.8)')
c	enddo
c	close(7)


c	endif
!        n_f=0
!       n_f=0
!      do i=1,ntotalb
!        x_i=x(1,i)
!            y_i=x(2,i)
!            ip_f(i)=1
!            ip0=0
!            ip1=0
!           coor_1=0
!	coor_2=0
!	coor_3=0
!	coor_4=0
!
!      coor_x_1=0
!	coor_x_2=0
!	coor_x_3=0
!	coor_x_4=0
!
!            do j=1,rn(i)
! 	      nj=rpn(j,i)
!          if(nj.le.ntotalb)then     
!	      x_j=x(1,nj)
!            y_j=x(2,nj)   
!        dis_i=sqrt((x(1,i)-x(1,nj))**2+(x(2,i)-x(2,nj))**2)  
!       if(dis_i.le.1.75*dx.and.dis_i.gt.0.0)then        
!	!dis_min_r!!!!!!!!!!!!!!!!!!!!!
!       
!c      call kernel(x_i,y_i,x_j,y_j,w,dwdx,dwdy,r0_h)
!
!	sin_sita=(x(1,nj)-x(1,i))/dis_i
!	cos_sita=(x(2,nj)-x(2,i))/dis_i
!	
!	if(cos_sita.gt.0.0.and.sin_sita.ge.0.0)then
!	coor_1=1  !First quadrant
!      
!      endif
!	if(cos_sita.le.0.0.and.sin_sita.gt.0.0)then
!	coor_2=1  !Fourh quadrant
!	endif
!	if(cos_sita.lt.0.0.and.sin_sita.le.0.0)then
!	coor_3=1  !Third quadrant
!	endif  
!	if(cos_sita.ge.0.0.and.sin_sita.lt.0.0)then
!	coor_4=1  !Second quadrant
!	endif      
!
! !     if(j0.gt.ntotal)then
!	!ip_s=ip_s
! !     elseif(j0.le.ntotal.and.j0.ne.i.and.dis_i.le.s_rad_group*dx)then
! !     ip_s=ip_s+1
! !     endif
!	if((x(1,nj)-x(1,i)).gt.0.0
!     &.and.abs(x(2,nj)-x(2,i)).le.(x(1,nj)-x(1,i)))then    
!	coor_x_1=1
!	endif
!	if((x(2,nj)-x(2,i)).gt.0.0
!     &.and.abs(x(1,nj)-x(1,i)).le.(x(2,nj)-x(2,i)))then    
!	coor_x_2=1
!	endif 
!	if((x(1,i)-x(1,nj)).gt.0.0
!     &.and.abs(x(2,nj)-x(2,i)).le.(x(1,i)-x(1,nj)))then     
!	coor_x_3=1	
!	endif
!	if((x(2,i)-x(2,nj)).gt.0.0
!     &.and.abs(x(1,nj)-x(1,i)).le.(x(2,i)-x(2,nj)))then      
!	coor_x_4=1 
!      endif 
!      
!      endif
!      endif
!      enddo
!      
!      ip0=coor_1+coor_2+coor_3+coor_4
!      ip1=coor_x_1+coor_x_2+coor_x_3+coor_x_4
!
!       if(ip0.lt.4.or.ip1.lt.4)then
!         ip_f(i)=0
!         n_f=n_f+1
!         fre_num(n_f)=i
!      !else if((ip0.lt.4.or.ip1.lt.4).and.ljsxyb(i).gt.0.1)then
!      !   ip_f(i)=0
!      !   n_f=n_f+1
!      !   fre_num(n_f)=i
!      endif
!       ! do j=1,rn(i)
! 	     ! nj=rpn(j,i)
!       !   if(nj.gt.ntotalb)then     
!	      !x_j=x(1,nj)
!       !     y_j=x(2,nj)   
!       ! dis_i=sqrt((x(1,i)-x(1,nj))**2+(x(2,i)-x(2,nj))**2)  
!       !if(dis_i.le.1.75*dx.and.dis_i.gt.0.0)then  
!       !    jc(i)=0
!       !endif
!       !endif
!       !    enddo  
!            enddo
!            do i=1,ntotalb
!                 jc(i)=1
!                if(ip_f(i)==0)then
!                  x_i=x(1,i)
!            y_i=x(2,i)   
!       do j=1,rn(i)
! 	      nj=rpn(j,i)
!          if(nj.gt.ntotalb)then     
!	      x_j=x(1,nj)
!            y_j=x(2,nj)   
!        dis_i=sqrt((x(1,i)-x(1,nj))**2+(x(2,i)-x(2,nj))**2)  
!       if(dis_i.le.2.0*dx.and.dis_i.gt.0.0)then  
!           jc(i)=0
!       endif
!       endif
!       enddo  
!       endif
!       enddo
      !write(*,*)n_f


	!call interforce(ti)
       !if(t.eq.1.or.mod(t,17).eq.0)then
       call interforce1(ti)
       !endif
c	write(*,*)t,'interforce'
c	write(*,'(''dvx='',e14.5,''dvy='',e14.5)')dvx(1,67),dvx(2,67)
c      open(1,file='initial_distribution.dat')
	
c	write(1,*)'TITLE="INITITAL DISTRIBUTION"'
c	write(1,*)'VARIABLES="X" "Y"'
c	write(1,*)'ZONE I=5 J=10 F=POINT'
c	do i=1,ntotal+nvirt1+nvirt2
c      write(1,'(e14.8,'','',e14.8,e14.6)')x(1,i),x(2,i),hsml(i) 
c	enddo
c	close(1)
	end