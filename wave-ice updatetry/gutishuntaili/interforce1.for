      subroutine interforce1(ti)
      implicit none 
	include 'menu.inc'
      double precision ti
       integer m,n,nj,tj,zj,zj1,xbs,xbx,xby,xbz,bjs,bys,byx,bzs,bzx
       integer nj1,nj2,bfj,szj,gzj,zj2,zj3,nfs,nfs1,bjs1
      double precision x_i,y_i,xi(dim),xj(dim),vx_i,vy_i,gx
     &,x_j,y_j,vx_j,vy_j,h_ij,rho_ij,c_ij
     &,v_ij_x_ij,phai_ij,hsml0,pai_ij
	double precision w,dwdx,dwdy,p_i,p_j,rho_i,rho_j
      
      double precision sr_i1,sr_i2,sr_i3,sr_i4,tr_i1,tr_i2,tr_i3
     &,tr_i4,ssr_i1,ssr_i2,ssr_i3,ssr_i4,sita,c,s,wij,sg_i1,
     & sg_i2,sg_i3,sg_i4,r_gi1,r_gi2,r_gi3,r_gi4,rgylx,rgyly
     & ,tssr_i1,tssr_i2,tssr_i3,tssr_i4
      double precision I1,J2,fl,xvx1,xvx2,x_ij,y_ij
     & ,pff_i1,pff_i2,pff_i3,pff_i4,psr_i1,psr_i2,psr_i3,psr_i4
     &, tpff_i1,tpff_i2,tpff_i3,tpff_i4
      double precision sts1,sts2,sts3,sts4,pp,bl,wij1
     &,tsts1,tsts2,tsts3,tsts4,ys,rx,rho_i2,rho_nj2,vw1,vw2,betla
      double precision r_1,r_2,L_11,L_12,L_21,L_22,JZ,
     & dwd(dim),cor(dim),mz,dlxi,dlxj,dlyi,dlyj,ldis,ldis1
      
       double precision del,def,ljsxybl,ljdxpyb,massj,eb,eb1,massi
      double precision  n00,c11,c12,c22,y1,y12,y2,y22,n11,n22,zbij,bij,
     & r11,r22,dmin,dis,fdeltx,fdelty,n0,ra,tpff,fij,a1,b1,c1,a2,b2,c2,
     & jdx,jdy,dd2,dd1,ww1,ww2,tdl1,tdl2,dl1,dl2,v_fraction,v_c,rhob
c	n03=0.0
      
c	ys=3.0e+6 !2.0e+6
     
     
       !  jiao=pi/5.8
       !jzj=0.0!jiao/4.

c        arf=tan(jiao)/(9.+12*tan(jiao)**2)**0.5! ! 2*sin(jiao)/3.**0.5/(3.-sin(jiao))
c      a_rf=2*sin(jzj)/3.**0.5/(3.-sin(jzj))! !sin(jiao)/3.**0.5/(3.-sin(jiao))
c      sgm=3./(9.+12.*tan(jiao)**2)**0.5 ! !sin(jiao)/3.**0.5/(3.-sin(jiao))6*cos(jiao)/3.**0.5/(3.-sin(jiao))

c     &  *6*cos(jiao)/3.**0.5/(3.-sin(jiao))
c      write(*,*)qfqd
    
      !k=0
!       
!      do i=ntotalb+1,ntotal
!c	if(bj(i)==0)then
!       ssr_i1=0.0
!       ssr_i2=0.0
!       ssr_i3=0.0
!       ssr_i4=0.0
!       !tr_i1=0.0
!       !tr_i2=0.0
!       !tr_i3=0.0
!       !tr_i4=0.0
!        
!       !sts1=sts(1,i)
!       !sts2=sts(2,i)
!       !sts3=sts(3,i)
!       !sts4=sts(4,i)
!       !  tsts1=tsts(1,i)
!       !tsts2=tsts(2,i)
!       !tsts3=tsts(3,i)
!       !tsts4=tsts(4,i)
!	      x_i=x(1,i)
!            y_i=x(2,i)
!             vx_i=vx(1,i)
!	      vy_i=vx(2,i)
!	     L_11=0.0
!       L_12=0.0
!      L_21=0.0
!       L_22=0.0 
!             do j=1,rn1(i)
! 	      nj=rpn1(j,i)
!          
!            if(nj.gt.ntotalb)then
!	      x_j=x(1,nj)
!            y_j=x(2,nj)                
!      !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
!      !hsml0=1/2.*(hsml+hsmls)
!      !      else
!        hsml0=1/2.*(hsmls+hsmls)
!        !endif
!	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
!  
!      r_1=x_j-x_i
!      r_2=y_j-y_i
!           
!c      r(k)=x(k,j0)
!
!c      dwd(1)=(dwdx-fz(1)/fm)/fm
!c	dwd(2)=(dwdy-fz(2)/fm)/fm 
!     
!      L_11=L_11+r_1*dwdx
!       L_12=L_12+r_1*dwdy
!      L_21=L_21+r_2*dwdx
!       L_22=L_22+r_2*dwdy 
!       endif
!             enddo
!C      L=.i.L
!             JZ=L_11*L_22-L_12*L_21
!             if(jz.ne.0.0)then
!         NL(1,1)=L_22/(JZ*dx1*dx1)
!       NL(1,2)=-L_12/(JZ*dx1*dx1)
!      NL(2,1)=-L_21/(JZ*dx1*dx1)
!       NL(2,2)=L_11/(JZ*dx1*dx1)
!             
!       else
!              NL(1,1)=0.0
!       NL(1,2)=0.0
!       NL(2,1)=-0.0
!       NL(2,2)=0.0    
!       endif
!C       call inv(L,dim,Q,is,js)
!            
!            
!	    do j=1,rn1(i)
! 	      nj=rpn1(j,i)
!           if(nj.gt.ntotalb)then
!            
!	      x_j=x(1,nj)
!            y_j=x(2,nj)
!             vx_j=vx(1,nj)
!	      vy_j=vx(2,nj)
!           
!      !  if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
!      !hsml0=1/2.*(hsml+hsmls)
!      !      else
!        hsml0=1/2.*(hsmls+hsmls)
!        !endif
!	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
!      dwd(1)=dwdx
!	dwd(2)=dwdy
!C	cor=matmul(NL,dwd) 
!C      write(*,*) cor(1),cor(2)
!      !cor(1)=dwdx
!      !cor(2)=dwdy  
!      cor(1)=NL(1,1)*dwd(1)+NL(1,2)*dwd(2)
!      cor(2)=NL(2,1)*dwd(1)+NL(2,2)*dwd(2)      
!      ssr_i1=ssr_i1+(vx_j-vx_i)*cor(1)-2/3.*((vx_j-vx_i)*cor(1)
!     & +(vy_j-vy_i)*cor(2)) !xx
!      ssr_i2=ssr_i2+(vy_j-vy_i)*cor(2)-2/3.*((vx_j-vx_i)*cor(1)
!     & +(vy_j-vy_i)*cor(2)) !xx !yy
!      ssr_i3=ssr_i3+0.5*((vx_j-vx_i)*cor(2)+(vy_j-vy_i)*cor(1))               !xy
!	ssr_i4=ssr_i3 !sr(4,i)+0.5*mass(nj)/rho(nj)*((vx(1,nj)-vx(1,i))*
!        ! &		  dwdx(2)+(vx(2,nj)-vx(2,i))*dwdx(1))               !yx 
!      endif
!          enddo
!      ssr_i1=ssr_i1*dx1*dx1
!      ssr_i2=ssr_i2*dx1*dx1
!      ssr_i3=ssr_i3*dx1*dx1             
!	ssr_i4=ssr_i3  
!      
!               
!      
!    
!     
!       sts(1,i)=mius*ssr_i1
!      sts(2,i)=mius*ssr_i2 !yy
!	sts(3,i)=mius*ssr_i3 !xy
!	sts(4,i)=mius*ssr_i4 !yx
!c      if(mod(ti,10000*dt).eq.0)then
!c      write(*,*) ljsxybl
!c      endif      
!c      if(cs(i)==1)then
!c      pyb(1,i)=pyb(1,i)+psr_i1*dt
!c      pyb(2,i)=pyb(2,i)+psr_i2*dt
!c      pyb(3,i)=pyb(3,i)+psr_i3*dt
!c      pyb(4,i)=pyb(4,i)+psr_i4*dt
!c      ljsxyb(i)=ljsxyb(i)+(ljsxybl)*dt
!c      elseif(cs(i)==0)then
!c         pyb(1,i)=pyb(1,i)+0.0
!c      pyb(2,i)=pyb(2,i)+0.0
!c      pyb(3,i)=pyb(3,i)+0.0
!c      pyb(4,i)=pyb(4,i)+0.0
!c      endif
!c      if(cs(i)==1)then
!c        dxpyb=dxpyb+(pyb(1,i)**2+pyb(2,i)**2+pyb(3,i)**2+pyb(4,i)  !+(ljsxyb(i)
!c     &**2)**0.5
!c      endif
!c      k=k+cs(i)
!      enddo
!      
!
!      do i=1,ntotalb 
!
!      fs(i)=0
!   
!      
!               x_i=x(1,i)
!	      y_i=x(2,i)
!            !if(ip_f(i)==0)then
! 
!       bys=0
!       byx=0
!       bzs=0
!       bzx=0
!       xbs=0
!       xbx=0
!       xbz=0
!       xby=0
!       bjs=0
!       do j=1,rn1(i)
!       nj=rpn1(j,i)
!c        do k=1,dim
!c        xj(k)=x(k,nj)
!c        enddo
!           if(nj.le.ntotalb)then
!       x_j=x(1,nj)
!       y_j=x(2,nj)
!       !z_j=x(3,nj)
!       if(sqrt((x_i-x_j)**2+(y_i-y_j)**2).le.1.8*dx)then
!c       vx_j=vx(1,nj)
!c      vy_j=vx(2,nj)
!      
!c       vx_i=vx(1,i)
!c       vy_i=vx(2,i)
!       r_1=x_i-x_j
!       r_2=y_i-y_j
!       !r_3=z_i-z_j
!       
!       if(r_1.lt.0.0.and.r_2.le.0.0)then
!           bys=1
!       else if(r_1.le.0.0.and.r_2.gt.0.0)then
!           byx=1
!       else if(r_1.ge.0.0.and.r_2.lt.0.0)then
!           bzs=1
!       else if(r_1.gt.0.0.and.r_2.ge.0.0)then
!           bzx=1
!       endif
!           
!           
!      if((x_j-x_i).gt.0.0
!     &.and.1.0*abs(y_j-y_i).le.(x_j-x_i))then
!	xby=1
!	endif
!	if((y_j-y_i).gt.0.0
!     &.and.1.0*abs(x_j-x_i).le.(y_j-y_i))then
!	xbs=1
!	endif 
!	if((x_i-x_j).gt.0.0
!     &.and.1.0*abs(y_j-y_i).le.(x_i-x_j))then
!	xbz=1
!	endif
!	if((y_i-y_j).gt.0.0
!     &.and.1.0*abs(x_j-x_i).le.(y_i-y_j))then   
!	xbx=1
!      endif 
!       endif
!        endif
!         if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
!       x_j=x(1,nj)
!       y_j=x(2,nj)  
!       if(sqrt((x_i-x_j)**2+(y_i-y_j)**2).le.3.0*dx)then
!           bjs=1
!           endif
!           endif
!       enddo
!       
!      if(((bys+byx+bzs+bzx).le.3.or.(xbs+xbx+xbz+xby).le.3).and.bjs==1)
!     & then
!           fs(i)=1
!       endif
!     !!  if(((xbs+xbx+xbz+xby).le.2).and.bjs==1).and.bjs==1
!     !!& then
!     !!      zbj(i)=2
!     !!  endif
!      enddo

!      !$omp parallel private(i,j,nj,x_i,y_i,x_j,y_j)!if(ntotalb>2000) 
!     & private( bys,byx,bzs,bzx,xbs,xbx,xbz,xby,bjs,bjs1,r_1,r_2)
!     & private(jz,w,dwdx,dwdy,hsml0,cor,c_ij)
!     & private(pp,n0)
!     !!& private(tssr_i1,tssr_i2,tssr_i3,tssr_i4,tssr_i5,tssr_i6)  
!     !!& private(tssr_i7,tssr_i8,tssr_i9)
!           
!      !$omp do
!      
!      
!       do i=ntotalb+1,ntotal 
!
!      fs(i)=0
!   
!       if(x(1,i).ge.10.0.and.x(1,i).le.16.8)then
!               x_i=x(1,i)
!	      y_i=x(2,i)
!            !if(ip_f(i)==0)then
! 
!       bys=0
!       byx=0
!       bzs=0
!       bzx=0
!       xbs=0
!       xbx=0
!       xbz=0
!       xby=0
!       bjs=0
!       bjs1=0
!       do j=1,rn1(i)
!       nj=rpn1(j,i)
!c        do k=1,dim
!c        xj(k)=x(k,nj)
!c        enddo
!      if(nj.gt.ntotalb.and.nj.le.(ntotal))then
!       x_j=x(1,nj)
!       y_j=x(2,nj)
!        if(sqrt((x_i-x_j)**2+(y_i-y_j)**2).le.1.8*dx1)then
!       !z_j=x(3,nj)
!       
!c       vx_j=vx(1,nj)
!c      vy_j=vx(2,nj)
!      
!c       vx_i=vx(1,i)
!c       vy_i=vx(2,i)
!       r_1=x_i-x_j
!       r_2=y_i-y_j
!       !r_3=z_i-z_j
!       
!       if(r_1.lt.0.0.and.r_2.le.0.0)then
!           bys=1
!       else if(r_1.le.0.0.and.r_2.gt.0.0)then
!           byx=1
!       else if(r_1.ge.0.0.and.r_2.lt.0.0)then
!           bzs=1
!       else if(r_1.gt.0.0.and.r_2.ge.0.0)then
!           bzx=1
!       endif
!           
!           
!      if((x_j-x_i).gt.0.0
!     &.and.1.0*abs(y_j-y_i).le.(x_j-x_i))then
!	xby=1
!	endif
!	if((y_j-y_i).gt.0.0
!     &.and.1.0*abs(x_j-x_i).le.(y_j-y_i))then
!	xbs=1
!	endif 
!	if((x_i-x_j).gt.0.0
!     &.and.1.0*abs(y_j-y_i).le.(x_i-x_j))then
!	xbz=1
!	endif
!	if((y_i-y_j).gt.0.0
!     &.and.1.0*abs(x_j-x_i).le.(y_i-y_j))then   
!	xbx=1
!      endif 
!       endif
!       endif
!         if(nj.le.ntotalb)then
!       x_j=x(1,nj)
!       y_j=x(2,nj)  
!       if(sqrt((x_i-x_j)**2+(y_i-y_j)**2).le.3.0*dx1)then
!           bjs=1
!       endif
!      if(sqrt((x_i-x_j)**2+(y_i-y_j)**2).le.1.1*dx1)then
!           bjs1=1
!           endif
!         endif
!        
!       enddo
!       !
!      if((((bys+byx+bzs+bzx).le.3.or.(xbs+xbx+xbz+xby).le.3).and.bjs==1)
!     & .or.bjs1==1) then
!           fs(i)=1
!       endif
!     !!  if(((xbs+xbx+xbz+xby).le.2).and.bjs==1)
!     !!& then
!     !!      zbj(i)=2
!     !!  endif
!       endif
!       enddo
!         !$omp end do 
!      
!      
!      !$omp end parallel 

      
 !   	do i=ntotalb+1,ntotal
	!  !p_n=p_n+mass*9.81*x(2,i)
 !       !rho_t=rho_t+rho(i)
 !         
	!  v_c=0
 !        fs(i)=0
 !        if(x(1,i).ge.10.0.and.x(1,i).le.16.8)then
	!  do j=1,rn1(i)
 !          nj=rpn1(j,i)
 !           !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
 !           if(nj.gt.ntotalb.and.nj.le.(ntotal))then !+nvirt1+nvirt2
 !       v_c=v_c+dx*dy
 !       endif
	!  enddo
 !       !p_n=p_n+0.5*mass*(vx(1,i)**2+vx(2,i)**2)
	!v_fraction=v_c/(pi*(s_rad1*dx)**2)
 !
	!if(v_fraction.le.0.75)then
 !     fs(i)=1
	!!i_s_f=i_s_f+1
	!!free_surface(i_s_f)=i
 !     endif
 !     endif
 !       enddo
 !   
 !       	do i=1,ntotalb
	!  !p_n=p_n+mass*9.81*x(2,i)
 !       !rho_t=rho_t+rho(i)
	!  v_c=0
 !        fs(i)=0
	!  do j=1,rn(i)
 !          nj=rpn(j,i)
 !           !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
 !           if(nj.le.ntotalb)then !+nvirt1+nvirt2
 !       v_c=v_c+dx*dy
 !       endif
	!  enddo
 !       !p_n=p_n+0.5*mass*(vx(1,i)**2+vx(2,i)**2)
	!v_fraction=v_c/(pi*(s_rad*dx)**2)
 !
	!if(v_fraction.le.0.65)then
 !     fs(i)=1
	!!i_s_f=i_s_f+1
	!!free_surface(i_s_f)=i
	!endif
	!enddo

       mz=0.0 
      !$omp parallel private(i,j,m,nj1,nfs,nfs1,zj,zj1,zj2,zj3,x_i,y_i)!if(ntotalb>2000) 
     & private(x_j,y_j,p_i,p_j,rho_i,massi,nj,rho_j,massj,x_ij,y_ij)
     & private(jz,w,dwdx,dwdy,hsml0,cor,c_ij,dis,dlxi,dlyi,dlxj,dlyj)
     & private(pp,n0,xvx1,xvx2,dd2,dd1,ww1,ww2,tdl1,tdl2,dl1,dl2)
     & private(ldis,ldis1,vw1,vw2,betla,rhob)
     !!& private(tssr_i1,tssr_i2,tssr_i3,tssr_i4,tssr_i5,tssr_i6)  
     !!& private(tssr_i7,tssr_i8,tssr_i9)
           
      !$omp do schedule(static)   
       do i=1,ntotalb 
       zj=0         
       k=0
      pp=0.0
      xvx1=0.0
      xvx2=0.0
      vw1=0.0
      vw2=0.0
      rhob=0.0
      !sts1=0.0
      !sts2=0.0
      !sts3=0.0
      !sts4=0.0
      n0=0.0  
 !     c11=0.0
	!c22=0.0
	!n11=0.0
	!n22=0.0
	!zbij=0.0
 !     r11=0.0
	!r22=0.0
 !     fij=0.0
 !     wij=0.0
 !     bfj=0
 !     zjbj(i)=0
      
               x_i=x(1,i)
	      y_i=x(2,i)

!          do j=1,rn1(i)
!       nj=rpn1(j,i)
!c        do k=1,dim
!c        xj(k)=x(k,nj)
!c        enddo
!           if(nj.le.ntotalb)then
!       x_j=x(1,nj)
!       y_j=x(2,nj)
!       !z_j=x(3,nj)
!       
!       dis=sqrt((x(1,ntotal+i)-sjx1(1,j))**2+(
!     & x(3,ntotal+i)-sjx1(3,j))**2)
!       if(dis.le.ldis)then
!           ldis=dis
!           vzj=j
!       endif
!          enddo     
      !      nfs=0
      !do m=1,rn1(i)
      !nj1=rpn1(m,i)
      !     if(nj1.le.ntotalb.and.fs(nj1)==1)then 
      !    nfs=nfs+1
      !     endif
      !     enddo
           

           do j=1,rn1(i)
 	      nj=rpn1(j,i)
            !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
      if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
	      x_j=x(1,nj)
            y_j=x(2,nj) 
            
        
                                
       !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
      hsml0=1/2.*(hsmls+hsmls)
      !else
      !  hsml0=1/2.*(hsmls+hsml)
      !  endif
      !call kernelw(x(1,zj),x(2,zj),x_j,y_j,w,dwdx,dwdy,hsml0)
      !wij=w
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0) 
     !!  if(sqrt((x_i-x_j)**2+
     !!& (y_i-y_j)**2).le.dx)then  
     !!             dis=sqrt((x_i-x_j)**2+(y_i-y_j)**2)
     !! dl1=0.5*dis
     !! dl2=0.5*dis
     !! dlxi=dl2*(x(1,nj)-x(1,i))/dis
     !! dlyi=dl2*(x(2,nj)-x(2,i))/dis
     !! dlxj=dl1*(x(1,nj)-x(1,i))/dis
     !! dlyj=dl1*(x(2,nj)-x(2,i))/dis
     !!  !write(*,*) i,x(2,i),dl2/dis,dl1/dis,dis,abs(y_j-y_i)
     !!  pp=pp+p(nj)*w+((g+dvx(2,i)))*rho01*!+dvx(2,i)g++ddvx(2,i)+ddvx(2,i)
     !!&(dlyj+rho0/rho01*dlyi)*w+rho01*(dlxj+rho0/rho01*dlxi)*
     !!&(dvx(1,i)) *w         !-dvx(1,nj)
     !!!!   pp=pp+p(nj)*w+((g+rho0/rho01*dvx(2,i)))*rho01*(x(2,nj)-!+dvx(2,i)g++ddvx(2,i)+ddvx(2,i)
     !!
     !!       betla=1.0+dl2/dl1
     !!       
     !! !betla=min(5.,1.0+dl2/dl1)
     !! !write(*,*) betla,w
     !!!!  if(sqrt((x(1,zj)-x_ij)**2+(x(2,zj)-y_ij)**2).le.
     !!!!&  sqrt((x(1,zj1)-x_ij)**2+(x(2,zj1)-y_ij)**2))then
     !! vw1=vw1+ betla*vx(1,i)*w
     !! vw2=vw2+ betla*vx(2,i)*w
     !! ! else
     !! ! vw1=vw1+ betla*vx(1,zj1)*w
     !! !vw2=vw2+ betla*vx(2,zj1)*w
     !! ! endif
     !! !  vw1=vw1+ betla*vx(1,i)*w
     !! !vw2=vw2+ betla*vx(2,i)*w
     !! xvx1=xvx1+(1.-betla)*vx(1,nj)*w!*mass(nj)/rho(nj)
     !! xvx2=xvx2+(1.-betla)*vx(2,nj)*w!*mass(nj)/rho(nj)
     !! !xvx3=xvx3+vx(3,nj)*w*mass(nj)/rho(nj)
     !!       n0=n0+w
     !!    else
     !!      
           
      
     !! nfs1=0
     !! do m=1,rn1(nj)
     !! nj1=rpn1(m,nj)
     !!  if(nj1.gt.ntotalb.and.nj1.le.(ntotal).and.fs(nj1)==1)then           
     !! nfs1=nfs1+1
     !!      endif
     !! enddo  
     !! 
     !! 
     !!  if(nfs.le.1.and.fs(i)==1)then
     !!   zj=i
     !!   zj1=i
     !!   
     !! elseif(nfs.le.1.and.fs(i).ne.1)then
     !!     zj=i 
     !!   
     !!  do m=1,rn1(i)
     !! nj1=rpn1(m,i)
     !!  if(nj1.le.ntotalb.and.fs(nj1)==1.and.nj1.ne.i)then 
     !!  !z_j=x(3,nj)
     !!      zj1=nj1
     !! endif
     !!  enddo 
     !!   elseif(nfs.ge.2)then
     !!     
     !! ldis=5.*dx1
     !!     ldis1=5.*dx1
     !!  do m=1,rn1(i)
     !! nj1=rpn1(m,i)
     !!      if(nj1.le.ntotalb.and.fs(nj1)==1)then 
     !!          
     !!  !z_j=x(3,nj) 
     !!  dis=sqrt(((x_i+x_j)/2.-x(1,nj1))**2+
     !!& ((y_i+y_j)/2.-x(2,nj1))**2)
     !!  if(dis.le.ldis)then
     !!      ldis=dis
     !!      zj=nj1
     !! endif
     !! endif
     !!  enddo 
     !!  !!write(*,*)"h",zbj(1),zbj(4),zbj(9)
     !!  do n=1,rn1(i)
     !!  nj1=rpn1(n,i)
     !!  !!write(*,*)"h",i,nj2
     !! if(nj1.le.ntotalb.and.fs(nj1)==1.and.nj1.ne.zj)then 
     !! !if(nj2.le.ntotalb.and.zbj(nj2)==1.and.nj2.ne.zj2)then 
     !!!!  !z_j=x(3,nj) 
     !!  dis=sqrt(((x_i+x_j)/2.-x(1,nj1))**2+
     !!& ((y_i+y_j)/2.-x(2,nj1))**2)
     !!  if(dis.le.ldis1.and.dis.gt.ldis)then
     !!      ldis1=dis
     !!      zj1=nj1
     !!  endif 
     !!  endif
     !!  enddo   
     !!  
     !! endif
     !! 
     !! if(nfs1.le.1.and.fs(nj)==1)then
     !!   zj2=nj
     !!   zj3=nj
     !!  elseif(nfs1.le.1.and.fs(nj).ne.1)then
     !!     zj2=nj 
     !!
     !!  do m=1,rn1(i)
     !! nj1=rpn1(m,i)
     !!  if(nj1.le.ntotalb.and.fs(nj1)==1.and.nj1.ne.nj)then 
     !!  !z_j=x(3,nj)
     !!      zj3=nj1
     !! endif
     !!  enddo 
     !!   elseif(nfs1.ge.2)then
     !! 
     !!  ldis=5.*dx1
     !!     ldis1=5.*dx1
     !!  do m=1,rn1(nj)
     !! nj1=rpn1(m,nj)
     !!  if(nj1.gt.ntotalb.and.nj1.le.(ntotal).and.fs(nj1)==1)then 
     !!          
     !!  !z_j=x(3,nj) 
     !!  dis=sqrt(((x_i+x_j)/2.-x(1,nj1))**2+
     !!& ((y_i+y_j)/2.-x(2,nj1))**2)
     !!  if(dis.le.ldis)then
     !!      ldis=dis
     !!      zj2=nj1
     !! endif
     !! endif
     !!  enddo 
     !!  !write(*,*)"h",zbj(1),zbj(4),zbj(9)
     !!  do n=1,rn1(nj)
     !!  nj1=rpn1(n,nj)
     !!  !write(*,*)"h",i,nj2
     !! if(nj1.gt.ntotalb.and.nj1.le.(ntotal).and.fs(nj1)==1
     !!&      .and.nj1.ne.zj2)then 
     !! !if(nj2.le.ntotalb.and.zbj(nj2)==1.and.nj2.ne.zj2)then 
     !!!!  !z_j=x(3,nj) 
     !!  dis=sqrt(((x_i+x_j)/2.-x(1,nj1))**2+
     !!& ((y_i+y_j)/2.-x(2,nj1))**2)
     !!  if(dis.le.ldis1.and.dis.gt.ldis)then
     !!      ldis1=dis
     !!      zj3=nj1
     !!  endif 
     !!  endif
     !! enddo   
     !!  endif
     !!!!
     !!!! dd2=0.0
     !!!! dd1=0.0
     !!!! ww1=0.0
     !!!! ww2=0.0
     !!  
     !!!!  do m=1,rn1(nj)
     !!!! nj1=rpn1(m,nj)
     !!!!  if(fs(nj1)==1)then
     !!!!  if(sqrt((x_i-x(1,nj1))**2+(y_i-x(2,nj1))**2).le.s_rad1*dx)then
     !!!! dd1=dd1+sqrt((x_j-x(1,nj1))**2+(y_j-x(2,nj1))**2)
     !!!! ww1=ww1+1.0
     !!!!  endif）
     !!!!  endif
     !!!! enddo
     !!!! 
     !!!!   do m=1,rn1(i)
     !!!! nj1=rpn1(m,i)
     !!!! if(fs(nj1)==1)then
     !!!!  if(sqrt((x_j-x(1,nj1))**2+(y_j-x(2,nj1))**2).le.s_rad1*dx)then
     !!!! dd2=dd2+sqrt((x_i-x(1,nj1))**2+(y_i-x(2,nj1))**2)
     !!!! ww2=ww2+1.0
     !!!! endif
     !!!! endif
     !!!!   enddo
     !!!!   
     !!!!   
     !!!! tdl1=dd1/ww1
     !!!! tdl2=dd2/ww2
     !! 
     !! dis=sqrt((x_i-x_j)**2+(y_i-y_j)**2)
     !!  x_ij=(x(1,zj)+x(1,zj1))/4.+(x(1,zj2)+x(1,zj3))/4.
     !! y_ij=(x(2,zj)+x(2,zj1))/4.+(x(2,zj2)+x(2,zj3))/4.
     !!   tdl1=sqrt((x_j-x_ij)**2+(y_j-
     !!& y_ij)**2) 
     !!  tdl2=sqrt((x_i-x_ij)**2+(y_i-
     !!& y_ij)**2) 
     !!  
     !! 
     !! dl1=tdl1*dis/(tdl1+tdl2)
     !! dl2=tdl2*dis/(tdl1+tdl2)
     !! dlxi=dl2*(x(1,nj)-x(1,i))/dis
     !! dlyi=dl2*(x(2,nj)-x(2,i))/dis
     !! dlxj=dl1*(x(1,nj)-x(1,i))/dis
     !! dlyj=dl1*(x(2,nj)-x(2,i))/dis
       !write(*,*) i,x(2,i),dl2/dis,dl1/dis,dis,abs(y_j-y_i)
       !!  pp=pp+p(nj)*w+((g+dvx(2,i)))*rho01*!+dvx(2,i)g++ddvx(2,i)+ddvx(2,i)
     !!&(dlyj+rho0/rho01*dlyi)*w+rho01*(dlxj+rho0/rho01*dlxi)*(dvx(1,i))
     !!& *w         !-dvx(1,nj)
        pp=pp+p(nj)*w+((g+rho0/rho01*dvx(2,i)))*rho01*(x(2,nj)-!+dvx(2,i)g++ddvx(2,i)+ddvx(2,i)
     &x(2,i))*w+rho01*(x(1,nj)-x(1,i))*(rho0/rho01*dvx(1,i))
     & *w         !-dvx(1,nj)   
     !!     
     !!        !pp=pp+p(nj)
     !!       ! sts1=sts1+sts(1,nj)*w
     !!       !sts2=sts2+sts(2,nj)*w
     !!       !sts3=sts3+sts(3,nj)*w
     !!       !sts4=sts4+sts(4,nj)*w
      xvx1=xvx1+vx(1,nj)*w!*mass(nj)/rho(nj)
      xvx2=xvx2+vx(2,nj)*w!*mass(nj)/rho(nj)
      !xvx3=xvx3+vx(3,nj)*w*mass(nj)/rho(nj)
      !!rhob=rhob+rho(nj)*w
      
            n0=n0+w
            endif
           enddo
      !write(*,*) i,zj,w,n0,bij,zbij,fij,pp
       !p(i)=pp  !/n0
 !       !p(i)=-(tsts(1,i)+tsts(2,i))/3.
 !      !if(i==2)then
 !      !    p(2)=pp/n0
 !      !           write(*,*) p(2),pp,k
 !      !           
 !      !           endif
	  if(n0.ne.0.0)then
            !write(*,*) "h1"
	p(i)=pp/n0
      bxvx(1,i)=2.*vx(1,i)-xvx1/n0
      bxvx(2,i)=2.*vx(2,i)-xvx2/n0
      !!rhox(i)=rhob/n0
      !p(i)=pp/k
       !if(i==1)then
       !           write(*,*) p(1),pp,k
       !           endif
      !ists(1,i)=sts1/n0
      !ists(2,i)=sts2/n0
      !ists(3,i)=sts3/n0
      !ists(4,i)=sts4/n0
            !else
            ! p(i)=0.0  
            endif
            !if(p(i).lt.0.0)then
            !    p(i)=0.0
            !    endif
       !endif
      !       if(ip_f(i)==0)then
      !      mz=mz+mass(i)
      !endif
      
       rhox(i)=rho01*(1+omiga*p(i)/
     &(c01*c01)/rho01)**(1./omiga)
      !write(*,*) xrho(1)
      enddo
        !$omp end do 
      
      
      !$omp end parallel 

      
      
       !$omp parallel private(i,j,x_i,y_i,vx_i,vy_i)!if(ntotalb>2000) 
     & private(L_11,L_12,L_21,L_22,p_i,p_j)
     & private(rho_i,massi,x_j,y_j,nj,rho_j,massj)  
     & private(r_1,r_2,jz,w,dwdx,dwdy,hsml0,cor,c_ij)
     & private(fdeltx,fdelty,eb1)
     !!& private(tssr_i1,tssr_i2,tssr_i3,tssr_i4,tssr_i5,tssr_i6)  
     !!& private(tssr_i7,tssr_i8,tssr_i9)
           
      !$omp do schedule(static)
            	  !!!delt密度修正项111111！！！111
      do i=ntotalb+1,ntotal!ntotalb+
c           if(bj(i)==0)then
           !dpx(i)=0.0
           ! dpy(i)=0.0

           rho_i=rho(i)
            x_i=x(1,i)
	      y_i=x(2,i)
	      vx_i=vx(1,i)
	      vy_i=vx(2,i)
c	rho(i)=rho(i)/n03(i)
             L_11=0.0
       L_12=0.0
      L_21=0.0
       L_22=0.0
             do j=1,rn1(i)
 	      nj=rpn1(j,i)
      !!if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
	      x_j=x(1,nj)
            y_j=x(2,nj)                
            
      !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
      !hsml0=1/2.*(hsmls+hsmls)
      ! call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
      !else
         hsml0=hsmls
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)  
            !endif    
  
      r_1=x_j-x_i
      r_2=y_j-y_i
           
c      r(k)=x(k,j0)

c      dwd(1)=(dwdx-fz(1)/fm)/fm
c	dwd(2)=(dwdy-fz(2)/fm)/fm 
     
      L_11=L_11+r_1*dwdx
       L_12=L_12+r_1*dwdy
      L_21=L_21+r_2*dwdx
       L_22=L_22+r_2*dwdy 
        !endif
             enddo
C      L=.i.L
             JZ=L_11*L_22-L_12*L_21
            if(jz.ne.0.0)then
      NL1(1,i)=L_22/(JZ*dx1*dx1)
      NL1(2,i)=-L_12/(JZ*dx1*dx1)
      NL2(1,i)=-L_21/(JZ*dx1*dx1)
      NL2(2,i)=L_11/(JZ*dx1*dx1)
             
       else
      NL1(1,i)=1.0
      NL1(2,i)=0.0
      NL2(1,i)=0.0
      NL2(2,i)=1.0
       endif
      enddo
            !$omp end do 
      
      !$omp end parallel
      
      !$omp parallel private(i,j,x_i,y_i,vx_i,vy_i)!if(ntotalb>2000) 
     & private(L_11,L_12,L_21,L_22,p_i,p_j)
     & private(rho_i,massi,x_j,y_j,nj,rho_j,massj)  
     & private(r_1,r_2,jz,w,dwdx,dwdy,hsml0,cor,c_ij)
     & private(fdeltx,fdelty,eb1)
     !!& private(tssr_i1,tssr_i2,tssr_i3,tssr_i4,tssr_i5,tssr_i6)  
     !!& private(tssr_i7,tssr_i8,tssr_i9)
           
      !$omp do schedule(static)
      	  !!!delt密度修正项111111！！！111
      do i=ntotalb+1,ntotal!ntotalb+
c           if(bj(i)==0)then
           dpx(i)=0.0
            dpy(i)=0.0

           rho_i=rho(i)
            x_i=x(1,i)
	      y_i=x(2,i)
	      vx_i=vx(1,i)
	      vy_i=vx(2,i)
c	rho(i)=rho(i)/n03(i)
            
       do j=1,rn1(i)
              rgylx=0.0
           rgyly=0.0
 	      nj=rpn1(j,i)
      if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then    
        x_j=x(1,nj)
	      y_j=x(2,nj)
        rho_j=rho(nj)     
       if(nj.le.ntotalb)then
            rho_j=rhox(nj)
             else
              rho_j=rho(nj) 
             endif
      !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
      !hsml0=1/2.*(hsmls+hsmls)
      ! call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
      !else
         hsml0=hsmls
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)  
            !endif      
	     
	    ! dwd(1)=dwdx
	    !dwd(2)=dwdy
C	cor=matmul(NL,dwd) 
C      write(*,*) cor(1),cor(2)
      !cor(1)=dwdx
      !cor(2)=dwdy  
      cor(1)=NL1(1,i)*dwdx+NL1(2,i)*dwdy
      cor(2)=NL2(1,i)*dwdx+NL2(2,i)*dwdy  
      
      
      dpx(i)=dpx(i)+(rho_j-rho_i)*cor(1)*dx1*dx1
      dpy(i)=dpy(i)+(rho_j-rho_i)*cor(2)*dx1*dx1
        endif
            enddo
      
       
       enddo
            !$omp end do 
      
      !$omp end parallel
    
      
       x_j=0.0+dx1
       y_j=0.0
       hsml0=hsmls
      call kernel_f(x_i,y_i,x_j,y_j,w,hsml0)
      wij=w
      dyl=0.0
      cdyl=0.0
      dyva=0.0
        !$omp parallel private(i,j,sts1,sts2,sts3,sts4)!if(ntotalb>2000) 
     & private(p_i,p_j,rho_i,rho_j,rho_i2)
     & private(x_i,y_i,vx_i,vy_i,massi,v_ij_x_ij,phai_ij) 
     & private(x_j,y_j,vx_j,vy_j,nj,rho_nj2,massj,pai_ij)  
     & private(tsts1,tsts2,tsts3,tsts4)
     & private(h_ij,rho_ij,w,dwdx,dwdy,hsml0,cor,c_ij)
     & private(fdeltx,fdelty,eb1)
     !!& private(tssr_i1,tssr_i2,tssr_i3,tssr_i4,tssr_i5,tssr_i6)  
     !!& private(tssr_i7,tssr_i8,tssr_i9)
           
      !$omp do schedule(static)
      
      
      
       do i=ntotalb+1,ntotal
           tj=0
c           if(bj(i)==0)then
           dvx(1,i)=0.0
            dvx(2,i)=0.0

            drho(i)=0.0
            va(1,i)=0.0  
           va(2,i)=0.0  
           p_i=p(i)
           rho_i=rho(i)
         sts1=sts(1,i)
       sts2=sts(2,i)
       sts3=sts(3,i)
       sts4=sts(4,i)
       !pp=p(i)
c            n01=0.0
c            n02=0.0
            x_i=x(1,i)
	      y_i=x(2,i)
	      vx_i=vx(1,i)
	      vy_i=vx(2,i)
c	rho(i)=rho(i)/n03(i)
             L_11=0.0
       L_12=0.0
      L_21=0.0
       L_22=0.0
!             do j=1,rn1(i)
! 	      nj=rpn1(j,i)
!            
!	      x_j=x(1,nj)
!            y_j=x(2,nj)                
!       if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
!      hsml0=1/2.*(hsml+hsmls)
!            else
!        hsml0=1/2.*(hsmls+hsmls)
!        endif
!	call kernel_df(x_i,y_i,x_j,y_j,dwdx,dwdy,hsml0)
!  
!      r_1=x_j-x_i
!      r_2=y_j-y_i
!           
!c      r(k)=x(k,j0)
!
!c      dwd(1)=(dwdx-fz(1)/fm)/fm
!c	dwd(2)=(dwdy-fz(2)/fm)/fm 
!     
!      L_11=L_11+r_1*dwdx
!       L_12=L_12+r_1*dwdy
!      L_21=L_21+r_2*dwdx
!       L_22=L_22+r_2*dwdy 
!
!             enddo
!C      L=.i.L
!             JZ=L_11*L_22-L_12*L_21
!            if(jz.ne.0.0)then
!         NL(1,1)=L_22/(JZ*dx*dx)
!       NL(1,2)=-L_12/(JZ*dx*dx)
!      NL(2,1)=-L_21/(JZ*dx*dx)
!       NL(2,2)=L_11/(JZ*dx*dx)
!             
!       else
!              NL(1,1)=0.0
!       NL(1,2)=0.0
!       NL(2,1)=-0.0
!       NL(2,2)=0.0    
!       endif
       
      
       
       
C       write(*,*) JZ
C       call inv(L,dim,Q,is,js)
          rho_i2=1.0/(rho_i*rho_i)     
	   do j=1,rn1(i)
              rgylx=0.0
           rgyly=0.0
           
 	      nj=rpn1(j,i)
            
            p_j=p(nj)
             if(nj.le.ntotalb)then
            rho_j=rho(nj)
             else
              rho_j=rho(nj) 
             endif
             
              !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
            !fp1(j,i)=0.0
            !fp2(j,i)=0.0
            !endif
          rho_nj2=1.0/(rho_j*rho_j) 
      !calculating the artificial viscosity
	      x_j=x(1,nj)
	      y_j=x(2,nj)
            
       !      if(nj.le.ntotalb)then
	      !vx_j=bxvx(1,nj)
	      !vy_j=bxvx(2,nj)
       !     else
           vx_j=vx(1,nj)
	      vy_j=vx(2,nj)
            !endif
            
           if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
     
            h_ij=1/2.*(hsmls+hsmls)
	      rho_ij=0.5*(rho_i+rho_j)
            
	     c_ij=0.5*(c01+c01)
                 
	      v_ij_x_ij=(vx_i-bxvx(1,nj))*(x_i-x_j)+(vy_i-bxvx(2,nj))*(y_i-y_j)
            
	      phai_ij=h_ij*(v_ij_x_ij)/
     &              ((x_i-x_j)**2+(y_i-y_j)**2+(0.1*h_ij)**2)
            if(v_ij_x_ij.lt.0.0)then
      !pai_ij=(-alpha_phai1*c_ij*phai_ij+belda_phai1*phai_ij**2)/rho_ij 
      pai_ij=(-alpha_phai1*c_ij*phai_ij+belda_phai1*phai_ij**2)/rho_ij          
     !!     pai_ij=(-(alpha_phai+alpha_phai1)/2.*c_ij*phai_ij+(belda_phai+
     !!&	  belda_phai1)/2.*phai_ij**2)/rho_ij
            else
	      pai_ij=0.0
            endif
      hsml0=1/2.*(hsmls+hsmls)
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
	     
	    ! dwd(1)=dwdx
	    !dwd(2)=dwdy
C	cor=matmul(NL,dwd) 
C      write(*,*) cor(1),cor(2)
      !!cor(1)=dwdx
      !!cor(2)=dwdy  
       cor(1)=NL1(1,i)*dwdx+NL1(2,i)*dwdy
      cor(2)=NL2(1,i)*dwdx+NL2(2,i)*dwdy 
     !!! if(((x_i-x_j)**2+(y_i-y_j)**2).ne.0.0)then
     !!! fdeltx=2.*(rhox(nj)-rho(i))*(x_j-x_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     !!!& -(dpx(i)+dpx(nj))
     !!! fdelty=2.*(rhox(nj)-rho(i))*(y_j-y_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     !!!& -(dpy(i)+dpy(nj))
     !!! else
     !!!     fdeltx=0.0
     !!!     fdelty=0.0
     !!! endif
      fdeltx=0.0
          fdelty=0.0
      elseif(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
        h_ij=hsmls
	      rho_ij=0.5*(rho_i+rho_j)
            
	      c_ij=0.5*(c01+c01)
                 
	      v_ij_x_ij=(vx_i-vx_j)*(x_i-x_j)+(vy_i-vy_j)*(y_i-y_j)
            
	      phai_ij=h_ij*(v_ij_x_ij)/
     &    ((x_i-x_j)**2+(y_i-y_j)**2+(0.1*h_ij)**2)
            if(v_ij_x_ij.lt.0.0)then
		  pai_ij=(-alpha_phai1*c_ij*phai_ij+belda_phai1*phai_ij**2)
     &	  /rho_ij
            else
	      pai_ij=0.0
            endif
               hsml0=hsmls
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
	     
	    ! dwd(1)=dwdx
	    !dwd(2)=dwdy
C	cor=matmul(NL,dwd) 
C      write(*,*) cor(1),cor(2)
      !!cor(1)=dwdx
      !!cor(2)=dwdy  
      cor(1)=NL1(1,i)*dwdx+NL1(2,i)*dwdy
      cor(2)=NL2(1,i)*dwdx+NL2(2,i)*dwdy  
      
     
      if(((x_i-x_j)**2+(y_i-y_j)**2).ne.0.0)then
      fdeltx=2.*(rho_j-rho_i)*(x_j-x_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     & -(dpx(i)+dpx(nj))
      fdelty=2.*(rho_j-rho_i)*(y_j-y_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     & -(dpy(i)+dpy(nj))
      else
          fdeltx=0.0
          fdelty=0.0
      endif
      endif

                 
      !cor(1)=NL(1,1)*dwd(1)+NL(1,2)*dwd(2)
      !cor(2)=NL(2,1)*dwd(1)+NL(2,2)*dwd(2)
           
      !cor(1)=NL(1,1)*dwd(1)+NL(1,2)*dwd(2)
      !cor(2)=NL(2,1)*dwd(1)+NL(2,2)*dwd(2)
   
      !if(nj.le.ntotalb)then
      !massj=rho_j*dx*dx
      !massi=rho_i*dx*dx
      
      massj=mass(nj)
      massi=mass(i)
      !else
      !    massj=mass(nj)
      !    endif
C      write(*,*) cor(1),cor(2),DWDX,DWDY
c      n01=n01+dx*dx*(x_j-x_i)*cor(1)
c           n02=n02+dx*dx*(y_j-y_i)*cor(2)
c	    if(mod(i,y_n).eq.0)then
c	      drho(i)=0.
c	    else
	    
c	    endif
c         if(art_vis.eq.1)then
c	dvx(1,i)=dvx(1,i)-((dl(1,i)/rho(i)**2+dl(1,nj)/rho(nj)**2)*dwdx(1)
c     &+(dl(3,i)/rho(i)**2+dl(3,nj)/rho(nj)**2)*dwdx(2))
c     &*mass(nj)-2*pai_ij*dwdx(1)*mass(nj)
c	 vx(2,i)=dvx(2,i)-((dl(4,i)/rho(i)**2+dl(4,nj)/rho(nj)**2)*dwdx(1)
c     &+(dl(2,i)/rho(i)**2+dl(2,nj)/rho(nj)**2)*dwdx(2))
c     &*mass(nj)-2*pai_ij*dwdx(2)*mass(nj)
     
c         if((x_i.lt.23.*dx.and.x_j.lt.23.*dx).or.(x_i.gt.42.*dx.and
c     &.x_j.gt.42.*dx))then
c      if(sqrt((x_i-x_j)**2+(y_i-y_j)**2).lt.dx.and.nj.le.ntotal)then  !
  !             if(wij.ne.0.0)then
  !    rgylx=(rg(1,i)+rg(1,nj))*(w/wij)**4*massj*cor(1)+(rg(3,i)
  !   &	+rg(3,nj))*(w/wij)**4*massj*cor(2)
	 !rgyly=(rg(2,i)+rg(2,nj))*(w/wij)**4*massj*cor(2)+(rg(4,i)
  !   &	+rg(4,nj))*(w/wij)**4*massj*cor(1)
  !             endif
c          endif
c          endif
c          else 
c           rgylx=0.0
c           rgyly=0.0
c           endif
             eb1=epsino1
                 if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
             !if(sqrt((x_i-x_j)**2+(y_i-y_j)**2).le.2.6*dx1)then         
                     !tj=tj+1
                     !jh(tj,i)=nj
             !write(*,*) p(nj)
                    !if(mod((i-ntotalb),15).eq.0)then
          !p(nj)=-(tsts(1,nj)+tsts(2,nj))/2.
          !pai_ij=0.0
          eb1=0.002!(epsino1+epsino1)/2.
     !!      if(ip_f(nj)==0.and.jc(nj)==0)then
     !! fp1(tj,i)=(-(p_j*rho_nj2+p_i*rho_i2)*cor(1)*massj
     !!& -pai_ij*cor(1)*massj)*massi
     !! fp2(tj,i)=(-(p_j*rho_nj2+p_i*rho_i2)*cor(2)*massj
     !!&-pai_ij*cor(2)*massj)*massi
     !! if(i==5990.and.tj==3)then
     !! write(*,*) "h", fp1(3,5990),nj
     !! endif
     !!      endif
 !     dvx(1,i)=dvx(1,i)-(-(tsts(1,nj)+tsts(2,nj))/2.*rho_nj2+p(i)
 !    &*rho_i2)*cor(1)*massj-pai_ij*cor(1)*massj
 !    & +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
 !    &)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
 !    &*massj  !+rgylx
	!dvx(2,i)=dvx(2,i)-(-(tsts(1,nj)+tsts(2,nj))/2.*rho_nj2+p(i)*rho_i2)
 !    &*cor(2)*massj-pai_ij*cor(2)*massj 
 !    &+((sts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
 !    &+((sts2)*rho_i2+(sts(2,nj))*rho_nj2)*cor(2))
 !    &*massj !+rgyly
             drho(i)=drho(i)+dx*dy*rho(i)*
     &	  ((vx_i-vx_j)*dwdx+(vy_i-vy_j)*dwdy)  
     &+deltf*hsml0*c_ij*(fdeltx*dwdx+fdelty*dwdy)*dx*dx
          
           dvx(1,i)=dvx(1,i)-
     &	  (p(nj)+p(i))*cor(1)*massj/(rho_j*rho01)
     &      -pai_ij*cor(1)*massj
     !!&         +(((sts1)*rho_i2+(ists(1,nj))*rho_nj2
     !!&)*cor(1)+(sts3*rho_i2+ists(3,nj)*rho_nj2)*cor(2))
     !!&*massj  !+rgylx   
	      dvx(2,i)=dvx(2,i)-
     &	  (p(nj)+p(i))*cor(2)*massj/(rho_j*rho01) 
     &      -pai_ij*cor(2)*massj 
     !!&       +((sts4*rho_i2+ists(4,nj)*rho_nj2)*cor(1)
     !!&+((sts2)*rho_i2+(ists(2,nj))*rho_nj2)*cor(2))
     !!&	*massj !+rgyly       
     !!       
            
     !! dyl(1,nj)=dyl(1,nj)+(-(p(nj)+p(i))*dwdx*massj/(rho_j*rho01)
     !!&    )*massi  ! -pai_ij*dwdx*massj
     !! dyl(2,nj)=dyl(2,nj)+(-(p(nj)+p(i))*dwdy*massj/(rho_j*rho01) 
     !!&  )*massi  !-pai_ij*dwdy*massj
     !! 
     !! cdyl(1,nj)=cdyl(1,nj)+(-(p(nj)+p(i))*dwdx*massj/(rho_j*rho01)
     !!&  )*massi  ! -pai_ij*dwdx*massj
     !! cdyl(2,nj)=cdyl(2,nj)+(-(p(nj)+p(i))*dwdy*massj/(rho_j*rho01) 
     !!&   )*massi  !-pai_ij*dwdy*massj
 !     dvx(1,i)=dvx(1,i)-(p_j*rho_nj2+p_i*rho_i2)*cor(1)*massj
 !    &-pai_ij*cor(1)*massj
 !    !! +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
 !    !!&)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
 !    !!&*massj  !+rgylx
	!dvx(2,i)=dvx(2,i)-(p_j*rho_nj2+p_i*rho_i2)*cor(2)*massj
 !    &-pai_ij*cor(2)*massj
 !    !! +((sts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
 !    !!&+((sts2)*rho_i2+(sts(2,nj))*rho_nj2)*cor(2))
 !    !!&	*massj !+rgyly
!c      if(nj.gt.ntotalb)then
      va(1,i)=va(1,i)-2.*eb1*massj/(rho_i+rho_j)
     &	  *(vx_i-vx(1,nj))*w
      va(2,i)=va(2,i)-2.*eb1*massj/(rho_i+rho_j)
     &	  *(vy_i-vx(2,nj))*w 
     !! dyva(1,nj)=dyva(1,nj)-2.*eb1*massj/(rho_i+rho_j)
     !!&	  *(vx_i-vx_j)*w*mass(i)
     !! dyva(2,nj)=dyva(2,nj)-2.*eb1*massj/(rho_i+rho_j)
     !!&	  *(vy_i-vy_j)*w*mass(i)
           !endif
           
           else if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
               eb1=epsino1
          !      p(nj)=0.0
          !pai_ij=0.0   
          !cor(1)=0.0
          !cor(2)=0.0
          !w=0.0
          !          endif
                    !endif
             drho(i)=drho(i)+mass(nj)*
     &	  ((vx_i-vx_j)*dwdx+(vy_i-vy_j)*dwdy)  
     &+deltf*hsml0*c_ij*(fdeltx*dwdx+fdelty*dwdy)*dx1*dx1    
              
                 dvx(1,i)=dvx(1,i)-
     &	  (p(nj)+p(i))*cor(1)*massj/(rho_j*rho01)
     &      -pai_ij*cor(1)*massj
     !!&         +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
     !!&)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
     !!&*massj  !+rgylx   
                 
	      dvx(2,i)=dvx(2,i)-
     &	  (p(nj)+p(i))*cor(2)*massj/(rho_j*rho01) 
     &      -pai_ij*cor(2)*massj  
     !!&       +((sts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
     !!&+((sts2)*rho_i2+(sts(2,nj))*rho_nj2)*cor(2))
     !!&	*massj !+rgyly
                    
 !     dvx(1,i)=dvx(1,i)-(p_j*rho_nj2+p_i*rho_i2)*cor(1)*massj
 !    &-pai_ij*cor(1)*massj
 !    !! +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
 !    !!&)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
 !    !!&*massj  !+rgylx
	!dvx(2,i)=dvx(2,i)-(p_j*rho_nj2+p_i*rho_i2)*cor(2)*massj
 !    &-pai_ij*cor(2)*massj
      !+((sts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
     !!&+((sts2)*rho_i2+(sts(2,nj))*rho_nj2)*cor(2))
     !!&	*massj !+rgyly
!c      if(nj.gt.ntotalb)then
      va(1,i)=va(1,i)-2.*eb1*massj/(rho_i+rho_j)
     &	  *(vx_i-vx_j)*w
      va(2,i)=va(2,i)-2.*eb1*massj/(rho_i+rho_j)
     &	  *(vy_i-vy_j)*w
c      endif
       
c      if(nj.gt.ntotal)then
c      jsd(i)=jsd(i) -(p(nj)+p(i))*cor(2)*mass/(rho(nj)*rho(i))
c     &-pai_ij*cor(2)*mass+((sts(4,i)/rho(i)**2+sts(4,nj)/rho(nj)**2
c     &)*cor(1)+(sts(2,i)/rho(i)**2+sts(2,nj)/rho(nj)**2)*cor(2))*mass !+rgyly	
c      endif
		  endif
           
      enddo
      
      !gs(i)=tj
c        write(55,*) p(11),sts(4,i),sts(2,i),dwdx(2),pai_ij,dvx(2,11) 
c          if(n01.ne.0.0.and.n02.ne.0.0)then  
c          dvx(1,i)=dvx(1,i)/n01
c          dvx(2,i)=dvx(2,i)/n02
c          endif
          !if(i.gt.ntotalb)then
	     dvx(2,i)=dvx(2,i)-9.81
           !endif
	!if(i.le.ntotalb)then
 !         if(mod(i,10).eq.0)then
 !     dvx(2,i)=dvx(2,i)-10./(rho(i)*dx*dx)
 !     endif
 !          endif
                    enddo
      !
          !$omp end do 
      
      !$omp end parallel           
                    
         do i=1,ntotalb
           tj=0
c           if(bj(i)==0)then
           !dvx(1,i)=0.0
           ! dvx(2,i)=0.0
           !
           ! drho(i)=0.0
           ! va(1,i)=0.0  
           !va(2,i)=0.0  
           p_i=p(i)
           rho_i=rho(i)
         sts1=sts(1,i)
       sts2=sts(2,i)
       sts3=sts(3,i)
       sts4=sts(4,i)
       !pp=p(i)
c            n01=0.0
c            n02=0.0
            x_i=x(1,i)
	      y_i=x(2,i)
	      vx_i=vx(1,i)
	      vy_i=vx(2,i)
c	rho(i)=rho(i)/n03(i)
             L_11=0.0
       L_12=0.0
      L_21=0.0
       L_22=0.0
!             do j=1,rn1(i)
! 	      nj=rpn1(j,i)
!            
!	      x_j=x(1,nj)
!            y_j=x(2,nj)                
!       if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
!      hsml0=1/2.*(hsml+hsmls)
!            else
!        hsml0=1/2.*(hsmls+hsmls)
!        endif
!	call kernel_df(x_i,y_i,x_j,y_j,dwdx,dwdy,hsml0)
!  
!      r_1=x_j-x_i
!      r_2=y_j-y_i
!           
!c      r(k)=x(k,j0)
!
!c      dwd(1)=(dwdx-fz(1)/fm)/fm
!c	dwd(2)=(dwdy-fz(2)/fm)/fm 
!     
!      L_11=L_11+r_1*dwdx
!       L_12=L_12+r_1*dwdy
!      L_21=L_21+r_2*dwdx
!       L_22=L_22+r_2*dwdy 
!
!             enddo
!C      L=.i.L
!             JZ=L_11*L_22-L_12*L_21
!            if(jz.ne.0.0)then
!         NL(1,1)=L_22/(JZ*dx*dx)
!       NL(1,2)=-L_12/(JZ*dx*dx)
!      NL(2,1)=-L_21/(JZ*dx*dx)
!       NL(2,2)=L_11/(JZ*dx*dx)
!             
!       else
!              NL(1,1)=0.0
!       NL(1,2)=0.0
!       NL(2,1)=-0.0
!       NL(2,2)=0.0    
!       endif
       
      
       
       
C       write(*,*) JZ
C       call inv(L,dim,Q,is,js)
          rho_i2=1.0/(rho_i*rho_i)     
	   do j=1,rn1(i)
              rgylx=0.0
           rgyly=0.0
           
 	      nj=rpn1(j,i)
          if(nj.gt.ntotalb)then  
            p_j=p(nj)
            ! if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
            !rho_j=rho(nj)
            ! else
              rho_j=rho(nj) 
             !endif
             
              !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
            !fp1(j,i)=0.0
            !fp2(j,i)=0.0
            !endif
          rho_nj2=1.0/(rho_j*rho_j) 
      !calculating the artificial viscosity
	      x_j=x(1,nj)
	      y_j=x(2,nj)
            
            
	      vx_j=vx(1,nj)
	      vy_j=vx(2,nj)
           
             
            h_ij=1/2.*(hsmls+hsmls)
	      rho_ij=0.5*(rho_i+rho_j)
            
	       c_ij=0.5*(c01+c01)
                 
	      v_ij_x_ij=(vx_i-bxvx(1,nj))*(x_i-x_j)+(vy_i-bxvx(2,nj))*(y_i-y_j)
            
      v_ij_x_ij=(bxvx(1,i)-vx_j)*(x_i-x_j)+(bxvx(2,i)-vy_j)*(y_i-y_j)
            
	      phai_ij=h_ij*(v_ij_x_ij)/
     &              ((x_i-x_j)**2+(y_i-y_j)**2+(0.1*h_ij)**2)
            if(v_ij_x_ij.lt.0.0)then
      pai_ij=(-alpha_phai1*c_ij*phai_ij+belda_phai1*phai_ij**2)/rho_ij 
      !pai_ij=(-0.005*c_ij*phai_ij+belda_phai1*phai_ij**2)/rho_ij          
     !!     pai_ij=(-(alpha_phai+alpha_phai1)/2.*c_ij*phai_ij+(belda_phai+
     !!&	  belda_phai1)/2.*phai_ij**2)/rho_ij
            else
	      pai_ij=0.0
            endif
      hsml0=1/2.*(hsmls+hsmls)
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
	     
	    ! dwd(1)=dwdx
	    !dwd(2)=dwdy
C	cor=matmul(NL,dwd) 
C      write(*,*) cor(1),cor(2)
      !!cor(1)=dwdx
      !!cor(2)=dwdy  
      cor(1)=NL1(1,i)*dwdx+NL1(2,i)*dwdy
      cor(2)=NL2(1,i)*dwdx+NL2(2,i)*dwdy  
     !! if(((x_i-x_j)**2+(y_i-y_j)**2).ne.0.0)then
     !! fdeltx=2.*(rho(nj)-rho(i))*(x_j-x_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     !!& -(dpx(i)+dpx(nj))
     !! fdelty=2.*(rho(nj)-rho(i))*(y_j-y_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     !!& -(dpy(i)+dpy(nj))
     !! else
     !!     fdeltx=0.0
     !!     fdelty=0.0
     !! endif
      fdeltx=0.0
          fdelty=0.0
   
      
      
                 
      !cor(1)=NL(1,1)*dwd(1)+NL(1,2)*dwd(2)
      !cor(2)=NL(2,1)*dwd(1)+NL(2,2)*dwd(2)
           
      !cor(1)=NL(1,1)*dwd(1)+NL(1,2)*dwd(2)
      !cor(2)=NL(2,1)*dwd(1)+NL(2,2)*dwd(2)
   
      !if(nj.le.ntotalb)then
      !massj=rho_j*dx*dx
      !massi=rho_i*dx*dx
      
      massj=mass(nj)
      massi=mass(i)
      !else
      !    massj=mass(nj)
      !    endif
C      write(*,*) cor(1),cor(2),DWDX,DWDY
c      n01=n01+dx*dx*(x_j-x_i)*cor(1)
c           n02=n02+dx*dx*(y_j-y_i)*cor(2)
c	    if(mod(i,y_n).eq.0)then
c	      drho(i)=0.
c	    else
	    
c	    endif
c         if(art_vis.eq.1)then
c	dvx(1,i)=dvx(1,i)-((dl(1,i)/rho(i)**2+dl(1,nj)/rho(nj)**2)*dwdx(1)
c     &+(dl(3,i)/rho(i)**2+dl(3,nj)/rho(nj)**2)*dwdx(2))
c     &*mass(nj)-2*pai_ij*dwdx(1)*mass(nj)
c	 vx(2,i)=dvx(2,i)-((dl(4,i)/rho(i)**2+dl(4,nj)/rho(nj)**2)*dwdx(1)
c     &+(dl(2,i)/rho(i)**2+dl(2,nj)/rho(nj)**2)*dwdx(2))
c     &*mass(nj)-2*pai_ij*dwdx(2)*mass(nj)
     
c         if((x_i.lt.23.*dx.and.x_j.lt.23.*dx).or.(x_i.gt.42.*dx.and
c     &.x_j.gt.42.*dx))then
c      if(sqrt((x_i-x_j)**2+(y_i-y_j)**2).lt.dx.and.nj.le.ntotal)then  !
  !             if(wij.ne.0.0)then
  !    rgylx=(rg(1,i)+rg(1,nj))*(w/wij)**4*massj*cor(1)+(rg(3,i)
  !   &	+rg(3,nj))*(w/wij)**4*massj*cor(2)
	 !rgyly=(rg(2,i)+rg(2,nj))*(w/wij)**4*massj*cor(2)+(rg(4,i)
  !   &	+rg(4,nj))*(w/wij)**4*massj*cor(1)
  !             endif
c          endif
c          endif
c          else 
c           rgylx=0.0
c           rgyly=0.0
c           endif
             eb1=epsino1
                 
             !if(sqrt((x_i-x_j)**2+(y_i-y_j)**2).le.2.6*dx1)then         
                     !tj=tj+1
                     !jh(tj,i)=nj
             !write(*,*) p(nj)
                    !if(mod((i-ntotalb),15).eq.0)then
          !p(nj)=-(tsts(1,nj)+tsts(2,nj))/2.
          !pai_ij=0.0
          eb1=0.002!epsino1+epsino1)/2.
     !!      if(ip_f(nj)==0.and.jc(nj)==0)then
     !! fp1(tj,i)=(-(p_j*rho_nj2+p_i*rho_i2)*cor(1)*massj
     !!& -pai_ij*cor(1)*massj)*massi
     !! fp2(tj,i)=(-(p_j*rho_nj2+p_i*rho_i2)*cor(2)*massj
     !!&-pai_ij*cor(2)*massj)*massi
     !! if(i==5990.and.tj==3)then
     !! write(*,*) "h", fp1(3,5990),nj
     !! endif
     !!      endif
 !     dvx(1,i)=dvx(1,i)-(-(tsts(1,nj)+tsts(2,nj))/2.*rho_nj2+p(i)
 !    &*rho_i2)*cor(1)*massj-pai_ij*cor(1)*massj
 !    & +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
 !    &)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
 !    &*massj  !+rgylx
	!dvx(2,i)=dvx(2,i)-(-(tsts(1,nj)+tsts(2,nj))/2.*rho_nj2+p(i)*rho_i2)
 !    &*cor(2)*massj-pai_ij*cor(2)*massj 
 !    &+((sts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
 !    &+((sts2)*rho_i2+(sts(2,nj))*rho_nj2)*cor(2))
 !    &*massj !+rgyly
         !!       
            
      dyl(1,i)=dyl(1,i)+(-(p(nj)+p(i))*cor(1)*massj/(rho_i*rho01)
     &  -pai_ij*cor(1)*massj   )*massi 
      dyl(2,i)=dyl(2,i)+(-(p(nj)+p(i))*cor(2)*massj/(rho_i*rho01) 
     &  -pai_ij*cor(2)*massj )*massi 
      
      cdyl(1,nj)=cdyl(1,nj)+(-(p(nj)+p(i))*cor(1)*massj/(rho_j*rho01)
     & -pai_ij*cor(1)*massj )*massi   
      cdyl(2,nj)=cdyl(2,nj)+(-(p(nj)+p(i))*cor(2)*massj/(rho_j*rho01) 
     & -pai_ij*cor(2)*massj  )*massi  
 !     dvx(1,i)=dvx(1,i)-(p_j*rho_nj2+p_i*rho_i2)*cor(1)*massj
 !    &-pai_ij*cor(1)*massj
 !    !! +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
 !    !!&)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
 !    !!&*massj  !+rgylx
	!dvx(2,i)=dvx(2,i)-(p_j*rho_nj2+p_i*rho_i2)*cor(2)*massj
 !    &-pai_ij*cor(2)*massj
 !    !! +((sts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
 !    !!&+((sts2)*rho_i2+(sts(2,nj))*rho_nj2)*cor(2))
 !    !!&	*massj !+rgyly
!c      if(nj.gt.ntotalb)then
     !! va(1,i)=va(1,i)-2.*eb1*massj/(rho_i+rho_j)
     !!&	  *(vx_i-vx_j)*w
     !! va(2,i)=va(2,i)-2.*eb1*massj/(rho_i+rho_j)
     !!&	  *(vy_i-vy_j)*w 
      dyva(1,i)=dyva(1,i)-2.*eb1*massj/(rho_i+rho_j)
     &	  *(vx_i-vx_j)*w*mass(i)
      dyva(2,i)=dyva(2,i)-2.*eb1*massj/(rho_i+rho_j)
     &	  *(vy_i-vy_j)*w*mass(i)
           !endif
           

		  endif
           
      enddo
      
      !gs(i)=tj
c        write(55,*) p(11),sts(4,i),sts(2,i),dwdx(2),pai_ij,dvx(2,11) 
c          if(n01.ne.0.0.and.n02.ne.0.0)then  
c          dvx(1,i)=dvx(1,i)/n01
c          dvx(2,i)=dvx(2,i)/n02
c          endif
          !if(i.gt.ntotalb)then
	     !dvx(2,i)=dvx(2,i)-9.81
           !endif
	!if(i.le.ntotalb)then
 !         if(mod(i,10).eq.0)then
 !     dvx(2,i)=dvx(2,i)-10./(rho(i)*dx*dx)
 !     endif
 !          endif
                    enddo
	end