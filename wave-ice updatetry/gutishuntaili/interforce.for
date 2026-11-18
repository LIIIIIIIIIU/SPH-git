      subroutine interforce(t,ti)
      implicit none 
	include 'menu.inc'
      double precision ti
       integer nj,m,t,nj1
      double precision x_i,y_i,xi(dim),xj(dim),vx_i,vy_i,gx
     &,x_j,y_j,vx_j,vy_j,h_ij,rho_ij,c_ij,h_ij1,c_ij1,phai_ij1
     &,v_ij_x_ij,phai_ij,hsml0,pai_ij,pai_ij1,p_i,p_j
	double precision w,dwdx,dwdy,JG_J2,sr_i1i2,rho_i,rho_j
      
      double precision sr_i1,sr_i2,sr_i3,sr_i4,tr_i1,tr_i2,tr_i3
     &,tr_i4,ssr_i1,ssr_i2,ssr_i3,ssr_i4,sita,c,s,wij,sg_i1,
     & sg_i2,sg_i3,sg_i4,r_gi1,r_gi2,r_gi3,r_gi4,rgylx,rgyly
     & ,tssr_i1,tssr_i2,tssr_i3,tssr_i4,dis,dis1,r0
      double precision I1,J2,fl,rg1,rg2,rg3,rg4
     & ,pff_i1,pff_i2,pff_i3,pff_i4,psr_i1,psr_i2,psr_i3,psr_i4
     &, tpff_i1,tpff_i2,tpff_i3,tpff_i4
      double precision sts1,sts2,sts3,sts4,pp,bl,wij1
     &,tsts1,tsts2,tsts3,tsts4,ys,rx,rho_i2,rho_nj2
      double precision r_1,r_2,L_11,L_12,L_21,L_22,JZ,
     & dwd(dim),cor(dim)
      
      
       double precision del,def,ljsxybl,ljdxpyb,massj,eb,eb1,massi
       double precision  n00,c11,c12,c22,y1,y12,y2,y22,y3,y32,y4,y42,r1,
     & r12,r2,r22,fdeltx,fdelty,n0,ra,tpff,r01
c	n03=0.0
      
c	ys=3.0e+6 !2.0e+6
     
     
       !  jiao=pi/5.8
       !jzj=0.0!jiao/4.

c        arf=tan(jiao)/(9.+12*tan(jiao)**2)**0.5! ! 2*sin(jiao)/3.**0.5/(3.-sin(jiao))
c      a_rf=2*sin(jzj)/3.**0.5/(3.-sin(jzj))! !sin(jiao)/3.**0.5/(3.-sin(jiao))
c      sgm=3./(9.+12.*tan(jiao)**2)**0.5 ! !sin(jiao)/3.**0.5/(3.-sin(jiao))6*cos(jiao)/3.**0.5/(3.-sin(jiao))

c     &  *6*cos(jiao)/3.**0.5/(3.-sin(jiao))
c      write(*,*)qfqd
    
      k=0
	
    
      
       !$omp parallel private(i,j,k,nj1,sr_i1,sr_i2,sr_i3,sr_i4)!if(ntotalb>2000) 
     & private(sr_i1i2,tr_i1,tr_i2,tr_i3,tr_i4,dis,dis1,r0,r01)
     & private(x_i,y_i,vx_i,vy_i,x_j,y_j,vx_j,vy_j,nj) 
     & private(y1,y12,y2,y22,c11,c12,c22,n00)  
     & private(tsts1,tsts2,tsts3,tsts4)
     & private(w,dwdx,dwdy,hsml0,cor)
     & private(ys,I1,J2,JG_J2,fl,tpff,ra,rx)
     & private(tssr_i1,tssr_i2,tssr_i3,tssr_i4)  
     & private(sts1,sts2,sts3,sts4,massj,rho_j)
     & private(sita,c,s,sg_i1,sg_i2,r_gi1,r_gi2)
     & private(L_11,L_12,L_21,L_22,r_1,r_2,jz)   
      !$omp do schedule(static)
    
       do i=1,ntotalb
c	if(bj(i)==0)then
           sx(i)=0
       sr_i1=0.0
       sr_i2=0.0
       sr_i3=0.0
       sr_i4=0.0
       tr_i1=0.0
       tr_i2=0.0
       tr_i3=0.0
       tr_i4=0.0
        c11=0.0
	c22=0.0
	c12=0.0
	y1=0.0
	y12=0.0
       y2=0.0
	y22=0.0
      
       sts1=sts(1,i)
       sts2=sts(2,i)
       sts3=sts(3,i)
       sts4=sts(4,i)
         tsts1=tsts(1,i)
       tsts2=tsts(2,i)
       tsts3=tsts(3,i)
       tsts4=tsts(4,i)
       
	      x_i=x(1,i)
            y_i=x(2,i)
             vx_i=vx(1,i)
	      vy_i=vx(2,i)
       L_11=0.0
       L_12=0.0
      L_21=0.0
       L_22=0.0
       
      do j=1,rn(i)
 	      nj=rpn(j,i)
             if(nj.le.ntotalb)then
	      x_j=x(1,nj)
            y_j=x(2,nj)  
      !       if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
      !hsml0=1/2.*(hsml+hsmls)
      !      else
        hsml0=1/2.*(hsml+hsml)
        !endif
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
  
      r_1=x_j-x_i
      r_2=y_j-y_i
           
c      r(k)=x(k,j0)

c      dwd(1)=(dwdx-fz(1)/fm)/fm
c	dwd(2)=(dwdy-fz(2)/fm)/fm 
     
      L_11=L_11+r_1*dwdx
       L_12=L_12+r_1*dwdy
      L_21=L_21+r_2*dwdx
       L_22=L_22+r_2*dwdy 
         endif
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
       endif!        
       
C       call inv(L,dim,Q,is,js)
            
            
	    do j=1,rn(i)
 	      nj=rpn(j,i)
         if(nj.le.ntotalb)then 
            
	      x_j=x(1,nj)
            y_j=x(2,nj)
             vx_j=vx(1,nj)
	      vy_j=vx(2,nj)
           
      hsml0=hsml
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
 !     dwd(1)=dwdx
	!dwd(2)=dwdy
C	cor=matmul(NL,dwd) 
C      write(*,*) cor(1),cor(2)
      !!cor(1)=dwdx
      !!cor(2)=dwdy  
       cor(1)=NL1(1,i)*dwdx+NL1(2,i)*dwdy
      cor(2)=NL2(1,i)*dwdx+NL2(2,i)*dwdy 
      
        !if(sx(i)==3.and.sx(nj)==4.and.sh(i)==1.and.sh(nj)==1)then
        !  if((x_j-x_i).gt.1.5*dx)then
        !    cor(1)=0.0
        !    cor(2)=0.0
        !    w=0.0
        !  endif
        !    endif
        !if(sx(i)==4.and.sx(nj)==3.and.sh(i)==1.and.sh(nj)==1)then
        !  if((x_j-x_i).lt.-1.5*dx)then
        !    cor(1)=0.0
        !    cor(2)=0.0
        !    w=0.0353
        !  endif
        !  endif 
       IF(NJ.NE.I.and.i.gt.50.and.nj.gt.50.and.ljsxyb(i).gt.0.0.and
     &     .ljsxyb(nj).gt.0.0)THEN
         dis=sqrt((x(1,i)-x(1,nj))**2
     &+(x(2,i)-x(2,nj))**2)
      dis1=sqrt((xc(1,i)-xc(1,nj))**2+(xc(2,i)-xc(2,nj))**2)
    
         if(dis1.ge.1.2*dx)then
          r0=0.45*dx/dis1
          r01=0.7*dx/dis1
      else
           r0=0.5*dx/dis1
           r01=0.75*dx/dis1
      endif
      if((dis-dis1)/dis1.ge.(r01))then
             
              if(x(1,i).lt.x(1,nj))then     
      sx(i)=3
      sx(nj)=4
      else
        sx(i)=4
      sx(nj)=3
      endif
          cor(1)=0.0
            cor(2)=0.0
            w=0.0
       endif
     
      
      
!!      if((dis-dis1)/dis1.ge.(r0))then
!!              do k=1,rn(i)
!!c           k=k+1
!! 	      nj1=rpn(k,i)
!!            if(nj1.gt.ntotalb.and.nj1.le.ntotal)then  
!!          !x_j=x(1,nj1)
!!          !  y_j=x(2,nj1)
!!            !z_j=x(3,nj1)                    
!!      if(sqrt(((x(1,nj)+x(1,i))/2.-x(1,nj1))**2+((x(2,nj)+x(2,i))/2.
!!     & -x(2,nj1))**2).le.3.0*dx)then     
!!              if(x(1,i).lt.x(1,nj))then     
!!      sx(i)=1
!!      sx(nj)=2
!!      else
!!        sx(i)=2
!!      sx(nj)=1
!!      endif
!!          cor(1)=0.0
!!            cor(2)=0.0
!!            w=0.0
!!       endif
!!      endif
!!      enddo
!!              endif 
        ENDIF      
      
 !     y1=y1+(vx_j-vx_i)*(x_j-x_i)*w    
 !      y12=y12+(vx_j-vx_i)*(y_j-y_i)*w
 !      y2=y2+(vy_j-vy_i)*(x_j-x_i)*w
 !      y22=y22+(vy_j-vy_i)*(y_j-y_i)*w
 !                
	!c11=c11+(x_j-x_i)**2*w
 !     c22=c22+(y_j-y_i)**2*w
 !
	!c12=c12+(x_j-x_i)*(y_j-y_i)*w  
      sr_i1=sr_i1+(vx_j-vx_i)*cor(1) !xx
      sr_i2=sr_i2+(vy_j-vy_i)*cor(2) !yy
      sr_i3=sr_i3+0.5*((vx_j-vx_i)*cor(2)+(vy_j-vy_i)*cor(1))               !xy
	sr_i4=sr_i3 !sr(4,i)+0.5*mass(nj)/rho(nj)*((vx(1,nj)-vx(1,i))*
        ! &		  dwdx(2)+(vx(2,nj)-vx(2,i))*dwdx(1))               !yx 
      tr_i1=0.0 !xx
      tr_i2=0.0 !yy
      tr_i3=tr_i3+0.5*((vx_j-vx_i)*cor(2)-(vy_j-vy_i)*cor(1))               !xy
      tr_i4=-tr_i3 !tr(4,i)+0.5*mass(nj)/rho(nj)*((vx(2,nj)-vx(2,i))*   
          !&		  dwdx(1)-(vx(1,nj)-vx(1,i))*dwdx(2))  !sr(2,i)     !yx 
c	n03(i)=n03(i)+mass/rho(nj)*w
      endif
          enddo
 !!       n00=c11*c22-c12*c12
 !!     if(n00.ne.0.0)then
 !!     sr_i1=((c22*y1-c12*y12))/n00
 !!     sr_i2=((c11*y22-c12*y2))/n00            
 !!     sr_i3=0.5*((c11*y12-c12*y1)+(c22*y2-c12*y22))/n00
 !!     sr_i4=sr_i3 
 !!     tr_i1=0.0
 !!     tr_i2=0.0           
 !!     tr_i3=0.5*((c11*y12-c12*y1)-(c22*y2-c12*y22))/n00
 !!     tr_i4=-tr_i3 
	!!else
 !!     sr_i1=0.0 !-(c22*y1-c12*y2)/n00/rho0
 !!     sr_i2=0.0 !-(c11*y2-c12*y1)/n00/rho0
 !!     sr_i3=0.0 !-(c22*y1-c12*y2)/n00/rho0
 !!     sr_i4=0.0
 !!      tr_i1=0.0
 !!     tr_i2=0.0           
 !!     tr_i3=0.0
 !!     tr_i4=0.0
	!!endif   
      
          
          
          ybl(1,i)=sr_i1
          ybl(2,i)=sr_i2
          ybl(3,i)=sr_i3
          ybl(4,i)=sr_i4
          
      sr_i1=sr_i1*dx*dx
      sr_i2=sr_i2*dx*dx
      sr_i3=sr_i3*dx*dx             
	sr_i4=sr_i3
      tr_i1=0.0 !xx
      tr_i2=0.0 !yy
      tr_i3=tr_i3*dx*dx                !xy
      tr_i4=-tr_i3 
          
      yb(1,i)=yb(1,i)+sr_i1*dt !xx
      yb(2,i)=yb(2,i)+sr_i2*dt !yy
	yb(3,i)=yb(3,i)+sr_i3*dt !xy
	yb(4,i)=yb(4,i)+sr_i4*dt !yx
      ! I1=(tsts1+tsts2)
      ! J2=(0.5*(sts1**2+sts2**2+sts3**2+sts4**2))
      !&**0.5
      
      ys=qfqd(i)*sgm !2.0e+6   ! 6.*cos(jiao)/3.**0.5/(3.-sin(jiao)) 3./(9.+12.*tan(jiao)**2)**0.5
        I1=(tsts1+tsts2)
      J2=(0.5*(sts1**2+sts2**2+sts3**2+sts4**2))
     &**0.5
      sr_i1i2=sr_i1+sr_i2
      if(J2.ge.(ys-arf*I1/2.))then
          JG_J2=JG/J2
      cs(i)=1
       fl=sts1*sr_i1+sts2*sr_i2+sts3*sr_i3+sts4*sr_i4
      tpff=((arf*KG*sr_i1i2+JG_J2*fl
     &)/(KG*arf*a_rf+JG))
      if(tpff.le.0.0)then
          tpff=0.0
          endif
c      tpff_i2=(arf*KG*(sr_i1+sr_i2)+(JG/J2)*fl
c     &)/(KG*arf*a_rf+JG)
c      tpff_i3=(arf*KG*(sr_i1+sr_i2)+(JG/J2)*fl
c     &)/(KG*arf*a_rf+JG)
c      tpff_i4=(arf*KG*(sr_i1+sr_i2)+(JG/J2)*fl
c     &)/(KG*arf*a_rf+JG)
!      ssr_i1=2.*JG*(sr_i1-1./3.*(sr_i1+sr_i2))-tpff*(JG_J2*sts1)!+sts1*tr_i1+sts3*tr_i3+
!c     & sts1*tr_i1+sts4*tr_i3  !+   !xx
!      ssr_i2=2.*JG*(sr_i2-1./3.*(sr_i1+sr_i2))-tpff*(JG_J2*sts2)!+sts4*tr_i4+sts2*tr_i2+
!c     &sts3*tr_i4+sts2*tr_i2!yy
!	ssr_i3=2.*JG*sr_i3-tpff*(JG_J2*sts3)!+sts1*tr_i4+sts3*tr_i2+sts3*tr_i1+sts2
!c     &*tr_i3!xy
!	ssr_i4=2.*JG*sr_i4-tpff*(JG_J2*sts4)!+sts4*tr_i1+sts2*tr_i3+sts1*tr_i4+sts4
         tssr_i1=2.*JG*(sr_i1-1./3.*sr_i1i2)+KG*(
     &sr_i1i2)-tpff*(KG*a_rf+JG_J2*sts1)+tsts1*tr_i1+tsts3*tr_i3+
     &tsts1*tr_i1+tsts4*tr_i3!xx   !xx
      tssr_i2=2.*JG*(sr_i2-1./3.*sr_i1i2)+KG*(
     & sr_i1i2)-tpff*(KG*a_rf+JG_J2*sts2)+tsts4*tr_i4+tsts2*tr_i2+
     &tsts3*tr_i4+tsts2*tr_i2!yy !yy
	tssr_i3=2.*JG*sr_i3-tpff*(JG_J2*sts3)+tsts1*tr_i4+tsts3*tr_i2+
     &tsts3*tr_i1+tsts2*tr_i3!xy
	tssr_i4=2.*JG*sr_i4-tpff*(JG_J2*sts4)+tsts4*tr_i1+tsts2*tr_i3+
     & tsts1*tr_i4+tsts4*tr_i2
      
      
      ljsxyb(i)=ljsxyb(i)+tpff*sgm*dt
      
      !if(i==421)then
      !write(*,*) ljsxyb(421)
      !endif
      else
        tpff=0.0
c      tpff_i2=0.0
c      tpff_i3=0.0
c      tpff_i4=0.0 
c        write(*,*) "h"  
         tssr_i1=2.*JG*(sr_i1-1./3.*sr_i1i2)+KG*(
     &sr_i1i2)+tsts1*tr_i1+tsts3*tr_i3+
     &tsts1*tr_i1+tsts4*tr_i3!xx   !xx
      tssr_i2=2.*JG*(sr_i2-1./3.*sr_i1i2)+KG*(
     & sr_i1i2)+tsts4*tr_i4+tsts2*tr_i2+
     &tsts3*tr_i4+tsts2*tr_i2!yy !yy
	tssr_i3=2.*JG*sr_i3+tsts1*tr_i4+tsts3*tr_i2+
     &tsts3*tr_i1+tsts2*tr_i3!xy
	tssr_i4=2.*JG*sr_i4+tsts4*tr_i1+tsts2*tr_i3+
     & tsts1*tr_i4+tsts4*tr_i2
!        ssr_i1=2.*JG*(sr_i1-1./3.*(sr_i1+sr_i2))!+sts1*tr_i1+sts3*tr_i3+
!c     & sts1*tr_i1+sts4*tr_i3  !+   !xx
!      ssr_i2=2.*JG*(sr_i2-1./3.*(sr_i1+sr_i2))!+sts4*tr_i4+sts2*tr_i2+
!c     &sts3*tr_i4+sts2*tr_i2!yy
!	ssr_i3=2.*JG*sr_i3!+sts1*tr_i4+sts3*tr_i2+sts3*tr_i1+sts2
!c     &*tr_i3!xy
!	ssr_i4=2.*JG*sr_i4!+sts4*tr_i1+sts2*tr_i3+sts1*tr_i4+sts4
      endif
c      write(*,*) tpff
!      ssr_i1=2.*JG*(sr_i1-1./3.*(sr_i1+sr_i2))!+sts1*tr_i1+sts3*tr_i3+
!c     & sts1*tr_i1+sts4*tr_i3  !+   !xx
!      ssr_i2=2.*JG*(sr_i2-1./3.*(sr_i1+sr_i2))!+sts4*tr_i4+sts2*tr_i2+
!c     &sts3*tr_i4+sts2*tr_i2!yy
!	ssr_i3=2.*JG*sr_i3!+sts1*tr_i4+sts3*tr_i2+sts3*tr_i1+sts2
!c     &*tr_i3!xy
!	ssr_i4=2.*JG*sr_i4!+sts4*tr_i1+sts2*tr_i3+sts1*tr_i4+sts4
!c     &*tr_i2!yx
    
c     &*tr_i2!yx
    
c       ssr(1,i)=ssr_i1
c        ssr(2,i)=ssr_i2
c        ssr(3,i)=ssr_i3
c        ssr(4,i)=ssr_i4
c        tssr(1,i)=tssr_i1
c        tssr(2,i)=tssr_i2
c        tssr(3,i)=tssr_i3
c        tssr(4,i)=tssr_i4
 !     sts(1,i)=sts(1,i)+ssr_i1*dt !xx
 !     sts(2,i)=sts(2,i)+ssr_i2*dt !yy
	!sts(3,i)=sts(3,i)+ssr_i3*dt !xy
	!sts(4,i)=sts(4,i)+ssr_i4*dt !yx
      tsts(1,i)=tsts(1,i)+tssr_i1*dt !xx
      tsts(2,i)=tsts(2,i)+tssr_i2*dt !yy
	tsts(3,i)=tsts(3,i)+tssr_i3*dt !xy
	tsts(4,i)=tsts(4,i)+tssr_i4*dt !yx
       sts(1,i)=tsts(1,i)-(tsts(1,i)+tsts(2,i))/2. !xx
      sts(2,i)=tsts(2,i)-(tsts(1,i)+tsts(2,i))/2. !yy
	sts(3,i)=tsts(3,i) !xy
	sts(4,i)=tsts(4,i) !yx
c      if(mod(ti,10000*dt).eq.0)then
c      write(*,*) ljsxybl
c      endif      
c      if(cs(i)==1)then
c      pyb(1,i)=pyb(1,i)+psr_i1*dt
c      pyb(2,i)=pyb(2,i)+psr_i2*dt
c      pyb(3,i)=pyb(3,i)+psr_i3*dt
c      pyb(4,i)=pyb(4,i)+psr_i4*dt
c      ljsxyb(i)=ljsxyb(i)+(ljsxybl)*dt
c      elseif(cs(i)==0)then
c         pyb(1,i)=pyb(1,i)+0.0
c      pyb(2,i)=pyb(2,i)+0.0
c      pyb(3,i)=pyb(3,i)+0.0
c      pyb(4,i)=pyb(4,i)+0.0
c      endif
c      if(cs(i)==1)then
c        dxpyb=dxpyb+(pyb(1,i)**2+pyb(2,i)**2+pyb(3,i)**2+pyb(4,i)  !+(ljsxyb(i)
c     &**2)**0.5
c      endif
c      k=k+cs(i)
      
      !!!!1 stress correction and cohesion softening
       sts1=sts(1,i)
       sts2=sts(2,i)
       sts3=sts(3,i)
       sts4=sts(4,i)  
       tsts1=tsts(1,i)
       tsts2=tsts(2,i)
       tsts3=tsts(3,i)
       tsts4=tsts(4,i)
c      dxyl(i)=(1.5*(sts1**2+sts2**2+sts3**2+sts4**2))**0.5    
c      I1=(sts1+sts2-2.*p(i))
      I1=(tsts1+tsts2)
      J2=(0.5*(sts1**2+sts2**2+sts3**2+sts4**2))
     &**0.5
      
       !IF(i.gt.50)THEN
          if(qfqd(i).gt.1./20.*q0) then
c           qfqd=q0-E*dx**1.1/1.75*ljdxpyb
         rx=1.2*E/10.
           qfqd(i)=q0-rx*ljsxyb(i)
            zjqd(i)=qfqd(i) 
c            zjpyb=ljdxpyb
            zjpyb(i)=ljsxyb(i)
      endif
!      if(qfqd(i).gt.1./50.*q0.and.qfqd(i).lt.1./3.*q0)
!     & then
!c      write(*,*) zjqd,1./3.*q0
!      rx=1.2*E/30.
!C            qfqd=zjqd-E*dx**1.1/1.75*1./3*(ljdxpyb-zjpyb)
!       qfqd(i)=zjqd(i)-rx*(ljsxyb(i)-zjpyb(i))   
!      endif
       if(qfqd(i).le.1./20.*q0) then
           rx=0.0
            qfqd(i)=qfqd(i)
        sh(i)=1
      endif 
       if(qfqd(i).le.0.0) then
           rx=0.0
            qfqd(i)=0.0
        sh(i)=2
      endif 
      
     
       ys=qfqd(i)*sgm
     
      if(J2.gt.(ys-arf*I1/2.))then
           if((ys-arf*I1/2.).lt.0.0)then
          tsts(3,i)=tsts3   
      tsts(4,i)=tsts4
       tsts(1,i)=tsts1-1./2.*(I1-2.*ys/arf) ! ra*tsts(1,i) 
      tsts(2,i)=tsts2-1./2.*(I1-2.*ys/arf) ! ra*tsts(2,i) 
c           endif
         else
           
 !         sts(1,i)=tsts(1,i)-(tsts(1,i)+tsts(2,i))/2. !xx
 !     sts(2,i)=tsts(2,i)-(tsts(1,i)+tsts(2,i))/2. !yy
	!sts(3,i)=tsts(3,i) !xy
	!sts(4,i)=tsts(4,i) !yx 
 !          sts1=sts(1,i)
 !      sts2=sts(2,i)
 !      sts3=sts(3,i)
 !      sts4=sts(4,i)  
 !
 !          I1=(tsts(1,i)+tsts(2,i))
 !     J2=(0.5*(sts1**2+sts2**2+sts3**2+sts4**2))
 !    &**0.5
 !         if(J2.gt.(ys-arf*I1/3.))then
c           cs(i)=cs(i)+1
              if(J2.ne.0.0)then
      ra=(ys-arf*I1/2.)/J2
      sts(1,i)=ra*sts1
      sts(2,i)=ra*sts2  
      sts(3,i)=ra*sts3
      sts(4,i)=ra*sts4
      tsts(3,i)=ra*tsts3    
      tsts(4,i)=ra*tsts4
       tsts(1,i)=sts(1,i)+1./2.*I1!  ra*tsts(1,i)  !
      tsts(2,i)=sts(2,i)+1./2.*I1  ! ra*tsts(2,i) !
c       if(mod(ti,10000*dt).eq.0)then
c       write(*,*)sts(1,i)+(tsts(1,i)+tsts(2,i))/2.,tsts(1,i),sts(2,i)+
c     &(tsts(1,i)+tsts(2,i))/2.,tsts(2,i)
      endif
      endif
      endif
      
      !!Endif
      
      
      !!!! artifical stress
       tsts1=tsts(1,i)
       tsts2=tsts(2,i)
       tsts3=tsts(3,i)
       tsts4=tsts(4,i)
       !pp=p(i)
	if((tsts1-tsts2).gt.0.0)then    
	sita=(atan(2.*tsts3/(tsts1-tsts2
     &  ))/2.) !+0.000001*rho0*(co/2.)**2
      else if((tsts1-tsts2).lt.0.0)then 
      sita=(atan(2.*tsts3/(tsts1-tsts2
     &  ))/2.) !-0.000001*rho0*(co/2.)**2
      else if ((tsts1-tsts2).eq.0.0)then
       sita=pi/4.
      endif
c      write(*,*)sita
	c=cos(sita)
	s=sin(sita)
	sg_i1=c*c*(tsts1)+2*s*c*tsts3+s*s*(tsts2)
      sg_i2=s*s*(tsts1)+2*s*c*tsts3+c*c*(tsts2) 
      if(sg_i1.gt.0.0)then
      r_gi1=-0.3*sg_i1/rho(i)**2
	else
	r_gi1=0.0
	endif
      if(sg_i2.gt.0.0)then
      r_gi2=-0.3*sg_i2/rho(i)**2
	else
	r_gi2=0.0
	endif
      rg(1,i)=c*c*r_gi1+s*s*r_gi2
      rg(2,i)=s*s*r_gi1+c*c*r_gi2
      rg(3,i)=c*s*(r_gi1-r_gi2)
	rg(4,i)=-rg(3,i)
c      if((tsts(i,1)).gt.0.0)then
c      rg(1,i)=-0.3*(tsts(1,i))/rho(i)**2
c	else
c	rg(1,i)=0.0
c	endif
c      if((tsts(i,2)).gt.0.0)then
c      rg(2,i)=-0.3*(tsts(2,i))/rho(i)**2
c	else
c	rg(2,i)=0.0
c	endif
c      if(sts(i,3).gt.0.0)then
c      rg(3,i)=-0.3*tsts(3,i)/rho(i)**2
c	else
c	rg(3,i)=0.0
c      endif
c      if(tsts(i,4).gt.0.0)then
c      rg(4,i)=-0.3*tsts(4,i)/rho(i)**2
c	else
c	rg(4,i)=0.0
c	endif
      
      enddo
        !$omp end do 
      
      
      !$omp end parallel 
      !if(mod(t,200).eq.0)then
      ! write(*,*) tsts(1,1),tsts(1,2)
      ! endif
        
!       do i=1,ntotalb 
!                
!       k=0
!      pp=0.0
!      !dpxx=0.0
!      !dpyy=0.0
!!      tsts1=0.0
!!      tsts2=0.0
!!      tsts3=0.0
!!      tsts4=0.0
!      n0=0.0        
!               x_i=x(1,i)
!	      y_i=x(2,i)
!        do j=1,rn(i)
! 	      nj=rpn(j,i)
!            !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
!            if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
!	      x_j=x(1,nj)
!            y_j=x(2,nj)                
!       !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
!      hsml0=1/2.*(hsml+hsml)
!      !else
!      !  hsml0=1/2.*(hsmls+hsml)
!      !  endif
!	call kernel_f(x_i,y_i,x_j,y_j,w,hsml0)
!  
!      
!      !if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
!      ! p(nj)=-(tsts(1,nj)+tsts(2,nj))/2. 
!      ! endif
!               k=k+1    
!      pp=pp+p(nj)*w+(g+dvx(2,i))*rho(nj)*(x(2,nj)-x(2,i))*w         !+dvx(2,i)+(dvx(2,i)-dvx(2,nj))
!     &+rho(nj)*(x(1,nj)-x(1,i))*(dvx(1,i))*w
!      !dpxx=dpxx+dpx(nj)*w
!      !dpyy=dpyy+dpy(nj)*w
!      
!             !pp=pp+p(nj)
!!             tsts1=tsts1+tsts(1,nj)*w
!!            tsts2=tsts2+tsts(2,nj)*w
!!            tsts3=tsts3+tsts(3,nj)*w
!!            tsts4=tsts4+tsts(4,nj)*w
!            n0=n0+w
!            endif
!        enddo
!       
!        !p(i)=-(tsts(1,i)+tsts(2,i))/3.
!       !if(i==2)then
!       !    p(2)=pp/n0
!       !           write(*,*) p(2),pp,k
!       !           
!       !           endif
!	  if(n0.ne.0)then
!            !write(*,*) "h1"
!	p(i)=pp/n0
!      !dpx(i)=dpxx/n0
!      !dpy(i)=dpyy/n0
!      !p(i)=pp/k
!       !if(i==1)then
!       !           write(*,*) p(1),pp,k
!       !           endif
!!      tsts(1,ntotal+nvirt1+nvirt2+nvirt3)=tsts1/n0
!!      tsts(2,ntotal+nvirt1+nvirt2+nvirt3)=tsts2/n0
!!      tsts(3,ntotal+nvirt1+nvirt2+nvirt3)=tsts3/n0
!!      tsts(4,ntotal+nvirt1+nvirt2+nvirt3)=tsts4/n0
!      endif
!!       
!      enddo
      

      

      
      
 !     do i=ntotalb+1,ntotal
 !       k=0
 !     pp=0.0
 !     tsts1=0.0
 !     tsts2=0.0
 !     tsts3=0.0
 !     tsts4=0.0
 !     n0=0.0 
 !            x_i=x(1,i)
 !           y_i=x(2,i)  
 !           do j=1,rn(i)
 !	      nj=rpn(j,i)
 !            !if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
 !           if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
	!      x_j=x(1,nj)
 !           y_j=x(2,nj)                
 !     ! if(nj.le.ntotalb.or.nj.gt.(ntotal+nvirt1+nvirt2))then
 !     !hsml0=1/2.*(hsml+hsmls)
 !     !      else
 !       hsml0=1/2.*(hsml+hsml)
 !       !endif
	!call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
 !        !     if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
 !        ! tsts(1,nj)=-p(nj)+sts(1,nj)
 !        !tsts(2,nj)=-p(nj)+sts(2,nj)
 !        !tsts(3,nj)=sts(3,nj)
 !        !tsts(4,nj)=sts(4,nj) 
 !        !        endif
 !                
 !           pp=pp+rho(nj)*w
 !                k=k+1
 !           ! tsts1=tsts1+tsts(1,nj)*w
 !           !tsts2=tsts2+tsts(2,nj)*w
 !           !tsts3=tsts3+tsts(3,nj)*w
 !           !tsts4=tsts4+tsts(4,nj)*w
 !           !tsts1=tsts1+tsts(1,nj)
 !           !tsts2=tsts2+tsts(2,nj)
 !           !tsts3=tsts3+tsts(3,nj)
 !           !tsts4=tsts4+tsts(4,nj)
 !           n0=n0+w
 !           endif
 !           !write(*,*) k,rn1(i)
 !           enddo
 !        !   tsts(1,i)=-p(i)+sts(1,i)
 !        !tsts(2,i)=-p(i)+sts(2,i)
 !        !tsts(3,i)=sts(3,i)
 !        !tsts(4,i)=sts(4,i)
	!  if(n0.ne.0)then
	!!p(nj)=pp/n0
 !     !ztsts(1,i)=tsts1/n0
 !     !ztsts(2,i)=tsts2/n0
 !     !ztsts(3,i)=tsts3/n0
 !     !ztsts(4,i)=tsts4/n0
 !     !tsts(1,i)=tsts1/k
 !     !tsts(2,i)=tsts2/k
 !     !tsts(3,i)=tsts3/k
 !     !tsts(4,i)=tsts4/k
 !     ! tsts(1,i)=tsts1/k
 !     !tsts(2,i)=tsts2/k
 !     !tsts(3,i)=tsts3/k
 !     !tsts(4,i)=tsts4/k
 !     !
 !            !if(i==5250)then
 !            !     write(*,*) ztsts(1,5250),tsts1,k
 !                 !endif
 !           xrho(i)=pp/n0
 !     endif
 !      
 !     enddo
        
      !do i=ntotalb+1,ntotal
      !tsts(1,i)=ztsts(1,i)
      !tsts(2,i)=ztsts(2,i)
      !tsts(3,i)=ztsts(3,i)
      !tsts(4,i)=ztsts(4,i)
      !enddo
       !write(*,*) ztsts(1,5250),"h"
       
	 !x_i=0.0
  !     y_i=0.0
  !      x_j=0.0
  !     y_j=0.0
  !     hsml0=hsml
  !    call kernel_f(x_i,y_i,x_j,y_j,w,hsml0)
  !    wij1=w 
      
      x_i=0.0
       y_i=0.0
      
       x_j=0.0+dx
       y_j=0.0
       hsml0=hsml
      call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
      wij=w
c      write(*,*) wij1/wij
           
      !$omp parallel private(i,j,k,nj1,sts1,sts2,sts3,sts4)!if(ntotalb>2000) 
     & private(p_i,p_j,rho_i,rho_j,rho_i2,rg1,rg2,rg3,rg4,dis,dis1,r0)
     & private(x_i,y_i,vx_i,vy_i,massi,v_ij_x_ij,phai_ij,pai_ij,c_ij) 
     & private(x_j,y_j,vx_j,vy_j,nj,rho_nj2,massj,h_ij,rho_ij,r01)  
     & private(tsts1,tsts2,tsts3,tsts4,eb)
     & private(w,dwdx,dwdy,hsml0,cor,rgylx,rgyly)
     & private(L_11,L_12,L_21,L_22,r_1,r_2,jz)   
     !!& private(tssr_i1,tssr_i2,tssr_i3,tssr_i4,tssr_i5,tssr_i6)  
     !!& private(tssr_i7,tssr_i8,tssr_i9)
      
      !$omp do schedule(dynamic,50)   
      
	 do i=1,ntotalb
c           if(bj(i)==0)then
         sx(i)=0  
           dvx(1,i)=0.0
            dvx(2,i)=0.0
         ddvx(1,i)=0.0
            ddvx(2,i)=0.0
            drho(i)=0.0
            va(1,i)=0.0  
           va(2,i)=0.0 
           p_i=p(i)
           rho_i=rho(i)
         tsts1=tsts(1,i)
       tsts2=tsts(2,i)
       tsts3=tsts(3,i)
       tsts4=tsts(4,i)
        sts1=sts(1,i)
       sts2=sts(2,i)
       sts3=sts(3,i)
       sts4=sts(4,i)
      rg1=rg(1,i)
       rg2=rg(2,i)
       rg3=rg(3,i)
       rg4=rg(4,i)
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
       
      do j=1,rn(i)
 	      nj=rpn(j,i)
             !!if(nj.le.ntotalb)then
	      x_j=x(1,nj)
            y_j=x(2,nj)  
      !       if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
      !hsml0=1/2.*(hsml+hsmls)
      !      else
        hsml0=1/2.*(hsml+hsml)
        !endif
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
  
      r_1=x_j-x_i
      r_2=y_j-y_i
           
c      r(k)=x(k,j0)

c      dwd(1)=(dwdx-fz(1)/fm)/fm
c	dwd(2)=(dwdy-fz(2)/fm)/fm 
     
      L_11=L_11+r_1*dwdx
       L_12=L_12+r_1*dwdy
      L_21=L_21+r_2*dwdx
       L_22=L_22+r_2*dwdy 
         !!endif
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
       endif!        
       
C       write(*,*) JZ
C       call inv(L,dim,Q,is,js)
          rho_i2=1.0/(rho_i*rho_i)     
	   do j=1,rn(i)
              rgylx=0.0
           rgyly=0.0
 	      nj=rpn(j,i)
            p_j=p(nj)
            rho_j=rho(nj)
          rho_nj2=1.0/(rho_j*rho_j) 
      !calculating the artificial viscosity
	      x_j=x(1,nj)
	      y_j=x(2,nj)
            
            
	      vx_j=vx(1,nj)
	      vy_j=vx(2,nj)
             massj=mass(nj)
          massi=mass(i)
            
           if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
             h_ij=1/2.*(hsml+hsml)
              !h_ij1=1/2.*(hsmls+hsmls)
	      rho_ij=0.5*(rho_i+rho_j)
            
	      c_ij=0.5*(co+co)
            !c_ij1=0.5*(c01+c01) 
            
	      v_ij_x_ij=(vx_i-vx_j)*(x_i-x_j)+(vy_i-vy_j)*(y_i-y_j)
            
	      phai_ij=h_ij*(v_ij_x_ij)/
     &              ((x_i-x_j)**2+(y_i-y_j)**2+(0.1*h_ij)**2)
     !!       phai_ij1=h_ij1*(v_ij_x_ij)/
     !!&              ((x_i-x_j)**2+(y_i-y_j)**2+(0.1*h_ij1)**2)
            if(v_ij_x_ij.lt.0.0)then
      !pai_ij=(-alpha_phai1*c_ij*phai_ij+belda_phai1*phai_ij**2)/rho_ij 
      pai_ij=(-0.05*c_ij*phai_ij+belda_phai1*phai_ij**2)/rho_ij 
     !! pai_ij1=(-alpha_phai1*c_ij1*phai_ij1+belda_phai1*phai_ij**2)
     !!&/rho_ij
		  !pai_ij=(-(alpha_phai+alpha_phai1)/2.*c_ij*phai_ij+(belda_phai+
    ! &	  belda_phai1)/2.*phai_ij**2)/rho_ij
            else
	      pai_ij=0.0
            !pai_ij1=0.0
            endif
            
             hsml0=1/2.*(hsml+hsml)
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
	     
	    ! dwd(1)=dwdx
	    !dwd(2)=dwdy
C	cor=matmul(NL,dwd) 
C      write(*,*) cor(1),cor(2)
      cor(1)=dwdx
      cor(2)=dwdy  
      !cor(1)=NL(1,1)*dwd(1)+NL(1,2)*dwd(2)
      !cor(2)=NL(2,1)*dwd(1)+NL(2,2)*dwd(2)   
     !!      if(((x_i-x_j)**2+(y_i-y_j)**2).ne.0.0)then
     !! fdeltx=2.*(rho(nj)-rho(i))*(x_j-x_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     !!& -(dpx(i)+dpx(nj))
     !! fdelty=2.*(rho(nj)-rho(i))*(y_j-y_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     !!& -(dpy(i)+dpy(nj))
     !! else
     !!     fdeltx=0.0
     !!     fdelty=0.0
     !! endif     
      else
         h_ij=hsml
	      rho_ij=0.5*(rho_i+rho_j)
            
	      c_ij=0.5*(co+co)
                 
	      v_ij_x_ij=(vx_i-vx_j)*(x_i-x_j)+(vy_i-vy_j)*(y_i-y_j)
            
	      phai_ij=h_ij*(v_ij_x_ij)/
     &              ((x_i-x_j)**2+(y_i-y_j)**2+(0.1*h_ij)**2)
            if(v_ij_x_ij.lt.0.0)then
		  pai_ij=(-alpha_phai*c_ij*phai_ij+belda_phai*phai_ij**2)
     &	  /rho_ij
            else
	      pai_ij=0.0
            endif
      hsml0=1/2.*(hsml+hsml)
	call kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
	     
	    ! dwd(1)=dwdx
	    !dwd(2)=dwdy
C	cor=matmul(NL,dwd) 
C      write(*,*) cor(1),cor(2)
      !!cor(1)=dwdx
      !!cor(2)=dwdy  
       cor(1)=NL1(1,i)*dwdx+NL1(2,i)*dwdy
      cor(2)=NL2(1,i)*dwdx+NL2(2,i)*dwdy 
      
     !!  if(((x_i-x_j)**2+(y_i-y_j)**2).ne.0.0)then
     !! fdeltx=2.*(rho(nj)-rho(i))*(x_j-x_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     !!& -(dpx(i)+dpx(nj))
     !! fdelty=2.*(rho(nj)-rho(i))*(y_j-y_i)/((x_i-x_j)**2+(y_i-y_j)**2)
     !!& -(dpy(i)+dpy(nj))
     !! else
     !!     fdeltx=0.0
     !!     fdelty=0.0
     !! endif
       if(wij.ne.0.0)then
      rgylx=(rg1+rg(1,nj))*(w/wij)**yt*massj*cor(1)+(rg3
     &	+rg(3,nj))*(w/wij)**yt*massj*cor(2)
      rgyly=(rg2+rg(2,nj))*(w/wij)**yt*massj*cor(2)+(rg4
     &	+rg(4,nj))*(w/wij)**yt*massj*cor(1)
               endif
      endif
          

      
          
      !cor(1)=NL(1,1)*dwd(1)+NL(1,2)*dwd(2)
      !cor(2)=NL(2,1)*dwd(1)+NL(2,2)*dwd(2)

      
      ! if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
      !     massj=mass(nj)
      !
      !else
        !massj=rho_j*dx*dx
      !massi=rho_i*dx*dx
      
    
          !endif
C      write(*,*) cor(1),cor(2),DWDX,DWDY
c      n01=n01+dx*dx*(x_j-x_i)*cor(1)
c           n02=n02+dx*dx*(y_j-y_i)*cor(2)
c	    if(mod(i,y_n).eq.0)then
c	      drho(i)=0.
c	    else
	    !!  drho(i)=drho(i)+massj*
     !!&	  ((vx_i-vx_j)*cor(1)+(vy_i-vy_j)*cor(2))
     !!&	   +deltf*hsml0*c_ij*(fdeltx*cor(1)+fdelty*cor(2))*dx*dx
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
              
c          endif
c          endif
c          else 
c           rgylx=0.0
c           rgyly=0.0
c           endif
           
       if(nj.gt.ntotalb.and.nj.le.(ntotal+nvirt1+nvirt2))then
             eb=0.005 !(epsino1+epsino1)/2.
     !! if((mod(i,6).eq.1.and.x(2,nj).lt.x(2,i)).or.(i.le.6.and.x(1,nj).lt
     !!&.x(1,i)).or.(i.ge.716.and.x(1,nj).gt.x(1,i)).or.(mod(i,6).eq.0
     !!& .and.x(2,nj).gt.x(2,i)).or.(ljsxyb(i).ge.0.5.and.(x(2,nj).ge.
     !!&  0.6)))then!x(2,i).and.x(2,nj).le.x(2,i)
          !if(ip_f(i)==0.and.jc(i)==0)then
         !tsts(1,nj)=-p(nj)+sts(1,nj)
         !tsts(2,nj)=-p(nj)+sts(2,nj)
         !tsts(3,nj)=sts(3,nj)
         !tsts(4,nj)=sts(4,nj)
         ! tsts(1,nj)=0.0
         !tsts(2,nj)=0.0
         !tsts(3,nj)=0.0
         !tsts(4,nj)=0.0
         !cor(1)=0.0
         ! cor(2)=0.0
         !pai_ij=0.0
 !     dvx(1,i)=dvx(1,i)-(-(tsts(1,i)+tsts(2,i))/2.*rho_nj2+p(nj)
 !    &*rho_i2)*cor(1)*massj-pai_ij*cor(1)*massj
 !    & +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
 !    &)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
 !    &*massj  !+rgylx
	!dvx(2,i)=dvx(2,i)-(-(tsts(1,i)+tsts(2,i))/2.*rho_nj2+p(nj)*rho_i2)
 !    &*cor(2)*massj-pai_ij*cor(2)*massj 
 !    &+((sts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
 !    &+((sts2)*rho_i2+(sts(2,nj))*rho_nj2)*cor(2))
 !    &*massj !+rgyly
 !      dvx(1,i)=dvx(1,i)-pai_ij*cor(1)*massj
 !    &+(((tsts1)*rho_i2+(p(nj))*rho_nj2
 !    &)*cor(1)+(tsts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
 !    &*massj !+rgylx
	!dvx(2,i)=dvx(2,i)-pai_ij*cor(2)*massj
 !    &+((tsts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
 !    &+((tsts2)*rho_i2+(p(nj))*rho_nj2)*cor(2))
 !    &	*massj !+rgyly
 !         hsml0=hsmls
	!call kernel_df(x_j,y_j,x_i,y_i,dwdx,dwdy,hsml0)
	!     
	!     dwd(1)=dwdx
	!    dwd(2)=dwdy
 !     fpp1=(-(p(nj)*rho_nj2+p(i)*rho_i2)*dwd(1)*massi
 !    & -pai_ij1*dwd(1)*massi)*massj
 !     fpp2=(-(p(nj)*rho_nj2+p(i)*rho_i2)*dwd(2)*massi
 !    &-pai_ij1*dwd(2)*massi)*massj  
      !!       drho(i)=drho(i)+dx*dy*rho(i)*
     !!&	  ((vx_i-vx_j)*cor(1)+(vy_i-vy_j)*cor(2))  
     !!!!        
     !!      
     !!      if(ip_f(i)==0)then     
     !!      dvx(1,i)=dvx(1,i)-
     !!&	  (p(nj)+p(i))*dwdx*massj/(rho(nj)*rho0)
     !!&      -pai_ij*dwdx*massj
     !!!!&         +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
     !!!!&)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
     !!!!&*massj  !+rgylx        
     !!      
     !!      
	    !!  dvx(2,i)=dvx(2,i)-
     !!&	  (p(nj)+p(i))*dwdy*massj/(rho(nj)*rho0) 
     !!&      -pai_ij*dwdy*massj 
     !!!!&       +((sts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
     !!!!&+((sts2)*rho_i2+(sts(2,nj))*rho_nj2)*cor(2))
     !!!!&	*massj !+rgyly   
     !!       
     !!       ddvx(1,i)=ddvx(1,i)-
     !!&	  (p(nj)+p(i))*dwdx*massj/(rho(nj)*rho0)
     !!!!&      -pai_ij*dwdx*massj
     !!!!&         +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
     !!!!&)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
     !!!!&*massj  !+rgylx        
	    !!  ddvx(2,i)=ddvx(2,i)-
     !!&	  (p(nj)+p(i))*dwdy*massj/(rho(nj)*rho0)  
     !!!!&      -pai_ij*dwdy*massj        
     !!       endif
     !!        dli(1,i)=dli(1,i)-
     !!&	  (p(nj)+p(i))*dwdx*massj/(rho(nj)*rho0)
     !!!!&      -pai_ij*dwdx*massj
     !!!!&         +(((sts1)*rho_i2+(sts(1,nj))*rho_nj2
     !!!!&)*cor(1)+(sts3*rho_i2+sts(3,nj)*rho_nj2)*cor(2))
     !!!!&*massj  !+rgylx        
	    !!  dli(2,i)=dli(2,i)-
     !!&	  (p(nj)+p(i))*dwdy*massj/(rho(nj)*rho0) 
     !!&      -pai_ij*dwdy*massj 
     !!&       +((sts4*rho_i2+sts(4,nj)*rho_nj2)*cor(1)
     !!&+((sts2)*rho_i2+(sts(2,nj))*rho_nj2)*cor(2))
     !!&	*massj !+rgyly        
            
     !! dvx(1,i)=dvx(1,i)+(-(p(nj)*rho_nj2+p(i)*rho_i2)*cor(1)*massj
     !!& -pai_ij*cor(1)*massj)
     !! dvx(2,i)=dvx(2,i)+(-(p(nj)*rho_nj2+p(i)*rho_i2)*cor(2)*massj
     !!&-pai_ij*cor(2)*massj)
     !! 
     !! ddvx(1,i)=ddvx(1,i)+(-(p(nj)*rho_nj2+p(i)*rho_i2)*cor(1)*massj)
     !!& -pai_ij*cor(1)*massj
     !! ddvx(2,i)=ddvx(2,i)+(-(p(nj)*rho_nj2+p(i)*rho_i2)*cor(2)*massj)
     !!&-pai_ij*cor(2)*massj    
 !        do m=1,gs(nj)
 !            if(i==jh(m,nj))then
 !               if(nj==5990.and.m==3)then
 !     write(*,*) "l",fp1(3,5990),i,m
 !     endif      
 !                
 !        dvx(1,i)=dvx(1,i)
 !    &-fp1(m,nj)/massi-pai_ij*cor(1)*massj
 !    !!    +(((tsts1)*rho_i2+(tsts(1,nj))*rho_nj2
 !    !!&)*cor(1)+(tsts3*rho_i2+tsts(3,nj)*rho_nj2)*cor(2))
 !    !!&*massj 
	!dvx(2,i)=dvx(2,i)
 !    & -fp2(m,nj)/massi-pai_ij*cor(2)*massj
 !    !!&+((tsts4*rho_i2+tsts(4,nj)*rho_nj2)*cor(1)
 !    !!&+((tsts2)*rho_i2+(tsts(2,nj))*rho_nj2)*cor(2))
 !    !!&	*massj
 !            endif
 !            enddo
         !    else if(nj.le.ntotalb)then
         !    tsts(1,nj)=0.0
         !tsts(2,nj)=0.0
         !tsts(3,nj)=0.0
         !tsts(4,nj)=0.0
         !pai_ij=0.0 
         !cor(1)=0.0
         ! cor(2)=0.0
         ! w=0.0
         !endif
            !if(ip_f(i)==0)then
     !!     va(1,i)=va(1,i)-2.*eb*massj/(rho_i+rho_j)
     !!&	  *(vx_i-vx_j)*w
     !! va(2,i)=va(2,i)-2.*eb*massj/(rho_i+rho_j)
     !!&	  *(vy_i-vy_j)*w
         !endif
       else if(nj.le.ntotalb)then
           eb=epsino
          ! if(sx(i)==3.and.sx(nj)==4.and.sh(i)==1.and.sh(nj)==1)then
          !if((x(1,nj)-x(1,i)).gt.1.5*dx)then
          !  cor(1)=0.0
          !  cor(2)=0.0
          !  rgylx=0.0
          !  rgyly=0.0
          !  eb=0.0
          !endif
          !  endif
          !if(sx(i)==4.and.sx(nj)==3.and.sh(i)==1.and.sh(nj)==1)then
          !if((x(1,nj)-x(1,i)).lt.-1.5*dx)then
          !  cor(1)=0.0
          !  cor(2)=0.0
          !  rgylx=0.0
          !  rgyly=0.0
          !  eb=0.0
          !endif
          !endif  
        IF(NJ.NE.I.and.i.gt.50.and.nj.gt.50.and.ljsxyb(i).gt.0.0.and
     &     .ljsxyb(nj).gt.0.0)THEN
               dis=sqrt((x(1,i)-x(1,nj))**2
     &+(x(2,i)-x(2,nj))**2)
      dis1=sqrt((xc(1,i)-xc(1,nj))**2+(xc(2,i)-xc(2,nj))**2)
      
       if(dis1.ge.1.2*dx)then
          r0=0.45*dx/dis1
          r01=0.7*dx/dis1
      else
           r0=0.5*dx/dis1
           r01=0.75*dx/dis1
      endif
      if((dis-dis1)/dis1.ge.(r01))then
             
              if(x(1,i).lt.x(1,nj))then     
      sx(i)=3
      sx(nj)=4
      else
        sx(i)=4
      sx(nj)=3
      endif
           cor(1)=0.0
            cor(2)=0.0
            rgylx=0.0
            rgyly=0.0
            eb=0.0
       endif
     
      
      
!!      if((dis-dis1)/dis1.ge.(r0))then
!!              do k=1,rn(i)
!!c           k=k+1
!! 	      nj1=rpn(k,i)
!!            if(nj1.gt.ntotalb.and.nj1.le.ntotal)then  
!!          !x_j=x(1,nj1)
!!          !  y_j=x(2,nj1)
!!          !  z_j=x(3,nj1)                    
!!      if(sqrt(((x(1,nj)+x(1,i))/2.-x(1,nj1))**2+((x(2,nj)+x(2,i))/2.
!!     & -x(2,nj1))**2).le.3.0*dx)then     
!!              if(x(1,i).lt.x(1,nj))then     
!!      sx(i)=1
!!      sx(nj)=2
!!      else
!!        sx(i)=2
!!      sx(nj)=1
!!      endif
!!            cor(1)=0.0
!!            cor(2)=0.0
!!            rgylx=0.0
!!            rgyly=0.0
!!            eb=0.0
!!       endif
!!      endif
!!      enddo
!!          endif        
       ENDIF    
           
          drho(i)=drho(i)+mass(nj)*
     &	  ((vx_i-vx_j)*cor(1)+(vy_i-vy_j)*cor(2))    
             
             
       dvx(1,i)=dvx(1,i)-pai_ij*cor(1)*massj
     &+(((tsts1)*rho_i2+(tsts(1,nj))*rho_nj2
     &)*cor(1)+(tsts3*rho_i2+tsts(3,nj)*rho_nj2)*cor(2))
     &*massj +rgylx
	dvx(2,i)=dvx(2,i)-pai_ij*cor(2)*massj
     &+((tsts4*rho_i2+tsts(4,nj)*rho_nj2)*cor(1)
     &+((tsts2)*rho_i2+(tsts(2,nj))*rho_nj2)*cor(2))
     &	*massj +rgyly
      
     !!   ddvx(1,i)=ddvx(1,i)-pai_ij*cor(1)*massj
     !!&+(((tsts1)*rho_i2+(tsts(1,nj))*rho_nj2
     !!&)*cor(1)+(tsts3*rho_i2+tsts(3,nj)*rho_nj2)*cor(2))
     !!&*massj +rgylx 
	    !!  ddvx(2,i)=ddvx(2,i)-pai_ij*cor(2)*massj
     !!&+((tsts4*rho_i2+tsts(4,nj)*rho_nj2)*cor(1)
     !!&+((tsts2)*rho_i2+(tsts(2,nj))*rho_nj2)*cor(2))
     !!&	*massj +rgyly
      
c      if(nj.le.ntotalb)then
       va(1,i)=va(1,i)-2.*eb*massj/(rho_i+rho_j)
     &	  *(vx_i-vx_j)*w
      va(2,i)=va(2,i)-2.*eb*massj/(rho_i+rho_j)
     &	  *(vy_i-vy_j)*w
            endif
      
c      if(nj.gt.ntotal)then
c      jsd(i)=jsd(i) -(p(nj)+p(i))*cor(2)*mass/(rho(nj)*rho(i))
c     &-pai_ij*cor(2)*mass+((sts(4,i)/rho(i)**2+sts(4,nj)/rho(nj)**2
c     &)*cor(1)+(sts(2,i)/rho(i)**2+sts(2,nj)/rho(nj)**2)*cor(2))*mass !+rgyly	
c      endif
		  !endif
           
      enddo
      !if(mod(i,7).eq.1)then !.or.i.le.7.or.i.ge.344
      !!dvx(1,i)=dvx(1,i)+p(i)/mass(i)
      !dvx(2,i)=dvx(2,i)+p(i)/mass(i)
      !endif
c        write(55,*) p(11),sts(4,i),sts(2,i),dwdx(2),pai_ij,dvx(2,11) 
c          if(n01.ne.0.0.and.n02.ne.0.0)then  
c          dvx(1,i)=dvx(1,i)/n01
c          dvx(2,i)=dvx(2,i)/n02
c          endif
          !if(i.gt.ntotalb)then
	     dvx(2,i)=dvx(2,i)-9.81
           ddvx(2,i)=ddvx(2,i)-9.81
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
           !if(ip_f(i)==0)then
           if(rhox(i)==0.0)then
                rhox(i)=rho0
               endif
          dvx(1,i)=dvx(1,i)+dyl(1,i)/mass(i) 
            dvx(2,i)=dvx(2,i)+dyl(2,i)/mass(i)   
           ddvx(1,i)=ddvx(1,i)+dyl(1,i)/mass(i) 
            ddvx(2,i)=ddvx(2,i)+dyl(2,i)/mass(i)
             va(1,i)=va(1,i)+dyva(1,i)/mass(i) 
            va(2,i)=va(2,i)+dyva(2,i)/mass(i)  
            !endif
      enddo
      
      !!  li1=0.0
      !!li2=0.0
      !do i=1,ntotalb
      !ddvx(1,i)=ddvx(1,i)+dvx(1,i)
      !ddvx(2,i)=ddvx(2,i)+dvx(2,i)
      !if(mod(t,33)==2)then
      ! ddvx(1,i)=0.0
      !ddvx(2,i)=0.0
      !endif
      !
      !enddo
          !write(*,*) dvx(1,6),ddvx(1,6)
          
          
      !!do i=1,ntotalb
      !!    if(ip_f(i)==0)then    
      !!    !nj1=nj1+ddvx(1,i)
      !!    !nj2=nj2+dou(2,i)
      !!    li1=li1+dli(1,i)*mass(i)
      !!    li2=li2+dli(2,i)*mass(i)
      !!    endif
      !!enddo
     !!  do i=1,ntotalb
     !!      !if(ip_f(i)==0)then    
     !! ddvx(1,i)=li1/mz!+nj1/(1./12.*mz*(xl**2+yl**2))*
     !!!!& (x(1,i)-x_zx)!/massi!-pai_ij*cor(1)*massj
     !! ddvx(2,i)=li2/mz-9.81!+nj2/(1./12.*mz*(xl**2+yl**2))*
     !!!!& (x(2,i)-y_zx) !/massi!-pai_ij*cor(2)*massj
     !! !endif
     !! 
     !! !dvx(1,i)=dvx(1,i)+li1/mz
     !! !dvx(2,i)=dvx(2,i)+li2/mz
     !! enddo

	end
      
      
      



      