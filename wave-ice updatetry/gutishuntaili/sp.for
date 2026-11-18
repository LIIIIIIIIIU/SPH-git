       subroutine sp(it,ntotal1)
      
      use rn_grid
      implicit none
      
	include 'menu.inc'

      integer ntotal1,it
     

      integer nix,niy,in_g,n_i,n_j

	double precision xi,yi,xj,yj
	
	integer i0,ip,ii1,j1,ii2,i00,j00,
     &j0,jj,c_j,ip1,temp,j01

      double precision dis_0,dis_1,dis


        do i=0,nx_g+2
      
        do j=0,ny_g+2
        
        mesh(0,j,i)=0
      
        enddo
        
        enddo

      
	
	do i=1,ntotal1
        if(x(1,i).ge.5.0.and.x(1,i).le.11.0.and.x(2,i).ge.0.3)then  
      rn(i)=0
      rn2(i)=0
      
      xi=x(1,i)
      yi=x(2,i)
      

      
      
      n_i=floor ((xi-wmx)/(nd*dx1)+0.1)+1
      n_j=floor ((yi-py_0)/(nd*dx1)+0.1)+1
     
      !!!!!!!!!Cut off the mistake information!!!!!!
      if(n_i.lt.1)n_i=1
      if(n_i.gt.nx_g+1)n_i=nx_g+1
      
      if(n_j.lt.1)n_j=1
      if(n_j.gt.ny_g+1)n_j=ny_g+1
      

      
      
      mesh(0,n_j,n_i)=mesh(0,n_j,n_i)+1      !某网格中的粒子数
      mesh(mesh(0,n_j,n_i),n_j,n_i)=i       !粒子在网格中的号
      endif

      enddo
      
      
      !$omp parallel  private(i,j,i0,i00,xi,yi,ip,ip1)
     &private(j0,j00,xj,yj,dis)
      !$omp do schedule(static)
      !!!!!Search neighbour particles!!!!!!
      do i=1,nx_g+1
      
          do j=1,ny_g+1
       
           do i0=1,mesh(0,j,i)
       
           i00=mesh(i0,j,i)
           
      if(x(1,i00).ge.5.0.and.x(1,i00).le.11.0.and.x(2,i00).ge.0.3)then
!           if(it.ge.1095.and.i00.eq.9001)then
!           write(*,*)
!           endif
           ip=0
           ip1=0
           xi=x(1,i00)
           yi=x(2,i00)
           
           
           !!!!!!!!!!Center!!!!!!!!
           do j0=1,mesh(0,j,i)
           
           j00=mesh(j0,j,i)  !!
       
           xj=x(1,j00)
           yj=x(2,j00)
           
           dis=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)
           
         
            if( dis.le.s_rad*s_rad*dx*dx)then
		    ip1=ip1+1
		    rpn(ip1,i00)=j00
            endif
            
             if( dis.le.s_rad1*s_rad1*dx1*dx1)then
		    ip=ip+1
		    rpn2(ip,i00)=j00
		    endif
            
		    enddo
        !!!!!!!!!!!!!E!!!!!!!!!!!!!!!!!!!
          
           do j0=1,mesh(0,j,i+1) !!
           
           j00=mesh(j0,j,i+1)  !!
       
           xj=x(1,j00)
           yj=x(2,j00)
           
          dis=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)
           
          
            if( dis.le.s_rad*s_rad*dx*dx)then
		    ip1=ip1+1
		    rpn(ip1,i00)=j00
            endif
            
             if( dis.le.s_rad1*s_rad1*dx1*dx1)then
		    ip=ip+1
		    rpn2(ip,i00)=j00
		    endif
            
		    enddo
            
        !!!!!!!!!!!!!ES!!!!!!!!!!!!!!!!!!
        
           do j0=1,mesh(0,j-1,i+1) !!
           
           j00=mesh(j0,j-1,i+1)  !!
       
           xj=x(1,j00)
           yj=x(2,j00)
           
           dis=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)
           
          
            if( dis.le.s_rad*s_rad*dx*dx)then
		    ip1=ip1+1
		    rpn(ip1,i00)=j00
            endif
            
             if( dis.le.s_rad1*s_rad1*dx1*dx1)then
		    ip=ip+1
		    rpn2(ip,i00)=j00
		    endif
		    enddo
         
        
        
        !!!!!!!!!!!!S!!!!!!!!!!!!!!!!!!!!
        
           do j0=1,mesh(0,j-1,i) !!
            
           j00=mesh(j0,j-1,i)  !!
       
           xj=x(1,j00)
           yj=x(2,j00)
           
           dis=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)
           
         
            if( dis.le.s_rad*s_rad*dx*dx)then
		    ip1=ip1+1
		    rpn(ip1,i00)=j00
            endif
            
             if( dis.le.s_rad1*s_rad1*dx1*dx1)then
		    ip=ip+1
		    rpn2(ip,i00)=j00
		    endif
		    enddo
          
        
        
        !!!!!!!!!!!!WS!!!!!!!!!!!!!!!!!!!
        
           do j0=1,mesh(0,j-1,i-1) !!
           
           j00=mesh(j0,j-1,i-1)  !!
       
           xj=x(1,j00)
           yj=x(2,j00)
           
           dis=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)
           
         
            if( dis.le.s_rad*s_rad*dx*dx)then
		    ip1=ip1+1
		    rpn(ip1,i00)=j00
            endif
            
             if( dis.le.s_rad1*s_rad1*dx1*dx1)then
		    ip=ip+1
		    rpn2(ip,i00)=j00
		    endif
		    enddo
          
        
        
        !!!!!!!!!!!!W!!!!!!!!!!!!!!!!!!!!
         
           do j0=1,mesh(0,j,i-1) !!
           
           j00=mesh(j0,j,i-1)  !!
       
           xj=x(1,j00)
           yj=x(2,j00)
           
           dis=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)
           
         
            if( dis.le.s_rad*s_rad*dx*dx)then
		    ip1=ip1+1
		    rpn(ip1,i00)=j00
            endif
            
             if( dis.le.s_rad1*s_rad1*dx1*dx1)then
		    ip=ip+1
		    rpn2(ip,i00)=j00
		    endif
		    enddo
          
        
        
        !!!!!!!!!!!!WN!!!!!!!!!!!!!!!!!!!
         
           do j0=1,mesh(0,j+1,i-1) !!
            
           j00=mesh(j0,j+1,i-1)   !!
       
           xj=x(1,j00)
           yj=x(2,j00)
           
           dis=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)
           
         
            if( dis.le.s_rad*s_rad*dx*dx)then
		    ip1=ip1+1
		    rpn(ip1,i00)=j00
            endif
            
             if( dis.le.s_rad1*s_rad1*dx1*dx1)then
		    ip=ip+1
		    rpn2(ip,i00)=j00
		    endif
		    enddo
        
        
        
        
        !!!!!!!!!!!!N!!!!!!!!!!!!!!!!!!!!
          
           do j0=1,mesh(0,j+1,i) !!
           
           j00=mesh(j0,j+1,i)  !!
       
           xj=x(1,j00)
           yj=x(2,j00)
           
           dis=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)
           
        
            if( dis.le.s_rad*s_rad*dx*dx)then
		    ip1=ip1+1
		    rpn(ip1,i00)=j00
            endif
            
             if( dis.le.s_rad1*s_rad1*dx1*dx1)then
		    ip=ip+1
		    rpn2(ip,i00)=j00
		    endif
		    enddo
          
        
        
        !!!!!!!!!!!!EN!!!!!!!!!!!!!!!!!!!
        
           do j0=1,mesh(0,j+1,i+1) !!
           
           j00=mesh(j0,j+1,i+1)  !!
       
           xj=x(1,j00)
           yj=x(2,j00)
           
           dis=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)
           
         
            if( dis.le.s_rad*s_rad*dx*dx)then
		    ip1=ip1+1
		    rpn(ip1,i00)=j00
            endif
            
             if( dis.le.s_rad1*s_rad1*dx1*dx1)then
		    ip=ip+1
		    rpn2(ip,i00)=j00
		    endif
		    enddo
          
       
       
          rn(i00)=ip1
          rn2(i00)=ip
          
                endif
          enddo
         
      enddo
      

      enddo
	
	 !$omp end do
       !$omp end parallel
	
	
	!!!subroutine bubble_sort(a,n)
       !implicit none
!       do i=1,ntotal+nvirt0
!       
!          do j=rn(i),1,-1
!          
!          
!           do j1=1,j
!          
!           if(rpn(j,i).lt.rpn(j1,i))then
!           
!           temp=rpn(j,i)
!           rpn(j,i)=rpn(j1,i)
!           rpn(j1,i)=temp
!           endif
!           enddo
!          enddo
!        enddo
!       
!      do i=1,ntotal+nvirt0
!       
!          do j=rn2(i),1,-1
!          
!          
!           do j1=1,j
!          
!           if(rpn2(j,i).lt.rpn2(j1,i))then
!           
!           temp=rpn2(j,i)
!           rpn2(j,i)=rpn2(j1,i)
!           rpn2(j1,i)=temp
!           endif
!           enddo
!          enddo
!        enddo
	
	
	end