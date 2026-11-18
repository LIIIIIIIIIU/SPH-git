        subroutine kernel(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
	implicit none
	include 'menu.inc'
      !kernel subroutine subprogram: it is used for calculating the
	!relue, devirative relue and second order devirative relue of the kernel function  
      integer d   
	double precision x_i,y_i,x_j,y_j,hsml0,w,factor
     &,q,dx_i,dy_i,dwdx,dwdy
      
	!typei=1: cubic spline function
	!typei=2: quintic spline function

	
	dx_i = x_i-x_j
	dy_i = y_i-y_j
      
	
      

	q = sqrt((x_i-x_j)**2+(y_i-y_j)**2)/hsml0
      

        
	     factor = 15.e0/(7.e0*pi*hsml0*hsml0)
        if (q.ge.0.and.q.le.1.e0) then
          w = factor * (2./3. - q*q + q**3 / 2.)
        dwdx = factor * (-2.+3./2.*q)*dx_i/(hsml0**2)
		dwdy = factor * (-2.+3./2.*q)*dy_i/(hsml0**2)
   
        else if (q.gt.1.e0.and.q.le.2) then          
       w = factor * 1.e0/6.e0 * (2.-q)**3   
       
       dwdx= -factor * 1./2.*(2.-q)**2*(dx_i/(hsml0*hsml0*q)) 
	 dwdy= -factor * 1./2.*(2.-q)**2*(dy_i/(hsml0*hsml0*q))  
	  else
	  w=0.
          dwdx= 0.
            dwdy= 0.  
       
        endif
      


      end subroutine
      
       subroutine kernel3(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
	implicit none
	include 'menu.inc'
      !kernel subroutine subprogram: it is used for calculating the
	!relue, devirative relue and second order devirative relue of the kernel function  
      integer d   
	double precision x_i,y_i,x_j,y_j,hsml0,w,factor
     &,q,dx_i,dy_i,dwdx,dwdy
      
      dx_i = x_i-x_j
	dy_i = y_i-y_j
      
	q = sqrt((x_i-x_j)**2+(y_i-y_j)**2)/hsml0
      
       factor = 7.e0 / (478.e0*pi*hsml0*hsml0)
            
	if(q.eq.0) then
         w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
      
	   dwdx = factor * (-5.*(3-q)**4+30.*(2-q)**4
     &   -75.*(1-q)**4)/hsml0*dx_i
         dwdy = factor * (-5.*(3-q)**4+30.*(2-q)**4
     &   -75.*(1-q)**4)/hsml0*dy_i

	else if(q.gt.0.and.q.le.1)then
         w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
      
	   dwdx = factor * (-5.*(3-q)**4+30.*(2-q)**4
     &   -75.*(1-q)**4)* (dx_i/(hsml0**2*q)) 
         dwdy = factor * (-5.*(3-q)**4+30.*(2-q)**4
     &   -75.*(1-q)**4)* (dy_i/(hsml0**2*q)) 
          
	else if(q.gt.1.and.q.le.2) then
          w = factor * ( (3-q)**5 - 6*(2-q)**5 )

          dwdx= factor * (-5.*(3-q)**4+30.*(2-q)**4)
     &   * (dx_i/(hsml0**2*q)) 
          dwdy= factor * (-5.*(3-q)**4+30.*(2-q)**4)
     &   * (dy_i/(hsml0**2*q)) 

        else if(q.gt.2.and.q.le.3) then
          w = factor * (3-q)**5 
          dwdx= factor * (-5.*(3-q)**4) * (dx_i/(hsml0**2*q)) 
          dwdy= factor * (-5.*(3-q)**4) * (dy_i/(hsml0**2*q)) 
        else   
	    w = 0.
          dwdx= 0.
          dwdy= 0.
        
        endif           
      
      end subroutine
      
      
      subroutine kernel_f(x_i,y_i,x_j,y_j,w,hsml0)
	implicit none
	include 'menu.inc'
      !kernel subroutine subprogram: it is used for calculating the
	!relue, devirative relue and second order devirative relue of the kernel function  
      integer d   
	double precision x_i,y_i,x_j,y_j,hsml0,w,factor
     &,q,dx_i,dy_i
      
	!typei=1: cubic spline function
	!typei=2: quintic spline function

	
	dx_i = x_i-x_j
	dy_i = y_i-y_j
      
	
      

	q = sqrt((x_i-x_j)**2+(y_i-y_j)**2)/hsml0
      

        
	     factor = 15.e0/(7.e0*pi*hsml0*hsml0)
        if (q.ge.0.and.q.le.1.e0) then
          w = factor * (2./3. - q*q + q**3 / 2.)
       
   
        else if (q.gt.1.e0.and.q.le.2) then          
       w = factor * 1.e0/6.e0 * (2.-q)**3           
	  else
	  w=0.
         
       
        endif
      


      end subroutine
      
      subroutine kernel_df(x_i,y_i,x_j,y_j,dwdx,dwdy,hsml0)
	implicit none
	include 'menu.inc'
      !kernel subroutine subprogram: it is used for calculating the
	!relue, devirative relue and second order devirative relue of the kernel function  
      integer d   
	double precision x_i,y_i,x_j,y_j,hsml0,dwdx,dwdy,factor
     &,q,dx_i,dy_i
      
	!typei=1: cubic spline function
	!typei=2: quintic spline function

	
	dx_i = x_i-x_j
	dy_i = y_i-y_j
      
	
      

	q = sqrt((x_i-x_j)**2+(y_i-y_j)**2)/hsml0


         ! factor = 1.e0/hsml
	     factor = 15.e0/(7.e0*pi*hsml0*hsml0)
        if (q.ge.0.and.q.le.1.e0) then
         
          dwdx = factor * (-2.+3./2.*q)*dx_i/(hsml0**2)
		dwdy = factor * (-2.+3./2.*q)*dy_i/(hsml0**2)
   
        else if (q.gt.1.e0.and.q.le.2) then          
       
       dwdx= -factor * 1./2.*(2.-q)**2*(dx_i/(hsml0*hsml0*q)) 
	 dwdy= -factor * 1./2.*(2.-q)**2*(dy_i/(hsml0*hsml0*q))      
             
	  else
	 
            dwdx= 0.
            dwdy= 0.        
        endif
      

      end subroutine
      
      
      subroutine kernelw(x_i,y_i,x_j,y_j,w,dwdx,dwdy,hsml0)
	implicit none
	include 'menu.inc'
      !kernel subroutine subprogram: it is used for calculating the
	!relue, devirative relue and second order devirative relue of the kernel function  
      integer d   
	double precision x_i,y_i,x_j,y_j,hsml0,w,factor
     &,q,dx_i,dy_i,dwdx,dwdy
      
      dx_i = x_i-x_j
	dy_i = y_i-y_j
      
	q = sqrt((x_i-x_j)**2+(y_i-y_j)**2)/hsml0
      
       factor = 7.e0 / (4.e0*pi*hsml0*hsml0)
            
	if(q.eq.0) then
         w = factor * (1-q/2.)**4*(2.*q+1) 
      
	   dwdx = factor * (-2.*(1-q/2.)**3*(2.*q+1)+2.*(1-q/2.)**4
     &   )/hsml0*dx_i
         dwdy = factor * (-2.*(1-q/2.)**3*(2.*q+1)+2.*(1-q/2.)**4
     &   )/hsml0*dy_i

	else if(q.gt.0.and.q.le.2)then
          w = factor * (1-q/2.)**4*(2.*q+1)
      
	   dwdx = factor * (-2.*(1-q/2.)**3*(2.*q+1)+2.*(1-q/2.)**4
     &   )* (dx_i/(hsml0**2*q)) 
         dwdy = factor * (-2.*(1-q/2.)**3*(2.*q+1)+2.*(1-q/2.)**4
     &   )* (dy_i/(hsml0**2*q)) 
          
        else   
	    w = 0.
          dwdx= 0.
          dwdy= 0.
        
        endif           
      
      end subroutine