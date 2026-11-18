      subroutine s_grid
      
	use rn_grid

	implicit none

	include 'menu.inc'

      

	double precision x_i,y_i,x_j,y_j

	nx_g=int(x0/(nd*dx1)+0.1*dx1)
      ny_g=int(y0/(nd*dx1)+0.1*dx1)
 



      k=0
	do i=1,nx_g+2

	    do j=1,ny_g+2
      k=k+1  
	   enddo

	enddo


	ngrid=k
 

c	allocate(rn_g(1:ngrid),rpn_g(1:10,1:ngrid))
c      allocate(tpnioc(1:ngrid),pgioc(1:maxnr*maxnr_g,1:ngrid))
c	allocate(x_g(1:dim,1:ngrid))
      allocate(mesh(0:maxnr,0:(ny_g+2),0:(nx_g+2)))

c      k=0
	
	
c	do j=1,ny_g+2
	    
		
		
c		do i=1,nx_g+2
c      k=k+1

c      x_g(1,k)=px_0-0.5*nd*dx+(i-1)*nd*dx
c	x_g(2,k)=py_0-0.5*nd*dx+(j-1)*nd*dx

c	   enddo

c	enddo



c	do i=1,ngrid

c      x_i=x_g(1,i)
c	y_i=x_g(2,i)
c	k=0

c	     do j=1,ngrid

c	x_j=x_g(1,j)
c	y_j=x_g(2,j)

c	dis=sqrt((x_i-x_j)**2+(y_i-y_j)**2)
c      if(dis.le.1.5*nd*dx)then
c      k=k+1
c      rpn_g(k,i)=j
c	endif
     
c	     enddo
c      rn_g(i)=k
c      enddo




	end subroutine
