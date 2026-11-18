	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		integer maxn,maxnr,maxnr_g
	parameter (maxn=450000,maxnr=80,maxnr_g=10)
	
	integer dim
	parameter (dim=2)

	double precision omiga,rho0,pi,E,miu
	parameter (omiga=7.0,rho0=900.,pi=3.1415926,E=1.3865e+8,miu=0.33)
      
	double precision dt,g,rho01,mius
	parameter (dt=0.000005,g=9.81,rho01=1000.,mius=0.1)
	  integer screan_p,sa_n
	parameter (screan_p=50000,sa_n=10000)

	real s_rad,s_rad1
	parameter (s_rad=3.1,s_rad1=3.1)
      
	integer nd,bs
	parameter (nd=3,bs=45)

	real vs,bj,np       
	parameter(vs=-0.00,bj=0.022,np=6.0)

	real alpha_phai,belda_phai,epsino,epsino1
	parameter (alpha_phai=1.0,belda_phai=1.0,epsino=0.1,epsino1=0.002)
	real alpha_phai1,belda_phai1,deltf
	parameter (alpha_phai1=0.2,belda_phai1=0.0,deltf=0.1)

       integer art_vis,kernel_type
	   parameter (art_vis=2,kernel_type=1)

	    real jiao,jzj,yt
		parameter (jiao=pi/5.,jzj=pi/15.,yt=4.0)
	real(kind=8),parameter::h=1.0,lf=16.,ld=3.        !!!ld

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Variable!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer i,j,k,k1,nc
	double precision x(dim,maxn),vx(dim,maxn),dvx(dim,maxn),rho(maxn)
     &,drho(maxn),mass(maxn),hsml,va(dim,maxn),hsmls,
     &sts(4,maxn),p(maxn),rhox(maxn),bxvx(dim,3000),xc(dim,3000)
	 double precision yb(4,maxn),ddvx(dim,3000),dx1,dy1
	 
	  double precision tsts(4,maxn),dyl(dim,3000),cdyl(dim,3000)
	  double precision qfqd(maxn),ljsxyb(maxn),dyva(dim,3000)  !!!??
	  integer x_n,y_n,ntotal,nvirt1,nvirt2,nvirt3,ntotals,ntotalb
	  double precision px_0,py_0,xl,yl,x0,y0,px_1,py_1,dx,dy,co
	  double precision xsl,ysl,c01,dpx(maxn),dpy(maxn),ybl(4,3000)
	  double precision x_0(dim,maxn),k1vx(dim,maxn),k1drho(maxn)
	  double precision k1dvx(dim,maxn),rho_0(maxn),vx_0(dim,maxn)
	  double precision NL1(dim,maxn),NL2(dim,maxn),v(maxn)
	  double precision zjqd(maxn),zjpyb(maxn),rg(4,maxn)
	  double precision vq(maxn),vf(maxn),pz(maxn),free_surface(7000)
      integer xs_n,ys_n,n_f,fre_num(3000),jc(3000),sh(maxn),zjbj(3000)
        double precision arf,a_rf,sgm,JG,KG,q0,wmv,wma,wmx
       integer rn(maxn),rpn(maxnr,maxn),rn1(maxn),rpn1(maxnr,maxn)
	   integer rn2(maxn),rpn2(maxnr,maxn),wq(maxn),zbj(3000),fs(maxn)
	   integer ngrid,nx_g,ny_g,current,cs(maxn),sx(maxn),ip_f(3000)    

	   common//x,vx,dvx,rho,drho,mass,hsml,hsmls,va,co,yb,c01,dyva
	   common//vx_0,rho_0,x_0,k1vx,k1drho,k1dvx,v,vq,vf,pz,free_surface
	  common/pohuai/sts,tsts,qfqd,ljsxyb,p,dyl,cdyl,ddvx,dx1,dy1,zjqd,zjpyb

	  common/scale/px_0,py_0,xl,yl,xsl,ysl,x0,y0,px_1,py_1,dx,dy,bxvx
	  common/grid/x_n,y_n,xs_n,ys_n,wq,zbj,zjbj
	  common/rn/rn,rpn,rn1,rpn1,rn2,rpn2,fs
	  common/ntotal/ntotal,nvirt1,nvirt2,nvirt3,ntotals,ntotalb
      common/rn_g/ngrid,nx_g,ny_g 
	  common/time/current,cs,nc,sx,ip_f,n_f,fre_num,jc 
	  common/bottom/arf,a_rf,sgm,JG,KG,q0,wmv,wma,wmx,rhox,xc
	  common/delta/dpx,dpy,NL1,NL2,rg,sh,ybl