*FROM THE GORING (1978),
*XP(T)=H/K*[TANH(X(T))+TANH(K/D*LAMDA)],X(T)=K/D*(C*T-XP(T)-LAMDA)
*K=SQRT(3H/4D),C=SQRT(G*(D+H))
*UP(T)=C*H/D*(1/(COSH(X(T))**2.+H/D)
*LAMDA=L/K*D, L=ARCCOSH(1/SQRT(EZ)),EZ=0.002
       SUBROUTINE wavemaker(t)
      implicit none
	include 'menu.inc'
      double precision t,wd,wh,cc,S,period,tx,tx1,cerr,err,ck,ck1,tu
     &	,ta,xt,ck2,wmx1,xt1,dt1
	
c	g=9.81
      dt1=dt*bs
      cerr=0.0001
      wh=0.10
	wd=ysl
	cc=sqrt(g*(wd+wh))
      ck=sqrt(4.*wh/(3.*wd))
	ck1=sqrt(3.*wh/(4.*wd**3))
	ck2=sqrt(3.*wh/(4.*wd))
	s=2.*wd*ck
      period=2./cc*sqrt(4.*wd**3/(3.*wh))*(3.8+wh/wd)
      tx1=0.0
	tx=wmx+wmv*dt1
c      xt=ck1*(cc*t-wmx)
c      wmx1=ck*wd*tanh(xt)
100      xt=ck1*(cc*(t-period/2.)-(tx-s/2.))
c      wmx1=ck*wd*(tanh(xt)+tanh(ck1*(cc*period/2.-s/2.)))
      tx1=ck*wd*(tanh(xt)+0.999)
c      xt=ck2/wd*(cc*t-wmx-xl)
c      wmx1=wh/ck2*(tanh(xt)+tanh(ck2/wd*xl))
c      xt=ck2/wd*(cc*t-tx-xl)
      err=ABS(tx1-tx)
	tx=tx1
	IF(err.gt.cerr)GOTO 100
     
      tu=cc*wh/(wd*cosh(-xt)**2+wh)
       
      ta=2*wh*cc*wd*COSH(-xt)*SINH(-xt)*(-ck1*cc+ck1*tu)/
     &	(wd*cosh(-xt)**2+wh)**2
c      ta=-2*wh*cc*wd*COSH(ck1*(tx-cc*t))*SINH(ck1*(tx-cc*t))*
c     &	(-ck1*cc+ck1*tu)/(wd*cosh(ck1*(tx-cc*t))**2+wh)**2 +ck1*tu 
      wma=ta
	wmv=tu
	wmx=tx+0.1
      RETURN 
	end




      !SUBROUTINE wavemaker(WMA,WMV,WMX,TTT,HEIGHT,H0,DT,G)
      !implicit none
	!include 'menu.inc'
      !real(kind=4) :: WMA,WMV,WMX,TTT,HEIGHT,H0,DT,G
	!real(kind=4) :: EZ,CL,CK,CK1,C,T0,SS,CLAMDA,CERR,
     !!&XPI,XP,XT,ERR,UP,VP
      
	!HEIGHT=0.45
	!H0=1.0
	!dt=0.0001
      !EZ=0.002
	!CL=LOG((1./(EZ**0.5))*(1+SQRT(1.-EZ)))
c      CL=ARCCOSH(1./EZ**0.5)
	!CK=SQRT(3.*HEIGHT/(4.*H0))
      !CK=SQRT(4.*HEIGHT/(3.*H0))
	!CK1=SQRT(3.*HEIGHT/(4.*H0**3))
	!C=SQRT(g*(H0+HEIGHT))

	!CLAMDA=CL/CK*H0
	!CERR=0.0001
      !XPI=0.0
	!XP=WMX+WMV*DT
      !T0=2*(3.8+HEIGHT/H0)/(C*CK1)
	!SS=SQRT(16*HEIGHT/(3*H0))

      !!!100   !XT=CK/H0*(C*TTT-XP-CLAMDA)
      !XPI=CK*H0*(1+TANH(CK1)*(C*(TTT-T0/2)-(XPI-SS/2)))
      !XPI=HEIGHT/CK*(TANH(XT)+TANH(CK*CLAMDA/H0))
	!ERR=ABS(XPI-XP)
	!XP=XPI
	!IF(ERR.GT.CERR)GOTO 100


      !UP=HEIGHT*C/(H0*(COSH(XT)**2.+HEIGHT/H0))
	!VP=2*HEIGHT*C*COSH(XT)*SINH(XT)*CK*C
	!VP=-SQRT(3.)*HEIGHT**1.5*(HEIGHT+H0)*COSH(XT)**3.*SINH(XT)
	!VP=VP/(H0*(COSH(XT)**2.+HEIGHT/H0))**2

	!WMA=VP  !ACCELERATE
	!WMV=UP  !VELOCITY
	!WMX=XP  !DISPLACEMENT

	!RETURN 
	!END
