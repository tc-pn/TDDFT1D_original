!-------------------------
c small program to compute 
c the 1D EOS 
c-------------------------

       program EOS1D
       
       IMPLICIT NONE
       INTEGER Nrho,I
       parameter (Nrho=50)
       REAL*8 drho, rho
       parameter (drho=0.05)
       REAL*8 n0,e0,ee
       
       REAL*8 c2,beta,cbeta
       REAL*8 PI,xk0
       REAL*8 MASS,HBAR
       
       MASS=938.91897D0
       HBAR=197.32705D0 


!--------cte       
       PI = ACOS(-1.)
       xk0 = HBAR**2./(2.D0*MASS)*PI*PI/12.
               
       print*,'xk0',xk0        
       e0   = -16.     ! MeV  
       e0   = -1.1*xk0     ! MeV 
       n0   = 0.5    ! fm-1
       beta = 3. ! no unit
!       beta = 2.+1./3. ! no unit
!       beta = 2.+1./6. ! no unit
       
       cbeta = -(e0 + xk0*n0*n0)
       cbeta = cbeta/(beta-1.)/n0**beta
       c2    = -2.*xk0*n0*n0-cbeta*beta*n0**beta
       c2    = c2/n0  
       
       rho = -drho 
       do I=1,Nrho
          rho = rho + drho
          ee = xk0*rho**2.+c2*rho+cbeta*rho**beta
          write(10,*) rho,ee
       enddo
       
       print*,'---beta---',beta
       print*,'c2=',c2
       print*,'cbeta=',cbeta
           
       
       end