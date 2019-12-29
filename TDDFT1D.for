c------------------------------------------------------------
c     Wave-packet evolution using 3 types of evolution (Ipro) from best to worth
c       Ipro =0 -- RK2 algorithm        -->evol_RK4 
c       Ipro =1 -- RK2 algorithm        -->evol_RK2 
c       Ipro =2 -- Leap frog algorithm  -->evol_LF
c       Ipro =3 -- Crank Nicholson      -->ecol_CN  ! worse--to avoid for non-occupied states 
c------------------------------------------------------------
c   Included
c       -Absorbing boundary conditions
c------------------------------------------------------------
      subroutine evol_mf

      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      INTEGER     I,J,IT,K,L
      COMPLEX*16  wave_sp_1(0:Noeud,Nmax_sp)
      COMPLEX*16  wave_sp_2(0:Noeud,Nmax_sp)
      REAL*8      prob_in,norm,prob_in_theo
      REAL*8      delta_t2
      REAL*8      time    
      REAL*8      Etot
      INTEGER     Ipro
      parameter  (Ipro = 0) 
       

      delta_t2 = 0.5*delta_t

c      open(unit=27,file='decr.dat',status='unknown')


!      DO I = 0,Noeud
!	    pot_imag(I) = 0. ! no imaginary
!        pot_one(I)  = 0. 
!      enddo

      time = -delta_t 
      
      do IT = 1,ntime
         time = time +delta_t 

         IF(Ipro .EQ. 0) THEN
            call evol_RK4(delta_t) 
         ENDIF
         IF(Ipro .EQ. 1) THEN
            call evol_RK2(delta_t) 
         ENDIF
            
         If((Ipro .EQ. 2) .OR. (Ipro .EQ. 3)) THEN
            
! mean-field calculation + one-body Hamiltonian     
             call dens_one(wave_sp)
	         call MEAN_FIELD
             call create_hprop !contain external+mf+absorption
         
!propagate up to dt/2
             DO j = 1,nmax_sp
                DO i = 0,noeud        
                   wave_sp_2(i,j) = wave_sp(i,j)
                ENDDO
             ENDDO
c         call evol_tddft(delta_t2,wave_sp,wave_sp_2) ! a bit less stable numerically
             IF(Ipro .EQ. 2) THEN
                call evol_LF(delta_t2,wave_sp,wave_sp_2) 
             ENDIF
             IF(Ipro .EQ. 3) THEN
                call evol_CN(delta_t2,wave_sp,wave_sp_2)
             ENDIF
             
         
c----------------------
c new Hamiltonian at t+dt/2
c-----------------------         
             call dens_one(wave_sp_2)
	         call MEAN_FIELD
	         call create_hprop
c--------------------------------
c- propagation from t to t+dt
c--------------------------------	     
c         call evol_tddft(delta_t,wave_sp,wave_sp_2) ! a bit less stable numerically
             IF(Ipro .EQ. 2) THEN
                call evol_LF(delta_t,wave_sp,wave_sp_2)
             ENDIF
             IF(Ipro .EQ. 3) THEN
                call evol_CN(delta_t,wave_sp,wave_sp_2)
             ENDIF
             DO j = 1,nmax_sp
                DO i = 0,noeud
                   wave_sp  (i,j) = wave_sp_2(i,j)
                ENDDO
             ENDDO       
          ENDIF   ! end Ipro=2 or 3
          
c-------------------------------
c Absorption at boundary
c-------------------------------
         call evol_IMAG
c-------------------------------
c Some output --
c-------------------------------      
         if(MOD(IT-1,Nprint) .EQ. 0) then   
            print*,'Itime',IT 
            call Output(time)
            call compute_emf(etot)
            call test_ortho(wave_sp)
            do I=0,Noeud
               write(100+NINT(float(IT-1)/float(Nprint)),*) 
     &           pas*float(I), densite(I)
            enddo          
         endif   
         
      enddo ! End of time loop
      
      call dens_one(wave_sp)
!      do I=0,noeud
!         write(30,*) float(I)*pas,densite(I)
!      ENDDO
      DO I=0,Noeud
	     write(30,*) SNGL(I*pas),SNGL(densite(I)),
     &	            SNGL(HamI(I)),SNGL(MF(I))
      ENDDO
	
      end

c-----------------------------------------------------
c     COMPLEX Hamiltonian for propagation            -
c-----------------------------------------------------
      subroutine create_hprop

      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      INTEGER    I

! only diagonal part is computed at once. 
! off-diagonal part directly implemented in propagation
      do I = 0,noeud
         HamI(I) = pot_one(I)+MF(I) 
      enddo


      end


      subroutine evol_LF(pas_temps,wave_sp_1,wave_sp_2)	

      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      REAL*8        pas_temps
      complex*16    wave_sp_1(0:Noeud,Nmax_sp)
      complex*16    wave_sp_2(0:Noeud,Nmax_sp)
c coefficient pour l'application de l'hamiltonien
      COMPLEX*16    ci
      integer       I,J

      ci   = DCMPLX(0.D0,-pas_temps/HBAR)
      do J = 1,nmax_sp
	     wave_sp_2(0,J)=        wave_sp_1(0,J)
     &      +ci*(2.*H1+HamI(0))*wave_sp_1(0,J)
     &      -ci*H1             *wave_sp_1(1,J)
	 
            do I = 1,noeud-1
	           wave_sp_2(I,J)=        wave_sp_1(I,J) 
     &            +ci*(2.*H1+HamI(I))*wave_sp_1(I,J)
     &            -ci*H1             *wave_sp_1(I+1,J)
     &            -ci*H1             *wave_sp_1(I-1,J)   
            ENDDO
	        wave_sp_2(noeud,J)=         wave_sp_1(noeud,J)
     &          +ci*(2.*H1+HamI(noeud))*wave_sp_1(noeud,J)
     &          -ci*H1                 *wave_sp_1(noeud-1,J) 
	  ENDDO

      end
c----------------------------------------------------
c Alternative (unitary) evolution/Crank-Nicholson ---
c----------------------------------------------------
      
      subroutine evol_CN(pas_temps,wave_sp_1,wave_sp_2)	 ! does not conserve well ortho

      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      REAL*8        pas_temps,DU
      complex*16    wave_sp_1(0:Noeud,Nmax_sp)
      complex*16    wave_sp_2(0:Noeud,Nmax_sp)
      integer       K,J,I
      COMPLEX*16    AX(0:Noeud),BX(0:Noeud),CX(0:Noeud)
      COMPLEX*16    AX2(0:Noeud),BX2(0:Noeud),CX2(0:Noeud)
      
      COMPLEX*16    DX(0:Noeud) 
      
      DU = pas_temps/HBAR
c----------------------------------
c     compute  (1-iHam*DT/2/hbar)*wave
c----------------------------------            
      do k=0,noeud
         AX(k)=CMPLX(0.D0,0.5D0*H1*DU)         
         BX(k)=CMPLX(1.D0,-DU*(H1+0.5D0*HamI(k)))
         CX(k)=AX(k)
c-----         
         AX2(k)=CMPLX(0.D0,-0.5D0*H1*DU)         
         BX2(k)=CMPLX(1.D0,+DU*(H1+0.5D0*HamI(k)))
         CX2(k)=AX2(k)
      enddo
      
      do I=1,Nmax_sp  
        
         DX(0)=BX(0)*wave_sp_1(0,I)
     &        +CX(0)*wave_sp_1(1,I)
         do k=1,(noeud-1)
           DX(k)=AX(k)*wave_sp_1(k-1,I)+BX(k)*wave_sp_1(k,I)
     &          +CX(k)*wave_sp_1(k+1,I)
         enddo
         DX(noeud)=AX(noeud)*wave_sp_1(noeud-1,I)
     &            +BX(noeud)*wave_sp_1(noeud,I)

         
c*******************************************************
c  compute (1+iH*DTV/(2*HBAR))*Y=X --Gauss method
c*******************************************************
         CX2(0)=CX2(0)/BX2(0)
         DX(0)=DX(0)/BX2(0)
         do k=1,noeud
            CX2(k)=CX2(k)/(BX2(k)-AX2(k)*CX2(k-1))
            DX(k)=(DX(k)-AX2(k)*DX(k-1))/(BX2(k)-AX2(k)*CX2(k-1))
         enddo
         wave_sp_2(noeud,I) = DX(noeud)
         do j=1,noeud                       
             k=noeud-j
             wave_sp_2(k,I)=DX(k)-CX2(k)*wave_sp_2(k+1,I)
         enddo          

       enddo   ! end of loops on levels 
                  
       end      

c-------------------------------------
c uses the RK2 algorithm to solve the problem 
c-------------------------------------      
      subroutine evol_RK2(dt) 
            
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'
      complex*16   k1(0:Noeud,Nmax_sp)
      complex*16   k2(0:Noeud,Nmax_sp)
      complex*16   wave_sp_1(0:Noeud,Nmax_sp)
      complex*16   ci
      REAL*8       dt
      INTEGER      I,J
       
      ci   = DCMPLX(0.D0,-1./HBAR)
      
c-----mf at time t

      call dens_one(wave_sp)
	  call MEAN_FIELD
c-------compute k1	           
      DO J = 1,Nmax_sp
	     k1(0,J)=ci*(2.*H1+pot_one(0)+MF(0))*wave_sp(0,J)
     &          -ci*H1                      *wave_sp(1,J)
	 
         DO I = 1,noeud-1
	        k1(I,J)=
     &        +ci*(2.*H1+pot_one(I)+MF(I))*wave_sp(I,J)
     &        -ci*H1             *wave_sp(I+1,J)
     &        -ci*H1             *wave_sp(I-1,J)   
         ENDDO
	     k1(noeud,J)=
     &          +ci*(2.*H1+pot_one(noeud)+MF(noeud))*wave_sp(noeud,J)
     &          -ci*H1                 *wave_sp(noeud-1,J) 
         DO I=0,noeud
            wave_sp_1(I,J) = wave_sp(I,J)+dt*k1(I,J)
         ENDDO 
	  ENDDO      
c------new mean-field at t+dt	  
      call dens_one(wave_sp_1)
	  call MEAN_FIELD
c-------compute k2	           
      DO J = 1,Nmax_sp
	     k2(0,J)=ci*(2.*H1+pot_one(0)+MF(0))*wave_sp_1(0,J)
     &          -ci*H1                      *wave_sp_1(1,J)
	 
         DO I = 1,noeud-1
	        k2(I,J)=
     &        +ci*(2.*H1+pot_one(I)+MF(I))*wave_sp_1(I,J)
     &        -ci*H1             *wave_sp_1(I+1,J)
     &        -ci*H1             *wave_sp_1(I-1,J)   
         ENDDO
	     k2(noeud,J)=
     &          +ci*(2.*H1+pot_one(noeud)+MF(noeud))*wave_sp_1(noeud,J)
     &          -ci*H1                 *wave_sp_1(noeud-1,J)          
c
c----final wf at t+dt
         DO I=0,noeud
            wave_sp(I,J) = wave_sp(I,J)+dt/2.*(k1(I,J)+k2(I,J))
         ENDDO 
	  ENDDO      
      
      end
      
c-------------------------------------
c uses the Runge-Kutta 4 algorithm to solve the problem 
c NB: could be compacted easily--NOT DONE
c-------------------------------------      
      subroutine evol_RK4(dt) 
            
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'
      complex*16   k1(0:Noeud,Nmax_sp)
      complex*16   k2(0:Noeud,Nmax_sp)
      complex*16   k3(0:Noeud,Nmax_sp)
      complex*16   k4(0:Noeud,Nmax_sp)
      complex*16   wave_sp_1(0:Noeud,Nmax_sp)
      complex*16   ci
      REAL*8       dt
      INTEGER      I,J
       
      ci   = DCMPLX(0.D0,-1./HBAR)
      
c-----mf at time t
      call dens_one(wave_sp)
	  call MEAN_FIELD
c-------compute k1	           
      DO J = 1,Nmax_sp
	     k1(0,J)=ci*(2.*H1+pot_one(0)+MF(0))*wave_sp(0,J)
     &          -ci*H1                      *wave_sp(1,J)
	 
         DO I = 1,noeud-1
	        k1(I,J)=
     &        +ci*(2.*H1+pot_one(I)+MF(I))*wave_sp(I,J)
     &        -ci*H1             *wave_sp(I+1,J)
     &        -ci*H1             *wave_sp(I-1,J)   
         ENDDO
	     k1(noeud,J)=
     &          +ci*(2.*H1+pot_one(noeud)+MF(noeud))*wave_sp(noeud,J)
     &          -ci*H1                 *wave_sp(noeud-1,J) 
         DO I=0,noeud
            wave_sp_1(I,J) = wave_sp(I,J)+0.5*dt*k1(I,J)
         ENDDO 
	  ENDDO      
c------new mean-field at t+dt/2 with wf+1/2*dt*k1	  
      call dens_one(wave_sp_1)
	  call MEAN_FIELD 
c-------compute k2	           
      DO J = 1,Nmax_sp
	     k2(0,J)=ci*(2.*H1+pot_one(0)+MF(0))*wave_sp_1(0,J)
     &          -ci*H1                      *wave_sp_1(1,J)
	 
         DO I = 1,noeud-1
	        k2(I,J)=
     &        +ci*(2.*H1+pot_one(I)+MF(I))*wave_sp_1(I,J)
     &        -ci*H1             *wave_sp_1(I+1,J)
     &        -ci*H1             *wave_sp_1(I-1,J)   
         ENDDO
	     k2(noeud,J)=
     &          +ci*(2.*H1+pot_one(noeud)+MF(noeud))*wave_sp_1(noeud,J)
     &          -ci*H1                 *wave_sp_1(noeud-1,J) 
         DO I=0,noeud
            wave_sp_1(I,J) = wave_sp(I,J)+0.5*dt*k2(I,J)
         ENDDO 
     
	  ENDDO 
c------new mean-field at t+dt/2 with wf+1/2 dt k2	  
      call dens_one(wave_sp_1)
	  call MEAN_FIELD    
c-------compute k3	           
      DO J = 1,Nmax_sp
	     k3(0,J)=ci*(2.*H1+pot_one(0)+MF(0))*wave_sp_1(0,J)
     &          -ci*H1                      *wave_sp_1(1,J)
	 
         DO I = 1,noeud-1
	        k3(I,J)=
     &        +ci*(2.*H1+pot_one(I)+MF(I))*wave_sp_1(I,J)
     &        -ci*H1             *wave_sp_1(I+1,J)
     &        -ci*H1             *wave_sp_1(I-1,J)   
         ENDDO
	     k3(noeud,J)=
     &          +ci*(2.*H1+pot_one(noeud)+MF(noeud))*wave_sp_1(noeud,J)
     &          -ci*H1                 *wave_sp_1(noeud-1,J) 
         DO I=0,noeud
            wave_sp_1(I,J) = wave_sp(I,J)+dt*k3(I,J)
         ENDDO 
	  ENDDO 
c------new mean-field at t+dt with wf+dt k3	  
      call dens_one(wave_sp_1)
	  call MEAN_FIELD 
	    
c-------compute k4	           
      DO J = 1,Nmax_sp
	     k4(0,J)=ci*(2.*H1+pot_one(0)+MF(0))*wave_sp_1(0,J)
     &          -ci*H1                      *wave_sp_1(1,J)
	 
         DO I = 1,noeud-1
	        k4(I,J)=
     &        +ci*(2.*H1+pot_one(I)+MF(I))*wave_sp_1(I,J)
     &        -ci*H1             *wave_sp_1(I+1,J)
     &        -ci*H1             *wave_sp_1(I-1,J)   
         ENDDO
	     k4(noeud,J)=
     &          +ci*(2.*H1+pot_one(noeud)+MF(noeud))*wave_sp_1(noeud,J)
     &          -ci*H1                 *wave_sp_1(noeud-1,J) 
! final value of the WF
         DO I=0,noeud
            wave_sp(I,J) = wave_sp(I,J)
     &                   + dt/6.*(k1(I,J)+2.*k2(I,J)+2.*k3(I,J)+k4(I,J))
         ENDDO 
	  ENDDO    	
	     
      
      end

      subroutine test_ortho(wave_sp_1)	

      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'
      complex*16    wave_sp_1(0:Noeud,Nmax_sp)
      integer       I,J,K
      REAL*8        xnormr,xnormi

      do J = 1,nmax_sp
         DO I=J,nmax_sp
            xnormr = 0. 
            xnormi = 0. 
            do K=0,Noeud
               xnormr=xnormr
     &          + DREAL(CONJG(wave_sp_1(K,J))*wave_sp_1(K,I))
               xnormi=xnormi
     &          + DIMAG(CONJG(wave_sp_1(K,J))*wave_sp_1(K,I))
            enddo
            xnormr=xnormr*pas
            xnormi=xnormi*pas
            if( I.EQ. J) THEN
              if(abs(xnormr-1.) .GT. 0.0001) THEN
                print*,'pb norm-r-D',I,xnormr
              endif
              if(abs(xnormi) .GT. 0.0001) THEN
                print*,'pb norm-i-D',I,xnormi
              endif
            else
              if(abs(xnormr) .GT. 0.0001) THEN
                print*,'pb norm-r-D2',I,J,xnormr
              endif
              if(abs(xnormi) .GT. 0.0001) THEN
                print*,'pb norm-i-D2',I,J,xnormi
              endif
            endif  
         ENDDO
      ENDDO  
      
      end
c---------------------------------------------------
c Absorbing boundary 
c---------------------------------------------------      
      subroutine evol_IMAG 
      
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'
      
      INTEGER I,J 
      REAL*8  xim, expo
   
      do I=0,Noeud 
         xim = - delta_t*ABS(pot_imag(I))/HBAR
         if(ABS(xim) .LE. 83.) THEN
            expo = DEXP(xim)
         endif
         do J = 1,Nmax_sp
            wave_sp(I,J) = wave_sp(I,J)*expo
         enddo
      enddo
      end
      
      subroutine Output(time)
      
c-----------------------------------------
c  |----------|-----|-----|----------|
c  -xc     -boxdim  0   +boxdim     +xc 
c-----------------------------------------      
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'
      
      REAL*8 time  
      integer       I
      REAL*8      meanx,flucx,xevap,XX,xc
      
      xc = pas*float(I_centre)
      call dens_one(wave_sp)
 
      meanx = 0.
      flucx = 0.
      xevap = 0.
      do I = 0,noeud    
         XX = pas*float(I) -xc
         meanx = meanx + XX*densite(I)
         flucx = flucx + XX*XX*densite(I) 
c--------------evaporation
         if(XX .LE. box_dim) xevap = xevap + densite(I)           
      enddo
      meanx = meanx*pas
      flucx =dsqrt(flucx*pas - meanx**2.)
      xevap = float(N_part) -xevap * pas 
      
      write(50,*) time,meanx,flucx,xevap
      
      end   
      
      