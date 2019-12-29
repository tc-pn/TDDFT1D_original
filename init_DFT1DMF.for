!----------------------
! init_mf  -- Initialize the self-consistent mean-field  
! N.B. static is solved by direct diagonalization in a box smaller than total size

!-----------------------
! init_pot -- Initialize different potentials 
!           - external potential Real and Imaginary 
!- 			-pot_one -external one-body potential
!-Iteractive procedure to get the ground stte
!- 		     Idiag=1 -> direct diagonalisation (Not fast, to use with less than 500 mesh points)          
!- 		     Idiag=0 -> Iterative procedure without diagonalization (fast, can be used with more mesh points)
!                       Concergence for (approx.) more than 20 particles delicate 
!---------------------------------------------------------
      subroutine init_mf
                        
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

c boucle,boucle,boucle sur le temps 
      INTEGER     I,J,IT,N,N0
c hamiltonien a t<0

      REAL*8      XMU
      REAL*8      lambda_BCS
!      INTEGER     Ibcs 
      REAL*8      etot_mf
      INTEGER     Idiag
      parameter (Idiag=0) ! if Idiag-1 -> direct diag, Idiag=0 -> diag by hand (allow biggest mesh)

      REAL*8    dd0,om0,xi,yy,expo,xnorm
      REAL*8    Herm(0:Nmax_sp-1)
      REAL*8    EHF(Nmax_sp)
  
c-----------------------------------------------------
c     Some constant
c-----------------------------------------------------
      MASS=938.91897D0
      HBAR=197.32705D0 
      HB  = HBAR**2./(2.D0*MASS)
      H1 =  HB/pas**2
      H2 =  2*H1
c-----------------------------------------------------
c     initialization of the wave-packet              -
c     with single particle part of the hamiltonian   -
c-----------------------------------------------------

      XMU = 0.5 ! damping factor to slow down the convergence
!      XMU = 0.8 ! damping factor to slow down the convergence
      
!--return different potentials 
      call ini_pot  ! init all potentials 
      
! Initial guess for the potential 
      DO I=0,Noeud
         Upot(I)       = pot_ini(I)
      ENDDO
            
      If(Idiag .EQ. 1) THEN      
         call create_h
         print *,"h created"
         call diag_hamilt ! Return the lowest wave-functions
         print *,'h diagonalised' 
      ELSE IF(Idiag .EQ. 0) Then
! analytical form of 1D oscillators      
         om0 = sqrt(lambda_ini/MASS) ! Not MeV units
         dd0 = sqrt(MASS*om0/HBAR)
c------- Initialize H.O. wave-function by hand Phi = Hn*Gaus      
         DO I=0,Noeud
            xi = (pas*float(I-I_centre))*dd0 !/2.
            Herm(0) = 1. 
            Herm(1) = 2.*xi 
            yy = -xi**2./2.d0
            if(abs(yy) .LE. 83.D0) THEN
              expo = dexp(yy)
            else
              expo = 0.
            endif    
            wave_sp(I,1) = expo*Herm(0)
            wave_sp(I,2) = expo*Herm(1)
            DO N=2,Nmax_sp-1
               Herm(N) =2.D0*(xi*Herm(N-1) - float(N-1)*Herm(N-2)) 
               wave_sp(I,N+1) = expo*Herm(N)
            ENDDO 
         ENDDO             
c-----------------
c normalization 
c-----------------         
         DO N=1,Nmax_sp
            xnorm = 0.         
            DO I=0,Noeud
               xnorm = xnorm+ABS(wave_sp(I,N))**2. 
            ENDDO
            xnorm= 1./sqrt(xnorm*pas)                   
            DO I=0,Noeud
               wave_sp(I,N) = wave_sp(I,N)*xnorm
            ENDDO
         ENDDO
      Endif      

! initially no pairing

      If(Ibcs .eq. 1) then
!         lambda_BCS = 0.5*(Ener( N_part/2)+Ener(N_part/2+1))
!         call BCS_init(lambda_BCS)
      else
        do I=1,N_part/2
           ni_sp(I) = 1.
        enddo
      endif   

      DO IT = 1, Iimag

	   
	    call dens_one(wave_sp)
        call MEAN_FIELD 
        
        if(IT .EQ. 1) then
          do J = 0,Noeud
             Upot(J)=(MF(J) + pot_one(J))
          enddo
c          stop
        else  
          do J = 0,Noeud
             Upot(J)=(1.-XMU)*(MF(J) + pot_one(J)) + XMU*Upot(J)
          enddo
        endif
c--------------------------------
c direct diagonalization
c--------------------------------

        
        IF(Idiag .EQ. 1) THEN
           call create_h

           call diag_hamilt
        
        ELSE IF(Idiag .EQ. 0) Then
c--------------------------------
c GAUSS INVERSION METHOD
c--------------------------------      
            call GAUSS_wf(IT,EHF)   
            call ORDERING(EHF)   
            call ORTHO          
        ENDIF
        
        IF(MOD(IT-1,10) .EQ. 0) THEN
    
           print *,'-------------------------------'
           print *,'----IT=',IT,'----'
           print *,'-------------------------------'    
           DO I=1,Nmax_sp
             if(Ener(I) .LT. 0) then
               print*,'E(',I,')=', Ener(I)
             endif  
           ENDDO
         ENDIF
         
c        lambda_BCS = 0.5*( Ener( N_part/2) + Ener( N_part/2+1) )
c        DO I=1,Nmax_sp
c  	  Delta(I) = Delta_0
c        ENDDO
c        IF(Ibcs .eq. 1) then
c           call BCS_init(lambda_BCS)
c        ELSE
c           call interaction_orbital ! direct interaction -> V_inter(I,J)
c           call compute_emf(etot_mf)
c        ENDIF
      ENDDO ! end loops on imaginary time
      
      call dens_one(wave_sp)
      
      open(unit=23,file='denspotini.dat',status='unknown')
      
      DO I=0,Noeud
	     write(23,*) SNGL(I*pas),SNGL(densite(I)),
     &	             SNGL(Upot(I)),SNGL(MF(I))
      ENDDO
      close(23)
      end
!-----------------------------------------      
! initialisation of all potential 
! pot 	   --> external potential for initialisation
! pot_cplx --> Wood-Saxon + imag pot 
!-----------------------------------------      
      subroutine ini_pot 
      
      IMPLICIT NONE 
      
      INCLUDE 'TDDFT1D.COMMON'

      integer     I,J
c constante pour le potentiel imaginaire harmonique
      real*8      Value
      real*8      a_imag,b_imag,c_imag
! external potential for initialisation 
     
      I_centre = (Noeud+1)/2        
      DO I = 0,Noeud
	      pot_ini(I) = 0.5*lambda_ini*(pas*float(I_centre-I))**2.  ! harmonic potential case
	      pot_one(I) = 0.5*lambda    *(pas*float(I_centre-I))**2.    ! harmonic potential case
     &          - V_0/(1+exp((pas*abs(I_centre-I)-r_0)/a))
      enddo 

      DO I=1,Nmax_sp
	      Delta(I) = 0! Delta_0
      ENDDO
      
      DO I = 0,Noeud
  	      write(26,*)I*pas,pot_one(I),pot_ini(I)
      enddo
! external potential for propagation 
! for the moment Wood-saxon+Imaginary external potential     

      print*,'rangeimag=',rangeimag,vimag
      
      c_imag = vimag
      a_imag = c_imag/(rangeimag)**2
      b_imag = -2.D0*pas*a_imag*rangeimag
      a_imag = a_imag*pas**2

c	print*,"potentiel imaginaire",a_imag,b_imag,c_imag
      rangeimag_I = int(rangeimag/pas)
      
      DO I = 0,rangeimag_I
	    pot_imag(I) = vimag*(float(rangeimag_I-I))/float(rangeimag_I)
      enddo
      DO I = rangeimag_I,Noeud-rangeimag_I
	     pot_imag(I) = 0.D0
      enddo
      DO I = Noeud-rangeimag_I,noeud
        pot_imag(I)= 
     &      vimag*(float(-Noeud+rangeimag_I+I))/float(rangeimag_I)
      enddo
! wood saxon + complex part 
      DO I = 0,Noeud
	     pot_cplx(I) = 
     &   dcmplx(- V_0/(1+exp((pas*abs(I_centre-I)-r_0)/a )),pot_imag(I))
      enddo
    
      open(unit=25,file='pot_ext.dat',status='unknown')
      DO I=0,Noeud
	    write(25,*)I*pas,real(pot_cplx(I)),imag(pot_cplx(I))
      ENDDO
      close(25)
       
      end 
    
c-----------------------------------------------------
c     Create the hamiltonion from Upot               -
c-----------------------------------------------------
      subroutine create_h

      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'
      INTEGER    I,J

      do J = 0,noeud
         do I = 0,noeud
            Hamilt(i,j)=0.D0
         enddo
      enddo

      do i = 0,noeud
         Hamilt(i,i) = 2.D0*H1 + Upot(i)
      enddo
      do i = 0,noeud-1
         HAMILT(i,i+1) = -H1
      enddo
      do i = 1,noeud
         HAMILT(i,i-1) = -H1
      enddo

      end      
c-----------------------------------------------------
c     Diagonalize the Hamilt -> Ener,wave_sp         -
c-----------------------------------------------------
      subroutine diag_hamilt
                                               
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      real*8      A_diag(Noeud+1,Noeud+1)
      real*8      Vec(Noeud+1,Noeud+1),E(Noeud+1)
      integer     I,J,IFAIL

      EXTERNAL F02ABF

      IFAIL = 1

      DO J = 0,Noeud
         Ener(J+1) = 0.
         E = 0.
         DO I = 0,Noeud
            A_diag(J+1,I+1) = Hamilt(J,I)
            Vec(J+1,I+1) = 0.
         ENDDO
      ENDDO

      call F02ABF(A_diag,Noeud+1,Noeud+1,Ener,Vec,Noeud+1,E,IFAIL)

      if(IFAIL .NE. 0) print*,'erreur de diagonalisation'

c------------------------------
c--------  eigenvectors -------
c------------------------------

      DO J=1,Nmax_sp
         DO I = 0,Noeud
            wave_sp(I,J) = DCMPLX(Vec(I+1,J)/sqrt(pas),0.D0)
         ENDDO
      ENDDO

      end
c-----------------------------
c Gauss inversion method
c-----------------------------

      subroutine GAUSS_wf(IT,EHF)
                                    
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'
      INTEGER       IT
      integer       I,J,K,N
      REAL*8        AX(0:Noeud),BX(0:Noeud),CX(0:Noeud)
      REAL*8        DX(0:Noeud),CY(0:Noeud),DY(0:Noeud)
      REAL*8        E,xnorm
      REAL*8        EHF(Nmax_sp)
      REAL*8        WAVE0(0:Noeud),SOL(0:Noeud)
      
      
c***********
c  boucle sur les niveaux
c***********   
      do N=1,Nmax_sp
         do k=0,Noeud  
           WAVE0(K) = wave_sp(K,N)
           
           AX(k)=  -H1
           BX(k)= 2.*H1+Upot(K)
           CX(k)=  -H1
         enddo
         DX(0)=BX(0)*wave0(0)+CX(1)*wave0(1)
         do k=1,Noeud-1
           DX(k)=AX(k)*wave0(K-1)+BX(k)*wave0(K)
     &          +CX(k)*wave0(K+1)
         enddo
         DX(Noeud)= AX(Noeud)*wave0(Noeud-1)
     &             +BX(Noeud)*wave0(Noeud)
         E=0.
         do k=0,Noeud
           E=E+DX(k)*wave0(K)
         enddo
         E=E*pas
!         print*,'E=',E
         if(IT.eq.1) then
               EHF(N)=E            
         endif
         Ener(N)=E  
         EHF(N) =E        
c         print*,'ener',n,EHF(n)
c***********
c  New wave-functions from (E-H)Y=X
c*********** 
         do k=0,Noeud
            AX(k)=H1
            BX(k)=EHF(N)-BX(k)
            CX(k)=H1            
            DX(k)=wave0(k)
            SOL(k)=0.
         enddo
         CY(0)=CX(0)/BX(0)
         DY(0)=DX(0)/BX(0)
         do J=0,Noeud-1
            CY(j+1)=CX(j+1)/(BX(j+1)-AX(j+1)*CY(j))
            DY(j+1)=(DX(j+1)-AX(j+1)*DY(j))/(BX(j+1)-AX(j+1)*CY(j))
         enddo
         SOL(Noeud)=DY(Noeud)
         do j=1,Noeud                       
             k=Noeud-j
            SOL(k)=DY(k)-CY(k)*SOL(k+1)
         enddo 
         do k=0,Noeud
            wave_sp(K,N)=SOL(K)
         enddo  
      enddo   ! end of loops on levels   
!      stop 
    
         DO N=1,Nmax_sp
            xnorm = 0.         
            DO I=0,Noeud
               xnorm = xnorm+ABS(wave_sp(I,N))**2. 
            ENDDO
            xnorm= 1./sqrt(xnorm*pas)                   
            DO I=0,Noeud
               wave_sp(I,N) = wave_sp(I,N)*xnorm
            ENDDO
         ENDDO
      end     
      
      subroutine dens_one(wave_sp_1)
                                   
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'
      COMPLEX*16  wave_sp_1(0:Noeud,Nmax_sp)
      
      INTEGER I,J
	  do I = 0,noeud
	     densite(I) = 0.D0
	  enddo

      do J = 1,nmax_sp
	     do I = 0,noeud
	       densite(I)=densite(I)+2.*ni_sp(j)*CDABS(wave_sp_1(i,j))**2
	     enddo
	  enddo
	   
	  end
	  
	  
	  
c--------------------------------------------
c Mean-field Hamiltonian !! skyrme like 
c--------------------------------------------	  
      subroutine MEAN_FIELD
  
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      INTEGER I,J

      DO I = 0,Noeud
         MF(I) = 2.*c2*densite(I)
     &         +(beta+1.)*cbeta*densite(I)**beta
      ENDDO

      end
c---------------------------------------------
c compute the HF+BCS energy
c---------------------------------------------

      subroutine  compute_emf(energie_tot)

      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      INTEGER     I,J
      REAL*8      energie_tot,e_c,e_p,e_mf,e_cor


!      call interaction_orbital
      
      e_c = 0.D0    ! kinetic term 
      e_p = 0.D0    ! potential-external
      e_cor = 0.D0  ! correlation part

      DO I = 1,Nmax_sp
c	energie = energie + (Ener(I)-lambda_BCS)*2*v(I)**2
c	e0 = e0 + (Ener(I)-lambda_BCS)*2*v(I)**2
c	print*,energy

! kinetic contribution 
        e_c = e_c+2.*ni_sp(I)*(2.*H1* CDABS(wave_sp(0,I))**2.D0
     &           - H1 * DCONJG(wave_sp(0,I))*wave_sp(1,I))
        DO J = 1,Noeud-1
          e_c = e_c+2.*ni_sp(I)*(2.*H1* CDABS(wave_sp(J,I))**2.D0
     &           - H1 * DCONJG(wave_sp(J,I))*
     &                      (wave_sp(J-1,I) + wave_sp(J+1,I)))
        ENDDO
        
        e_c = e_c+2.**ni_sp(I)*(2.*H1 * CDABS(wave_sp(noeud,I))**2.D0
     &             - H1 * DCONJG(wave_sp(noeud,I))*wave_sp(noeud-1,I))

! potential part--external field contribution 
        DO J = 0,Noeud
	       e_p = e_p + 2.*ni_sp(I)*pot_one(J)
     &                       *CDABS(wave_sp(J,I))**2.D0
        ENDDO
        
!  mean-field +pairing part in the canonical basis 
      enddo
      e_p = e_p*pas
      e_c = e_c*pas
      
      e_mf = 0.D0   ! mean-field energy
      DO I = 0,noeud
	     e_mf = e_mf+c2*densite(I)**2.+cbeta*densite(I)**(beta+1.)
      enddo
      e_mf = e_mf*pas 
      
        
      
	  energie_tot = e_c + e_p + e_mf 
	  
	  print*,'e_c',SNGL(e_c),'e_p',SNGL(e_p),
     &     'e_mf',SNGL(e_mf),'Etot',SNGL(energie_tot)

      end      
c---------------------------------
c re-order the sp wave-function 
c from lowest energy to highest
c---------------------------------
      
      subroutine ORDERING(EHF)
      
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      INTEGER I,J,K
      REAL*8  Enermax 
      REAL*8  EE(Nmax_sp),EHF(Nmax_sp),EHF0(Nmax_sp)
      INTEGER Nlow(Nmax_sp)
      INTEGER Istore 
      complex*16   wave_sp_1(0:Noeud,Nmax_sp)
      
      
      DO I=1,Nmax_sp 
         Nlow(I) = 0
         EE(I)   = Ener(I)
         if(Ener(I) .LE. 0) then
         endif
      ENDDO
      do I=1,Nmax_sp
         Enermax = 1000.
         Do J=1,Nmax_sp
            IF(EE(J) .LE. Enermax) THEN
               Enermax = EE(J)
               Istore  = J 
            ENDIF
         Enddo
         Nlow(I) = Istore
         EE(Istore) = 10000.
      enddo    
c-------
c reorder wave-functions      
      DO I=1,Nmax_sp
         Istore = Nlow(I)
         
         EE(I)  = Ener(Istore)
         EHF0(I)= Ener(Istore) ! EHF(Istore) 
         DO K=0,Noeud
            wave_sp_1(K,I) = wave_sp(K,Istore)
         ENDDO 
      ENDDO   
      DO I=1,Nmax_sp
         Ener(I) = EE(I)
         EHF(I)  = EHF0(I)      
         DO K=0,Noeud
            wave_sp(K,I) = wave_sp_1(K,I)
         ENDDO   
      ENDDO
! Note that the first levels have automatically ni=1      
      end
!     orthonormalization of wave-functions
      
      subroutine ORTHO
      
      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      INTEGER I,J,K
      COMPLEX*16   xover
      REAL*8       xnorm
      
      DO J=2,Nmax_sp
         DO I=1,J-1
!           compute the overlap         
            xover = 0
            DO K=0,Noeud
               xover = xover+wave_sp(K,J)*CONJG(wave_sp(K,I))
            ENDDO
            xover = xover*pas
            DO K=0,Noeud
              wave_sp(K,J) = wave_sp(K,J)- wave_sp(K,I)*xover
            ENDDO
         ENDDO
         xnorm = 0.
                
         DO K=0,Noeud
            xnorm = xnorm+ABS(wave_sp(K,J))**2. 
         ENDDO
         xnorm= 1./sqrt(xnorm*pas)                   
         DO K=0,Noeud
            wave_sp(K,J) = wave_sp(K,J)*xnorm
         ENDDO
!  checking overlap

         DO I=1,J-1
!           compute the overlap         
            xover = 0
            DO K=0,Noeud
               xover = xover+wave_sp(K,J)*CONJG(wave_sp(K,I))
            ENDDO
            xover = xover*pas
            DO K=0,Noeud
              wave_sp(K,J) = wave_sp(K,J)- wave_sp(K,I)*xover
            ENDDO
            if(ABS(xover) .GT. 0.000001) THEN
                print*,'problem ortho',abs(xover)           
            endif
         ENDDO
         
                 
      ENDDO
      
      end
      
      
      
      