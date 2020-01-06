!------------------------------------------------------------------------------------
! This program solves the TDDFT problem of N particles in 1D
! on a mesh with 
! -a mean-field potential is of the form: EMF-c2*rho^2+cbeta*rho^(beta+1)
! -with an external initial harmonic field -lambda
! -with a WS external field - pot_one v_0,a,r_0
! -with imaginary external boundary conditions -rangeimag,vimag -- Not yet implemented
!-----------------------------------------------------------------------------------


      program mainTDDFT1D

      IMPLICIT NONE

      INCLUDE 'TDDFT1D.COMMON'

      integer     I,J
c indication du nbr de noeud pour trace dans le .in
      integer     noeud_check,Nmax_sp_check
c occupation des orbitals

      read(5,*) Ibcs
      print*,'',Ibcs
      read(5,*) pas,delta_t
      read(5,*) Ntime , NPRINT
      read(5,*) name_run
      read(5,*) lambda_ini
      read(5,*) lambda
      read(5,*) v_0,a,r_0
      read(5,*) rangeimag,vimag
      read(5,*) box_dim 
c      read(5,*) noeud_check
c      read(5,*) Nmax_sp_check
      read(5,*) N_part
      read(5,*) c2,cbeta,beta ! SKYRME like mean-field E=c2*rho^2+cbeta*rho^(beta+1)
      read(5,*) Iimag
      read(5,*) Delta_0
      
!--------------------some dimension checking-------------
!      if ( noeud_check .ne. noeud ) then
!	     print*,"erreur noeud"
!	     stop
!      endif
!      if ( Nmax_sp_check .ne. Nmax_sp ) then
!	      print*,"erreur Nmax_sp"
!	      stop
!      endif
      if ( 2*Nmax_sp .lt. N_part ) then
	     print*,"Error : Nmax_sp too small"
	     stop
      endif
c---------------------
c Initialize mean-field 
c Ibcs = 0 -> no pairing
c Ibcs = 1 -> with pairing  !! DO NOT CONVERGE YET !!
c---------------------
 
       call init_mf
        
       call evol_mf
      
      end 





