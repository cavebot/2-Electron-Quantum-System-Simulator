!#####################################################################
!#
!#    L. A.A. Nikolopoulos
!#
!#
!#    Calculates the phase and normalization constant for the 
!#    continuum states of a 1-e atomic systems. 
!#    Normalization is done according the Burgess article (1963)  
!#     
!#    A. Burgess, PROC. PHYS. SOC. 1963, VOL. 81, pg. 442
!#   ' The determination of Phases and Amplitudes of Wave Functions'
!#
!#    out:
!               * Run once for all calculated partial waves l = 0 
!               * Gets  input data from command line (angular momentum)
!               * produces out/hnorm-l.out, out/hphase-l.out
!   
!
!######################################################################

PROGRAM HNORM

  USE PRECISION, ONLY : dpk
  USE param, ONLY:lmin, lmax
  USE param, ONLY:input
  USE units

  IMPLICIT NONE

  INTEGER NOUT, NINP, NHWF1EFILE, NFILE
  INTEGER NSTATES, NPOINTS, L
  REAL(DPK) E_THRESHOLD, ZEFF, A_POL, V_POL

  REAL(DPK), DIMENSION( :),  POINTER :: energy  ! eigenenergies
  REAL(DPK), DIMENSION(:,:), POINTER :: pl,dpl  ! radial wavefunctions
  REAL(DPK), DIMENSION( :),  POINTER :: r       ! grid points    
  CHARACTER*100 ARGV

!........................................
  INTERFACE
     SUBROUTINE READ_WF1E(L, NHWF1EFILE, PL, DPL, ENERGY, ZEFF, R, NSTATES, NPOINTS)
       USE PRECISION, ONLY : DPK
       IMPLICIT NONE
       !
       INTEGER                            :: L
       INTEGER                            :: nhwf1efile
       REAL(DPK), DIMENSION(:,:), POINTER :: pl
       REAL(DPK), DIMENSION(:,:), POINTER :: dpl 
       REAL(DPK), DIMENSION(:),   POINTER :: energy  ! eigenenergies
       REAL(DPK), INTENT(OUT)             :: zeff
       REAL(DPK), DIMENSION( :),  POINTER :: r       ! grid points    
       INTEGER,   INTENT(OUT)             :: nstates
       INTEGER,   INTENT(OUT)             :: npoints

     END SUBROUTINE READ_WF1E
     !
     SUBROUTINE NORM(L, energy, pl, r, nstates, npoints, zeff, E_THRESHOLD, a, vp)
       !
       USE PRECISION, ONLY: DPK
       USE units, ONLY: M_PI
       USE ioroutines, ONLY: hnormfile, hnormasciifile, hphaseasciifile
       USE wf_1e, ONLY: p_asymptotic, delta, phi, zeta
       !
       IMPLICIT NONE
       !
       INTEGER                                 :: L
       REAL(DPK), DIMENSION(NSTATES)           :: ENERGY
       REAL(DPK), DIMENSION(NSTATES, NPOINTS)  :: PL
       REAL(DPK), DIMENSION(NPOINTS)           :: R
       INTEGER                                 :: NSTATES
       INTEGER                                 :: NPOINTS
       REAL(DPK)                               :: ZEFF
       REAL(DPK)                               :: E_THRESHOLD
       REAL(DPK)                               :: A
       REAL(DPK)                               :: VP
     END SUBROUTINE NORM
END INTERFACE




  
!........................................

  NHWF1EFILE = 12
  NFILE      = 13
  NINP       = 14
  NOUT       = 16

!.....
      

! GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT

  CALL GETARG(1, ARGV)
  READ(ARGV,*) L

  E_THRESHOLD = 0 
  A_POL = 0 
  V_POL = 0 
  
!.....

  IF(L<0) THEN ! calculate only for L = abs(l)-1

     LMIN = ABS(L) 
     LMAX = ABS(L)

  ENDIF

!............



!......      P_nl(r_i) i = 1, 2,... npoints
 

     WRITE(*,*) "###########################################"

     ! read unnormalized wf1e data and
     !... Find Normalization & Phase  A_L(E_I), PHASE_L(E_I), I = 1, NSTATES

     CALL read_wf1e(l, nhwf1efile, pl, dpl, energy, zeff, r, nstates, npoints)
     CALL norm(l, energy, pl, r, nstates, npoints, zeff, e_threshold, a_pol, v_pol)

!..............

     DEALLOCATE( ENERGY )
     DEALLOCATE(  PL )
     DEALLOCATE( DPL )
     DEALLOCATE(  R )
  

END PROGRAM hnorm

!###########################################################
SUBROUTINE READ_WF1E(L, NHWF1EFILE, PL, DPL, ENERGY, ZEFF, R, NSTATES, NPOINTS)
  
  USE UNITS,     ONLY : DPK
  USE IOROUTINES, ONLY: HWF1EFILE

  IMPLICIT NONE

  INTEGER L
  INTEGER NHWF1EFILE
  REAL(DPK), DIMENSION(:,:), POINTER :: pl     !  P_kl(r_j)
  REAL(DPK), DIMENSION(:,:), POINTER :: dpl    !  dP_kl(r_j)/dr
  REAL(DPK), DIMENSION(:),   POINTER :: energy  !  E_i
  REAL(DPK), DIMENSION(:),   POINTER :: r       !  r_j
  INTEGER,   INTENT(OUT)             :: NSTATES
  INTEGER,   INTENT(OUT)             :: NPOINTS
  REAL(DPK), INTENT(OUT)             :: ZEFF
  !
  !  REAL(DPK), DIMENSION(:):: dr       !  dr_j

  REAL(DPK) HH       !  dr_j
  INTEGER I, J
!...............


     CALL HWF1EFILE(NHWF1EFILE, L)

     READ(NHWF1EFILE,'(I5)') NSTATES
     READ(NHWF1EFILE,'(I5,1X,G25.10)') NPOINTS, HH
     READ(NHWF1EFILE,*) ZEFF

     WRITE(*,*) '#                                      L = ', L
     WRITE(*,*) '#                                NSTATES = ', NSTATES
     WRITE(*,*) '#                                NPOINTS = ', NPOINTS
     WRITE(*,*) '#                                      H = ', HH
     WRITE(*,*) '#                                   ZEFF = ', ZEFF

     ALLOCATE( ENERGY( NSTATES ) )
     ALLOCATE(  PL(NSTATES, NPOINTS) )
     ALLOCATE( DPL(NSTATES, NPOINTS) )
     ALLOCATE(   R(NPOINTS) )

     DO I = 1, NSTATES
     
        READ(NHWF1EFILE,*) ENERGY(I) 
        WRITE(*,*) '#                                  EN(I) = ', ENERGY(I)

        DO J = 1, NPOINTS
        READ(NHWF1EFILE,*) R(J), PL(I,J)
        ENDDO
     ENDDO

     CLOSE(NHWF1EFILE)
     
   END SUBROUTINE READ_WF1E
!###########################################################

SUBROUTINE NORM(L, energy, pl, r, nstates, npoints, zeff, E_THRESHOLD, a, vp)

  USE units, ONLY: M_PI
  USE precision, ONLY: DPK
  USE ioroutines, ONLY: hnormfile, hnormasciifile, hphaseasciifile
  USE wf_1e, ONLY: p_asymptotic, delta, phi, zeta!, zeta_bound !, pscoul

  IMPLICIT NONE
  INTEGER I,J
  INTEGER NOUT, NNORMASCII, NNORM, NPHASEASCII, NRATIOASCII
  
  INTEGER NPOINTS, L, NSTATES
  REAL(DPK) E_THRESHOLD, A, VP , ZEFF

  REAL(DPK), DIMENSION(300)               :: DEL
  REAL(DPK), DIMENSION(NPOINTS)           :: P, R
  REAL(DPK), DIMENSION(NSTATES)           :: ENERGY
  REAL(DPK), DIMENSION(NSTATES, NPOINTS)  :: PL

  REAL(DPK)  K_E, DE
  REAL(DPK)  A_DE
  REAL(DPK)  G_1, PHI_1, A_1, AN_1, W_1, PHI_1_RMAX
  REAL(DPK)  G_2, PHI_2, A_2, AN_2, W_2, A_1_RMAX

  INTEGER ID1

  INTEGER  P0, P_STEP
  INTEGER  P1, POINT_1

  !............................
  
  
  NOUT        = 16
  nnormascii  = 17
  nphaseascii = 18
  NNORM       = 19
  nratioascii = 20

  !....................  start matching

  P0     = 100
  P_STEP = 1
  P1     = 5


  IF( P0 > NPOINTS.OR.P1 > NPOINTS) THEN 

     WRITE(*,*) '# ERROR IN SELECTING MATCHING POINTS '
     WRITE(*,*) '#  P_STEP <P1 < PO < NPOINTS         '
     WRITE(*,*) '#                           NPOINTS =', NPOINTS
     WRITE(*,*) '#                                PO =', P0
     WRITE(*,*) '#                                P1 =', P1
     
     STOP
     
  ENDIF

  !...................

  CALL       HNORMFILE( NNORM,       L ) 
  CALL HNORMASCIIFILE(  NNORMASCII,  L )
  CALL HPHASEASCIIFILE( NPHASEASCII, L )

  OPEN(NRATIOASCII, FILE='out/hratio.out')

  WRITE(NNORM) L, NSTATES 

     WRITE(*,*)'#                 RADIAL WAVE         L  = ', L
     WRITE(*,*)'#                 NOF STATES   NSTATES_L = ', NSTATES

  
  DO I = 1, NSTATES


     IF(( ENERGY(I) - E_THRESHOLD ) <= 0.0_dpk) THEN 

        w_2   = 0.0_dpk
        an_2  = 1.0_dpk


!........  bound states 

        w_1      = 0.0_dpk
        a_1_rmax = 0.0_dpk
        an_1     = 1.0_dpk
        phi_1_rmax = 0.0_dpk

!........

     ELSE

!
        K_E = SQRT( 2 * ( ENERGY(I) - E_THRESHOLD) )


        P = PL(I,:)

!!!        P = PL(:,I) when read as pl(r,energy)


        !#
        !#   P(r) ---> A * sqrt(1/pi*k) * sqrt( k / zeta) * sin( \phi + \delta )
        !#
        !#       w1 : short range (scattering) phase shift


        !#     method No 1              (phi_1, w_1, a_1)
        

        !#   Fitting to the last P0 points with P_STEP
        !#   for getting the short-range scattering phase shift
        !#   
        !#   w1 = del
        !#  


        

     id1 = 0
     DO  j = npoints - p0, npoints - 2,  p_step 

        id1 = id1 + 1

        del(id1) = delta(r(j-1),r(j),p(j-1),p(j),l, k_e, zeff, a, vp)

     ENDDO
            


     !     w_1 = SUM(del(1:id1))/id1

     w_1 = 0.0_dpk
     DO  j = 1, id1
        w_1 = w_1 + del(j)
     ENDDO

     w_1 = w_1/id1     !short range scattering phase shift

!...............


        point_1 = npoints - p1

        phi_1 = phi(  r( point_1 ), k_e, l, zeff, a, vp)
        a_1   = zeta( r( point_1 ), k_e, l, zeff, vp)

! 
!   An is determined such a way as  :
!
!
!                  ______
!   A_k * P_k  =  /2/pi*k) * sqrt( k / z(r) ) * sin[ phi + delta ] ===>
!
!                   ________  
!                  /    2
!              =  / -------- *   sin[ phi + delta ] ===>
!                /  pi * z(r)
!              
!                      ___________                                  
!                     /2/(pi*z(R)) * sin[ phi(R) + delta(R) ]      
!        1./A_k  =  --------------------------------------- 
!                              P_k(R)                           
!
!
!   note that    z(R) --> k
!
!
        
        g_1 = SQRT( 2.0_dpk/ ( M_PI * a_1 ) ) * SIN ( phi_1 + w_1 ) 

        an_1 = g_1 / p( point_1 )



!#     method No 2              (phi_2, w_2, a_2)


         w_2  = delta( r(point_1 - 1 ), r(point_1), p(point_1 - 1), p(point_1)    &
              &   ,l, k_e, zeff, a, vp )

         phi_2 = phi(  r( point_1 ), k_E, l, zeff, a, vp)
         a_2   = zeta( r( point_1 ), K_E, l, zeff, vp)


         g_2  = SQRT(2.0_dpk / ( M_PI*a_2) ) * SIN( phi_2 + w_2 )


         an_2 = g_2 / p(point_1) 

!.....
         IF(w_1.LT.0.0_dpk)  w_1 = w_1 + 2 * M_PI
         IF(w_2.LT.0.0_dpk)  w_2 = w_2 + 2 * M_PI


!#     method No3 : energy differences       an
         
         IF(energy(i-1).LT.E_THRESHOLD) THEN
         
            DE = 2.0_dpk * ABS( energy(i+1) - energy(i) )

         ELSE IF(I.EQ.NSTATES) THEN

            DE = 2.0_dpk * ABS( energy(i) - energy(i-1) )

         ELSE

            DE =  ABS( energy(i+1) - energy(i-1) )
         ENDIF


            a_de = SQRT(2.0_dpk / DE )
!.....
!            
         !phase at rmax

            phi_1_rmax = phi(  r( npoints), k_e, l, zeff, a, vp)
            a_1_rmax   = zeta( r( npoints), k_e, l, zeff, vp)

         ENDIF

         IF(w_1.GT.M_PI) w_1= w_1 - 2.0_dpk * M_PI

         WRITE(nnormascii, '(6G20.10)') energy(i), w_1,w_2,an_1, an_2, a_de
         WRITE(nphaseascii,'(5G20.10)') energy(i), w_1, w_1/M_PI, phi_1_rmax + w_1, a_1_rmax

         WRITE(nnorm)  i, k_e, w_1, an_1, a_1_rmax, phi_1_rmax + w_1
      ENDDO

      OPEN(30, file="out/wf.dat")
      DO J = 1, NPOINTS
         WRITE(30,*) R(J), AN_1 * P(J)
      ENDDO
      CLOSE(30)

      CLOSE(NNORM)
      CLOSE(NNORMASCII)
      CLOSE(NPHASEASCII)
      CLOSE(NRATIOASCII)


      RETURN

      
 END SUBROUTINE NORM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*EOF







