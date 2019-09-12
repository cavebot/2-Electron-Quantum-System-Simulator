!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
MODULE bs_h1e

IMPLICIT NONE
!
CONTAINS
  SUBROUTINE SOLVE_BANDED_DIAG_NAG77(L, H, B, EN, CE)
  
!  USE nag_f77_f_chapter
    USE PRECISION,  ONLY : DPK
    USE PARAM,      ONLY : KB, NDIM
    USE UTILS,      ONLY : print_mx
    !......................

    IMPLICIT NONE
    INTEGER I, J
    INTEGER N, K 
    !
  INTEGER                                           :: L
  REAL(DPK), DIMENSION( ndim, ndim),  INTENT(IN)    :: H
  REAL(DPK), DIMENSION( ndim, ndim),  INTENT(IN)    :: B
  REAL(DPK), DIMENSION( ndim, ndim),  INTENT(INOUT) :: CE
  REAL(DPK), DIMENSION( ndim),        INTENT(INOUT) :: EN
  ! DSBGV/LAPACK/BANDED CALCULATION               
  REAL(DPK), DIMENSION( kb, ndim)                   :: hb, bb     
  REAL(DPK), DIMENSION( ndim,  ndim)                :: z
  REAL(DPK), DIMENSION( 3*ndim )                    :: w
  !
  INTEGER INFO
  INTEGER NOUT
!..............................
  WRITE(*,*) '# subroutine:: solve_banded_diag_nag77 in.'
  WRITE(*,'(a60,i10)') 'ndim = ', ndim

  nout = 12

  K   = KB - 1


  !#  get the symmetric banded format

  OPEN(nout, file='out/subsolve_fxd.out')

  !  CALL print_mx(ndim,h,'h','f')
  !  CALL print_mx(ndim,b,'b','f')


  ! transform the (n-1)x(n-1) 'free' matrices to (n-2)x(n-2) 'fxd' matrices 
  ! equivalent to setting P(R) = 0 <=> or excluding the B_n spline.
  !
      hb = 0.0_dpk
      bb = 0.0_dpk
      DO j = 1, ndim
         DO i = MAX(1, j - k), j
            hb( K + 1 + I - J, J ) =   h(I, J) 
            bb( K + 1 + I - J, J ) =   b(I, J) 
         ENDDO
      ENDDO
      !      CALL print_mx(ndim,hb,'hb','f')
      !      CALL print_mx(ndim,bb,'bbb','f')

     
      ! go on and diagonalize


      !
      !  DSBGV/LAPACK :     A * C = e * B * C     
      !   
      ! where, A,B banded real symmetric
      !
      !   EV NORMALIZED AS : Z^T * C * Z = 1
      !   ev(i) <====> z(j,i)
      !

      info   = 0
      CALL DSBGV('V', 'U', ndim, k, k, hb, kb, bb, kb, en, z, ndim, w, info)

     
      IF(INFO.NE.0) THEN 
         WRITE(*,'(a60,i10)') 'error in  DSBGV/LAPACK subroutine,  ierror = ', info 
         WRITE(*,'(a60)') 'solve_banded_diag_nag77  exit.'
         STOP
      ELSE         
         DO  i = 1 , ndim
            WRITE(nout,"(I4,3X,3E20.10)") i, en(i)
         ENDDO

         ! ASSIGN Z(NB-2,NB-2) ---> C(NB-2,NB-1)

         CE = TRANSPOSE(Z)         ! z(nb-2,nb-2) --> c(nb-2,nb-1)

         !
         ! get the solutions with dP(0)/dr > 0
         !

         DO i = 1, SIZE(z,dim=1) 
            IF(ce(i,2).LT.0.0_dpk) THEN
               ce(i,:) = -ce(i,:)
            END IF
         ENDDO
         

!         CE(:, NB-1 )  = 0.0D+00
                       
!         DO I = 1, INDEX_B
!            WRITE(nout, *) '#           COEFF(I) = ', I
!            DO  J = 1, NB/2              
!               WRITE(NOUT, "(I4,1X,E25.14,2X,I4,1X,E25.14)") &
!                    & J, CE(I,J), J+NB/2, CE(I,J+NB/2)
!            ENDDO
!         ENDDO
         
      ENDIF
      
      CLOSE(NOUT)
      WRITE(*,*) '# subroutine:: solve_banded_diag_nag77 out.'     
      RETURN
    END SUBROUTINE SOLVE_BANDED_DIAG_NAG77
    !S
    !S
    !S
    SUBROUTINE ENERGY_SPECTRUM(l, en, mode)
      !
      USE param
      USE DATA, ONLY: write_v
      !
      IMPLICIT NONE
      !
      INTEGER                                   L
      REAL(DPK), DIMENSION(:), INTENT(INOUT) :: EN
      CHARACTER(LEN=*)                          MODE
      ! L/
      INTEGER I, NS

!.....................
      WRITE(*,*) '# subroutine::  energy_spectrum in. '

     
      ns = SIZE(en)


      WRITE(*,'(a60,i10)') ' nof b-splines basis ns = ', ns
      
      IF (MODE =='OUT') THEN 

         DO  I = L + 1 , NS
            
            IF(I.LE.NBS + L) THEN 

               EN(I) = - 0.5D+00 * (ZNUC/I)**2

            ELSE IF(I.LE.NS) THEN 
               
               EN(I) =   EN_1E  +    DE_1E * ( I - NBS - 1 - L) 
     
            ENDIF

         ENDDO

         CALL write_v(l, en, "en1e-")

      ELSE IF (MODE =='INOUT') THEN      ! do not store

         DO  I = 1 , NS
            
            EN(I) =   EN_1E  +    DE_1E * ( I - 1) 
            
         ENDDO

      ELSE IF (MODE == 'IN') THEN      
         CALL write_v(l, en, "en1e-")
      ENDIF

      WRITE(*,*) '# subroutine::  energy_spectrum out. '
    END SUBROUTINE ENERGY_SPECTRUM
    !S
    !S
    !S
    
  END MODULE bs_h1e


!
!
!

    !###########################################
!!!!!!!! main subroutine for free-boundary diagonalization (1)
!!%SUBROUTINE SOLVE_BANDED_LIN_NAG77(L, H, B, C0, EN, CE)
!!%  
!!%!  USE nag_f77_f_chapter
!!%  USE param, ONLY : DPK, METHOD, SPECTRUM
!!%  USE param, ONLY : en_1e, de_1e
!!%  USE param, ONLY : nbs, ncs 
!!%  USE param, ONLY : nb, kb, rmax, znuc
!!%  USE utils, ONLY : check, print_mx
!!%  USE units, ONLY : M_PI, M_PI_2
!!%
!!%  IMPLICIT NONE
!!%
!!%  INTEGER I, IN, J, INFO, INDEX_B 
!!%  INTEGER K
!!%
!!%!...........
!!%
!!%  INTEGER  NN, NS, L
!!%  
!!%  REAL(DPK), DIMENSION(:,:),    INTENT(IN) :: H,B
!!%  REAL(DPK), DIMENSION(:,:),    INTENT(IN) :: c0
!!%
!!%  REAL(DPK), DIMENSION(NB-1, 1)            :: c
!!%
!!%  REAL(DPK), DIMENSION(:),   INTENT(inout) :: en
!!%  REAL(DPK), DIMENSION(:,:), INTENT(inout) :: ce
!!%
!!%! linear equation solver 
!!%
!!%  REAL(DPK) ENERGY
!!%  REAL(DPK), DIMENSION(3*KB-2, NB-1)       :: A  ! linear
!!%  INTEGER,               DIMENSION(NB-1)   :: iw
!!%  REAL(DPK)  NORM, K_E, KR, C_NB
!!%  INTEGER NOUT
!!%  INTEGER dim, dim_1, dim_2
!!%  integer er
!!%
!!%  REAL(DPK) PSCOUL
!!%!..............................
!!%
!!%  WRITE(*,*) '# h1e::solve_banded_lin_nag77:       call.'
!!%
!!%  NOUT     = 12
!!%
!!%
!!%  K   = KB - 1             ! number of sub-diagonials (upper diagonals)
!!%  NN  = NB - 1             ! dimension of coefficient vector
!!%
!!%
!!%
!!%!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES
!!%
!!%
!!%  OPEN(NOUT, FILE='out/subsolve_lin.out')
!!%
!!%  WRITE(nout, *) '#       LINEAR SOLUTION : '
!!%  WRITE(nout, *) '#'
!!%  WRITE(nout, *) '#                          K + 1 = KB     = ', K + 1
!!%  WRITE(nout, *) '#                            NN  = NB - 2 = ', NN
!!%  WRITE(nout, *) '#'
!!%  WRITE(nout, *) "# STATE        E_FREE"
!!%
!!%
!!%!.........................................
!!%  
!!%  dim = size(en)
!!%  dim_1 = SIZE(ce,1)
!!%  dim_2 = SIZE(ce,2)
!!%
!!%  call check( dim, dim_1, er)
!!%  call check( dim_2,nb-1, er)
!!%  
!!%  
!!%  DO IN = 1, dim
!!%
!!%
!!%     DO J = 1, NN
!!%        DO I = MAX(1, J - K), MIN(NN, J + K)
!!%           A( 2*K + I - J + 1, J) = H(I,J) - EN(IN) * B(I,J)
!!%        ENDDO
!!%     ENDDO
!!%
!!%!..........................................
!!%
!!%     WRITE(nout, *)  '# linear:: input energy was en(i) = ', IN, EN(IN)
!!%
!!%     ! FACTORIZE FIRST
!!%
!!%     IW = 0
!!%
!!%
!!%     CALL F07BDF(NN, NN, K, K, A, 3*KB-2, IW, INFO) 
!!%
!!%     IF(INFO.NE.0) THEN 
!!%        CALL X04CEF(NN, NN, K, 2*K, A, NN,'Details of factorization', INFO)
!!%     ELSE
!!%             
!!%        c = c0
!!%        !             C        = 0.0D+00
!!%        !             C(NN, 1) = 1.0D+00
!!%
!!%        CALL F07BEF('N', NN, K, K, 1, A, 3*KB-2, IW, C, NN, INFO)
!!%     ENDIF
!!%
!!%
!!%     IF(INFO.NE.0) THEN 
!!%        WRITE(nout,*) '# solve_banded_lin_nag77: error in factorization subroutine f07bef& 
!!%                      &  info = ', info                     
!!%        STOP
!!%     ELSE
!!%        WRITE(nout, *) '#   INPUT     ENERGY      E = ', EN(IN)
!!%                          
!!%        NORM = DOT_PRODUCT( C(:,1), MATMUL( B, C(:, 1) ))
!!%
!!%        CE(IN,:) = C(:,1) /SQRT(NORM)
!!%     ENDIF
!!%     
!!%     WRITE(NOUT, *) IN, EN(IN) 
!!%
!!%!?     IF(IN.EQ.30.AND.L.EQ.0) THEN
!!%!        DO  J = 1, (NB-1)/2
!!%!           WRITE(NOUT, "(I4,1X,E25.14,2X,I4,1X,E25.14)") J, CE(IN,J), J+NB/2, CE(IN, J+NB/2)
!!%!        ENDDO
!!%!     ENDIF
!!%     
!!%  ENDDO
!!%
!!%
!!%  CLOSE(NOUT)
!!%
!!%  RETURN
!!%       
!!%END SUBROUTINE SOLVE_BANDED_LIN_NAG77
!!%!EOF
!###############################################          DALGARNO-LEWIS METHOD 
! main subroutine for free-boundary diagonalization (2)
!!%SUBROUTINE SOLVE_BANDED_LIN_NAG77_DL(L, H, B, E_C, COEFF_C, E_B, COEFF_B, CHOICE)
!!%  
!!%!  USE nag_f77_f_chapter
!!%  USE param, ONLY : DPK, METHOD
!!%  USE param, ONLY : en_1e, de_1e
!!%  USE param, ONLY : ncs, nbs
!!%  USE param, ONLY : nb, kb, rmax, znuc
!!%  USE ioroutines, ONLY: d1efile
!!%  USE units, ONLY : M_PI, M_PI_2
!!%
!!%  IMPLICIT NONE
!!%
!!%  INTEGER I, IN, J, INFO, INDEX_B 
!!%  INTEGER K
!!%
!!%!...........
!!%
!!%  INTEGER  NN, NSTATES, L
!!%  
!!%  REAL(DPK), DIMENSION(NB-1, NB-1), INTENT(IN) :: H,B
!!%  REAL(DPK), DIMENSION(NB-1, 1):: c
!!%  REAL(DPK), DIMENSION(NCS), INTENT(OUT) :: e_c
!!%  REAL(DPK), DIMENSION(NCS, NB-1), INTENT(OUT):: coeff_c
!!%  REAL(DPK), ALLOCATABLE,DIMENSION(:) :: E_FREE
!!%!...  for bound states
!!%  REAL(DPK), DIMENSION(NBS), INTENT(OUT) :: e_b
!!%  REAL(DPK), DIMENSION(NBS, NB-1), INTENT(OUT):: coeff_b
!!%  REAL(DPK), ALLOCATABLE,DIMENSION(:) :: E_HYD
!!%!.... for bound states
!!%  REAL(DPK), DIMENSION(NB-2) :: ER
!!%
!!%
!!%
!!%! linear equation solver 
!!%
!!%  REAL(DPK) ENERGY
!!%  REAL(DPK), DIMENSION(3*KB-2, NB-1) :: hh  ! linear
!!%
!!%  REAL(DPK)  NORM, K_E, KR, C_NB
!!%  INTEGER, DIMENSION(NB-1)   :: iw
!!%  INTEGER NOUT, ND1EFILE
!!%  INTEGER CHOICE    ! 0 --> BOUND + CONTINUUM , 1 == CONTINUUM SPECTRUM
!!%
!!%  REAL(DPK) PSCOUL
!!%!..............................
!!%
!!%  WRITE(*,*) '# h1e::solve_banded_lin_nag77:       call.'
!!%
!!%  ND1EFILE = 11
!!%  NOUT     = 12
!!%
!!%
!!%!..............................
!!%
!!%
!!%  K   = KB - 1             ! number of sub-diagonials (upper diagonals)
!!%  NN  = NB - 1             ! dimension of coefficient vector
!!%
!!%
!!%!............................
!!%
!!%
!!%!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES
!!%
!!%
!!%  OPEN(NOUT, FILE='out/subsolve_lin.out')
!!%
!!%  WRITE(nout, *) '#       LINEAR SOLUTION : '
!!%  WRITE(nout, *) '#'
!!%  WRITE(nout, *) '#                          K + 1 = KB     = ', K + 1
!!%  WRITE(nout, *) '#                            NN  = NB - 2 = ', NN
!!%  WRITE(nout, *) '#'
!!%  WRITE(nout, *) "# STATE        E_FREE"
!!%
!!%
!!%!................................
!!%
!!%
!!%
!!%
!!%
!!%!......
!!%
!!%  CALL D1EFILE(ND1EFILE, L) 
!!%
!!%  DO I = 1, NB-2
!!%
!!%     READ(ND1EFILE) ER(I)
!!%  ENDDO
!!%
!!%  READ(ND1EFILE) INDEX_B
!!%      
!!%  CLOSE(ND1EFILE)
!!%
!!%
!!%!.................................
!!%
!!%
!!%  IF(DISCRETE.EQ.'FX'.AND.NCS.GT.( NB - 2 - INDEX_B ) ) THEN
!!%
!!%     WRITE(*,*) '#                N_CONTINUUM <= NB-2'
!!%     WRITE(*,*) '#                            N_C  = ', NCS
!!%     WRITE(*,*) '#                         N_C_FXD = ', NB - 2 - INDEX_B
!!%     WRITE(*,*) '#'
!!%
!!%     NCS = NB - 2 - INDEX_B
!!%!     STOP
!!%  END IF
!!%
!!%
!!%  ALLOCATE( E_FREE(NCS) ) 
!!%
!!%
!!%      
!!%  IF(CHOICE.EQ.0) THEN 
!!%
!!%     ALLOCATE( E_HYD( NBS + L ) ) 
!!% 
!!%     DO  I = L + 1 , NBS + L
!!%            
!!%        E_HYD(I) = - 0.5D+00 * (ZNUC/I)**2
!!%     ENDDO
!!%
!!%
!!%!         LINEAR MATRICES 
!!%!         BOUND  STATES
!!%!
!!%
!!%     DO IN = 1, NBS
!!%
!!%        IF(DISCRETE.EQ.'FX') THEN
!!%             
!!%           E_B(IN) = ER(IN + INDEX_B)
!!%        ELSE
!!%             
!!%           E_B(IN) = E_HYD(IN +L)
!!%        ENDIF
!!%
!!%
!!%!......... make the matrix
!!%
!!%
!!%        DO J = 1, NN
!!%           DO I = MAX(1, J - K), MIN(NN, J + K)
!!%              
!!%              HH( 2*K + I - J + 1, J) = H(I,J) - E_B(IN) * B(I,J)
!!%                  
!!%           ENDDO
!!%        ENDDO
!!%
!!%!..........................................
!!%
!!%        WRITE(NOUT, *)  '#   INPUT ENERGY  E = ', IN, E_B(IN)
!!%     
!!%
!!%
!!%        ! FACTORIZE FIRST
!!%
!!%        IW = 0.0D+00
!!%
!!%
!!%        CALL F07BDF(NN, NN, K, K, HH, 3*KB-2, IW, INFO) 
!!%
!!%        
!!%        IF(INFO.NE.0) THEN 
!!%             
!!%           WRITE(nout,*) '# ERROR IN FACTORIZATION SUBROUTINE '
!!%           WRITE(nout,*) '#'
!!%           WRITE(nout,*) '#                     INFO = ', INFO
!!%           WRITE(nout,*) '#'
!!%               
!!%
!!%           CALL X04CEF(NN, NN, K, 2*K, HH, NN,'Details of factorization', INFO)
!!%
!!%        ELSE
!!%               
!!%             
!!%           C        = 0.0D+00
!!%           C(NN, 1) = 1.0D+00
!!%
!!%           CALL F07BEF('N', NN, K, K, 1, HH, 3*KB-2, IW, C, NN, INFO)
!!%             
!!%        ENDIF
!!%
!!%
!!%        IF(INFO.NE.0) THEN 
!!%
!!%           WRITE(nout,*) '#         ERROR IN   F07BEF          '
!!%           WRITE(nout,*) '#'
!!%           WRITE(nout,*) '#                     INFO = ', INFO
!!%           WRITE(nout,*) '#'
!!%               
!!%           STOP
!!%               
!!%        ELSE
!!%
!!%           WRITE(nout, *) '#   INPUT     ENERGY      E = ', E_B(IN)
!!%           WRITE(nout, *) '#                     NORM  = ', NORM
!!%           WRITE(nout, *) '#                COEFF_B(I) = '
!!%        ENDIF
!!%          
!!%        WRITE(NOUT, *) IN, E_B(IN) 
!!%
!!%        IF(IN.EQ.2.AND.L.EQ.0) THEN
!!%           DO  J = 1, (NB-1)/2
!!%             
!!%              WRITE(NOUT, "(I4,1X,E25.14,2X,I4,1X,E25.14)") J,      COEFF_B(IN, J )     &
!!%                   & ,J+NB/2, COEFF_B(IN, J+NB/2)
!!%           ENDDO
!!%        ENDIF
!!%            
!!%     ENDDO
!!%       
!!%  ENDIF
!!%
!!%!  CONTINUUM STATES
!!%!.........................................
!!%
!!%
!!%  DO IN = 1, NCS
!!%
!!%
!!%     IF(DISCRETE.EQ.'FX') THEN
!!%             
!!%        E_C(IN) = ER(IN + INDEX_B)
!!%
!!%     ELSE
!!%
!!%        E_C(IN) =   EN_1E  +    DE_1E * ( IN - 1 ) 
!!%
!!%     ENDIF
!!%
!!%
!!%
!!%
!!%!.........
!!%
!!%!!$          K_E = SQRT( 2 * E_C(IN))
!!%!!$          KR = K_E * RMAX 
!!%!!$
!!%!!$
!!%!!$          C_NB =  SQRT( 2.0D+00/(M_PI*K_E) ) &
!!%!!$               & * SIN(  ( K_E + ZNUC * LOG(2 * KR) / KR )*RMAX - M_PI_2 * L & 
!!%!!$               & + PSCOUL(L, ZNUC, K_E) )
!!%!!$
!!%!!$          IF(L.EQ.1) THEN
!!%!!$
!!%!!$             WRITE(15,*) K_E, C_NB, IN
!!%!!$          ENDIF
!!%!!$
!!%!!$!.........
!!%
!!%
!!%     DO J = 1, NN
!!%        DO I = MAX(1, J - K), MIN(NN, J + K)
!!%
!!%           HH( 2*K + I - J + 1, J) = H(I,J) - E_C(IN) * B(I,J)
!!%
!!%        ENDDO
!!%     ENDDO
!!%
!!%!..........................................
!!%
!!%     WRITE(NOUT, *)  '#   INPUT ENERGY  E = ', IN, E_C(IN)
!!%     
!!%
!!%
!!%     ! FACTORIZE FIRST
!!%
!!%     IW = 0.0D+00
!!%
!!%
!!%          CALL F07BDF(NN, NN, K, K, HH, 3*KB-2, IW, INFO) 
!!%
!!%
!!%          IF(INFO.NE.0) THEN 
!!%
!!%             WRITE(nout,*) '# ERROR IN FACTORIZATION SUBROUTINE '
!!%             WRITE(nout,*) '#'
!!%             WRITE(nout,*) '#                     INFO = ', INFO
!!%             WRITE(nout,*) '#'
!!%
!!%
!!%             CALL X04CEF(NN, NN, K, 2*K, HH, NN,'Details of factorization', INFO)
!!%
!!%          ELSE
!!%
!!%             
!!%             C        = 0.0D+00
!!%             C(NN, 1) = 1.0D+00
!!%
!!%             CALL F07BEF('N', NN, K, K, 1, HH, 3*KB-2, IW, C, NN,INFO)
!!%
!!%          ENDIF
!!%
!!%
!!%          IF(INFO.NE.0) THEN 
!!%             
!!%             WRITE(nout,*) '#         ERROR IN   F07BEF          '
!!%             WRITE(nout,*) '#'
!!%             WRITE(nout,*) '#                     INFO = ', INFO
!!%             WRITE(nout,*) '#'
!!%
!!%             STOP
!!%
!!%          ELSE
!!%
!!%             
!!%             WRITE(nout, *) '#   INPUT     ENERGY      E = ', E_C(IN)
!!%             WRITE(nout, *) '#'
!!%             
!!%                          
!!%             NORM = DOT_PRODUCT( C(:,1), MATMUL( B, C(:, 1) ))
!!%
!!%
!!%
!!%             COEFF_C(IN,:) = C(:,1) /SQRT(NORM) 
!!%
!!%             !             WRITE(*,*) E_C(IN), K_E, KR, C_NB,  2 * ZNUC * LOG(2 *KR)/KR
!!%             !             COEFF_C(IN,:) = COEFF_C(IN,:) *  C_NB / COEFF_C(IN,NB-1)
!!%             !             WRITE(*,*) IN, SQRT(2*E_C(IN)), C_NB,  C(NB-1,1), COEFF_C(IN,NB-1)
!!%             !             STOP
!!%
!!%
!!%             WRITE(nout, *) '#           NORM     = ', NORM
!!%             WRITE(nout, *) '#           COEFF_C(I) = '
!!%             
!!%             
!!%          ENDIF
!!%
!!%          WRITE(NOUT, *) IN, E_C(IN) 
!!%
!!%          IF(IN.EQ.30.AND.L.EQ.0) THEN
!!%             DO  J = 1, (NB-1)/2
!!%               
!!%                WRITE(NOUT, "(I4,1X,E25.14,2X,I4,1X,E25.14)") J, COEFF_C(IN,J), J+NB/2, COEFF_C(IN, J+NB/2)
!!%             ENDDO
!!%          ENDIF
!!%
!!%       ENDDO
!!%
!!%
!!%       CLOSE(NOUT)
!!%
!!%       RETURN
!!%       
!!%     END SUBROUTINE SOLVE_BANDED_LIN_NAG77_DL
!!%!###################################################### 
!!%SUBROUTINE SOLVE_BBANDED_INVERSE_NAG77(HB, BBB, ENERGY_BOUND, C)
!!%  
!!%!  USE nag_f77_f_chapter
!!%  USE param, ONLY : dpk
!!%  USE param, ONLY : ncs, nbs
!!%  USE param, ONLY : nb, kb, rmax, znuc
!!%  IMPLICIT NONE
!!%
!!%  INTEGER I,II, J, INFO 
!!%  INTEGER N, K, NSTATES 
!!%
!!%  INTEGER  NN
!!%
!!%  REAL(DPK), DIMENSION(NB-1,NB-1), INTENT(IN) :: HB, BBB
!!%  REAL(DPK), DIMENSION(NBS,NB-1),INTENT(inout) :: c
!!%
!!%! diagonalization (fixed boundary conditions) 
!!%
!!%  REAL(DPK), DIMENSION(KB, NB-2) :: a, b         !  fxd diagonalization 
!!%  REAL(DPK), DIMENSION(NB-2) :: er
!!%
!!%  REAL(DPK), DIMENSION(2*KB-1, NB-1) :: aa, bb    ! inverse  iteration
!!%  REAL(DPK), DIMENSION(NB-2), INTENT(OUT):: energy_bound
!!%
!!%  REAL(DPK), DIMENSION(30)   :: d
!!%  REAL(DPK), DIMENSION( (NB-1)*(KB+1) )   :: w
!!%
!!%  REAL(DPK) REL_ERROR, NORM, ENERGY
!!%  INTEGER, DIMENSION(NB-1)   :: iw
!!%  INTEGER LW
!!%  INTEGER NOUT
!!%
!!%  REAL(DPK), ALLOCATABLE,DIMENSION(:) :: E_FREE, E_HYD
!!%!..............................
!!%
!!%  WRITE(*,*) '# h1e::solve_bbanded_inverse_nag77:       call.'
!!%
!!%  nout = 12
!!%
!!%  
!!%!..............................
!!%
!!%
!!%      K   = KB - 1
!!%      N   = NB - 2
!!%      NN  = NB - 1
!!%
!!%
!!%      ALLOCATE(E_FREE(N)) 
!!%      ALLOCATE(E_HYD(N)) 
!!%
!!%
!!%!............................
!!%
!!%!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES
!!%
!!%
!!%      OPEN(NOUT, FILE='out/subsolve_inv.out')
!!%
!!%      WRITE(nout, *) '#        DIMENSION OF BANDED MATRICES  (K+1) X N '
!!%      WRITE(nout, *) '#'
!!%      WRITE(nout, *) '#                                K + 1 = KB     = ', K + 1
!!%      WRITE(nout, *) '#                                   N  = NB - 2 = ', N
!!%      WRITE(nout, *) '#                                  NN  = NB - 2 = ', NN
!!%      WRITE(nout, *) '#'
!!%
!!%
!!%      WRITE(nout, *) '#  FIXED BOUNDARY DIAGONALIZATION : '
!!%
!!%
!!%      A = 0.0D+00
!!%      B = 0.0D+00
!!%
!!%      DO J = 1, N
!!%         DO I = MAX(1, J - K), J
!!%
!!%            A( K + 1 + I - J, J ) =   HB(I, J) 
!!%            B( K + 1 + I - J, J ) =   BBB(I, J) 
!!%
!!%         ENDDO
!!%      ENDDO
!!%
!!%
!!%!     
!!%      !   SOLVES THE BANDED FORM OF THE SCHONDIGER HAMILTONIAN EQNS  
!!%      !   DIAGONALIZATION IS PERFORMED
!!%      !   
!!%      !                  ONLY EIGENVALUES
!!%      !
!!%
!!%
!!%      INFO   = 0
!!%      LW     = MAX( N, ( 3*KB + KB ) * ( KB + KB + 1) )
!!%
!!%      CALL F02FHF(N, K, A, KB, K, B, KB, ER, W, LW, INFO)
!!%      
!!%      IF(INFO.NE.0) THEN 
!!%
!!%        WRITE(*,*) '# ERROR IN  F02FHF/NAG_77/M20  SUBROUTINE '
!!%        WRITE(*,*) '#'
!!%        WRITE(*,*) '#                     INFO = ', INFO
!!%        WRITE(*,*) '#'
!!%
!!%        STOP
!!%
!!%     ENDIF
!!%
!!%     
!!%
!!%      WRITE(nout,*) " STATE        E_DISCRETE           E_FREE             E_HYDROGEN"
!!% 
!!%       DO  I = 1, N
!!% 
!!%          E_FREE(I) = 0.5D+00* ( I * 3.141592654D+00/ RMAX )**2
!!%          E_HYD(I) = - 0.5D+00 * (ZNUC/I)**2
!!%
!!%          IF(ER(I)<0.0D+00) THEN
!!%
!!%             WRITE(nout,"(I4,3X,3E20.10)") I, ER(I), E_HYD(I), E_HYD(I)
!!%          ELSE
!!%
!!%             WRITE(nout,"(I4,3X,3E20.10)") I, ER(I), E_FREE(I), E_HYD(I)
!!%          ENDIF
!!%
!!%       ENDDO
!!%
!!%
!!%       WHERE(ER < 0.0D+00)  ENERGY_BOUND = ER
!!%
!!%
!!%!.......................    INVERSE ITERATION
!!%!!$!
!!%!   NOW WE GET THE ENERGY EIGENVALUES FROM THE PREVIOUS DIAGONALIZATION 
!!%!   AND WE ASK THE SOLUTION AT THE CONTINUUM ( E(I) > 0 ) WITH 
!!%!   INVERSE ITERATION METHOD
!!%!
!!%!
!!%
!!%
!!%
!!%!.......................
!!%
!!%       DO I = 1, NBS 
!!%
!!%
!!%       DO J = 1, NN
!!%          DO II = MAX(1, J - K), MIN(NN, J + K)
!!%
!!%             AA( K + 1 + II - J, J ) =   HB(J, II) 
!!%             BB( K + 1 + II - J, J ) =   BBB(J, II) 
!!%
!!%          ENDDO
!!%       ENDDO
!!%
!!%          ENERGY  = E_HYD(I)   
!!%
!!%          WRITE(nout,*)  '#   INPUT ENERGY  E = ', ENERGY
!!%       
!!%
!!%!          C(I,:)  = 0.0D+00
!!%          INFO      = 0
!!%          REL_ERROR = 0.0D+00
!!%          LW        = NN * ( KB + 1 )
!!%          IW        = 0.0D+00
!!%          D(1)      = 1.0D+00 
!!%     
!!%          CALL F02SDF(NN, KB, KB, AA, 2*K+1, BB, 2*K+1, .FALSE., REL_ERROR, ENERGY, C(I,:),&
!!%               &                                   D, IW, W, LW, INFO)
!!%
!!%
!!%          IF(INFO.NE.0) THEN 
!!%
!!%             WRITE(nout,*) '# ERROR IN  F02SDF/NAG_77/M20 OR LINEAR SUBROUTINE '
!!%             WRITE(nout,*) '#'
!!%             WRITE(nout,*) '#                     INFO = ', INFO
!!%             WRITE(nout,*) '#                     D(1) = ', D(1)
!!%             WRITE(nout,*) '#'
!!%        
!!%             STOP
!!%             
!!%     ELSE
!!%        
!!%        WRITE(nout, *) '#    INVERSE ITERATION : '
!!%        WRITE(nout, *) '#'
!!%        WRITE(nout, *) '#          NBS  = ', NBS
!!%        WRITE(nout, *) '#      NCS  = ', NCS
!!%        WRITE(nout, *) '#   INPUT     ENERGY     E = ', ENERGY
!!%        WRITE(nout, *) '#   CORRECTED ENERGY E_INV = ', ENERGY + D(30)
!!%        WRITE(nout, *) '#                     DE/E = ', ABS(D(30) / ENERGY) 
!!%            
!!%
!!%        NORM = 0.0D+00
!!%
!!%        NORM = DOT_PRODUCT(C(I,:), MATMUL(B, C(I,:) ) )
!!%          
!!%
!!%        C(I,:) = C(I,:) / SQRT(NORM)
!!%
!!%
!!%        WRITE(nout, *) '#           NORM     = ', NORM
!!%        WRITE(nout, *) '#           COEFF_C(I) = ', I
!!%
!!%        DO  J = 1, NB/2
!!%           WRITE(NOUT, "(I4,1X,E25.14,2X,I4,1X,E25.14)") J, C(I,J), J+NB/2, C(I,J+NB/2)
!!%        ENDDO
!!%
!!%
!!%     ENDIF
!!%
!!%     ENDDO
!!%
!!%     CLOSE(NOUT)
!!%
!!%     RETURN
!!%
!!%   END SUBROUTINE SOLVE_BBANDED_INVERSE_NAG77
!!%
!!%!###################################################### 
!!%SUBROUTINE SOLVE_BBANDED_LIN_NAG77(AB, BB, E_C, C)
!!%  
!!%
!!%!  USE nag_f77_f_chapter
!!%  USE param, ONLY : dpk
!!%  USE param, ONLY : nbs, ncs
!!%  USE param, ONLY : nb, kb, rmax, znuc
!!%  IMPLICIT NONE
!!%
!!%
!!%  INTEGER I, J, INFO 
!!%  INTEGER K
!!%
!!%
!!%!...........
!!%
!!%  INTEGER  NN
!!%  
!!%  REAL(DPK), DIMENSION(NB-1,NB-1), INTENT(IN) :: aB,bB
!!%  REAL(DPK), DIMENSION(NB-1,1):: c
!!%  REAL(DPK), DIMENSION(NCS), INTENT(OUT) :: e_c
!!%
!!%
!!%! linear equation solver 
!!%
!!%  REAL(DPK) ENERGY
!!%  REAL(DPK), DIMENSION(3*KB-2, NB-1) :: a, b, h  ! linear
!!%  REAL(DPK)  NORM
!!%  INTEGER, DIMENSION(NB-1)   :: iw
!!%  INTEGER NOUT
!!%
!!%  REAL(DPK), ALLOCATABLE,DIMENSION(:) :: E_FREE, E_HYD
!!%!..............................
!!%
!!%  WRITE(*,*) '# h1e::solve_bbanded_lin_nag77:       call.'
!!%
!!%  nout = 12
!!%
!!%
!!%!..............................
!!%
!!%
!!%      K   = KB - 1    ! number of sub-diagonials (upper diagonals)
!!%      NN  = NB - 1    ! dimension of coefficient vector
!!%
!!%
!!%
!!%      ALLOCATE( E_FREE(NB-2) ) 
!!%      ALLOCATE( E_HYD(NB-2)  ) 
!!%
!!%!............................
!!%
!!%!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES
!!%
!!%
!!%      OPEN(NOUT, FILE='out/subsolve_lin.out')
!!%
!!%      WRITE(nout, *) '#        DIMENSION OF BANDED MATRICES  (K+1) X N '
!!%      WRITE(nout, *) '#'
!!%      WRITE(nout, *) '#                                K + 1 = KB     = ', K + 1
!!%      WRITE(nout, *) '#                                  NN  = NB - 2 = ', NN
!!%      WRITE(nout, *) '#'
!!%      WRITE(nout, *) "# STATE        E_FREE"
!!% 
!!%
!!%       DO  I = 1, NB-2
!!% 
!!%          E_FREE(I) = 0.5D+00* ( I * 3.141592654D+00/ RMAX )**2
!!%          E_HYD(I) = - 0.5D+00 * (ZNUC/I)**2
!!%
!!%
!!%          WRITE(nout,"(I4,3X,2E20.10)") I, E_HYD(I), E_HYD(I)
!!%
!!%
!!%       ENDDO
!!%
!!%!! LINEAR MATRICES 
!!%
!!%!.........................................
!!%
!!%
!!%
!!%          DO J = 1, NN
!!%             DO I = MAX(1, J - K), MIN(NN, J + K)
!!%
!!%                A( K + K + I - J + 1, J ) =   AB(I, J) 
!!%                B( K + K + I - J + 1, J ) =   BB(I, J) 
!!%
!!%             ENDDO
!!%          ENDDO
!!%
!!%!..........................................
!!%!     NSTATES = 1
!!%
!!%
!!%          E_C(1)  = 0.001 !E_HYD(NBS) 
!!%
!!%
!!%          WRITE(nout, *)  '#   INPUT ENERGY  E = ', E_C(1)
!!%     
!!%          
!!%!          DO I = 1, 3*KB - 2
!!%!             DO J = 1, NB - 1
!!%
!!%                A = A - E_C(1) * B
!!%
!!%!             ENDDO
!!%!       ENDDO
!!%
!!%          IW = 0.0D+00
!!%
!!%
!!%        ! FACTORIZE FIRST
!!%
!!%          CALL F07BDF(NN, NN, K, K, A, 3*KB-2, IW, INFO) 
!!%
!!%
!!%
!!%          IF(INFO.NE.0) THEN 
!!%
!!%             WRITE(nout,*) '# ERROR IN FACTORIZATION SUBROUTINE '
!!%             WRITE(nout,*) '#'
!!%             WRITE(nout,*) '#                     INFO = ', INFO
!!%             WRITE(nout,*) '#'
!!%
!!%
!!%             CALL X04CEF(NN, NN, K, 2*K, A, NN,'Details of factorization', INFO)
!!%
!!%
!!%          ELSE
!!%
!!%             
!!%             C        = 0.0D+00
!!%             C(NN, 1) = 1.0D+00
!!%
!!%             CALL F07BEF('N', NN, K, K, 1, A, 3*KB-2, IW, C, NN,INFO)
!!%
!!%          ENDIF
!!%
!!%
!!%
!!%      IF(INFO.NE.0) THEN 
!!%
!!%         WRITE(nout,*) '#         ERROR IN   F07BEF          '
!!%         WRITE(nout,*) '#'
!!%         WRITE(nout,*) '#                     INFO = ', INFO
!!%         WRITE(nout,*) '#'
!!%
!!%         STOP
!!%
!!%      ELSE
!!%
!!%
!!%
!!%         WRITE(nout, *) '#    LINEAR SOLUTION : '
!!%         WRITE(nout, *) '#'
!!%         WRITE(nout, *) '#   INPUT     ENERGY      E = ', E_C(1)
!!%         WRITE(nout, *) '#'
!!%         WRITE(nout, *) '#        NCS  = ', NCS
!!%
!!%
!!%
!!%        
!!%         NORM = DOT_PRODUCT( C(:,1), MATMUL( BB, C(:, 1) ))
!!%          
!!%         C(:,1) = C(:,1) / SQRT(NORM)
!!%
!!%
!!%         WRITE(nout, *) '#           NORM     = ', NORM
!!%         WRITE(nout, *) '#           COEFF(I) = '
!!%         
!!%         DO  I = 1, NN/2
!!%
!!%            WRITE(nout, "(I4,1X,E25.14,2X,I4,1X,E25.14)") I, C(I,1), I+NB/2, C(I+NB/2,1)
!!%         ENDDO
!!%
!!%
!!%      ENDIF
!!%
!!%
!!%      CLOSE(NOUT)
!!%
!!%      RETURN
!!%
!!%       
!!%    END SUBROUTINE SOLVE_BBANDED_LIN_NAG77
!!%
!!%!###################################################### 
!!%SUBROUTINE SOLVE_BANDED_INVERSE_NAG77(L, H, B, ENERGY_BOUND, C, CHOICE)
!!%  
!!%!  USE nag_f77_f_chapter
!!%  USE PRECISION, ONLY : DPK
!!%  USE PARAM,     ONLY : nbs
!!%  USE PARAM,     ONLY : nb, kb, rmax, znuc
!!%  USE IOROUTINES, ONLY: d1efile
!!%!......................
!!%
!!%  IMPLICIT NONE
!!%
!!%
!!%  INTEGER I,II, J, INFO 
!!%  INTEGER N, K 
!!%
!!%  INTEGER  NN, L
!!%
!!%  REAL(DPK), DIMENSION(NB-1,NB-1), INTENT(IN) :: H, B
!!%  REAL(DPK), DIMENSION(NBS, NB-1),INTENT(inout) :: c
!!%
!!%  REAL(DPK), DIMENSION(KB, NB-2) :: hb, bb         !  fxd diagonalization 
!!%  REAL(DPK), DIMENSION(NB-2) :: er
!!%
!!%  REAL(DPK), DIMENSION(2*KB-1, NB-1) :: ainv, binv    ! inverse  iteration
!!%  REAL(DPK), DIMENSION(NBS + L), INTENT(OUT):: energy_bound
!!%
!!%  REAL(DPK), DIMENSION(30)   :: D
!!%  REAL(DPK), DIMENSION( (NB-1)*(KB+1) )   :: w
!!%
!!%  REAL(DPK) REL_ERROR, NORM, ENERGY
!!%  INTEGER, DIMENSION(NB-1)   :: iw
!!%  INTEGER LW, INDEX_B
!!%  INTEGER NOUT, ND1EFILE
!!%  INTEGER CHOICE
!!%  REAL(DPK), ALLOCATABLE,DIMENSION(:) :: E_FREE, E_HYD
!!%
!!%!..............................
!!%  WRITE(*, *) '#'
!!%  WRITE(*,*) '# h1e::solve_banded_inverse_nag77:    call.'
!!%
!!%  ND1EFILE = 11
!!%  NOUT     = 12
!!%
!!%  
!!%!..............................
!!%
!!%      K   = KB - 1
!!%      N   = NB - 2
!!%      NN  = NB - 1
!!%
!!%      ALLOCATE(E_FREE( N )) 
!!%      ALLOCATE(E_HYD( NBS + L ) ) 
!!%
!!%!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES
!!%
!!%      OPEN(NOUT, FILE='out/subsolve_inv.out')
!!%
!!%      WRITE(*, *) '# h1e::solve_banded_inverse_nag77: banded matrices dimensions: (kb+1) x (nb-2) '
!!%      WRITE(*, *) '# h1e::solve_banded_inverse_nag77:     bandwidth k = ', k + 1
!!%      WRITE(*, *) '# h1e::solve_banded_inverse_nag77:     n  = nb - 2 = ', n
!!%      WRITE(*, *) '# h1e::solve_banded_inverse_nag77:    nn  = nb - 1 = ', nn
!!%      WRITE(*, *) '# h1e::solve_banded_inverse_nag77:    fxd boundary diagonalization : '
!!%
!!%
!!%      HB = 0.0D+00
!!%      BB = 0.0D+00
!!%
!!%      DO J = 1, N
!!%         DO I = MAX(1, J - K), J
!!%
!!%            HB( K + 1 + I - J, J ) =   H(I, J) 
!!%            BB( K + 1 + I - J, J ) =   B(I, J) 
!!%
!!%         ENDDO
!!%      ENDDO
!!%
!!%
!!%!     
!!%      !!   SOLVES THE BANDED FORM OF THE SCHRONDIGER HAMILTONIAN EQNS  
!!%      !!   DIAGONALIZATION IS PERFORMED
!!%      !!   
!!%      !!                  ONLY EIGENVALUES
!!%      !!
!!%
!!%
!!%      INFO   = 0
!!%      LW     = MAX( N, ( 3*KB + KB ) * ( KB + KB + 1) )
!!%
!!%      CALL F02FHF(N, K, HB, KB, K, BB, KB, ER, W, LW, INFO)
!!%
!!%      WRITE(*, *) '# h1e::solve_banded_inverse_nag77:    fxd boundary diagonalization done.'
!!%!................................
!!%
!!%      CALL D1EFILE(ND1EFILE, L) 
!!%      
!!%      INDEX_B = 0
!!%      DO I = 1, NB - 2
!!%
!!%         IF(ER(I).LT.0.0D+00) INDEX_B = INDEX_B + 1
!!%
!!%         WRITE(ND1EFILE) ER(I)
!!%      ENDDO
!!%
!!%      WRITE(ND1EFILE) INDEX_B
!!%
!!%      CLOSE(ND1EFILE)
!!%
!!%!................................
!!%      
!!%      IF(INFO.NE.0) THEN 
!!%
!!%         WRITE(*,*) '# h1e::solve_banded_inverse_nag77: error in F02FHF/NAG_77/M20 .stop '
!!%         WRITE(*,*) '# h1e::solve_banded_inverse_nag77:                           info = ', info
!!%         WRITE(*,*) '#'
!!%         STOP
!!%
!!%      ENDIF
!!%
!!%!............................
!!%      
!!%
!!%      WRITE(*,*) " STATE        E_DISCRETE          E_FREE             E_HYDROGEN"
!!% 
!!%       DO  I = L + 1 , NBS + L
!!% 
!!%          E_FREE(I) = 0.5D+00* ( I * 3.141592654D+00/ RMAX )**2
!!%          E_HYD(I) = - 0.5D+00 * (ZNUC/I)**2
!!%
!!%          IF(ER(I)<0.0D+00) THEN
!!%
!!%             WRITE(*,"(I4,3X,3E20.10)") I, ER(I), E_HYD(I), E_HYD(I)
!!%          ELSE
!!%
!!%             WRITE(*,"(I4,3X,3E20.10)") I, ER(I), E_FREE(I), E_HYD(I)
!!%          ENDIF
!!%
!!%       ENDDO
!!%
!!%
!!%!       WHERE(ER < 0.0D+00)  ENERGY_BOUND = ER
!!%!       WHERE(ER < 0.0D+00)  ENERGY_BOUND = E_HYD
!!%
!!%
!!%!.......................    INVERSE ITERATION
!!%
!!%
!!%       WRITE(*, *) '# h1e::solve_banded_inverse_nag77:    inverse iteration call.'
!!%
!!%       IF( CHOICE == 1) THEN 
!!%
!!%          WRITE(*, *) '# h1e::solve_banded_inverse_nag77:                    l = ',l
!!%
!!%          DO I = 1, NBS 
!!%
!!%             
!!%             DO J = 1, NN
!!%                DO II = MAX(1, J - K), MIN(NN, J + K)
!!%                   
!!%                   AINV( K + 1 + II - J, J ) =   H(J, II) 
!!%                   BINV( K + 1 + II - J, J ) =   B(J, II) 
!!%
!!%                ENDDO
!!%             ENDDO
!!%
!!%!             ENERGY  = E_HYD( I + L )
!!%
!!%!! 28.10.2005LAAN
!!%!! note ( inverse iteration method extremely unstable!
!!%!! for instance: for energy = -0.5              wf comes out o.k.
!!%!!                   energy = -0.5000000000013  wf comes out bad  
!!%! 
!!%             ENERGY  = ER( I + L )
!!%!              ENERGY = -0.501D+00
!!%
!!%             ENERGY_BOUND(I)  = ENERGY
!!%
!!%!             WRITE(*, *) '# h1e::solve_banded_inverse_nag77:      nbs = ', i
!!%!             WRITE(*, *) '# h1e::solve_banded_inverse_nag77:      hen = ',e_hyd(i+l)
!!%!             WRITE(*, *) '# h1e::solve_banded_inverse_nag77:      den = ',er(i+l)
!!%!             WRITE(*, *) '# h1e::solve_banded_inverse_nag77:       en = ',energy
!!%
!!%          
!!%             INFO      = 0
!!%             REL_ERROR = 0.0D+00
!!%             LW        = NN * ( KB + 1 )
!!%             IW        = 0.0D+00
!!%             D(1)      = 1.0D+00 
!!%             
!!%             CALL F02SDF(NN, KB, KB, AINV, 2*K+1, BINV, 2*K+1    &
!!%                  &          , .FALSE., REL_ERROR, ENERGY, C(I,:)    &
!!%                  &          , D, IW, W, LW, INFO)
!!%          
!!%          
!!%             IF(INFO.NE.0) THEN 
!!%                WRITE(*,*) '# h1e::solve_banded_inverse_nag77: error in F02SDF/NAG_77/M20, stop '
!!%                WRITE(*,*) '# h1e::solve_banded_inverse_nag77:                           info = ', info
!!%                WRITE(*,*) '# h1e::solve_banded_inverse_nag77:                           D(1) = ', d(1)
!!%                WRITE(*,*) '#'
!!%                STOP
!!%             ELSE
!!%!                WRITE(*, *) '# h1e::solve_banded_inverse_nag77: nof bound_states   nbs = ', nbs
!!%!                WRITE(*, *) '# h1e::solve_banded_inverse_nag77: input energy        en = ', energy
!!%!                WRITE(*, *) '# h1e::solve_banded_inverse_nag77: correct. energy  e_inv = ', ENERGY + D(30)
!!%!                WRITE(*, *) '# h1e::solve_banded_inverse_nag77:                   de/e = ',ABS(D(30) / ENERGY) 
!!%             
!!%                NORM = 0.0D+00
!!%             
!!%                NORM = DOT_PRODUCT(C(I,:), MATMUL(B, C(I,:) ) )
!!%                          
!!%                C(I,:) = C(I,:) / SQRT(NORM)
!!%             
!!%!                WRITE(*, *) '# h1e::solve_banded_inverse_nag77:   norm   = ', norm
!!%!                WRITE(*, *) '# h1e::solve_banded_inverse_nag77: coeff(i) = ', i
!!%             
!!%!             DO  J = 1, NB/2
!!%!                WRITE(*, "(I4,1X,E25.14,2X,I4,1X,E25.14)") J, C(I,J), J+NB/2, C(I,J+NB/2)
!!%!             ENDDO
!!%             
!!%             
!!%             ENDIF
!!%          
!!%          ENDDO
!!%       
!!%       ENDIF
!!%    
!!%       CLOSE(NOUT)
!!%    
!!%       RETURN
!!%    
!!%     END SUBROUTINE SOLVE_BANDED_INVERSE_NAG77
!!%      FUNCTION pscoul(l, qz, qk)
!!%
!!%      USE units, ONLY: M_PI
!!%      IMPLICIT COMPLEX*16(a-h,r-z), REAL*8(o-q)
!!%           
!!%      q1 = l + 1
!!%
!!%      ak = CMPLX(0.0D+00, qz / qk )
!!%      
!!%      ga = q1 - ak
!!%
!!%      cc = log(ga) + log(1.d0 + ga) + log(2.d0 + ga) + log(3.d0 + ga)
!!%      cc = cc + log(4.d0+ga)+ log(5.d0+ga)+log(6.d0+ga)+log(7.d0+ga)
!!%
!!%      z = ga + 8.0D+00
!!%      
!!%      gam=-cc+(z-0.5d0)* log(z)-z+0.5d0* log(2.d0*M_PI)+1.d0/(12.d0*z)
!!%
!!%      gam=gam-1.d0/(360.d0*z**3)+1.d0/(1260.d0*z**5)-1.d0/(1680.d0*z**7)
!!%
!!%      pscoul= aimag(gam)
!!%
!!%      RETURN
!!%    END FUNCTION pscoul
!###################################################### 
