MODULE w2e

PUBLIC diagonalize, read_cfg_file

CONTAINS

SUBROUTINE diagonalize(nh2eb, nhx, h, u)
  !
  USE PRECISION, only:dpk
  !
  IMPLICIT NONE 
  !     
  INTEGER,             INTENT(in) :: nh2eb
  INTEGER,             INTENT(in) :: nhx
  REAL(dpk),  DIMENSION(nhx,nhx)  :: h
  REAL(dpk),  DIMENSION(nhx,nhx)  :: u
  !
  REAL(dpk),  DIMENSION(nhx)      :: c
  !
  !external subroutine entries
  REAL(dpk), DIMENSION(5*nhx)     :: iwork
  REAL(dpk), DIMENSION(10*nhx)    :: work
  INTEGER,   DIMENSION(nhx)       :: ifail
  INTEGER                         :: lwork, il, iu, info, m
  REAL(dpk)                       :: abstol, vl,vu
  !
  INTEGER                         :: j,k

  !exe!


  !**  lapack routine
  
  lwork  = SIZE(work)
  vl     = 0.0_dpk
  vu     = 0.0_dpk
  il     = 1
  iu     = nhx
  abstol = 0.0_dpk


  PRINT*, "this is a test", nhx
  
!diagonalize

  CALL dsyevx('V','I','U', nhx, h, nhx, vl, vu, il, iu, abstol, m, &
       & c, u, nhx, work, lwork, iwork, ifail, info)

! and save

!nhx = size(c)

  WRITE(nh2eb) SIZE(c)
  save_eigen_values_coe:DO  k = 1, SIZE(c)
     WRITE(nh2eb)  c(k) 
     WRITE(nh2eb) (u(j,k), j = 1, SIZE(c) )
  ENDDO save_eigen_values_coe

  DO k = 1, nhx 
     WRITE(*,'(i6,2x,2E18.8)') k, c(k)/2.,c(k)*27.211396_dpk/2.0_dpk
  ENDDO
!  CLOSE(nh2eb)
  
!!%  DO k = 1, nhx
!!%
!!%     k1 = 0
!!%     sumProb10 = 0.0d+00
!!%
!!%     DO  i = 1, 10       
!!%        p(i) = 0.0D+00         
!!%        DO  j = 1, nd(i)
!!%            
!!%           k1 = k1 + 1
!!%
!!%           IF(k1.GT.nhmx)  CYCLE
!!%
!!%           p(i) = p(i) +  u(k1, k) * u(k1, k)
!!%
!!%        ENDDO
!!%
!!%        !         write(*,*) i, p(i)
!!%
!!%         sumProb10  = sumProb10  + p(i)
!!%        
!!%      ENDDO
!!%
!!%      WRITE(*,8) k, en(k)*0.5D+00,( p(i), i = 1, 6), sumProb10
!!%   ENDDO

  RETURN
END SUBROUTINE diagonalize

!#######################################################################
      
SUBROUTINE read_cfg_file(l, s, n1, l1, l2, n2_min, n2_max, idcs, ncs, ndim)
        !
        IMPLICIT NONE
        !
        INTEGER,                INTENT(in)  :: l        ! total angular momentum
        INTEGER,               INTENT(out)  :: s        ! total spin
        INTEGER, DIMENSION(:),     POINTER  :: n1, l1   !
        INTEGER, DIMENSION(:),     POINTER  :: l2, n2_min, n2_max
        INTEGER, DIMENSION(:),     POINTER  :: idcs     ! 0 or 1
        INTEGER,               INTENT(out)  :: ncs      ! nof channels included
        INTEGER,               INTENT(out)  :: ndim
        !
        INTEGER                             :: n_total_l, inpl
        INTEGER                             :: k
        !

        integer                             :: ncfg
        !
        !EXE!
        

        ncfg = 1

        !


        CALL cfile(ncfg,"inp","cfg",l)

        READ(ncfg, *) n_total_l
        READ(ncfg, *) inpl, s
        READ(ncfg, *) ncs

        WRITE(*,*) n_total_l, inpl, s, ncs

        check_partial_wave:IF( L.NE.inpl ) THEN
           WRITE(*,*) '# read_cfg: inconsistent cfg file'
           WRITE(*,*) '#                            L = ', l
           WRITE(*,*) '#                         inpl = ', inpl
           STOP
        END IF check_partial_wave
        

        
        ALLOCATE(    n1(ncs) )
        ALLOCATE(    l1(ncs) )
        ALLOCATE(    l2(ncs) )
        ALLOCATE(n2_min(ncs) )
        ALLOCATE(n2_max(ncs) )
        ALLOCATE(  idcs(ncs) )


        idcs = 0
        ndim = 0
        read_cfg:DO k = 1, ncs
           READ(ncfg,*) n1(k), l1(k), l2(k), n2_min(k), n2_max(k)!, idcs(k)
           
           ndim = ndim + n2_max(k) - n2_min(k) + 1
        ENDDO read_cfg
        

        CLOSE(ncfg)
        

        WRITE(*,*) '# read_cfg::       hamiltonian matrix dimension   ndim = ', ndim
        WRITE(*,*) '# read_cfg::        nof channels included          ncs = ', ncs        
        

        RETURN

      END SUBROUTINE read_cfg_file

    END MODULE w2e
!eof
!
