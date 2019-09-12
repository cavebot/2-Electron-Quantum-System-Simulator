!qub/24102008/laan/:
!  --prints selected eigenvalues/eigenvectors from the h2eb.f90 output data
!  -- transformed to f90 from old bspci2e/cpc/wf2e.f progran 
!

PROGRAM w2eb
  !
  USE precision, ONLY:dpk
  USE units,     only:enau
  !
  IMPLICIT NONE
  !
  INTEGER                                :: ncs, nhx
  REAL(dpk)                              :: en_2e               !
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: c_w2e               !nhx
  INTEGER,   ALLOCATABLE, DIMENSION(:)   :: nhf, lhf, ll        !ncs
  INTEGER,   ALLOCATABLE, DIMENSION(:)   :: nllmin, nllmax, nd  !ncs
  !
  INTEGER                                :: l2e, l, ls
  INTEGER                                :: n2e
  !
  INTEGER                                :: ie, k 
  INTEGER                                :: i_e2e_threshold  = 0
  REAL(dpk)                              :: e2e_threshold    = 0.0_dpk
  LOGICAL                                :: print_out        = .false.
  !
  !i/o!
  CHARACTER(len=100)                     :: argv 
  INTEGER                                :: nh2eb = 1
  INTEGER                                :: nw2eb = 2
  
  !


  CALL GETARG(1, ARGV)
  READ(ARGV,*) l2e
  CALL GETARG(2, ARGV)
  READ(ARGV,*) n2e

  !
  CALL hfile(nh2eb,"dat","h2eb","bin",l2e)       !  CALL d2efile(9, lo)

  READ(nh2eb)   l, ls

  IF(l.NE.l2e) THEN
     WRITE(*,*) '# w2eb:: inconsistent input and data file'
     WRITE(*,*) '# w2eb::                input       l2e = ', l2e
     WRITE(*,*) '# w2eb::                datafile      l = ', l
     STOP
  ENDIF
  
  READ(nh2eb)   ncs
  
  ALLOCATE(    nhf(ncs) )
  ALLOCATE(    lhf(ncs) )
  ALLOCATE(     ll(ncs) )
  ALLOCATE( nllmin(ncs) ) 
  ALLOCATE( nllmax(ncs) )
  ALLOCATE(     nd(ncs) )

  READ(nh2eb) ( nhf(k),    k = 1, ncs)
  READ(nh2eb) ( lhf(k),    k = 1, ncs)
  READ(nh2eb) ( ll(k),     k = 1, ncs)
  READ(nh2eb) ( nllmin(k), k = 1, ncs)
  READ(nh2eb) ( nllmax(k), k = 1, ncs)
  READ(nh2eb) ( nd(k),     k = 1, ncs)
  READ(nh2eb)   nhx

  ALLOCATE(c_w2e(nhx))      

  WRITE(*,'(a10,30I3)') '# lhf = ', lhf
  WRITE(*,'(a10,30I3)') '#  ll = ', ll


  CALL hfile(nw2eb,"dat","w2eb","bin",l)       !  CALL d2efile(9, lo)

  WRITE(nw2eb) l,ls
  WRITE(nw2eb) ncs
  WRITE(nw2eb) ( nhf(k),    k = 1, ncs)
  WRITE(nw2eb) ( lhf(k),    k = 1, ncs)
  WRITE(nw2eb) ( ll(k),     k = 1, ncs)
  WRITE(nw2eb) ( nllmin(k), k = 1, ncs)
  WRITE(nw2eb) ( nllmax(k), k = 1, ncs)
  WRITE(nw2eb) ( nd(k),     k = 1, ncs)
  WRITE(nw2eb)   nhx
  WRITE(nw2eb)   n2e

  !
  WRITE(*,*) '# w2e::    2-e partial wave l2e = ', l2e
  WRITE(*,*) '# w2e::    nof 2-e states   n2e = ', n2e
  WRITE(*,*) '# w2e::    2e matrix dim    nhx = ', nhx
  WRITE(*,*) '#'
  WRITE(*,'(a5,a25,a25)') '# n','E(n)/a.u.', 'E(n)/eV'

  IF(n2e < 0) THEN
     n2e = - n2e
     print_out = .true.
  ENDIF

  read_and_write:DO  ie = 1, n2e

     READ(nh2eb) en_2e
     READ(nh2eb) (c_w2e(k), k = 1, SIZE(c_w2e))


     IF(en_2e > e2e_threshold.and.i_e2e_threshold.eq.0) i_e2e_threshold = ie

     print_max_state:IF(ie==n2e.OR.ie==1.or.ie==i_e2e_threshold) THEN 
        WRITE(*,'(i5,3E25.12,a10)') ie, en_2e/2.0_dpk,enau*en_2e/2.0_dpk
     ENDIF print_max_state
     IF(print_out)         CYCLE
        
     WRITE(nw2eb) en_2e
     WRITE(nw2eb) (c_w2e(k), k = 1, SIZE(c_w2e))
      
  ENDDO read_and_write
  
  CLOSE(nh2eb)
  CLOSE(nw2eb)
  
END PROGRAM w2eb
!EOF
