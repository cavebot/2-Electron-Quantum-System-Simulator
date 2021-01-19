!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
MODULE param
  !
  USE PRECISION, ONLY : DPK
  USE io, only:ninp
  !
  IMPLICIT NONE
  !
  PUBLIC
  !
  INTEGER, PARAMETER                  :: kx = 15     ! should conform with fxd code
  INTEGER, PARAMETER                  :: np = 2000   ! should conform with fxd code


!..........  h1e.inp file   ...........

  REAL(DPK)                           :: znuc, mass
  CHARACTER(len=20)                   :: potential
  REAL(DPK)                           :: r_m              ! internuclear distance
  INTEGER                             :: nl
  !

  !
  REAL(dpk)                           :: rmin, rmax, xi, rs
  INTEGER                             :: no, idbsp, nrp
  INTEGER                             :: nb, kb, ib
  INTEGER                             :: first_spline, last_spline
  INTEGER                             :: ndim, nldim
  !
  INTEGER                             :: l1,l2
  !
  CHARACTER(LEN=3)                    :: method, spectrum
  REAL(dpk)                           :: en_1e, de_1e
  INTEGER                             :: nbs,  ncs
  !hartree-fock
  INTEGER                             :: ihf, imodel
  INTEGER                             :: ncore(3), no_s,  no_p, no_d
  INTEGER                             :: ntc, nsp
  INTEGER                             :: lcore
  REAL(dpk)                           :: ap, rc
  REAL(dpk), DIMENSION(:),ALLOCATABLE :: zk, ek
  INTEGER,   DIMENSION(:),ALLOCATABLE :: lk
!  REAL(dpk), DIMENSION(:),ALLOCATABLE :: ri, dri
  REAL(dpk)                           :: crt, di, df
  INTEGER                             :: itrmx, id

CONTAINS
!!!xxxxxxxxxxxxxxxxxxxxx

  SUBROUTINE input_read_write(v_0, z_n, r_b, n_b)
    !
    IMPLICIT NONE
    !
    REAL(dpk)              :: v_0
    REAL(dpk)              :: z_n
    REAL(dpk)              :: r_b
    INTEGER                :: n_b
    !
    CHARACTER(len=20)      :: h1efile
    !

    ! read general file

!!!...............
    OPEN(ninp,file ='inp/bsci2e.inp',status='old')
! 'hydrogenic' atomic parameters
    READ(ninp, *) h1efile               ! input file for h1e/w1e/d1e/r12/h2e/d2e
    READ(ninp, *) mass                  ! Atomic number
    READ(ninp, *) potential             ! system (h, h2p, model)
! basis parameters
    READ(ninp, *) kb                    ! order of B-splines
    READ(ninp, *) rs                    ! rs = first non-zero point)(sine)
    READ(ninp, *) idbsp                 ! grid (0=sin,1=lin)
    READ(ninp, *) no                    ! interval (rmin, rmax), no: nof grid points
    READ(ninp, *) rmin                  ! rmin = first non-zero point for sine-grid (1,..,no), (idbsp=0)
    READ(ninp, *) nrp                   ! exponent for exp-like knots (idbsp = 1)
!spectrum information
    READ(ninp, *) first_spline
    READ(ninp, *) last_spline           ! 1st(last)_spline(0/1) is (excluded/included),
    READ(ninp, *) method
    READ(ninp, *) spectrum              ! d(diag), dl(diag-linear) l(linear/inverse)  (fxd,free,mxd)
    READ(ninp, *) en_1e
    READ(ninp, *) de_1e
    READ(ninp, *) ncs
    READ(ninp, *) nbs
    !
    READ(ninp, *) ap                    ! ap potential = 'hn','model'
    READ(ninp, *) rc                    ! ap potential = 'hn','model'
    !hartree-fock
    READ(ninp,*) ihf

    hf_parameters:IF(ihf.EQ.1) THEN

       READ(ninp, *) ncore           ! ns, nop, nod
       READ(ninp, *) lcore           ! outer-shell core angular momentum

       ALLOCATE(zk(1:SUM(ncore)) )

       READ(NINP, *) zk              ! effective charges
       READ(NINP, *) crt             !
       READ(NINP, *) di, df, id      !
       READ(ninp, *) itrmx           ! HF iteration for core orbitals
    ENDIF hf_parameters

    !       READ(NINP, *) ( alp1(i + lmin - 1), i = 1, nw )
    !       READ(NINP, *) (  r01(i + lmin - 1), i = 1, nw )

    CLOSE(ninp)


    ap = v_0

    !    ALLOCATE(  ri(no) )
    !    ALLOCATE( dri(no) )


    !
    !  write specialized input file for system, box radius and basis size
    !

    IF(kb.GT.kx) THEN
       PRINT*, "# hard-coded parameter kx in bs1e_parameter.f90  kx = ", kx
       PRINT*, "#                      kb in bspci2e.inp         kb = ", kb
       PRINT*, "#                                    must be kb <= kx "
    ENDIF

    IF(np.NE.no) THEN
       PRINT*, "# hard-coded parameter 'np' in bs1e_parameter.f90  np = ", np
       PRINT*, "#                      'no' in bspci2e.inp         no = ", no
       PRINT*, "#                                         must be equal "
       STOP
    ENDIF


    h1efile = "inp/"//trim(adjustl(h1efile))
    OPEN(ninp,file =h1efile)
! 'hydrogenic' atomic parameters
    WRITE(ninp, '(2f15.5,10x,a10)')           z_n, mass,                 '# zn, mass'
    WRITE(ninp, '(a10,2f10.5,10x,a25)')       potential, ap, rc,         '# model potential, ap, rc'
! basis parameters
    WRITE(ninp, '(2I15,10x,a8)')              n_b, kb,                   '# nb, kb'
    WRITE(ninp, '(1E10.3,F10.5,I10,10x,a14)') rmin, r_b, no,             '# rmin, rb, no'
    WRITE(ninp, '(i10,1E10.3,i10,10x,a54)')   idbsp, rs, nrp,            '# grid (0=sin,1=lin), rs = first non-zero point)(sine)'
!spectrum information
    WRITE(ninp, '(2I15,10x,a45)')             first_spline, last_spline, '# 1st(last)_spline(0/1) is (excluded/included)'
    WRITE(ninp, '(2a15,10x,a60)')             method, spectrum,  '# d(diag), dl(diag-linear) l(linear/inverse)  (fxd,free,mxd)'
    WRITE(ninp, '(2f10.5,2I5,10x,a18)')       en_1e, de_1e, ncs, nbs,    '# e0, de, ncs, nbs'
    WRITE(ninp, '(1I15,10x,a50)')             ihf,                       '# ihf = 1, perform a HF calculation'

    hf_write_parameters:IF(ihf.EQ.1) THEN
       WRITE(ninp, '(1I15,10x,a45)')          lcore,                     '# lcore (outer-shell ang. mom)'
       WRITE(ninp, '(3I10,10x,a45)')          ncore,                     '# ncore = nos,nop,nod, lcore (outer-shell ang. mom)'
       WRITE(ninp, '(3f10.5,7x,a45)')         zk,                        '# zk: effective charges for shell orbitals'
       WRITE(ninp, '(1E10.3,2f10.5,1x,a33)')  crt, di, df,               '# HF control parameters '
       WRITE(ninp, '(2I15,7x,a14)')           id, itrmx,                 '# id, itrmx'
    ENDIF hf_write_parameters

    WRITE(ninp, '(a40)') '###############################################'
    WRITE(ninp, '(a1)') '#'
    WRITE(ninp, '(a40)') '#  input file for running bsci2e programs created'
    WRITE(ninp, '(a1)') '#'
    CLOSE(ninp)
    PRINT*, '#  output file = ', h1efile




  END SUBROUTINE input_read_write


!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  SUBROUTINE input

    IMPLICIT NONE

!    INTEGER  :: i

!!!...............
    OPEN(ninp,file ='inp/h1e.inp',status='old')
! 'hydrogenic' atomic parameters
    READ(ninp, *) znuc, mass                ! Atomic number
    READ(ninp, *) potential, ap, rc         ! model potential,
! basis parameters
    READ(ninp, *) nb, kb                    ! nof B-splines, order of B-splines
    READ(ninp, *) rmin, rmax, no            ! interval (rmin, rmax), no: nof knot points
    READ(ninp, *) idbsp, rs, nrp            ! grid (0=sin,1=lin), rs = first non-zero point)(sine)
    !    READ(ninp, *) rs                   !                     rs = 0, nrp (only used in exp knots)
    !spectrum information
    READ(ninp, *) first_spline, last_spline ! 1st(last)_spline(0/1) is (excluded/included),
    READ(ninp, *) method, spectrum          ! d(diag), dl(diag-linear) l(linear/inverse)  (fxd,free,mxd)
    READ(ninp, *) en_1e, de_1e, ncs, nbs    !
    READ(ninp, *) ihf

    hf_write_parameters:IF(ihf.EQ.1) THEN

       READ(ninp, *)  lcore
       READ(ninp, *)  ncore
       ALLOCATE(zk(1:SUM(ncore)))
       READ(ninp, '(3f10.5)')  zk
       READ(ninp, *)  crt, di, df
       READ(ninp, *)  id, itrmx
       !


       no_s = ncore(1)
       no_p = ncore(2)
       no_d = ncore(3)
       nsp = no_s + no_p
       ntc = SUM(ncore)



       ALLOCATE( ek(1:SUM(ncore) ) )
       ALLOCATE( lk(1:SUM(ncore) ) )


        lk(       1:  no_s ) = 0
        lk( no_s +1:   nsp ) = 1
        lk( nsp + 1:   ntc ) = 2


     ENDIF hf_write_parameters


    CLOSE(ninp)

    ndim = nb - 2 + first_spline + last_spline

    xi = 0.0_dpk


!.............................. output the values
! atomic system
    WRITE(*,'(a45,G15.2)') ' # param::input:      atomic number      z = ', znuc
! basis
    WRITE(*, '(a45,G15.2)')' # param::input:             box radius  R = ', rmax
    WRITE(*, '(a45,I4)')   ' # param::input:      b-splines number  nb = ', nb
    WRITE(*, '(a45,I4)')   ' # param::input:      b-splines order   kb = ', kb
    WRITE(*, '(a45,G15.2)')' # param::input:      first knot point  rs = ', rs
!!%    IF(idbsp==0) THEN
!!%       WRITE(*, *)  '# param::input:             knot sequence   sine ', idbsp
!!%    ELSEIF(idbsp==1) THEN
!!%       WRITE(*, *)  '# param::input:             knot sequence linear ', idbsp
!!%    ENDIF
! spectrum information
    IF(method=='l'.OR.method=='dl') THEN
       WRITE(*,'(a45,a2)')  ' # param::input:           solution method = ', method
       WRITE(*,'(a45,a3)')  ' # param::input:                  spectrum = ', spectrum
       WRITE(*,'(a45,G15.2)')' # param::input:       initial energy   e0 = ', en_1e
       WRITE(*,'(a45,G15.2)')' # param::input:       step    energy   de = ', de_1e
       WRITE(*,'(a45,I4)')  ' # param::input:  nof bound states     nbs = ', nbs
       WRITE(*,'(a45,I4)')  ' # param::input:  nof free  states     ncs = ', ncs
    ELSE

       WRITE(*,'(a45,a2)')  ' # param::input:                fxd method = ', method
       WRITE(*,'(a45,I2)')  ' # param::input:              first_spline = ', first_spline
       WRITE(*,'(a45,I2)')  ' # param::input:               last_spline = ', last_spline
       WRITE(*, '(a45,I4)') ' # param::input:    matrix dimension  ndim = ', ndim
    ENDIF
    !potential information
    WRITE(*,'(a45,G15.2)') ' # param::input:                      mass = ', mass
    !hartree-fock parameters
    IF(ihf.EQ.1) THEN
       WRITE(*, '(a45,3I4)')   ' # param::input:                    lcore = ', lcore
       WRITE(*, '(a45,3I4)')   ' # param::input:              nos,nop,nod = ', ncore
       WRITE(*, '(a45,3I4)')   ' # param::input:                       lk = ', lk
       WRITE(*,'(a45,3f10.5)') ' # param::input:  effective charges    zk = ', zk

! basis
       WRITE(*, '(a45,3G15.2,2X,I5)')' # param::input:   crt, di, df, id  = ', crt, di, df, id
       WRITE(*, '(a45,3I4)')         ' # param::input:              itrmx = ', itrmx
    endif


  END SUBROUTINE input

END MODULE param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
