!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE param

  USE PRECISION, ONLY : DPK
  IMPLICIT NONE
  PUBLIC

  INTEGER, PARAMETER :: kx = 15, np = 1200
  INTEGER            :: ntc, nsp, id
!
!..........  h1e.inp file
!
  CHARACTER(LEN=10)  :: problem
  REAL(DPK)          :: znuc 
  REAL(DPK)          :: za, ma, mass
  REAL(DPK)          :: zb, mb
  REAL(DPK)          :: r_m                ! internuclear distance
  character(len=20)  :: potential
  REAL(dpk)          :: rmin, rmax, xi, rs
  INTEGER            :: no, idbsp, nrp 
  INTEGER            :: nos,  nop,     nod
  INTEGER            :: ihf, imodel, idl 
  INTEGER            :: lmin, lcore, itrmx
  INTEGER            :: l1,l2
  REAL(dpk)          :: ap, rc
  REAL(dpk), DIMENSION(:),ALLOCATABLE :: zk
  REAL(dpk), dimension(3):: ncore_orbitals
  INTEGER            :: NB, KB, IB
  CHARACTER(LEN=6)   :: method, spectrum, code
  REAL(dpk)          :: en_1e, de_1e
  INTEGER            :: nbs, ncs
  REAL(dpk)          :: b_f                          ! magnetic field
  INTEGER            :: m_l,  l0                   ! azimuthial q. number, (0(1),2(3),4(5),....),
  INTEGER            :: first_spline, last_spline
  INTEGER            :: ndim, nl, nldim
  CHARACTER(LEN=3)   :: print_out


CONTAINS

!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  SUBROUTINE input(h2p_input_file)

    IMPLICIT NONE
    CHARACTER(len=100)  :: h2p_input_file
    INTEGER             :: ninp

!!!...............

    ninp = 9
!    OPEN(ninp,file ='inp/h2p.inp',status='old')
    OPEN(ninp, file =h2p_input_file,status='old')

! problem type: atomic, molecular, magnetic, qdot


! 'hydrogenic' atomic parameters
    READ(ninp, *) znuc, ma                 ! Atomic number, mass
    READ(ninp, *) potential, ap, rc      ! potential, static polarizability, cutoff radius
! basis parameters
    READ(ninp, *) nb, kb                 ! number of B-splines, order of B-splines 
    READ(ninp, *) rmin, rmax, no         ! interval (rmin, rmax), no number of knot points
    READ(ninp, *) idbsp, rs, nrp         ! knot sequence selects the grid 0 = sin, 1=exp
!    READ(ninp, *) idbsp, xi, nrp        ! knot sequence selects the grid 0 = sin, 1=exp
!    READ(ninp, *) rs                    ! sine-like first non-zero knot points
!spectrum information
    READ(ninp, *) first_spline, last_spline ! 1st(last)_spline(0/1) is (excluded/included),
    READ(ninp, *) method, spectrum       ! (lin, diag, inv), (fxd,free,mxd)
    READ(ninp, *) en_1e, de_1e, ncs, nbs 
    READ(ninp, *) print_out, code
!problem type 
    READ(ninp, *) problem             ! problem type (1e, 2e, h2p, 1eb, qdot)  
    READ(ninp, *) r_m, m_l            ! molecular hydrogen  

    za = znuc
    zb = znuc
    mb = ma                           ! diatomic potential       
    mass = ma

!    IF(problem == "1e".or.problem=='2e') THEN 
!       nl   = 1
!       m_l  = 0 
!    ELSE IF(problem =="1eb") THEN
!       READ(ninp, *) b_f, nl, m_l        ! magnetic field (in a.u.)
!    ELSE IF(problem =="h2p") THEN
!    ENDIF

    CLOSE(ninp)
    
    

    ndim = nb - 2 + first_spline + last_spline

    xi = 0.0_dpk


    !    nldim = nl*ndim

    !    ndim = nb - 1 
    !    IF(method=="d") THEN 
    !       ndim = nb - 2
    !    ENDIF
    !    IF(method=="dl") THEN 
    !       ndim = nb - 2
    !    ENDIF
    
    




!.............................. output the values
    WRITE(*,'("# param::input:               problem type = ", a5)'), problem

! atomic system
    WRITE(*,'(a45,G8.2)') ' # param::input:        atomic number,  za = ', za
    WRITE(*,'(a45,G8.2)') ' # param::input:          atomic mass,  ma = ', ma

    IF(problem == 'h2p') THEN
       WRITE(*,'(a45,G8.2)') ' # param::input:        atomic number,  zb = ', zb
       WRITE(*,'(a45,G8.2)') ' # param::input:          atomic mass,  ma = ', mb
    ELSEIF(problem == '1eb') THEN
       WRITE(*,'(a45,G8.2)') ' # param::input:       magnetic field  b_f = ', b_f
    ENDIF

    WRITE(*,'(a45,I2)')   ' # param::input:  azimuthial q. number m_l = ', m_l
    WRITE(*,'(a45,I2)')   ' # param::input:                        l0 = ', l0
! basis
    WRITE(*, '(a45,G8.2)')' # param::input:             box radius  R = ', rmax
    WRITE(*, '(a45,I4)')  ' # param::input:      b-splines number  nb = ', nb
    WRITE(*, '(a45,I4)')  ' # param::input:      b-splines order   kb = ', kb
    WRITE(*, '(a45,G8.2)')' # param::input:      first knot point  rs = ', rs
    IF(idbsp==0) THEN 
    WRITE(*,'(a45,a5)')   ' # param::input:      knot sequence  idbsp = ','sine'
    ELSEIF(idbsp==1) THEN
    WRITE(*,'(a45,a5)')   ' # param::input:      knot sequence  idbsp = ','exp'
    ENDIF
! spectrum information
    IF(method=='l'.OR.method=='dl') THEN 
       WRITE(*,'(a45,a2)')  ' # param::input:           solution method = ', method
       WRITE(*,'(a45,a3)')  ' # param::input:                  spectrum = ', spectrum
       WRITE(*,'(a45,G8.2)')' # param::input:       initial energy   e0 = ', en_1e
       WRITE(*,'(a45,G8.2)')' # param::input:       step    energy   de = ', de_1e
       WRITE(*,'(a45,I4)')  ' # param::input:  nof bound states     nbs = ', nbs
       WRITE(*,'(a45,I4)')  ' # param::input:  nof free  states     ncs = ', ncs
    ELSE
       WRITE(*,'(a45,a2)')  ' # param::input:           solution method = ', method

       WRITE(*, *)  '# param::input:         fixed BC are used' 
    ENDIF
!potential information

    IF(ma==1) THEN
       WRITE(*,*)    '# param::input: hydrogenic radial functions calculation'
    ELSE IF(ma==2) THEN
       WRITE(*,*) '# param::input: positronium radial functions calculated.'
    ELSE 
       WRITE(*,*) '# param::input: other radial functions produced. Maybe vibrational functions?'
    ENDIF

    WRITE(*,'(a45,a3)')   ' # param::input:                     print = ', print_out

  END SUBROUTINE input

END MODULE param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



