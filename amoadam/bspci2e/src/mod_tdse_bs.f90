!
MODULE deriv_tdse_fxd
  USE PRECISION
  IMPLICIT NONE
  REAL(DPK), POINTER, DIMENSION(:) :: yderiv
END MODULE deriv_tdse_fxd

!...

MODULE parameter_tdse_fxd
  !
  USE PRECISION
  !
  IMPLICIT NONE
  !
  PUBLIC
  !common to all tdse*.f90
  CHARACTER( LEN =  30 )                   :: atomic_name   !input file
  INTEGER                                  :: lmax          !input file
  INTEGER                                  :: nmax          !input file
  REAL(dpk)                                :: en_ground     !input file
  REAL(dpk)                                :: en_ion_1      !input file
  REAL(dpk)                                :: en_cut        !input file
  CHARACTER(len=6)                         :: gauge         !input file
  CHARACTER(len=6)                         :: bin_type      !input file
  LOGICAL                                  :: hohg          !input file
  LOGICAL                                  :: ac            !input file
  INTEGER                                  :: irestart      !input file
  REAL(dpk)                                :: tol           !input file
!  REAL(dpk)                                :: td            !ac delay time: command line
  !
  INTEGER                                  :: neq, ntot
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:) :: dzr, dr, da
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: en
  REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: dag
  INTEGER,   ALLOCATABLE, DIMENSION(:)     :: n, nsum

  !
  !fcn 
  INTEGER                                  :: neval
  !
  !

CONTAINS
  !
  SUBROUTINE input_tdse_fxd
    !
    use io
    !
    IMPLICIT NONE
    !

!    CALL getarg(1, argv)                ! peak intensity
!    READ(argv,*)   td                   ! in SI 


    OPEN(ninp,file="tinp/tdse_bs.inp")
    READ(ninp, * )  atomic_name              ! he,mg,,,,
    READ(ninp, * )  bin_type                 ! use netCDF (cdf) or f90(f90) binaries
    READ(ninp, * )  Lmax                     ! L = 0,1,2,. lmax
    READ(ninp, * )  nmax                     ! Max no of States   in L = 0,1,2,. lmax+1-1
    READ(ninp, * )  en_ion_1                 ! E_g(+) (a.u.)
    READ(ninp, * )  en_ground                ! E_g    (a.u.)
    READ(ninp, * )  en_cut                   ! E_c    (a.u.)
    READ(ninp, * )  gauge                    ! gauge used 'l', 'v'
    READ(ninp, * )  hohg                     ! true or false 
    READ(ninp, * )  irestart                 ! restart
    READ(ninp, * )  tol                      ! tolerance
  CLOSE(ninp)
  !
  END SUBROUTINE input_tdse_fxd
  !
  !
  !
  SUBROUTINE output_tdse_fxd
    !
    use io
    !
    IMPLICIT NONE
    !
    !    OPEN(nout,file="tout/tdse_bs.out")
    WRITE(*, '(a1,a60,a10)'   )'&', ' system = ', atomic_name     ! he,mg,,,,
    WRITE(*, '(a1,a60,i10)'   )'&', '   Lmax = ', Lmax            ! L = 0,1,2,. lmax
    WRITE(*, '(a1,a60,i10)'   )'&', '   nmax = ', nmax            ! Max no of States   in L = 0,1,2,. lmax+1-1
    WRITE(*, '(a1,a60,E15.5)' )'&', '     E+ = ', en_ion_1        ! E_g(+) (a.u.)
    WRITE(*, '(a1,a60,E15.5)' )'&', '      E = ', en_ground       ! E_g    (a.u.)
    WRITE(*, '(a1,a60,E15.5)' )'&', '     Ec = ', en_cut          ! E_c    (a.u.)
    WRITE(*, '(a1,a60,a10)'   )'&', '  gauge = ', gauge           ! gauge used 'l', 'v'
    WRITE(*, '(a1,a60,i10)'   )'&', '  start = ', irestart        ! restart
    WRITE(*, '(a1,a60,E15.8)' )'&', '    tol = ', tol             ! tolerance
!    WRITE(*, '(a1,a60,E15.8)' )'&', '     td = ', td              ! tolerance
!    WRITE(nout, '(a60,i10)'   )  hohg                    ! true or false 
  CLOSE(ninp)
END SUBROUTINE output_tdse_fxd
  !
  !
  !

SUBROUTINE make_space
    !exe!
    !LCL
    INTEGER       :: neq_max
    !exe!

    neq_max = 2 * nmax * ( lmax + 1 ) 

    ALLOCATE(    n( 1:lmax + 1  )                )
    ALLOCATE( nsum( 1:lmax + 1  )                )
    ALLOCATE(  dag( 1:neq_max   )                )
    ALLOCATE(   en( 1:nmax, 1:lmax + 1 )         )

    !
    ! matrix elements
    !

    IF(gauge=="v") THEN

       ALLOCATE(  dzr( 1:nmax, 1:nmax, 1:lmax )  )

       IF(hohg) THEN 
          ALLOCATE(  dr( 1:nmax, 1:nmax, 1:lmax )  )
          !ALLOCATE(  da( 1:nmax, 1:nmax, 1:lmax )  )
       ENDIF


    ELSE IF (gauge=="l") THEN

       ALLOCATE(  dr( 1:nmax, 1:nmax, 1:lmax )     )

       IF(hohg) THEN 
          ALLOCATE(  dzr( 1:nmax, 1:nmax, 1:lmax )  )
          !ALLOCATE(  da( 1:nmax, 1:nmax, 1:lmax )  )
       ENDIF

    ELSE IF  (gauge=="lv") THEN

       ALLOCATE(  dr( 1:nmax, 1:nmax, 1:lmax )     ) 

       IF(hohg) THEN 
          ALLOCATE(  dzr( 1:nmax, 1:nmax, 1:lmax )  )
          !ALLOCATE(  da( 1:nmax, 1:nmax, 1:lmax )  )
       ENDIF

    ENDIF

    
  END SUBROUTINE make_space

  !
  ! reads the field-free atom dynamical quantities: 
  !
  !    E_{nL} = < nL| H_0 | nL >        Energies (Ryd) 
  ! 
  !    V_{nL,mL+1} = <nl| V | m L+1>    dipole matrix elements
  !
  !    form =       length (l), 
  !               velocity (v), 
  !           acceleration (a)
  !
  !

  SUBROUTINE read_system(form)
    !
    IMPLICIT NONE
    !ARG!    
    CHARACTER(len=6)            :: form
    !LCL
    INTEGER                     :: mode
    INTEGER                     :: li,lf
    INTEGER                     :: ni,nf
    INTEGER                     :: l
    !    
    INTEGER                     :: nbin = 1
    !exe

    WRITE(*,*) '# opening atomic files:'
     
    en = 0.0_dpk 
    
    IF(form.EQ.'v')     dzr = 0.0_dpk
    IF(form.EQ.'l')     dr  = 0.0_dpk
    IF(form.EQ.'lv')    dr  = 0.0_dpk
    IF(form.EQ.'a')     da  = 0.0_dpk

    open_dmx_files:DO l = 1, lmax


       CALL dmxfile(nbin, "dat",TRIM(atomic_name),"bin", form, l-1, l) 

       READ(nbin) mode
       READ(nbin) li, lf, n(l), n(l+1)
       
       IF(  n(l) > nmax )   n(l) = nmax
       IF(n(l+1) > nmax ) n(l+1) = nmax
       
       WRITE(*,'(a1,a60,i6)')  '&', 'li,lf = ', li,  lf
       WRITE(*,'(a1,a60,i6)')  '&', '   ni = ', n(l) 
       WRITE(*,'(a1,a60,i6)')  '&', '   nf = ', n(l+1) 
       
       READ(nbin) ( en(ni, li + 1 ), ni = 1, n(l)  )
       READ(nbin) ( en(nf, lf + 1 ), nf = 1, n(l + 1) )
       
       DO  ni = 1, n(l)
          IF(form.EQ.'v') THEN 
             READ(nbin) (dzr(nf,ni,l), nf = 1, n(l+1))
          ELSE IF(form.EQ.'l') THEN 
             READ(nbin) (dr(nf,ni,l),  nf = 1, n(l+1))
             !READ(nbin) (dzr(nf,ni,l),  nf = 1, n(l+1))
          ELSE IF(form.EQ.'lv') THEN 
             READ(nbin) (dr(nf,ni,l),  nf = 1, n(l+1))
             !READ(nbin) (dzr(nf,ni,l),  nf = 1, n(l+1))
             !      ELSE IF(form.EQ.'a') THEN 
             !READ(nbin) (da(nf,ni,l),  nf = 1, n(l+1))
          ENDIF         
       END DO
       CLOSE(nbin)
       
       WRITE(*,'(a1,a60,i6)') '&','sizeof  e_i = ', SIZE( en(:,l  ) )
       WRITE(*,'(a1,a60,i6)') '&','sizeof  e_f = ', SIZE( en(:,l+1) )


       
    ENDDO open_dmx_files

    WRITE(*,*) '& atomic files read done.:'    
    
    !    en(1,1) = 2.0_dpk * en_ground ! set to the real ground state (in Ryd)

    WRITE(*,*)'& energies: (in a.u.)'
    WRITE(*,*)'& dipoles: (ni = 1), nf = 1,2,.. '
    DO l = 0, lmax
       WRITE(*,*)'& li = ', l
       WRITE(*,'(8E20.8)') (  en( ni, l+1  ) / 2.0_dpk,  ni = 1, MIN(nmax,4))

          IF(form.EQ.'v') THEN 
             WRITE(*,*) 'ok', l
             DO nf = 1, MIN(nmax,4) 
                WRITE(*,'(8E20.8)')  dzr( nf, 1, l+1 )
!                WRITE(*,'(8E20.8)') ( dzr( nf, 1, l ),            nf = 1, MIN(nmax,4))
             ENDDO
          ELSE IF(form.EQ.'l') THEN 
             WRITE(*,'(8E20.8)') ( dr( nf, 1, l+1 ),            nf = 1, MIN(nmax,4))
!            WRITE(*,'(8E20.8)') ( dzr( nf, 1, l ),            nf = 1, MIN(nmax,4))
          ELSE IF(form.EQ.'lv') THEN 
             WRITE(*,'(8E20.8)') ( dr( nf, 1, l+1 ),            nf = 1, MIN(nmax,4))
!            WRITE(*,'(8E20.8)') ( dzr( nf, 1, l ),            nf = 1, MIN(nmax,4))
          ELSE IF(form.EQ.'a') THEN 
             WRITE(*,'(8E20.8)') ( da( nf, 1, l+1 ),            nf = 1, MIN(nmax,4))
          ENDIF         
    ENDDO
       
  END SUBROUTINE read_system

  !
  !
  !

  SUBROUTINE read_system_ncFile(form)
    !
    USE netcdf         ! netcdf module
    USE ncFile         ! my module
    USE PRECISION
    USE units
    USE io
    !
    IMPLICIT NONE
    !
    !ARG!    
    CHARACTER(len=6)                        :: form         !input file
    !
    REAL(dpk), ALLOCATABLE,  DIMENSION(:,:) :: d_if       ! ni,nf    d(l,l+1)
    REAL(dpk), ALLOCATABLE,    DIMENSION(:) :: e_i        ! e_i(n_i) 
    REAL(dpk), ALLOCATABLE,    DIMENSION(:) :: e_f        ! e_f(n_f)
    !
    CHARACTER( LEN = 100 )                  :: filename
    CHARACTER( LEN =  15 )                  :: sl1,sl2
    INTEGER                                 :: l
    !EXE!

    WRITE(*,*) '& opening atomic files:', atomic_name

    
    ALLOCATE( e_i( 1:nmax))
    ALLOCATE( e_f( 1:nmax))
    ALLOCATE(d_if( 1:nmax, 1:nmax))

    loop_over_pws:DO l = 0, lmax-1

       WRITE(sl1,'(I15)') l
       WRITE(sl2,'(I15)') l+1

     filename  = "dat/"//trim(adjustl(atomic_name)) & 
                       //TRIM(ADJUSTL(sl1))//TRIM(ADJUSTL(sl2)) &
                       //TRIM(ADJUSTL(form))                    &  
                       //".nc"
 
     IF(atomic_name == "he") THEN
        CALL read_atomic_ncFile(filename, l, e_i, e_f, d_if)
     ELSE IF(atomic_name == "h") THEN
        CALL read_atomic_ncFile_1e(filename, l, e_i, e_f, d_if)
     ENDIF
     
     n(l+1)           = nmax
     n(l+2 )          = nmax
     en(1:nmax,l+1)   = e_i(1:nmax)
     en(1:nmax,l+2)   = e_f(1:nmax)
     IF     (form.EQ.'v') THEN 
        dzr(1:nmax,1:nmax, l+1) = d_if(1:nmax, 1:nmax)
     ELSE IF(form.EQ.'l') THEN 
        dr(1:nmax,1:nmax, l+1) = d_if(1:nmax, 1:nmax)
     ELSE IF(form.EQ.'a') THEN 
        da(1:nmax,1:nmax, l+1) = d_if(1:nmax, 1:nmax)
     ENDIF

  END DO loop_over_pws  

  !initial version
  IF(atomic_name == "he") THEN
     en = 2.0_dpk * en ! transform to Ryd
  ENDIF
  !

  !  en(1,1) = 2.0_dpk * en_ground ! set to the real ground state (in Ryd)

  DEALLOCATE( e_i)
  DEALLOCATE( e_f)
  DEALLOCATE(d_if)

  !
  END SUBROUTINE read_system_ncFile
!
!
!
!
END MODULE PARAMETER_TDSE_FXD
!eof
