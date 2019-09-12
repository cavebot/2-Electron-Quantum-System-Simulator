!
! contains: write_ncFile, read_ncFile (prototypes)
!           check
!        

! write_atomic_ncFile(filename, system, gauge, l, e_i, e_f, d_if)    
! read_atomic_ncFile(filename, l, e_i, e_f, d_if)
! read_atomic_ncFile_1e(filename, l, e_i, e_f, d_if)
! open_tdse_ncFile(system, gauge, ncid)
! save_pulse_ncFile(ncid, intensity, frequency, duration)
! save_en_ncFile(ncid, n_l, e_nl)
! save_ct_ncFile(ncid, t, c_nl, t_name, c_name)
! close_tdse_ncFile(ncid)
!
!
! prototype functions
!
! write_ncFile(datafile, dim_s, charact, iscalar, rscalar, dim_v, ivector, rvector, rmatrix)
! read_ncFile(datafile, charact, iscalar, rscalar, ivector, rvector, rmatrix)
! check(status)



MODULE ncFile
  !
  USE PRECISION, ONLY: dpk
  !
  IMPLICIT NONE
  PUBLIC  write_ncFile

CONTAINS
  !
  !
  !
  SUBROUTINE write_atomic_ncFile(filename, system, gauge, l, e_i, e_f, d_if)    
    !
    USE netcdf
    !
    IMPLICIT NONE
    !
    CHARACTER (len = 100),          PARAMETER :: units='a.u.'      !dipole form
    !variable_name
    CHARACTER (len = *), PARAMETER ::        l_i_name = "l_i"
    CHARACTER (len = *), PARAMETER ::        l_f_name = "l_f"
    CHARACTER (len = *), PARAMETER ::        e_i_name = "e_i"
    CHARACTER (len = *), PARAMETER ::        e_f_name = "e_f"
    CHARACTER (len = *), PARAMETER ::       d_if_name = "d_if"
    !dimension_name
    CHARACTER (len = *), PARAMETER ::    e_i_dim_name =    "e_i_dim"
    CHARACTER (len = *), PARAMETER ::    e_f_dim_name =    "e_f_dim"
    CHARACTER (len = *), PARAMETER :: d_if_dim_1_name = "d_if_dim_1"
    CHARACTER (len = *), PARAMETER :: d_if_dim_2_name = "d_if_dim_2"

    ! attributes
    !    CHARACTER (len = *), PARAMETER ::      att_units  = "units"
    !    CHARACTER (len = *), PARAMETER ::  e_i_att_units  = "a.u."
    !    CHARACTER (len = *), PARAMETER ::  e_f_att_units  = "a.u."
    !    CHARACTER (len = *), PARAMETER :: d_if_att_units  = "a.u."

    CHARACTER (len = *), PARAMETER ::      att_state  = "state"
    CHARACTER (len = *), PARAMETER ::  l_i_att_state  = "initial"
    CHARACTER (len = *), PARAMETER ::  l_f_att_state  = "final"

    ! file location
    CHARACTER (len = 100),               INTENT(in) :: filename   !ncFile location
    ! variables stored 
    CHARACTER (len = 30),                INTENT(in) :: system     !atomic system
    CHARACTER (len = 30),                INTENT(in) :: gauge      !dipole form
    !
    INTEGER,                             INTENT(in) :: l
    REAL(dpk),   DIMENSION(:),  POINTER, INTENT(in) :: e_i
    REAL(dpk),   DIMENSION(:),  POINTER, INTENT(in) :: e_f
    REAL(dpk), DIMENSION(:,:),  POINTER, INTENT(in) :: d_if

    !file id
    INTEGER                                       :: ncid        
    !variable_id
    INTEGER                                       :: l_i_id
    INTEGER                                       :: l_f_id
    INTEGER                                       :: e_i_id
    INTEGER                                       :: e_f_id
    INTEGER                                       :: d_if_id
    !variable_dim
    INTEGER                                       ::  e_i_dim_id
    INTEGER                                       ::  e_f_dim_id
    INTEGER, DIMENSION(2)                         :: d_if_dim_id
    !
    !
    INTEGER                                       :: dim_i, dim_f
    !
    

    !set dimensions
    dim_i = SIZE(e_i)
    dim_f = SIZE(e_f)

!    dim_s = SIZE(system)
!    dim_g = SIZE(gauge)
    
    CALL check( nf90_create(trim(filename), nf90_clobber, ncid) )      ! Create the file. 

    WRITE(*,*) '# create_atomic_ncFile:: file opened: ', filename

    !define dimensions.
    CALL check( nf90_def_dim( ncid,    e_i_dim_name,   dim_i,     e_i_dim_id  ) )    ! 
    CALL check( nf90_def_dim( ncid,    e_f_dim_name,   dim_f,     e_f_dim_id  ) )    !
    CALL check( nf90_def_dim( ncid, d_if_dim_1_name,   dim_i,  d_if_dim_id(1) ) )    ! 
    CALL check( nf90_def_dim( ncid, d_if_dim_2_name,   dim_f,  d_if_dim_id(2) ) )    !

    
    WRITE(*,*) '# read_ncFile:: dimensions given'

    !define variables
    CALL check( nf90_def_var(ncid,    l_i_name,  nf90_int,                   l_i_id  ) )
    CALL check( nf90_def_var(ncid,    l_f_name,  nf90_int,                   l_f_id  ) )
    CALL check( nf90_def_var(ncid,    e_i_name, nf90_real,    e_i_dim_id,    e_i_id  ) )
    CALL check( nf90_def_var(ncid,    e_f_name, nf90_real,    e_f_dim_id,    e_f_id  ) )
    CALL check( nf90_def_var(ncid,   d_if_name, nf90_real,   d_if_dim_id,   d_if_id  ) )
    
    WRITE(*,*) '# read_ncFile:: variables defined'

    !attributes
    !    CALL check( nf90_put_att(ncid, charact_id, att_units, charact_att_units) )
    !    CALL check( nf90_put_att(ncid,  e_i_id, att_units,  e_i_att_units) )
    !    CALL check( nf90_put_att(ncid,  e_f_id, att_units,  e_f_att_units) )
    !    CALL check( nf90_put_att(ncid, d_if_id, att_units, d_if_att_units) )
    !
    CALL check( nf90_put_att(ncid, l_i_id, att_state, l_i_att_state) )
    CALL check( nf90_put_att(ncid, l_f_id, att_state, l_f_att_state) )

    
    WRITE(*,*) '# read_ncFile:: attributes given '

    ! Global attributes
    CALL check(nf90_put_att(ncid, nf90_global, "atomic_ncFile", &
          &   "created by l.aa. nikolopoulos OCT2008QUB"))
    CALL check(nf90_put_att(ncid, nf90_global, "system", system ))
    CALL check(nf90_put_att(ncid, nf90_global,  "gauge",  gauge ))
    CALL check(nf90_put_att(ncid, nf90_global,  "units",  units ))

    WRITE(*,*) '# read_ncFile:: global identification set '

    CALL check( nf90_enddef(ncid) )         ! end of defining mode
    
    
    ! now write in ncFile
    


    CALL check( nf90_put_var( ncid,  l_i_id, l   ) )     ! l_i
    CALL check( nf90_put_var( ncid,  l_f_id, l+1 ) )     ! l_f
    CALL check( nf90_put_var( ncid,  e_i_id, e_i ) )     ! e_i
    CALL check( nf90_put_var( ncid,  e_f_id, e_f ) )     ! e_f
    CALL check( nf90_put_var( ncid, d_if_id, d_if) )     ! d_if

    ! close the ncFile

    CALL check( nf90_close(ncid) )                            
    
    PRINT *,"# write_netCDF done!"
    
  END SUBROUTINE write_atomic_ncFile
  !
  !

  !prototype
  SUBROUTINE read_atomic_ncFile(filename, l, e_i, e_f, d_if)
    !
    USE netcdf
    !
    IMPLICIT NONE    
    !variable_name
    CHARACTER (len = *), PARAMETER ::        l_i_name = "l_i"
    CHARACTER (len = *), PARAMETER ::        l_f_name = "l_f"
    CHARACTER (len = *), PARAMETER ::        e_i_name = "e_i"
    CHARACTER (len = *), PARAMETER ::        e_f_name = "e_f"
    CHARACTER (len = *), PARAMETER ::       d_if_name = "d_if"

    !dimension_name
    CHARACTER (len = *), PARAMETER ::    e_i_dim_name =    "e_i_dim"
    CHARACTER (len = *), PARAMETER ::    e_f_dim_name =    "e_f_dim"
    CHARACTER (len = *), PARAMETER :: d_if_dim_1_name = "d_if_dim_1"
    CHARACTER (len = *), PARAMETER :: d_if_dim_2_name = "d_if_dim_2"

    ! file location
    CHARACTER (len = 100),                   INTENT(in) :: filename   !ncFile location

    ! variables stored 
    INTEGER,                                 INTENT(in) :: l
    REAL(dpk), ALLOCATABLE,   DIMENSION(:), INTENT(out) :: e_i
    REAL(dpk), ALLOCATABLE,   DIMENSION(:), INTENT(out) :: e_f
    REAL(dpk), ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: d_if
    !LCL
    INTEGER                                             :: l_i,l_f
    
    !file id
    INTEGER                                         :: ncid        
    !variable_id 
    INTEGER                                         :: l_i_id
    INTEGER                                         :: l_f_id

    INTEGER                                         :: e_i_id
    INTEGER                                         :: e_f_id
    INTEGER                                         :: d_if_id
    !variable_dim
    INTEGER                                         ::  e_i_dim_id
    INTEGER                                         ::  e_f_dim_id
    INTEGER, DIMENSION(2)                           :: d_if_dim_id
    !
    INTEGER                                         :: dim_i, dim_f
    INTEGER                                         :: ndims_in,  nvars_in
    INTEGER                                         :: ngatts_in, unlimdimid_in
    !exe!

    !
 
    WRITE(*,*) 'check ', trim(filename)
    ! open ncFile
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )
    ! inquire info
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )
    ! read dims
    CALL check( nf90_inq_dimid(ncid,    e_i_dim_name,    e_i_dim_id  ) )
    CALL check( nf90_inq_dimid(ncid,    e_f_dim_name,    e_f_dim_id  ) )
    CALL check( nf90_inq_dimid(ncid, d_if_dim_1_name, d_if_dim_id(1) ) )
    CALL check( nf90_inq_dimid(ncid, d_if_dim_2_name, d_if_dim_id(2) ) )
    ! now read  ncFile
    ! inquire variable_id from variable_name
    CALL check( nf90_inq_varid( ncid,  l_i_name,  l_i_id )  )     ! read var scalar
    CALL check( nf90_inq_varid( ncid,  l_f_name,  l_f_id )  )     ! read var scalar
    CALL check( nf90_inq_varid( ncid,  e_i_name,  e_i_id )  )     ! read var vector
    CALL check( nf90_inq_varid( ncid,  e_f_name,  e_f_id )  )     ! read var vector
    CALL check( nf90_inq_varid( ncid, d_if_name, d_if_id )  )     ! read var matrix
    !
    CALL check( nf90_inquire_dimension(ncid,    e_i_dim_id,  len=dim_i) )
    CALL check( nf90_inquire_dimension(ncid,    e_f_dim_id,  len=dim_f) )
    CALL check( nf90_inquire_dimension(ncid, d_if_dim_id(1), len=dim_i) )
    CALL check( nf90_inquire_dimension(ncid, d_if_dim_id(2), len=dim_f) )
    !
    PRINT*, '&  e_i:      dim_i = ', dim_i
    PRINT*, '&  e_f:      dim_f = ', dim_f
    !
    ALLOCATE(  e_i( dim_i) )
    ALLOCATE(  e_f( dim_f) )
    ALLOCATE( d_if( dim_i, dim_f))

    ! inquire a
    !    CALL check( nf90_get_att(ncid, charact_id, att_units, charact_att_units))     ! character
    !
    !

    CALL check( nf90_get_var( ncid,  l_i_id,   l_i  )  )
    CALL check( nf90_get_var( ncid,  l_f_id,   l_f  )  )

    IF((l_i.NE.l).OR.(l_f.NE.l+1)) THEN 
       WRITE(*,*) '# read_atomic_file::         inconsistent file'
       WRITE(*,*) '# read_atomic_file::                      l = ',l
       WRITE(*,*) '# read_atomic_file::                    l_i = ',l_i
       WRITE(*,*) '# read_atomic_file::                    l_f = ',l_f
    ENDIF
       
    CALL check( nf90_get_var( ncid,  e_i_id,   e_i  )  )
    CALL check( nf90_get_var( ncid,  e_f_id,   e_f  )  )
    CALL check( nf90_get_var( ncid, d_if_id,   d_if )  )

    !
    CALL check( nf90_close(ncid) )                            ! close the file.
    
    PRINT *,"# read_ncFile done!", l_i, l_f    

  END SUBROUTINE read_atomic_ncFile
  !
  ! get data from 1-e code (c++)
  !
  SUBROUTINE read_atomic_ncFile_1e(filename, l, e_i, e_f, d_if)
    !
    USE netcdf
    !
    IMPLICIT NONE
!    CHARACTER (len = 100), PARAMETER          :: units='a.u.'      !dipole form    
    !variable_name
    CHARACTER (len = *), PARAMETER ::        l_i_name = "l_i"
    CHARACTER (len = *), PARAMETER ::        e_i_name = "e_i"
    CHARACTER (len = *), PARAMETER ::        e_f_name = "e_f"
    CHARACTER (len = *), PARAMETER ::       d_if_name = "d_if"
    !dimension_name
    CHARACTER (len = *), PARAMETER ::    e_i_dim_name =    "e_i_dim"
    CHARACTER (len = *), PARAMETER ::    e_f_dim_name =    "e_f_dim"
    
    !   CHARACTER (len = *), PARAMETER :: d_if_dim_1_name = "d_if_dim_1"
    !   CHARACTER (len = *), PARAMETER :: d_if_dim_2_name = "d_if_dim_2"    
    ! attributes
    !    CHARACTER (len = *), PARAMETER ::      att_units  = "units"
    !    CHARACTER (len = *), PARAMETER ::  e_i_att_units  = "a.u."
    !    CHARACTER (len = *), PARAMETER ::  e_f_att_units  = "a.u."
    !    CHARACTER (len = *), PARAMETER :: d_if_att_units  = "a.u."

!    CHARACTER (len = *), PARAMETER ::      att_state  = "state"
!    CHARACTER (len = *), PARAMETER ::  l_i_att_state  = "initial"


    ! file location
    CHARACTER (len = 100),                   INTENT(in) :: filename   !ncFile location

    ! variables stored 
    INTEGER,                                 INTENT(in) :: l
    !
    INTEGER                                             :: l_i
    REAL(dpk), ALLOCATABLE,   DIMENSION(:), INTENT(out) :: e_i
    REAL(dpk), ALLOCATABLE,   DIMENSION(:), INTENT(out) :: e_f
    REAL(dpk), ALLOCATABLE, DIMENSION(:,:), INTENT(out) :: d_if
    !LCL

    !file id
    INTEGER                                         :: ncid        
    !variable_id 
    INTEGER                                         :: l_i_id
    INTEGER                                         :: e_i_id
    INTEGER                                         :: e_f_id
    INTEGER                                         :: d_if_id
    !variable_dim
    INTEGER                                         ::  e_i_dim_id
    INTEGER                                         ::  e_f_dim_id
    !  INTEGER, DIMENSION(2)                           :: d_if_dim_id
    !
    INTEGER                                         :: dim_i, dim_f
    INTEGER                                         :: ndims_in,  nvars_in
    INTEGER                                         :: ngatts_in, unlimdimid_in
    !exe!

    !
 
    WRITE(*,*) 'check 1e: ', trim(filename)


    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                    ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) ) ! inquire info
    ! read dimension_id from dimension_name
    CALL check( nf90_inq_dimid(ncid,    e_i_dim_name,    e_i_dim_id  ) )
    CALL check( nf90_inq_dimid(ncid,    e_f_dim_name,    e_f_dim_id  ) )

!    CALL check( nf90_inq_dimid(ncid, d_if_dim_1_name, d_if_dim_id(1) ) )
!    CALL check( nf90_inq_dimid(ncid, d_if_dim_2_name, d_if_dim_id(2) ) )

    ! now read  ncFile
    ! inquire variable_id from variable_name
    CALL check( nf90_inq_varid( ncid,  l_i_name,  l_i_id )  )     ! read var scalar
    CALL check( nf90_inq_varid( ncid,  e_i_name,  e_i_id )  )     ! read var vector
    CALL check( nf90_inq_varid( ncid,  e_f_name,  e_f_id )  )     ! read var vector
    CALL check( nf90_inq_varid( ncid, d_if_name, d_if_id )  )     ! read var matrix
    !inquire dimensions of variables
    CALL check( nf90_inquire_dimension(ncid, e_i_dim_id,  len=dim_i) )
    CALL check( nf90_inquire_dimension(ncid, e_f_dim_id,  len=dim_f) )

    !   CALL check( nf90_inquire_dimension(ncid, d_if_dim_id(1), len=dim_i) )
    !   CALL check( nf90_inquire_dimension(ncid, d_if_dim_id(2), len=dim_f) )


    PRINT*, '&  e_i:      dim_i = ', dim_i
    PRINT*, '&  e_f:      dim_f = ', dim_f
 
    ALLOCATE(  e_i( dim_i) )
    ALLOCATE(  e_f( dim_f) )
    ALLOCATE( d_if( dim_i, dim_f))


    CALL check( nf90_get_var( ncid,  l_i_id,   l_i  ) )

    IF((l_i.NE.l)) THEN 
       WRITE(*,*) '# read_atomic_file::         inconsistent file'
       WRITE(*,*) '# read_atomic_file::                      l = ',l
       WRITE(*,*) '# read_atomic_file::                    l_i = ',l_i
    ENDIF
       
    CALL check( nf90_get_var( ncid,  e_i_id,   e_i  )  )
    CALL check( nf90_get_var( ncid,  e_f_id,   e_f  )  )
    CALL check( nf90_get_var( ncid, d_if_id,   d_if )  )

    !
    CALL check( nf90_close(ncid) )                            ! close the file.
    
    PRINT *,"# read_ncFile done!", l_i

  END SUBROUTINE read_atomic_ncFile_1e
  !
  !
  ! routines for the tdse_bs program
  !
  !
  SUBROUTINE open_tdse_ncFile(system, gauge, ncid)
    !
    USE netCDF
    !
    IMPLICIT NONE
    !    
    CHARACTER (len = 100),          PARAMETER ::     units='a.u.'      !dipole form
    !
    CHARACTER (len = 30),          INTENT(in) :: system     !atomic system
    CHARACTER (len = 6),           INTENT(in) :: gauge      !dipole form
    INTEGER,                      INTENT(out) :: ncid
    !
    CHARACTER (len = 100)                     :: filename    !ncFile location    
    !exe!

    filename  = "tdat/"//TRIM(ADJUSTL(system))//TRIM(ADJUSTL(gauge))//".nc"

    CALL check( nf90_create(trim(filename), nf90_clobber, ncid) )  !Create the file.

    WRITE(*,'(a60,a30)') '# read_ncFile::  ncFile created: ', TRIM(filename) 

    ! Global attributes
    CALL check(nf90_put_att(ncid, nf90_global, "tdse_ncFile", &
          &   "created by l.aa. nikolopoulos NOV2008QUB"))
    CALL check(nf90_put_att(ncid, nf90_global,  "file",   filename ))
    CALL check(nf90_put_att(ncid, nf90_global, "system",    system ))
    CALL check(nf90_put_att(ncid, nf90_global,  "gauge",     gauge ))
    CALL check(nf90_put_att(ncid, nf90_global,  "units",     units ))

    WRITE(*,*) '# read_ncFile:: global identification set '

    !
  END SUBROUTINE open_tdse_ncFile


  
  SUBROUTINE save_pulse_ncFile(ncid, intensity, frequency, duration)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !variable_name
    CHARACTER (len = *),         PARAMETER :: i_name = "intensity"     !scalar (real)
    CHARACTER (len = *),         PARAMETER :: w_name = "frequency"     !scalar (real)
    CHARACTER (len = *),         PARAMETER :: d_name = "duration"      !scalar (real)    
    !
    INTEGER,                    INTENT(in) :: ncid       !nc file id
    REAL(dpk),                  INTENT(in) :: intensity
    REAL(dpk),                  INTENT(in) :: frequency
    REAL(dpk),                  INTENT(in) :: duration
    !variable_id
    INTEGER                                :: i_id
    INTEGER                                :: w_id
    INTEGER                                :: d_id
    !

    WRITE(*,*) '# write_en_ncFile:: file opened:'

    !define variables
    CALL check( nf90_def_var(ncid, i_name,  nf90_real, i_id  ) ) !intensity
    CALL check( nf90_def_var(ncid, w_name,  nf90_real, w_id  ) ) !frequency
    CALL check( nf90_def_var(ncid, d_name,  nf90_real, d_id  ) ) !duration
    CALL check( nf90_enddef(ncid) )                               ! end of define
    WRITE(*,*) '# read_ncFile:: variables defined'
        
    ! now write in ncFile    
    CALL check( nf90_put_var( ncid,   i_id, intensity ) )     ! 
    CALL check( nf90_put_var( ncid,   w_id, frequency ) )     ! 
    CALL check( nf90_put_var( ncid,   d_id, duration  ) )     !

    ! close the ncFile
    !    CALL check( nf90_close(ncid) )
    
    PRINT *,"# pulse parameters, stored!"
    
  END SUBROUTINE save_pulse_ncFile
  !

  SUBROUTINE save_en_ncFile(ncid, n_l, e_nl)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !

    !variable_name
    CHARACTER (len = *), PARAMETER ::     l_name = "maximum_partial_wave" !scalar (int)
    CHARACTER (len = *), PARAMETER ::     n_name = "energy_index"         !vector (real)
    CHARACTER (len = *), PARAMETER ::     e_name = "energy_matrix"        !matrix (real)
    !dimension_name
    CHARACTER (len = *), PARAMETER ::   n_dim_name =  "l_dim"
    CHARACTER (len = *), PARAMETER :: e_dim_1_name =  "e_dim_1"
    CHARACTER (len = *), PARAMETER :: e_dim_2_name =  "e_dim_2"
    ! variables stored 
    INTEGER,                         INTENT(in) :: ncid

    INTEGER,   DIMENSION(:),         INTENT(in) :: n_l
    REAL(dpk), DIMENSION(:,:),       INTENT(in) :: e_nl
    INTEGER                                     :: lmax
    !variable_id
    INTEGER                                     :: l_id
    INTEGER                                     :: n_id
    INTEGER                                     :: e_id
      !variable_dim
    INTEGER                                     :: n_dim_id
    INTEGER, DIMENSION(2)                       :: e_dim_id
    !
    INTEGER                                     :: dim_l, dim_e
    
    ! 
    !
      
    lmax = SIZE(n_l) - 1

    !set dimensions
    dim_l = SIZE(n_l)
    dim_e = SIZE(e_nl,dim=1) 


    WRITE(*,*) '# save_en_ncFile:: file opened:'

    !define new variables
    
    CALL check( nf90_redef(ncid) )
    !define dimensions.
    CALL check( nf90_def_dim( ncid,    n_dim_name,   dim_l,  n_dim_id    ) )    !
    CALL check( nf90_def_dim( ncid,  e_dim_1_name,   dim_e,  e_dim_id(1) ) )    ! 
    CALL check( nf90_def_dim( ncid,  e_dim_2_name,   dim_l,  e_dim_id(2) ) )    !    
    WRITE(*,*) '# read_ncFile:: dimensions given'
    !define variables
    CALL check( nf90_def_var(ncid,   l_name, nf90_int,                l_id  ) )
    CALL check( nf90_def_var(ncid,   n_name, nf90_int,    n_dim_id,   n_id  ) )
    CALL check( nf90_def_var(ncid,   e_name, nf90_real,   e_dim_id,   e_id  ) )    
    CALL check( nf90_enddef(ncid) )         ! end of defining mode
    
    WRITE(*,*) '# read_ncFile:: new added variables defined'

    
    ! now write in ncFile
    
    CALL check( nf90_put_var( ncid,  l_id, lmax  ) )     ! l_i
    CALL check( nf90_put_var( ncid,  n_id, n_l   ) )     ! e_f
    CALL check( nf90_put_var( ncid,  e_id, e_nl  ) )     ! d_if

    ! close the ncFile
    !    CALL check( nf90_close(ncid) )            
    
    PRINT *,"# wf2e, field-free data saved!"

    
  END SUBROUTINE save_en_ncFile
  !
  !

  !
  !
  !
  !
  SUBROUTINE save_ct_ncFile(ncid, t, c_nl, t_name, c_name)
    !
      USE netcdf
    !
      IMPLICIT NONE
      !

      !variable_name

      !dimension_name
      ! variables stored 
      INTEGER,                         INTENT(in) :: ncid
      REAL(dpk),   DIMENSION(:),       INTENT(in) :: t
      REAL(dpk), DIMENSION(:,:),       INTENT(in) :: c_nl
      CHARACTER (len=2),               INTENT(in) :: c_name     !matrix (real)
      CHARACTER (len=2),               INTENT(in) :: t_name     !matrix (real)
      !variable_id
      INTEGER                                     :: t_id
      INTEGER                                     :: c_id
      !variable_dim
      INTEGER                                     :: t_dim_id
      INTEGER, DIMENSION(2)                       :: c_dim_id
      !
      INTEGER                                     :: dim_t, dim_c
      CHARACTER (len = 6)                         :: t_dim_name 
      CHARACTER (len = 8)                         :: c_dim_1_name 
      CHARACTER (len = 8)                         :: c_dim_2_name 

    !

        t_dim_name =  t_name//"_dim"
      c_dim_1_name =  c_name//"_dim_1"
      c_dim_2_name =  c_name//"_dim_2"
      
      
    ! 
      dim_t = SIZE(t)
      dim_c = SIZE(c_nl,dim=2)

    !

    !set dimensions

    WRITE(*,*) '# write_ct_ncFile:: file opened: ', dim_t, dim_c

    CALL check( nf90_redef(ncid))    !define mode.
    CALL check( nf90_def_dim( ncid,    t_dim_name,   dim_t,  t_dim_id    ) )    !
    CALL check( nf90_def_dim( ncid,  c_dim_1_name,   dim_t,  c_dim_id(1) ) )    ! 
    CALL check( nf90_def_dim( ncid,  c_dim_2_name,   dim_c,  c_dim_id(2) ) )    !
    !define variables
    CALL check( nf90_def_var(ncid,   t_name, nf90_real,   t_dim_id,   t_id  ) )
    CALL check( nf90_def_var(ncid,   c_name, nf90_real,   c_dim_id,   c_id  ) )    
    WRITE(*,*) '# read_ncFile:: variables defined'

    CALL check( nf90_enddef(ncid) )         ! end of defining mode
    
    
    ! now write in ncFile    
    CALL check( nf90_put_var( ncid,  t_id, t    ) )    ! t(i)
    CALL check( nf90_put_var( ncid,  c_id, c_nl ) )    ! c_nl(t)

    ! close the ncFile
    !    CALL check( nf90_close(ncid) )            
    
    PRINT *,"# wf2e, field-free data saved!"    
  END SUBROUTINE save_ct_ncFile
  !
  !
  ! close the tdse netcdf file
  !
  !
  SUBROUTINE close_tdse_ncFile(ncid)
    USE netCDF
    INTEGER,      INTENT(in) :: ncid
    CALL check( nf90_close(ncid) ) 
  END SUBROUTINE close_tdse_ncFile
  !





  !
  !
  !
  !prototype
  !
  !

  SUBROUTINE write_ncFile(datafile, dim_s, charact, iscalar, rscalar, dim_v, ivector, rvector, rmatrix)
    
    !
    USE netcdf
    !
    IMPLICIT NONE
    !variable_name
    CHARACTER (len = *), PARAMETER :: charact_name = "charact"
    CHARACTER (len = *), PARAMETER :: iscalar_name = "iscalar"
    CHARACTER (len = *), PARAMETER :: rscalar_name = "rscalar"
    CHARACTER (len = *), PARAMETER :: ivector_name = "ivector"
    CHARACTER (len = *), PARAMETER :: rvector_name = "rvector"
    CHARACTER (len = *), PARAMETER :: rmatrix_name = "rmatrix"
    !dimension_name
    CHARACTER (len = *), PARAMETER ::   charact_dim_name = "charact_dim"
    CHARACTER (len = *), PARAMETER ::   ivector_dim_name = "ivector_dim"
    CHARACTER (len = *), PARAMETER ::   rvector_dim_name = "rvector_dim"
    CHARACTER (len = *), PARAMETER :: rmatrix_dim_1_name = "rmatrix_dim_1"
    CHARACTER (len = *), PARAMETER :: rmatrix_dim_2_name = "rmatrix_dim_2"

    ! file location
    CHARACTER (len = *),             INTENT(in) :: datafile 
    ! variables stored in file
    INTEGER,                           INTENT(in) :: dim_s
    CHARACTER(len=dim_s),              INTENT(in) :: charact
    INTEGER,                           INTENT(in) :: iscalar
    REAL(dpk),                         INTENT(in) :: rscalar
    INTEGER,                           INTENT(in) :: dim_v
    INTEGER, DIMENSION(dim_v),         INTENT(in) :: ivector
    REAL(dpk), DIMENSION(dim_v),       INTENT(in) :: rvector
    REAL(dpk), DIMENSION(dim_v,dim_v), INTENT(in) :: rmatrix
    !variable_id
    INTEGER                                   :: ncid
    INTEGER                                   :: charact_id
    INTEGER                                   :: iscalar_id, rscalar_id
    INTEGER                                   :: ivector_id, rvector_id
    INTEGER                                   :: rmatrix_id
    !variable_dim
    INTEGER                                   :: charact_dim_id
    INTEGER                                   :: ivector_dim_id,  rvector_dim_id
    INTEGER, DIMENSION(2)                     :: rmatrix_dim_id
    !
    ! attributes
    CHARACTER (len = *),            PARAMETER ::           att_units  = "units"
    CHARACTER (len = *),            PARAMETER ::   charact_att_units  = "ch_units"
    CHARACTER (len = *),            PARAMETER ::   iscalar_att_units  = "is_units"
    CHARACTER (len = *),            PARAMETER ::   rscalar_att_units  = "rs_units"
    CHARACTER (len = *),            PARAMETER ::   ivector_att_units  = "iv_units"
    CHARACTER (len = *),            PARAMETER ::   rvector_att_units  = "rv_units"
    CHARACTER (len = *),            PARAMETER ::   rmatrix_att_units  = "rm_units"
    !
    INTEGER                                   :: dim_1, dim_2
    !
    

    !set dimensions

    dim_1 = dim_v
    dim_2 = dim_v

    
    CALL check( nf90_create(trim(datafile), nf90_clobber, ncid) )      ! Create the file. 

    WRITE(*,*) '# read_ncFile:: file created'
    !call check(nf90_create(path = trim(datafile),cmode = nf90_clobber,ncid = ncFileID))

    !define dimensions.


    CALL check( nf90_def_dim(ncid,   ivector_dim_name,   dim_v,   ivector_dim_id  ) )    !iv
    CALL check( nf90_def_dim(ncid,   rvector_dim_name,   dim_v,   rvector_dim_id  ) )    !rv
    CALL check( nf90_def_dim(ncid,   charact_dim_name,   dim_s,   charact_dim_id  ) )    !char
    CALL check( nf90_def_dim(ncid, rmatrix_dim_1_name,   dim_1, rmatrix_dim_id(1) ) )    ! row
    CALL check( nf90_def_dim(ncid, rmatrix_dim_2_name,   dim_2, rmatrix_dim_id(2) ) )    ! column

    
    WRITE(*,*) '# read_ncFile:: dimensions given'

    !define variables
    CALL check( nf90_def_var(ncid, iscalar_name, nf90_int,  iscalar_id)              )
    CALL check( nf90_def_var(ncid, rscalar_name, nf90_real, rscalar_id)              )
    CALL check( nf90_def_var(ncid, charact_name, nf90_char, charact_dim_id, charact_id) )
    CALL check( nf90_def_var(ncid, ivector_name, nf90_int,  ivector_dim_id, ivector_id) )
    CALL check( nf90_def_var(ncid, rvector_name, nf90_real, rvector_dim_id, rvector_id) )
    CALL check( nf90_def_var(ncid, rmatrix_name, nf90_real, rmatrix_dim_id, rmatrix_id) )
    
    WRITE(*,*) '# read_ncFile:: variables defined'
    !attributes
    CALL check( nf90_put_att(ncid, charact_id, att_units, charact_att_units) )
    CALL check( nf90_put_att(ncid, iscalar_id, att_units, iscalar_att_units) )
    CALL check( nf90_put_att(ncid, rscalar_id, att_units, rscalar_att_units) )
    CALL check( nf90_put_att(ncid, ivector_id, att_units, ivector_att_units) )
    CALL check( nf90_put_att(ncid, rvector_id, att_units, rvector_att_units) )
    CALL check( nf90_put_att(ncid, rmatrix_id, att_units, rmatrix_att_units) )
    
    WRITE(*,*) '# read_ncFile:: attributes given '

    ! Global attributes
    CALL check(nf90_put_att(ncid, nf90_global, "nc_File", &
          &   "created by l.aa. nikolopoulos OCT2008QUB"))
    CALL check(nf90_put_att(ncid, nf90_global, "dipole", &
         & "bs_fxd_2e"))
    
    WRITE(*,*) '# read_ncFile:: global identification set '

    CALL check( nf90_enddef(ncid) )         ! end of defining mode
    
    
    ! now write in ncFile
    
    CALL check( nf90_put_var(ncid, charact_id, charact) )     ! integer scalar
    CALL check( nf90_put_var(ncid, iscalar_id, iscalar) )     ! integer scalar
    CALL check( nf90_put_var(ncid, rscalar_id, rscalar) )     ! real    scalar
    CALL check( nf90_put_var(ncid, ivector_id, ivector) )     ! integer vector
    CALL check( nf90_put_var(ncid, rvector_id, rvector) )     ! real    vector
    CALL check( nf90_put_var(ncid, rmatrix_id, rmatrix) )     ! real    matrix 

    ! close the ncFile

    CALL check( nf90_close(ncid) )                            
    
    PRINT *,"# write_netCDF done!"    
  END SUBROUTINE write_ncFile
  !
  !
  !prototype
  SUBROUTINE read_ncFile(datafile, charact, iscalar, rscalar, ivector, rvector, rmatrix)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !variable_name  in ncFile 
    CHARACTER (len = *),      PARAMETER :: charact_name = "charact"
    CHARACTER (len = *),      PARAMETER :: iscalar_name = "iscalar"
    CHARACTER (len = *),      PARAMETER :: rscalar_name = "rscalar"
    CHARACTER (len = *),      PARAMETER :: ivector_name = "ivector"
    CHARACTER (len = *),      PARAMETER :: rvector_name = "rvector"
    CHARACTER (len = *),      PARAMETER :: rmatrix_name = "rmatrix"
    !
    !dimension_name
    CHARACTER (len = *), PARAMETER ::   charact_dim_name = "charact_dim"
    CHARACTER (len = *), PARAMETER ::   ivector_dim_name = "ivector_dim"
    CHARACTER (len = *), PARAMETER ::   rvector_dim_name = "rvector_dim"
    CHARACTER (len = *), PARAMETER :: rmatrix_dim_1_name = "rmatrix_dim_1"
    CHARACTER (len = *), PARAMETER :: rmatrix_dim_2_name = "rmatrix_dim_2"
    !actual variables 
    CHARACTER (len = *),                  INTENT(in)  :: datafile 
    CHARACTER (len = *),                  INTENT(out) :: charact
    INTEGER,                              INTENT(out) :: iscalar
    REAL(dpk),                              INTENT(out) :: rscalar
    INTEGER, DIMENSION(:),   POINTER,     INTENT(out) :: ivector
    REAL(dpk), DIMENSION(:),   POINTER,     INTENT(out) :: rvector
    REAL(dpk), DIMENSION(:,:), POINTER,     INTENT(out) :: rmatrix
    !variable_dim
    INTEGER                                          :: ivector_dim_id,  rvector_dim_id
    INTEGER                                          :: charact_dim_id
    INTEGER, DIMENSION(2)                            :: rmatrix_dim_id
    !variable_id ( to be identified via the variable_name)
    INTEGER                                          :: ncid
    INTEGER                                          :: charact_id
    INTEGER                                          :: iscalar_id, rscalar_id
    INTEGER                                          :: ivector_id, rvector_id
    INTEGER                                          :: rmatrix_id
    !variable_att (to be identified via the variable id and the att_units) 
!    CHARACTER (len = *),                   PARAMETER ::  att_units  = "units"
!    CHARACTER (len = *)                              ::  charact_att_units
!    CHARACTER (len = *)                              ::  iscalar_att_units, rscalar_att_units
!    CHARACTER (len = *)                              ::  ivector_att_units, rvector_att_units
!    CHARACTER (len = *)                              ::  rmatrix_att_units
    !
    INTEGER                                           :: dim_v, dim_s, dim_1, dim_2
    INTEGER :: ndims_in, nvars_in, ngatts_in, unlimdimid_in
    !exe!

    !
 
    ! open ncFile

    CALL check( nf90_open(TRIM(datafile), nf90_nowrite, ncid) )

    ! inquire info
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )


    ! read dims
    CALL check( nf90_inq_dimid(ncid,   ivector_dim_name,   ivector_dim_id  ) )
    CALL check( nf90_inq_dimid(ncid,   rvector_dim_name,   rvector_dim_id  ) )
    CALL check( nf90_inq_dimid(ncid,   charact_dim_name,   charact_dim_id  ) )
    CALL check( nf90_inq_dimid(ncid, rmatrix_dim_1_name, rmatrix_dim_id(1) ) )
    CALL check( nf90_inq_dimid(ncid, rmatrix_dim_2_name, rmatrix_dim_id(2) ) )




    ! now read  ncFile
    ! inquire variable_id from variable_name
    CALL check( nf90_inq_varid(ncid, iscalar_name, iscalar_id )  )     ! integer scalar
    CALL check( nf90_inq_varid(ncid, rscalar_name, rscalar_id )  )     ! real scalat
    CALL check( nf90_inq_varid(ncid, charact_name, charact_id )  )     ! character
    CALL check( nf90_inq_varid(ncid, ivector_name, ivector_id )  )     ! character
    CALL check( nf90_inq_varid(ncid, rvector_name, rvector_id )  )     ! character
    CALL check( nf90_inq_varid(ncid, rmatrix_name, rmatrix_id )  )     ! character

    CALL check( nf90_inquire_dimension(ncid,    ivector_dim_id, len=dim_v) )
    CALL check( nf90_inquire_dimension(ncid,    rvector_dim_id, len=dim_v) )
    CALL check( nf90_inquire_dimension(ncid,    charact_dim_id, len=dim_s) )
    CALL check( nf90_inquire_dimension(ncid, rmatrix_dim_id(1), len=dim_1) )
    CALL check( nf90_inquire_dimension(ncid, rmatrix_dim_id(2), len=dim_2) )

!    CALL check( nf90_inquire_dimension(ncid, rmatrix_id, len=rmatrix_dim) )

    PRINT*, 'dim_v = ', dim_v
    PRINT*, 'dim_s = ', dim_s
    PRINT*, 'dim_1 = ', dim_1
    PRINT*, 'dim_2 = ', dim_2


    ALLOCATE( ivector(dim_v))
    ALLOCATE( rvector(dim_v))
    ALLOCATE( rmatrix(dim_1,dim_2))
    ! inquire a
    !    CALL check( nf90_get_att  (ncid, charact_id, att_units, charact_att_units))     ! character
    !
    !
    CALL check( nf90_get_var  (ncid, iscalar_id,   iscalar  )  )
    CALL check( nf90_get_var  (ncid, rscalar_id,   rscalar  )  )
    CALL check( nf90_get_var  (ncid, charact_id,   charact  )  )
    CALL check( nf90_get_var  (ncid, ivector_id,   ivector  )  )
    CALL check( nf90_get_var  (ncid, rvector_id,   rvector  )  )
    CALL check( nf90_get_var  (ncid, rmatrix_id,   rmatrix  )  )



    !    CALL check( nf90_put_var(ncid, iscalar_id, iscalar) )     ! integer scalar
    !    CALL check( nf90_put_var(ncid, rscalar_id, rscalar) )     ! real    scalar
    !    CALL check( nf90_put_var(ncid, ivector_id, ivector) )     ! integer vector
    !    CALL check( nf90_put_var(ncid, rvector_id, rvector) )     ! real    vector
    !    CALL check( nf90_put_var(ncid, rmatrix_id, rmatrix) )     ! real    matrix        
    !    CALL check( nf90_close(ncid) )                            ! close the file.
    
    PRINT *,"# read_ncFile done!"    

  END SUBROUTINE read_ncFile




  !check netCDF based calls
  SUBROUTINE check(status)
    !
    use netcdf 
    !
    implicit none
    !
    INTEGER, INTENT ( in) :: status
    !
    IF(status /= nf90_noerr) THEN 
       PRINT *, TRIM(nf90_strerror(status))
      STOP "Stopped"
   END IF
   !
  END SUBROUTINE check
  
END MODULE ncFile

