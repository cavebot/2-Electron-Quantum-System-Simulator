!
! contains: open_ncFile, close_ncFile
!
!           write_rs_ncFile(filename, r, r_name)  : write REAL    scalar
!           write_is_ncFile(filename, r, r_name)  : write INTEGER scalar
!           write_rv_ncFile(filename, r, r_name)  : write REAL    vector
!           write_iv_ncFile(filename, r, r_name)  : write INTEGER vector
!           write_rm_ncFile(filename, r, r_name)  : write REAL    matrix
!           write_im_ncFile(filename, r, r_name)  : write INTEGER scalar
!
!
!           read_rs_ncFile(filename, r, r_name)   : ! as above but 'read'
!           read_is_ncFile(filename, r, r_name)
!           read_rv_ncFile(filename, r, r_name)
!           read_iv_ncFile(filename, r, r_name)
!           read_rm_ncFile(filename, r, r_name)
!           read_im_ncFile(filename, r, r_name)
!                                       ! when the dimensions names are needed
!           read_rv_dim_ncFile(filename, v, v_name, vd_name_1)
!           read_iv_dim_ncFile(filename, v, v_name, vd_name_1)
!           read_rm_dim_ncFile(filename, mx, mx_name, mxd_name_1, mxd_name_2)
!           read_im_dim_ncFile(filename, mx, mx_name, mxd_name_1, mxd_name_2)



! prototype functions
!
! write_ncFile(datafile, dim_s, charact, iscalar, rscalar, dim_v, ivector, rvector, rmatrix)
! read_ncFile(datafile, charact, iscalar, rscalar, ivector, rvector, rmatrix)
! check(status)



MODULE ncFile_new
  !
  !
  USE PRECISION, ONLY: dpk
  !
  IMPLICIT NONE
  !
  PUBLIC  create_ncFile,     close_ncFile
  PUBLIC  write_is_ncFile, read_is_ncFile
  PUBLIC  write_rs_ncFile, read_rs_ncFile 
  PUBLIC  write_iv_ncFile, read_iv_ncFile
  PUBLIC  write_rv_ncFile, read_rv_ncFile
  PUBLIC  write_im_ncFile, read_im_ncFile
  PUBLIC  write_rm_ncFile, read_rm_ncFile 

  CHARACTER(len=2), PARAMETER :: ncFile_start = "nc"
  CHARACTER(len=3), PARAMETER :: ncFile_end   = ".nc"
  CHARACTER(len=6), PARAMETER :: ncFile_logo  = "netCDF"
  

CONTAINS
  !
  !
  ! creat a ncFile with the name 'filename'
  ! 
  ! 
  SUBROUTINE create_ncFile(filename)
    !
    USE netCDF
    !
    IMPLICIT NONE
    !    
    CHARACTER (len = 30), INTENT(in) :: filename
    !
    INTEGER                          :: ncid
    !

        
    CALL check( nf90_create(TRIM(filename), nf90_clobber, ncid) )              !Create the file.
    CALL check( nf90_put_att(ncid, nf90_global,   "file",     filename ) )     ! Global attributes
    CALL check( nf90_put_att(ncid, nf90_global, "ncFile", "LAANDEC2016") )    
    CALL check( nf90_close( ncid ) )                                           ! close the file.    
    !
    
    PRINT*,'# create_ncFile::', TRIM(filename)
    !
    RETURN
    
    !    
  END SUBROUTINE create_ncFile
  !
  !
  ! close a ncFile with the name 'filename'
  !
  !
  SUBROUTINE close_ncFile(filename)
    !
    USE netCDF
    !
    IMPLICIT NONE
    !    
    CHARACTER (len = 30), INTENT(in) :: filename
    !
    INTEGER                          :: ncid
    !exe

        
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )     ! open ncFile!define new real variable
    CALL check( nf90_close( ncid ) )                               ! close the file.
    !
    PRINT*,'# close_ncFile::', TRIM(filename)
    !
    return    
    !    
  END SUBROUTINE close_ncFile


  ! write scalar
  !
  !


  SUBROUTINE write_rs_ncFile(filename, s, s_name )
    !
    USE netcdf
    !
    IMPLICIT NONE
    !

    !variable_name
    !    INTEGER,              INTENT(in) :: ncid
    CHARACTER (len = 30),INTENT(in)  :: filename 
    REAL(dpk),            INTENT(in) :: s
    CHARACTER (len = *),  INTENT(in) :: s_name
    !
    INTEGER                          :: s_id            !variable_id
    INTEGER                          :: ncid
    !

    CALL check( nf90_open(TRIM(filename), nf90_write, ncid) )     ! open ncFile!define new real variable
    CALL check(   nf90_redef(ncid))
    CALL check( nf90_def_var(ncid, s_name, nf90_real, s_id  ) ) !
    CALL check(  nf90_enddef(ncid)                            )   ! end of defining mode
    CALL check( nf90_put_var(ncid, s_id,  s                 ) ) ! save 'r'
    CALL check( nf90_close  (ncid)                            )   ! close the file.


    
    WRITE(*,*) '# write_rs_ncFile::', trim(s_name), ' saved.'

    RETURN
    !    
  END SUBROUTINE write_rs_ncFile
  !
  !
  !
  !
  SUBROUTINE read_rs_ncFile(filename, s, s_name )
    !
    USE netcdf
    !
    IMPLICIT NONE
    !

    !variable_name
    CHARACTER (len = 30), INTENT(in)  :: filename        !ncFile location
    REAL(dpk),            INTENT(out) :: s
    CHARACTER (len = *),  INTENT(in)  :: s_name
    !
    INTEGER                           :: s_id            !variable_id
    INTEGER                           :: ncid
    !
    INTEGER                           :: ndims_in,  nvars_in
    INTEGER                           :: ngatts_in, unlimdimid_in
    !exe!

    !
 

    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  ! inquire info
    CALL check( nf90_inq_varid( ncid,  s_name, s_id )  )      ! read real 'rs_name' id 
    CALL check(   nf90_get_var( ncid,  s_id,   s    )  )      ! read 'rs'
    CALL check(     nf90_close(ncid) )                          ! close ncFile


    
    WRITE(*,*) '# read_rs_ncFile::', TRIM(s_name), ' read.'

    RETURN
    
  END SUBROUTINE read_rs_ncFile
  !
  !
  ! write integer scalar
  !
  !
  SUBROUTINE write_is_ncFile(filename, s, s_name )
    !
    USE netcdf
    !
    IMPLICIT NONE
    !

    !variable_name
    !    INTEGER,              INTENT(in) :: ncid
    CHARACTER (len = 30),INTENT(in)  :: filename 
    INTEGER,              INTENT(in) :: s
    CHARACTER (len = *),  INTENT(in) :: s_name
    !
    INTEGER                          :: s_id            !variable_id
    INTEGER                          :: ncid
    !
    !    CHARACTER (len = 100)            :: filename    
    !
    
    CALL check( nf90_open(TRIM(filename), nf90_write, ncid) )                     ! open ncFile
    CALL check(   nf90_redef(ncid))                              !define new real variable
    CALL check( nf90_def_var(ncid, s_name, nf90_int, s_id     ) ) !
    CALL check(  nf90_enddef(ncid)                              ) ! end of defining mode
    CALL check( nf90_put_var(ncid, s_id,  s                   ) ) ! save 'r'
    CALL check( nf90_close(ncid) )                                ! close the file.
    !

    WRITE(*,*) '# write_is_ncFile::', TRIM(s_name), ' saved.'
    
    RETURN
    !    
  END SUBROUTINE write_is_ncFile
  !
  !
  !  read INTEGER scalar variable
  !
  !
  SUBROUTINE read_is_ncFile(filename, s, s_name )
    !
    USE netcdf
    !
    IMPLICIT NONE
    !

    !variable_name
    CHARACTER (len = 30), INTENT(in)  :: filename        !ncFile location
    INTEGER,              INTENT(out) :: s
    CHARACTER (len = *),  INTENT(in)  :: s_name
    !
    INTEGER                           :: s_id            !variable_id
    !
    INTEGER                           :: ndims_in,  nvars_in
    INTEGER                           :: ngatts_in, unlimdimid_in
    INTEGER                           :: ncid
    !exe!

    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  ! inquire info
    CALL check( nf90_inq_varid( ncid,  s_name, s_id )  )                            ! read real 's_name' id 
    CALL check(   nf90_get_var( ncid,  s_id,   s    )  )                            ! read 'is'
    CALL check( nf90_close(ncid) )                                                  ! close
    !
    
    WRITE(*,*) '# read_is_ncFile::', TRIM(s_name), ' read.'


    RETURN
    
  END SUBROUTINE read_is_ncFile

  !
  !
  !
  
  SUBROUTINE write_rv_ncFile(filename, v, v_name)
    !
    USE netcdf
    !
    IMPLICIT NONE
    ! variable_name
    CHARACTER (len = 30),              INTENT(in)  :: filename        !ncFile location
    REAL(dpk), DIMENSION(:),  POINTER,  INTENT(in) :: v
    CHARACTER (len = 100),              INTENT(in) :: v_name
    !
    CHARACTER (len = *), PARAMETER                 :: v_dim_1_name="rv_dim_1"    !variable dimension 1 name
    INTEGER                                        :: v_id              !variable_id  
    INTEGER                                        :: v_dim_id          !variable_dim id
    !
    INTEGER                                        :: ncid
    INTEGER                                        :: dim_1    
    ! 

    !get 1d array dimension
    
    dim_1 = SIZE(v)    

    
    CALL check( nf90_open(TRIM(filename), nf90_write, ncid) )                   ! open ncFile
    CALL check( nf90_redef(ncid) )                                               !sof define new variables
    CALL check( nf90_def_dim( ncid, v_dim_1_name,     dim_1, v_dim_id       ) )  !'v_name' dim variable(s)   
    CALL check( nf90_def_var( ncid,       v_name, nf90_real, v_dim_id, v_id ) )  !'v_name' variable   
    CALL check( nf90_enddef (ncid) )                                             ! eof define        
    CALL check( nf90_put_var( ncid,  v_id, v))                                   ! put
    CALL check( nf90_close(ncid) )                                               ! close ncFile

    
    WRITE(*,*) '# write_rv_ncFile::', TRIM(v_name), ' saved.'
    !
    RETURN
    !    
  END SUBROUTINE write_rv_ncFile
  !
  !
  !
  !   read 1d array data
  !
  !
  !
  
  SUBROUTINE read_rv_ncFile(filename, v, v_name)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !

    CHARACTER (len = 30),              INTENT(in)     :: filename        !ncFile location
    REAL(dpk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: v
    CHARACTER (len = *),               INTENT(in)     :: v_name          ! variable name
    !
    CHARACTER (len = *), PARAMETER          :: v_dim_1_name="rv_dim_1"    !variable dimension 1 name
    INTEGER                                 :: v_id            !variable id
    INTEGER                                 :: v_dim_id        !variable_dim id
    !
    INTEGER                                 :: dim_1
    !
    INTEGER                                 :: ndims_in,  nvars_in
    INTEGER                                 :: ngatts_in, unlimdimid_in
    INTEGER                                 :: ncid

    !

    !    v_dim_1_name = "e_dim_1" ! mx_name//"_dim_1"

    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  !inquire info
    CALL check( nf90_inq_varid( ncid,       v_name, v_id     ) )       !access id of variable 'v_name'    
    CALL check( nf90_inq_dimid( ncid, v_dim_1_name, v_dim_id ) )       !access ids of 'v_name' dimensions
    CALL check( nf90_inquire_dimension(ncid, v_dim_id, len=dim_1) ) !get dim 1 of vector 'v_name' 
    ALLOCATE( v(dim_1) )                                               !allocate storage for 'v_name' 
    CALL check( nf90_get_var( ncid, v_id,  v )  )                      ! read 'v'
    CALL check( nf90_close(ncid) )                                     ! close file
    !
    PRINT *,"# read_rv_ncFile::", TRIM(v_name), " read."      
    !
    RETURN
    !
  END SUBROUTINE read_rv_ncFile
  !
  ! read real vector variable 'v_name' with 'vd_name_1' name of its dimension
  !
  !
    SUBROUTINE read_rv_dim_ncFile(filename, v, v_name, vd_name_1)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !
    CHARACTER (len = 30),                 INTENT(in)  :: filename        !ncFile location
    REAL(dpk), DIMENSION(:), ALLOCATABLE, INTENT(out) :: v
    CHARACTER (len = *),                  INTENT(in)  :: v_name          ! variable name
    CHARACTER (len = *),                  INTENT(in)  :: vd_name_1       ! variable name
    !
    CHARACTER (len = 30)                    :: v_dim_1_name
    !
    !    CHARACTER (len = *), PARAMETER                    :: v_dim_1_name="iv_dim_1"    !variable dimension 1 name
    INTEGER                                 :: v_id            !variable id
    INTEGER                                 :: v_dim_id        !variable_dim id
    !
    INTEGER                                 :: dim_1
    !
    INTEGER                                 :: ndims_in,  nvars_in
    INTEGER                                 :: ngatts_in, unlimdimid_in
    INTEGER                                 :: ncid             !ncFile id
    !

    v_dim_1_name = TRIM(vd_name_1)
    
    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  !inquire info
    !
    CALL check( nf90_inq_varid( ncid,       v_name, v_id     ) )       !access id of variable 'v_name'    
    CALL check( nf90_inq_dimid( ncid, v_dim_1_name, v_dim_id ) )       !access ids of 'v_name' dimensions
    CALL check( nf90_inquire_dimension(ncid, v_dim_id, len=dim_1) )    !get dim 1 of vector 'v_name' 
    ALLOCATE( v(dim_1) )                                               !allocate storage for 'v_name' 
    CALL check( nf90_get_var( ncid, v_id,  v )  )                      ! read 'v'
    CALL check( nf90_close(ncid) )                                     ! close file
    !
    PRINT *,"# read_iv_ncFile::", TRIM(v_name), " read."      
    !
    RETURN
    !
  END SUBROUTINE read_rv_dim_ncFile

  !
  !
  !  write INTEGER vector
  !
  !
  
  SUBROUTINE write_iv_ncFile(filename, v, v_name)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !

    !variable_name
    CHARACTER (len = 30),              INTENT(in)  :: filename        !ncFile location
    INTEGER, DIMENSION(:),     POINTER, INTENT(in) :: v
    CHARACTER (len = *),                INTENT(in) :: v_name
    !

    CHARACTER (len = *), PARAMETER          :: v_dim_1_name="iv_dim_1"    !variable dimension 1 name
    INTEGER                                        :: v_id              !variable_id  
    INTEGER                                        :: v_dim_id          !variable_dim id
    !  
    INTEGER                                        :: dim_1
    INTEGER                                        :: ncid
    !exe
    

    
    dim_1 = SIZE(v)        !get 'v' dimensions
    
    !
    CALL check( nf90_open(TRIM(filename), nf90_write, ncid) )                  ! open ncFile
    CALL check( nf90_redef(ncid) )                                               !sof define new variables
    CALL check( nf90_def_dim( ncid, v_dim_1_name,     dim_1, v_dim_id       ) )  !'v_name' dim variable(s)   
    CALL check( nf90_def_var(ncid,        v_name, nf90_int,  v_dim_id, v_id ) )  !'v_name' variable   
    CALL check( nf90_enddef(ncid) )                                              ! eof define        
    CALL check( nf90_put_var( ncid,  v_id, v))                                   ! put
    CALL check( nf90_close(ncid) )                                               ! close ncFile
    !
    
    WRITE(*,*) '# write_iv_ncFile::', TRIM(v_name), ' saved.'
    !
    RETURN
    
  END SUBROUTINE write_iv_ncFile
  !
  !
  !
  !   read 1d array data
  !
  !
  !
  
  SUBROUTINE read_iv_ncFile(filename, v, v_name)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !
    CHARACTER (len = 30),                 INTENT(in)  :: filename        !ncFile location
    INTEGER,   DIMENSION(:), ALLOCATABLE, INTENT(out) :: v
    CHARACTER (len = *),                  INTENT(in)  :: v_name          ! variable name
    !
    CHARACTER (len = *), PARAMETER                    :: v_dim_1_name="iv_dim_1"    !variable dimension 1 name
    INTEGER                                 :: v_id            !variable id
    INTEGER                                 :: v_dim_id        !variable_dim id
    !
    INTEGER                                 :: dim_1
    !
    INTEGER                                 :: ndims_in,  nvars_in
    INTEGER                                 :: ngatts_in, unlimdimid_in
    INTEGER                                 :: ncid             !ncFile id
    !

    
    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  !inquire info
    !
    CALL check( nf90_inq_varid( ncid,       v_name, v_id     ) )       !access id of variable 'v_name'    
    CALL check( nf90_inq_dimid( ncid, v_dim_1_name, v_dim_id ) )       !access ids of 'v_name' dimensions
    CALL check( nf90_inquire_dimension(ncid, v_dim_id, len=dim_1) )    !get dim 1 of vector 'v_name' 
    ALLOCATE( v(dim_1) )                                               !allocate storage for 'v_name' 
    CALL check( nf90_get_var( ncid, v_id,  v )  )                      ! read 'v'
    CALL check( nf90_close(ncid) )                                     ! close file
    !
    PRINT *,"# read_iv_ncFile::", TRIM(v_name), " read."      
    !
    RETURN
    !
  END SUBROUTINE read_iv_ncFile
  !
  !
  ! read INTEGER vector variable 'v_name' with 'vd_name_1' name of its dimension
  !
  SUBROUTINE read_iv_dim_ncFile(filename, v, v_name, vd_name_1)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !
    CHARACTER (len = 30),                 INTENT(in)  :: filename        !ncFile location
    INTEGER,   DIMENSION(:), ALLOCATABLE, INTENT(out) :: v
    CHARACTER (len = *),                  INTENT(in)  :: v_name          ! variable name
    CHARACTER (len = *),                  INTENT(in)  :: vd_name_1       ! variable name
    !
    CHARACTER (len = 30)                    :: v_dim_1_name
    !
    !    CHARACTER (len = *), PARAMETER                    :: v_dim_1_name="iv_dim_1"    !variable dimension 1 name
    INTEGER                                 :: v_id            !variable id
    INTEGER                                 :: v_dim_id        !variable_dim id
    !
    INTEGER                                 :: dim_1
    !
    INTEGER                                 :: ndims_in,  nvars_in
    INTEGER                                 :: ngatts_in, unlimdimid_in
    INTEGER                                 :: ncid             !ncFile id
    !

    v_dim_1_name = TRIM(vd_name_1)
    
    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  !inquire info
    !
    CALL check( nf90_inq_varid( ncid,       v_name, v_id     ) )       !access id of variable 'v_name'    
    CALL check( nf90_inq_dimid( ncid, v_dim_1_name, v_dim_id ) )       !access ids of 'v_name' dimensions
    CALL check( nf90_inquire_dimension(ncid, v_dim_id, len=dim_1) )    !get dim 1 of vector 'v_name' 
    ALLOCATE( v(dim_1) )                                               !allocate storage for 'v_name' 
    CALL check( nf90_get_var( ncid, v_id,  v )  )                      ! read 'v'
    CALL check( nf90_close(ncid) )                                     ! close file
    !
    PRINT *,"# read_iv_ncFile::", TRIM(v_name), " read."      
    !
    RETURN
    !
  END SUBROUTINE read_iv_dim_ncFile
  !

  
  !
  !
  !
  ! write REAL matrix data
  !
  !
  !
  !
  
  SUBROUTINE write_rm_ncFile(filename, mx, mx_name)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !

    !variable_name
    CHARACTER (len = 30),              INTENT(in)  :: filename        !ncFile location
    REAL(dpk), DIMENSION(:,:), POINTER, INTENT(in) :: mx
    CHARACTER (len = *),                INTENT(in) :: mx_name
    !
    CHARACTER (len = *), PARAMETER                 :: mx_dim_1_name="rmx_dim_1"    !variable dimension 1 name
    CHARACTER (len = *), PARAMETER                 :: mx_dim_2_name="rmx_dim_2"    !variable dimension 1 name
    INTEGER                                        :: mx_id              !variable_id  
    INTEGER, DIMENSION(2)                          :: mx_dim_id          !variable_dim id
    !  
    INTEGER                                        :: dim_1, dim_2
    INTEGER                                        :: ncid             !ncFile id
    !

    !get 'mx' dimensions
    
    dim_1 = SIZE(mx, dim = 1)
    dim_2 = SIZE(mx, dim = 2) 

    !
    CALL check( nf90_open(TRIM(filename), nf90_write, ncid) )                     ! open ncFile
    CALL check( nf90_redef(ncid) )                                               !sof define new variables
    CALL check( nf90_def_dim( ncid,  mx_dim_1_name,   dim_1,  mx_dim_id(1) ) )   !'mx_name' dim variable(s)   
    CALL check( nf90_def_dim( ncid,  mx_dim_2_name,   dim_2,  mx_dim_id(2) ) )   ! 
    CALL check( nf90_def_var(ncid,   mx_name, nf90_real,  mx_dim_id,  mx_id  ) ) !'mx_name' variable   
    CALL check( nf90_enddef(ncid) )                                              ! eof define        
    CALL check( nf90_put_var( ncid,  mx_id, mx))                                 ! put
    CALL check( nf90_close(ncid) )                                               ! close ncFile
    !    
    WRITE(*,*) '# write_rm_ncFile::', TRIM(mx_name), ' saved.'
    !
    RETURN
    
  END SUBROUTINE write_rm_ncFile
  !
  !
  !
  !
  SUBROUTINE read_rm_ncFile(filename, mx, mx_name)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !
    CHARACTER (len = 30),                   INTENT(in)  :: filename        !ncFile location
    REAL(dpk), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: mx
    CHARACTER (len = *),                    INTENT(in)  :: mx_name          ! variable name
    !
    CHARACTER (len = *), PARAMETER          :: mx_dim_1_name="rmx_dim_1"    !variable dimension 1 name
    CHARACTER (len = *), PARAMETER          :: mx_dim_2_name="rmx_dim_2"    !variable dimension 1 name
    INTEGER                                 :: mx_id            ! variable id
    INTEGER, DIMENSION(2)                   :: mx_dim_id        ! variable_dim id
    !
    INTEGER                                 :: dim_1, dim_2    
    !
    INTEGER                                 :: ndims_in,  nvars_in
    INTEGER                                 :: ngatts_in, unlimdimid_in
    INTEGER                                 :: ncid                       !ncFile id
    !exe
    
    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  !inquire info
    CALL check( nf90_inq_varid( ncid,       mx_name, mx_id        ) )    !access id of variable 'mx_name'    
    CALL check( nf90_inq_dimid( ncid, mx_dim_1_name, mx_dim_id(1) ) )    !access ids of 'mx_name' dimensions
    CALL check( nf90_inq_dimid( ncid, mx_dim_2_name, mx_dim_id(2) ) )    !        
    CALL check( nf90_inquire_dimension(ncid, mx_dim_id(1),  len=dim_1) ) !get dim 1 of matrix 'mx_name' 
    CALL check( nf90_inquire_dimension(ncid, mx_dim_id(2),  len=dim_2) ) !get dim 2 of matrix 'mx_name'    
    ALLOCATE( mx(dim_1,dim_2) )                                          !allocate storage 'mx_name' matrix
    CALL check( nf90_get_var( ncid, mx_id,  mx )  )                      !read 'mx'
    CALL check( nf90_close(ncid) )                                       !close file
    !
    PRINT *,"# read_rm_ncFile::", TRIM(mx_name), " read."      
    !
    RETURN
    !
  END SUBROUTINE read_rm_ncFile
  !
  !
  ! read REAL 2D array variable 'mx_name' with 'mxd_name_1, mxd_name_2' names of its dimensions
  !
  !
    SUBROUTINE read_rm_dim_ncFile(filename, mx, mx_name, mxd_name_1, mxd_name_2)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !
    CHARACTER (len = 30),                   INTENT(in)  :: filename        !ncFile location
    REAL(dpk), DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: mx
    CHARACTER (len = *),                    INTENT(in)  :: mx_name          ! variable name
    CHARACTER (len = *),                    INTENT(in)  :: mxd_name_1       ! variable dimensions 1 name
    CHARACTER (len = *),                    INTENT(in)  :: mxd_name_2       ! variable dimensions 2 name
    !
    CHARACTER (len = 30)                    :: mx_dim_1_name  !variable dimension 1 name
    CHARACTER (len = 30)                    :: mx_dim_2_name         !variable dimension 1 name
    INTEGER                                 :: mx_id            ! variable id
    INTEGER, DIMENSION(2)                   :: mx_dim_id        ! variable_dim id
    !
    INTEGER                                 :: dim_1, dim_2    
    !
    INTEGER                                 :: ndims_in,  nvars_in
    INTEGER                                 :: ngatts_in, unlimdimid_in
    INTEGER                                 :: ncid                       !ncFile id
    !exe

    mx_dim_1_name = TRIM(mxd_name_1)
    mx_dim_2_name = TRIM(mxd_name_2)
    
    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  !inquire info
    CALL check( nf90_inq_varid( ncid,       mx_name, mx_id        ) )    !access id of variable 'mx_name'    
    CALL check( nf90_inq_dimid( ncid, mx_dim_1_name, mx_dim_id(1) ) )    !access ids of 'mx_name' dimensions
    CALL check( nf90_inq_dimid( ncid, mx_dim_2_name, mx_dim_id(2) ) )    !        
    CALL check( nf90_inquire_dimension(ncid, mx_dim_id(1),  len=dim_1) ) !get dim 1 of matrix 'mx_name' 
    CALL check( nf90_inquire_dimension(ncid, mx_dim_id(2),  len=dim_2) ) !get dim 2 of matrix 'mx_name'    
    ALLOCATE( mx(dim_1,dim_2) )                                          !allocate storage 'mx_name' matrix
    CALL check( nf90_get_var( ncid, mx_id,  mx )  )                      !read 'mx'
    CALL check( nf90_close(ncid) )                                       !close file
    !
    PRINT *,"# read_rm_ncFile::", TRIM(mx_name), " read."      
    !
    RETURN
    !
  END SUBROUTINE read_rm_dim_ncFile


  !
  !
  !
  ! write INTEGER matrix data
  !
  !
  !
  !
  
  SUBROUTINE write_im_ncFile(filename, mx, mx_name)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !

    !variable_name
    CHARACTER (len = 30),              INTENT(in)  :: filename        !ncFile location
    INTEGER,   DIMENSION(:,:), POINTER, INTENT(in) :: mx
    CHARACTER (len = *),                INTENT(in) :: mx_name
    !
    CHARACTER (len = *), PARAMETER                 :: mx_dim_1_name="imx_dim_1"    !variable dimension 1 name
    CHARACTER (len = *), PARAMETER                 :: mx_dim_2_name="imx_dim_2"    !variable dimension 1 name    
    INTEGER                                        :: mx_id              !variable_id  
    INTEGER, DIMENSION(2)                          :: mx_dim_id          !variable_dim id
    !  
    INTEGER                                        :: dim_1, dim_2
    INTEGER                                        :: ncid             !ncFile id
    !

    !get 'mx' dimensions
    
    dim_1 = SIZE(mx, dim = 1)
    dim_2 = SIZE(mx, dim = 2) 
    
    !
    CALL check( nf90_open(TRIM(filename), nf90_write, ncid) )                     ! open ncFile
    CALL check( nf90_redef(ncid) )                                               !sof define new variables
    CALL check( nf90_def_dim( ncid,  mx_dim_1_name,   dim_1,  mx_dim_id(1) ) )   !'mx_name' dim variable(s)   
    CALL check( nf90_def_dim( ncid,  mx_dim_2_name,   dim_2,  mx_dim_id(2) ) )   ! 
    CALL check( nf90_def_var(ncid,   mx_name, nf90_INT, mx_dim_id,  mx_id  ) ) !'mx_name' variable   
    CALL check( nf90_enddef(ncid) )                                              ! eof define        
    CALL check( nf90_put_var( ncid,  mx_id, mx))                                 ! put
    CALL check( nf90_close(ncid) )                                               ! close ncFile
    !

    WRITE(*,*) '# write_im_ncFile::', TRIM(mx_name), ' saved.'

    !
    RETURN
    !
  END SUBROUTINE write_im_ncFile
  !
  !
  !  read INTEGER matrix from ncFile
  !
  !
  SUBROUTINE read_Im_ncFile(filename, mx, mx_name)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !
    CHARACTER (len = 30),                 INTENT(in)  :: filename        !ncFile location
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: mx
    CHARACTER (len = *),                  INTENT(in)  :: mx_name          ! variable name
    !
    CHARACTER (len = *), PARAMETER                 :: mx_dim_1_name="imx_dim_1"    !variable dimension 1 name
    CHARACTER (len = *), PARAMETER                 :: mx_dim_2_name="imx_dim_2"    !variable dimension 1 name        
    INTEGER                                 :: mx_id            ! variable id
    INTEGER, DIMENSION(2)                   :: mx_dim_id        ! variable_dim id
    !
    INTEGER                                 :: dim_1, dim_2    
    !
    INTEGER                                 :: ndims_in,  nvars_in
    INTEGER                                 :: ngatts_in, unlimdimid_in
    INTEGER                                 :: ncid                       !ncFile id
    !exe


    WRITE(*,*) '# read_en_ncFile:: file opened:'
    
    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  !inquire info
    CALL check( nf90_inq_varid( ncid,       mx_name, mx_id        ) )    !access id of variable 'mx_name'    
    CALL check( nf90_inq_dimid( ncid, mx_dim_1_name, mx_dim_id(1) ) )    !access ids of 'mx_name' dimensions
    CALL check( nf90_inq_dimid( ncid, mx_dim_2_name, mx_dim_id(2) ) )    !        
    CALL check( nf90_inquire_dimension(ncid, mx_dim_id(1),  len=dim_1) ) !get dim 1 of matrix 'mx_name' 
    CALL check( nf90_inquire_dimension(ncid, mx_dim_id(2),  len=dim_2) ) !get dim 2 of matrix 'mx_name'    
    ALLOCATE( mx(dim_1,dim_2) )                                          !allocate storage 'mx_name' matrix
    CALL check( nf90_get_var( ncid, mx_id,  mx )  )                      !read 'mx'
    CALL check( nf90_close(ncid) )                                       !close file
    !
    PRINT *,"# read_im_ncFile::", TRIM(mx_name), " read."      
    !
    RETURN
    !
  END SUBROUTINE read_Im_ncFile
  !
  !
  ! read INTEGER 2D array variable 'mx_name' with 'mxd_name_1, mxd_name_2' names of its dimensions
  !
  !
  SUBROUTINE read_im_dim_ncFile(filename, mx, mx_name, mxd_name_1, mxd_name_2)
    !
    USE netcdf
    !
    IMPLICIT NONE
    !
    CHARACTER (len = 30),                   INTENT(in)  :: filename        !ncFile location
    INTEGER,   DIMENSION(:,:), ALLOCATABLE, INTENT(out) :: mx
    CHARACTER (len = *),                    INTENT(in)  :: mx_name          ! variable name
    CHARACTER (len = *),                    INTENT(in)  :: mxd_name_1       ! variable dimensions 1 name
    CHARACTER (len = *),                    INTENT(in)  :: mxd_name_2       ! variable dimensions 2 name
    !
    CHARACTER (len = 30)                    :: mx_dim_1_name  !variable dimension 1 name
    CHARACTER (len = 30)                    :: mx_dim_2_name         !variable dimension 1 name
    INTEGER                                 :: mx_id            ! variable id
    INTEGER, DIMENSION(2)                   :: mx_dim_id        ! variable_dim id
    !
    INTEGER                                 :: dim_1, dim_2    
    !
    INTEGER                                 :: ndims_in,  nvars_in
    INTEGER                                 :: ngatts_in, unlimdimid_in
    INTEGER                                 :: ncid                       !ncFile id
    !exe

    mx_dim_1_name = TRIM(mxd_name_1)
    mx_dim_2_name = TRIM(mxd_name_2)
    
    !
    CALL check( nf90_open(TRIM(filename), nf90_nowrite, ncid) )                     ! open ncFile
    CALL check( nf90_inquire(ncid, ndims_in, nvars_in, ngatts_in, unlimdimid_in) )  !inquire info
    CALL check( nf90_inq_varid( ncid,       mx_name, mx_id        ) )    !access id of variable 'mx_name'    
    CALL check( nf90_inq_dimid( ncid, mx_dim_1_name, mx_dim_id(1) ) )    !access ids of 'mx_name' dimensions
    CALL check( nf90_inq_dimid( ncid, mx_dim_2_name, mx_dim_id(2) ) )    !        
    CALL check( nf90_inquire_dimension(ncid, mx_dim_id(1),  len=dim_1) ) !get dim 1 of matrix 'mx_name' 
    CALL check( nf90_inquire_dimension(ncid, mx_dim_id(2),  len=dim_2) ) !get dim 2 of matrix 'mx_name'    
    ALLOCATE( mx(dim_1,dim_2) )                                          !allocate storage 'mx_name' matrix
    CALL check( nf90_get_var( ncid, mx_id,  mx )  )                      !read 'mx'
    CALL check( nf90_close(ncid) )                                       !close file
    !
    PRINT *,"# read_rm_ncFile::", TRIM(mx_name), " read."      
    !
    RETURN
    !
  END SUBROUTINE read_im_dim_ncFile


  
  !
  !
  !  NOT USED
  !
  ! end of save and read of 0-2d arrays (0d) scalars, 1d (vectors), 2d matrices
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
    REAL(dpk),                            INTENT(out) :: rscalar
    INTEGER, DIMENSION(:),     POINTER,   INTENT(out) :: ivector
    REAL(dpk), DIMENSION(:),   POINTER,   INTENT(out) :: rvector
    REAL(dpk), DIMENSION(:,:), POINTER,   INTENT(out) :: rmatrix
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
  
END MODULE ncFile_new

