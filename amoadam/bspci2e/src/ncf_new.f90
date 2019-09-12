!
!  test new netCDF module for I/O
!
!

PROGRAM ncf_new
  !
  use netcdf
  USE PRECISION
  USE ncFile_new

  !
  IMPLICIT NONE
  !
  !
  INTEGER                                :: i
  INTEGER                                :: ncid
  CHARACTER( LEN =  30 )                 :: base_name
  CHARACTER( LEN = 100 )                 :: filename
  !real in
  REAL(dpk)                              :: rsin
  REAL(dpk), DIMENSION(:),       POINTER :: rvin
  REAL(dpk), DIMENSION(:,:),     POINTER :: rmin
  !int in
  INTEGER                                :: isin
  INTEGER, DIMENSION(:),         POINTER :: ivin
  INTEGER, DIMENSION(:,:),       POINTER :: imin
  !real out
  REAL(dpk)                              :: rsout
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: rvout
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: rmout
  !int out
  INTEGER                                :: isout
  INTEGER,   ALLOCATABLE, DIMENSION(:)   :: ivout
  INTEGER,   ALLOCATABLE, DIMENSION(:,:) :: imout
  !
  CHARACTER( LEN = 100)                  :: rs_name
  CHARACTER( LEN = 100)                  :: rv_name
  CHARACTER( LEN = 100)                  :: rm_name
  CHARACTER( LEN = 100)                  :: is_name
  CHARACTER( LEN = 100)                  :: iv_name
  CHARACTER( LEN = 100)                  :: im_name
  !
  INTEGER                                :: dim_1, dim_2
  
  !
  
  base_name = "ncFile-tst.nc"
  filename  = "tdat/"//TRIM(ADJUSTL(base_name))


  dim_1 = 4
  dim_2 = 4
  
  !
  rs_name   = "real 0d array"
  rv_name   = "real 1d array"
  rm_name   = "real 2d array"
  is_name   = "int  0d array"
  iv_name   = "int  1d array"
  im_name   = "int  2d array"

  !
  ALLOCATE(rvin(dim_1))
  ALLOCATE(ivin(dim_1))
  ALLOCATE(rmin(dim_1,dim_2))
  ALLOCATE(imin(dim_1,dim_2))

  
  rsin = 0.0_dpk
  rvin = 1.0_dpk
  rmin = 2.0_dpk
  
  isin = 0
  ivin = 1
  imin = 2.0_dpk
  

  PRINT*, "   write : "
  PRINT*, "    rsin = ", rsin
  PRINT*, "    rvin = "
  WRITE(*,'(10E15.4))') rvin
  PRINT*, "    rmin = "
  DO i = 1, SIZE(rmin,dim=1)
     WRITE(*,'(10E15.4))') rmin(i,:)
  ENDDO
  PRINT*, "    isin = ", isin
  PRINT*, "    ivin = "
  WRITE(*,'(10I5))') ivin
  PRINT*, "    imin = "
  DO i = 1, SIZE(imin,dim=1)
     WRITE(*,'(10I10))') imin(i,:)
  ENDDO


  
  CALL create_ncFile(filename)          ! open file
  
  !write in the ncFile named 'filename'
  CALL write_rs_ncFile(filename, rsin,  rs_name)  !save real scalar as 'rs_name'
  CALL write_rv_ncFile(filename, rvin,  rv_name)  !save real vector as 'rv_name'
  CALL write_rm_ncFile(filename, rmin,  rm_name)  !save real matrix as 'rm_name'
  CALL write_is_ncFile(filename, isin,  is_name)  !save int scalar  as 'is_name'
  CALL write_iv_ncFile(filename, ivin,  iv_name)  !save int vector  as 'iv_name'
  CALL write_im_ncFile(filename, imin,  im_name)  !save int matrix  as 'im_name'
  
  !read from the ncFile named 'filename'
  CALL read_rs_ncFile(filename, rsout, rs_name)  !read real scalar saved as 'rs_name'
  CALL read_rv_ncFile(filename, rvout, rv_name)  !read real vector saved as 'rv_name'
  CALL read_rm_ncFile(filename, rmout, rm_name)  !read real matrix saved as 'rm_name'
  CALL read_is_ncFile(filename, isout, is_name)  !read int scalar  saved as 'is_name'
  CALL read_iv_ncFile(filename, ivout, iv_name)  !read int vector  saved as 'iv_name'
  CALL read_im_ncFile(filename, imout, im_name)  !read int matrix  saved as 'im_name'

  !tst now if the values are read ok.
  
  PRINT*, "  rsout = ", rsout
  PRINT*, "  rvout = "
  WRITE(*,'(10E15.4))') rvout
  PRINT*, "  rmout = "
    DO i = 1, SIZE(rmout,dim=1)
     WRITE(*,'(10E15.4))') rmout(i,:)
  ENDDO  
  PRINT*, "  isout = ", isout
  PRINT*, "  ivout = "
  WRITE(*,'(10I10))') ivout
  PRINT*, "  imout = "
  DO i = 1, SIZE(imout,dim=1)
     WRITE(*,'(10I10))') imout(i,:)
  ENDDO
  
  
  RETURN
     
END PROGRAM ncf_new




