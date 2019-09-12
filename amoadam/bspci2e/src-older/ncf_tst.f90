

PROGRAM ncf_tst
  !
  use netcdf
  USE ncFile

  !
  IMPLICIT NONE
  !
  CHARACTER( LEN =  30 )                   :: gauge
  !
  REAL(dpk), DIMENSION(:,:),       POINTER :: d_if       ! ni,nf    d(l,l+1)
  REAL(dpk), DIMENSION(:),         POINTER :: e_i        ! e_i(n_i) 
  REAL(dpk), DIMENSION(:),         POINTER :: e_f        ! e_f(n_f)  !
  INTEGER                                  ::l
  !
  CHARACTER( LEN =  30 )                  :: filename
  CHARACTER( LEN =  30 )                  :: atomic_name


     atomic_name  = "dat_1"
     atomic_name  = "dat_2"
     !

     CALL write_atomic_ncFile(filename, atomic_name, gauge, l, e_i, e_f, d_if)

      !      CALL write_ncFile(datafile, nchar, s_in, i_in, r_in, ndim, iv_in, rv_in, rm_in)
      !      CALL  read_ncFile(datafile, s_out, i_out, r_out, iv_out, rv_out, rm_out)
 
  

   END PROGRAM ncf_tst




