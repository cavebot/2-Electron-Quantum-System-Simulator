

PROGRAM create_atomic_ncFile
  !
  use netcdf
  USE ncFile
  USE PRECISION
  USE units
  USE PARAMETER
  USE deriv
  USE io

  !
  IMPLICIT NONE
  !
  INTEGER                                  :: l_max
  !
  REAL(dpk), DIMENSION(:,:),       POINTER :: d_if       ! ni,nf    d(l,l+1)
  REAL(dpk), DIMENSION(:),         POINTER :: e_i        ! e_i(n_i) 
  REAL(dpk), DIMENSION(:),         POINTER :: e_f        ! e_f(n_f)
  REAL(dpk)                                :: n_i, n_f
  INTEGER                                  :: l_i, l_f   !
  !
  INTEGER                                  :: i,j, k,l
  !
  CHARACTER(LEN=*)                         :: filename
  CHARACTER( LEN = 30  )                   :: datafile
  CHARACTER( LEN = 15  )                   :: sl1,sl2
  CHARACTER( LEN = 6  )                    :: gauge
  !EXE!
    !
    !........
    

  CALL getarg(1, argv)               ! 
  READ(argv,*) l_max                 ! max partial wave
  CALL getarg(2, argv)               ! 
  READ(argv,*) smode                 ! gauge
  

  filename = 'he'

  DO l = 0, lmax-1


     CALL dmxfile(nbin, "dat", "dmx2ebb","bin", gauge, l, l+1) 

     !CALL DMX2EFILE(17, i-1, i )

         READ(nbin) mode
         READ(nbin) l_i, l_f, n_i, n_f

         ALLOCATE(  e_i(n_i)     )
         ALLOCATE(  e_f(n_f)     )
         ALLOCATE( d_if(n_f,n_i) )

         READ(nbin) ( e_i(ki), ki = 1, size(e_i)  )
         READ(nbin) ( e_f(kf), kf = 1, size(e_f)  )

         read_dmx2e:DO  ki = 1, n_i
            READ(nbin) (dzr(kf,ki), kf = 1, n_f)
         END DO read_dmx2e


         CLOSE(nbin)

         WRITE(*,'(a7,i3,a1,i3,a6)') 'dipole(', i, ',', i + 1,') read'

         ! write atomic ncFile

         WRITE(sl1,'(I15)') l
         WRITE(sl2,'(I15)') l+1

     datafile  = "dat/"//FILENAME & 
                       //TRIM(ADJUSTL(sl1))//TRIM(ADJUSTL(sl2)) &
                       //"."//TRIM(ADJUSTL(gauge))              &  
                       //".dat"



     CALL write_atomic_ncFile(datafile, e_i, e_f, d_if)

     DEALLOCATE(e_i,e_f)
     DEALLOCATE(d_if)

  ENDDO



      !      CALL write_ncFile(datafile, nchar, s_in, i_in, r_in, ndim, iv_in, rv_in, rm_in)
      !      CALL  read_ncFile(datafile, s_out, i_out, r_out, iv_out, rv_out, rm_out)
  

  

END PROGRAM create_atomic_ncFile




