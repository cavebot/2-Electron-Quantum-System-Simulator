!
! LAAN: accepts as an input a fortran dmx file and turns into netCDF ready 
!       to use it in tdse calculations.
!
!

PROGRAM ncf
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
  INTEGER                                  :: n_max
  INTEGER                                  :: l_min, l_max
  CHARACTER( LEN =  30 )                   :: gauge
  !
  REAL(dpk), DIMENSION(:,:),       POINTER :: d_if       ! ni,nf    d(l,l+1)
  REAL(dpk), DIMENSION(:),         POINTER :: e_i        ! e_i(n_i) 
  REAL(dpk), DIMENSION(:),         POINTER :: e_f        ! e_f(n_f)
  INTEGER                                  :: n_i, n_f
  INTEGER                                  :: l_i, l_f   !
  !
  INTEGER                                  :: ki,kf,l
  !
  CHARACTER( LEN = 100 )                  :: filename
  CHARACTER( LEN =  30 )                  :: atomic_name
  CHARACTER( LEN =  15 )                  :: sl1,sl2
  INTEGER                                 :: mode
  !EXE!
    !
    !........
    

  CALL getarg(1, argv)               ! 
  READ(argv,*) l_min                 ! max partial wave
  CALL getarg(2, argv)               ! 
  READ(argv,*) l_max                 ! max partial wave
  CALL getarg(3, argv)               ! 
  READ(argv,*) gauge                 ! gauge
  CALL getarg(4, argv)               ! 
  READ(argv,*) n_max                 ! nof states taken
!  CALL getarg(5, argv)               ! 
!  READ(argv,*) atomic_name           ! gauge

  atomic_name = 'he'

  DO l = l_min, l_max-1

     CALL dmxfile(nbin, "dat", "he","bin", gauge, l, l+1) 

         READ(nbin) mode
         READ(nbin) l_i, l_f, n_i, n_f

         IF(n_max.NE.0) THEN 
            n_i = n_max
            n_f = n_max
            WRITE(*,'(a1,a60,i5)') '#','sizeof symmetries changed to: ', n_max
         ENDIF
         
         ALLOCATE(  e_i(n_i)     )
         ALLOCATE(  e_f(n_f)     )
         ALLOCATE( d_if(n_f,n_i) )



         READ(nbin) ( e_i(ki), ki = 1, size(e_i)  )
         READ(nbin) ( e_f(kf), kf = 1, size(e_f)  )
         read_dmx2e:DO  ki = 1, size(e_i)
            READ(nbin) (d_if(kf,ki), kf = 1, SIZE(e_f))
         END DO read_dmx2e

         CLOSE(nbin)

         
         WRITE(*,'(a1,a60,i6)') '#','sizeof  e_i = ', SIZE(e_i)
         WRITE(*,'(a1,a60,i6)') '#','sizeof  e_f = ', SIZE(e_f)

         e_i = e_i*0.5_dpk
         e_f = e_f*0.5_dpk
         !         IF(l == 0) e_i(1) = -2.903

         ! write atomic ncFile

         WRITE(sl1,'(I15)') l
         WRITE(sl2,'(I15)') l+1

     filename  = "dat/"//trim(adjustl(atomic_name)) & 
                       //TRIM(ADJUSTL(sl1))//TRIM(ADJUSTL(sl2)) &
                       //TRIM(ADJUSTL(gauge))                   &  
                       //".nc"
     
     !

     CALL write_atomic_ncFile(filename, atomic_name, gauge, l, e_i, e_f, d_if)

     !
     DEALLOCATE(e_i,e_f)
     DEALLOCATE(d_if)

  ENDDO



      !      CALL write_ncFile(datafile, nchar, s_in, i_in, r_in, ndim, iv_in, rv_in, rm_in)
      !      CALL  read_ncFile(datafile, s_out, i_out, r_out, iv_out, rv_out, rm_out)
  

  

END PROGRAM ncf




