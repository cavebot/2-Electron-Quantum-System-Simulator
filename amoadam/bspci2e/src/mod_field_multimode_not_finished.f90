MODULE pulse
  !
  USE PRECISION
  !
  IMPLICIT NONE
  !
  PUBLIC
  !pulse!

  REAL(dpk)                           :: e0
  REAL(dpk)                           :: w0
  INTEGER                             :: cycles
  REAL(dpk)                           :: phi
  CHARACTER(len=10)                   :: pulseType
  REAL(dpk)                           :: tau
  REAL(dpk)                           :: dw_fwhm, dw_step
  REAL(dpk)                           :: w_min, w_max
  INTEGER                             :: n_modes
  REAL(dpk), ALLOCATABLE, DIMENSION(:):: e_i
  REAL(dpk), ALLOCATABLE, DIMENSION(:):: w_i
  REAL(dpk), ALLOCATABLE, DIMENSION(:):: p_i
  REAL(dpk), ALLOCATABLE, DIMENSION(:):: dw_i

  !!
  !!  REAL(DPK) FUNCTION E_T( T, PULSETYPE, E0, TAU, OMEGA, PHI)
  !!

CONTAINS


  SUBROUTINE READ_FIELD
    !
    USE units, ONLY:i0_td, enau, m_pi, t_fs
    USE io,    only:ninp
    !    
    IMPLICIT NONE
    !
    CHARACTER(LEN=25)      :: field_file = "tinp/pulse.inp"
    REAL(dpk)              :: i0                   ! e0 = sqrt(i0): electric field
    REAL(dpk)              :: it_av
    INTEGER                :: i,j
    !
    WRITE(*,*) '# opening field file :', field_file
    !


    ! CALL getarg(1, argv)              ! peak intensity
    ! READ(argv,*)   i0                 ! in SI 
    ! CALL getarg(2, argv)              ! photon frequency
    ! READ(argv,*)   omega              ! in SI
    ! CALL getarg(3, argv)              ! nof cycles
    ! READ(argv,*)   cycles             !  

   
    OPEN( ninp, file=TRIM(field_file),status='old')
    READ( ninp, * ) n_modes
    READ( ninp, * ) pulsetype ! pulse shape ('agauss', 'acos2' 'acos2flat', 'asin2', 'esin2','composite' )
    READ( ninp, * ) i0
    READ( ninp, * ) w0
    

    e0  = SQRT(i0/i0_td)                   ! E_o : W/cm^2 -->     (a.u.)
    w0  = w0/enau                          ! Wph :     eV -->     (a.u.)

    ALLOCATE( e_i(ABS(n_modes)) )
    ALLOCATE( w_i(ABS(n_modes)) )
    ALLOCATE( p_i(ABS(n_modes)) )
    ALLOCATE(dw_i(ABS(n_modes)) )


    IF(n_modes.EQ.1) THEN

       READ( ninp, * ) phi
       READ( ninp, * ) cycles 
       CLOSE(ninp)
       tau = cycles * (2.0*m_pi/w0)           ! T_L : pulse duration (a.u.)
       
       e_i(1) = e0
       w_i(1) = w0
       p_i(1) = phi
       
       !    CLOSE(ninp)
       
       ! peak intensity in W/cm^2 
       ! central field frequency in eV
       
       
    ELSE    IF( n_modes > 1 ) THEN 
       
       read_field_components:DO i = 1, n_modes
          READ(ninp,*)  e_i(i), w_i(i), p_i(i) 
       ENDDO read_field_components
       CLOSE(ninp)

    ELSE IF(n_modes < 0) THEN

       n_modes = - n_modes
       
       READ(ninp,*) w_min, w_max, dw_fwhm
       CLOSE(ninp)

       dw_step = ABS(w_max - w_min)/(n_modes-1)

       make_w_components:DO i = 1, n_modes
          w_i(i) = w_min + (i-1) * dw_step 
       ENDDO make_w_components

       w0 = w_i((n_modes-1)/2) ! should be    = (w_max + w_min)/2.0_dpk  

       make_amplitude_components:DO i = 1, n_modes
          e_i(i) = EXP( - ( (w_i(i) - w0)/dw_fwhm)**2)
       ENDDO make_amplitude_components
       
       make_phase_components:DO i = 1, n_modes
          p_i(i) = 0.0_dpk
       ENDDO make_phase_components

             
    ENDIF

!    CLOSE(ninp)


       !       p_i = 0.0_dpk
       !       w_i   = w_i * 10.0_dpk * t_fs              ! in a.u. 


    ! time average of the field 

    it_av = SUM(e_i**2)

    dw_i(1) = 0.0_dpk 
    calculate_field_norm:DO i = 2, n_modes - 1
       dw_i(i) = ABS(w_i(i+1) - w_i(i)) 
       
       !       it_av = it_av + dw_i * e_i(i)**2       
    ENDDO calculate_field_norm
    dw_i(n_modes) = ABS(w_i(n_modes) - w_i(n_modes-1) ) 


    it_av = DOT_PRODUCT(dw_i,e_i**2)        !intended to be:  it_av = it_av + dw_i * e_i(i)**2       

       !

       !
!       renormalize_field:DO i = 1, n_modes-1
!          e_i(i) = (e_i(i)*sqrt(dw_i)/ SQRT(it_av) ) * e0
!       ENDDO renormalize_field
!       e_i(n_modes) = (e_i(n_modes)*sqrt(dw_i)/ SQRT(it_av) ) * e0

     ! renormalize
       e_i = (e_i*SQRT(dw_i)/ SQRT(it_av) ) * e0


       !       e_i = (e_i/ sqrt(it_av) ) * e0


       WRITE(*,*) '& it_av = ', it_av    

!       DO i = 1, n_modes
!          WRITE(*,*) w_i(i), e_i(i)
!       ENDDO
       
    !IF(pulsetype=='composite_file') THEN       
       !    ENDIF




    !
    WRITE(*,*) '#  field file read done.'
    !
  END SUBROUTINE READ_FIELD
!
!
!

  SUBROUTINE SAVE_FIELD(ti, nout_T)
    !
    use io
    use units
    !    
    implicit none
    !
    REAL(dpk)              :: ti                     ! initial time
    INTEGER                :: nout_T                 ! nof output
    !
    REAL(dpk)              :: t, tp, td
    REAL(dpk)              :: i0 
    INTEGER                :: iout
    CHARACTER(LEN=25)      :: field_file  = "tout/pulse.out"


    !
    i0 = e0**2 * i0_td ! W/cm^2
    Tp = tau/cycles

    OPEN( nout, file=trim(field_file))
    WRITE( nout,  '(a30,2E25.10)') '&                 E0 = ', e0                ! a.u.
    WRITE( nout,  '(a30,2E25.10)') '&                  w = ', w0,w0*enau  ! a.u., SI
    WRITE( nout,  '(a30,i5)'     ) '&             cycles = ', cycles           
    WRITE( nout,  '(a30,E25.5)'  ) '&                phi = ', phi
    WRITE( nout,  '(a30,a25)'    ) '& pulsetype    shape = ', pulsetype 
    WRITE( nout,  '(a30,2E25.10)') '& intensity       i0 = ', e0**2, i0         ! a.u., SI
    WRITE( nout,  '(a30,2E25.10)') '& pulse duration tau = ', tau
    WRITE( nout,  '(a30,2E25.10)') '& initial time    ti = ', ti
    WRITE( nout,  '(a30, E25.10)') '& pulse period    tp = ', tp                !a.u.
    WRITE( nout,  '(a1,3a25)'    ) '&','t(au)', 'A(t)', 'E(t)'


    td = 20.

    DO iout = -2*cycles*nout_T, 2*cycles*nout_T
       t = ti + iout * tp/nout_T
       WRITE(nout,'(3E25.14)') t, A_t(t), A_t(t+td)
!       WRITE(nout,'(3E25.14)') t, A_t(t), E_t(t)
    ENDDO


    CLOSE(nout)
    WRITE(*,*) "# pulse::save_field:          i_end =  ", iout    
    WRITE(*,*) "# pulse::save_field:          t_end =  ", t    
    WRITE(*,*) "# pulse::save_field:   pulse saved at ", field_file
  END SUBROUTINE SAVE_FIELD




!
!--------------------------------------------------------------------
!
!     ..
!
!  Purpose
!  =======
!
!  POTENTIAL_AMPLI returns the value at time T of the vector potential A( T )
!  coresponding to the physical field E( T ) and defined as :
!
!     E( T ) = - D / DT A( T )  
!  
!  with 
!
!     A( T ) = E0 / OMEGA * F( T, TAU ) * COS( OMEGA * T + PHI )
!
!  where F( T, TAU ) is the envelope of the pulse, OMEGA the photon energy
!  and PHI the initial phase.
!

!    REAL(DPK) FUNCTION A_T( T, PULSETYPE, E0, TAU, OMEGA, PHI )
!      REAL(DPK)             E0, TAU, PHI, OMEGA 
!      CHARACTER(DPK)        PULSETYPE


  REAL(DPK) FUNCTION A_T(T)

      !
      USE PRECISION
      USE UNITS
      !
      IMPLICIT NONE
      !
      REAL(DPK)             T 
      !
      REAL(DPK)     SHAPE, SHAPE_DER
      REAL(DPK)     OMEGA_P
      REAL(DPK)     C_1, C_2
      REAL(DPK)     T1, T2
      REAL(DPK)     a_i
      INTEGER       i
      !EX!
        
      !     .. Executable Statement ..

      SHAPE     = 0.0_dpk
      SHAPE_DER = 0.0_dpk
    

      
      composite_pulse: IF(pulsetype.EQ.'composite') THEN 

         zero_the_pulse:IF ( DABS( 2.D0 * T ).LT.10*TAU ) THEN

            a_t = 0.0_dpk
            build_the_pulse:DO i = 1, n_modes
               
               a_i = 1.0_dpk !e_i(i) / w_i(i)
               
               a_t = a_t +  a_i * COS( w_i(i)* t + p_i(i) )
            ENDDO build_the_pulse
            
         ENDIF zero_the_pulse

         RETURN

      ENDIF composite_pulse



      !
      !     Gaussian Pulse ( Tau = Half-width at 1/e)
      !

      gauss_pulse:IF ( PULSETYPE.EQ.'agauss' ) THEN
         SHAPE = DEXP( - 0.5D0 * ( T / TAU ) * ( T / TAU ) )
      ENDIF Gauss_pulse


      !
      !     Cosine-square pulse ( Tau = Total width )
      !

      cos2_pulse:IF ( PULSETYPE.EQ.'acos2' ) THEN
         IF ( DABS( 2.D0 * T ).LT.TAU ) THEN
            SHAPE = DCOS( M_PI * T /TAU ) ** 2
         ENDIF   
      ENDIF cos2_pulse


      !
      !     Sine-squared pulse ( Tau = Total width )
      !


      sin2_pulse:IF ( PULSETYPE.EQ.'asin2' ) THEN
         IF ( DABS( T ).LT.TAU ) THEN
            SHAPE = DSIN( M_PI * T /TAU ) ** 2
         ENDIF   
      ENDIF sin2_pulse
      


 !     stop

      !
      !     Flat pulse with CosSqr turn on and off :
      !
      !        - CosSqr turn on to -T1 with a time duration T2 
      !          ( half the duration of an entire CosSqr pulse ) 
      !        - Constant amplitude ( = 1 ) from -T1 to T1
      !        - CosSqr turn off from T1 with a time duration T2 
      !
      !     T1 and T2 are coded in TAU as the following :
      !        - T1 = integer part of TAU / 100 ( 2 decimals )
      !        - T2 = decimal part of TAU * 10000 ( Maximum 9999 u.a.)
      !
      !     Example : TAU = 12504.002084 => T1 = 125.04, T2 = 20.84
      !

      cos2flat_pulse:IF (PULSETYPE.EQ.'acos2flat') THEN

	 T1 = DINT( TAU ) / 1.D2
	 T2 = ( TAU - DINT( TAU ) ) * 1.D4
         T2 = 2.D0 * T2
	 SHAPE = 0.D0

         IF ( DABS( T ).LT.( T1 + T2/2.D0 ) ) THEN
            IF ( DABS( T ).LT.T1 ) THEN
               SHAPE = 1.D0
            ELSE
               SHAPE = DCOS( M_PI * ( DABS( T ) - T1 ) / T2 ) ** 2
            ENDIF
         ENDIF
      ENDIF cos2flat_pulse
      
      !



      IF (PULSETYPE.EQ.'esin2') THEN

         C_2 = cycles / DBLE( cycles**2 - 1) 
         C_1 = cycles * C_2

         OMEGA_P = 2 * M_PI/TAU

         SHAPE     = 0.0_dpk      
         SHAPE_DER = 1.0_dpk - C_1 
         IF ( DABS( T ).LT.TAU ) THEN
            SHAPE_DER = 1.0_DPK - C_1 * DCOS( OMEGA_P * T )
            SHAPE =           C_2 * DSIN( OMEGA_P * T )
         ENDIF

         !         WRITE(*,*) c_1, c_2, cycles
         !         WRITE(*,*) omega_p,tau
         !         WRITE(*,*) shape, shape_der

         A_T = - 0.5_dpk * E0 *  ( SHAPE * DCOS( w0 * T + PHI ) & 
              &         +      SHAPE_DER * DSIN( w0 * T + PHI )  )/w0


      ELSE

         A_T = E0 * SHAPE * DCOS( w0 * T + PHI )/w0
         
      ENDIF

      !
      RETURN


    END FUNCTION A_T

  !     ..
  !
  !  Purpose
  !  =======
  !
  !  E_T returns the value of an electromagnitc field E( T )
  !  defined as :
  !
  !     E( T ) = - D / DT A( T )
  !
  !  where A( T ) is defined as    
  !
  !  A( T ) = E0 / OMEGA * F( T, TAU ) * COS( OMEGA * T + PHI ) 
  !  
  !  and where F( T, TAU ) represents the shape of the pulse and is specified by
  !  PULSETYPE.
  !
  !  Available Shapes:
  !
  !     - Gauss  (Ok)
  !     - CosSqr (Ok)
  !     - CosSqF (Ok)
  !

  !  REAL(DPK) FUNCTION E_T( T, PULSETYPE, E0, TAU, OMEGA, PHI)
  !REAL(DPK)         :: E0, TAU, PHI, OMEGA 
  !CHARACTER(LEN=6)     ::PULSETYPE

  REAL(DPK) FUNCTION E_T(T)

    USE PRECISION
    USE UNITS
    !
    IMPLICIT NONE
    !
    REAL(DPK)             T 
    !
    REAL(DPK)     SHAPE, SHAPE_DERIV, T1, T2
    INTEGER       i
    !     .. Executable Statement ..

    !
    
      SHAPE = 0.D0
      SHAPE_DERIV = 0.D0

      
      composite_pulse: IF(pulsetype.EQ.'composite') THEN 

         e_t = 0.0_dpk
         build_the_pulse:DO i = 1, n_modes
            e_t = e_t + e_i(i) * cos( w_i(i)* t + p_i(i) )
         ENDDO build_the_pulse
         RETURN

      ENDIF composite_pulse




      !     Gaussian Pulse ( Tau = Half-width at 1/e )


      gauss_pulse:IF ( PULSETYPE.EQ.'agauss' ) THEN

	 SHAPE = DEXP( - 0.5D0 * ( T / TAU ) * ( T / TAU ) )
	 SHAPE_DERIV = - T / ( TAU * TAU ) * SHAPE

      ENDIF gauss_pulse


      !     Cosine-square pulse ( Tau = Total width )

      cos2_pulse:IF ( PULSETYPE.EQ.'acos2' ) THEN          ! -ta/2 < t < ta/2

         IF ( DABS( 2.D0 * T ).LT.TAU ) THEN

            SHAPE       = DCOS( M_PI * T / TAU ) ** 2
            SHAPE_DERIV = - (M_PI / TAU ) * DSIN( 2.0D+00* M_PI * T / TAU )

         ENDIF

      ENDIF cos2_pulse

!     Sine-square pulse ( Tau = Total width )

      sin2_pulse:IF ( PULSETYPE.EQ.'asin2' ) THEN          ! 0 < t < ta

         IF ( DABS( T ).LT.TAU ) THEN
            SHAPE       = DSIN( M_PI * T / TAU ) ** 2
            SHAPE_DERIV = (M_PI / TAU ) * DSIN( 2.0D+00* M_PI * T / TAU )
         ENDIF   
      ENDIF sin2_pulse



!
!     Flat pulse with CosSqr turn on and off :
!
!        - CosSqr turn on to -T1 with a time duration T2 
!          ( half the duration of an entire CosSqr pulse ) 
!        - Constant amplitude ( = E0 ) from -T1 to T1
!        - CosSqr turn from T1 with a time duration T2 
!
!     T1 and T2 are coded in TAU as the following :
!        - T1 = integer part of TAU / 100 ( 2 decimals )
!        - T2 = decimal part of TAU * 10000 ( Maximum 9999 u.a.)
!
!     Example : TAU = 12504.002084 => T1 = 125.04, T2 = 20.84
!

      cos2flat_pulse:IF (PULSETYPE.EQ.'acos2flat') THEN

	 T1 = DINT( TAU ) / 1.D2

	 T2 = ( TAU - DINT( TAU ) ) * 1.D4

         T2 = 2.D0 * T2


         SHAPE = 0.D0
         SHAPE_Deriv = 0.D0 

         IF ( DABS( T ).LT.( T1 + T2 / 2.D0 ) ) THEN


            IF ( DABS( T ).LT.T1 ) THEN

               SHAPE = 1.D0

            ELSE

               SHAPE = DCOS( M_PI * ( DABS( T ) - T1 ) / T2 ) ** 2
               
               IF ( T.LT.-T1 ) THEN

                  SHAPE_DERIV = - 2.D0 * M_PI / T2               &
                              *  DCOS( M_PI * ( T + T1 ) / T2 )  &
                              * DSIN( M_PI * ( T + T1 ) / T2 )
               ENDIF

               IF ( T.GT.T1 ) THEN

                  SHAPE_DERIV = - 2.D0 * M_PI / T2               &
                              * DCOS( M_PI * ( T - T1 ) / T2 )   &
                              * DSIN( M_PI * ( T - T1 ) / T2 )
               ENDIF

            ENDIF
         ENDIF

      ENDIF cos2flat_pulse
!

      IF (PULSETYPE.EQ.'esin2') THEN

         IF ( DABS( T ).LT.TAU )  SHAPE = 0.5_dpk * ( 1 - DCOS( 2* M_PI * T/TAU ) )
         

         E_T = E0 * SHAPE * DCOS( w0 * T + PHI )

      ELSE

         E_T = E0 * (                SHAPE * DSIN( w0 * T + PHI )    &
                     - SHAPE_DERIV * DCOS( w0 * T + PHI )/w0 )

      ENDIF


      RETURN

    END FUNCTION E_T

    !     ..
    !
    !  Purpose
    !  =======
    !
    !  E_T returns the value of an electromagnitc field E( T )
    !  defined as :
    !
    !     E( T ) = - D / DT A( T )
    !
    !  where A( T ) is defined as    
    !
    !  E( T ) = E0 * F( T, TAU ) * COS( OMEGA * T + PHI ) 
    !  
    !    F(T,TAU) = SIN(M_PI*T/TAU)**2 = 1/2 *( 1 - COS(OMEGA_P * T ) ) 
    !
    !        OMEGA_P = 2*PI/TAU
    !
    !  and where F( T, TAU ) represents the shape of the pulse and is specified by
    !  PULSETYPE.
    

  END MODULE pulse
!EOF!
