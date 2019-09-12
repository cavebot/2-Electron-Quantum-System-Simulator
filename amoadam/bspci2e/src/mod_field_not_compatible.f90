MODULE pulse
  !
  USE PRECISION
  !
  IMPLICIT NONE
  !
  PUBLIC
  !pulse!

  REAL(dpk)                 :: E0
  REAL(dpk)                 :: omega
  REAL(dpk)                 :: tau
  REAL(dpk)                 :: phi
  CHARACTER(len=10)         :: pulseType

  !!
  !!  REAL(DPK) FUNCTION E_T( T, PULSETYPE, E0, TAU, OMEGA, PHI)
  !!

CONTAINS
  
  SUBROUTINE READ_FIELD
    !
    USE units, ONLY:i0_td, enau, m_pi
    USE io,    only:ninp
    !    
    IMPLICIT NONE
    !
    CHARACTER(LEN=25)      :: field_file = 'tinp/pulse.inp'
    REAL(dpk)              :: i0   ! e0 = sqrt(i0): electric field
    INTEGER                :: cycles
    !
    WRITE(*,*) '# opening field file :', field_file

    !
    OPEN( ninp, file=TRIM(field_file),status='old')
    READ( ninp, * ) i0                   ! peak intensity in W/cm^2  
    READ( ninp, * ) omega                ! field frequency in eV
    READ( ninp, * ) cycles        
    READ( ninp, * ) phi
    READ( ninp, * ) pulsetype ! pulse shape ('gauss', 'cos2' 'cos2flat', 'sin2')
    CLOSE(ninp)

    e0    = SQRT(i0/i0_td)                      ! E_o : W/cm^2 -->     (a.u.)
    omega = omega/enau                          ! Wph :     eV -->     (a.u.)
    tau   = cycles * (2.0*m_pi/omega)           ! T_L : pulse duration (a.u.)
    
    !
    WRITE(*,*) '#  field file read done.'
  END SUBROUTINE READ_FIELD




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

  REAL(DPK) FUNCTION E_T(T)

    USE PRECISION
    USE UNITS
    !
    IMPLICIT NONE
    !
    REAL(DPK)             T, E0, TAU, PHI, OMEGA 
    CHARACTER(LEN=6)      PULSETYPE
    REAL(DPK)             SHAPE, SHAPE_DERIV, T1, T2
    !     .. Executable Statement ..
    !
    
      SHAPE = 0.D0
      SHAPE_DERIV = 0.D0



      !     Gaussian Pulse ( Tau = Half-width at 1/e )


      gauss_pulse:IF ( PULSETYPE.EQ.'gauss' ) THEN

	 SHAPE = DEXP( - 0.5D0 * ( T / TAU ) * ( T / TAU ) )
	 SHAPE_DERIV = - T / ( TAU * TAU ) * SHAPE

      ENDIF gauss_pulse


      !     Cosine-square pulse ( Tau = Total width )



      cos2_pulse:IF ( PULSETYPE.EQ.'cos2' ) THEN          ! -ta/2 < t < ta/2

         IF ( DABS( 2.D0 * T ).LT.TAU ) THEN
            SHAPE       =                   DCOS(          M_PI * T / TAU ) ** 2
            SHAPE_DERIV = - (M_PI / TAU ) * DSIN( 2.0D+00* M_PI * T / TAU )
         ENDIF

      ENDIF cos2_pulse

!     Sine-square pulse ( Tau = Total width )

      sin2_pulse:IF ( PULSETYPE.EQ.'sin2' ) THEN          ! 0 < t < ta

         IF ( DABS( T ).LT.TAU ) THEN
            SHAPE       =                 DSIN(          M_PI * T / TAU ) ** 2
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

      cos2flat_pulse:IF (PULSETYPE.EQ.'cos2flat') THEN

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
      E_T = E0 * ( SHAPE * DSIN( OMEGA * T + PHI )    &
           - SHAPE_DERIV / OMEGA * DCOS( OMEGA * T + PHI ) )


!      E_T = E0 *  SHAPE * DSIN( OMEGA * T + PHI )    
!

      RETURN

    END FUNCTION E_T

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

    REAL(DPK) FUNCTION A_T(T)

      !
      USE PRECISION
      USE UNITS
      !
      IMPLICIT NONE
      !
      REAL(DPK)             T, E0, TAU, PHI, OMEGA 
      CHARACTER(DPK)        PULSETYPE
      REAL(DPK)             SHAPE, T1, T2
        
      !     .. Executable Statement ..

      SHAPE = 0.D0


      !     Gaussian Pulse ( Tau = Half-width at 1/e)


      gauss_pulse:IF ( PULSETYPE.EQ.'gauss' ) THEN

         SHAPE = DEXP( - 0.5D0 * ( T / TAU ) * ( T / TAU ) )
      ENDIF Gauss_pulse


      !     Cosine-squared pulse ( Tau = Total width )

      cos2_pulse:IF ( PULSETYPE.EQ.'cos2' ) THEN
         IF ( DABS( 2.D0 * T ).LT.TAU ) THEN
            SHAPE = DCOS( M_PI * T /TAU ) ** 2
         ENDIF   
      ENDIF cos2_pulse


!     sine-squared pulse ( Tau = Total width )

      sin2_pulse:IF ( PULSETYPE.EQ.'sin2' ) THEN
         IF ( DABS( T ).LT.TAU ) THEN
            SHAPE = DSIN( M_PI * T /TAU ) ** 2
         ENDIF
      ENDIF sin2_pulse


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

      cos2flat_pulse:IF (PULSETYPE.EQ.'cos2flat') THEN

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
      A_T = E0 * SHAPE / OMEGA * DCOS( OMEGA * T + PHI )
      !
      RETURN
    END FUNCTION A_T

  END MODULE pulse
!EOF!
