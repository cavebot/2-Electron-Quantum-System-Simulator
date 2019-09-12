C#
C#
C#  Purpose
C#  =======
C#
C#  Evaluate the harmonic spectrum produced by an atom iluminated by an
C#  laser field. This spectrum is calculated using the Fourier 
C#  transform of the interferometric ionization yield
C#  
C#     author     :  Eric Cormier
C#      
C#     2001.12.05 : modified to accomodate Spline interpolations using 
C#                  NAG library,
C#     2004.01.02 : modified to perfom DFT for time-delayed yield 
C#     
C#

      PROGRAM FOURIER

      IMPLICIT NONE
      INTEGER        MaxValues, MaxPts, Nsplines, NWRK
      PARAMETER    ( MaxValues = 16800, MaxPts = 8192 )
      PARAMETER    ( Nsplines = MaxValues + 4, NWRK= 6 * MaxValues + 16)
      REAL*8         PI
      PARAMETER    ( PI = 3.141592654D+00 )
      REAL*8         ENAU
      PARAMETER    ( ENAU = 27.211396181D+00 )


      REAL*8             MaxFrequency, Tmin, Tmax, Delta, PhotonEnergy,
     $                   MaxEnergy, Width, Period, NbtimesDuration,
     $                   Duration
      CHARACTER*20       FileDipole, FileHarmonics, FileOutput

      INTEGER            NPow, i, NbPoints, NbValues, BreakLength
      INTEGER            NbPeriodWin 

      REAL*8             Time( MaxValues ) 
      REAL*8             at( MaxValues ),vt( MaxValues ),rt( MaxValues )
      REAL*8             e_t( MaxValues ),a_t( MaxValues )

      REAL*8                X( MaxPts    ),    Dip( MaxPts    )
      REAL*8                W( MaxPts )
      REAL*8             Spline(Nsplines), Break(Nsplines), WRK(NWRK)
      INTEGER            IFAIL
      COMPLEX*16         DipFFT( MaxPts )
      COMMON             /WORKSP/  RWKSP
      REAL*8             RWKSP(900000)
      REAL*8             DCSDER, FUNC
      EXTERNAL           FFT, FUNC

C......... NAG  ROUTINES

      EXTERNAL           E01BAF, E02BBF


C.....................................

C#
C#     Read the input parameters
C#


      OPEN( UNIT = 1, FILE ='tinp/tdse_hhg.inp', STATUS = 'UNKNOWN' )
      READ( 1, FMT = '( A20 )' ) FileDipole
      READ( 1, FMT = '( A20 )' ) FileHarmonics
      READ( 1, FMT = '( A20 )' ) FileOutput
      READ( 1, * )               MaxEnergy
      READ( 1, * )               PhotonEnergy
      READ( 1, * )               NbPeriodWin
      READ( 1, * )               NbtimesDuration
      CLOSE(1)

C.................................



      OPEN( UNIT = 2, FILE ='tout/tdse_hhg.out', STATUS = 'UNKNOWN' )

C#
C#     Read the function to be transformed. 
C#     data in time domain:
C#     t, f(t)
C#
C#

      OPEN( UNIT = 1, FILE = FileDipole, STATUS = 'OLD' )

      i = 1
C 1    READ( 1, *, END = 2 ) time( i ), dipole(i), acc1(i)

  1    READ( 1, *, END = 2 ) time(i), at(i),vt(i),rt(i),a_t(i),e_t(i)



      i = i + 1

      IF ( i.GE.MaxValues ) THEN

         WRITE( 6, * )' tdse_hhg : Number of input values too large. '
         WRITE( 6, * )' MaxValues = ',MaxValues
         WRITE( 6, * )' Program stopped '
         STOP

      ENDIF   

      GO TO 1

 2    NbValues = i - 1

      CLOSE( 1 )

C...............................


      OPEN( UNIT = 2, FILE ='tdat/observables.dat',STATUS = 'UNKNOWN' )
        WRITE(2,'(a1,a60)') '&','wf calculated in vel. gauge'
        WRITE(2,'(a1,8a20)')
     & '&','t','at=dv/dt-ZqE(t)',
     &      'vt=<psi|v|psi>+ZA(t)',
     &      'rt=<psi|r|psi>',
     &      'A(t)',
     &      'E(t)=-dA(t)/dt)'
      DO i = 1, NbValues
        at(i) = at(i) - e_t(i)
        vt(i) = vt(i) + a_t(i)
        WRITE(2,'(6E20.8)') time(i),at(i),vt(i),rt(i),a_t(i),e_t(i)
      ENDDO
      CLOSE(2) 


C#
C#     Interpolate real and imaginary at particular values 
C#

      Duration     = Time( NbValues ) - Time( 1 )
      Tmin         = Time( 1 ) - NbtimesDuration * Duration / 2.0D+00
      Tmax         = Time( NbValues ) + NbtimesDuration * Duration/2.D0
      MaxFrequency = MaxEnergy / ( 2.0D+00 * PI )
      Delta        = 0.5D+00 / MaxFrequency
      NbPoints     = IDINT( ( Tmax - Tmin ) / Delta )

   
C#
C#     Calculate the nearest power of two
C#


      NPow         = IDINT( DLOG( DFLOAT( NbPoints ))/DLOG( 2.D0 )) + 1
      NbPoints     = 2**NPow
      Delta        = ( Tmax - Tmin ) / DFLOAT( NbPoints - 1 )
      MaxFrequency = 0.5D+00 / Delta

      DO i = 1, NbPoints
         X( i ) = Tmin + ( i - 1 ) * Delta
      ENDDO


      WRITE( *, * ) '            Starting time : ', Tmin
      WRITE( *, * ) '            Ending time   : ', Tmax
      WRITE( *, * ) 'Number of sampling points : ', NbPoints 
      WRITE( *, * ) '     Maximum frequency    : ', MaxFrequency 
      WRITE( *, * ) ' Maximum energy (a.u.)    : ', MaxFrequency*2.D0*PI
      WRITE( *, * ) '        Frequency spacing : ', Delta
      WRITE( *, * ) 'Window cut (Nb of cycles) : ', NbPeriodWin 
      CLOSE( 2 )


C#
C#     Multiply the input data by a window function
C#

      Period = 2.D0 * PI / PhotonEnergy 

      if(NbPeriodWin.ne.0) then 

         Width = DBLE( NbPeriodWin ) * Period

         write(*,*) '# window is implemented. '  
         WRITE(*,* )'#               cw = ', NbPeriodWin 
         write(*,*) '#               tw = ', width  
         write(*,*) '#                T = ', period
         
         DO i = 1, NbValues
                
           IF(TIME(I+1).LT.TIME(I)) THEN
              WRITE(*,*) i+1,i, time(i+1),time(i)
              cycle
           ENDIF
     
           at( i ) = at( i )
     s             * FUNC(Time( i ),Width,Time( 1 ),Time( NbValues ) )
         ENDDO
      
      endif
C.........................................



C      CALL SplineNat( Time, DipoleRe, NbValues, X, Re, NbPoints,
C     &                MaxValues, MaxPts )      
C      CALL SplineNat( Time, DipoleIm, NbValues, X, Im, NbPoints,
C     &                MaxValues, MaxPts )   


C.......................................


C#
C#         DIPOLE(I)  INTERPOLATION
C#


      IFAIL = 0 

      CALL E01BAF(NbValues, Time, at ,Break, Spline, Nsplines, 
     1     WRK, NWRK, IFAIL)


      DO i = 1, NbPoints


         Dip( i ) = 0.0D+00

         IF ( ( X( i ).GT.Time( 1        ) ).AND.
     c           ( X( i ).LT.Time( NbValues ) ) ) THEN


        CALL E02BBF( NbValues+4, Break, Spline, X( i ), Dip( i ), IFAIL)
C     1           IFAIL) 

         ENDIF
      ENDDO



C...
C...
C...
C...
C...

      OPEN( UNIT = 22, FILE ='tdat/dip_window.dat', STATUS = 'UNKNOWN' )
      DO i = 1, NbPoints
         DipFFT( i ) = DCMPLX( Dip( i ), 0.0D+00 )
         WRITE( 22, * ) X( i ), Dip( i )
      ENDDO
      CLOSE(22) 



C#
C#     Calculate the Fast Fourier transform 
C#

      CALL FFT( DipFFT, NbPoints, 1 )



C...........................................

C#
C#     Rescale the frequency in harmonics number
C#

      DO i = 1, NbPoints / 2 + 1
         W( i ) = ( i - 1 ) / ( Delta * NbPoints ) * 2.0D+00 * PI / 
     $            PhotonEnergy 

      ENDDO

C
      DO i = NbPoints / 2 + 2, NbPoints

         W( i ) = ( i - NbPoints - 1) / ( Delta * NbPoints ) * 2.0D+00 * 
     $            PI / PhotonEnergy

      ENDDO 


C...........................................

C#
C#     Renormalise the spectrum
C#         

      DO i = 2, NbPoints
         DipFFT( i ) = DipFFT( i ) * Delta / DSQRT( 2.D0 * PI )
      ENDDO


C#
C#     Write on disk the spectrum as a function of the harmonic number
C#

      OPEN( UNIT = 1, FILE =FileHarmonics, STATUS = 'UNKNOWN' )
      DO i = 2, NbPoints / 2 + 1
	 WRITE( 1, * ) W( i ), CDABS( DipFFT( i ) )**2.0D+00
      ENDDO
      CLOSE(1)
      


      END
C#
C####################################################################
C#
      SUBROUTINE FFT(data, nn, isign)

      INTEGER isign,nn
      DOUBLE PRECISION data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif

        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then

        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0

        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=dble(wr)*data(j)-dble(wi)*data(j+1)
            tempi=dble(wr)*data(j+1)+dble(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
C#######################################################################

*
* User defined function for the window
*
C#        * func(t, tw, t0, tf )
C#
C#         ( cos[ pi* ( t- (t0+wd) ) ] / 2*wd )^2             t0  < t < t0 + tw 
C# w(t) =  1.0D+00                                        t  + tw < t < tf - tw 
C#         ( cos[ pi* ( t- (tf-wd) ) ] / 2*wd )^2         tf - tw < t < tf      
C#

      REAL*8 FUNCTION FUNC( X, W, LB, RB )
C#
      REAL*8             PI
      PARAMETER        ( PI = 3.1415926549D0 ) 
C#
      REAL*8             X, W, LB, RB
C#LCL
      REAL*8             Left, Right, VALUE
      REAL*8             DCOS
      INTRINSIC          DCOS      
*
* Cosine Squared function 
*
       Left = LB + W
      Right = RB - W
*
      VALUE = 0.D0 
      IF ( ( X.GE.LB ).AND.( X.LE.Left ) ) THEN
         VALUE = DCOS( PI * ( X - LEFT) / ( 2.D0 * W ) )**2.D0
      ENDIF  
      IF ( ( X.GT.Left ).AND.( X.LT.Right ) ) VALUE = 1.D0 
      IF ( ( X.GE.Right ).AND.( X.LE.RB ) ) THEN
         VALUE = DCOS( PI * ( X - RIGHT ) / ( 2.D0 * W ) )**2.D0
      ENDIF  
*
      FUNC  = VALUE
*
      RETURN
      END
C#######################################################################
