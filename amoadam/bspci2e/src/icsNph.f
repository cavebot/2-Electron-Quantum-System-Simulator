C#                                                                       
C#  prepare data files to be used by the csNph.f program
C#   --
C#   -- requires date from the d2eb.f  program
C#   -- requires data from the  n2e.f  program
C#
C#
      PROGRAM CSNPH_I

      IMPLICIT NONE 
      INTEGER MaxDim, MaxL, MaxFile, MODE
      PARAMETER( MaxDim = 2000, MaxL = 4, MaxFile = MaxL )

      INTEGER NBIN, NINP

      INTEGER i, l, n, nn1, nn2, j, lmax, lmin, NDimension
      INTEGER NbInit, NbFin, Linit, Lfin,nbphot, inp_mode

      REAL*8 de,norm

       LOGICAL IncludeAngul
       CHARACTER*20 DipoleFile, enFile, OutFile

       REAL*8 en_dip( MaxDim, 0:MaxL ),en_nrm( MaxDim, 0:MaxL )
       REAL*8 an( MaxDim, 0:MaxL ), phase(MaxDim,0:MaxL)
       REAL*8 p_sign( MaxDim, 0:MaxL )
       REAL*8 dmx( MaxDim, MaxDim, 0:MaxL)
c........................................

       CHARACTER*100 ARGV
       CHARACTER*6  GAUGE, SYSTEM
c       CHARACTER*1 SMODE
C.........................................

      NINP = 9
      NBIN = 1

C.............................



C GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT

      CALL GETARG(1, ARGV)
      READ(ARGV,*) NBPHOT
      CALL GETARG(2, ARGV)
      READ(ARGV,*) GAUGE
      CALL GETARG(3, ARGV)
      READ(ARGV,*) SYSTEM



c      IF(SMODE.EQ."v") THEN              
c         INP_MODE = 0
c      ELSE IF(SMODE.EQ."l") THEN         
c         INP_MODE = 1
c      ELSE  
c         INP_MODE = 1
c         WRITE(*,*) '# csNph_l::   only the inputs l or v are available'
c         WRITE(*,*) '# csNph_l::        set default mode, length gauge '
c      ENDIF
      
      WRITE(*,*) '# csNph_l::        gauge =  ', gauge
      
C...........................................

C#  to be changed for the general case Nph > 4

       DipoleFile = 'dat/dmx.dat'
          enFile  = 'dat/en.dat'
      OutFile     = 'out/icsNph.out'

C.........................................

      WRITE(*,*) "#  ICSNPH PROGRAM EXECUTION "
      WRITE(*,*) "#"
      WRITE(*,*) "# icsNph ::                system = ", system
      WRITE(*,*) "# icsNph ::         # of  photons = ", nbphot
      WRITE(*,*) "# icsNph ::                 gauge = ", gauge
      WRITE(*,*) "#"
      WRITE(*,*) "# icsNph :: 2-e dipole matrix elements file:"
      WRITE(*,*) "#"
      WRITE(6,*) "# icsNph ::           energy file = ", enFile
      WRITE(6,*) "# icsNph ::           dipole file = ", DipoleFile

C
      NDimension = 0
      LMin = MaxL
      LMax = 0

      
      DO L = 0, NBPHOT-1

        CALL dmxfile(nbin, "dat", system, "bin", gauge, l, l+1)

!         CALL DMX2EFILE(NBIN,L,INP_MODE)

         READ( NBIN ) mode
         READ( NBIN ) Linit, Lfin, Nbinit, Nbfin

         NDimension = MAX( NDimension, Nbinit, Nbfin )
               LMin = MIN( LMin, Linit, Lfin )
               LMax = MAX( LMax, Linit, Lfin )

         WRITE(*,*), Linit, Lfin, Nbinit, Nbfin, LMin, LMax


         READ( NBIN ) ( en_dip( i, Linit ), i = 1, Nbinit )
         READ( NBIN ) ( en_dip( i, Lfin ),  i = 1, Nbfin )


         IF ( Linit.LT.Lfin ) THEN

            DO i = 1, Nbinit
               READ( NBIN ) ( dmx( i, j, Lfin ), j = 1, Nbfin )
            ENDDO   

         ELSE
            DO i = 1, Nbinit
               READ( NBIN ) ( dmx( j, i, Lfin ), j = 1, Nbfin )
            ENDDO      
         ENDIF

         CLOSE( NBIN )
      ENDDO

      WRITE(*,*) "# csNph_i::             dmx2e data file read. "
C#   read normalization factor and phase for each 2e state

      DO L = 0, NBPHOT 
 
        CALL hfile(nbin,  "dat","n2e" ,"bin"  , l)       

         READ( nbin ) Nbfin, Lfin

         if(Lfin.ne.L) then
            WRITE(*,*) "# csNph_i:: norm data file is not       "
            WRITE(*,*) "# csNph_i:: the appropriate one.        "
            WRITE(*,*) "# csNph_i::                      Lfin = ",Lfin
            WRITE(*,*) "# csNph_i::                      Lfin = ",L
         else
            WRITE(*,*) "# csNph_i:: norm data  for          L = ",Lfin
            WRITE(*,*) "# csNph_i:: norm data       NSTATES_L = ",nbfin
         endif

         DO i = 1, Nbfin
            READ( nbin ) en_nrm(i,Lfin), phase(i, Lfin), an( i, Lfin )
         ENDDO

         CLOSE( nbin )

      ENDDO
      

C      IncludeAngul = .true.
      IncludeAngul = .true.



C#   save data in dat/dmx2e.dat


      OPEN( nbin, FILE = DipoleFile, FORM='UNFORMATTED' )

      WRITE( *, * ) mode
      WRITE( *, * ) Lmin, Lmax, NDimension


      WRITE( nbin ) mode
      WRITE( nbin ) Lmin, Lmax, NDimension, IncludeAngul
C#
      DO l = Lmin, Lmax - 1

         DO nn1 = 1, NDimension
            DO nn2 = 1, NDimension

               WRITE( nbin ) dmx( nn1, nn2, l + 1 )               

            ENDDO
         ENDDO
  
      ENDDO

      CLOSE( nbin )

C#   save energy, phase and normalization phactor


      OPEN( 16, FILE = OutFile)
      WRITE(*,*) Lmin, Lmax, NDimension 

      OPEN( nbin, FILE = enFile, FORM='UNFORMATTED' )
      WRITE( nbin ) Lmin, Lmax, NDimension


C      Phase = 0.D+00

      DO l = Lmin, Lmax 
         DO n = 1, NDimension

C#    'density of states' normalization    sqrt(2/DE(a.u.)

            de = en_dip(n+1, l) - en_dip(n-1,l) 

            IF(N.EQ.1.OR.N.EQ.NDimension)   THEN

C# to avoid  en(n-1,l), en(n+1,l) out of bounds.

               norm = 1.0D+00 
            ELSE
!               norm = 2.0D+00/sqrt(de)
               norm = an(n,l)
            ENDIF
 
            WRITE(16,*) 0.5D+00*en_dip(n,l),norm,phase(n,l)
            WRITE(nbin) 0.5D+00*en_dip(n,l),norm,phase(n,l),an( n,l)


!            WRITE( nbin ) 0.5D+00*en_dip(n,l ),an( n,l),phase(n,l),norm
         ENDDO
      ENDDO

      CLOSE( nbin )
C
      END
C###########################################################################
