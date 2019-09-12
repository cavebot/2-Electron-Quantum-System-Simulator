C----- PROGRAM modified to save memory  94.5.2 --------------------
C****   NCS   : maximum number of configuration series    |nll' >
C****   NHX   : Total number of configuration             |nln'l' >
C****   NS-2  : Maximum number of 1e orbitals in dmx1e.inp files 
C****   NWF2E : 
C****   NL1E  : ( <= NL )
      PROGRAM DMX2E
C -----------------------------------------------------------------
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

      INCLUDE "parameter.2e.inc"
C      PARAMETER(NS=802, NL1E=5, NHX=2500, NCS=30, NWF2E=600)

      DIMENSION VZ(NL1E-1,NS-2,NS-2),DZ(NS-2,NS-2), DMXZ(NHX,NHX),ANG(4)
      DIMENSION CIIN(NWF2E,NHX),CFIN(NWF2E,NHX),CI(NHX),CF(NHX),
     1  ENI(NWF2E),ENF(NWF2E),DMX(NWF2E),TL1(NHX)
      DIMENSION NHFI(NCS),LHFI(NCS),LLI(NCS),NMINI(NCS),NMXI(NCS),
     1  NDI(NCS),NHFF(NCS),LHFF(NCS),LLF(NCS),NMINF(NCS),NMXF(NCS),
     2  NDF(NCS),ID(4)
      COMMON/ANG/J1,J2,LS
      character*16 gauge
      CHARACTER*100 ARGV
C.........................................

      CALL GETARG(1, ARGV)
      READ(ARGV,*) L
      CALL GETARG(2, ARGV)
      READ(ARGV,*) GAUGE

      WRITE(*,*) '# d2eb     initial wave        L =', INPL
      WRITE(*,*) '# d2eb     gauge               G =', gauge

C      open(9,file='inp/dmx2e.inp')
C      read(9,*) L, gauge
C      close(9)


C************************************************************************
C*
C*       Read        radial 1e-dme :   
C*                                          _  
C*       d(nl,n^l^) = < nl | d | n^l^ > = _/ dr Pnl(r) d/dr( Pn^l(r)) 
C*
C*       d_r = vz / dw( Rydberg )
C*
C*       therefore 
C*       vz          are in   Rydberg scale
C*
C*       vz(l,n,n^) == < n l | d | n^ l + 1 > 
C*                                                    l.n.n.   12/03/2000
C*************************************************************************  

C# length
      inpmode = 1
      if(gauge.eq.'v')  inpmode = 0


      DO  I = 1, NL1E - 1


         call dmx1efile(3, I - 1, inpmode) 

         READ(3) MODE
         READ(3) LB, LZ, NWFB, NWFZ
         READ(3) ((VZ(I,IC,IR), IC=1,NWFZ), IR = 1,NWFB)

         if(inpmode.ne.mode) then 

            write(*,*)'Incosistent gauge in dmx1e file and input gauge.'
            stop
         endif

         CLOSE(3)

      ENDDO
      

C#
C#   initial angular symmetry
C#

      LOI = l

      call wf2efile(1, loi) 

      READ(1)    LOI, LSI
      READ(1)    NCMXI
      READ(1)  ( NHFI(K), K = 1, NCMXI)
      READ(1)  ( LHFI(K), K = 1, NCMXI)
      READ(1)  ( LLI(K),  K = 1, NCMXI)
      READ(1)  ( NMINI(K),K = 1, NCMXI)
      READ(1)  ( NMXI(K), K = 1, NCMXI)
      READ(1)  ( NDI(K),  K = 1, NCMXI)
      READ(1)    NHMXI
      READ(1)    NSI


C#
C#   final angular symmetry
C#

      LOF = l + 1 

      call wf2efile(2,lof) 

      READ(2)    LOF, LSF
      READ(2)    NCMXF
      READ(2)  ( NHFF(K), K = 1, NCMXF)
      READ(2)  ( LHFF(K), K = 1, NCMXF)
      READ(2)  ( LLF(K),  K = 1, NCMXF)
      READ(2)  ( NMINF(K),K = 1, NCMXF)
      READ(2)  ( NMXF(K), K = 1, NCMXF)
      READ(2)  ( NDF(K), K = 1, NCMXF)
      READ(2)    NHMXF
      READ(2)    NSF

      
C..........................................

      open(16,file='out/dmx2e.out')

      call dmx2efile(17, loi, mode)


      J1 = LOF
      J2 = LOI

      DO  KR = 1,NHX
         DO   KC = 1,NHX

            DMXZ(KR,KC) = 0.0D+00
         ENDDO
      ENDDO

      IF(LSI .NE. LSF) STOP
      LS = LSI
      IF(LS .EQ. 0) PHASE =  1.0D0
      IF(LS .EQ. 1) PHASE = -1.0D0

C.......................................................................
      KF = 1
      DO 100 IR = 1,NCMXF
      KI = 1
      DO 101 IC = 1,NCMXI

         DO NR = 1, NS-2
            DO  NC = 1, NS-2            
               DZ(NR,NC) = 0.0D+00
            ENDDO
      ENDDO

      CALL ANGFC(LHFF(IR),LLF(IR),LHFI(IC),LLI(IC),LOF,LOI,ID,ANG)

C-------------------------------------------------------------------
      DO 20 I = 1,4

         IF(ID(I) .EQ. 0) GO TO 20
         IF(I .GT. 1) GO TO 25

         DO 30 NR = 1,NDF(IR)
            MR = NMINF(IR)+NR-1
            DO 30 NC = 1,NDI(IC)
               MC = NMINI(IC)+NC - 1
               IF(MC .NE. MR) GO TO 30

               IF(2*INT(ID(I)/2).eq.ID(I)) then

                  IDN = ID(I)/2

         IF(MODE.EQ.0) THEN 
            DZ(NR,NC) = DZ(NR,NC)+ANG(I)*(-VZ(IDN,NHFI(IC),NHFF(IR)))
         ELSE
            DZ(NR,NC) = DZ(NR,NC)+ANG(I)*(VZ(IDN,NHFI(IC),NHFF(IR)))
         ENDIF

      ELSE
         IDN = (ID(I)+1)/2
         DZ(NR,NC) = DZ(NR,NC)+ANG(I)*VZ(IDN,NHFF(IR),NHFI(IC))
      END IF

   30 CONTINUE
      GO TO 20
C.......................................................................
   25 IF(I .GT. 2) GO TO 26
      DO 31 NR=1,NDF(IR)
      MR=NMINF(IR)+NR-1
      IF(MR .NE. NHFI(IC)) GO TO 31
      DO 32 NC=1,NDI(IC)

      MC = NMINI(IC) + NC - 1

      IF(2*INT(ID(I)/2).EQ.ID(I)) THEN
 
         IDN=ID(I)/2

         IF(MODE.EQ.0) THEN 
            DZ(NR,NC)=DZ(NR,NC)+PHASE*ANG(I)*(-VZ(IDN,MC,NHFF(IR)))
         ELSE
            DZ(NR,NC)=DZ(NR,NC)+PHASE*ANG(I)*(VZ(IDN,MC,NHFF(IR)))
         ENDIF

      ELSE
         IDN=(ID(I)+1)/2
         DZ(NR,NC)=DZ(NR,NC)+PHASE*ANG(I)*VZ(IDN,NHFF(IR),MC)
      END IF

   32 continue
   31 CONTINUE
      GO TO 20
C ---------------------------------------------------------------
   26 IF(I .GT. 3) GO TO 27

      DO 33 NC=1,NDI(IC)
      MC=NMINI(IC)+NC-1
      IF(MC .NE. NHFF(IR)) GO TO 33

      DO 34 NR=1,NDF(IR)
         MR = NMINF(IR) + NR - 1

      IF(2*INT(ID(I)/2).EQ.ID(I)) THEN

         IDN=ID(I)/2

         IF(MODE.EQ.0) THEN 
            DZ(NR,NC) =DZ(NR,NC) + PHASE*ANG(I)*(-VZ(IDN,NHFI(IC),MR))
         ELSE
            DZ(NR,NC) =DZ(NR,NC) + PHASE*ANG(I)*(VZ(IDN,NHFI(IC),MR))
         ENDIF
      
      ELSE

      IDN=(ID(I)+1)/2
      DZ(NR,NC) = DZ(NR,NC) + PHASE*ANG(I)*VZ(IDN,MR,NHFI(IC))

      END IF

   34 CONTINUE
   33 CONTINUE
      GO TO 20
C --------------------------------------------------------------
   27 IF(NHFF(IR) .NE. NHFI(IC)) GO TO 20
      DO 35 NR=1,NDF(IR)
      MR=NMINF(IR)+NR-1
      DO 35 NC=1,NDI(IC)
      MC = NMINI(IC) + NC - 1

      IF(2*INT(ID(I)/2).EQ.ID(I)) THEN
         IDN = ID(I)/2
         IF(MODE.EQ.0) THEN 
            DZ(NR,NC) = DZ(NR,NC) + ANG(I) * (-VZ(IDN,MC,MR))
         ELSE
            DZ(NR,NC) = DZ(NR,NC) + ANG(I) * (VZ(IDN,MC,MR))
         ENDIF

      ELSE
         IDN = (ID(I) + 1) / 2
         DZ(NR,NC) = DZ(NR,NC) + ANG(I)*VZ(IDN,MR,MC)
      END IF
 35   CONTINUE
 20   CONTINUE
C-------------------------------------------------------------------
      IF(IR.NE.1 .OR. IC.NE.1) GO TO 95

   95 CONTINUE 

      CALL TRNSMX(DMXZ,DZ,KF,KI,NDF(IR),NDI(IC))

      KI = KI + NDI(IC)

  101 CONTINUE
      KF = KF + NDF(IR)

  100 CONTINUE
C-------------------------------------------------------------------
      DSQ2 = DSQRT(2.0D0)

      IPI  = LHFI(1) + LLI(1)
      IP   = MOD(IPI,2)
      IT   = 1
      IF(IP.EQ.1) GO TO 65

      DO 60 IC = 1,NCMXI
        IF(LHFI(IC) .NE. LLI(IC))   GO TO 60
        IF(NMINI(IC) .GT. NHFI(IC)) GO TO 60
        DO 61 NR = 1, NHX
 61       DMXZ(NR,IT) = DMXZ(NR,IT) / DSQ2
 60       IT = IT + NDI(IC)
          GO TO 70

   65 CONTINUE
      DO 62 IR = 1,NCMXF
      IF(LHFF(IR) .NE. LLF(IR)) GO TO 62
      IF(NMINF(IR) .GT. NHFF(IR)) GO TO 62
      DO 63 NC = 1,NHX
   63 DMXZ(IT,NC) = DMXZ(IT,NC)/DSQ2
   62 IT = IT + NDF(IR)

   70 CONTINUE

C...................................

      DO 45 NSC = 1, NSI
      READ(1) ENI(NSC)
 45   READ(1) (CIIN(NSC,K), K = 1, NHMXI)

C...................................

      DO 46 NSR = 1, NSF
      READ(2) ENF(NSR)
 46   READ(2) (CFIN(NSR,K), K=1,NHMXF)

C....................................
      CLOSE(1)
      CLOSE(2)

C*********************************************************
C*                  UNIT  17
C*********************************************************
C*    Write in binary mode  OSLL+1.DAT L = 0,1,2 
C*    
C*    Storing Form : 
C* 
C*    Lin, Lfin, Nin, Nfin
C*    En Initial
C*    En Final states 
C*    DME( Lin -> Lfin)
C* 
C*    Note:
C*    These dataFiles are read from reados.f programm
C********************************************************

      WRITE(17) MODE
      WRITE(17) LOI,LOF,NSI,NSF
      WRITE(17) (ENI(K),K=1,NSI)
      WRITE(17) (ENF(K),K=1,NSF)


C************* Unit 16 ***********************************
C*    Write in ascii mode  OSLL+1.DAT L = 0,1,2 
C*    in similar form as in binary mode
C*    For testing purposes
c
   55   format(4i5)

        WRITE(16,55) MODE
        WRITE(16,55) LOI,LOF,NSI,NSF
        WRITE(16,5) (ENI(NSC), NSC = 1, NSI )
        WRITE(16,5) (ENF(NSR), NSR = 1, NSF )
C********************************************************

      DO 180 NSC = 1, NSI


      DO  KC = 1, NHMXI

      CI(KC) = CIIN(NSC,KC)

      ENDDO
C#
C#    Make the matrix-vector multiplication
C#

      CALL DGEMV('N',NHMXF,NHMXI,1.0d0,DMXZ,NHX,CI,1,0.0d0,TL1,1)

      DO NSR = 1, NSF

      DO KR  = 1, NHMXF
       CF(KR) = CFIN(NSR,KR)
      ENDDO

C#
C#    Make the vector-vector multiplication
C#

      DMX(NSR)  = DDOT(NHMXF,CF,1,TL1, 1)

      ENDDO

      WRITE(17)   (DMX(NSR), NSR = 1, NSF)
      WRITE(16,5) (DMX(NSR), NSR = 1, NSF)

 180   CONTINUE

       close(16)
       close(17) 

C.......................................................................
C*    Format Statements
 577  FORMAT(5A16)  
    1 FORMAT(/2X,'THE INITIAL STATE OF THE TRANSITION --')
    2 FORMAT(/2X,'THE FINAL STATE OF THE TRANSITION --')
    3 FORMAT(2X,'TOTAL L AND TOTAL S=',2I3,4x,'# of Config. is ',i4)
    4 FORMAT(/2X,'ENERGY EIGENVALUES OF INITIAL & FINAL STATES')
    5 FORMAT(2X,1P8E15.7)
c    6 FORMAT(/2X,'OS/(2*LOF+1): Length - top & Velocity - bottom'/
    6 FORMAT(/2X,'OS Values: Length - top & Velocity - bottom.   {LF is
     & included. For Absorp. LF=2*LOF+1. For Emiss.("-" sign) LF=2*LOI+1
     1 }' /2X,'ROW - FINAL STATES & COLUMN - INITIAL STATES')
    7 FORMAT(2X)
    8 FORMAT(/2X,'LHFF,LLF,LHFI,LLI,ID(I),ANG(I) --'/)
    9 FORMAT(2X,4I2,3X,4I2,1P4E12.4)
   11 FORMAT(/2X,'EXCITATION ENERGY IN RYD.')
C.......................................................................

       END

C########################################################################
C#EOF












