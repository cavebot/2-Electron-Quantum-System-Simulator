
C#
C* 
C*             CPC/ SUBMITTED VERSION  
C*             BSPCI2E PACKAGE : 1.0
C*
C*             PROGRAM NUMBER 7          DMX2E
C*
C#LAAN2002
C#######################################################################
      PROGRAM d2eb

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      INCLUDE "parameter.2e.inc"
      PARAMETER ( DSQ2  = 1.414 213 562 4 D+00 )
      PARAMETER ( PI    = 3.141 592 653 589 793 D+00 )
      PARAMETER ( ENAU  = 27.211 396 181  D+00 )
      PARAMETER ( ALPHA = 1.0D+00/137.035 989 561 D+00 )

      integer l1e_max, l1e_max_i, l1e_max_f
      
      DIMENSION VZ(NL,NS-2,NS-2),DZ(NS-2,NS-2),DMXZ(NHX,NHX),ANG(4)
      DIMENSION CIIN(NWF2E,NHX),CFIN(NWF2E,NHX),CI(NHX),CF(NHX)
      DIMENSION ENI(NWF2E),ENF(NWF2E),DMX(NWF2E),TL1(NHX)
      DIMENSION NHFI(NCS),LHFI(NCS),LLI(NCS),NMINI(NCS),NMXI(NCS)
      DIMENSION  NDI(NCS),NHFF(NCS),LHFF(NCS),LLF(NCS),NMINF(NCS)
      DIMENSION NMXF(NCS),NDF(NCS),ID(4)
      COMMON/ANG/J1,J2,LS
      CHARACTER*6  GAUGE, SYSTEM
      INTEGER NINP, NOUT
      integer nd1e, nwf2e_i, nwf2e_f, nd2e
      CHARACTER*100 ARGV
C.EXE


      nd1e       = 9
      nwf2e_i    = 10
      nwf2e_f    = 11
      nd2e       = 12

      NINP       = 13
      NOUT       = 16

C.....................

      CALL GETARG(1, ARGV)
      READ(ARGV,*) L             ! 0,1,...
      CALL GETARG(2, ARGV) 
      READ(ARGV,*) gauge         ! l,v
      CALL GETARG(3, ARGV)
      READ(ARGV,*) system        ! he0,li1,be2

      nini = 1

      WRITE(*,*) '& d2eb     initial wave        L =', L
      WRITE(*,*) '& d2eb     gauge               G =', gauge
      WRITE(*,*) '& d2eb                   system  =', system


c      OPEN(NINP, FILE='inp/dmx2e.inp')
c      READ(NINP, *) L, GAUGE
c      READ(NINP, *) NINI
c      READ(NINP, *) NFIN
c      CLOSE(NINP)

C.......

      

C#   INITIAL SYMMETRY : L


      LOI = L

      CALL hfile(nwf2e_i,"dat","w2eb","bin",loi)              !initial states

!      CALL WF2EFILE(1, LOI) 

      READ(unit=nwf2e_i)    LOI, LSI
      READ(unit=nwf2e_i)    NCMXI
      READ(unit=nwf2e_i)  ( NHFI(K), K = 1, NCMXI)
      READ(unit=nwf2e_i)  ( LHFI(K), K = 1, NCMXI)
      READ(unit=nwf2e_i)  ( LLI(K),  K = 1, NCMXI)
      READ(unit=nwf2e_i)  ( NMINI(K),K = 1, NCMXI)
      READ(unit=nwf2e_i)  ( NMXI(K), K = 1, NCMXI)
      READ(unit=nwf2e_i)  ( NDI(K),  K = 1, NCMXI)
      READ(unit=nwf2e_i)    NHMXI
      READ(unit=nwf2e_i)    NSI



      WRITE(*,'(a10,30I3)') '# lhfi = ', (lhfi(k), k=1,ncmxi)
      WRITE(*,'(a10,30I3)') '#  lli = ', (lli(k),  k=1,ncmxi)


C#   FINAL SYMMETRY : L + 1

      LOF = L + 1 
      CALL hfile(nwf2e_f,"dat","w2eb","bin",lof)   !final states

      READ(unit=nwf2e_f)    LOF, LSF
      READ(unit=nwf2e_f)    NCMXF
      READ(unit=nwf2e_f)  ( NHFF(K), K = 1, NCMXF)
      READ(unit=nwf2e_f)  ( LHFF(K), K = 1, NCMXF)
      READ(unit=nwf2e_f)  ( LLF(K),  K = 1, NCMXF)
      READ(unit=nwf2e_f)  ( NMINF(K),K = 1, NCMXF)
      READ(unit=nwf2e_f)  ( NMXF(K), K = 1, NCMXF)
      READ(unit=nwf2e_f)  ( NDF(K),  K = 1, NCMXF)
      READ(unit=nwf2e_f)    NHMXF
      READ(unit=nwf2e_f)    NSF
      

      WRITE(*,'(a10,30I3)') '# lhfi = ', (lhff(k), k=1,ncmxf)
      WRITE(*,'(a10,30I3)') '#  lli = ', (llf(k),  k=1,ncmxf)


C*
C*      read 1-electron radial dmx (in ryd units)  
C*                                          _  


C# determine maximum 1-e q. angular number l_max = ?


C# search initial LOI symmetry

!      write(*,*) '&        lli = '
!      print*, lli(1)
   
      l1e_max_i = lli(1)
      do k = 2, ncmxi     
        if(lli(k).gt.l1e_max_i) l1e_max_i = lli(k)
      enddo

      do k = 2, ncmxi
        if(lhfi(k).gt.l1e_max_i) l1e_max_i = lhfi(k)
      enddo

      write(*,*) '&         Li = ', loi
      write(*,*) '&  l1e_max_i = ', l1e_max_i


C# search final LOF = LOI + 1 symmetry

      l1e_max_f = llf(1)
      do k = 2, ncmxf        
        if(llf(k).gt.l1e_max_f) l1e_max_f = llf(k)
      enddo

      do k = 2, ncmxf
        if(lhff(k).gt.l1e_max_f) l1e_max_f = lhff(k)
      enddo

      write(*,*) '&         Lf = ', lof
      write(*,*) '&  l1e_max_f = ', l1e_max_f

C      if(abs(l1e_max_f-l1e_max_i).ne.1) then
C         write(*,*)'& d2eb::  wrong 2e configuration files'
C         write(*,*)'& d2eb:: l1e_max_f-l1e_max_i = ',l1e_max_f-l1e_max_i
C         write(*,*)'& d2eb:: exiting.'
C         stop
C      else
C         l1e_max = max(l1e_max_i, l1e_max_f)
C      endif

         l1e_max = max(l1e_max_i, l1e_max_f)

      
      IF(GAUGE.EQ.'v')    INPMODE = 0
      IF(GAUGE.EQ.'l')    INPMODE = 1 
      IF(GAUGE.EQ.'a')    INPMODE = 2 


      DO  l = 1, l1e_max 
        CALL dfile(nd1e,l-1,l,'d1e-',gauge)
         READ(nd1e) MODE
         READ(nd1e) LB, LZ, NWFB, NWFZ
         READ(nd1e) ((VZ(l,IC,IR), IC=1,NWFZ), IR = 1,NWFB)

         IF(INPMODE.NE.MODE) THEN 
            WRITE(*,*) 'd2eb: inconsistent gauge mode in d1e files'
            write(*,*) 'd2eb:      mode = ', mode
            write(*,*) 'd2eb:   inpmode = ', inpmode
            STOP
         ENDIF

         CLOSE(nd1e)

         write(*,*) l, vz(l,1,1)
       ENDDO



C#
C#  calculate matrix element
C#


      OPEN(NOUT,FILE='out/dmx2e.out')


!     dmx file!
      CALL dmxfile(nd2e, "dat", system,"bin", gauge, loi, lof) 

      J1 = LOF
      J2 = LOI

      DO  KR = 1,NHX
         DO   KC = 1,NHX
            DMXZ(KR,KC) = 0.0D+00
         ENDDO
      ENDDO

      IF(LSI.NE.LSF) THEN 

         WRITE(*,*) ' ERROR IN INPUT DATA FILES : '
         WRITE(*,*) ' INITIAL SPIN DIFFERENT THAN FINAL SPIN : '
         WRITE(*,*) ' LSI = ', LSI
         WRITE(*,*) ' LSF = ', LSF
         WRITE(*,*) ' d2eb aborts. '

         STOP
      ENDIF

      LS = LSI


      phase = 1.0 - 2.0 * lsi

C# singlet
C      IF(LSI.EQ.0) PHASE =  1.0D+00
C# triplet
C      IF(LSI.EQ.1) PHASE = -1.0D+00

C.....

      KF = 1
      DO 100 IR = 1, NCMXF
         KI = 1
         DO 101 IC = 1,NCMXI

            DO NR = 1, NS - 2
               DO  NC = 1, NS - 2            
                  DZ(NR,NC) = 0.0D+00
               ENDDO
             ENDDO

             
          CALL ANGFC(LHFF(IR),LLF(IR),LHFI(IC),LLI(IC),LOF,LOI,ID,ANG)

C          write(*,'(6I5)') ir,ic, id

C.........................................

          
          DO 20 I = 1, 4

            IF(ID(I) .EQ. 0) GO TO 20
            IF(I .GT. 1)     GO TO 25

         DO 30 NR = 1, NDF(IR)

            MR = NMINF(IR) + NR - 1

            DO 30 NC = 1, NDI(IC)

               MC = NMINI(IC) + NC - 1

               IF(MC .NE. MR) GO TO 30

               IF(2*INT(ID(I)/2).EQ.ID(I)) THEN
C#
C#  lf-li < 0
C#
                  IDN = ID(I)/2

                  IF(MODE.EQ.0) THEN
                     DZ(NR,NC) = DZ(NR,NC)
     &                         + ANG(I) * ( -VZ(IDN,NHFI(IC),NHFF(IR) ))
                  ELSE
                     DZ(NR,NC) = DZ(NR,NC)
     &                         + ANG(I) * ( VZ(IDN,NHFI(IC),NHFF(IR) ) )
                  ENDIF
               ELSE

                  IDN = (ID(I)+1)/2
                  DZ(NR,NC) = DZ(NR,NC)+ANG(I)*VZ(IDN,NHFF(IR),NHFI(IC))
               END IF

 30         CONTINUE

            GO TO 20

C.........................................

 25         IF(I .GT. 2) GO TO 26

      DO 31 NR = 1, NDF(IR)

         MR = NMINF(IR) + NR - 1

         IF(MR .NE. NHFI(IC)) GO TO 31

         DO 32 NC = 1, NDI(IC)
            
            MC = NMINI(IC) + NC - 1

            IF(2*INT(ID(I)/2).EQ.ID(I)) THEN
 
               IDN=ID(I)/2

               IF(MODE.EQ.0) THEN 

                  DZ(NR,NC) = DZ(NR,NC) 
     &                      + PHASE * ANG(I) * ( - VZ(IDN,MC,NHFF(IR) ))

         ELSE

            DZ(NR,NC)=DZ(NR,NC)+PHASE*ANG(I)*(VZ(IDN,MC,NHFF(IR)))

         ENDIF

      ELSE

         IDN = ( ID(I) + 1 ) /2
         DZ(NR,NC) = DZ(NR,NC) + PHASE * ANG(I) * VZ(IDN,NHFF(IR),MC)
      END IF

 32   CONTINUE
 31   CONTINUE

      GO TO 20

C.................................................

   26 IF(I .GT. 3) GO TO 27

      DO 33 NC=1,NDI(IC)

         MC = NMINI(IC) + NC - 1

      IF(MC .NE. NHFF(IR)) GO TO 33

      DO 34 NR = 1, NDF(IR)

         MR = NMINF(IR) + NR - 1

      IF(2*INT(ID(I)/2).EQ.ID(I)) THEN

         IDN = ID(I) / 2

         IF(MODE.EQ.0) THEN 
            DZ(NR,NC) =DZ(NR,NC) + PHASE*ANG(I)*(-VZ(IDN,NHFI(IC),MR))
         ELSE
            DZ(NR,NC) =DZ(NR,NC) + PHASE*ANG(I)*(VZ(IDN,NHFI(IC),MR))
         ENDIF
      
      ELSE

         IDN = ( ID(I) + 1 ) / 2


         DZ(NR,NC) = DZ(NR,NC) + PHASE * ANG(I) * VZ( IDN,MR,NHFI(IC))

C         write(*,*) nr,nc, VZ( IDN,MR,NHFI(IC))

      END IF

 34   CONTINUE
 33   CONTINUE

      GO TO 20

C.................................................

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

C.................................................

      IF(IR.NE.1 .OR. IC.NE.1) GO TO 95

 95   CONTINUE 


      CALL TRNSMX( DMXZ,NHX, DZ, NS,KF, KI, NDF(IR), NDI(IC))

      KI = KI + NDI(IC)

 101  CONTINUE
      KF = KF + NDF(IR)

 100  CONTINUE

C.................................................


C      DSQ2 = SQRT(2.0D+00)
 
      IPI  = LHFI(1) + LLI(1)
      IP   = MOD(IPI, 2)
      IT   = 1

      IF(IP.EQ.1) GO TO 65

      DO 60 IC = 1, NCMXI

         IF(LHFI(IC) .NE. LLI(IC))   GO TO 60
         IF(NMINI(IC) .GT. NHFI(IC)) GO TO 60

         DO  NR = 1, NHX
            DMXZ(NR,IT) = DMXZ(NR,IT) / DSQ2
         ENDDO

 60      IT = IT + NDI(IC)

         GO TO 70

 65      CONTINUE

         DO 62 IR = 1, NCMXF

            IF(LHFF(IR).NE.LLF(IR))   GO TO 62
            IF(NMINF(IR).GT.NHFF(IR)) GO TO 62

            DO  NC = 1,NHX

               DMXZ(IT,NC) = DMXZ(IT,NC)/DSQ2
            ENDDO

 62            IT = IT + NDF(IR)

 70            CONTINUE


C..... initial states [E,C(E)]

      DO  NSC = 1, NSI

         READ(unit=nwf2e_i)   ENI(NSC)
         READ(unit=nwf2e_i) ( CIIN(NSC,K), K = 1, NHMXI)

      ENDDO


C..... final states [E,C(E)]

      DO  NSR = 1, NSF

         READ(unit=nwf2e_f)  ENF(NSR)
         READ(unit=nwf2e_f) (CFIN(NSR,K), K=1,NHMXF)

      ENDDO

C.

      CLOSE(nwf2e_i)
      CLOSE(nwf2e_f)

C     *
C     *             STORE AS 
C     *    LIN, LFIN, NIN, NFIN
C     *    E_IN (I)  ENERGIES OF INITIAL STATES
C     *    E_FIN(I)  ENERGIES OF FINAL   STATES
C     *    DME( L_IN  ----> L_FIN)
C     *


      WRITE(nd2e)   MODE
      WRITE(nd2e)   LOI, LOF, NSI, NSF
      WRITE(nd2e)   ( ENI(K), K = 1, NSI )
      WRITE(nd2e)   ( ENF(K), K = 1, NSF )

C.......................................................

 55   FORMAT(4I5)

         IF(NINI.EQ.0) THEN
            WRITE(16, 55) MODE
            WRITE(16, 55) LOI,LOF,NSI,NSF
            WRITE(16, 5) (ENI(NSC), NSC = 1, NSI )
            WRITE(16, 5) (ENF(NSR), NSR = 1, NSF )
         ENDIF
C........................................................




         DO  NSC = 1, NSI


            DO  KC = 1, NHMXI
               CI(KC) = CIIN(NSC,KC)
            ENDDO


C     #    MAKE THE MATRIX-VECTOR MULTIPLICATION


            CALL DGEMV('N',NHMXF,NHMXI,1.0D0,DMXZ,NHX,CI,1,0.0D0,TL1,1)

            DO NSR = 1, NSF


               DO KR  = 1, NHMXF
                  CF(KR) = CFIN(NSR, KR)                  
               ENDDO
            
C     #    MAKE THE VECTOR-VECTOR MULTIPLICATION

               DMX(NSR)  = DDOT( NHMXF, CF, 1, TL1, 1 )

            ENDDO

C     #    STORE IN BINARY FILE THE 2-E DIPOLE MATRIX ELEMENTS

            WRITE(nd2e)    (DMX(NSR), NSR = 1, NSF)

            IF((NSC.EQ.1).AND.(NSR.EQ.1)) THEN
               DO I = 1, NSR
                 WRITE(16, '(2e20.10,1X,i5)') ENF(NSR), DMX(NSR), NSR
               ENDDO
            ENDIF
            
      
         IF(NSC.EQ.NINI) THEN

            WRITE(*,*) '           DE(N - ', nini,' )           DMX2E '
            DO NSR = 1, NSF-1
               WRITE(*, '(4E20.10)') ABS(ENF(NSR) - ENI(NINI)),
     &                DMX(NSR)**2/abs(enf(nsr-1)-enf(nsr+1)),
     &                               DMX(NSR)**2,
     &                               DMX(NSR)
            ENDDO
            
         ENDIF

      ENDDO

       CLOSE(16)
       CLOSE(nd2e) 

C*    FORMAT STATEMENTS
C..................................................

 5     FORMAT(2X, 1P8E15.7 )
C..................................................

       END
C########################################################################
C#EOF
