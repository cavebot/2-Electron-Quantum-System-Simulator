C#    PROGRAM TO CALCULATE THE H-MATRIX FOR CA-LIKE ATOMS
C#    BSP HF WF ARE USED FOR ALL CONFIGURATIONS (NORMAL PARITY)



      PROGRAM R12
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.2e.inc"
C      PARAMETER(NP = 2000, NS = 802, NHX = 2500, NCS = 30, NL=15 )
      DIMENSION PR1(NP),PC1(NP),PR(NS,NP),PC(NS,NP),YK(NP)
      DIMENSION VX(NS,NS),EHF(NL,NS)
      DIMENSION NHF(NCS),LHF(NCS),LL(NCS),NLLMIN(NCS)
      DIMENSION NLLMAX(NCS),NOLL(NCS),IS(NCS)
      DIMENSION nd(ncs) 
C     
      COMMON/HX/HMX(NS,NHX)
      COMMON/RARR/RPW(NL,NP)
      COMMON/RRO/RHO(6),alpha,beta
      COMMON/RH/H,R(NP),DR(NP)
      CHARACTER*100 ARGV
C.................................................



C      Q2 = DSQRT(2.0D+00)


C calculate for symmetry L


      CALL GETARG(1, ARGV)
      READ(ARGV,*) INPL

      WRITE(*,*) '# v2eb::    partial wave  L =', INPL
      WRITE(*,*) '# v2eb::    parameters:'
      WRITE(*,*) '# v2eb::          ns = ', ns
      WRITE(*,*) '# v2eb::          np = ', np
      WRITE(*,*) '# v2eb::          nl = ', nl
      WRITE(*,*) '# v2eb::          lf = ', lf
      WRITE(*,*) '# v2eb::         nhx = ', nhx
      WRITE(*,*) '# v2eb::         ncs = ', ncs




      
C#    alpha :: dielectronic polarization potential: Vd
C#    beta  :: mass polarization potential        : Vm

      OPEN(9,FILE='inp/r12.inp')
      READ(9,*) alpha, beta
      READ(9,*) (rho(i), i = 1, 5)
      CLOSE(9)

      OPEN(16,FILE='out/r12.out')

C#     Read 1-e angular momenta and energies (wf2e-en.inp)


      OPEN(15,FILE='dat/en1e.dat')

      READ(15,3) LLMIN,LLMAX

      DO  LP = LLMIN, LLMAX

         DO  IE = 1, NS
            EHF(LP,IE) = 0.0D+00
         ENDDO 

         READ(15,3) NMX, NVALENCE
         DO  IE = 1, NMX
            READ(15,'(E25.14)') EHF(LP,IE)
         ENDDO      

      ENDDO

      CLOSE(15)

      write(*,*) ' read en1e.dat'

c      CALL read_v(l, en, "en1e-")
c


C*    prepare for the grid. (open 'grid.inp') 

      CALL rin(r, dr, h, no, idr)
      close(10)

      write(*,*) ' read inp/r12.inp, grid constructed.'


      DO   KK = 1, NL
        DO  J = 1, NO
          RPW(KK,J) = R(J)**KK             
        ENDDO
      ENDDO



C#   read configurations from inp/cfg-L.inp file 


      NCFG = 15
      CALL CFGFILE(NCFG, INPL) 


      READ(NCFG,*) NSYMX

C      DO 200 NSY = 1, NSYMX


         CALL CFIN(NCFG, NHF, LHF, LL, NLLMIN, NLLMAX, NOLL,
     1        IS, LO, ISL, NCSMX, NCMX)

         itry = 1

C         READ( NCFG, 3 ) itry

         IF(ncsmx.gt.ncs) THEN
            WRITE(*,*)' Error:: ncsmx or itry exceeds maximum!Exiting..'
            STOP
         ENDIF
 
C         DO  it = 1, itry
C            READ(NCFG,7) ( nd(k), k = 1, ncsmx ),iouti,ioutf
C       ENDDO

C#
      if(isl.eq.1) then
         write(*,*) ' SINGLET SYMMETRY'
      else if(isl.eq.3)   then
         write(*,*) ' TRIPLET SYMMETRY'
      else
         write(*,*) ' Allowed values of isl in r12.inp are 1 (Singlet),
     1        or 3 (triplet) '
      endif

      if(inpl.ne.lo) then
         write(*,*) ' Input L and L from cfg file differ. Exiting'
         stop
      endif


C     # 


      WRITE(16,6) LO, NCMX


      call r12file(3, LO)

 
      write(*,*)' # configuration series,                 ncs = ',ncsmx
      write(*,*)' # two-electron  uncorrelated states,    n2e = ',ncmx


      KR = 1
      DO IR = 1,  NCSMX


      LLR   = LL(IR)  + 1
      LLHFR = LHF(IR) + 1
      NHFIR = NHF(IR)


C#
C#     n1,    l1,  l2, n2_min, n2_max
C#     nhfir, lhf, ll, nllmin, nllmax
C#

C#
C#  store P_{n1,l1}(r_j)    r_j : grid points
C#


      CALL WFHFIN( PR1, LHF(IR), NO, NHFIR)

C#
C#  store P_{n2,l2}(r_j)    n2_min < n_2 < n2_max, r_j : grid points
C#



      CALL WFNLIN( PR, YK, LL(IR), NO, NLLMIN(IR), NLLMAX(IR))

      KC = KR

      DO  I = 1, NS
         DO  K = 1, NHX
            hmx(i,k) = 0.0D+00
         ENDDO
      ENDDO


      DO  IC = IR, NCSMX

         LLC   = LL(IC) + 1
         LLHFC = LHF(IC) + 1

C*


         CALL WFHFIN( PC1, LHF(IC), NO, NHF(IC))
         CALL WFNLIN( PC, YK, LL(IC), NO, NLLMIN(IC), NLLMAX(IC))

C#   vx(ns,ns)
         call submx(PR1, PR, PC1, PC,R,DR,VX,ISL,NHF,NLLMIN,NOLL,LHF,LL,
     &        YK, EHF, IR, IC, LO, NO, H, IDR)

C# vx(ns,ns) ---> hmx(ns,nhx)

         CALL TRNSMX(HMX, VX, 1, KC, NOLL(IR), NOLL(IC))
         KC = KC + NOLL(IC)
       enddo

c#     upper part     noll = n2_max - n2_min + 1
       write(*,*) ir, noll(ir)
       KR = KR - 1
       DO  NR = 1, NOLL(IR) 

         KR = KR + 1

         WRITE(3) ( HMX( NR, NC), NC = KR, NCMX )

       enddo

       KR = KR + 1

      enddo

      CLOSE(3)
      CLOSE(4)

C  200 CONTINUE

      STOP
C================================================================
C#
C#                    format statements
C#

 1    FORMAT(4D15.7)
 2    FORMAT(2X,1P9D14.6)
 3    FORMAT(6I5)
 4    FORMAT(1PD16.8)
 5    FORMAT(/2X,'ROW -- NL & LL = ',3I3,',  COLUMN -- NL & LL = ',3I3)
 6    FORMAT(2X,'L & S =',2I4/)
 7    FORMAT(10i4)
 566  FORMAT(5A16)

C#
C#
C=================================================================
      END

C*EOF
