C* 
C*             CPC/ SUBMITTED VERSION    1.0 / F77
C*             BSPCI2E PACKAGE : 
C*
C*             PROGRAM NUMBER 1       D1E
C*     
C*
C*LAAN
C*######################################################################
      PROGRAM D1E

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

      INCLUDE "prm.1e.inc"

      DIMENSION AMAT(NS, NS),BMAT(NS,NS), ER(NS)
      DIMENSION AD(NS, NS), AE(NS,NS)
      DIMENSION POD(NL,NP),EN(NL,NS)
      DIMENSION T(NS+NK)
      DIMENSION R(NP),P(NL,NP),VD(NP),FF(NP),DR(NP)
      DIMENSION RS(NL),YK(NP),ZK(NP),WY(NP),WZ(NP)
      DIMENSION BR(NP),BC(NP), DTM(NL),N2(NL)
      DIMENSION NCORE(3)
      INTEGER   NO_S, NO_P, NO_D, NSP, NTC
      INTEGER   NINP, NOUT, NBIN
      DIMENSION ALP1(NL), R01(NL) 

C...................................................................

        WRITE(*,*) 'I GOT THIS FAR'
        WRITE(*,*) 'I GOT THIS EVEN FURTHER'

C...................................................................


      NINP  = 14
      NOUT  = 16
      NBIN  = 1


C...................................................................

      OPEN(NINP, FILE='inp/d1e.inp', STATUS='OLD')

      READ(NINP, *) ZNUC, RMIN, RMAX
      READ(NINP, *) NPO, REDMASS
      READ(NINP, *) ( NCORE(I),  I = 1, 3 )
      READ(NINP, *) ISL, NW
      READ(NINP, *) LMIN, LCORE, LMAX, ITRMX
      READ(NINP, *) ( ALP1(I + LMIN - 1), I = 1, NW )
      READ(NINP, *) ( R01(I  + LMIN - 1), I = 1, NW )
      READ(NINP, *) ( ZK(I),  I = 1, NCORE(1) + NCORE(2) + NCORE(3) )
      READ(NINP, *) ( N2(I),  I = 1, NW ), K
      READ(NINP, *) ( RS(I),  I = 1, NW )
      READ(NINP, *)  CRT, DI, DF, ID

      CLOSE(NINP)

C...................................................................
 
      OPEN(NOUT, FILE='out/d1e.out')

      WRITE(*,*)'#######################################'
      WRITE(*,*)' '
      WRITE(*, *)'NUMBER OF          NBSP = ', N2(1)
      WRITE(*, *)'ORDER BSPL.           K = ', K
      WRITE(*, *)'BOX  RADIUS           R = ', RMAX
      WRITE(*, *)'FIRST POINT         FKN = ', RS(1)    
      WRITE(*, *)'#     OF S,P,D ORBITALS = ',NCORE(1),NCORE(2),NCORE(3)

      IF((ISL.NE.1).AND.(ITRMX.NE.1)) THEN
         WRITE(*,*)'HARTREE-FOCK     LCORE = ', LCORE
      ELSE 
         WRITE(*,*)'HYD - LIKE FUNCTIONS OR EXISTING CORE COEFF USED'
      ENDIF

      IF(REDMASS.EQ.1) THEN

         WRITE(*,*)'HYD - LIKE FUNCTIONS ARE CACLULATED'

      ELSE IF(REDMASS.EQ.2) THEN

         WRITE(*,*)'PS - LIKE FUNCTIONS ARE CALCULATED'
      ELSE 

         WRITE(*,*)'NOT ACCEPTABLE VALUE FOR PARAM REDMASS. EXITING...'

         STOP
      ENDIF


      IF( (ALP1(LMIN).EQ.0.D0) ) THEN
         WRITE(*,*)' NO CORE -POLARIZATION '
      ELSE
         WRITE(*,*)' CORE - POLARIZATION IS INCLUDED'
      ENDIF

      WRITE(*,*)'#######################################'


C............................................................................
      

      NO_S = NCORE(1) 
      NO_P = NCORE(2)
      NO_D = NCORE(3) 
            

      NSP = NO_S + NO_P
      NTC =  NSP + NO_D


C............................................................................


C#      
C#     SOME CHECKS FIRST
C#      


      DO  I = LMIN, LMAX - LMIN + 1

         IF(N2(I).GT.NS.OR.N2(I).LT.(K+2)) THEN

            WRITE(*, *) ' NK + 2 < N2  < K+2'
            WRITE(*, *) '  N2 = ', N2(I)
            WRITE(*, *) '  NS = ',  NS
            WRITE(*, *) ' K+2 = ',  K + 2

            GOTO 999

         ENDIF
      ENDDO


      IF(K.GT.NK.OR.K.LT.2) THEN

         WRITE(*, *)' PROBLEM WITH K !! K=',K
         GOTO 999

      ENDIF

      IF(DI.GT.1.0D+00 .OR. DF.GT.1.0D+00) THEN

         WRITE(*, *)' DI OR DF GREATER THAN 1  '
         GOTO 999

      ENDIF

C...............................

C#
C#  GENERATE KNOT POINTS. DATA PROVIDED FROM 'RIN.DIN'
C#

      CALL RIN(R, DR, H, NPO, IDR)


      OPEN(26, FILE ='dat/knot.dat')

      DO II = 1, NPO

         WRITE(26,*) II,'  ', R(II)      

      ENDDO
 
      CLOSE(26)

C...............................

C#
C#       PREPARE H-LIKE CORE WF'S
C#

      IF (ISL.EQ.1) THEN

         DO  I = 1, NO_S
            
            IF ( R(1).NE.0.0D+00) THEN 

               P(I, 1 ) = HWF( I, 0, ZK(I), R(1), ER(I) )
            ELSE

               P(I, 1 ) = 0.0D+00
            ENDIF
                        
            DO J = 2, NPO

               P(I, J)  = HWF( I, 0, ZK(I), R(J), ER(I) )

            ENDDO
         ENDDO

c         IF(NO_S.EQ.NTC) GO TO 45

         IF(NO_S.EQ.NTC) then
c DONOTHING            
         ELSE
           DO  I = NO_S + 1, NSP

              IF (R(1).NE.0.0D+00) THEN 
                 
                 P(I, 1) = HWF( I, 1, ZK(I), R(1), ER(I) )

              ELSE
                    
                 P(I, 1) = 0.0D+00

              ENDIF
                 
              DO  J = 2, NPO

                 P(I, J) = HWF( I, 1, ZK(I), R(J), ER(I) )
                 
               ENDDO

             ENDDO
          ENDIF

c           IF((NSP).EQ.NTC) GO TO 45

        IF((NSP).EQ.NTC) THEN
c DO NOTHING           
        ELSE

           DO  I = NSP + 1, NTC

              IF (R(1).NE.0.0D+00) THEN  
                 
                 P(I,1) = HWF(I, 2, ZK(I), R(1), ER(I) )
              ELSE
                    
                 P(I,1) = 0.0D+00
              ENDIF

              DO  J = 2, NPO

                 P(I, J) = HWF(I, 2, ZK(I), R(J), ER(I))

              ENDDO
              
           ENDDO
            
        ENDIF

      ENDIF
c 45   ENDIF


C         DO I = 1, NTC
C              DO  J = 1, NOP

C                 WRITE(28+I, *) R(J), P(I, J) 

C              ENDDO
C           ENDDO

C#
C#  READ EXISTING COEF. & GENERATE CORE WF'S.
C#
  
C#       
C# this is not a HF calculation
C#
         IF(ISL.NE.1) THEN

         NI = 1
         NF = NO_S

         DO 15 L = 1, LCORE

            LANG = L - 1

C........ OPEN CORE WF FILES

            CALL COREFILE(NBIN, LANG) 

C            CALL D1EFILE(1, LANG) 


            READ(1) RS1, RMX, N1, K1, LIN


            CALL MKGRID( N1, K1, RMX, RS1, T )


            DO NPP = 1, N1 - 2

               READ(NBIN) EN(L, NPP), ( AMAT(J, NPP),J = 1, N1 - 2)
            ENDDO



            WRITE(NOUT, * ) ' L = ', LANG,'  INITIALLY INPUT ENERGIES ='
            WRITE(NOUT, 8 ) ( 0.5 * EN(L,I),I = N1 - 2,1,-1)


C#    SUBROUTINE COREWF() TAKES CMAT, R, DR, N, K, NI, NF, NS, NOP, H,
C#    RETURNS CORE WF'S P(I,J)(( I = NI, NF ),(J = 1, NOP) ).


            CALL COREWF(AMAT,NS,P,NL,FF,R,DR,T,NPO,N1,K1,NI,NF,H )

C            CALL COREWF( AMAT, P, FF, R, DR, T, N1, K1, NI, NF, NOP, H )

            NI = NF + 1

            IF(L.EQ.1) NF = NF + NO_P
            IF(L.EQ.2) NF = NF + NO_D

            CLOSE(1)

 15      CONTINUE
      ENDIF

C#       THIS LOOP IMPROVES COEF. FOR CORE ORBITS FOR ITRMX TIMES.
C#       D1 = DI + DFLOAT(ITRY) / DFLOAT(ID) * (DF - DI)

C#
C#     L > LCORE   ---> DO NOT CALCULATE HF ORBITALS
C#       

      IF (LMIN.GT.LCORE) GOTO 201


      DO 200 ITRY = 1, ITRMX

         D1  = DI
         IF(ITRY.GT.ID) D1 = DF

         WRITE(*,*)'  #####################################'
         WRITE(*,*)'  ITERATION NUMBER : ', ITRY
         WRITE(NOUT,*)'  #####################################'
         WRITE(NOUT,*)'  ITERATION NUMBER : ', ITRY

         NI = 1
         NF = NO_S
         DTMA = 0.0D+00

         DO 100  LL = 1, LCORE

C         WRITE(*,*)'....................................... L = ', LL-1

            DR0 = RS(LL)
            N   = N2(LL)

            LANG = LL - 1

            CALL D1EFILE(NBIN, LANG)

C            WRITE(NOUT,3) ZNUC,H,RMIN,DR0,RMAX,N,K,NPO,IDR,LANG,D1
            
C#
C#      GETEIG(...) TAKES LANG,N,K,DR0,RMAX,NOP. NCORE,N,M,R,DR,P,IDR,
C#                        
C#      RETURNS COEF. CMAT & ENERGIES ER. 
C#      OTHER PARAMETERS ONLY PROVIDE STORAGE SPACE 
C#


            CALL GETEIG(LANG, ALP1, R01, ZNUC, REDMASS, N, K, DR0, RMAX,
     &          NPO, NCORE, H, T, R, DR, P, NL,VD,FF,YK,ZK,WY,WZ,BR,BC, 
     &           AMAT, BMAT, AD, AE, NS, ER, IDR)


               
               DO  I = NI, NF
                  DO  J = 1, NPO                     
                     POD(I,J) = P(I,J)
                  ENDDO
               ENDDO


               CALL COREWF(AMAT,NS,P,NL,FF,R,DR,T,NPO,N,K,NI,NF,H)

C               CALL COREWF( AMAT, P, FF, R, DR, T, N, K, NI, NF, NOP, H)
               

C     *          ORN = |< WF_N | WF_N > | ^2   ( SHOULD BE == 1)

               DO 22 I = NI, NF

                  FF(1) = 0.0D+00
                  DO  J = 2, NPO

                     P(I,J) = D1 * P(I,J) + (1.D0 - D1) * POD(I,J)
                     FF(J)  = P(I,J) * P(I,J) * R(J) / DR(J)

                  ENDDO

                  ORN = RINT(FF, 1, NPO, 14, H )

                  SQ  = DSQRT(ORN)


                     IF(ITRMX.NE.1) THEN

                        CALL WRITENORM(I, LL, ORN) 
                     
                     ENDIF

                     DO 22 J = 1, NPO

 22                   P(I,J) = P(I,J) / SQ 

                      
C***  FROM HERE TO NEXT COMMENT LINE ARE AUXILIARY. 

                      IF( (ITRY.EQ.1) .AND. (N1.NE.N) ) GOTO 52
         
                      DLTEM = 0.0D+00

                      DO 51 IE = 1, N - 2
                              
                         DLTE = ABS( ER(IE) - EN(LL,IE) )

                         IF(DLTEM.LT.DLTE) THEN

                            KIE   = ( N - 2 ) - IE + 1
                            DELT  = ER(IE) - EN(LL, IE)
                            DLTEM = DLTE

                         ENDIF

 51                   CONTINUE

 52                   DO 53 IE = 1, N - 2

 53                      EN(LL,IE) = ER(IE)
                              
                         DTM(LL) = 0.0D+00
                         
                         DO I = NI, NF
                            DO J = 1, NPO

                               DLT = ABS(P(I,J) - POD(I,J))

                               IF(DTM(LL).LT.DLT) DTM(LL) = DLT

                            ENDDO
                         ENDDO

C**   ...........LOG
C                         WRITE(NOUT, *) KIE, DELT, DTM(LL)
                         WRITE(NOUT, * ) ' ENERGIES  '
                         WRITE(NOUT, 8 ) (0.5 * ER(J),J = N-2, 1, -1)

                         IF(DTMA.LT.DTM(LL))  DTMA = DTM(LL)

                         IF((DTM(LL).LT.CRT).OR.(ITRY.EQ.ITRMX)) THEN
  
                            
                            WRITE(1) DR0, RMAX, N, K, LANG

C***  WRITE IN BINARIES IN INCREASING EIGENVALUE ORDER

                            DO I = 1, N - 2

                               WRITE(NBIN) ER(I),(AMAT(J,I), J = 1, N-2)              
                            ENDDO
                                    
                         ENDIF

                         IF(DTMA.LT.CRT.AND.LL.EQ.LCORE) GOTO 201
                         NI = NF + 1
                        
                         IF(LANG.EQ.0) NF = NF + NO_P
                         IF(LANG.EQ.1) NF = NF + NO_D

                         CLOSE(NBIN)

 100               CONTINUE
 200    CONTINUE



C#
C#       IF  LCORE > LMAX     -->   LEAVE
C#

 201    IF(LCORE.GE.LMAX) GOTO 999

        IF(LMIN.GT.LCORE) LMN = LMIN
        IF(LMIN.LE.LCORE) LMN = LCORE + 1
                           
C**   THIS LOOP USES COEF FOR CORE TO CALCULATE COEF. FOR NON-CORE ORBITS.
C**    LL > LCORE

        DO 205 LL = LMN, LMAX

           LANG = LL - 1
           DR0  = RS( LL - LMN + 1 )
           N    = N2( LL - LMN + 1 )

           CALL D1EFILE( NBIN, LANG)  

           WRITE(*,*) ' LL ', LL
           WRITE(*,*) ' LMIN = ', LMN
           WRITE(*,*) ' LMAX = ', LMAX
           WRITE(*,*) ' R(',LANG,') = ', DR0
           WRITE(*,*) ' N(',LANG,') = ', N
           WRITE(*,*) '**********************'
          
           WRITE(NOUT,*) ZNUC,H,RMIN,DR0,RMAX,N,K,NPO,IDR,LANG

           CALL GETEIG(LANG, ALP1, R01, ZNUC, REDMASS, N, K, DR0, RMAX,
     1          NPO, NCORE, H, T, R, DR, P, NL, VD,FF,YK,ZK,WY,WZ,BR,BC,
     2          AMAT, BMAT, AD, AE, NS, ER, IDR)

             WRITE(NOUT,*)' EN = '
             WRITE(NOUT,8) ( 0.5 * ER(J),J = N - 2, 1, -1 )

C     #  WRITE IN BINARIES IN INCREASING EIGENVALUE ORDER

           WRITE(NBIN) DR0, RMAX, N, K, LANG
           
           DO I = 1, N - 2

              WRITE(NBIN) ER(I),( AMAT(J,I), J = 1, N - 2 )

           ENDDO

           CLOSE(NBIN)
        
 205    CONTINUE

 999    STOP

C...................................................................

 3      FORMAT(/2X,' FRACTION OF NEW WF USED: D1 = ',E9.2) 
 5      FORMAT(10I5)
 7      FORMAT(5E13.4)
 8      FORMAT(5D20.12)

 566  FORMAT(5A16)
C...................................................................

        END
C############################################################

      SUBROUTINE WRITENORM(I, LL, ORN) 
      INTEGER LL, I
      DOUBLE PRECISION ORN
C.............................

      IF(LL.EQ.1.AND.I.EQ.1) THEN 
         WRITE(*,*) ' <P_1S | P_1S >     = ', ORN
      ELSE IF (LL.EQ.1.AND.I.EQ.2) THEN 
         WRITE(*,*) ' <P_2S | P_2S >     = ', ORN
      ELSE IF (LL.EQ.1.AND.I.EQ.3) THEN 
         WRITE(*,*) ' <P_3S | P_3S >     = ', ORN
      ELSE IF (LL.EQ.2.AND.I.EQ.3) THEN 
         WRITE(*,*) ' <P_2P | P_2P >     = ',  ORN
      ELSE IF (LL.EQ.2.AND.I.EQ.4) THEN 
         WRITE(*,*) ' <P_2P | P_2P >     = ', ORN
      ELSE IF (LL.EQ.2.AND.I.EQ.5) THEN 
         WRITE(*,*) ' <P_3P | P_3P >     = ', ORN
      ELSE
         WRITE(*,*) ' ORBITAL  : ', I
         WRITE(*,*) ' NORM     = ', ORN
      ENDIF
      RETURN
      END
C#############################################################
C#  EOF 












