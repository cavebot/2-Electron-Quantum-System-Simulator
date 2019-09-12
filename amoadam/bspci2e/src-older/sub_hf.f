C* 
C*             CPC/ SUBMITTED VERSION    1.0 / F77
C*             BSPCI2E PACKAGE : 
C*
C*             PROGRAM NUMBER         SUBD1E
C*     
C*             SUBROUTINES            FUNCTIONS
C*
C*             GETEIG                 Q, QQ                                    
C*             COREWF                 HWF 
C*             VHFBSP 
C*             RKPREP
C*             DIRECTB
C*             XCHB
C*             YKFCTB
C*             YKEXPB
C*             YKSINB
C*             SOLVEMAT
C*             SETMAT
C*
C*LAAN    
C*######################################################################
      SUBROUTINE GETEIG(LANG, ALP1, R01, ZNUC, REDMASS, N, K, DR0, RMAX, 
     1     NPO,NCORE, H, T, R, DR, P, NL, VD, FF, YK, ZK, WY, WZ,BR,BC, 
     2     AMAT, BMAT, AD, AE, NS, ER, IDR)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION AMAT(NS, *), BMAT(NS, *)
      DIMENSION   AE(NS, *),   AD(NS, *)
      DIMENSION    NCORE(3),    P(NL, *)
      DIMENSION       ER(*) 
      DIMENSION     ALP1(*),  R01(*)
      DIMENSION  T(*), R(*), DR(*)
      DIMENSION VD(*), FF(*), YK(*), ZK(*)
      DIMENSION WY(*), WZ(*), BR(*), BC(*) 

C......... 

C#
C#       CREATE THE GRID KNOT SEQUENCE
C#
C#       CALCULATE THE DIRECT AND EXCHANGE POTENTIAL FOR THE FROZEN CORE 
C#
C#       SET THE ( N-2 X N - 2 ) Schrondinger Equation 
C#
C#          AMAT --> B-splines HAMILTONIAN MATRIX
C#          BMAT --> B-splines OVERLAP MATRIX
C#
C#         SOLVE THE  SCHRODINGER MATRIX EQNS ON B SPLINES 
C#                     A * X = E * B * X 


!      CALL MKGRID( N, K, RMAX, DR0, T)

      CALL VHFBSP( T, R, DR, P, NL, VD, FF, YK, ZK, WY, WZ, BR,
     1             BC, AD, AE, NS, H, N, K, NPO, LANG, NCORE, IDR)

!      CALL SETMAT( AMAT, BMAT, AD, AE, NS, N, K, T, 
!     1             ALP1, R01, ZNUC, LANG, REDMASS)

!      CALL SOLVEMAT( AMAT, BMAT, NS, ER, N-2)


      RETURN
      END
C#######################################################################
      SUBROUTINE COREWF(C, NS, P, NL, F, R, DR, T, NPO, N, K, NI, NF, H)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION C(NS,*), P(NL,*)
      DIMENSION F(*), R(*), DR(*), T(*)
      DIMENSION COE(NS), WORK(3*K) 
C..............................................................

      IN = N - 1
      COE(1) = 0.0D+00
      COE(N) = 0.0D+00

      DO  NPP = NI, NF

         IN = IN - 1

         DO  I = 1, N - 2

            COE( I + 1 ) = C(I, IN)

         ENDDO

C     #         D'BOOR ROUTINE
        
         P(NPP, 1 ) = 0.0D+00
         DO  J = 2, NPO

            INV = 1 

            P(NPP, J) = DBVALU(T, COE, N, K, 0, R(J), INV, WORK)

         ENDDO

C..... ?
        IF(P(NPP, 6).GT.0.0D+00) GOTO 14

        DO  J = 1, NPO
           P(NPP,J) = - P(NPP,J)

        ENDDO

 14     CONTINUE

        F(1) = 0.0D+00

        DO  J = 2, NPO

           F(J) = P(NPP,J) * P(NPP,J) * R(J) / DR(J)
        ENDDO

        OTH = RINT( F, 1, NPO, 14, H)

        SQ2 = DSQRT(OTH)

        IF(SQ2.EQ.0.0D+00) THEN 

           WRITE(*,*) ' COREWF : PROBLEM FOR SQ2 = ', SQ2 
           STOP           
        ENDIF

        DO J = 1, NPO
           
           P(NPP,J) = P(NPP,J) / SQ2
        ENDDO

      ENDDO

      RETURN
      END
C#######################################################################
      SUBROUTINE VHFBSP(T,R,DR,P,NL,VD,FF,YK,ZK,WY,WZ,BR,BC,AD,AE,NS,
     1                  H, N, K, NPO, LANG, NCORE, IDR)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION AD(NS, *), AE(NS, *), BF(NS-2, NPO)
      DIMENSION  P(NL, *) 
      DIMENSION  C(NS),  WORK(3*K)
      DIMENSION T(*),  R(*),  DR(*), VD(*), FF(*), NCORE(*)
      DIMENSION YK(*),  ZK(*), WY(*), WZ(*), BR(*), BC(*)

C......................................


C    4 FORMAT(1X,1P9E14.6)

C........................

      CALL DIRECTB( P, NL, VD, YK, ZK, WY, WZ, DR, R, FF, H, NPO, NCORE,
     1 IDR)

C.......................


      DO 50 NC = 1, N - 2


         DO  NN = 1, N

            C(NN)   = 0.0D+00
         ENDDO


            C( NC + 1 ) = 1.0D+00
            DO  J = 1, NPO

C#      SLATEC  D'BOOR

               INV  = 1
               BF(NC, J) = DBVALU( T, C, N, K, 0, R(J), INV, WORK )

            ENDDO

 50         CONTINUE

C......................

            DO 100 KR = 1, N - 2

               DO  J = 1, NPO

                  BR(J) = BF( KR, J)
               ENDDO

C#
C#  EXCHANGE POTENTIAL KERNEL
C#
               CALL XCHB( P, NL, BR, YK, ZK, WY, WZ, DR, R, FF, H, LANG,
     1             NPO, NCORE, IDR )


               
 10            DO 101 KC = 1, KR


C#          CALCULATE AND STORE DIRECT POTENTIAL IN AD  
                  DO  J = 1, NPO

                     BC(J) = BF(KC, J )

                  ENDDO

                  IF(DR(1).GT.1.0D-15) THEN

                     FF(1) = - BR(1) * BC(1) * VD(1) / DR(1)
                  ENDIF

                  DO  J = 2, NPO

                     FF(J) = - BR(J) * BC(J) * VD(J) / DR(J)
                  ENDDO



                     AD(KR, KC) = RINT( FF, 1, NPO, 14, H)
 
                     AD(KC, KR) = AD(KR, KC)


C#          CALCULATE AND  TORE EXCHANGE POTENTIAL IN AD

                     DO J = 2, NPO

                        FF(J) = YK(J) * BC(J) / DR(J)
                        
                     ENDDO

                     AE(KR,KC) = RINT( FF, 1, NPO, 14, H)

                     AE(KC,KR) = AE(KR, KC)

                     
 101              CONTINUE


                  
 100  CONTINUE
      
      RETURN
      END
C#######################################################################
      SUBROUTINE RKPREP(R, DRKR, RKK,K, NPO)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION DRKR(*),RKK(*),R(*)
C.........................................

      KR = K 
      KK = K + 1

      IF(KR.EQ.0) GO TO 30

      DRKR(1) = 0.0D+00

      IF(R(1).NE.0.0D0) DRKR(1) = 1.0D+00 / R(1)**KR

      DO  J = 2, NPO

         DRKR(J) = 1.0D+00 / R(J)**KR
      ENDDO

 30   CONTINUE

      DO J = 1, NPO

         RKK(J)   =  R(J)**KK
       ENDDO

      RETURN

      END
C#######################################################################
      SUBROUTINE DIRECTB( P, NL, FT, YK, ZK, WY, WZ, DR, R, F, H, 
     1                   NPO, NCORE, IDR)


      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION P(NL, *), NCORE(*)
      DIMENSION YK(*), ZK(*), WY(*), WZ(*)
      DIMENSION FT(*), F(*),  DR(*),  R(*)
      DIMENSION DRKR(NPO), RKK(NPO)

C.....................................

      NO_S = NCORE(1)
      NO_P = NCORE(2)
      NO_D = NCORE(3) 

C.....................................
      
C*    (1S^2)(2S^2)(2P^6)(3S^2)(3P^6)(3D^10)

      QS = 2.0D+00
      QP = 6.0D+00
      QD = 10.0D+00

C#
C#   CALCULATE   R^(K+1) , 1/R^K AT GRID POINTS
C#  
C#                 K = 0
C#


      CALL RKPREP(R, DRKR, RKK, 0, NPO)


C**    NO CORE  HE, H- 


      DO  J  = 1, NPO

         F(J)  = P(1,J)
         FT(J) = 0.0D+00

      ENDDO

      IF(NO_S.EQ.0) RETURN


C*    ( 1S^2 )  CORE 
C*      S = 0   CORE STATES

      CALL YKFCTB(F,F,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,0,NPO,IDR)

      DO 11 J  = 1, NPO
 11      FT(J) = FT(J) + 2.D0 * QS * YK(J)

      IF(NO_S.EQ.1) RETURN
C********************************************************
C*
C*    ( 1S^2 2S^2 2P^6)  CORE 
C*
C********************************************************

C*
C*      S = 0   CORE STATES
C*

      DO 12 J = 1, NPO
 12      F(J) = P(2,J)

      CALL YKFCTB(F,F,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,0,NPO,IDR)

      DO 13 J = 1, NPO
         F(J)  = P( NO_S + 1, J)
 13      FT(J) = FT(J) + 2.D0 * QS * YK(J)


C*
C*      P = 1   CORE STATES
C*

      CALL YKFCTB(F,F,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,0,NPO,IDR)

      DO 14 J = 1, NPO
 14      FT(J)= FT(J) + 2.D0 * QP *YK(J)

      IF(NO_S.EQ.2) RETURN

C********************************************************
C*
C*    ( 1S^2 2S^2 2P^6 3S^2 3P^6 3D^10)  CORE 
C*
C********************************************************
      DO 15 J = 1, NPO
 15      F(J) = P(3,J)

C*
C*      S = 0   CORE STATES
C*

      CALL YKFCTB(F,F,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,0,NPO,IDR)

      DO 16 J = 1, NPO
         F(J) =  P(NO_S + 2, J )
 16      FT(J) = FT(J) + 2.D0 * QS *YK(J)


C*
C*      P = 1   CORE STATES
C*

      CALL YKFCTB(F,F,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,0,NPO,IDR)

      DO 17 J = 1, NPO
 17      FT(J) = FT(J) + 2.D0 * QP * YK(J)

      IF(NO_S.EQ.3) RETURN

C*
C*      S = 0   CORE STATES
C*


      DO 18 J = 1, NPO
 18      F(J) = P(4,J)

      CALL YKFCTB(F,F,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,0,NPO,IDR)

      DO 19 J = 1, NPO
         F(J)  = P(NO_S+3,J)
 19      FT(J) = FT(J) + 2.D0* QS * YK(J)

C*
C*      P = 1   CORE STATES
C*

      CALL YKFCTB(F,F,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,0,NPO,IDR)


      DO 20 J=1,NPO
         F(J) = P(NO_S+NO_P+1,J)
 20      FT(J) = FT(J) + 2.D0 * QP * YK(J)


C*
C*      D = 2   CORE STATES
C*

      CALL YKFCTB(F,F,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,0,NPO,IDR)

      DO 21 J = 1, NPO
 21      FT(J) = FT(J) + 2.D0* QD *YK(J)

      IF(NO_S.EQ.4) RETURN

      END
C#######################################################################
      SUBROUTINE XCHB(P,NL,PR,YK,ZK,WY,WZ,DR,R,F,H,L,NPO,NCORE,IDR)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION P(NL,*), NCORE(*)       
      DIMENSION PR(*), F(*),  DR(*), R(*), YK(*)
      DIMENSION ZK(*), WY(*), WZ(*)
      DIMENSION DRKR(NPO), RKK(NPO), Y1(NPO)

C....................................................

      NO_S = NCORE(1) 
      NO_P = NCORE(2)
      NO_D = NCORE(3) 

C....................................................

      FS   = 2.0D+00                       /  DBLE( 2*L + 1 )
      FPM  = 3.0D+00 * FS  * DBLE( L )     /  DBLE( 2*L - 1 )
      FPP  = 3.0D+00 * FS  * DBLE( L + 1 ) /  DBLE( 2*L + 3 )
      FDM  = 2.5D+00 * FPM * DBLE( L - 1 ) /  DBLE( 2*L - 3 )
      FD   = 5.0D+00 * FPM * DBLE( L + 1 ) /  DBLE( 2*L + 3 ) / 3.0D+00
      FDP  = 2.5D+00 * FPP * DBLE( L + 2 ) /  DBLE( 2*L + 5 )

      LPM  = L - 1
      LPP  = L + 1
      LDM  = L - 2
      LDP  = L + 2

C#    PROBLEM HERE FOR GNU COMPILER
C       WRITE(*,*) ' LDP = ', LDP
C      WRITE(*,*) ' LPP = ', LPP
C      WRITE(*,*) ' LDM = ', LDM
C      WRITE(*,*) ' LPM = ', LPM

      IF(NO_S.EQ.0) RETURN

      DO  J = 1, NPO

         F(J) = P(1,J)
         YK(J) = 0.0D+00

      ENDDO

C*
C*   CALCULATE   R^(K+1) , 1/R^K AT GRID POINTS
C*  

      CALL RKPREP( R, DRKR, RKK, L, NPO)      
      CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NPO,IDR)

      IF(NO_S.EQ.0)  RETURN

      DO  J = 1, NPO
         YK(J) = YK(J) + FS * Y1(J) * F(J)
      ENDDO

      IF(NO_S.EQ.1) RETURN

      DO  J = 1, NPO
         F(J) = P(2,J)
      ENDDO


      CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NPO,IDR)


      DO  J = 1, NPO
         F(J) = P(NO_S+1,J)
         YK(J) = YK(J) + FS * Y1(J) * P(2,J)
      ENDDO


      CALL RKPREP(R,DRKR,RKK,LPP,NPO)
      CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPP,NPO,IDR)

      DO  J = 1, NPO
         YK(J) = YK(J) + FPP * Y1(J) * F(J)
      ENDDO

       IF(LPM.GE.0) THEN

          CALL RKPREP(R,DRKR,RKK,LPM,NPO)
          CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPM,NPO,IDR)

          DO  J = 1,NPO
             YK(J) = YK(J) + FPM * Y1(J) * F(J)
          ENDDO

       ENDIF

      IF(NO_S.EQ.2) RETURN


      DO  J = 1, NPO
         F(J) = P(3,J)
      ENDDO

      CALL RKPREP(R,DRKR,RKK,L,NPO)
      CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NPO,IDR)

      DO  J = 1, NPO
         F(J) = P(NO_S + 2,J)
         YK(J) = YK(J) + FS * Y1(J) * P(3,J)
      ENDDO

      CALL RKPREP(R,DRKR,RKK,LPP,NPO)
      CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPP,NPO,IDR)

      DO  J = 1, NPO

         YK(J) = YK(J) + FPP * Y1(J) * F(J)
      ENDDO

       IF(LPM.GE.0) THEN

          CALL RKPREP(R,DRKR,RKK,LPM,NPO)
          CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPM,NPO,IDR)

          DO J = 1, NPO
             YK(J) = YK(J) + FPM * Y1(J) * F(J)
          ENDDO

       ENDIF

      IF(NO_S.EQ.3) RETURN

      DO  J = 1, NPO
         F(J) = P(4,J)
      ENDDO

      CALL RKPREP(R,DRKR,RKK,L,NPO)
      CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NPO,IDR)

      DO J = 1, NPO
         F(J) = P(NO_S+3,J)
         YK(J) = YK(J) + FS * Y1(J) * P(4,J)
      ENDDO
      
      CALL RKPREP(R,DRKR,RKK,LPP,NPO)
      CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPP,NPO,IDR)

      DO  J = 1, NPO

         YK(J) = YK(J) + FPP * Y1(J) * F(J)
      ENDDO

      IF(LPM.GE.0) THEN

         CALL RKPREP(R,DRKR,RKK,LPM,NPO)
         CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPM,NPO,IDR)

          DO  J = 1, NPO
             YK(J) = YK(J) + FPM * Y1(J) * F(J)
          ENDDO

       ENDIF

       DO  J = 1,NPO
          F(J) = P(NO_S + NO_P + 1, J)
       ENDDO
 
       CALL RKPREP(R,DRKR,RKK,LDP,NPO)
       CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LDP,NPO,IDR)


       DO  J = 1, NPO
          YK(J) = YK(J) + FDP * Y1(J) * F(J)
       ENDDO

       
       CALL RKPREP(R,DRKR,RKK,L,NPO)
       CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NPO,IDR)

       DO  J = 1, NPO
          YK(J) = YK(J) + FD * Y1(J) * F(J)
       ENDDO

       IF(LDM.GE.0) THEN

          CALL RKPREP(R, DRKR, RKK, LDM, NPO)
          CALL YKFCTB(PR,F,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LDM,NPO,IDR)

          DO  J = 1, NPO

             YK(J) = YK(J) + FDM * Y1(J) * F(J)
          ENDDO

       ENDIF

       RETURN
       END
C#######################################################################
      SUBROUTINE YKFCTB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NPO,
     1     IDR)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION P(*),  PP(*), YK(*),ZK(*)
      DIMENSION WY(*), WZ(*), R(*), DR(*), DRKR(*), RKK(*)

      IF(IDR.EQ.0) CALL YKSINB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NPO)

      IF(IDR.EQ.1) CALL YKEXPB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NPO)

      IF(IDR.GT.1) STOP

      RETURN
      END
C#######################################################################
      SUBROUTINE YKEXPB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NPO)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION P(*),PP(*),YK(*),ZK(*)
      DIMENSION WY(*),WZ(*),R(*),DR(*),DRKR(*),RKK(*)


      KR = K 
      KK = K + 1

      IF(KR.EQ.0) GO TO 30

      DO  J = 2, NPO

         WZ(J) = P(J) * PP(J) * DRKR(J)/DR(J)
         WY(J) = P(J) * PP(J) * RKK(J)/DR(J)
      ENDDO

      CALL YINT(WY,WZ,YK,ZK,NPO,H)

      DO  J = 2, NPO

         YK(J) = YK(J) * DRKR(J) + ZK(J) * RKK(J)
      ENDDO

      RETURN

 30   DO 40 J=2,NPO
         WZ(J)=P(J)*PP(J)/DR(J)
 40      WY(J)=P(J)*PP(J)*RKK(J)/DR(J)

         CALL YINT(WY,WZ,YK,ZK,NPO,H)

         DO  J = 2, NPO
            YK(J) = YK(J) + ZK(J) * RKK(J)
         ENDDO

            RETURN
            
            END
C#######################################################################
      SUBROUTINE YKSINB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NPO)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION P(*),PP(*),YK(*),ZK(*)
      DIMENSION WY(*),WZ(*),R(*),DR(*),DRKR(*),RKK(*)


      KR    = K 
      KK    = K + 1
      WZ(1) = 0.0D+00
      WY(1) = 0.0D+00

C*
C*    FOR K = 0  NO INTEGRATION FOR Y_1 IS NEEDED
C*

      IF(KR .EQ. 0) GO TO 30

C*
C*   PREPARE INTEGRAL ARGUMENT Y_1, Y_2
C*       
C*     Y_1  ~    PI(R) * PJ(R) * R^K        
C*       
C*     Y_2  ~    PI(R) * PJ(R) / R^(K+1)           
C*

      DO  J = 2, NPO

         WY(J) = P(J) * PP(J) * RKK(J) / DR(J)
         WZ(J) = P(J) * PP(J) * DRKR(J) / DR(J)

      ENDDO

C*
C*    INTEGRATE NOW 
C*

      CALL YINT(WY, WZ, YK, ZK, NPO, H)


      YK(1)   =  0.0D+00

         DO J = 2, NPO

           YK(J) = YK(J) * DRKR(J) + ZK(J) * RKK(J)  

         ENDDO

      RETURN
C*
C*    IF K = 0  R^(K+1) = R    AND 1/R^K = 1    
C*

   30  DO  J = 2, NPO

          WY(J)  = P(J) * PP(J) * RKK(J) / DR(J)
          WZ(J) =  P(J) * PP(J) / DR(J)

      ENDDO
C*    INTEGRATE NOW

      CALL YINT(WY,WZ,YK,ZK,NPO,H)

      YK(1) = 0.0D+00
         DO J = 2, NPO

           YK(J) = YK(J) + ZK(J) * RKK(J) 

         ENDDO


      RETURN
      END
C#######################################################################
      SUBROUTINE SOLVEMAT( A, B, NS, ER, N)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DOUBLE PRECISION A(NS, *), B(NS, *)
      DOUBLE PRECISION ER(*)
      DOUBLE PRECISION WORK(NS*N)
      CHARACTER UPLO, JOBZ
      INTEGER ITYPE, LWORK, INFO
C..................................

      ITYPE = 1
      LWORK =  NS * N
      UPLO  = "U"
      JOBZ  = "V"
      LDZ   = NS
      INFO  = 0 
      
C#    DSYGV / LAPACK

      CALL DSYGV(ITYPE,JOBZ,UPLO, N, A, NS, B, NS, ER, WORK,LWORK,INFO)


      IF(INFO.NE.0) THEN 

        WRITE(*,*) 'ERROR IN DSYGV LAPACK SUBROUTINE '
        WRITE(*,*) 'INFO = ', INFO
        WRITE(*,*) 'OPTIMAL LWORK IS ', WORK(1) 

        STOP
      ENDIF


      DO  I = 1, N

         ER(I) = - ER(I)
      ENDDO


      RETURN
      END
C#######################################################################
c$$$      SUBROUTINE SETMAT( A, B, AD, AE, NS, N, K, T,
c$$$     1             ALP1, R01, ZNUC, LANG, REDMASS)
c$$$
c$$$      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
c$$$      DIMENSION  A(NS, *),  B(NS, *)
c$$$      DIMENSION AE(NS, *), AD(NS, *)
c$$$      DIMENSION   ALP1(*), R01(*)
c$$$      DIMENSION T(*)
c$$$      DIMENSION DB(K, K), XG(K), WG(K)
c$$$      DIMENSION WORK((K+1)*(K+2)/2)
c$$$C.................................................
c$$$
c$$$*
c$$$*  INITIALIZE  A AND B
c$$$*
c$$$      DO  J = 1, N - 2
c$$$        DO   I = 1, N - 2
c$$$
c$$$          A( I, J ) = AD(I,J) + AE(I,J)
c$$$          B( I, J ) = 0.0D+00
c$$$
c$$$          ENDDO
c$$$       ENDDO
c$$$
c$$$
c$$$      CALL GAUSS(K, XG, WG)
c$$$
c$$$*  SET UP NON ZERO LOOPS  I = 2 ... N-1  ,J = I,I+K-1
c$$$
c$$$      DO 500 I = 2, N - 1
c$$$
c$$$         JHI = MIN0(I+K-1,N-1)
c$$$
c$$$         DO 400 J = I,JHI
c$$$            LOW = MAX0(K,J)
c$$$            LHI = MIN0(I+K-1,N)
c$$$*        SUM OVER INTERMEDIATE SEGMENTS
c$$$            DO 300 L = LOW,LHI
c$$$               SEG = 0.D0
c$$$               SEH = 0.D0
c$$$               DL = (T(L+1)-T(L))
c$$$               X0 = T(L)
c$$$C           SUM OVER GAUSSIAN WEIGHTS
c$$$               DO 200 M = 1, K
c$$$
c$$$                  XM = DL * XG(M) + X0
c$$$
c$$$C#
c$$$C#               DE BOOR'S SOUBROUTINE
c$$$C#               SLATEC LIBRARY
c$$$C#                  
c$$$
c$$$                  CALL DBSPVD(T, K, 2, XM, L, K, DB, WORK)
c$$$
c$$$
c$$$                  FM = DB(I-L+K,1)*DB(J-L+K,1)
c$$$
c$$$                  GM = Q(XM, ALP1, R01, ZNUC, LANG, REDMASS) * FM
c$$$
c$$$                  HM = DB(I-L+K,2)*DB(J-L+K,2)
c$$$
c$$$                  SEG = SEG + WG(M)*FM
c$$$
c$$$C*
c$$$C*                 REDMASS = 1./MASS_REDUCED 
c$$$C*
c$$$C*
c$$$C*                  == 1  HYDROGEN    - LIKE
c$$$C*                  == 2  POSITRONIUM - LIKE
c$$$C*
c$$$C*
c$$$
c$$$                  SEH = SEH + WG(M)* REDMASS * ( GM - HM )
c$$$
c$$$ 200           CONTINUE
c$$$
c$$$               A(I-1,J-1) = A(I-1,J-1) + DL*SEH
c$$$
c$$$               B(I-1,J-1) = B(I-1,J-1) + DL*SEG
c$$$
c$$$ 300        CONTINUE
c$$$ 400     CONTINUE
c$$$ 500  CONTINUE
c$$$
c$$$*   STORE REMAINDER
c$$$      DO 600 I = 1, N - 2
c$$$         DO 550 J = 1, I
c$$$            
c$$$            B(I,J) = B(J,I)
c$$$            A(I,J) = A(J,I)
c$$$
c$$$ 550     CONTINUE
c$$$ 600  CONTINUE
c$$$
c$$$      RETURN
c$$$       END
C#######################################################################
      FUNCTION Q(X, ALP1, R01, ZNUC, LANG, REDMASS)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION ALP1(*), R01(*)
C..............................

      QL = DBLE(LANG*(LANG+1))

      IF(ALP1(LANG+1).EQ.0.0D+00) THEN
         Q  = (2.0D+00 * ( ZNUC/REDMASS )  - QL/X)/X
      ELSE

         IF(R01(LANG+1).EQ.0.0D+00) THEN
            WRITE(*,*)'CUT-OFF PARAMETER EQUALS ZERO R01 = ',R01(LANG+1)  
         ENDIF

         Q  = (2.0D+00 * ( ZNUC/REDMASS )  - QL/X)/X
         Q  = Q + ALP1(LANG+1)*(1.0D+00-DEXP(-(X/R01(LANG+1))**6))/X**4

      ENDIF


      RETURN
      END
C#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C# Improved core-polarization model potential (for Ca)
C# C. Laughlin Physica Scripta, 45, 238, 1992
C#
C#   alp(1) = a_d
C#   alp(2) = a_q - 6*b1
C#   alp(3) = r_a 
C#
      FUNCTION QQ(X, ALP, R0, ZNUC, L, REDMASS)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION ALP(*), R0(*)
      real*8 x, redmass,znuc
      real*8 w6,w8
      real*8 a0,a1,a2,th
      integer l
C..............................

      QL = DBLE(L*(L +1))

      IF(ALP(L + 1).EQ.0.0D+00) THEN

         QQ  = (2.0D+00 * ( ZNUC/REDMASS )  - QL/X)/X

      ELSE

         w6 = 1.0D+00 - exp(-( x/alp(3)**6 ))
         w8 = 1.0D+00 - exp(-( x/alp(3)**8 ))
         vp_4 =  alp(1) * w6/x**4
         vp_6 =  alp(2) * w8/x**6 
         
         
         call model_parameters(znuc,l,a0,a1,a2,th)
         
         vl = ( a0 + a1*x + a2*x*x ) * exp(-th*x) 
         
         IF( R0(L+1).EQ.0.0D+00) THEN

            WRITE(*,*)'CUT-OFF PARAMETER EQUALS ZERO R01 = ',R0(L+1)  
         ENDIF

         QQ  = (2.0D+00 * ( ZNUC/REDMASS )  - QL/X) / X
         QQ  = QQ + vp_4 + vp_6

      ENDIF

!         write(*,*) " a1,a2 = ", alp(1),alp(2),alp(3),alp(4),alp(5)
!         stop

      RETURN
      END
C#
C#  parameters used in the hf-pol potential by C. Laughlin
C#  Physica Scripta Vol. 45,2380245,1992
C#
C#
      subroutine model_parameters(znuc,l,a0,a1,a2,th)
      implicit none
      real*8 znuc
      integer l
      real*8 a0,a1,a2,th

      if(znuc.eq.12) then 
C. mg
         if(l.eq.0) then 
            a0 = -3.12953040d+00
            a1 = 1.15200880d+00
            a2 = 0.0d+00 
            th = 1.6d+00            
         else if(l.eq.1) then
            a0 = -2.19728520d+00
            a1 = 0.60062794d+00
            a2 = 0.0d+00 
            th = 1.5d+00            
         else if(l.eq.2) then
            a0 = -4.09339770d+00
            a1 = 1.44441530d+00
            a2 = 0.0d+00 
            th = 2.2d+00            
         else if(l.eq.3) then
            a0 = -0.13821573d+00
            a1 =  0.04547868d+00
            a2 = 0.0d+00 
            th = 1.6d+00
         else if(l.eq.4) then
            a0 = -0.03410672d+00
            a1 = 0.01260802d+00
            a2 = 0.0d+00 
            th = 1.1d+00
         endif
      else if(znuc.eq.20) then
C     . ca
         if(l.eq.0) then 
            a0 = -3.21619190d+00
            a1 = -1.86020650d+00
            a2 = 0.7d+00
            th = 1.8d+00
         else if(l.eq.1) then
            a0 = -3.87022490d+00
            a1 = 0.27150462d+00
            a2 = 0.05d+00 
            th = 1.8d+00            
         else if(l.eq.2) then
            a0 = 15.533477d+00
            a1 = -30.077412d+00
            a2 = 0.1d+00 
            th = 2.9d+00
         else if(l.eq.3) then
            a0 = -6.46541250d+00
            a1 =  1.172720040d+00
            a2 = 0.2d+00 
            th = 1.8d+00
         else if(l.eq.4) then
            a0 = -8.95074570d+00
            a1 = 4.41149130d+00
            a2 = -0.5d+00 
            th = 1.8d+00
         endif
      else if(znuc.eq.11) then
C     . na
         if(l.eq.0) then
            a0 = -3.56981940d+00
            a1 = 0.97882157d+00
            a2 = 0.0d+00 
            th = 1.8d+00
         else if(l.eq.1) then
            a0 = -2.39201830d+00
            a1 = -0.29097386d+00
            a2 = 0.0d+00 
            th = 1.9d+00            
         else if(l.eq.2) then            
            a0 = -2.39657490d+00
            a1 = 0.66112582d+00
            a2 = 0.0d+00 
            th = 1.8d+00
         else if(l.eq.3) then
            a0 = 0.11340957d+00
            a1 = 0.006850348d+00
            a2 = 0.0d+00 
            th = 1.2d+00
         else if(l.eq.4) then
            a0 =  0.66604094d+00
            a1 = -0.08410312d+00
            a2 = 0.0d+00 
            th = 0.7d+00
         endif         
      else if(znuc.eq.19) then
C     . k
         if(l.eq.0) then
            a0 = -1.93873490d+00
            a1 =  0.65255104d+00
            a2 = -0.055d+00 
            th = 0.8d+00
         else if(l.eq.1) then
            a0 = -0.41538107d+00
            a1 = -0.55613230d+00
            a2 = 0.1d+00 
            th = 0.8d+00            
         else if(l.eq.2) then            
            a0 = -1.291847d+00
            a1 =  0.42647161d+00
            a2 = -0.032d+00 
            th =  0.6d+00
         else if(l.eq.3) then
            a0 = 5.85015660d+00
            a1 = -1.33154540d+00
            a2 = 0.007d+00 
            th = 0.6d+00
         else if(l.eq.4) then
            a0 = -0.33270253d+00
            a1 = -8.41885410d+00
            a2 = 0.5d+00 
            th = 0.6d+00
         endif         
      endif

      end subroutine
C#####################################################################
C#####################################################################
C     #EOF





