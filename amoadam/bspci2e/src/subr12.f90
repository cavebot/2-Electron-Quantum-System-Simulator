!
!
!
!
      SUBROUTINE SUBMX(P1,PR,P3,PC,R,DR,VX,ISL,NHF,NLLMIN,NOLL,LHF,LL,&
     &     YK, EHF, IR, IC, LO, NPO, H, IDR, RHO, ALPHA,BETA,NS,NL)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

      DIMENSION  PR(NS, *),   PC(NS, *),   EHF(NL, *)
      DIMENSION      P1(*),       P3(*),        YK(*) 
      DIMENSION       R(*),        DR(*)
      DIMENSION  NHF(*), NLLMIN(*), NOLL(*), LHF(*), LL(*) 
      DIMENSION  RHO(*)
      DIMENSION  VX(NS, NS),   YE(NS, NPO)
      DIMENSION      P2(NPO),      P4(NPO)
      DIMENSION    ANGD(NL),     KAD(NL)
      DIMENSION    ANGE(NL),     KAE(NL)
      DIMENSION      YD(NPO)

!..........................................................

!..........................................................
!#
!#        V_12(AB,CD) = < N1 L1, N2 L2 | 1 / R_12 | N3 L3, N4 L4 > 
!#  
!#
!#                     ___
!#       V_12(AB,CD) = \     A_K * R_K(AB,CD) +- B_K * R_K(AB;DC)
!#                     /__
!#



!    ANGULAR COEFFICIENT OF 'DIRECT/EXCHANGE' TERM :   a_k, b_k
!    and                            phase factor   :   pha

      CALL ANGFK(LHF(IR),LL(IR),LHF(IC),LL(IC),LO,ANGD,KAD,KNOD)
      CALL ANGFK(LHF(IR),LL(IR),LL(IC),LHF(IC),LO,ANGE,KAE,KNOE)

      PHA =(-1)**( IABS( LO + LHF(IC) - LL(IC) ) )


!*   VX = V_12(AB,CD)  

      vx = 0.0_dpk


!.. r_k(ab,cd) :  direct term



!*   CALCULATE Y_K(J_2) FOR THE 'DIRECT' TERM AND SUM OVER 'K'
!*   ELECTRON '1' COORDINATES INTEGRATION




      yd = 0.0_dpk      
      direct_term:DO  K = 1, knod

         CALL SET_YK( P1, P3, YK, R, DR, KAD(K), NPO, H, IDR)

         yd(:) = yd(:) + angd(k) * yk(:) 

!         DO  J = 1, NPO
!            YD(J) = YD(J) + ANGD(K) * YK(J)
!         ENDDO

      ENDDO direct_term


!... r_k(ab,dc) : exchange term


!*   CALCULATE Y_K(J_2) FOR THE 'EXCHANGE' TERM AND SUM OVER 'K'
!*   ELECTRON '1' COORDINATES INTEGRATION


      exchange_term:DO  N = 1, NOLL(IC)
         
         ye(n,:) = 0.0_dpk       ! 1:npo
         p4(:) = pc(n,:)         ! 1:npo


         DO  K = 1, knoe
         
            CALL set_yk(P1, P4, YK, R, DR, KAE(K), NPO, H, IDR)

            ye(n,:) = ye(n,:) + ange(k) * yk(:) 

!            DO  J = 1, NPO
!               YE(N,J) = YE(N,J) + ANGE(K) * YK(J)
!            ENDDO

         ENDDO
      ENDDO excahnge_term


!#
!#   v_12(ab,cd) =  < A_12 * (N1 L1, N2 L2 ) | V(1,2) | A_12 *( N3 L3, N4 L4)>
!#


      DO  NR = 1, NOLL(IR)             !# orbital    b = | n2, l2 >
                                          !#           < AB | == < N1 L1 N2 L2 | 
            p2(:) = pr(nr,:)              ! 1:npo

            DO  NC = 1, NOLL(IC)          ! orbital    D = | N4, L4 >
                                          !#| CD > == | N3 L3 N4 L4  > 


               p4(:) = pc(nc,:)        ! 1:npo

!'d' term
            YK(1) = 0.0_dpk

            yk = yd * p2 * p4 / dr        ! 2:npo

            IF(dr(1).NE.0.0_dpk) THEN 
               yk(1) = yd(1) * p2(1) * p4(1) / dr(1)
            ENDIF


            vd = 2.0_dpk * rint( yk, 1, npo, 14, h)
!'x' term


            YK(1) = 0.0_dpk

            YK(:) = YE(NC,:) * P2(:) * P3(:) / DR(:)     !2:npo

            IF(dr(1).NE.0.0_dpk) THEN 
               yk(1) = ye(nc,1) * p2(1) * p3(1) / dr(1)
            ENDIF
            

            ve = pha * 2.0_dpk * rint(yk, 1, npo, 14, h) 


            !ls = 0 (singlet), 1= (triplet)

            vx(nr,nc) = vd + (1-2*ls) * ve


            ! for equivalent cfg normalize properly

         
            CALL Q2TST1(LHF, LL, NHF, NLLMIN, IR, IC, NR, NC, VX(NR,NC))


         ENDDO
      ENDDO

      
      !     #    ADD THE NON-INTERACTING PART OF THE CI MATRIX :
      !     #
      !     #               H = H(E_1) + H(E_2) + 1/R_12
      !     #
      
      LLR   = LL(IR)  + 1
      LLHFR = LHF(IR) + 1
      NHFIR = NHF(IR)
      
      DO  KK = 1, NOLL(IR)
         NK = KK + NLLMIN(IR) - 1
         VX(KK, KK) = VX(KK, KK) + EHF(LLHFR, NHFIR) + EHF(LLR, NK)
      ENDDO

      RETURN
    END SUBROUTINE SUBMX
    !#######################################################################


      SUBROUTINE Q2TST1(LHF,LL,NHF,NLLMIN,IR,IC,NR,NC,VXS)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION LHF(*), LL(*), NHF(*), NLLMIN(*)

!.........................................................

      IF(NC.GT.1.AND.NR.GT.1) RETURN

      IF(NHF(IR).EQ.NLLMIN(IR).AND.LHF(IR).EQ.LL(IR)) THEN

         IF(NR.GT.1) GOTO 10

         VXS = VXS / DSQRT(2.0D+00)

 10   ENDIF

      IF(NHF(IC).EQ.NLLMIN(IC).AND.LHF(IC).EQ.LL(IC)) THEN

         IF(NC.GT.1) GOTO 11

         VXS = VXS/DSQRT(2.0D+00)

11    ENDIF
      
      RETURN
      END
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      SUBROUTINE WFHFIN( P, L, NPO, NHF)
        !
        IMPLICIT DOUBLEPRECISION(A-H,O-Z)
        DIMENSION P(*)
        !.........
        
      CALL WF1EFILE(2, L) 
      
      DO  K = 1, NHF         
         READ(2) (P(J), J = 1, NPO)
      ENDDO
      CLOSE(2)
      
      RETURN
      END
      !#######################################################
      SUBROUTINE WFNLIN( PP, NS, P, L, NPO, NMIN, NMX)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

      DIMENSION PP(NS, *), P(*)
      !...........................................

      CALL WF1EFILE(2, L) 

      KMX  = NMX - NMIN + 1
      KMIN = NMIN - 1

      IF(KMIN.EQ.0) GO TO 20

      DO  K = 1, KMIN

         READ(2) (P(J), J = 1, NPO)

      ENDDO

 20   DO  K = 1, KMX

         READ(2) (P(J), J = 1, NPO)

         DO  J = 1, NPO
            PP(K,J) = P(J)
         ENDDO
      
      ENDDO

      CLOSE(2)

      RETURN
      END
      !
      ! subroutine
      !
      SUBROUTINE TRNSMX(v, v_bs, ns, ihr, ihc, ir, ic)
        !exe!
        REAL(dpk) DIMENSION(ns,nhx) :: h
        REAL(dpk) DIMENSION(ns,ns)  :: v_bs
        INTEGER                     :: ns, ihr, ihc,ir,ic
        !exe!
 

      IRF = IHR + IR - 1
      ICF = IHC + IC - 1

      DO  K = IHR, IRF
         DO  KK = IHC, ICF

            KR = K  - IHR + 1
            KC = KK - IHC + 1

            h(k, kk) = v_bs(kr, kc)
         ENDDO
      ENDDO
      
      RETURN
    END SUBROUTINE TRNSMX
    !eof
