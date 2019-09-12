C###################################################################
      SUBROUTINE SUBMX(P1,PR,P3,PC,R,DR,VX,ISL,NHF,NLLMIN,NOLL,LHF,LL,YK
     1     ,EHF,IR,IC,LO,NO,H,IDR)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.2e.inc"
C     PARAMETER(NS=802, NP=2000, NHX=2500, NL=15)
C     np,ns,nl      
      DIMENSION P1(NP),PR(NS,NP),P3(NP),PC(NS,NP),R(NP),DR(NP)
      DIMENSION VX(NS,NS),NHF(NS),NOLL(1),LHF(1),ANGD(NL),KAD(NL)
      DIMENSION ANGE(NL),KAE(NL),YD(NP),YE(NS,NP),YK(1),P2(NP)
      DIMENSION P4(NP),LL(1),NLLMIN(1), EHF(NL,NS)
      COMMON/RARR/RPW(NL,NP)
      COMMON/RRO/rho(6),alpha,beta
C*  
C*    calculate coefficient of 'direct'/excahnge term and phase factor

      CALL ANGFK(LHF(IR),LL(IR),LHF(IC),LL(IC),LO,ANGD,KAD,KNOD)
      CALL ANGFK(LHF(IR),LL(IR),LL(IC),LHF(IC),LO,ANGE,KAE,KNOE)
      PHA =(-1)**( IABS( LO + LHF(IC) - LL(IC) ) )


C*    zero the relevant matrix before starting calculation
 
      DO  NR = 1, NS
         DO  NC = 1, NS
            VX( NR, NC ) = 0.0D+00
         ENDDO
      ENDDO

C*
C*   calculate Y_k(j_2) for the 'direct' term and sum over 'k'
C*   electron '1' coordinates integration
C*

      DO  J = 1, NO
         YD(J) = 0.0D+00
      ENDDO

      
      DO 10 K = 1, KNOD

         CALL SET_YK( P1, P3, YK, R, DR, KAD(K), NO, H, IDR)

         DO  J = 1, NO

            YD(J) = YD(J) + ANGD(K) * YK(J)
         ENDDO

 10   CONTINUE

C*
C*   calculate Y_k(j_2) for the 'exchange' term and sum over 'k'
C*   electron '1' coordinates integration
C*

      DO 20 N = 1, NOLL(IC)
         

         DO J = 1, NO

            YE(N,J) = 0.0D+00
            P4(J)   = PC(N, J)

         ENDDO


         DO 20 K = 1, KNOE
         
            CALL SET_YK(P1, P4, YK, R, DR, KAE(K), NO, H, IDR)

            DO  J = 1, NO

               YE(N,J) = YE(N,J) + ANGE(K) * YK(J)
            ENDDO

 20      CONTINUE

C##############################################


         DO 40 NR = 1, NOLL(IR)


C...........            state    b = | n2, l2 >
C#
C#
C#           < ab | == < n1 l1 n2 l2 | 
C#
C#

            DO  J = 1, NO

               P2(J) = PR(NR, J)

            ENDDO



         DO 40 NC = 1, NOLL(IC)


C...........            state    d = | n4, l4 >
C#
C#
C#           | cd > == | n3 l3 n4 l4  > 
C#
C#


            DO  J = 1, NO

               P4(J) = PC(NC, J)
            ENDDO


C#
C#         V_12(ab,cd) =  < A_12 * (n1 l1, n2 l2 ) | V(1,2) | A_12 *( n3 l3, n4 l4) >
C#



C     *   calculate Y_k for the 'direct' term and sum over 'k'
C     *   electron '2' coordinates integration
C*

            YK(1) = 0.0D+00
            IF(DR(1).NE.0.D0) YK(1) = YD(1) * P2(1) * P4(1) / DR(1)

            DO  J = 2, NO

               YK(J) = YD(J) * P2(J) * P4(J) / DR(J)
            ENDDO

            VD = RINT( YK, 1, NO, 14, H)
            VD = 2.0D+00 * VD

C     *   calculate Y_k for the 'exchange' term and sum over 'k'
C     *   electron '2' coordinates integration

            YK(1) = 0.0D+00

            IF(DR(1).NE.0.D0) YK(1) = YE(NC,1) * P2(1) * P3(1) / DR(1)
            
            DO  J = 2, NO

               YK(J) = YE(NC,J) * P2(J) * P3(J) / DR(J)
            ENDDO

            VE = RINT(YK, 1, NO, 14, H) 

            VE = PHA * 2.0D+00 * VE
 
C....................

C#
C#    vmd    :      Mass polarization - 'direct' part
C#    vme    :      Mass polarization - 'exchange' part
C#    vdd    :      dielectronic polarization  - 'direct' part 
C#    vde    :      dielectronic polarization  - 'exchange' part 
C#
C#
C#  vd,     ve  : 1/r_12                                ci
C#  vdd,    vde : - a_d * r_1 * r_2 * V(r_1) * V(r_2),  diel. pot.
C#  vmd,    vme : - p_1 * p_2                           Ma polarization pot.
C#

C...............................................
      
            LHF1 = LHF(IR) + 1
            LHF3 = LHF(IC) + 1

            LL2  = LL(IR)  + 1
            LL4  = LL(IC)  + 1

c            write(*,*) ' nr, nc = ', nr, nc 

            e1 = EHF( LHF1, NHF(IR))

c            if(lhf(ir).eq.0) then 
c               write(*,*) ir, lhf(ir), nhf(ir), e1
c            endif
               
            e2 = EHF( LL2, NR)
            e3 = EHF( LHF3, NHF(IC))
            e4 = EHF( LL4, NC )
            

            n1 = nhf(ir) + lhf(ir)
            l1 = lhf(ir)

            n2 = nr + ll(ir)
            l2 = ll(ir)

            n3 = nhf(ic) + lhf(ic)
            l3 = lhf(ic)

            n4 = nc + ll(ic)
            l4 = ll(ic)
            
C................................................
            
            
            VDD = 0.0D+00
            VDE = 0.0D+00
            
            RHO1 = RHO( LHF(IR) )
            RHO2 = RHO( LL(IR)  )

C     *   'direct term'
            
            do 55 k = 1, knod
               
               if(kad(k).ne.1) goto 55
               
               VDD = VDD - 2.0D+00 * alpha * angd(k)
     1              * VDR( P1, P3, RHO1, NO )
     1              * VDR( P2, P4, RHO2, NO )
               
 55         continue
            
            
            
C     *   'exchange term'
            
            do 56 k = 1, knoe
               
               if(kae(k) .ne. 1) goto 56
               
               VDE  =  VDE - 2.0D+00 * alpha * PHA * ange(k)
     1              * VDR( P1, P4, RHO1,  NO )
     1              * VDR( P2, P3, RHO2,  NO )
               
 56         continue
            
C     #    skip Vm 
            
C     beta= 0.0D+00
            
            
            
C     #
C     #                                      1    
C     #           V( 1s2p | 1s2p)  = +- ---------- | d_v_radial(1s,2p) |^2   
C     #                                 3( M + m) 
C     #
C     #              M      : nucleus mass
C     #              m      : electron mass
C     #              d_v_r  : radial matrix element in velocity gauge
C     #
C     #          For  Ps : 
C     #                        M = 1, m = 1, 
C     #
C     #               d_v_radial(a,b) = mu dE_ab * d_r_radial(a,b) 
C     #
C     #
C     #           Hydrogen     ---> d_r_radial( H / 1s,2p) = 1.29 a.u 
C     #           Positronium  ===> d_r_radila(Ps / 1s,2p) =     
C     #


C...............................................


            IF(BETA.EQ.1) THEN 

               BETA = 0.5D+00

            ENDIF
            
            VMD = 0.0D+00
            VME = 0.0D+00

C     *   'direct term'
            vvv = 0.0D+00

            DE13 = 0.5D+00 * ( EHF( LHF1, NHF(IR) ) - EHF(LHF3,NHF(IC)))
            DE24 = 0.5D+00 * ( EHF( LL2, NR ) -  EHF(LL4, NC))

      do 57 k = 1, knod

         if(kad(k).ne.1) goto 57

C         if(n1.eq.1.and.n2.eq.1.and.n3.eq.2.and.n4.eq.2.and.l1.eq.0
C     1 .and.l2.eq.0.and.l3.eq.1.and.l4.eq.1)       then

         if(e1.lt.0.0D+00.and.e2.lt.0.0D+00.and.e3.lt.0.0D+00
     1        .and.e4.lt.0.0D+00) then

            VMD = VMD -  2.0D+00 * angd(k)
     1           * BETA * DE13 * VMR( P1, P3, NO )
     1           * BETA * DE24 * VMR( P2, P4, NO )
            
            vvv =    2.0D+00 * angd(k)
     1           * BETA * DE13 * VMR( P1, P3, NO )
     1           * BETA * DE24 * VMR( P2, P4, NO )

          

C     IF(l3.eq.1.and.l4.eq.1) then

C     IF((n1+1).eq.n3.and.(n2+1).eq.n4)
            
            if(ABS(vvv).gt.1) then


            WRITE(*,*) ' LLHF3, LL4  == ', LHF3, LL4
            WRITE(*,*) ' LHF1,  LL2  == ', LHF1, LL2
            WRITE(*,*) '      IR, IC == ', IR, IC
            WRITE(*,*) '      NR, NC == ', NR, NC

            
               WRITE(*,*) ' xxxxxxxxxxxxxx  DIRECT PART xxxxxxxxxxxxxx'
               WRITE(*,*) nr,nc,ir,ic
               WRITE(*,'(a1,4i3,a4,4i3,a5,E25.14)') 
     1               '(',n1,l1,n2,l2,'|VMD|',n3,l3,n4,l4,' > = ',vvv
               WRITE(*,*) '          a(k),k = ', angd(k), k 
               WRITE(*,*) ' DE(1,3), D(2,4) = ', DE13, DE24
               WRITE(*,*) '   VMR(1,3),r,v  = ', VMR( P1, P3, NO ), 
     1              BETA * DE13 * VMR( P1, P3, NO )
               WRITE(*,*) '   VMR(2,4),r,v  = ', VMR( P2, P4, NO ),
     1              BETA * DE24 * VMR( P2, P4, NO )
               
            endif
            
         endif
      
c         write(*,*) 'vmd = ', vmd, 'xxxxxxxxx'
 57   continue


       
C*   'exchange term'

      DE14 = 0.5D+00 * ( EHF( LHF1, NHF(IR) ) - EHF( LL4, NC ) )
      DE23 = 0.5D+00 * ( EHF( LL2, NR ) -  EHF(LHF3, NHF(IC)) )

      do 58 k = 1, knoe

         if(kae(k).ne.1) goto 58


c         if(n1.eq.1.and.n2.eq.1.and.n3.eq.2.and.n4.eq.2.and.l1.eq.0
c     1 .and.l2.eq.0.and.l3.eq.1.and.l4.eq.1)       then

         
         if(e1.lt.0.0D+00.and.e2.lt.0.0D+00.and.e3.lt.0.0D+00
     1        .and.e4.lt.0.0D+00) then
         

         VME  =  VME - 2.0D+00 * PHA  * ange(k)      
     1                         * BETA * DE14 * VMR( P1, P4, NO)
     1                         * BETA * DE23 * VMR( P2, P3, NO)




         vvv  = 2.0D+00 * PHA  * ange(k)      
     1                         * BETA * DE14 * VMR( P1, P4, NO)
     1                         * BETA * DE23 * VMR( P2, P3, NO)


C         IF(l3.eq.1.and.l4.eq.1) then

c$$$         IF((n1+1).eq.n3.and.(n2+1).eq.n4) then
c$$$


         if(ABS(vvv).gt.1) then

            write(*,*) nr, ir, ll(ir),lhf(ir), nhf(ir), e1, e2
            write(*,*) nc, ic, ll(ic),lhf(ic), nhf(ic), e3, e4

 
            WRITE(*,*) ' LLHF3, LL4  == ', LHF3, LL4
            WRITE(*,*) ' LHF1,  LL2  == ', LHF1, LL2
            WRITE(*,*) '      IR, IC == ', IR, IC
            WRITE(*,*) '      NR, NC == ', NR, NC

         do j = 1, no

            write(40,*) r(j), p1(j), p2(j) 
            write(41,*) r(j), p3(j), p4(j) 

         enddo

            WRITE(*,*) ' xxxxxxxxxxxxxx  EXCHANGE PART xxxxxxxxxxxxxxxx'
            WRITE(*,*) nr,nc,ir,ic
            WRITE(*,'(a1,4i3,a4,4i3,a5,E25.14)') 
     1       '(',n1,l1,n2,l2,'|VME|',n3,l3,n4,l4,' > = ',  vvv
C            WRITE(*,*) '          a(k),k = ', ange(k), k 
C            WRITE(*,*) ' DE(1,4), D(2,3) = ', DE14, DE23
C            WRITE(*,*) '    VMR(1,4),r,v = ', VMR( P1, P4, NO),
C     1           BETA * DE14 * VMR( P1, P4, NO)
C            WRITE(*,*) '    VMR(2,3)r,v = ', VMR( P2, P3, NO ),
C     1           BETA * DE23 * VMR( P2, P3, NO)

         endif
      endif


c          write(*,*) 'vmd = ', vme, 'xxxxxxxxx'
 58   continue

      if(isl.eq.1) then

C#                  --1/r12-- -- v1*v2 --  -- r1*r2 --                      

         VX(NR,NC) = VD + VE  + VMD + VME + VDD + VDE

         if(n1.eq.1.and.n2.eq.1.and.n3.eq.2.and.n4.eq.2.and.l1.eq.0
     1 .and.l2.eq.0.and.l3.eq.1.and.l4.eq.1)  then

c$$$         IF((n1+1).eq.n3.and.(n2+1).eq.n4) then
c$$$
c$$$C         IF(l3.eq.1.and.l4.eq.1) then

c         WRITE(*,*) ' VD,   VE = ',  VD,  VE,  VD  + VE
c         WRITE(*,*) ' VDD, VDE = ',  VDD, VDE, VDD + VDE
c         WRITE(*,*) ' VMD, VME = ',  VMD, VME, VMD + VME
               WRITE(*,'(a1,4i3,a4,4i3,a5,3E25.14)') 
     1               '<',n1,l1,n2,l2,'|V|',n3,l3,n4,l4,' > = ',vx(nr,nc)
     1                                           , vd,ve
         WRITE(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

         endif

      else if(isl.eq.3) then 

         VX(NR,NC) = VD - VE + VMD - VME + VDD - VDE
      else
         write(*,*)'only singlet (1) and triplet (2) states are allowed'
      endif

C     #    Normalize properly when configurations are 'equivalent'

      CALL Q2TST1(LHF,LL,NHF,NLLMIN, IR, IC, NR, NC,VX(NR,NC))


 40   CONTINUE

      IF(IR .NE. IC) GO TO 48

C#
C#    Add the non-interacting part of the CI matrix :
C#
C#    H = e_1 + e_2 + V(aa)
C#
      LLR   = LL(IR)  + 1
      LLHFR = LHF(IR) + 1 
      NHFIR = NHF(IR)

      DO  KK = 1, NOLL(IR)
         
         NK = KK + NLLMIN(IR) - 1

         VX(KK, KK) = VX(KK, KK) + EHF(LLHFR, NHFIR) + EHF(LLR, NK)
      ENDDO
         
 48   CONTINUE

      RETURN
      END
C#########################################################
      SUBROUTINE Q2TST1(LHF,LL,NHF,NLLMIN,IR,IC,NR,NC,VXS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION LHF(1),LL(1),NHF(1),NLLMIN(1)
C.........................................................

      IF(NC.GT.1.AND.NR.GT.1) RETURN

C      if(nhf(ir).eq.nllmin(ir).and.lhf(ir).eq.ll(ir)) then
C         if(nr.gt.1) goto 10
C	vxs = vxs / dsqrt(2.0D+00)
C 10   endif

C      if(nhf(ic).eq.nllmin(ic).and.lhf(ic).eq.ll(ic)) then
C         if(nc .gt. 1) goto 11
C         vxs = vxs/dsqrt(2.0d+00)
C 11   endif


      if(nr.eq.1) then 
         if(nhf(ir).eq.nllmin(ir).and.lhf(ir).eq.ll(ir)) then
            vxs = vxs /dsqrt(2.0D+00)
         endif

      endif

      if(nc.eq.1) then 
         if(nhf(ic).eq.nllmin(ic).and.lhf(ic).eq.ll(ic)) then
            vxs = vxs/dsqrt(2.0d+00)
         endif         
      endif
      
      return
      end
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C.
C.         calculates the radial integral 
C.
C.            < n l| r | n' l' > 
C.
C.        where v_r the velocity operator
C.
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      FUNCTION VMR(PR, PC, NO)
      IMPLICIT REAL*8(a-h,o-z)
      INCLUDE "parameter.2e.inc"
C      PARAMETER(NP = 2000 )
      DIMENSION PR(NP),PC(NP),F(NP)
      COMMON/RH/H, R(NP), DR(NP)
C.................................................


      F(1) = 0.0D+00

      IF(DR(1).NE.0.0D0) F(1) = PR(1) * PC(1) * R(1)/DR(1)
      
      DELR = 0.0D+00

      IF(R(1).NE.0.0D+00) DELR = 0.5D+00 * F(1) * R(1)

      DO j = 2, NO

         F(J) = PR(J) * PC(J) * R(J)**2 / DR(j)

      END DO

      VMR = RINT( F, 1, NP, 14, H ) + DELR

      RETURN
      END
C################################################
c      function vdr(PR,LR,PC,LC,RHO,NO)
      function vdr(PR,PC,RHO,NO)
      implicit real*8(a-h,o-z)
      INCLUDE "parameter.2e.inc"
C      PARAMETER(NP=2000)
      DIMENSION PR(NP),PC(NP),VP(NP),F(NP)
      common/RH/H,R(NP),DR(NP)
C...............................................

      CALL VPOL(R, VP, RHO, H, NO)

      F(1) = 0.0D+00
      DO j = 2, NP

         F(J) = PR(J) * PC(J) * VP(J) / DR(j)
         
      END DO

      VDR = RINT( F, 1, NP, 14, H )

      RETURN
      END
C##########################################
      SUBROUTINE VPOL(R, VP, RHO, H, NO)
      INCLUDE "parameter.2e.inc"
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(NP),VP(NP)
C.........................................

      DO 20 J = 2, NP

         VP(J) = 1.0D+00
         A = (R(J) / RHO)**6

      IF(A.GT.1.0D+02) GO TO 20

      VP(J) = DSQRT( 1.D0 - DEXP(-A) )
 20   VP(J) = VP(J) / R(J)

      RETURN
      END
C############################################
      SUBROUTINE WFHFIN( P, L, NO, NHF)
      INCLUDE "parameter.2e.inc"
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(NP)
C     
      call wf1eb_file(2, l) 
      DO k = 1, nhf
        READ(2) (P(J), J = 1, NP)
      enddo
      CLOSE(2)
      
      RETURN
      END
C#############################################
      SUBROUTINE WFNLIN( PP, P, L, NO, NMIN, NMX)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.2e.inc"
C      PARAMETER(NP = 2000, NS = 802)
      DIMENSION PP(NS, NP), P(NP)
C...........................................



      call wf1eb_file(2, L) 


      KMX  = NMX - NMIN + 1
      KMIN = NMIN - 1

!      print*, kmx, nmx, nmin

      IF(KMIN.EQ.0) GO TO 20

      DO  K = 1, KMIN
         READ(2) (P(J), J = 1, NP)
      ENDDO

 20   DO  K = 1, KMX
         READ(2) (P(J), J = 1, NP)

         DO  J = 1, NO
            PP(K,J) = P(J)
         ENDDO
      
       ENDDO

       CLOSE(2)
   
       RETURN
       END
C#
C#
C#
      SUBROUTINE Q2TST(LHF,LL,NHF,NLLMIN,IR,IC,NR,NC,VXS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION LHF(1),LL(1),NHF(1),NLLMIN(1)
C.........................................................
      IF(NC.GT.1 .AND. NR.GT.1) RETURN
      IF(LHF(IR).NE.LL(IR) .AND. LHF(IC).NE.LL(IC)) RETURN
      IF(NHF(IR).NE.NLLMIN(IR) .AND. NHF(IC).NE.NLLMIN(IC)) RETURN
      IF(LHF(IR).NE.LL(IR) .AND. NHF(IC).NE.NLLMIN(IC)) RETURN
      IF(LHF(IC).NE.LL(IC) .AND. NHF(IR).NE.NLLMIN(IR)) RETURN
      IF(LHF(IR).NE.LL(IR) .AND. NHF(IC).EQ.NLLMIN(IC)) GO TO 10
      IF(LHF(IC).NE.LL(IC) .AND. NHF(IR).EQ.NLLMIN(IR)) GO TO 11
      IF(NHF(IR).EQ.NLLMIN(IR) .AND. NHF(IC).EQ.NLLMIN(IC)) GO TO 12
      STOP
 10   IF(NC .EQ. 1) GO TO 20
      RETURN
 11   IF(NR .EQ. 1) GO TO 20
      RETURN
 20   VXS = VXS / DSQRT(2.0D0)

      RETURN
 12   IF(NR .EQ. 1 .AND. NC .NE. 1) GO TO 20
      IF(NR .NE. 1 .AND. NC .EQ. 1) GO TO 20

      VXS = VXS/2.0D+00

      RETURN
      END
C##############################################
      SUBROUTINE TRNSMX(HMX,XM,IHR,IHC,IR,IC)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.2e.inc"
      DIMENSION hmx(ns,nhx),xm(ns,ns)
C..............................................
 
      IRF = IHR + IR - 1
      ICF = IHC + IC - 1

      DO  K = IHR, IRF
        DO  KK = IHC, ICF
          
          KR = K  - IHR + 1
          KC = KK - IHC + 1
          HMX(K,KK)=XM(KR,KC)
        enddo
      enddo
      RETURN
      END
C##########################################
      SUBROUTINE MXPRNT(XM,IR,IC)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.2e.inc"
C      PARAMETER(NS=802)
      DIMENSION XM(NS,NS)
C..........................................
    1 FORMAT(4X,1P9D14.5)
    2 FORMAT(2X)
C..........................................
      WRITE(16,2)
      WRITE(*,2)
      DO 10 K=1,IR
   10 WRITE(*,1) (XM(K,KK),KK=1,IC)
      RETURN
      END
C####################################################################
      SUBROUTINE CFIN(NCFG, NHF,LHF,LL,NMIN,NMX,NOL,IS,L,ISL,NCSMX,NCMX)
      DIMENSION NHF(1),LHF(1),LL(1),NMIN(1),NMX(1),NOL(1),IS(1)
C...............................................
    3 FORMAT(5I5)
    6 FORMAT(2X,'L=',I3/)
C...............................................

      READ(NCFG ,3) L, ISL
      READ(NCFG, 3) NCSMX

      NCMX = 0
      DO  K = 1, NCSMX

      READ(NCFG, 3) NHF(K),LHF(K),LL(K),NMIN(K),NMX(K)

      NOL(K) = NMX(K) - NMIN(K) + 1
      IS(K)  = NCMX + 1
      NCMX   = NCMX + NOL(K)
      ENDDO

      if(isl.eq.0.or.isl.eq.1) then 

         isl = 2*isl + 1 
      else

         isl = isl
      endif

      RETURN
      END
C####################################################################
      SUBROUTINE WF1EB_FILE(N,L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE WF1EFILE. EXITING...'
         STOP

      ELSE

C         WRITE(*,*) 'WF1E  :: DATA FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(N,FILE='dat/w1e-0.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/w1e-1.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/w1e-2.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/w1e-3.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/w1e-4.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/w1e-5.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/w1e-6.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/w1e-7.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/w1e-8.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/w1e-9.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/w1e-10.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/w1e-11.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/w1e-12.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/w1e-13.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/w1e-14.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
      END
C#EOF
