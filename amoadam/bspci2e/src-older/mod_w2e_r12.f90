MODULE w2e_r12
  
  USE PRECISION, only: dpk
  IMPLICIT NONE

  PUBLIC diagonalize, read_cfg_file, submx, wfhfin, wfnlin, trnsmx

CONTAINS

SUBROUTINE diagonalize(nh2eb, nhx, h, u)
  !
  !  USE PRECISION, only:dpk
  !
  IMPLICIT NONE 
  !     
  INTEGER,            INTENT(in) :: nh2eb
  INTEGER,             INTENT(in) :: nhx
  REAL(dpk),  DIMENSION(nhx,nhx)  :: h
  REAL(dpk),  DIMENSION(nhx,nhx)  :: u
  !
  REAL(dpk),  DIMENSION(nhx)      :: c
  !
  !external subroutine entries
  REAL(dpk), DIMENSION(5*nhx)     :: iwork
  REAL(dpk), DIMENSION(10*nhx)     :: work
  INTEGER,   DIMENSION(nhx)       :: ifail
  INTEGER                         :: lwork, il, iu, info, m
  REAL(dpk)                       :: abstol, vl,vu
  !
  INTEGER                         :: i,j,k

  !exe!


  !**  lapack routine
  
  lwork  = SIZE(work)
  vl     = 0.0D+00
  vu     = 0.0D+00
  il     = 1
  iu     = nhx
  abstol = 0.0D+00
  
!diagonalize

  CALL dsyevx('V','I','U', nhx, h, nhx, vl, vu, il, iu, abstol, m,&
       & c, u, nhx, work, lwork, iwork, ifail, info)

! and save

!nhx = size(c)

  WRITE(nh2eb) SIZE(c)
  save_eigen_values_coe:DO  k = 1, SIZE(c)
     WRITE(nh2eb)  c(k) 
     WRITE(nh2eb) (u(j,k), j = 1, SIZE(c) )
  ENDDO save_eigen_values_coe

  DO k = 1, nhx 
     WRITE(*,'(i6,2x,2E18.8)') k, c(k)/2.,c(k)*27.211396/2.
  ENDDO
!  CLOSE(nh2eb)
  
!!%  DO k = 1, nhx
!!%
!!%     k1 = 0
!!%     sumProb10 = 0.0d+00
!!%
!!%     DO  i = 1, 10       
!!%        p(i) = 0.0D+00         
!!%        DO  j = 1, nd(i)
!!%            
!!%           k1 = k1 + 1
!!%
!!%           IF(k1.GT.nhmx)  CYCLE
!!%
!!%           p(i) = p(i) +  u(k1, k) * u(k1, k)
!!%
!!%        ENDDO
!!%
!!%        !         write(*,*) i, p(i)
!!%
!!%         sumProb10  = sumProb10  + p(i)
!!%        
!!%      ENDDO
!!%
!!%      WRITE(*,8) k, en(k)*0.5D+00,( p(i), i = 1, 6), sumProb10
!!%   ENDDO

  RETURN
END SUBROUTINE diagonalize

!#######################################################################
      
SUBROUTINE read_cfg_file(l, s, n1, l1, l2, n2_min, n2_max, idcs, ncs, ndim)
        !
        IMPLICIT NONE
        !
        INTEGER,                INTENT(in)  :: l        ! total angular momentum
        INTEGER,               INTENT(out)  :: s        ! total spin
        INTEGER, DIMENSION(:),     POINTER  :: n1, l1   !
        INTEGER, DIMENSION(:),     POINTER  :: l2, n2_min, n2_max
        INTEGER, DIMENSION(:),     POINTER  :: idcs     ! 0 or 1
        INTEGER,               INTENT(out)  :: ncs      ! nof channels included
        INTEGER,               INTENT(out)  :: ndim
        !
        INTEGER                             :: n_total_l, inpl
        INTEGER                             :: k
        !

        integer                             :: ncfg
        !
        !EXE!
        

        ncfg = 1

        !


        CALL cfile(ncfg,"inp","cfg",l)

        READ(ncfg, *) n_total_l
        READ(ncfg, *) inpl, s
        READ(ncfg, *) ncs

        WRITE(*,*) n_total_l, inpl, s, ncs

        check_partial_wave:IF( L.NE.inpl ) THEN
           WRITE(*,*) '# read_cfg: inconsistent cfg file'
           WRITE(*,*) '#                            L = ', l
           WRITE(*,*) '#                         inpl = ', inpl
           STOP
        END IF check_partial_wave
        

        
        ALLOCATE(    n1(ncs) )
        ALLOCATE(    l1(ncs) )
        ALLOCATE(    l2(ncs) )
        ALLOCATE(n2_min(ncs) )
        ALLOCATE(n2_max(ncs) )
        ALLOCATE(  idcs(ncs) )


        idcs = 0
        ndim = 0
        read_cfg:DO k = 1, ncs
           READ(ncfg,*) n1(k), l1(k), l2(k), n2_min(k), n2_max(k)!, idcs(k)
           
           ndim = ndim + n2_max(k) - n2_min(k) + 1
        ENDDO read_cfg
        

        CLOSE(ncfg)
        

        WRITE(*,*) '# read_cfg::       hamiltonian matrix dimension   ndim = ', ndim
        WRITE(*,*) '# read_cfg::        nof channels included          ncs = ', ncs        
        

        RETURN

      END SUBROUTINE read_cfg_file
!
!
!  r12 subroutines
!
      SUBROUTINE SUBMX(P1,PR,P3,PC,R,DR,VX,ISL,NHF,NLLMIN,NOLL,LHF,LL,&
     &     YK, EHF, IR, IC, LO, NP, H, IDR, RHO, ALPHA,BETA,NS,NL,NCS)

        
!!%      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
!!%      DIMENSION  PR(NS, *),   PC(NS, *),   EHF(NL, *)
!!%      DIMENSION      P1(*),       P3(*),        YK(*) 
!!%      DIMENSION       R(*),        DR(*)
!!%      DIMENSION  NHF(*), NLLMIN(*), NOLL(*), LHF(*), LL(*) 
!!%      DIMENSION  RHO(*)
!!%      DIMENSION  VX(NS, NS),   YE(NS, NPO)
!!%      DIMENSION      P2(NPO),      P4(NPO)
!!%      DIMENSION    ANGD(NL),     KAD(NL)
!!%      DIMENSION    ANGE(NL),     KAE(NL)
!!%      DIMENSION      YD(NPO)

        !
        USE PRECISION, ONLY: dpk
        !
        IMPLICIT NONE
        !

        !

        INTEGER   :: ir,ic,lo, np, idr,ns,nl
        INTEGER   :: isl,ncs
        REAL(dpk) :: h, alpha, beta
        REAL(dpk) :: vd, ve, pha
        INTEGER   :: k, n, nr, nc, nk, kk 
        INTEGER   :: ls, llr, llhfr, nhfir
        INTEGER   :: knoe,  knod


        REAL(dpk), DIMENSION(:,:) :: pr(ns,np),   pc(ns,np)
        REAL(dpk), DIMENSION(:,:) :: vx(ns, ns)
        REAL(dpk), DIMENSION(:)   :: p1(np),p3(np), p2(np),p4(np)
        REAL(dpk), DIMENSION(:)   :: ye(ns, np), yk(np), yd(np)     
        REAL(dpk), DIMENSION(:)   :: r(np), dr(np), rho(nl)
        !1e
        REAL(dpk), DIMENSION(:,:) :: ehf(nl, ns)
        INTEGER,   DIMENSION(:)   :: nhf(ncs), nllmin(ncs), noll(ncs), lhf(ncs), ll(ncs) 
        REAL(dpk), DIMENSION(:)   :: angd(nl), kad(nl), ange(nl), kae(nl)
        

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

         CALL SET_YK( P1, P3, YK, R, DR, KAD(K), np, H, IDR)

         yd(:) = yd(:) + angd(k) * yk(:) 

!         DO  J = 1, NPO
!            YD(J) = YD(J) + ANGD(K) * YK(J)
!         ENDDO

      ENDDO direct_term


!... r_k(ab,dc) : exchange term


!*   CALCULATE Y_K(J_2) FOR THE 'EXCHANGE' TERM AND SUM OVER 'K'
!*   ELECTRON '1' COORDINATES INTEGRATION

      exchange_term:DO  N = 1, NOLL(IC)
         
         ye(n,:) = 0.0           ! 1:npo
         p4(:) = pc(n,:)         ! 1:npo


         DO  K = 1, knoe
         
            CALL set_yk(p1, p4, yk, r, dr, kae(k), np, h, idr)

            ye(n,:) = ye(n,:) + ange(k) * yk(:) 

!            DO  J = 1, NPO
!               YE(N,J) = YE(N,J) + ANGE(K) * YK(J)
!            ENDDO

         ENDDO
      ENDDO exchange_term


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
            YK(1) = 0.0

            yk = yd * p2 * p4 / dr        ! 2:npo

            IF(dr(1).NE.0.0_dpk) THEN 
               yk(1) = yd(1) * p2(1) * p4(1) / dr(1)
            ENDIF


            vd = 2.0_dpk * rint( yk, 1, np, 14, h)
!'x' term


            YK(1) = 0.0_dpk

            YK(:) = YE(NC,:) * P2(:) * P3(:) / DR(:)     !2:npo

            IF(dr(1).NE.0.0) THEN 
               yk(1) = ye(nc,1) * p2(1) * p3(1) / dr(1)
            ENDIF
            

            ve = pha * 2.0_dpk * rint(yk, 1, np, 14, h) 


            !ls = 0 (singlet), 1= (triplet)

            vx(nr,nc) = vd + (1-2*ls) * ve


            ! for equivalent cfg normalize properly

         
            CALL Q2TST1(lhf, ll, nhf, nllmin, ir, ic, nr, nc, vx(nr,nc),NCS)


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
    
    !######################

      SUBROUTINE Q2TST1(LHF,LL,NHF,NLLMIN,IR,IC,NR,NC,VXS,ncs)
        !
      INTEGER,   DIMENSION(:)   :: lhf(ncs), ll(ncs), nhf(ncs), nllmin(ncs)
      INTEGER                   :: ir,ic,nr,nc,ncs
      REAL(dpk)                 :: vxs

!      IMPLICIT REAL*8(A-H,O-Z)
      !      DIMENSION LHF(1),LL(1),NHF(1),NLLMIN(1)
!.........................................................
      
      IF(NC.GT.1.AND.NR.GT.1) RETURN
      
      IF(nr.EQ.1) THEN 
         IF(nhf(ir).EQ.nllmin(ir).AND.lhf(ir).EQ.ll(ir)) THEN
            vxs = vxs /dsqrt(2.0_dpk)
         ENDIF
         
      ENDIF
      
      IF(nc.EQ.1) THEN 
         IF(nhf(ic).EQ.nllmin(ic).AND.lhf(ic).EQ.ll(ic)) THEN
            vxs = vxs/dsqrt(2.0_dpk)
         ENDIF
      ENDIF
      
      RETURN
    END SUBROUTINE Q2TST1
!
!

!!%    SUBROUTINE Q2TST1(LHF,LL,NHF,NLLMIN,IR,IC,NR,NC,VXS,NCS)
!!%      
!!%!      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
!!%      !      DIMENSION LHF(*), LL(*), NHF(*), NLLMIN(*)
!!%
!!%      INTEGER,   DIMENSION(:)   :: lhf(ncs), ll(ncs), nhf(ncs), nllmin(ncs)
!!%      INTEGER                   :: ir,ic,nr,nc,ncs
!!%      REAL(dpk)                 :: vxs
!!%
!!%!.........................................................
!!%
!!%      IF(NC.GT.1.AND.NR.GT.1) RETURN
!!%
!!%      IF(NHF(IR).EQ.NLLMIN(IR).AND.LHF(IR).EQ.LL(IR)) THEN
!!%         
!!%         IF(NR.GT.1) GOTO 10
!!%         
!!%         
!!%         VXS = VXS / DSQRT(2.0_dpk)
!!%         
!!%10    ENDIF
!!%
!!%      IF(NHF(IC).EQ.NLLMIN(IC).AND.LHF(IC).EQ.LL(IC)) THEN
!!%
!!%         IF(NC.GT.1) GOTO 11
!!%
!!%         VXS = VXS/DSQRT(2.0_dpk)
!!%
!!%11    ENDIF
!!%      
!!%      RETURN
!!%    END SUBROUTINE Q2TST1


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE WFHFIN( P, L, NPO, NHF)
        !
        !        IMPLICIT DOUBLEPRECISION(A-H,O-Z)
        ! DIMENSION P(*)
        !
        USE PRECISION, only:dpk
        !
        IMPLICIT NONE        
        !
        REAL(dpk), DIMENSION(npo) :: p
        INTEGER, INTENT(in)       :: l, npo, nhf
        !
        INTEGER                 :: k,j

        !.........
        
        CALL WF1EFILE(2, L) 
      
        DO  K = 1, NHF         
           READ(2) (P(J), J = 1, NPO)
        ENDDO
        CLOSE(2)
        
        RETURN
      END SUBROUTINE WFHFIN
      !
      !sub
      !
      SUBROUTINE WFNLIN( PP, NS, P, L, NPO, NMIN, NMX)
        
!        IMPLICIT DOUBLEPRECISION(A-H,O-Z)       
!        DIMENSION PP(NS, *), P(*)
        USE PRECISION, only:dpk
        !
        IMPLICIT NONE        
        !
        REAL(dpk), DIMENSION(ns,npo) :: pp
        REAL(dpk), DIMENSION(npo)    :: p
        INTEGER, INTENT(in)          :: ns, l, npo, nmin, nmx
        !
        INTEGER                      :: k, kmx, kmin, j

        !...........................................
        
        CALL WF1EFILE(2, L) 
        
        KMX  = NMX - NMIN + 1
        KMIN = NMIN - 1
        
        ! check this structure against the f77 version
        IF(KMIN.EQ.0) THEN 
           
           DO  K = 1, KMX
              
              READ(2) (P(J), J = 1, NPO)
              
              DO  J = 1, NPO
                 PP(K,J) = P(J)
              ENDDO
           ENDDO

        ELSE
           
           DO  K = 1, KMIN           
              READ(2) (P(J), J = 1, NPO)           
           ENDDO
        
        ENDIF

        CLOSE(2)
        
        RETURN
        !
      END SUBROUTINE WFNLIN
      !
      ! subroutine
      !
      SUBROUTINE TRNSMX(v, v_bs, ns, nhx, ihr, ihc, ir, ic)
        !exe!
        !
        USE PRECISION, ONLY:dpk
        !
        REAL(dpk), DIMENSION(ns,nhx) :: v
        REAL(dpk), DIMENSION(ns,ns)  :: v_bs
        INTEGER                      :: ns, nhx, ihr, ihc,ir,ic
        INTEGER                      :: kr,kc,kk,k, icf, irf
        !exe!
 

        IRF = IHR + IR - 1
        ICF = IHC + IC - 1
        
        DO  K = IHR, IRF
           DO  KK = IHC, ICF

              KR = K  - IHR + 1
              KC = KK - IHC + 1
              
              v(k, kk) = v_bs(kr, kc)
           ENDDO
        ENDDO
        
        RETURN
      END SUBROUTINE TRNSMX


!##########
!
!  this program calculates the integral of the function f from point na
!  to point nb using a nq points quadrature ( nq is any integer between
!  1 and 14 ).  h is the grid size.
!                                      written by c. c. j. roothaan
!

      REAL(dpk) FUNCTION rint (f,na,nb,nq,h)
        !
        !        IMPLICIT REAL*8(a-h,o-z)
        !
        
        REAL(dpk),DIMENSION(:) :: c(105),c1(25),c2(80),d(14),f(1)       
        INTEGER   na, nb, nq 
        REAL(dpk) :: h
        INTEGER l, m, i, j, n 
        REAL(dpk) :: a 
        
        EQUIVALENCE (c1(1),c(1)),(c2(1),c(26))
        DATA c1/1.d0,2.d0,1.d0,23.d0,28.d0,9.d0,25.d0,20.d0,31.d0,8.d0,&
             &1413.d0,1586.d0,1104.d0,1902.d0,475.d0,1456.d0,1333.d0,1746.d0,&
             &944.d0,1982.d0,459.d0,119585.d0,130936.d0,89437.d0,177984.d0/
             DATA c2/54851.d0,176648.d0,36799.d0,122175.d0,11108.d1,156451.d0,&
             & 46912.d0,220509.d0,29336.d0,185153.d0,35584.d0,7200319.d0,&
             & 7783754.d0,5095890.d0,12489922.d0,-1020160.d0,16263486.d0,&
             &261166.d0,11532470.d0,2082753.d0,7305728.d0,6767167.d0,9516362.d0,&
             & 1053138.d0,18554050.d0,-7084288.d0,20306238.d0,-1471442.d0,&
             & 11965622.d0,2034625.d0,952327935.d0,1021256716.d0,636547389.d0,&
             & 1942518504.d0,-1065220914.d0,3897945600.d0,-2145575886.d0,&
             & 3373884696.d0,-454944189.d0,1637546484.d0,262747265.d0,&
             & 963053825.d0,896771060.d0,1299041091.d0,-196805736.d0,&
             & 3609224754.d0,-3398609664.d0,6231334350.d0,-3812282136.d0,&
             & 4207237821.d0,-732728564.d0,1693103359.d0,257696640.d0,&
             & 5206230892907.d0,5551687979302.d0,3283609164916.d0,&
             & 12465244770050.d0,-13155015007785.d0,39022895874876.d0,&
             & -41078125154304.d0,53315213499588.d0,-32865015189975.d0,&
             & 28323664941310.d0,-5605325192308.d0,9535909891802.d0,&
             & 1382741929621.d0,5252701747968.d0,4920175305323.d0,&
             & 7268021504806.d0,-3009613761932.d0,28198302087170.d0,&
             & -41474518178601.d0,76782233435964.d0,-78837462715392.d0,&
             & 81634716670404.d0,-48598072507095.d0,34616887868158.d0,&
             & -7321658717812.d0,9821965479386.d0,1360737653653.d0/
             DATA d/2.0d0,2.0d0,24.0d0,24.0d0,1440.0d0,1440.0d0,120960.0d0,&
                  &120960.0d0,7257600.d0,7257600.d0,958003200.d0,958003200.d0,&
                  &5230697472000.0d0,5230697472000.0d0/
             
             !..............................
             
             a = 0.0d0
             l = na
             m = nb
             i = nq*(nq+1)/2
             DO  j = 1, nq
                a = a+ c(i) * (f(l)+f(m))
                l = l+1
                m = m-1
                i = i-1
             ENDDO
             a=a/d(nq)
             
             DO  n = l, m
                a = a + f(n)
             ENDDO
             rint = a*h
             RETURN
           END FUNCTION rint
          !############################################################
                !#eof
                


              END MODULE w2e_r12
!eof
!
