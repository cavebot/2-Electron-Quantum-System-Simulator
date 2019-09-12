!#
!# 
!#             
!#             transforming the f77 version to f90. In progress.... 28/May/2013
!#
!#             
!#
!#LAAN
!*######################################################################
PROGRAM r12

  USE PRECISION, ONLY: dpk
  USE param,     ONLY: np,ns
  USE set_grid
  USE w2e_r12,   ONLY: read_cfg_file
  USE bs_frb_2e, ONLY: read_target_energies
  USE wf_1e,  only: read_partial_waves_l
  implicit NONE

  INTEGER                                :: ncs, nhx             
  REAL(dpk), ALLOCATABLE,DIMENSION(:,:)  :: hmx                  ! ns, nhx
  REAL(dpk), ALLOCATABLE,DIMENSION(:,:)  :: vx                   ! ns, ns
  INTEGER,   DIMENSION(:,:), POINTER     :: ehf                  ! nl,ns
!  REAL(DPK), DIMENSION(:),   POINTER     :: en                  !

  INTEGER,   DIMENSION(:),   POINTER     :: nhf, lhf             !ncs
  INTEGER,   DIMENSION(:),   POINTER     :: ll, nllmin, nllmax   !ncs
  INTEGER,   DIMENSION(:),   POINTER     :: is, noll             !ncs
  INTEGER,   DIMENSION(:),   POINTER     :: idcs                 !ncs
  

  
  !INTEGER, ALLOCATABLE,   DIMENSION(:)   :: is, nd, n            !ncs
  !INTEGER, ALLOCATABLE,   DIMENSION(:)    :: ih                   !nhx
  REAL(DPK), DIMENSION( :),  POINTER      :: r, dr                ! r_i, dr_i
  REAL(DPK), ALLOCATABLE, DIMENSION(:,:), :: pl1, pl2, pl3, pl4   ! p_l(r), p'_l(r)
  REAL(DPK), ALLOCATABLE, DIMENSION(:),   :: r1, r2, r3, r4       ! r_i, dr_i
  REAL(dpk), ALLOCATABLE, DIMENSION(:)    :: yk                   !np
  REAL(dpk), ALLOCATABLE, DIMENSION(:)    :: pr1, pr3             !np
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)  :: pr2, pr4             !ns,np



  REAL(dpk), ALLOCATABLE, DIMENSION(:)  :: rho !nl
  REAL(dpk)                             :: alpha, beta
  !
  INTEGER                                :: ncfg, nh2eb, nv2eb

  !
  !      DIMENSION HMX(NS,NHX), VX(NS, NS),EHF(NL, NS)
  !      DIMENSION PR1(NP),PC1(NP),PR(NS,NP),PC(NS,NP),YK(NP)
  !      DIMENSION NHF(NCS),LHF(NCS),LL(NCS),NLLMIN(NCS)
  !      DIMENSION NLLMAX(NCS),NOLL(NCS),IS(NCS)
  !    DIMENSION R(NP), DR(NP), RHO(NL)
  !    DOUBLE PRECISION ALPHA, BETA

      INTEGER NINP, NOUT, NBIN
      INTEGER inpl, nl

      !
!.............

      NBIN = 1
      NINP = 9
      NCFG = 15
      NOUT = 16
      nv2eb = 3
      
      nl = 5

!.............

      
      CALL GETARG(1, ARGV)
      READ(ARGV,*) inpl
      
      WRITE(*,*) '# h2eb::        partial wave L = ', inpl 


!..............................................


...

!C#    ALPHA :: DIELECTRONIC POLARIZATION POTENTIAL: VD
!C#    BETA  :: MASS POLARIZATION POTENTIAL        : VM


      ALLOCATE(rho(nl)) 

      alpha = 0.0_dpk
      beta  = 0.0_dpk
      rho   = 0.0_dpk

!

    CALL read_target_energies(ehf)                 ! read 1e energies (in Ryd)



!#    MAKE THE GRID
!    CALL r_grid
!    CALL RIN(R, DR, H, NPO, IDR)



    lo = inpl
    CALL read_cfg_file(lo, ls, nhf,lhf,ll, nllmin, nllmax, idcs, ncs, nhx)

    ALLOCATE(   n(ncs) )
    ALLOCATE(  nd(ncs) )
    ALLOCATE(  is(ncs) )
  
    n = nhf + lhf                  !principal q.number
    nd = nllmax - nllmin + 1
  
    ncmx = 0
    DO k = 1, ncs
       is(k) = ncmx + 1
       ncmx  = ncmx + nd(k)
    ENDDO

!.................................................


    WRITE(*, *) lo, ncs, ncmx,nhx
    
    ALLOCATE(yk(np))
    !
    ALLOCATE(pl1(ns,np))
    ALLOCATE(pl2(ns,np))
    ALLOCATE(pr1(np))
    ALLOCATE(pr2(ns,np))
    !
    ALLOCATE(pl3(ns,np))
    ALLOCATE(pl4(ns,np))
    ALLOCATE(pr3(np))
    ALLOCATE(pr4(ns,np))

    

    ALLOCATE( hmx(ns,nhx) )
    ALLOCATE( vx(ns,ns) )




    
    CALL r12file(nv2eb, lo)         ! read r12 data


    
    kr = 1
    outer_loop: DO  ir = 1,  ncs

       llr   =  ll(ir) + 1
       llhfr = lhf(ir) + 1
       nhfir = nhf(ir)
       

!#     N1,    L1,  L2, N2_MIN, N2_MAX
!#     NHFIR, LHF, LL, NLLMIN, NLLMAX


!#  STORE P_{ N1, L1}( R_J )    R_J : GRID POINTS
!#  STORE P_{N2,L2}(R_J)    N2_MIN < N_2 < N2_MAX, R_J : GRID POINTS

      call read_all_wf1e(lhf(ir), pl1, r1)
      call read_all_wf1e(ll(ir),  pl2, r2)         


      pr1(:) = pl1( nhf(ir), :)                            
      DO k = 1, nd(ir)  
         pr2(k,:) = pl2(nllmin(ir) - 1 + k, : )  
      ENDDO
      yk(:) = pl2(nllmax(ir), : )  

!? pr(1:nd(ir),:) = pl2(nllmin(ir):nllmax(ir), : )




!      CALL WFHFIN( PR1, LHF(IR), NP, NHFIR)
!      CALL WFNLIN( PR, NS, YK, LL(IR), NP, NLLMIN(IR), NLLMAX(IR))

      kc = kr

!    INTEGER llr, llhfr, nhfir,ir,ic, kr, kc, llc, llhfc

      hmx = 0.0_dpk
      inner_loop: DO  ic = ir, ncs


         llc   = ll(ic) + 1
         llhfc = lhf(ic) + 1



         CALL read_all_wf1e(lhf(ic), pl3, r1)
         CALL read_all_wf1e(ll(ic),  pl4, r2)         


         pr3(:) = pl3( nhf(ic), :)                            
         DO k = 1, nd(ic)  
            pr4(k,:) = pl4(nllmin(ic) - 1 + k, : )  
         ENDDO
         yk(:) = pl4(nllmax(ic), : )  


!         CALL WFHFIN( PC1, LHF(IC), NP, NHF(IC))
!         CALL WFNLIN( PC, NS, YK, LL(IC), NP, NLLMIN(IC), NLLMAX(IC))


         WRITE(*,*) '  ic = ',  ic

         CALL SUBMX(PR1, PR2, PR3, PR4,R,DR,VX,LS,NHF,NLLMIN,NOLL,LHF,LL, &
              &        YK, EHF, IR, IC, LO, NP, H, IDR,                &
              &        RHO, ALPHA, 0.0_dpk, ns, nl,ncs)



         CALL TRNSMX(hmx, vx, ns, nhx, 1, kc, noll(ir), noll(ic))


     kc = kc + noll(ic)

  ENDDO inner_loop
      
      kr = kr - 1
      save_r12_matrix_elements: DO  nr = 1, noll(ir)

         kr = kr + 1

         WRITE(nv2eb) ( hmx( nr, nc), nc = kr, ncmx)

      ENDDO save_r12_matrix_elements
      
      kr = kr + 1

   ENDDO outer_loop
   
   CLOSE(nv2eb)


 END PROGRAM r12
 !#######################################################################
 !*EOF
