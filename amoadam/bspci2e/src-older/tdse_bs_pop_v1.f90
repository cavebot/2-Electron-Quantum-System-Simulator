!
! version where it was assumed as icll to include e.g sp and ps case. Nevertheless this is obsolete since
! in configuration files only (sp) appears and not (ps) again. So the case where (k1l2;k2l1) is useless essentially.
!
!
! [1] PRA 68,013409, Laulan and Bachau,
! 'Correlation effects in two-photon single and Double ionization of helium'
!

!
! calculates ionization probabilities for the 2-e Hamiltonian problem
!
!
!  id/dt psi(r1,r2;t) = [ H(r1,r2) + V(r1,r2;t)] psi(r1,r2;t) 
!
!    H = h(r1) + h(r2) + 1/|r_1-r_2| = H_0 + V12
!
!    if |i> defined as the CI 2e states (E;LS)  (ignore here M_L,M_S)
!
!    H|i> = E |i>,  S_i |i><i| = 1 and <i|i'> = delta_ii' 
!
!    and |j> the zero-order 2e states (uncorrelated)    (n1l1;n2l2;E_0=e1+e2;LS)
!
!          H_0|j> = E_0 |j>,  S_i |j><j| = 1 and <j|j'> = delta_jj', E_0 = e_1 + e_2
!
!     with e1,e2 the eigenenergies of the 1-e SE:
!
!         h(r1)|e1> = e1 |e1> and h(r2) |e2> = e2|e2>
!
! then this program calculates two different ionization probabilities are calculated
!
! - 
!
!  (1) psi(t) = S_i c_i(t) |i>,   TDSE calculation provides c_i(t) due to the laser interaction
!x
!  with    i == (nL) 
!
!  |i> = S_{j} v_ij |j>,           CI calculation provides elements v_ij of the V12 interaction
!
!    j  = n1l1;n2l2 and Phi_j the zero-order  
! v_ij  = the  CI-mixing coefficients
! 
!  |j > = the zero-order (uncorrelated) 2e-hamiltonian
!
! we have two expansions for the TD wf:
!
!  psi(t) = S_i c_i(t) |i> = S_j v_j(t) |j>,   v_j = S_i c_i(t) v_ij
!
!
! then
!       P_j(t) = |<j|psi(t)>| = ...= |S_i c_i(t)v_ij|^2
!
!  probability that the system will be found at state |j>
!
!  with energy E= e_1 + e_2, 
!
! and
!
!      P_i(t)  = |<i|psi(t)>| = ...= |c_i(t)|^2
!
! probability that the system will be found in state |i>
!
! with energy E_n and angular momentum L.
!
!
!
!
!
PROGRAM tdse_pop_sdi
  !
  USE PRECISION
  USE parameter_tdse_fxd
  USE bs_frb_2e,      ONLY: read_target_energies
  USE io
  USE atom_2e!, ONLY: read_w2e_ci, read_w2e_ct
  !  USE atom_1e!, ONLY: read_w1e  ! to be implemented instead of read_target_energies

  IMPLICIT NONE
  !
  !
  INTEGER                                     :: k1, k2
  INTEGER                                     :: ie, ie_min, ie_max
  INTEGER                                     :: ic, icl, icll, ncl
  INTEGER                                     :: l,  le1,le2
  INTEGER                                     :: ll, lle1,lle2
  INTEGER                                     :: nll_max
  INTEGER,   ALLOCATABLE, DIMENSION(:,:)      :: chl               ! chl(0:lmax, 1:nll_max)
  INTEGER                                     :: n_ll
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: l1, l2            ! l1(1:nll), l2(1:n_ll)
  ! coefficients
  !field! 

  !
  !grid!

  !field-free atomic configurations
  TYPE(symmetry_ls2e), ALLOCATABLE, DIMENSION(:) :: w2e
  INTEGER                                        :: l1e_max
  REAL(dpk),   DIMENSION(:,:), POINTER           :: en1e                !nl,ns
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)         :: e1e                 !nl,ns
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)         :: rho_e1e             !nl,ns
  INTEGER                                        :: ne_1e, nl_1e
  !(summed over all 2e L = l1 + l2)
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:)    :: popl_e1e2          ! (n_ll, e1,e2)  
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: popl_e1, popl_e2   ! (n_ll,e1), (n_ll,e2)
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: pop_e1e2           ! pop_e1e2(e1,e2)
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_e1, pop_e2      ! pop_e1e2(e1,e2)
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_ll_e1, pop_ll_e2      ! pop_e1e2(e1,e2)
  !
  REAL(dpk)                                   :: sum_this, sum_this_1, sum_this_2  
  REAL(dpk)                                   :: pe, pop_icl
  COMPLEX(dpk)                                :: zsum
  !
  REAL(dpk)                                   :: pop_di_1, pop_di_2
  REAL(dpk)                                   :: e1,e2               !e1 = e1e(l1,n1), e2 = e1e(l2,n2)
  !
  INTEGER                                     :: ne_1e_cut       !ne_1e_cut = ne_1e - nof_states_to_exclude
  !
  REAL(dpk)                                   :: norm_0, norm_v  ! all population
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: popl_0, popl_v  !population in L-channel
  REAL(dpk)                                   :: pop_b_0, pop_s_0, pop_d_0
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: pop_Lll_0 !
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: popl_b_v,  popl_b_0 !
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: popl_s_v,  popl_s_0 
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: popl_d_v,  popl_d_0
  !var
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: ndi_l           ! pop with  E(ndi_l,l) = E++
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nsi_l           ! pop with  E(ndi_l,l) = E+
  INTEGER                                     :: n_l2e           ! nof CI states of symmetry l2e
  !I/O
  CHARACTER(LEN=25)                           :: pop_file         !i
  CHARACTER(LEN=25)                           :: ion_file         !o
  CHARACTER(len = 60 )                        :: data_file
  CHARACTER(len = 6  )                        :: sl,sl1
  LOGICAL                                     :: check_zero
  INTEGER                                     :: n_cols
  !
  ! now read
  !
  
  pop_file = "pes/pes_"
  ion_file = 'ion/pop_di.dat'
  !  coe_file = 'tdat/coe.dat'

  !
  CALL input_tdse_fxd              ! tinp/tdse_bs_fxd_2e.inp
  CALL output_tdse_fxd             ! tout/tdse_bs_fxd_2e.out
  CALL read_target_energies(en1e)  ! read 1e energies (in Ryd)
  CALL read_w2e_ct(w2e,lmax)       ! read TDSE wavefunction $\Psi(t) = \sum_i c_{i}(t) \Phi_{i}$
  DO l = 0, lmax
     CALL read_w2e_ci(w2e(l),l)    !> read coefs $\Phi_j = \sum_i v_{ij} \Phi^{(0)}_{i}$
  ENDDO

  
  
  !
  WRITE(*,'(a1)')'&'
  WRITE(*,'(a1,a60,G15.3)')'&',' input w1e, w2e read finished ' 
  WRITE(*,'(a1)')'&'
  WRITE(*,'(a1)')'&'

  !  
  ! some checks
  !  filetype = "cdf"
  !
        
     DO l = 0, lmax
        WRITE(*,'(a1,a60)') '&','----------------------------------'
        WRITE(*,'(a1,a60,I6)') '&','        ang.  symmetry      l2e = ', w2e(l)%l
        WRITE(*,'(a1,a60,I6)') '&','        spin  symmetry      s2e = ', w2e(l)%s
        WRITE(*,'(a1,a60,I6)') '&',' nof zero-2e-states         n2e = ', w2e(l)%nev
        WRITE(*,'(a1,a60,I6)') '&',' nof config. series         ncs = ', w2e(l)%nch
        WRITE(*,'(a1,a60)'   ) '&',' read configuration finished.'
        !!%        
     ENDDO

     !
     !  data processing now


     !
     ! work out the 1-e energies
     !          
     ! l1e_max : is the max angular channel of the 1e-states included in the cfg
     ! nl_1e-1 : the nof angular channels calculated for the 1e-states
     ! In a full basis l1e_max =lmax_1e but in a truncated basis might differ. 
     !
     !

     !
     !NOTE::    en1e pointer index from 1:nl_1e NOT 0:nl_1e
     !
     
     en1e = en1e * 0.5_dpk              ! convert energies in a.u.
     
     
     nl_1e = SIZE(en1e,dim=1)  ! e1e(l1_max, n1_max) 
     ne_1e = SIZE(en1e,dim=2)  ! e1e(l1_max, n1_max)

     
     
     ! for conformity reasons we rewrite
     ! en1e(1:nl_1e,1:ne_1e) --> e1e(0:nl_1e-1,1:ne_1e) -->



     !First find the largest angular momentum included in the cfg channels.

     l1e_max = 0
     DO l = 0, lmax   
        IF(SIZE(w2e(l)%l12,dim=1).GT.l1e_max) l1e_max = SIZE(w2e(l)%l12,dim=1)
     ENDDO
     l1e_max = l1e_max - 1   ! since l1e_max inside the loop returns the size of l12
     
     PRINT*,"&             l1e_max = ", l1e_max
     
     IF((nl_1e-1).LT.l1e_max) THEN
        PRINT*,"something wrong here:"
        PRINT*," max l1e calculated              ne_1e = ", nl_1e-1
        PRINT*," max l1e included in cfg files l1e_max = ", l1e_max
        STOP
     ENDIF

     
     ALLOCATE( e1e(0:l1e_max,1:ne_1e))

     
     DO l = 0, l1e_max
        DO ie = 1, ne_1e
           e1e(l,ie) = en1e(l+1,ie)
        ENDDO
     ENDDO
     DEALLOCATE(en1e)
     

  ALLOCATE(rho_e1e(0:l1e_max, 1:ne_1e) )  ! e1e(l1_max, n1_max)

  energy_density_1e:DO l = 0, l1e_max
     rho_e1e(l,1) =  1.0_dpk / ABS(e1e(l,2)-e1e(l,1))
     DO ie = 2, ne_1e-1
        rho_e1e(l,ie) = 2.0_dpk / ABS(e1e(l,ie-1)-e1e(l,ie+1))
     ENDDO
     rho_e1e(l,ne_1e) = 1.0_dpk / ABS(e1e(l,ne_1e-1)-e1e(l,ne_1e))     
  ENDDO energy_density_1e


  check_zero = ALL(rho_e1e > 0.0_dpk)
  IF(check_zero.EQV.(.FALSE.)) THEN
     PRINT*, "# tdse_bs_pop:: rho_e1e contains zero values. Stop."
     STOP
  ENDIF


  
  !
  !
  ! 2-electron work out now
  !
  !


     
     DO l = 0, lmax
        w2e(l)%e2e = 0.5_dpk * w2e(l)%e2e    !transform to a.u. ?  
     ENDDO

     
     ALLOCATE( ndi_l(0:lmax) )
     ALLOCATE( nsi_l(0:lmax) )


     ndi_l =  0
     nsi_l  =  0
     find_si_di_thresholds:DO l = 0, lmax
        
        n_l2e = w2e(l)%net

        si:DO  ie = 1, n_l2e     ! find SI threshold energy  for symmetry L        
           nsi_l(l) = ie
           !        IF(  w2e(l)%e2e(ie) > 0.0_dpk ) EXIT
           IF(  w2e(l)%e2e(ie) > en_ion_1 ) EXIT        
        ENDDO si
        
     
     di:DO  ie = 1, n_l2e     ! find DI threshold energy  for symmetry L       
        ndi_l(l) = ie       

        ! since 2-electron systemd by definition en_ion_2 = 0.0_dpk
        
        IF( ( w2e(l)%e2e(ie) ) > 0.0_dpk ) EXIT            
     ENDDO di
     WRITE(*,*) "tdse_sdi_pop:: si threshold              en_ion_1 = ", en_ion_1
     WRITE(*,*) "tdse_sdi_pop:: si threshold  E(", nsi_l(l), l, ") = ", w2e(l)%e2e(nsi_l(l))
     WRITE(*,*) "tdse_sdi_pop:: di threshold  E(", ndi_l(l), l, ") = ", w2e(l)%e2e(ndi_l(l))
     WRITE(*,*) "tdse_sdi_pop:: max energy    E(", n_l2e, l,    ") = ", w2e(l)%e2e(n_l2e)

     !     IF(nsi_l(l).EQ.1) THEN
     !        WRITE(*,*) "tdse_sdi_pop:: ie for SI threshold can't be ", nsi_l(l)
     !        STOP
     !     ENDIF
     
  ENDDO find_si_di_thresholds

  


  ! population based on  Eq.(11) of Ref [1]   (Laulan and Bachau)
  

  ALLOCATE(popl_0(0:lmax) )
  
  popl_0 = 0.0_dpk
  partial_probability_Lv:DO  l = 0, lmax

     PRINT*, '#          l2e = ',l
     PRINT*, '#          ncf = ',w2e(l)%ncf


     
     !     ie_min = nsi_l(l)
     ie_min = ndi_l(l) 
     WRITE(*,*) "tdse_sdi_pop::       ie_min(",l,") = ", ie_min
     

    ! ie_min = 1    ;
    ! IF(l.EQ.0) ie_min = 2   ! do not include population of GROUND state (L=0)
                             ! exclude from projection (see Eq.(11) of Ref [1])


     ALLOCATE( w2e(l)%pvt_0(1:w2e(l)%ncf ))
     
     
     ie_max = w2e(l)%net             
     population_of_cfg_a:  DO ic = 1, w2e(l)%ncf  ! cfg index  
        
        zsum = 0.0_dpk
           DO ie = ie_min, ie_max                              ! j = (n1l1;n2l2;LS)
              zsum = zsum + w2e(l)%cv(ie,ic) * w2e(l)%ct(ie)   !v_j(t) = S_i c_i(t)v_ij
           ENDDO
        
        w2e(l)%pvt_0(ic) = ABS(zsum)**2                  !p_j(t) = |v_j(t)|^2  
        
     ENDDO population_of_cfg_a
        
     popl_0(l)  = SUM( w2e(l)%pvt_0 )   ! sum over all cfgs (n1l1;n2l2) ionization in L channel

  ENDDO partial_probability_Lv
  
  !
  PRINT*, "# proceed with partial wave probabilities now"
  !
  
  !REMARKS:
  !
   ! probabilities in terms of the zero 2e-functions
  !
  ! \[
  !   \psi({\bf r}_1,{\bf r}_2;t) = \sum_{L;nln'l'} v_{L;nln'l'}(t) \Phi_{L;nln'l'}({\bf r}_1,{\bf r}_2)
  ! \]
  !with
  ! \[
  !  v_{L;nln'l'}(t) = \sum_{N} C_{NL}(t) 
  ! \]  
  !


  !
  !
  !  probabilities in terms of  (E = e1 + e2) and (n1,l1, n2,l2)
  !
  !

  nll_max = 0 
  DO l = 0, lmax
     nll_max = nll_max + w2e(l)%chl_n 
  ENDDO

  
  ALLOCATE( pop_Lll_0(0:lmax,1:nll_max) )
  ALLOCATE(   popl_b_0(0:lmax) )          !bound
  ALLOCATE(   popl_s_0(0:lmax) )          !single ionization
  ALLOCATE(   popl_d_0(0:lmax) )          !double ionization
  !
  pop_Lll_0 = 0.0_dpk 
  popl_b_0  = 0.0_dpk
  popl_s_0  = 0.0_dpk
  popl_d_0  = 0.0_dpk
  !
  !
  angular_symmetry_2e_zero_order: DO l = 0, lmax
     
     PRINT*, "l2e = ", l     


     ALLOCATE( w2e(l)%popd_ee_Lll( w2e(l)%chl_n, 1:ne_1e, 1:ne_1e) )



     w2e(l)%popd_ee_Lll = 0.0_dpk ;
     pop_b_0  = 0.0_dpk
     pop_s_0  = 0.0_dpk
     pop_d_0  = 0.0_dpk
     sum_all_configurations: DO ic = 1, w2e(l)%ncf  ! based on (L,ic) --->  e1,l1 and e2,l2

        !ic = (k1 l1 ; k2 l2)
        
        k1  = w2e(l)%n1(ic)    
        le1 = w2e(l)%l1(ic)
        e1  = e1e(le1, k1 )
        
        k2  = w2e(l)%n2(ic)  
        le2 = w2e(l)%l2(ic)
        e2  = e1e(le2, k2 )

        !note:
        ! according Eq. (11) of ref [1]
        ! population of (L,ic) from all states excluding ONLY the GROUND state

        !
        ! pop(L,e1,e2) = |S_i c_i(t) * v_ij|^2
        !
        
        ! the below is associated with   rho_e1e(le1,k1) * rho_e1e(le1,k2)  density of states

        
        pe   = w2e(l)%pvt_0(ic)       ! population p(L:e_1,e2) = |S_i c_i(t)v_ij|^2  

        
        ! ic <--> (e1 l1 ; e2 l2)
        pick_di_channels:IF( (e1>0).AND.(e2>0) )  THEN

           ! with the above the pop(ich,k1,k2) = pop_Ll1l2(k1,k2) EQ  pop_Ll2l1(k2,k1) 
           !              but   pop_Ll1l2(k1,k2) NEQ pop_Ll1l2(l2,k1)

           !
           ! k1 is associated with chl_1(icl) 
           ! k2 is associated with chl_2(icl) 
           !
           ! this means that for pop(icl,k1,k2) was multiplied by rho_(chl_1(icl),e1) 
           !
           
           ! find the proper icl for the (le1,le2)
           pick_l12_channel:DO icl = 1, w2e(l)%chl_n       !channels (l1,l2) of L symmetry

              pop_icl = rho_e1e(w2e(l)%chl_1(icl), k1) * rho_e1e(w2e(l)%chl_2(icl), k2) * pe
              
              IF( (w2e(l)%chl_1(icl).EQ.le1).AND.(w2e(l)%chl_2(icl).EQ.le2) ) THEN

                 !                 w2e(l)%popd_ee_Lll(icl,k1,k2) = rho_e1e(le1, k1) * rho_e1e(le2, k2) * pe

                 w2e(l)%popd_ee_Lll(icl,k1,k2) = pop_icl
                 
                 EXIT pick_l12_channel                 

              ELSE IF( (w2e(l)%chl_1(icl).EQ.le2).AND.(w2e(l)%chl_2(icl).EQ.le1) ) THEN

                 !                 w2e(l)%popd_ee_Lll(icl,k2,k1) = rho_e1e(le2, k1) * rho_e1e(le1, k2) * pe
                  
                 w2e(l)%popd_ee_Lll(icl,k2,k1) = pop_icl
                 
                 EXIT pick_l12_channel                  
              ELSE
                 CYCLE
              ENDIF
              
              ! should never reach this point if chl_1,chl_2 were detected correctly
              PRINT*, " something wrong (l, icl) = ", l, icl
              STOP
              
           ENDDO pick_l12_channel
           !
           !P_Ll1l2(k1,k2)           
           !           w2e(l)%popd_ee_Lll(icl,k1,k2) = rho_e1e(w2e(l)%chl_1(icl), k1) * rho_e1e(w2e(l)%chl_2(icl), k2) * pe
           
           !           w2e(l)%popd_ee_Lll(icl,k1,k2) = rho_e1e(le1, k1) * rho_e1e(le2, k2) * pe
           !           w2e(l)%popd_ee_Lll(icl,k2,k1) = w2e(l)%popd_ee_Lll(icl,k1,k2)
           
           pop_d_0 = pop_d_0 + pe
           pop_Lll_0(l,icl) =  pop_Lll_0(l,icl) + pe 
           
        ENDIF pick_di_channels
        
           
        si_channels:IF ( (e1<0).AND.(e2<0) ) THEN
           pop_b_0 = pop_b_0 + pe   ! (e1,e2 < 0)              
        ELSE
           pop_s_0 = pop_s_0 + pe    !  (e1<0 and e2>0 or e2 < 0 and e1 > 0)
        ENDIF si_channels

        
     ENDDO sum_all_configurations
     !     STOP      
     popl_b_0(l)   = pop_b_0  ! ionization in L-channel with e1<0 and e2<0     
     popl_s_0(l)   = pop_s_0  ! ionization in L-channel with e1<0 and e2>0 or vice versa
     popl_d_0(l)   = pop_d_0  ! ionization in L-channel with e1>0 and e2>0
     !
  ENDDO angular_symmetry_2e_zero_order

  norm_0 = SUM(popl_b_0 + popl_s_0 + popl_d_0)     ! all population in L-channgel bound + single + double
  ! sum over L-channels  (normalization to 1)



  
  
  !
  !  pop_d_e12_channels_0 = dP(e1,e2)/de1de2            (for fixed l,le1,le2)
  !
  !
  !-  Marginal probability for e1  
  !  pop_d_e1_channels_0 = I_de2 [ dP(e1,e2)/de1de2]    (for fixed l,le1,le2)
  !
  !
  !-  Marginal probability for e2  
  !  pop_d_e2_channels_0 = I_de1 [ dP(e1,e2)/de1de2]    (for fixed l,le1,le2)
  !
  !
  !
  !  pop_d_channels_0 = I_de_1_de2 [ dP(e1,e2)/de1de2]  (for fixed l,le1,le2)
  !
  !
  !  S_{l,le1,le2} = pop_double_ionization (e1>0, e2>)
  !
  

  ne_1e_cut = ne_1e - 11

  
  !
  ! now calculate P_{L,l1,l2}(k_1),
  !               P_{Ll1l2} = Sum_{k_1} { P_{L,l1,l2}(k_1) },
  !               P_L       = Sum_{l1l2}{P_Ll1l2} 
  !               P_DI      = Sum_L{P_L}


  !
  !  PLll'(ek_1)
  !
  pop_1D_e1:DO l = 0, lmax
     
     ALLOCATE(  w2e(l)%popd_e1_Lll( w2e(l)%chl_n, 1:ne_1e ) )
     ALLOCATE(     w2e(l)%popd_Lll_1( w2e(l)%chl_n ) )
     
     w2e(l)%popd_e1_Lll = 0.0_dpk ;
     w2e(l)%popd_Lll_1  = 0.0_dpk ;
     
     pop_Lll_e1:DO icl = 1, w2e(l)%chl_n
        le2 = w2e(l)%chl_2(icl)
        el1:DO k1 = 1, ne_1e_cut
           sum_this_1 = 0.0_dpk
           sum_over_e2:DO k2 = 1, ne_1e_cut
              sum_this_1 = sum_this_1 + w2e(l)%popd_ee_Lll(icl,k1,k2)/rho_e1e(le2,k2)
           ENDDO sum_over_e2
           w2e(l)%popd_e1_Lll(icl,k1) = sum_this_1
        ENDDO el1
        
     ENDDO pop_Lll_e1
  ENDDO pop_1D_e1  
  
  

  ! Calculate P_{L,l1,l2}(e1) by integrating over e2
  
  pop_e1_Lll:DO l = 0, lmax
     
     DO icl = 1, w2e(l)%chl_n
        le1 = w2e(l)%chl_1(icl)        
        sum_this = 0.0_dpk
        sum_over_ee1:DO k1 = 1, ne_1e_cut
           sum_this = sum_this + w2e(l)%popd_e1_Lll(icl, k1)/rho_e1e(le1,k1) 
        ENDDO sum_over_ee1
        w2e(l)%popd_Lll_1(icl) = sum_this
     ENDDO
     
     w2e(l)%popd_L_1 = SUM(w2e(l)%popd_Lll_1)            ! now calculate P_{L}   
  ENDDO pop_e1_Lll

  
  ! Calculate DI


  pop_di_1 = SUM(w2e(0:lmax)%popd_L_1)

  !
  !
  ! repeat the same for electron  '2'
  !
  ! now calculate P_{L,l1,l2}(k_2),
  !               P_{Ll1l2} = Sum_{k_2} { P_{L,l1,l2}(k_2) },
  !               P_L       = Sum_{l1l2}{P_Ll1l2} 
  !               P_DI      = Sum_L{P_L}
  !

  
  ! PLll'(k2)
  !
  pop_1D_e2:DO l = 0, lmax
     ALLOCATE(  w2e(l)%popd_e2_Lll( w2e(l)%chl_n, 1:ne_1e ) )
     ALLOCATE(     w2e(l)%popd_Lll_2( w2e(l)%chl_n ) )
     w2e(l)%popd_e2_Lll = 0.0_dpk ;
     w2e(l)%popd_Lll_2  = 0.0_dpk ;
     pop_Lll_e2:DO icl = 1, w2e(l)%chl_n
        le1 = w2e(l)%chl_1(icl)
        el2:DO k2 = 1, ne_1e_cut
           sum_this_2 = 0.0_dpk
           sum_over_e1:DO k1 = 1, ne_1e_cut   !sum over e1
              sum_this_2 = sum_this_2 + w2e(l)%popd_ee_Lll(icl,k1,k2)/rho_e1e(le1,k1)
           ENDDO sum_over_e1
           w2e(l)%popd_e2_Lll(icl,k2) = sum_this_2
        ENDDO el2
     ENDDO pop_Lll_e2
  ENDDO pop_1D_e2
  

  
  ! Calculate P_{L,l1,l2}(e1) by integrating over e2


  pop_e2_Lll:DO l = 0, lmax   
     DO icl = 1, w2e(l)%chl_n
        le2 = w2e(l)%chl_2(icl)
        sum_this = 0.0_dpk
        DO k2 = 1, ne_1e_cut
           sum_this = sum_this + w2e(l)%popd_e2_Lll(icl, k2)/rho_e1e(le2,k2) 
        ENDDO
        w2e(l)%popd_Lll_2(icl) = sum_this
     ENDDO
     w2e(l)%popd_L_2 = SUM(w2e(l)%popd_Lll_2)            ! now calculate P_{L}
  ENDDO pop_e2_Lll


  ! Calculate DI 
  pop_di_2 = SUM(w2e(0:lmax)%popd_L_2)


  
  !
  ! now sum over l channels
  !


  nll_max = 0 
  DO l = 0, lmax
     nll_max = nll_max + w2e(l)%chl_n 
  ENDDO
  
  ALLOCATE( chl(0:lmax, 1:nll_max) )

  chl = 0
  ncl = 1  
  GATHER_L1L2_FROM_ALL_L_SYMMETRIES:DO l = 0, lmax     
     
     IF(l.EQ.0) THEN
        
        DO icl = 1, w2e(l)%chl_n

           chl(l,icl) = ncl                ! w2e(l)%chl_2(icl)           
            
           ncl = ncl + 1
           
        ENDDO
           
     ELSE  ! now find common (l1,l2) in other symmetries
        
        channel:DO icl = 1, w2e(l)%chl_n
            
            le1 = w2e(l)%chl_1(icl)
            le2 = w2e(l)%chl_2(icl)
            
            check_for_duplications_on_all_ls:DO ll = 0, l-1
               check_for_duplications_inside_ll:DO icll = 1, w2e(ll)%chl_n
                  
                  lle1 = w2e(ll)%chl_1(icll)
                  lle2 = w2e(ll)%chl_2(icll)
                  
                  IF( (le1==lle1).AND.(le2==lle2)) THEN
                     chl(l,icl) = chl(ll,icll)
                     CYCLE channel
                  ENDIF
                  
                  IF( (le1==lle2).AND.(le2==lle1)) THEN
                     chl(l,icl) = chl(ll,icll)
                     CYCLE channel
                  ENDIF
                  
               ENDDO check_for_duplications_inside_ll
            ENDDO check_for_duplications_on_all_ls

            ! only if no other channel same as the current (l1,l2) was found to the other
            ! symmetries register as a NEW one
            
            chl(l,icl) = ncl                ! w2e(l)%chl_2(icl) 
            ncl = ncl + 1
         ENDDO channel
         
      ENDIF
      !
        
   ENDDO GATHER_L1L2_FROM_ALL_L_SYMMETRIES

   n_ll = MAXVAL(chl)
   
   ALLOCATE( l1( MAXVAL(chl) ))
   ALLOCATE( l2( MAXVAL(chl) ))
   
   DO l = 0, lmax
      DO icl = 1, w2e(l)%chl_n
         l1( chl(l,icl) ) = w2e(l)%chl_1(icl)
         l2( chl(l,icl) ) = w2e(l)%chl_2(icl)
         PRINT*, l, icl, chl(l,icl), l1(chl(l,icl)), l2(chl(l,icl)) 
      ENDDO
   ENDDO


   
!   PRINT*, " total number of partial channels in all symmetries: ", n_ll
!   DO icl = 1, n_ll
!      PRINT*, 'channel ', icl, ' --> (', l1(icl),l2(icl), ' )'
!   ENDDO


     
  ALLOCATE(popl_e1e2( 1:n_ll,  1:ne_1e, 1:ne_1e) )
  ALLOCATE(  popl_e1( 1:n_ll,  1:ne_1e) )
  ALLOCATE(  popl_e2( 1:n_ll,  1:ne_1e) )
  ALLOCATE( pop_e1e2( 1:ne_1e, 1:ne_1e) )
  ALLOCATE(   pop_e1( 1:ne_1e) )
  ALLOCATE(   pop_e2( 1:ne_1e) )
  !


  ! rho_l1l2(e1,e2)                    = popl_e1e2(icll,k1,k2)
  ! rho_l1l2(e1), rho_l1l2(e2)         = popl_e1(icll,k1), popl_e2(icll,k2)
  
  DO icll = 1, n_ll
     
     popl_e1e2(icll,:,:) = 0.0_dpk
     popl_e1(icll,:)     = 0.0_dpk
     popl_e2(icll,:)     = 0.0_dpk
     DO l = 0, lmax
        DO icl = 1, w2e(l)%chl_n
            
            IF( chl(l,icl) == icll)  THEN
               popl_e1e2(icll,:,:) = popl_e1e2(icll,:,:) +  w2e(l)%popd_ee_Lll(icl,:,:)
               popl_e1(icll,:)     = popl_e1(icll,:)     +  w2e(l)%popd_e1_Lll(icl,:)
               popl_e2(icll,:)     = popl_e2(icll,:)     +  w2e(l)%popd_e2_Lll(icl,:)
            ENDIF
            
         ENDDO
      ENDDO  
   ENDDO
   
   !
   ! now summed over all partial channels (l1,l2) --> P(e1,e2), P(e1), P(e2)
   ! (requires interpolation in the energies (e1,e2)
   !
   
   pop_e1e2 = SUM(popl_e1e2,dim = 1)     ! rho(e1,e2)   
   pop_e1   = SUM(popl_e1,  dim = 1)     ! rho(e1)
   pop_e2   = SUM(popl_e2,  dim = 1)     ! rho(e2)


   ALLOCATE(pop_ll_e1(1:n_ll),pop_ll_e2(1:n_ll))


   DO icll = 1, n_ll 
      le1 = l1(icll)
      sum_this = 0.0_dpk
      DO k1 = 1, ne_1e_cut
         sum_this = sum_this + popl_e1(icll,k1)/rho_e1e(le1,k1) 
      ENDDO
      pop_ll_e1(icll) = sum_this
   ENDDO
   !
   DO icll = 1, n_ll 
      le2 = l2(icll)
      sum_this = 0.0_dpk
      DO k2 = 1, ne_1e_cut
         sum_this = sum_this + popl_e2(icll,k2)/rho_e1e(le2,k2) 
      ENDDO
      pop_ll_e2(icll) = sum_this
   ENDDO

   !
   ! STORE RESULTS
   !
   !
   ! w2e(l)%popd_ee_Lll( icl, k1, k2)
   ! w2e(l)%popd_e1_Lll( icl, k1    ) 
   ! w2e(l)%popd_Lll_1 ( icl        ) 
   !



   !
   ! store info for the various (l1,l2) channels in L 
   !

        
   !   data_file = TRIM(pop_file)//"ll_yld.dat"
   data_file = "pes/yield_ll.dat"
   OPEN(nout,file=data_file)
   WRITE(nout,'(a60)') "..................................................................."
   WRITE(nout,'(1a5,11X,2a20)') "nll","total DI(e1)","total DI(e2)"
   WRITE(nout,'(3a5,1X,2a20)') "l1","l2","icll","pop_e1_ch", "pop_e2_ch"
   WRITE(nout,'(a60)') "..................................................................."
   WRITE(nout,'(1I5,11X,2ES20.8)') n_ll, SUM(pop_ll_e1), SUM(pop_ll_e2)
   DO icll = 1, n_ll
      WRITE(nout,'(3I5,1X, 2ES20.8)') l1(icll), l2(icll), icll, pop_ll_e1(icll), pop_ll_e2(icll)
   ENDDO
   WRITE(nout,'(a60)') "..................................................................."
   WRITE(nout,'(a40,a20)') "total DI(e1)","total DI(e2)"
   WRITE(nout,'(2a5,10X,2a20)') "l,","n_l", "pop_L_e1","pop_L_e2"
   WRITE(nout,'(4a5,2a20)')      "l1","l2","chl","icl", "pop_Lll_e1","pop_Lll_e2"
   WRITE(nout,'(a60)') "...................................................................."
   WRITE(nout,'(20X,2ES20.8)') pop_di_1, pop_di_2 
   DO l = 0, lmax
      WRITE(nout,'(2I5,10X,2ES20.8)') l, w2e(l)%chl_n,  w2e(l)%popd_L_1,   w2e(l)%popd_L_2
      DO icl = 1, w2e(l)%chl_n
         l1( chl(l,icl) ) = w2e(l)%chl_1(icl)
         l2( chl(l,icl) ) = w2e(l)%chl_2(icl)
         WRITE(nout,'(4I5,2ES20.8)') l1(chl(l,icl)), l2(chl(l,icl)), chl(l,icl), icl, &
              &                      w2e(l)%popd_Lll_1(icl), w2e(l)%popd_Lll_2(icl) 
      ENDDO
   ENDDO
   CLOSE(nout)
   !
   !
   
  
  DO l = 0, lmax
     PRINT*, 
     WRITE(*,'(i5,13X,3ES20.8)') l, w2e(l)%popd_L_1, w2e(l)%popd_L_2, SUM(pop_Lll_0(l,1:nll_max))
     DO icl = 1, w2e(l)%chl_n        
        WRITE(*,'(i5,a2,2i3,a5,3ES20.8)')&
             icl, " (", w2e(l)%chl_1(icl), w2e(l)%chl_2(icl), " ) = ", &
             & w2e(l)%popd_Lll_1(icl), w2e(l)%popd_Lll_2(icl), pop_Lll_0(l,icl)
     ENDDO
     PRINT*,"..............................."
  ENDDO
  PRINT*, pop_di_1, pop_di_2
  

  
  !
  !  pop_d_e12_channels_0 = dP(e1,e2)/de1de2            (for fixed l,le1,le2)
  !
  !
  !-  Marginal probability for e1  
  !  pop_d_e1_channels_0 = I_de2 [ dP(e1,e2)/de1de2]    (for fixed l,le1,le2)
  !
  !-  Marginal probability for e2  
  !  pop_d_e2_channels_0 = I_de1 [ dP(e1,e2)/de1de2]    (for fixed l,le1,le2)
  !
  !  pop_d_channels_0 = I_de_1_de2 [ dP(e1,e2)/de1de2]  (for fixed l,le1,le2)
  !
  !  S_{l,le1,le2} = pop_double_ionization (e1>0, e2>)
  !

  !
  !
  !
  !  PES 1D  P_Ll1l2(k1), P_Ll1l2(k2)
  !
  !          P_l1l2(k1) = Sum_{L}  P_Ll1l2(k1)
  !          P_l1l2(k2) = Sum_{L}  P_Ll1l2(k2) ! not stored at the moment. 
  !           
  !  icl == (l1,l2)
  !
  

  !
  ! electron '1'
  !
  
  save_PES_1D1_Ll1l2: DO l = 0, lmax   
     partial_channels_l:DO icl = 1, w2e(l)%chl_n
        
        WRITE(sl, '(I6)') l
        WRITE(sl1,'(I6)') icl
        
        data_file = TRIM(pop_file)//"1D1_LC"//TRIM(ADJUSTL(sl))//TRIM(ADJUSTL(sl1))//".dat"     
        OPEN(nout,file=data_file)
        
        le1 = w2e(l)%chl_1(icl)
        le2 = w2e(l)%chl_2(icl)

        n_cols = 2
        WRITE(nout,'(a1,4a5)') '&',"icl","l","l1", "l2" 
        WRITE(nout,'(a1,4I5)') '&',  icl, l, le1, le2
        WRITE(nout,'(1i5)')      n_cols
!        WRITE(nout,'(a1,1i5)') '&', n_cols
        DO k1 = 1, ne_1e_cut
           IF(e1e(le1,k1) > 0.0_dpk) THEN 
              WRITE(nout,'(3E15.6)') e1e(le1,k1), w2e(l)%popd_e1_Lll(icl,k1)
           ENDIF
        ENDDO
        CLOSE(nout)
        
     ENDDO partial_channels_l
  ENDDO save_PES_1D1_Ll1l2

  !
  !
  ! electron '2'
  ! 
  !
  
    save_PES_1D2_Ll1l2: DO l = 0, lmax   
       partial_channels_l_2:DO icl = 1, w2e(l)%chl_n
        
          WRITE(sl, '(I6)') l
          WRITE(sl1,'(I6)') icl
        
          data_file = TRIM(pop_file)//"1D2_LC"//TRIM(ADJUSTL(sl))//TRIM(ADJUSTL(sl1))//".dat"     
          OPEN(nout,file=data_file)
        
          le1 = w2e(l)%chl_1(icl)
          le2 = w2e(l)%chl_2(icl)
          n_cols = 2
          WRITE(nout,'(a1,4a5)') '&',"icl","l","l1", "l2" 
          WRITE(nout,'(a1,4I5)') '&', icl, l, le1, le2
          WRITE(nout,'   (1i5)')  n_cols
          !          WRITE(nout,'(a1,1i5)') '&', n_cols
          DO k2 = 1, ne_1e_cut
             IF(e1e(le2,k2) > 0.0_dpk) THEN 
                WRITE(nout,'(3E15.6)') e1e(le2,k2), w2e(l)%popd_e2_Lll(icl,k2)
             ENDIF
          ENDDO
          CLOSE(nout)          
       ENDDO partial_channels_l_2
    ENDDO save_PES_1D2_Ll1l2
  
  !
  !
  ! do the same for the (l1,l2) channels summed over L:
  !
  !


    !
    ! electron '1'
    !
    
  partial_channels_1:DO icll = 1, n_ll
        
     WRITE(sl1, '(I6)') icll
        
     data_file = TRIM(pop_file)//"1D1_C"//TRIM(ADJUSTL(sl1))//".dat"     
     OPEN(nout,file=data_file)

     le1 = l1(icll)
     le2 = l2(icll)
     
     n_cols = 2
     WRITE(nout,'(a1,3a5)') '&',"icll","l1","l2"
     WRITE(nout,'(a1,3I5)') '&', icll, le1, le2
     WRITE(nout,'(1i5)') n_cols
!     WRITE(nout,'(a1,1i5)') '&', n_cols
        
     DO k1 = 1, ne_1e_cut
        IF(e1e(le1,k1) > 0.0_dpk) THEN 
           WRITE(nout,'(3E15.6)') e1e(le1,k1), popl_e1(icll,k1)
        ENDIF
     ENDDO
     
     CLOSE(nout)
        
  ENDDO partial_channels_1

  !
  ! electron '2'
  !
 
  partial_channels_2:DO icll = 1, n_ll
        
     WRITE(sl1, '(I6)') icll
        
     data_file = TRIM(pop_file)//"1D2_C"//TRIM(ADJUSTL(sl1))//".dat"     
     OPEN(nout,file=data_file)

     le1 = l1(icll)
     le2 = l2(icll)
     
     n_cols = 2
     WRITE(nout,'(a1,3a5)') '&',"icll","l1","l2"
     WRITE(nout,'(a1,3I5)') '&', icll, le1, le2
     WRITE(nout,'(1i5)') n_cols
!     WRITE(nout,'(a1,1i5)') '&', n_cols

        
     DO k2 = 1, ne_1e_cut
        IF(e1e(le2,k2) > 0.0_dpk) THEN 
           WRITE(nout,'(3E15.6)') e1e(le2,k2), popl_e2(icll,k2)
        ENDIF
     ENDDO
     
     CLOSE(nout)
        
  ENDDO partial_channels_2
 
  

  !
  !  PES 2D  P_Ll1l2(k1,k2)
  !
  !          P_l1l2(k1,k2) = Sum_{L}  P_Ll1l2(k1,k2)

  !           
  !          icll == (l1,l2)
  !
  
  
  save_PES_2D_Ll1l2: DO l = 0, lmax   
     partial_channelsL:DO icl = 1, w2e(l)%chl_n
        
        WRITE(sl, '(I6)') l
        WRITE(sl1,'(I6)') icl

        
        data_file = TRIM(pop_file)//"2D_LC"//TRIM(ADJUSTL(sl))//TRIM(ADJUSTL(sl1))//".dat"     
        OPEN(nout,file=data_file)

                
        le1 = w2e(l)%chl_1(icl)
        le2 = w2e(l)%chl_2(icl)

        n_cols = 3
        WRITE(nout,'(a1,4a5)') '&',"icl","l","l1", "l2" 
        WRITE(nout,'(a1,4I5)') '&',  icl, l, le1, le2
        WRITE(nout,'(i5)'    )  n_cols
         
        save_e1e2_rho_L:DO k1 = 1, ne_1e_cut      !function values on the grid
           DO k2 = 1, ne_1e_cut
              IF( (e1e(le1,k1)>0).AND.(e1e(le2,k2)>0) )  THEN
                 WRITE(nout,'(3E16.8)') e1e(le1,k1),e1e(le2,k2),w2e(l)%popd_ee_Lll(icl,k1,k2)
              ENDIF
           ENDDO
        ENDDO save_e1e2_rho_L
        CLOSE(nout)
        
     ENDDO partial_channelsL
  ENDDO save_PES_2D_Ll1l2
  !  


  
  save_PES_2D_l1l2:DO icll = 1, n_ll
     !
     WRITE(sl1, '(I6)') icll
     !  
     data_file = TRIM(pop_file)//"2D_C"//TRIM(ADJUSTL(sl1))//".dat"     
     OPEN(nout,file=data_file)
     !
     le1 = l1(icll)
     le2 = l2(icll)
     !
     n_cols = 3
     WRITE(nout,'(a1,3a5)') '&', "icll","l1","l2"
     WRITE(nout,'(a1,3I5)') '&',  icll, le1, le2
     WRITE(nout,'(i5)'    ) n_cols
     !
     save_e1e2_rho:DO k1 = 1, ne_1e_cut      !function values on the grid
        DO k2 = 1, ne_1e_cut
           IF( (e1e(le1,k1)>0).AND.(e1e(le2,k2)>0) )  THEN
              WRITE(nout,'(3E16.8)') e1e(le1,k1), e1e(le2,k2), popl_e1e2(icll,k1,k2)
           ENDIF
           !           IF( (e1e(le1,k2)>0).AND.(e1e(le2,k1)>0) )  THEN
           !              WRITE(nout,'(3E16.8)') e1e(le1,k2), e1e(le1,k1), popl_e1e2(icll,k2,k1)
           !           ENDIF
        ENDDO
     ENDDO save_e1e2_rho     
     CLOSE(nout)
     !
  ENDDO save_PES_2D_l1l2


  
  !
  !
  !  
  !
  !> probabilities in terms of the CI 2e-functions  psi(t) = S_NL C_NL(t) Phi_NL(r_1,r_2)
  !
  !
  !
  !
  
  ALLOCATE(   popl_v(0:lmax) )  
  ALLOCATE(   popl_b_v(0:lmax) ) !bound
  ALLOCATE(   popl_s_v(0:lmax) ) !single ionization
  ALLOCATE(   popl_d_v(0:lmax) ) !double ionization

  
  popl_b_v = 0.0_dpk
  popl_s_v = 0.0_dpk
  popl_d_v = 0.0_dpk     
  populations_relative_to_ci_energies:DO l = 0, lmax
     !
     popl_b_v(l) = SUM(w2e(l)%pt(        1: nsi_l(l)-1) )   ! E_0 < E < E+
     popl_s_v(l) = SUM(w2e(l)%pt( nsi_l(l): ndi_l(l)-1) )   ! E+  < E < E++
     popl_d_v(l) = SUM(w2e(l)%pt( ndi_l(l): w2e(l)%net) )   ! E++ < E
     !   
  ENDDO populations_relative_to_ci_energies

  popl_v  = popl_b_v + popl_s_v + popl_d_v     ! all population in L-channgel bound + single + double
  norm_v = SUM(popl_v)                   ! sum over L-channels  (normalization to 1)

  ! printout results
    
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   bound population'
  PRINT*, 'x'
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                        E  < E+     pop_b_0 = ',   popl_b_0,  SUM(popl_b_0)
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                        E  < E+     pop_b_v = ',   popl_b_v,  SUM(popl_b_v)
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   single ionization'
  PRINT*, 'x'
  WRITE(*,'(A60,2X, 1P10E12.4)')  ' (e1>0 or e2>0) and  E+  < E < E++    pop_s_0 = ', popl_s_0,  SUM(popl_s_0)
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                     E+  < E < E++    pop_s_v = ', popl_s_v,  SUM(popl_s_v)
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   double ionization'
  PRINT*, 'x'
  WRITE(*,'(A60,2X, 1P10E12.4)')  ' (e1,e2 > 0) and           E++ < E    pop_d_0 = ', popl_d_0,  SUM(popl_d_0)
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                           E++ < E    pop_d_v = ', popl_d_v,  SUM(popl_d_v)
  PRINT*, 'x'  
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   total population'
  PRINT*, 'x'
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                                       norm_0 = ', norm_0, SUM(popl_0)  !`
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                                       norm_v = ', norm_v, SUM(popl_v)  !
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx '
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   double ionization'
  PRINT*, 'x'
  WRITE(*,'(A30,2X, 1P10E12.4)')  ' (e1,e2 > 0) and           E++ < E    pop_dd = '!, popl_d_0,   SUM(popl_d_0)
  PRINT*, 'x          L = 0,1,...,Lmax,       sum over L'
  PRINT*, 'x'
  PRINT*, 'x'  
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   (l1,l2) channels: 1=present, 0=absent'
  PRINT*, 'xxxx DI: uncorrelated basis projection'
  PRINT*, 'xxxx  l, l1, l2, pop_d(l,l1,l2), exists(1)'

  


END PROGRAM tdse_pop_sdi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!EOF
