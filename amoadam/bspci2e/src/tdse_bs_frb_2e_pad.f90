MODULE channel 
  USE PRECISION 
  TYPE channels
     !
     ! time-dependent output y(i)
     !
     DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: ct           !time-dependent coefficient for chann
     REAL(DPK), ALLOCATABLE, DIMENSION(:)      :: ct_en        !corresponding energy for ct
     INTEGER,   ALLOCATABLE, DIMENSION(:)      :: ch_index_en  !n
     INTEGER,   ALLOCATABLE, DIMENSION(:)      :: noch_en      !n
     INTEGER                                   :: ne           !nof channels for symmetry L
     !
     ! cfg channels information
     !
     INTEGER                                   :: nstates_L !1(free 2e state) or 0(bound 2e state)
     INTEGER,   ALLOCATABLE, DIMENSION(:)      :: n         !n
     INTEGER,   ALLOCATABLE, DIMENSION(:)      :: l         !l
     INTEGER,   ALLOCATABLE, DIMENSION(:)      :: lp        !l'
     INTEGER,   ALLOCATABLE, DIMENSION(:)      :: np_min    !n'_min
     INTEGER,   ALLOCATABLE, DIMENSION(:)      :: np_max    !n'_max
     INTEGER,   ALLOCATABLE, DIMENSION(:)      :: nfree     !1(free 2e state) or 0(bound 2e state)
     INTEGER                                   :: nch       !nof channels for symmetry L
     INTEGER                                   :: nchf      !nof free channels for symmetry L
     INTEGER                                   :: nchb      !nof bound channels for symmetry L
     INTEGER                                   :: nsi       !nof sionization channels for symmetry L
     INTEGER                                   :: ndi       !nof dionization channels for symmetry L
     !
     ! e1,e2 energies
     !
     REAL(DPK), ALLOCATABLE, DIMENSION(:)      :: en_t      !target energy e1 for channel (l,i)
     REAL(DPK), ALLOCATABLE, DIMENSION(:)      :: ek_2       !energy for e2 : e2 = ct_en - e1

     DOUBLE COMPLEX,ALLOCATABLE, DIMENSION(:,:,:):: dpsi_en   !dp_Single_Ionization at enery ie
     DOUBLE COMPLEX,ALLOCATABLE, DIMENSION(:,:,:):: dpdi_en   !dp_Single_Ionization at enery ie

  END TYPE channels
END MODULE channel
 
PROGRAM tdse_bs_mch_pad
   !
!   use par_tdse_bs_mch_pad
  USE PRECISION 
  USE units
  USE spherical_harmonics
  USE anglib
  USE io
  USE channel
  USE pad_utils

  !
  IMPLICIT NONE
  !
  !2e
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)         :: en
  INTEGER,   ALLOCATABLE, DIMENSION(:)           :: ne
  REAL(dpk), ALLOCATABLE, DIMENSION(:)           :: y
  REAL(dpk)                                      :: ecut
  INTEGER                                        :: nemax, nmax, ntot 
  INTEGER                                        :: lm, nl 
  INTEGER                                        :: nbsi, nbsf, li, lf, ntl, ch_index
  INTEGER                                        :: l_total, spin_total  
  REAL(dpk)                                      :: threshold_system  !(he++ energy axis)
  REAL(dpk)                                      :: z1,z2             ! effective SI and DI charges
  REAL(dpk)                                      :: E12, ken          ! k = sqrt(2E) !test
  REAL(dpk)                                      :: kx,ky,kz
  REAL(dpk)                                      :: ek2, ke2, kx_2,ky_2,kz_2
  REAL(dpk)                                      :: ek1, ke1, kx_1,ky_1,kz_1
  REAL(dpk)                                      :: clb,yr,yi,y2
  DOUBLE COMPLEX                                 :: al,all,tst
   !
  INTEGER,           ALLOCATABLE, DIMENSION(:,:) :: noch
  INTEGER,           ALLOCATABLE, DIMENSION(:)   :: nec
  INTEGER,           ALLOCATABLE, DIMENSION(:)   :: ns
  REAL(DPK),         ALLOCATABLE, DIMENSION(:,:) :: energy
  ! target data (he+)
  TYPE(channels), ALLOCATABLE, DIMENSION(:)      :: ch
  REAL(dpk)                                      :: en_target_i
  REAL(dpk)                                      :: delta_1, delta_2
  INTEGER                                        :: ich, k1
  REAL(dpk)                                      :: threshold_target  !(he+ energy axis)
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)         :: ent ! target energies
  INTEGER,   ALLOCATABLE, DIMENSION(:)           :: nbs ! nof bound target energies
  INTEGER                                        :: net ! nof target energies/pwave
  INTEGER                                        :: lmt, nlt, mt ! nof pwaves
  REAL(dpk)                                      :: eni ! initial target state energy (he+(1s))
  REAL(dpk)                                      :: en1,en2,en3 ! initial target state energy (he+(1s))
  !pad
  DOUBLE COMPLEX                                 :: zg
  !   DOUBLE COMPLEX,  ALLOCATABLE, DIMENSION(:,:)   :: y_lm
  DOUBLE COMPLEX                                 :: y_lm_sum, sph_arm_1, sph_arm_2 
  INTEGER                                        :: nlm
  REAL(dpk)                                      :: theta, phi, theta_p, phi_p
  REAL(dpk), ALLOCATABLE, DIMENSION(:)           :: th, ph
  REAL(dpk), ALLOCATABLE, DIMENSION(:)           :: th_p,ph_p
  INTEGER                                        :: ith_max, jph_max
  INTEGER                                        :: ith_p_max, jph_p_max, ip,jp 
  DOUBLE COMPLEX,ALLOCATABLE, DIMENSION(:,:)     :: dpsi      !dp_Single_Ionization (total signal)
  DOUBLE COMPLEX,ALLOCATABLE, DIMENSION(:,:)     :: dpdi      !dp_Double_Ionization (total signal)
  REAL(DPK),ALLOCATABLE, DIMENSION(:,:)          :: dpddi     !dp_Double_Ionization (total signal)
  REAL(dpk)                                      :: theta_step, phi_step
  INTEGER                                        :: m, m1,m2, l1, l2
  INTEGER                                        :: jdat, it
  LOGICAL                                        :: pad_di, pad_si 
  LOGICAL                                        :: check, debug
  !... tdse_
  REAL(DPK)                            tim, i0, omeg, tau
  INTEGER                              irespec, irestart
  !var
  INTEGER                i,j,l
  INTEGER                ie,il,in
  INTEGER                 L_INP, ie_inp
  CHARACTER(len=80)      cfg_file, dip_file
  CHARACTER(len=6)       s_l
  CHARACTER*25           coe_t, pad_file, ent_file
  INTEGER                ncoe, npad, nent, ncfg, ndip, ncheck
  INTEGER                nxy, nxz, nyz
  REAL(dpk)              zero
  !
  ! executable statements
  !
  
  ! i/o !
  
     CALL getarg(1, argv)
     READ(argv, *) zero
  !   CALL getarg(2, argv)
  !   READ(argv, *) pad_di
  
!  zero = 0.0D+00
  pad_si = .FALSE.
  pad_di = .TRUE.
  check  = .FALSE.
  debug  = .false.
  !
  npad  = 10
  ncoe  = 11
  nent  = 20
  ncfg  = 32
  ndip  = 33
  ncheck = 34

  !
  z1 = 2.0d+00

  !
  ent_file   ='dat/ent.dat'       !i/ target energies
  coe_t      ='tdat/coe.dat'      !i/ coefficients
  pad_file   ='pad/pad.dat'      !o/ output pad file 

  !i/read main input file

  OPEN(ninp,file="inp/tdse_bs_frb_2e.inp")
  READ(ninp,*)             nmax         ! Max no of States   in L = 0,1,2,. ltot-1
  READ(ninp,*)             nemax        ! Max no of Energies in L = 0,1,2,. ltot-1
  CLOSE(ninp)
  
  WRITE(*,*) "#                     nmax = ", nmax
  WRITE(*,*) "#                    nmaxL = ", nemax
  

  ! some settings


  threshold_target = 0.0D+0
  threshold_system = 2.0D+00
!  zi = (0,1)

!read coefficient data at final time.

   OPEN(ncoe,file=coe_t,form='unformatted',access='sequential')
   READ(ncoe) i0, omeg, tau, nl, ecut, tim
   WRITE(*,*) nl, ecut, tim
   lm = nl - 1
   WRITE(*,*) nmax,lm
   ALLOCATE( en(nmax,0:lm))
   ALLOCATE( ne(0:lm))
   READ(ncoe) (ne(l), l = 0, lm), ntot 
   WRITE(*,*) (ne(l), l=0, lm), ntot

   ALLOCATE( y(2*ntot))
   DO l = 0, lm
      READ(ncoe) ( en(ie,l), ie = 1, ne(l))        !relative to he+(1s)
   END DO
   READ(ncoe) (y(in), in = 1, 2*ntot)         ! save 
   CLOSE(ncoe)
!
   WRITE(*,*) "#                               i0 = ", i0
   WRITE(*,*) "#                             omeg = ", omeg
   WRITE(*,*) "#                              tau = ", tau
   WRITE(*,*) "#            nof symmetries     nl = ", nl
   WRITE(*,*) "#            max symmetry       lm = ", lm
   WRITE(*,*) "#                             Emax = ", ecut
   DO l = 0, nl-1
      WRITE(*,'( " #",25X,"   n(",i2,") = ",1X, i4 )') l, ne(l)
   ENDDO
   WRITE(*,*) "# tdse_pad_mch_pad::         ntot = ", ntot
   WRITE(*,*) "#       coefficient data file read. ", coe_t


   ! read target energies
   OPEN(nent,file=ent_file)

   READ(nent,*)  nlt, net

   lmt = nlt - 1

   ALLOCATE(ent(0:lmt,net))
   ALLOCATE(nbs(0:lmt))

   DO l = 0, lmt

      READ(nent,*) net

      nbs(l) = 0 
      DO i = 1, net
         READ(nent,*) ent(l,i)
         IF(ent(l,i).LT.0d+00) nbs(l) = nbs(l) + 1
      ENDDO
   ENDDO
   CLOSE(nent)

   ent = 0.5 * ent  ! move to a.u.

   eni = ent(0,1)   ! 1s
   en1 = eni
   en2 = ent(0,2)
   en3 = ent(0,3)

   WRITE(*,*) "#                nof p. waves  nlt = ", nlt
   WRITE(*,*) "#                max p. wave   ltm = ", nlt
   WRITE(*,*) "#                nof energies  net = ", net
   WRITE(*,*) "#                gst  energy   eni = ", eni
   WRITE(*,*) "#                ion threshold en1 = ", en1
   WRITE(*,*) "#                ion threshold en2 = ", en2
   WRITE(*,*) "#                ion threshold en3 = ", en3


   DO i = 0, nlt-1
      WRITE(*,'( " #",25X,"nbs(",i2,") = ",1X, i2 )') i, nbs(i)
   ENDDO
   WRITE(*,*) "#              target energies read. ", ent_file


   !
   !check coulombic phase shift

   make_checks:IF(check) THEN   
      OPEN(ncheck, FILE="pad/phase-shifts.dat")
   
      ken = 0.
      DO l = 1, lmt
         DO i = 1, 100
            ken = ken + DBLE(i)/10000.
            zg = CMPLX(DBLE(l + 1), - z1/ken)
            WRITE(ncheck,'(4E25.15)') ken, sigma_kl(l,z1,ken), AIMAG(gammln(zg)) 
         ENDDO
         WRITE(ncheck,*) "&"   !XMGRACE REQUIRED
      ENDDO
      WRITE(*,*) "#         Coulombic phase shift checked  s_kl = Im(ln(G(z)).  "
      
      CLOSE(ncheck)
   ENDIF make_checks

   
   ! read channel data: L,M_L, n,l,lp,np_min,np_max nfree  (inp/cfg-L.inp)
      
   ALLOCATE(ch(0:nl-1))

   DO L = 0, Lm 

      cfg_file   ="inp/cfg-"            !i/ channnel conffigurations 
      WRITE(s_l,'(i6)') L         
      cfg_file = trim(cfg_file)//trim(adjustl(s_l))//".inp"

      OPEN(ncfg,file=cfg_file)

      READ(ncfg,*) l_total, spin_total

      IF(l_total.NE.l) THEN
         WRITE(*,*) "#         error in input file : ", cfg_file
         WRITE(*,*) "#                           L = ", L
         WRITE(*,*) "#                     L_total = ", L_total
         STOP
      ENDIF

      READ(ncfg,*)  ch(L)%nch

      j = ch(L)%nch

      ALLOCATE( ch(L)%n(j)     )
      ALLOCATE( ch(L)%l(j)     )
      ALLOCATE( ch(L)%lp(j)    )
      ALLOCATE( ch(L)%np_min(j))
      ALLOCATE( ch(L)%np_max(j))
      ALLOCATE( ch(L)%nfree(j) )
      ALLOCATE( ch(L)%en_t(j)  )

      ch(L)%nsi  = 0
      ch(L)%ndi  = 0
      ch(L)%nchb = 0
      ch(L)%nchf = 0

      DO i = 1, j

         READ(ncfg,*) ch(L)%n(i),  ch(L)%l(i),                        &
&                     ch(L)%lp(i), ch(L)%np_min(i), ch(L)%np_max(i),  &
&                     ch(L)%nfree(i)


!         WRITE(*,*) l,i,  ent(l, ch(l)%n(i) ),  ch(l)%nfree(i)
         IF(ch(L)%nfree(i) == 1 ) THEN 

            ch(L)%en_t(i) = ent(L, ch(l)%n(i) ) 

            IF(ch(L)%en_t(i) < 0.0D+00) THEN 
               ch(L)%nsi = ch(L)%nsi + 1          ! nof single ionization channels
            ELSE
               ch(L)%ndi = ch(L)%ndi + 1          ! nof double ionization channels
            ENDIF
            ch(L)%nchf = ch(L)%nchf + 1           ! nof ionization channels
         ELSE
            ch(L)%nchb = ch(L)%nchb + 1           ! nof bound channels
         ENDIF
      ENDDO

      IF(L==0)  THEN 
         WRITE(*,'( " # channel_file,      l     nch,    nchf,   ndi,    nsi,   nchb")')
         WRITE(*,*)"#---------------------------------------------------------------"
      ENDIF
      WRITE(*,'(" # ",a15,2X,I3,5X,I3,5X,I3,5X,I3,5X,I3,5X,I3 )') cfg_file,l, &
                            ch(L)%nch, &
              &             ch(L)%nchf,& 
              &             ch(L)%ndi ,&
              &             ch(L)%nsi ,& 
              &             ch(L)%nchb 


      !check! nof channels = nof continuum + nof bound channels
      IF( ch(L)%nch.NE.(ch(L)%nchf + ch(L)%nchb) )  THEN 
         WRITE(*,*) "#         error in input file : ",  cfg_file
         WRITE(*,*) "#                        nch  = ",  ch(L)%nch
         WRITE(*,*) "#                 nchf + nchb = ",  ch(L)%nchf + ch(L)%nchb
         STOP      
      ENDIF
      ! check! nof continuum channels = SI + DI ionization channels
      IF( ch(l)%nchf.NE.( ch(l)%nsi  + ch(l)%ndi ) ) THEN
         WRITE(*,*) "#         error in input file : ", cfg_file
         WRITE(*,*) "#                       nchf = ",  ch(L)%nchf
         WRITE(*,*) "#                  nsi + ndi = ",  ch(L)%nsi  + ch(L)%ndi
         STOP
      ENDIF
   ENDDO
   WRITE(*,*)"#---------------------------------------------------------------"
   WRITE(*,*)"# target configuration files read."
   
   CLOSE(ncfg)


   ! check dipole matrix elements !
   ! decompose system state vector (y) to the various channels
   ! based on en(ie,L) and ch(L) inputs
   !

!!%   ALLOCATE(energy(nemax,0:lm) ) 
!!%   ALLOCATE(  noch(nemax,0:lm) )
!!%   ALLOCATE(   nec(0:lm) )
!!%   ALLOCATE(   ns(0:lm) )
!!%!   ALLOCATE(  nbs(0:lm) )
!!%
!!%
!!%   ntl = 0
!!%   DO l = 0, lm-1
!!%      li = l
!!%      lf = l+1
!!%      dip_file   ="dat/os"            !i/ channnel configurations 
!!%      WRITE(s_l,'(i6)') li         
!!%      dip_file = TRIM(dip_file)//TRIM(ADJUSTL(s_l))
!!%      WRITE(s_l,'(i6)') lf
!!%      dip_file = TRIM(dip_file)//TRIM(ADJUSTL(s_l))//".dat"
!!%      !         
!!%      !data
!!%      !01          de1 = 0.008   (E<0) in a.u.
!!%      !            de2 = 0.0175  (E>0)  
!!%      !12          de1 = 0.0175  (E<0)
!!%      !            de2 = 0.0175  (E>0)
!!%      !   
!!%      OPEN(ndip, file = dip_file)
!!%      READ(ndip,*) li,lf, ns(L), ns(L+1), nbsi, nbsf, nec(L), nec(L+1)
!!%      
!!%      IF(ns(l).GT.nmax) WRITE(*,*)    '# tdse_bs_mch::     ni > nmax,   ni = ', ns(L)
!!%      IF(ns(l+1).GT.nmax) WRITE(*,*)  '# tdse_bs_mch:: n(i+1) > nmax,   nf = ', ns(L+1)
!!%      
!!%      READ(ndip,*) ( energy(ie, L   ),  ie =        1, nec(L)   )        ! in Ryd
!!%      READ(ndip,*) ( energy(ie, L+1 ),  ie =        1, nec(L+1) )        ! in Ryd
!!%      READ(ndip,*) (   noch(ie, L   ),  ie = nbsi + 1, nec(L)   ) 
!!%      READ(ndip,*) (   noch(ie, L+1 ),  ie = nbsf + 1, nec(L+1) )
!!%      CLOSE(ndip)
!!%      
!!%      nbs(L) = nbsi
!!%      IF(l==lm-1) nbs(l+1) = nbsf
!!%      noch(1:nbsi, l)   = 1
!!%      noch(1:nbsf, l+1) = 1
!!%      !         write(*,*) lb, lz, n(i), n(i+1), nbs(i), nbs(i+1), ne(i), ne(i+1)
!!%      WRITE(*,*) "#-------------------------------------------------------"
!!%      WRITE(*,*) "#                               li = ", li
!!%      WRITE(*,*) "#                              nei = ", ns(L)
!!%      WRITE(*,*) "#                             nbsi = ", nbs(L)
!!%      WRITE(*,*) "#                              nec = ", nec(L)
!!%      WRITE(*,*) "#------"
!!%      WRITE(*,*) "#                               lf = ", lf
!!%      WRITE(*,*) "#                              nef = ", ns(L+1)
!!%      WRITE(*,*) "#                             nbsf = ", nbsf
!!%      IF(l==lm-1) THEN
!!%      WRITE(*,*) "#                             nbsf = ", nbs(L+1)
!!%   ENDIF
!!%      WRITE(*,*) "#                              nec = ", nec(L+1)
!!%   
!!%   END DO
!!%   WRITE(*,*) "#-------------------------------------------------------"
!!%      



   OPEN(nout,file="dat/yt.dat") 
   OPEN(ndat,file="pad/yt-channel.dat") 
   READ(nout,'(5e14.6)')         i0,  omeg, tau
   READ(nout,'(2e14.6,1X,2I5)') tim, ecut, nl, ntot

   lm = nl - 1
   ntl = 0 
   DO L = 0, lm

      READ(nout,'(2I8)')  ch(L)%ne, l_inp
      IF(L.NE.l_inp) THEN
         WRITE(*,*)"# input error in y(t).dat file:"
         WRITE(*,*)"#                           L = ", L 
         WRITE(*,*)"#                       L_INP = ", L_INP 
         STOP
      ENDIF

      ALLOCATE( ch(L)%ct         ( ch(L)%ne ) )
      ALLOCATE( ch(L)%ct_en      ( ch(L)%ne ) )
      ALLOCATE( ch(L)%ch_index_en( ch(L)%ne ) )
      ALLOCATE( ch(L)%noch_en    ( ch(L)%ne ) )
      ALLOCATE( ch(L)%ek_2       ( ch(L)%ne)  )

      DO ie = 1, ch(L)%ne

         READ(nout,'(4I8,2X,4e14.6)')      ie_inp                  &
                                       &,  l_inp                   &
                                       &, ch(L)%noch_en(ie)        &
                                       &, ch(L)%ch_index_en(ie)    &
                                       &, ch(L)%ct_en(ie)          &
                                       &, yr                       &
                                       &, yi &
                                       &, y2

         ch(L)%ct(ie) = CMPLX(yr,yi)


         WRITE(ndat,'(7I8,2X,4e14.6)')     ie                               & ! i-state
                                       &,  l                                & ! L symmetry
                                       &, ch(L)%noch_en(ie)                 & ! noch(i,L) 
                                       &, ch(L)%ch_index_en(ie)             & ! jth-channel at (i,L) state.
                                       &, ch(L)%n( ch(L)%ch_index_en(ie) )  & !   n(ie,j,l)
                                       &, ch(L)%l( ch(L)%ch_index_en(ie) )  & !   l(ie,j,l)
                                       &, ch(L)%lp( ch(L)%ch_index_en(ie) ) & !  lp(ie,j,l)
                                       &, ch(L)%ct_en(ie)                   & !  E(ie,L)
                                       &, REAL(ch(L)%ct(ie))                & !  re(coef(ie,L))
                                       &, AIMAG(ch(L)%ct(ie))               & !  re(coef(ie,L))
                                       &, ABS(ch(L)%ct(ie))**2                !  re(coef(ie,L))

         IF(L.NE.l_inp) THEN
            WRITE(*,*)"# input error in y(t).dat file: "
            WRITE(*,*)"#                           L = ",L 
            WRITE(*,*)"#                       L_INP = ",L_INP 
            STOP
         ENDIF

      ENDDO
   ENDDO
   
   CLOSE(nout)
   CLOSE(ndat)

   WRITE(*,*) "# wavefunction decomposition done."

                 
   !... theta, phi discretization

      ith_max = 200             ! theta discretization
      jph_max = 3               ! phi   discretization

      ALLOCATE(th(ith_max+1))
      ALLOCATE(ph(jph_max+1))

      !theta
      DO i = 1, ith_max+1
         th(i) = DBLE(i-1) * (2.0*m_pi)/ith_max
      ENDDO
      !phi
      DO j = 1, jph_max+1 
         ph(j)  = DBLE(j-1) * (2.0*m_pi)/jph_max
      ENDDO

   WRITE(*,*) "# theta and phi discretization grid done."

      !calculate pad and momentum distribution for single-ionized channels




   calculate_pad_si: IF(pad_si) THEN


      ALLOCATE(dpsi(ith_max+1,jph_max+1))

      clb     = 0.0d+00 
      al      = 0.0d+00
      delta_2 = 0.0D+00


      dpsi = 0.0D+00

      total_angular_momentum: DO L = 0, nl - 1        ! 2e-channels (L)

         ALLOCATE( ch(L)%dpsi_en( ith_max+1, jph_max+1, ch(L)%ne) )

            ch(L)%dpsi_en = 0.0d+00
            ch(L)%ek_2    = 0.0d+00


            total_energy: DO ie = 1, ch(L)%ne


               E12 = ch(L)%ct_en(ie)           ! total energy measured from he+(1s)
               
               IF(ch(L)%noch_en(ie) > 1) CYCLE     ! take only (1s) ionic core states
               IF(E12<0.0d+00)           CYCLE     ! exclude bound states

               k1 =   ch(L)%n(  ch(L)%ch_index_en(ie) )   
               l1 =   ch(L)%l(  ch(L)%ch_index_en(ie) )
               l2 =   ch(L)%lp( ch(L)%ch_index_en(ie) )

               ek1 =  ent(l1, k1) 
               ek2 =  E12 - (ek1-en1)
               ke2 = SQRT(2.0D+00*ek2) 
               ch(L)%ek_2(ie) =  ek2
               

               delta_2 = sigma_kl(l2, z2, ke2)
               all = zi**l2 * EXP(-zi*delta_2)


!                  ek1 =  ch(L)%en_t(ich)!ch(L)%en_t(ich)

               
               IF(debug) THEN                  
                  IF(l==1) WRITE(*,'(2I4,1X,5f12.4)') ie, l, ek1, ek2*27.211, ke2, delta_2, E12
               ENDIF

                  !x
                  !x IF( (ch(l)%nfree(ich) == 0)) CYCLE   ! ionization channels only
                  !x IF  ( ek1  > 0.0D+00       ) CYCLE   ! only single ionization channels
                  !x IF  ( E12  < (ek1-en1)     ) CYCLE   ! exclude states E < e_1  (closed channels)
                  !x
                  !x go ahead
                  !x

               cleb_m2: DO m2 = -l2, l2
                  m1 = -m2

                  clb = cleb( 2*l1, 2*m1, 2*l2, 2*m2, 2*l, 0 )
                  al = clb * all

                  theta_loop: DO i = 1, ith_max+1
                     phi_loop : DO j = 1, jph_max+1

                        sph_arm_2 = sphharm(l2, m2, th(i), ph(j))
                           

                        ch(L)%dpsi_en(i,j,ie) =  ch(L)%dpsi_en(i,j,ie) +  &
                                & al *  ch(L)%ct(ie) * CONJG( sph_arm_2 )

                        dpsi(i,j) = dpsi(i,j) +  ch(L)%dpsi_en(i,j,ie)
                           
!                           WRITE(ndat, '(4E20.8)') kx_2,ky_2,kz_2, ABS(ch(L)%dpsi(i,j,ie))**2

                           
                     ENDDO phi_loop
                  ENDDO theta_loop

               ENDDO cleb_m2

               
!!               WRITE(*,'(5I4,2X,3f12.6,2X,I4,1X,E12.6)') ie, L, k1, l1, l2, E12, ek1, ch(L)%ek2(ie), &
!!                    & ch_index,ABS(ch(l)%dpsi(1,1,ie))**2
    
            ENDDO Total_energy
         ENDDO total_angular_momentum

!!pad!!

         !theta
         OPEN(ndat, file='pad/theta.dat')
         store_theta : DO i = 1, ith_max+1
            WRITE(ndat, "(1E20.5)") th(i)
         ENDDO store_theta
         CLOSE(ndat)
         !phi
         OPEN(ndat, file='pad/phi.dat')
         store_phi: DO j = 1, jph_max+1
            WRITE(ndat,"(1E20.5)") ph(j)
         ENDDO store_phi
         CLOSE(ndat)
         !pad-s.dat

         OPEN(ndat, file='pad/pad-s.dat')
         store_pad_single : DO j = 1, jph_max+1
            WRITE(ndat, '(101E20.5)')( ABS(dpsi(i,j))**2, i = 1,ith_max+1)
         ENDDO store_pad_single
         CLOSE(ndat)



!!% momentum distributions



!!%         nxy = 1
!!%         nxz = 2
!!%         nyz = 3
!!%         OPEN(ndat, file='pad/pcolor.dat')
!!%         OPEN(nxy,  file='pad/kxy.dat')
!!%         OPEN(nxz,  file='pad/kxz.dat')
!!%         OPEN(nyz,  file='pad/kyz.dat')
!!%
!!%
!!%         L_wave: DO L = 0, nl - 1        ! 2e-channels (L)
!!%            WRITE(*,*)  L, ch(L)%ne
!!%            CH_wave:  DO ie = 1,  ch(L)%ne
!!%               
!!%               ke2 = SQRT( 2.0d+00* ch(L)%ek_2(ie) )
!!%               
!!%               pcolor_th: DO i = 1, ith_max
!!%                  pcolor_phi : DO j = 1, jph_max
!!%
!!%                     IF( ABS(ch(L)%dpsi_en(i,j,ie))**2 >= 1.e-14 ) THEN
!!%
!!%                     kx = ke2 * SIN(th(i)) * COS(ph(j))
!!%                     ky = ke2 * SIN(th(i)) * SIN(ph(j))
!!%                     kz = ke2 * COS(th(j))
!!%
!!%                     WRITE(ndat, '(4E15.6)')kx,ky,kz,ABS(ch(L)%dpsi_en(i,j,ie))**2
!!%
!!%
!!%!                     IF(debug) THEN
!!%!                        WRITE(ndat, '(2I3,1X,2f10.4,2X,3E15.6,1X,2I3)') ie, L, ch(L)%ek_2(ie)*enau,     &
!!%!                             &                                                 ch(L)%ct_en(ie)*enau,    &
!!%!                             &                                           ABS(ch(L)%dpsi_en(i,j,ie))**2, &
!!%!                             &                                                         th(i),ph(j),i,j
!!%!                     ENDIF
!!%                  ENDIF
!!%               ENDDO pcolor_phi
!!%            ENDDO pcolor_th
!!%
!!%
!!%               !xy  (theta = pi/2)
!!%               xy_ph : DO j = 1, jph_max 
!!%                  IF( ABS(ch(L)%dpsi_en(20,j,ie))**2 >= 1.e-14 ) THEN
!!%                     kx_2 = ke2 * COS(ph(j))*SIN(th(20))
!!%                     ky_2 = ke2 * SIN(ph(j))*SIN(th(20)) 
!!%                     WRITE(nxy, '(4E20.8)')  kx_2, ky_2, ABS(ch(L)%dpsi_en(20,j,ie))**2
!!%                  ENDIF
!!%               ENDDO xy_ph
!!%               
               !xz (phi = 0)
!!%               xz_th: DO i = 1, ith_max
!!%                  IF( ABS(ch(L)%dpsi_en(i,1,ie))**2 >= 1.e-14 ) THEN
!!%                     kx_2 = ke2 * SIN(th(i))
!!%                     kz_2 = ke2 * COS(th(i)) 
!!%                     WRITE(nxz, '(4E20.8)') kx_2,kz_2, ABS(ch(L)%dpsi_en(i,1,ie))**2
!!%                  ENDIF
!!%               ENDDO xz_th

               !yz (phi = pi/2)
!!%               yz_th: DO i = 1, ith_max
!!%                  IF( ABS(ch(L)%dpsi_en(i,jph_max/4,ie))**2 >= 1.e-14 ) THEN
!!%                     ky_2 = ke2 * SIN(th(i))
!!%                     kz_2 = ke2 * COS(th(i))
!!%                     WRITE(nyz, '(4E20.8)') ky_2,kz_2, ABS(ch(L)%dpsi_en(i,jph_max/4,ie))**2
!!%                  ENDIF
!!%               ENDDO yz_th               !
!!%            ENDDO CH_WAVE
!!%         ENDDO L_WAVE
!!%      
!!%         CLOSE(ndat)
!!%         CLOSE(nxy)
!!%         CLOSE(nxz)
!!%         CLOSE(nyz)
!!%
        

      ENDIF calculate_pad_si



      ! double  ionization photoelectron angular distributions (pad_di)


      calculate_pad_di: IF(pad_di) THEN
         
         theta_p = zero*M_PI
         phi_p   = 0

      !theta_p
!      DO i = 1, ith_max+1
!         th(i) = DBLE(i-1) * m_pi/ith_max
!      ENDDO
      !phi_p  
!      DO j = 1, jph_max+1 
!         ph(j)  = DBLE(j-1) * (2.0*m_pi)/jph_max
!      ENDDO

         ALLOCATE(dpdi(ith_max+1,jph_max+1))
         ALLOCATE(dpddi(ith_max+1,jph_max+1))
         dpdi = 0.0D+00
         dpddi = 0.0D+00

!    WRITE(*,*) "ie,   l,  k1, l1, l2,      E12,        ek1,        ek2,         ke1,       ke2"



         total_angular_momentum_di: DO L = 0, nl - 1

            ALLOCATE( ch(L)%dpdi_en( ith_max+1, jph_max+1, ch(L)%ne) ) ;  
            
            ch(L)%dpdi_en   = 0.0d+00
            ch(L)%ek_2      = 0.0d+00

            total_energy_di: DO ie = 1, ch(L)%ne


               E12 = ch(L)%ct_en(ie)                  ! total energy measured from he+(1s)
               IF(E12<0.0d+00)           CYCLE        ! exclude bound + single ionization states

               k1 =   ch(L)%n(  ch(L)%ch_index_en(ie) )
               l1 =   ch(L)%l(  ch(L)%ch_index_en(ie) )
               l2 =   ch(L)%lp( ch(L)%ch_index_en(ie) )

               ek1 =  ent(l1, k1)

               di:IF( ek1 > 0.0D+00 ) THEN

                  !               WRITE(*,*) 'ie,ek1, E12 = ', ie, ek1, E12
               

                  z1 = 2.0D+00
                  Z2 = 2.0D+00
                  ke1 = SQRT(2.0d+00 * ek1)                  ! electron '1'
                  delta_1 = sigma_kl(l1, z1, ke1)
               !
                  ek2     =  E12 - (ek1-en1)                 ! electron '2'
                  ke2     = SQRT(2.0D+00*ek2) 
                  delta_2 = sigma_kl(l2, z2, ke2)
                  
                  ch(L)%ek_2(ie) =  ek2
                  

                  all = zi**l2 * EXP(-zi*delta_2) * zi**l1 * EXP( -zi* delta_1 )


               IF(L.EQ.0.AND.ie.GE.882.AND.ie.LE.920) THEN               
                  WRITE(*,'(5I4,1X,7f12.4,1E15.7)')   ie,   l,  k1,  l1,  l2,      &
                                                   &      E12, ek1, ek2, ke1, ke2, & 
                                                   & delta_1, delta_2,             &
                                                   & ABS(ch(L)%ct(ie))**2
               ENDIF
               IF(L.EQ.2.AND.ie.GE.1090.AND.ie.LE.1131) THEN

                  WRITE(*,'(5I4,1X,7f12.4,1E15.7)') ie,   l,  k1,  l1,  l2,        &
                                              &         E12, ek1, ek2, ke1, ke2,   &
                                              &         delta_1, delta_2,          & 
                                              &         ABS(ch(L)%ct(ie))**2
               ENDIF

              cleb_m:DO m2 = -l2,l2
                 m1 = -m2
                     
                 clb = cleb( 2*l1, 2*m1, 2*l2, 2*m2, 2*L, 0 )

!                 WRITE(*,'(1G10.5,5I4)') clb, l1,l2,L,m1,m2
                 al = clb * all 

                 IF(clb==0.0d+00) CYCLE

!              WRITE(*,'(2I4,1X,8f12.4)') ie, l, ek1, ek2, ke1, ke2, delta_1, delta_2, E12, clb

                  theta_lp: DO i = 1, ith_max+1
                     phi_lp : DO j = 1, jph_max+1

!                        IF(th(i).GT.m_pi)    thi(i) = 2*m_pi - thi(i)

                        sph_arm_2 = sphharm(l2, m2,   th(i), ph(j) )
                        sph_arm_1 = sphharm(l1, m1, theta_p, phi_p )

!                        IF(th(i).GT.m_pi)    sph_arm_2 = (-1)**(l2)*sph_arm_2

                                                      
!                     dpdi(i,j) = dpdi(i,j) +  al * CONJG( sph_arm_1 )  &
!                                                 * CONJG( sph_arm_2 )  &
!                                                 * ch(L)%ct(ie)

                        
                        ch(L)%dpdi_en(i,j,ie) = ch(L)%dpdi_en(i,j,ie) +  al * CONJG( sph_arm_1 ) &
                                               * CONJG( sph_arm_2 ) &
                                               * ch(L)%ct(ie)


!                           dpdi(i,j) = dpdi(i,j) + ch(L)%dpdi_en(i,j,ie)

!L=0
!(1) 915,914,913
!(2) 903,904,905
!(3) 898,899,900
                       IF(L.EQ.0.and.ie.eq.915) THEN
                           dpdi(i,j) = dpdi(i,j) + ch(L)%dpdi_en(i,j,ie) !ss
                        ENDIF
                       IF(L.EQ.0.and.ie.eq.914) THEN
                           dpdi(i,j) = dpdi(i,j) + ch(L)%dpdi_en(i,j,ie) !pp
                        ENDIF
                        IF(L.EQ.0.AND.ie.eq.913) THEN
                           dpdi(i,j) = dpdi(i,j) + ch(L)%dpdi_en(i,j,ie) !dd
                        ENDIF 

!L=2
!(1)1131,1129
!(2)1123,1121

                        IF(L.EQ.2.AND.ie.eq.1129) THEN
                         dpdi(i,j) = dpdi(i,j) + ch(L)%dpdi_en(i,j,ie)   !sd
                        ENDIF
                        IF(L.EQ.2.AND.ie.eq.1130) THEN
                           dpdi(i,j) = dpdi(i,j) + ch(L)%dpdi_en(i,j,ie) !pf
                        ENDIF
                        IF(L.EQ.2.AND.ie.eq.1131) THEN
                           dpdi(i,j) = dpdi(i,j) - ch(L)%dpdi_en(i,j,ie) !pp
                        ENDIF 


!!%                        IF(L.EQ.3.AND.ie.eq.1115) THEN 
!!%                           dpdi(i,j) = dpdi(i,j) + ch(L)%dpdi_en(i,j,ie)
!!%                        ENDIF
!!%                        IF(L.EQ.3.AND.ie.eq.1116) THEN
!!%                           dpdi(i,j) = dpdi(i,j) + ch(L)%dpdi_en(i,j,ie)
!!%!                           dpddi(i,j) = dpddi(i,j) + ABS(ch(L)%dpdi_en(i,j,ie))**2
!!%                        ENDIF


                  ENDDO phi_lp
               ENDDO theta_lp
            ENDDO cleb_m
!
         ENDIF di
!
      END DO total_energy_di
   ENDDO total_angular_momentum_di


         !pad-d.dat
         OPEN(ndat, file='pad/pad-2d-d.dat')
         store_pad_double : DO i = 1, ith_max+1
            WRITE(ndat, '(2E20.5)')180*th(i)/m_pi ,ABS(dpdi(i,1))**2
         ENDDO store_pad_double

!         store_pad_double_negative_axis : DO i = ith_max+1, 1,-1
!            WRITE(ndat, '(2E20.5)') 360 - 180*th(i)/ M_PI , ABS(dpdi(i,1))**2
!         ENDDO store_pad_double_negative_axis

         CLOSE(ndat)

         !pad-d-polar.dat
         OPEN(ndat, file='pad/pad-2d-polar.dat')
         store_pad_double_polar : DO i = 1, ith_max+1
            WRITE(ndat, '(2E20.5)') SIN(th(i))*ABS(dpdi(i,1))**2, COS(th(i))*ABS(dpdi(i,1))**2
         ENDDO store_pad_double_polar

!         store_pad_double_polar_negative_axis : DO i = 1, ith_max
!            WRITE(ndat, '(2E20.5)') -SIN( th(ith_max+1-i) ) * ABS(dpdi(ith_max+1-i,1))**2, &
!                 COS( th(ith_max+1-i) ) * ABS(dpdi(ith_max+1-i,1))**2                                  
!         ENDDO store_pad_double_polar_negative_axis


         CLOSE(ndat)


         OPEN(ndat, file='pad/pad-3d-d.dat')
         store_pad_ddouble : DO j = 1, jph_max+1
            WRITE(ndat, '(101E20.5)')(ABS(dpdi(i,j))**2, i = 1, ith_max+1)
         ENDDO store_pad_ddouble
         CLOSE(ndat)
         

      END IF calculate_pad_di

      
! store pad for gnuplot


!!%!theta
!!%   OPEN(ndat, file='tmp/theta-info.dat')
!!%   DO i = 1, ith_max
!!%     WRITE(ndat, "(1E20.5,2X,1f                    12.5, 2X,I3)") th(i), th(i)*180/m_pi, i
!!%  ENDDO
!!%   CLOSE(ndat)
!!%!phi
!!%   OPEN(ndat, file='tmp/phi-info.dat')
!!%   DO j = 1, jph_max
!!%      WRITE(ndat,"(1E20.5,2X,1f12        .5,2X, I3)") ph(j),  ph(j)*360/(2.*m_pi), j
!!%   ENDDO
!!%   CLOSE(ndat)
!!%

         
!pad-double
!!%   OPEN(ndat, file='tmp/pad-d.dat')
!!%   store_pad_double : DO j = 1, jph_max+1
!!%      WRITE(ndat, '(201E20.5)')( ABS(dp_di(1,1,i,j))**2, i = 1,ith_max+1)
!!%   ENDDO store_pad_double
!!%   CLOSE(ndat)

!   DEALLOCATE(y_lm)
!   DEALLOCATE(en,y,ne)
!   DEALLOCATE(ent,nbs)
!   DEALLOCATE(ph,th)
!   DEALLOCATE(dp_di)   
   
   write(*,*) " # pad:: data stored in dat/theta.dat dat/phi.dat dat/pad.dat."
   write(*,*) " # pad::                                          nof theta = ", ith_max 
   write(*,*) " # pad::                                          nof   phi = ", jph_max 


! format statements

100 FORMAT ("# spse::   (TDSE_BS_MCH_PAD)")
 94 format ("# spse::              :                        system = ", A6)
204 format ("# spse::              :     internuclear distance R0  = ", f9.6)
104 format ("# spse::              : gs partial waves lmax, lmax_i = ", I6, I6)
114 format ("# spse::              :     box boundary         rmax = ", 1E15.6)
124 format ("# spse::              :     nof grid points        nr = ", I6)
134 format ("# spse::              :     nof time steps     nsteps = ", I7)
154 format ("# spse::              :     nuclear charge       znuc = ", f9.6)
164 format ("# spse::              :           grid step        dr = ", f9.6)
174 format ("# spse::              :           time step        dt = ", 1E15.6)
 91 format ("# spse::              :                     Real time = ", /)
101 format ("# spse::              :     Field     amplitude    E0 = ", 1E15.6, " a.u.   ")
201 format ("# spse::              :     Intensity amplitude    I0 = ", 1E15.6, " a.u. = ", 1E15.6, " W/cm^2")
111 format ("# spse::              :     Field frequency        wL = ", 1E15.6, " a.u. = ", 1E15.6, " eV." )
211 format ("# spse::              :     Pulse duration         Tp = ", 1E15.6, " a.u. = ", 1E15.6, " fs." )
112 format ("# spse::              :     Field shape    pulse type = ", A)
311 format ("# spse::              :     Optical cycle          TL = ", 1E15.6, " a.u. = ", 1E15.6, " fs." )
304 format ("# spse::              :     nof cycles        ncycles = ", I6)
404 format ("# spse::              :     nof      ncycles_no_field = ", I6)
411 format ("# spse::              :     Total duration       Tfin = ", 1E15.6, " a.u. = ", 1E15.6, " fs." )
121 format ("# spse::              :     Field/Molecule angle beta = ","(",f4.1," deg)")
!131 format ("# spse::              :                  Output files : ", A)
701 format ("# spse::              :              wf at file : ", A)
702 format ("# spse::              :            pad  at file : ", A)
!
!...........
!
END PROGRAM tdse_bs_mch_pad

 


