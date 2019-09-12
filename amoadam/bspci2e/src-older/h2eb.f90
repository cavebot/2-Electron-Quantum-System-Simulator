!############################################################
!
!
PROGRAM h2eb
  !
  USE precision, ONLY:dpk
  USE w2e, ONLY: read_cfg_file, diagonalize
  !
  IMPLICIT NONE
  !
  INTEGER                                :: ncs, nhx             
  REAL(dpk), ALLOCATABLE,DIMENSION(:,:)  :: h, u
  INTEGER,   DIMENSION(:),   POINTER     :: nhf, lhf             !ncs
  INTEGER,   DIMENSION(:),   POINTER     :: ll, nllmin, nllmax   !ncs
  INTEGER,   DIMENSION(:),   POINTER     :: idcs                 !ncs
  !
  INTEGER, ALLOCATABLE,   DIMENSION(:)   :: is, nd, n            !ncs
  INTEGER, ALLOCATABLE,   DIMENSION(:)   :: ih                   !nhx
  !
  INTEGER                                :: lo, inpl, ls
  INTEGER                                :: ncmx
  INTEGER                                :: ncfg, nh2eb, nv2eb
  !var!
  INTEGER                                :: k, kk, ir, ic
  INTEGER                                :: ihd, id 
  !i/o!
  CHARACTER(len=100)                     :: argv 
  !exe!

  !i/o files

  ncfg  = 1
  nh2eb = 2
  nv2eb = 3


  ! commamd line argument

  
  CALL GETARG(1, ARGV)
  READ(ARGV,*) inpl

  WRITE(*,*) '& h2eb::        partial wave L = ', inpl 
  
  ! read cfg-file for (l)
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



  !xxxxxxxxxxxxxx   read  r12 data


!  CALL hfile(nr12,"dat","v2eb","bin",lo)

  ALLOCATE( u(nhx,nhx) )

  CALL r12file(nv2eb, lo)         ! read r12 data

  
  read_v12_upper: DO  ir = 1, SIZE(u,dim=1)
     
     READ(nv2eb) ( u(ir,ic), ic = ir, SIZE(u,dim=2))
!     WRITE(*,*) u(1,4), u(1,5),u(1,6)
     
     fill_lower_part:DO  ic = ir, SIZE(u,dim=2)                   !  fill matrix since symmetric 

        u(ic,ir) = u(ir,ic)
        
     ENDDO fill_lower_part
  ENDDO read_v12_upper
  
  CLOSE(nv2eb)
  
  !xxxxxxxxxxxxxxxxxx     hamiltonian matrix


  ALLOCATE(    h(nhx,nhx) )
  
  CALL hfile(nh2eb,"dat","h2eb","bin",lo)       !  CALL d2efile(9, lo)


  
  WRITE(nh2eb) lo,ls
  WRITE(nh2eb) ncs
  WRITE(nh2eb) (    nhf(k), k = 1, ncs)
  WRITE(nh2eb) (    lhf(k), k = 1, ncs)
  WRITE(nh2eb) (     ll(k), k = 1, ncs)
  WRITE(nh2eb) ( nllmin(k), k = 1, ncs)
  WRITE(nh2eb) ( nllmax(k), k = 1, ncs)
  WRITE(nh2eb) (     nd(k), k = 1, ncs)

  !


  
  ALLOCATE(ih(nhx))

  ih = 0

  define_ihd:DO  k = 1, ncs

     IF(nd(k).EQ.0) CYCLE
        
     ihd = is(k) - 1
        
     DO  id = 1, nd(k)           
        ihd = ihd + 1
        ih(ihd) = 1        
     ENDDO

  ENDDO define_ihd



  
  !# exclude configurations from v2eb matrix

  k = 0
  include_cfg:DO  ir = 1, SIZE(u,dim=1)
               
     IF(ih(ir).EQ.0) CYCLE     ! exclude 

     k = k + 1 
     kk = 0
     DO  ic = 1, SIZE(u,dim=2)
        
        IF(ih(ic).EQ.0) CYCLE   ! exclude
        
        kk = kk + 1
        h(k, kk) = u(ir, ic)
        
     ENDDO
     
  ENDDO include_cfg
                        
  
  CALL diagonalize( nh2eb, nhx, h, u)       ! diagonalize now

  !
  CLOSE(nh2eb)
      

END PROGRAM h2eb

!


    !eof



!!%      subroutine cfin(nhf,lhf,ll,nmin,nmx,nol,is,l,ls,ncsmx,ncmx)
!!%      dimension nhf(1),lhf(1),ll(1),nmin(1),nmx(1),nol(1),is(1)
!!%    3 format(5i5)
!!%
!!%      read(15,3) l,ls
!!%      read(15,3) ncsmx
!!%
!!%      ncmx = 0
!!%      do  k = 1, ncsmx
!!%
!!%      read(15,3) nhf(k), lhf(k), ll(k), nmin(k), nmx(k)
!!%
!!%      nol(k) = nmx(k) - nmin(k) + 1
!!%      is(k)  = ncmx + 1
!!%      ncmx   = ncmx + nol(k)
!!%
!!%      end do
!!%
!!%      write(*,*) ' Total angular momentum      L = ', l
!!%      write(*,*) ' Total Spin                  S = ', ls
!!%      write(*,*) ' No of channel series    ncsmx = ', ncsmx
!!%      write(*,*) ' No of configurations     ncmx = ', ncmx
!!%      write(*,*) 
!!%
!!%      return
!!%      end
!#######################################################################
!!%C---------------------------------------------------------------
!!%C*  CALCULATE % CONTRIBUTION (PROBABILITY DENSITY) OF 
!!%C*  CONFIGURATION   |n1 l1 l2 > (SUM TO ALL N2) TO EIGENSTATES  
!!%C---------------------------------------------------------------
!!%C*  EXAMPLE :
!!%C*
!!%C*
!!%C*  In the CI procedure the eigenstates of atom with quantum 
!!%C*  number of LS is a sum over the 2e configurations  (n1 l1,n2 l2)
!!%C*  (allowed by the Clebsch-Gordan coefficients).
!!%C*  therefore:
!!%C*
!!%C*                 ____ 
!!%C*                 \ 
!!%C*   |He(gs):LS > = \     C(n1 l1,n2 l2)  * | LS : (n1 l1,n2 l2) > 
!!%C*                  /
!!%C*                 /___  
!!%C*          (to all n1 l1,n2 l2)  
!!%C*   
!!%C*     
!!%C*   We can calculate the contribution of the 1ss as :
!!%C*
!!%C*   1ss:      n1 = 1, l1 = s , l2 =s
!!%C*
!!%C*            ___  
!!%C*   p(1ss) = \   C(1s n2 s) * C(1s n2 s)
!!%C*            /__ 
!!%C*           (n2)
!!%C*          
!!%C*  Also
!!%C*  We can calculate the contribution of the 2ss as :
!!%C*
!!%C*   2ss:      n1 = 2, l1 = s , l2 = s
!!%C*            ___  
!!%C*
!!%C*   p(2ss) = \   C(2s n2 s) * C(2s n2 s)
!!%C*            /__ 
!!%C*           (n2)
!!%C*          
!!%C*
!!%C* 
!!%C*              Energy      p(1ss)          p(2ss)   ..... and so on 
!!%C*
!!%C*   HE(1S^2) :: xxx       p(1) * 100      p(2)*100   ..... p(i)*100 
!!%C*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!%C*     calculate for the first 10 configuration series
!!%C*  
!!%C*  for the  calculation of each of the | u(nln'l')|^2             
!!%C*            uncomment the 'c*1' comments 
!!%C* 
!!%C*

!*                output on "diagxxal.dat" :
!------------------------------------------------------------------
!*    Lo and Ls   ::  Orbital and Spin Angular Momenta
!*    NCSMX       ::  Number of Configuration Series
!*    NHF(K), LHF(K), LL(K), NLLMIN(K), NLLMAX(K),   K = 1, ncsmx
!*    ITRY        ::  Number of data sets in the file
!*    ND(K),      ::  Number of conf. included from each series
!*                     K=1,NCSMX 
!*    NHMX        ::  total number of config. in the eigenvectors
!*    NDTOTAL     ::  Number of states with output in the file, 
!*                    repeating energy eigenvalue and
!*                    energy eigenvector  (1 to nhmx)

!* EOF








