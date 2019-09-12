MODULE jl
  IMPLICIT NONE
  INTEGER j1, j2, ls
END MODULE jl

MODULE wf_files
  USE PRECISION, ONLY:dpk
  USE ioroutines
  
  IMPLICIT NONE
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:) :: vz
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:) :: vzba
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:) :: sp_a

CONTAINS
  !
  !
  !
  SUBROUTINE read_dp(nsx, ltotalMax,mode)
    !
    USE utils, ONLY: dfile, datafile
    !
    IMPLICIT NONE
    !arg!
    integer          :: nsx
    integer          :: ltotalmax
    integer          :: mode
    !loc!
    CHARACTER(len=6) :: gauge
    !
    INTEGER          :: m, n, li, lf, nwfi, nwff, ir, ic
    !io!    
    INTEGER         :: nout, nbin
    INTEGER         :: i,j

!!!......................
    WRITE(*, '(a60)') 'subroutine::read_dp in.'
    WRITE(*, '(a60)') 'dipole matrix elements for target states.'

    !
    nbin   = 2
    nout   = 16


    IF(mode==0) THEN
       gauge = 'v'
    ELSE IF(mode==1) THEN
       gauge = 'l'
    ELSE IF(mode==2) THEN
       gauge = 'a'
    ENDIF
       
       !
    WRITE(*,'(a60,i10)') ' lmax = ' , ltotalMax
    WRITE(*,'(a60,i10)') ' nof basis functions nb = ', nsx
    WRITE(*,'(a60,a10)') '                  gauge = ', gauge


    ALLOCATE(vz(ltotalMax, nsx, nsx))             ! <P_nl|T|P_ml+1>, l=0,...,lmax-1

    DO i = 1, ltotalMax  

!       CALL HDMX1EFILE(NDMX1E, i-1, mode)

       CALL dfile(nbin,i-1,i,'d1e-',gauge)
       READ(nbin) mode
       READ(nbin) li, lf, nwfi, nwff
       READ(nbin) ( ( vz(i,ir,ic), ir = 1, nwff), ic = 1, nwfi)       
       CLOSE(nbin)
         
       WRITE(*,'(a60,5i5)') '(mode,la,lb,na,nb) = ', mode, li,lf,nwfi,nwff
    END DO

    WRITE(*,'(a60)') 'dipoles  <P|T|P> read.'

    ALLOCATE( vzba(2*ltotalMax, nsx, nsx) )

    j = 0
    get_overlap_dipoles: DO i = 1, ltotalMax   

!       CALL RWFBFILE(nbin, i-1, i, mode)
       CALL dfile(nbin,i-1,i,'overlap-d-',gauge)

       READ(nbin) mode
       READ(nbin) li,lf,nwfi,nwff
       READ(nbin) ((vzba(i + j, ir, ic), ir = 1, nwff), ic = 1, nwfi)

       CLOSE(nbin)


!       CALL RWFBFILE(NRWFB, i, i-1, MODE)
       CALL dfile(nbin,i,i-1,'overlap-d-',gauge)

       READ(nbin) mode
       READ(nbin) li,lf,nwfi,nwff
       READ(nbin) ((vzba(i+j+1, ir, ic), ir = 1, nwff), ic = 1, nwfi)

       j = j + 1

       CLOSE(nbin)
       
       WRITE(*,'(a60,5i5)') '(mode,la,lb,na,nb) = ', mode, li,lf,nwfi,nwff
    END DO get_overlap_dipoles

    WRITE(*,'(a60)') 'overlap dipoles  <P|T|B> read.'



    ALLOCATE(sp_a(ltotalMax, nsx, nsx))          ! <P_nl|B_j> overlaps

    get_overlaps:DO i = 1, ltotalMax

!       CALL RBBFILE(nbin, i-1)
       CALL datafile(nbin,i-1,'overlap-1-') 

       READ(nbin) li,nwfi,nwff
       READ(nbin) ((sp_a(li+1,ir,ic), ir=1,nwff), ic=1,nwfi)
       CLOSE(nbin)

       WRITE(*,'(a60,3i5)') '(la,na,j) = ', li,nwfi,nwff
    END DO get_overlaps

    WRITE(*,'(a60)') 'overlap <P|B> read.'


    WRITE(*, '(a60)') 'subroutine::read_dp out.'
  END SUBROUTINE read_dp

END MODULE wf_files
!S
!S
!S
!S
MODULE configuration
  IMPLICIT NONE                       ! n1,   l1,   l2,  n2_min, n2_max , n2_max-n2_min+1
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhfi,lhfi,lli,nmini,nmxi, ndi
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhff,lhff,llf,nminf,nmxf, ndf
  INTEGER, ALLOCATABLE, DIMENSION(:) :: isf, nolf, idcs
CONTAINS
!sub
  SUBROUTINE config_space(ncsi, ncsf)  
    !
    IMPLICIT NONE
    !
    INTEGER              :: ncsi        ! nof channels for the initial state
    INTEGER              :: ncsf        ! nof channels for the final   state


    ALLOCATE( nhfi (ncsi) )           !initial state
    ALLOCATE( lhfi (ncsi) )
    ALLOCATE( lli  (ncsi) )
    ALLOCATE( nmini(ncsi) )
    ALLOCATE( nmxi (ncsi) )
    ALLOCATE( ndi  (ncsi) )
    ALLOCATE( nhff (ncsf) )            !final state
    ALLOCATE( lhff (ncsf) )
    ALLOCATE( llf  (ncsf) )
    ALLOCATE( nminf(ncsf) )
    ALLOCATE( nmxf (ncsf) )
    ALLOCATE( ndf  (ncsf) )    
    ALLOCATE( nolf (ncsf) )            ! 
    ALLOCATE( isf  (ncsf) )            ! singlet/triplet
    ALLOCATE( idcs (ncsf) )            ! fxd(0)/free(1) states
    !
  END SUBROUTINE config_space
  !
END MODULE configuration


MODULE channel_dipole
  !
  USE wf_files
  IMPLICIT NONE
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: dz
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: dz_b

CONTAINS
  !
  ! get_dipole_corrections
  !
  SUBROUTINE evaluate_channel_dipole_corrections(ir, ic, id, ang, phase, mode, &
                                                e_outer,w_outer,e_inner,w_inner)
    !
    USE PRECISION, ONLY: dpk
    !
    USE configuration
    !
    IMPLICIT NONE
    !
    INTEGER                   :: ir             ! final   channel |ir> = |P(ni,li) * B(ni',li') > 
    INTEGER                   :: ic             ! initial channel |ic> = |P(nj,lj) * B(nj',lj') >  
    INTEGER,    DIMENSION(4)  :: id             ! info for li-->li \pm 1
    REAL(dpk),  DIMENSION(4)  :: ang
    REAL(dpk)                 :: phase          
    INTEGER                   :: mode           ! gauge
    REAL(dpk), DIMENSION(:,:) :: e_outer
    REAL(dpk), DIMENSION(:,:) :: w_outer
    REAL(dpk), DIMENSION(:,:) :: e_inner
    REAL(dpk), DIMENSION(:,:) :: w_inner
    !
    INTEGER                  :: i, idn
    INTEGER                  :: nr, nc
    INTEGER                  :: mr, mc
    REAL(dpk)                :: de,sb, rho_b
    !
    INTEGER                    :: ni,nj,nip,njp
    INTEGER                    :: nbi,nbf, nbi_1st, nbf_1st
    INTEGER                    :: gg
    

    !
    !
    !    vzba(i+j,   j,n) = < P_n,(l, ) |D| B_j, (l+1)>
    !    vzba(i+j+1, j,n) = < P_n,(l+1),|D| B_j, (l)  >
    !
    !

!if  < a b |D| c d > = <P(na,la)P(nb,lb)|D|P(nc,lc)B(nd,ld)>

  gg = 1
  IF(mode==0) gg = -1

  !
  ! channel of initial state : ic
  ! channel of final   state : ir
  !
  !  d(nj,ni) = < P(njp,ljp) * P(nj,lj) |D|P(nip,lip) * P(ni,li) >
  !
  !
  !

  njp     = nhff(ir)      ! final state 'inner' orbital (ljp)
  nj      = mr            ! final state 'outer' orbital (lj)
  nbf_1st = nminf(ir)     ! 
  nbf     = ndf(ir)       ! basis dimension for the outer orbital
  !
  nip     = nhfi(ic)      ! final state 'inner' orbital (lip)
  ni      = mc            ! final state 'outer' orbital (li)
  nbi     = ndi(ic)       !
  nbi_1st = nmini(ic)     ! basis dimension for the outer orbital

  !
  ! 1->01
  ! 2->10
  ! 3->12
  ! 4->21
!  
  antisymmetrization_terms: DO i = 3, 4

     IF(id(i).NE.0) THEN ! |l-l'| = 1
        
        SELECT CASE(i)

        CASE(3)         ! li=ljp ---> <Pi|Pj'> <Pi'|t|Bj>

             !rho_b = (kb-1)/dt
 
             channel_initial_3: DO  nc = 1, nbi

                ni = nbi_1st + nc - 1

                overlap_3:IF(ni.EQ.njp) THEN     !<i|j'> = \delta(i-j')

                   nj = nbf_1st + nbf - 1

                   de = ehf(nip) - e2(nj)         !de = ehf(nip) - e2(nbf)
                   sb = wb_f(nj) * wb_i(nip)      !sb = wb_f(nbf) * wb_i(nip) 

                   
                   IF(mode.EQ.0) THEN
                      dz_b(nbf,nc)=  dz(nbf,nc) - phase * ang(i) * sb *rho_b/de
                   ELSE IF(mode.EQ.1) THEN
                      dz_b(nbf,nc)=  dz(nbf,nc) + phase * ang(i) * sb *rho_b
                   ENDIF

                ENDIF overlap_3

             ENDDO channel_initial_3
             !
             !
             !             
          CASE(4)          ! li'=lj' ---> <Pi'|Pj'> <Pi|t|Bj>


             overlap_4:IF(njp.EQ.nip) THEN 

                   nj = nbf_1st + nbf - 1

                   initial_channels_21:DO  nc = 1, nbi

                      ni = nbi_1st + nc - 1

                      de = ehf(ni) - e2(nj)         !de = ehf(nip) - e2(nbf)
                      sb = wb_f(nj) * wb_i(ni)      !sb = wb_f(nbf) * wb_i(nip) 

                   
                      IF(mode.EQ.0) THEN
                         dz_b(nbf,nc)=  dz(nbf,nc) - phase * ang(i) * sb *rho_b/de
                      ELSE IF(mode.EQ.1) THEN
                         dz_b(nbf,nc)=  dz(nbf,nc) + phase * ang(i) * sb *rho_b
                      ENDIF
                                           
                      dz_b(nr,nc) = dz_b(nr,nc) + ang(i) * sb * rho_b/de

                   ENDDO initial_channels_21
                ENDIF overlap_4


             END SELECT

          ENDIF
       
       ENDDO antisymmetrization_terms

     END SUBROUTINE evaluate_channel_dipole_corrections
!!!#######################################################################
                 
  SUBROUTINE evaluate_channel_pp_dipole_pp(ir, ic, id, ang, phase, mode) 
  !mod!
  USE PRECISION, ONLY: dpk
  USE configuration
  !
  IMPLICIT NONE
  !arg!
  INTEGER                    :: ir
  INTEGER                    :: ic
  INTEGER,      DIMENSION(4) :: id
  REAL(dpk),    DIMENSION(4) :: ang
  REAL(dpk)                  :: phase
  INTEGER                    :: mode
  INTEGER                    :: i
  INTEGER                    :: nr, nc, mr, mc, idn
  INTEGER                    :: ni,nj,nip,njp
  INTEGER                    :: nbi,nbf, nbi_1st, nbf_1st
  integer                    :: gg
  !exe!

  
  gg = 1
  IF(mode==0) gg = -1

  !
  ! channel of initial state : ic
  ! channel of final   state : ir
  !
  !  d(nj,ni) = < P(njp,ljp) * P(nj,lj) |D|P(nip,lip) * P(ni,li) >
  !
  !
  !

  njp     = nhff(ir)      ! final state 'inner' orbital (ljp)
  nj      = mr            ! final state 'outer' orbital (lj)
  nbf_1st = nminf(ir)     ! 
  nbf     = ndf(ir)       ! basis dimension for the outer orbital
  !
  nip     = nhfi(ic)      ! final state 'inner' orbital (lip)
  ni      = mc            ! final state 'outer' orbital (li)
  nbi     = ndi(ic)       !
  nbi_1st = nmini(ic)     ! basis dimension for the outer orbital

!
! 1->01
! 2->10
! 3->12
! 4->21
!  

  antisymmetrization_terms:DO  i = 1, 4

     IF(id(i).NE.0) THEN      ! |l-l'| = 1

        SELECT CASE(i)

        CASE(1)         ! li=lj ---> <i|j> <i'|t|j'>

           channel_final_1:DO nr = 1, nbf
              nj = nbf_1st + nr - 1
              channel_initial_1:DO nc = 1, nbi                 
                 ni = nbi_1st + nc - 1
                 
!,,,
                 overlap_1:IF(ni.EQ.nj) THEN    ! <i|j> = \delta(i,j) (overlap integral)
                    

                    !<j'|t|i'>  
                    IF(MOD(id(i),2).EQ.0) THEN      ! ljp = lip + 1  
                       
                       idn = id(i) / 2                       
                       dz(nr,nc)= dz(nr,nc)+ gg * ang(i) * vz(idn,nip,njp)
                       
                    ELSE !id(i) = 1,3,5,..          !lip = ljp + 1

                       idn=(id(i)+1)/2
                       dz(nr,nc)= dz(nr,nc)+ang(1) * vz(idn,njp,nip)

                    END IF

                 ENDIF overlap_1
                 !
              END DO channel_initial_1
           END DO channel_final_1
                            
           
        CASE(2) ! lj=lip ---> <i'|j> <i|t|j'>

           channel_final_2:DO nr = 1, nbf
              
              nj = nbf_1st + nr - 1
              
              overlap_2:IF(nip.EQ.nj) THEN       ! <i'|j> = delta(i'-j)
                 
                 channel_initial_2:DO nc = 1, nbi

                    ni = nbi_1st + nc - 1
                    
                    !<i|t|j'>
                    IF(MOD(id(i),2).EQ.0) THEN   ! ljp = li + 1
                       idn=id(i)/2
                       dz(nr,nc) = dz(nr,nc) + gg * phase * ang(i) * vz(idn,ni,njp)
                    ELSE                         ! li = ljp + 1
                       idn=(id(i)+1)/2
                       dz(nr,nc) = dz(nr,nc) +      phase * ang(i) * vz(idn,njp,ni)
                    END IF
                    
                    
                 END DO channel_initial_2
                 
              ENDIF overlap_2
              
           END DO channel_final_2
           !           

        CASE(3)       ! li=lj ---> <i|j'> <i'|t|j>

           channel_initial_3:DO nc = 1, nbi

              ni = nbi_1st + nc - 1

              overlap_3:IF(ni.EQ.njp) THEN !<i|j'> = \delta(i-j')
                 
                 channel_final_3:DO  nr = 1, nbf

                    nj = nbf_1st + nr - 1
                    
                    !<i'|t|j>
                    IF(MOD(id(i),2).EQ.0) THEN      ! lip = lj + 1
                       idn=id(i)/2
                       dz(nr,nc) = dz(nr,nc) + gg *phase*ang(i) * vz(idn,nip,nj)
                    ELSE                            ! lj = lip + 1
                       idn = (id(i)+1) / 2
                       dz(nr,nc) = dz(nr,nc) + phase*ang(3) * vz(idn,nj,nip)
                    END IF

                 ENDDO channel_final_3
                 
              ENDIF overlap_3
              
           ENDDO channel_initial_3
           

        CASE(4)          ! li'=lj' ---> <i'|j'> <i|t|j>
           overlap_4:IF(njp.EQ.nip) THEN  !<i'|j'> = delta(i'-j')
!                           

              !<i|t|j>
              channel_final_4:DO  nr = 1, nbf
                 nj = nbf_1st + nr - 1                 
                 channel_initial_4:DO   nc = 1, nbi
                    ni = nbi_1st + nc - 1

                    
                    IF(MOD(id(i),2).EQ.0) THEN    ! li = lj + 1
                       idn=id(i)/2
                       dz(nr,nc) = dz(nr,nc) + gg * ang(i)*vz(idn,ni,nj)
                    ELSE                          ! lj = li - 1
                       idn=(id(i)+1)/2
                       dz(nr,nc) = dz(nr,nc) + ang(i) * vz(idn,nj,ni)
                    END IF
                    
                 ENDDO channel_initial_4
              ENDDO channel_final_4              
           ENDIF overlap_4
           
        END SELECT
     ENDIF
     
  ENDDO antisymmetrization_terms 
  
  
END SUBROUTINE evaluate_channel_dipole

  !S
  !S
  !S
  SUBROUTINE evaluate_channel_pb_dipole_pp(ir, ic, id, ang, phase, mode) 
    !
    USE PRECISION, only: dpk
    !
    USE configuration
    !
    IMPLICIT NONE
    INTEGER                 :: ir
    INTEGER                 :: ic
    INTEGER,   DIMENSION(4) :: id          
    REAL(dpk), DIMENSION(4) :: ang
    REAL(dpk)               :: phase
    INTEGER                 :: mode
    !
    INTEGER                 :: i
    INTEGER                 :: nr, nc, mr, mc, idn
    INTEGER                 :: ni,nj,nip,njp
    INTEGER                 :: li, lip
    INTEGER                 :: nbi,nbf, nbi_1st, nbf_1st
    INTEGER                 :: gg
    !exe!

  
    gg = 1
    IF(mode==0) gg = -1

  !
  ! channel of initial state : ic
  ! channel of final   state : ir
  !
  !  d(nj,ni) = < P(njp,ljp) * P(nj,lj) |D|P(nip,lip) * B(ni,li) >
  !
  !
  !lj = lli(ic)
  !

  njp     = nhff(ir)      ! final state 'inner' orbital (ljp)
  nj      = mr            ! final state 'outer' orbital (lj)  
  nbf_1st = nminf(ir)     ! 
  nbf     = ndf(ir)       ! basis dimension for the outer orbital
  !
  nip     = nhfi(ic)      ! final state 'inner' orbital (lip)
  ni      = mc            ! final state 'outer' orbital (li)
  li      = lli(ic)       !
  lip     = lhfi(ic)       !
  nbi     = ndi(ic)       !
  nbi_1st = nmini(ic)     ! basis dimension for the outer orbital



!
! 1->01
! 2->10
! 3->12
! 4->21

    !

    antisymmetrization_terms: DO  i = 1, 4

       IF(id(i).NE.0) THEN    ! |l-l'| = 1

          SELECT CASE(i) 

          CASE(1)               ! li=lj ---> <Pi|Bj> <Pi'|t|Pj'>

             channel_final_1:DO nr = 1, nbf
                nj = nbf_1st + nr - 1
                channel_initial_1:DO nc = 1, nbi
                   ni = nbi_1st + nc - 1

                   IF(MOD(id(1),2).EQ.0) THEN ! ljp = lip + 1  
                      idn = id(1)/2
                      dz(nr,nc) = dz(nr,nc)+ gg * ang(1)* vz(idn,nip,njp) * sp_a(li+1,nj,ni)
                   ELSE                       ! lip = ljp + 1  
                      idn=(id(1)+1)/2
                      dz(nr,nc)= dz(nr,nc)+ang(1) * vz(idn,njp,nip) * sp_a(li+1,nj,ni)
                   END IF 
                   
                END DO channel_initial_1
             END DO channel_final_1
             

          CASE(2) ! lj=lip ---> <Pi'|Bj> <Pi|t|Pj'>

             channel_final_2: DO nr = 1, nbf
                nj = nbf_1st + nr - 1
                channel_initial_2: DO nc = 1, nbi
                   ni = nbi_1st + nc - 1


                   IF(MOD(id(i),2).EQ.0) THEN   ! ljp = li + 1
                      idn=id(2)/2
                      dz(nr,nc) = dz(nr,nc) + gg * phase * ang(i)* vz(idn,ni,njp) * sp_a(lip+1,nj,nip)
                   ELSE                         ! li = ljp + 1
                      idn=(id(i)+1)/2
                      dz(nr,nc) = dz(nr,nc) +      phase * ang(i)* vz(idn,njp,ni) * sp_a(lip+1,nj,nip)
                   END IF

                END DO channel_initial_2
             END DO channel_final_2
!
!
!
          CASE(3)        ! li=ljp ---> <Pi|Pj'> <Pi'|t|Bj>

             channel_initial:DO  nc = 1, nbi
                ni = nbi_1st + nc - 1
                overlap_3:IF(ni.EQ.njp) THEN       !<i|j'> = \delta(i-j')
                   channel_final:DO  nr = 1, nbf

                      nj = nbf_1st + nr - 1

                      idn = id(i)
                      dz(nr,nc)= dz(nr,nc)+ phase * ang(i) * vzba(idn,nj,nip)

                   ENDDO channel_final
                ENDIF overlap_3
             ENDDO channel_initial
             !
             !
             !
                CASE(4)          ! li'=lj' ---> <Pi'|Pj'> <Pi|t|Bj>

                   overlap_4:IF(njp.EQ.nip) THEN   !<i'|j'> = delta(i'-j')

                      channel_final_4:DO  nr=1,nbf
                         
                         nj = nbf_1st + nr - 1
                         channel_initial_4:DO  nc = 1, nbi
                            
                            ni = nbi_1st + nc - 1
                            idn = id(i)
                            dz(nr,nc) = dz(nr,nc)+ ang(i) * vzba(idn,nj,ni)

                         ENDDO
                      ENDDO channel_final_4
                   ENDIF overlap_4


                END SELECT
             
             ENDIF

          ENDDO antisymmetrization_terms

        END SUBROUTINE evaluate_channel_pb_dipole_pp
!####################################################
                       !
      END MODULE channel_dipole
                     !eof         
                     
