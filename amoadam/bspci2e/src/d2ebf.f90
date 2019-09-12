!!!   nhmx : nhmx > total number of configurations ( how much ?)
!!!   ncsi : # configuration series for the initial state
!!!   ncsf : # configuration series for the final   states
!!!............................

PROGRAM d2ebf
  !  use jl
  USE PRECISION, ONLY:dpk
  USE one_electron_data
  USE interchannel_dipole
  USE configuration
  USE utils, only: asciifile
!  USE ioroutines
  
  IMPLICIT NONE
  INTEGER j1, j2, ls
  COMMON/ang/j1,j2,ls
  !
  CHARACTER(len = 16), ALLOCATABLE, DIMENSION(:) :: inidat
  !
  !
  INTEGER                                :: i_fn
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: dmxz_b
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: e_inner, w_inner
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: e_outer, w_outer
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: e1i,e2i
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: e1f,e2f
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: w1i,w2i, w1f, w2f
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: e1
  INTEGER,      ALLOCATABLE,DIMENSION(:) :: n1,l1,l2 
  INTEGER                                :: i_in
  REAL(dpk)                              :: e_total
  INTEGER                                :: nd, dim_1
  !
  INTEGER                                :: ltotalMax
  !
  !

  INTEGER   lof, lsf, ncsmxf, loi,lsi, ncsmxi, ntotf
  INTEGER   nhmx, nbsp, ncsi, ncsf, nsi, nsf, ntoti, ne
  INTEGER   kf, ki, ir, ic, nr, nc, i, j, it, n
  INTEGER   LINITIAL, NMININ, NMAXI, LFINAL, NMINFIN, NMAXF, k
  INTEGER                                :: incx, incy
  INTEGER,  ALLOCATABLE, DIMENSION(:)    :: ndef
  INTEGER,   DIMENSION(4)                :: id
  REAL(dpk), DIMENSION(4)                :: ang
  REAL(dpk)                              :: alpha,beta 
  REAL(dpk)                              :: fli, flf, phase
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: dmxz, cifin
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: ciin, eni, enf
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: dpz, tt
  CHARACTER(len = 16)                    :: config
  CHARACTER(len = 16)                    :: cdiag
  CHARACTER(len = 16)                    :: outfile
  CHARACTER(len = 16)                    :: gauge
  CHARACTER(len = 100)                   :: argv 
  INTEGER                                :: icn
  INTEGER                                :: mode
  INTEGER                                :: NOUT, NWF2E, NH2E, NCFG, NHDMX2E
  integer                                :: ifile
  EXTERNAL                                  dgemv

! executable statements


      CALL GETARG(1, ARGV)
      READ(ARGV,*) LINITIAL         
      CALL GETARG(2, ARGV)
      READ(ARGV,*) LFINAL
      CALL GETARG(3, ARGV)
      READ(ARGV,*) NMAXF
      CALL GETARG(4, ARGV)
      READ(ARGV,*) GAUGE


      OPEN(9,file='inp/dmx2ebf.inp',status='old')
      READ(9,*) nminin,  nmaxi
      READ(9,*) nminfin
      CLOSE(9)


!.......

      NH2E    = 1
      NWF2E   = 3         
      NHDMX2E = 4
      NOUT    = 16
      NCFG    = 18

      ltotalmax = 4             ! It is connected with nm=5 (h1e)


!!!...........................


      OPEN(NOUT,FILE='out/d2ebf.out')

      WRITE(*,*) '# dmx2ebf::  bound-free dipole matrix elements calculation.'
      WRITE(*,*) '#'
      WRITE(*,*) '# dmx2ebf::                                        gauge = ', gauge
      WRITE(*,*) '# dmx2ebf::               initial state symmetry     lin = ', linitial
      WRITE(*,*) '# dmx2ebf::               final   state symmetry     lfn = ', lfinal

      MODE = 1
      IF(gauge.EQ.'V')    MODE = 0     



      ! get the number of channel series for initial and final states 

      !note:: implementation needs to change.

      nbsp = -1  
      CALL getnumberofchannels(ncfg, linitial, ncsi, nbsp)   !initial cfg-l.inp  file
      CALL getnumberofchannels(ncfg, lfinal,   ncsf, nbsp)   !final   HCFG-L.INP file
      CALL config_space(ncsi, ncsf)

      !
      CALL hfile(nwf2e,"dat","wf2e","bin",linitial)              !initial states

      READ(UNIT = NWF2E)   loi, lsi
      READ(UNIT = NWF2E)   ncsmxi
      READ(UNIT = NWF2E) ( nhfi(k),  k = 1, ncsmxi )
      READ(UNIT = NWF2E) ( lhfi(k),  k = 1, ncsmxi )
      READ(UNIT = NWF2E) ( lli(k),   k = 1, ncsmxi )
      READ(UNIT = NWF2E) ( nmini(k), k = 1, ncsmxi )
      READ(UNIT = NWF2E) ( nmxi(k),  k = 1, ncsmxi )
      READ(UNIT = NWF2E) ( ndi(k),   k = 1, ncsmxi )
      READ(UNIT = NWF2E)   ntoti
      READ(UNIT = NWF2E)   nsi
      
      CLOSE(NWF2E, STATUS='KEEP')


!      l1e_max = MAX(lhf)

      fli = 2 * loi + 1

      !#  some checks
      IF(LOI.NE.LINITIAL) THEN
         WRITE(NOUT,*) '# dmx2ebf::      wrong file is opened. '
         WRITE(NOUT,*) '# dmx2ebf::                    loi   = ', loi
         WRITE(NOUT,*) '# dmx2ebf::                 linitial = ', linitial
         STOP
      ENDIF

      IF(NMAXI.GT.NSI) THEN 
         WRITE(NOUT,*)'MAXIMUM NUMBER STATE OF INITIAL SYMMETRY CHANGED TO  nsi = ', nsi
         NMAXI = NSI
      ENDIF


      !#  final states cfg file

      lof = lfinal

      CALL cxfin(ncfg, nhff,lhff, llf, nminf, nmxf, nolf,&
                &isf, lof, lsf, ncsmxf, ntotf, ndf, idcs )

      
      FLF = 2*lof + 1


 !......... some checks
      IF(lsf.NE.lsi) THEN
         WRITE(*,*) '# dmx2ebf:: states of different parity symmetry'
         WRITE(*,*) '# dmx2ebf:;                              lsi = ',  lsi
         WRITE(*,*) '# dmx2ebf:;                              lsf = ',  lsf
         STOP
      ENDIF
      !
      j1 = lof             ! used by angls
      j2 = loi
      ls = lsi
      !

      phase = 1 - 2 * ls

      nhmx = MAX(ntoti, ntotf)

      WRITE(*,*) '# dmx2ebf::   max matrix dim,              nhmx = ', nhmx
      WRITE(*,*) '# dmx2ebf::   max cfg for initial state   ntoti = ', ntoti
      WRITE(*,*) '# dmx2ebf::   max cfg for final state     ntotf = ', ntotf
      WRITE(*,*) '# dmx2ebf::                               phase = ', phase


      !1-e dipole matrix elements

      !      
      ! vz(l,i,j):   1e-dipole matrix elements
      !


      ALLOCATE(   dz( nbsp, nbsp ) )               ! <pp|t|pb>
      ALLOCATE( dmxz( nhmx, nhmx ) )

      dmxz = 0.0_dpk


      CALL read_target_target_dipoles(nbsp, ltotalmax, mode)   !read <p|d|p>
      CALL read_target_basis_dipoles( nbsp, ltotalmax, mode)   !read <p|d|b>
      CALL read_target_basis_overlaps(nbsp, ltotalmax)         !read <p|b>

      WRITE(*,*) '# dmx2ebf:: one-electron data read.'
 
      !fill in dmxz

      kf = 1           
      dipole_1e_final_states: DO ir = 1, ncsmxf
         ki = 1         
         dipole_1e_initial_states: DO ic = 1, ncsmxi

            CALL angfc( lhff(ir),llf(ir),lhfi(ic),lli(ic),lof,loi,id,ang)

            WRITE(*,'(a60,2I5)')'(l_jp,l_j|d|l_ip,l_i) (J,I) = ', ic, ir 
            WRITE(*,'(5i5,E20.10)')llf(ir), lhff(ir), lli(ic), lhfi(ic),id(1), ang(1)
            WRITE(*,'(5i5,E20.10)')llf(ir), lhff(ir), lhfi(ic),lli(ic), id(2), ang(2)
            WRITE(*,'(5i5,E20.10)')lhff(ir), llf(ir), lli(ic), lhfi(ic),id(3), ang(3)
            WRITE(*,'(5i5,E20.10)')lhff(ir), llf(ir), lhfi(ic),lli(ic), id(4), ang(4)

         
            dz = 0.0_dpk
            SELECT CASE(idcs(ir))
            CASE(0)
               CALL evaluate_channel_pp_dipole_pp( ir, ic, id, ang, phase, mode)
            CASE(1)
               CALL evaluate_channel_pb_dipole_pp( ir, ic, id, ang, phase, mode)
            END SELECT
            WRITE(*,'(a60,E20.10)') ' d2e(J,I) = ', dz(ir,ic)

            CALL trnsmx(dmxz, nhmx, dz, nbsp, kf, ki, ndf(ir), ndi(ic))

            ki = ki + ndi(ic)

         ENDDO dipole_1e_initial_states
         kf = kf + ndf(ir)
      ENDDO dipole_1e_final_states

      CALL deallocate_one_electron_data() 


      !,,,,,,


      DEALLOCATE(dz)           !       


      ! normalize for equivalent orbitals by 1/sqrt(2)


      !initial states

      it = 1
      equivalent_orbitals_case:SELECT CASE( MOD( lhfi(1)+lli(1), 2) ) ! get initial state parity

      CASE(0)
         initial_channels: DO ic = 1, ncsmxi
            equiv_orbitals_initial:IF( (lhfi(ic).EQ.lli(ic)).AND.(nhfi(ic).EQ.nmini(ic)) ) THEN 

               dmxz(:,it) = dmxz(:,it) / dsqrt(2.0_dpk)

            ENDIF equiv_orbitals_initial
            it = it + ndi(ic) 
         END DO initial_channels
         
      CASE(1)              !,,,,,,,final states
         final_channels:DO ir = 1, ncsmxf

            SELECT CASE(idcs(ir))      ! only the fxd 2e-states
            CASE(0)
               equiv_orbitals_final:IF( (lhff(ir).EQ.llf(ir)).AND.(nhff(ir).EQ.nminf(ir)) ) THEN
                  dmxz(it,:) = dmxz(it,:)/dsqrt(2.0_dpk)
               ENDIF equiv_orbitals_final               
            END SELECT
            it = it + ndf(ir)

         END DO final_channels
         
      END SELECT equivalent_orbitals_case


      WRITE(*,'(a60)') 'renormalization of dipole for equivalent-orbitals done.'


      !allocations
      
      ALLOCATE(   eni(ntoti) )      
      ALLOCATE(  ciin(ntoti) )
      ALLOCATE(  cifin(ntotf, ncsmxf) )
      ALLOCATE(  tt(ntotf)   )

      !
      ! allocation for final symmetry
      !
      
      CALL hfile(nh2e,"dat","h2e","bin",lof)
      READ(NH2E) ne
      CLOSE(NH2E) 
      IF(NMAXF.GT.NE) THEN
         NMAXF = NE
         WRITE(NOUT,*) ' nmaxf changed to ne, nmaxf =', ne
      ENDIF
      ALLOCATE(enf(ne))
      ALLOCATE(ndef(ne))
      CLOSE(nh2e)

      !
      !end of allocation for final symmetry
      !



      CALL hfile(nwf2e,"dat","wf2e","bin",linitial)      !initial states

      READ(UNIT=NWF2E)   loi, lsi
      READ(UNIT=NWF2E)   ncsmxi
      READ(UNIT=NWF2E) ( nhfi(k),  k = 1, ncsmxi )
      READ(UNIT=NWF2E) ( lhfi(k),  k = 1, ncsmxi )
      READ(UNIT=NWF2E) ( lli(k),   k = 1, ncsmxi )
      READ(UNIT=NWF2E) ( nmini(k), k = 1, ncsmxi )
      READ(UNIT=NWF2E) ( nmxi(k),  k = 1, ncsmxi )
      READ(UNIT=NWF2E) ( ndi(k),   k = 1, ncsmxi )
      READ(UNIT=NWF2E)   ntoti
      READ(UNIT=NWF2E)   nsi

      ! dmx file
      CALL dmxfile(nhdmx2e, "dat", "dmx2ebf","bin", gauge, loi, lof) 

      WRITE(NOUT,*) '# gauge =  ', mode      !(1) vel
      WRITE(UNIT=NHDMX2E) mode
      WRITE(UNIT=NHDMX2E) NMININ, NMAXI
      WRITE(UNIT=NHDMX2E) NMINFIN, NMAXF


      !
      !,,,,,,,,,,,,,,,,
      !

!!%      CALL hfile(ifile,"dat","z2e","ascii",lof)     !final-symmetry amplitude file
!!%      !
!!%      READ(ifile, '(1e15.7,2x,i5)')  e_total, nd
!!%      !
!!%      ALLOCATE( n1( nd ) )
!!%      ALLOCATE( l1( nd ) )
!!%      ALLOCATE( l2( nd ) )
!!%      ALLOCATE( e1( nd ) )
!!%      ALLOCATE( e_outer( nd ) )
!!%      ALLOCATE( w_outer( nd, nd ) )
!!%      !
!!%      channel_q_numbers: DO i = 1, nd
!!%         READ(ifile, '(4i5,2E25.14)') icn, n1(i),l1(i),l2(i), e1, e_outer(i,i)
!!%      ENDDO channel_q_numbers
!!%      !
!!%      get_amplitudes: DO ic = 1, nd
!!%         READ(ifile, '(5e15.7)')   (w2(i,ic), i = 1, nd)     ! amplitude - matrix in wb-L.DAT
!!%      END DO get_amplitudes
!!%      CLOSE(ifile)
!!%
!!%      
!!%      WRITE(*,'(a60)')'read energies and amplitudes for inner orbitals'
!!%      WRITE(*,'(a60,i10)')' max angular q. number for inner orbitals lmax = ', dim_1e_basis
!!%     
!!%      DO icn = 1, ltotalmax + 1
!!%         CALL asciifile(ifile,lhfi(i_in),"wb-")
!!%         READ(ifile,*)  dim_1
!!%         IF(icn==0) THEN 
!!%            ALLOCATE(e_inner(dim_1,dim_1))
!!%            ALLOCATE(w_inner(dim_1,dim_1))            
!!%         ENDIF
!!%         DO i= 1, dim_1
!!%            READ(ifile,'(2E25.14)') e_inner(icn,i), w_inner(icn,i)
!!%         ENDDO
!!%      ENDDO
!!%      CLOSE(ifile)
!!%      
!!%      !
!!%      ALLOCATE(e1i(dim_1))
!!%      ALLOCATE(e2i(dim_1))
!!%      ALLOCATE(e1f(dim_1))
!!%      ALLOCATE(e2f(dim_1))
!!%      ALLOCATE(w1i(dim_1))
!!%      ALLOCATE(w2i(dim_1))
!!%      ALLOCATE(w2f(dim_1))
      !


      !
      !,,,,,,,,,,,,,,,,
      !      

      ALLOCATE(   dz_b( nbsp, nbsp ) )               ! <PP|T|PB>
      ALLOCATE( dmxz_b( nhmx, nhmx ) )                  

      dipole_2e_initial_states: DO ic = nminin, nmaxi            !

         READ(NWF2E)  eni(ic)
         READ(NWF2E)  (ciin(j), j = 1, ntoti)

         WRITE(*,*)  '  e_i = ', ic, eni(ic)
         WRITE(*,*)  '   nf = ', ne


         WRITE( UNIT=NHDMX2E) eni(ic)

         WRITE( NOUT, '(e14.6,2x,i5)') eni(ic)

         CALL hfile(nh2e,  "dat","h2e","bin",  lof)     !final-symmetry  h-file

         READ(NH2E) ne
            
         dipole_2e_final_states: DO ir = 1, nmaxf

            READ(NH2E) enf(ir), ndef(ir)
            READ(NH2E) ((cifin(n,j), n = 1, ntotf), j = 1, ndef(ir))
            WRITE(*, '(e14.6,2x,i5)') enf(ir), ndef(ir)


            WRITE(UNIT=NHDMX2E)      enf(ir), ndef(ir)

            !

            ALLOCATE( dpz( ndef(ir) ) )           ! D_g(j) = <g|D|f(j)>
            dpz = 0.0_dpk
            

            !...........,
            alpha = 1.0_dpk 
            beta  = 0.0_dpk
            incx  = 1
            incy  = 1
            channels_of_final_state: DO n = 1, ndef(ir)         !(solutions)


               !calculate correction terms since e_2(final) and c_2(final) is now known
               

!!%               e1i = 0.0_dpk
!!%               e2i = 0.0_dpk
!!%               e1f = 0.0_dpk
!!%               e2f = 0.0_dpk
!!%               w1i = 0.0_dpk
!!%               w2i = 0.0_dpk
!!%               w1f = 0.0_dpk
!!%               w2f = 0.0_dpk
!!%               dz_b = 0.0_dpk        !correction depends on e2 = E-e1
!!%               kf = 1
!!%               loop_over_channel_final_states: DO i_fn = 1, ncsmxf
!!%                  ki = 1               
!!%                  loop_over_cfg_initial_states: DO i_in = 1, ncsmxi
!!%
!!%                     ! < P(:,l1_f),P(:,l2_f)|D|P(:,l1_i),B(nb,l2_f)>
!!%
!!%                     IF(idcs(i_fn) == 1) THEN
!!%
!!%
!!%
!!%                        ! target channels /l1,l2> energies and amplitudes
!!%
!!%                        e1i = e_inner(lhfi(i_in) + 1, : )
!!%                        e2i = e_inner( lli(i_in) + 1, : )
!!%                        e1f = e_inner(lhff(i_fn) + 1, : )
!!%                        e2f = enf(ir) - e1f
!!%                        !
!!%                        w1i = w_inner(lhfi(i_in) + 1, : )
!!%                        w2i = w_inner( lli(i_in) + 1, : )
!!%                        w1f = w_inner(lhff(i_fn) + 1, : )
!!%                        w2f = w_outer(i_fn, n )
!!%                        !
!!%
!!%                        ! determine angular factors: ang(4)/phase(4)
!!%                        
!!%                        CALL angfc( lhff(i_fn),llf(i_fn),lhfi(i_in),lli(i_in),lof,loi,id,ang)
!!%                        CALL evaluate_channel_dipole_corrections( i_fn, i_in, id, ang, phase, mode,&
!!%                             e1i, e2i,e1f,e2f, w1i,w2i,w1f,w2f)
!!%                        CALL add_dipole_corrections(dmxz_b, nhmx, dz_b, nbsp, &
!!%                             kf, ki, ndf(i_fn), ndi(i_in))
!!%                     ENDIF
!!%                     ki = ki + ndi(i_in)
!!%                     
!!%                  ENDDO loop_over_cfg_initial_states
!!%                  kf = kf + ndf(i_fn)               
!!%               ENDDO loop_over_channel_final_states
               

               !               dmxz = dmxz + dmxz_b

                  !  calculate the final dipole

                  CALL dgemv('n', ntotf, ntoti, alpha, dmxz, nhmx, ciin, incx, beta, tt, incy)

                  dpz(n) = 0.0_dpk
                  dpz(n) = DOT_PRODUCT(cifin(:,n), tt)

               END DO channels_of_final_state


               WRITE(UNIT=NHDMX2E)    ( dpz(j), j = 1, ndef(ir) )
               WRITE(NOUT,'(4e14.6)') ( dpz(j), j = 1, ndef(ir) )

               DEALLOCATE(dpz)
            END DO dipole_2e_final_states

            CLOSE(NH2E)
         END DO dipole_2e_initial_states

      
!!!.................................
      
      DEALLOCATE(eni)
      DEALLOCATE(enf)
      DEALLOCATE(ciin)
      DEALLOCATE(cifin)
      DEALLOCATE(ndef)
      DEALLOCATE(tt)

      CLOSE(NHDMX2E)
      CLOSE(NWF2E)
      CLOSE(NH2E)
!!!.....................................

1     FORMAT(/2x,'the initial state of the transition --')
2     FORMAT(/2x,'the final state of the transition --')
3     FORMAT(2x,'total l and total s=',2i3,4x,'# of config. is ',i4)
4     FORMAT(/2x,'energy eigenvalues of initial & final states')
5     FORMAT(2x,8e15.7)
6     FORMAT(/2x,'os values: length - top & velocity - bottom.   {lf is&
           & included. for absorp. lf=2*lof+1. for emiss.("-" sign) lf=2*loi+1&
           & }' /2x,'row - final states & column - initial states')
7     FORMAT(2x)
8     FORMAT(/2x,'lhff,llf,lhfi,lli,id(i),ang(i) --'/)
9     FORMAT(2x,4i2,3x,4i2,4e12.4)
11    FORMAT(/2x,'excitation energy in ryd.')
!!!.................................
      
    END PROGRAM d2ebf
    !
    !
    !
    !
    SUBROUTINE trnsmx(hmx,nhmx,xm,nbsp,ihr,ihc,ir,ic)
      !
      USE PRECISION, ONLY:dpk
      !
      IMPLICIT NONE
      REAL(DPK),        DIMENSION(nhmx,nhmx) ::  hmx
      INTEGER                                :: nhmx
      REAL(DPK),        DIMENSION(nbsp,nbsp) ::  xm
      INTEGER                                :: nbsp
      INTEGER                                :: ihr
      INTEGER                                :: ihc
      INTEGER                                :: ir
      INTEGER                                :: ic
      !
      INTEGER                                :: irf, icf, kr, k, kc, kk
      !

      
      irf = ihr + ir - 1
      icf = ihc + ic - 1
      
!!!      WRITE(*,*) ihr,irf, ihc,icf
      DO  k = ihr, irf
         DO  kk = ihc, icf
            kr        = k - ihr + 1
            kc        = kk- ihc + 1
            hmx(k, kk) = xm(kr,kc)
         ENDDO
      ENDDO

    END SUBROUTINE trnsmx
    
    ! EOF
