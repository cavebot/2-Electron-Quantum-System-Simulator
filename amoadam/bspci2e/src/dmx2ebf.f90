!!!        
!!!         laan/iesl/2001
!!!
!!!         06.09.2001: 
!!!                    Reconstructing the code. Mainly I/O facilities 
!!!                    have been added for a more automatic running
!!!                    of the whole package   
!!!
!!!           note  (0) old name was   dpbf.f90 (dpbf_mod.f90,dpbf_sub77.f)  
!!!           note  (1) The file opened by the Routine HWF2EFILE   
!!!                        cannot be read at once. The file closes  
!!!                        during the execution of the programm and
!!!                        and message error comes that do not exist
!!!                        The correspondent programm dmx2e.f for the  
!!!                        fxd code (compiled with pgf77) not such   
!!!                        problem appears.
!!!                        Temporarily, i close and open again this
!!!                        file.
!!!          
!!!  
!!!
!!!  purpose:: Calculates un-normalized 2-e DME from a bound state | b > to 
!!!  a multichannel continuum (single or double) | c >
!!!  
!!!  | b > ----> | c >       
!!!
!!!  ! b >  2-e bound state calculated using the BSPCI2E code (fxd)
!!!  | c >  continuum state calculated on the mixed basis states (free) 
!!!
!!!
!!!   nhmx : nhmx > total number of configurations ( how much ?)
!!!   ncsi : # configuration series for the initial state
!!!   ncsf : # configuration series for the final   states
!!!............................

PROGRAM dmx2ebf
  !  use jl
  USE PRECISION, ONLY:dpk
  USE wf_files
  USE dz_value
  USE configuration
  USE ioroutines
  
  IMPLICIT NONE
  INTEGER j1, j2, ls
  COMMON/ang/j1,j2,ls
  !
  CHARACTER(len = 16)  config, cdiag, outfile, gauge
  CHARACTER(len = 16), ALLOCATABLE, DIMENSION(:) :: inidat
  !
  INTEGER   lof, lsf, ncsmxf, loi,lsi, ncsmxi, ntotf
  INTEGER   nhmx, nbsp, ncsi, ncsf, nsi, nsf, ntoti, ne
  INTEGER   kf, ki, ir, ic, nr, nc, i, j, it, n, ltotalMax
  INTEGER   LINITIAL, NMININ, NMAXI, LFINAL, NMINFIN, NMAXF, k
  INTEGER   MODE, NOUT, NWF2E, NH2E, NCFG, NHDMX2E
  INTEGER   incx, incy
  INTEGER,  ALLOCATABLE, DIMENSION(:)    :: ndef
  INTEGER,   DIMENSION(4)                :: id
  REAL(dpk), DIMENSION(4)                :: ang
  REAL(dpk)                                  alpha,beta 
  REAL(dpk)                                 fli, flf, phase
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: dmxz, cifin
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: ciin, eni, enf
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: dpz, tt
  CHARACTER(LEN=100) ARGV 
  EXTERNAL dgemv

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
      READ(9,*) NMININ,  NMAXI
      READ(9,*) NMINFIN 
      CLOSE(9)


!.......

      NH2E    = 1
      NWF2E   = 3         
      NHDMX2E = 4
      NOUT    = 16
      NCFG    = 18

      LTOTALMAX = 4             ! It is connected with nm=5 (h1e)


!!!...........................


      OPEN(NOUT,FILE='out/dmx2ebf.out')

      WRITE(*,*) '# dmx2ebf::  bound-free dipole matrix elements calculation.'
      WRITE(*,*) '#'
      WRITE(*,*) '# dmx2ebf::                                        gauge = ', gauge
      WRITE(*,*) '# dmx2ebf::               initial state symmetry     lin = ', linitial
      WRITE(*,*) '# dmx2ebf::               final   state symmetry     lfn = ', lfinal

      MODE = 1
      IF(gauge.EQ.'V')    MODE = 0     



      ! get the number of chanell series for initial and final states 

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


      ALLOCATE( dmxz( nhmx, nhmx ) )
      ALLOCATE(   dz( nbsp, nbsp ) )


      dmxz = 0.0_dpk
  
      CALL read_dp(nbsp, ltotalMax, mode)          ! read dmx1e


      WRITE(*,*) '# dmx2ebf:: dmx-1e read.'

 
      !fill in dmxz

      kf = 1           
      dipole_1e_final_states: DO ir = 1, ncsmxf
         
         ki = 1
         
         dipole_1e_initial_states: DO ic = 1, ncsmxi

            CALL angfc( lhff(ir),llf(ir),lhfi(ic),lli(ic),lof,loi,id,ang)
            
            SELECT CASE(idcs(ir))

            CASE(0)
               dz = 0.0_dpk              ! <P|T|P>
               CALL cal_nlnl( ir, ic, id, ang, phase, mode)                  
            CASE(1)
               dz = 0.0_dpk              ! <P|T|B>, <P|B>
               CALL cal_nlsp( ir, ic, id, ang, phase, mode)
            END SELECT

            !    put      dz( nb x nb ) matrix
            !
            !   at the right place in the large 
            !
            !     dmxz(nhmx,nhmx) matrix
            !

            CALL trnsmx(dmxz, nhmx, dz, nbsp, kf, ki, ndf(ir), ndi(ic))

            ki = ki + ndi(ic)

         ENDDO dipole_1e_initial_states

         kf = kf + ndf(ir)

      ENDDO dipole_1e_final_states

      DEALLOCATE(dz)   
      DEALLOCATE(vz)
      DEALLOCATE(vzba)
      DEALLOCATE(sp_a)
      

      ! normalize for equivalent orbitals by 1/sqrt(2)


!initial states

      it = 1
      SELECT CASE( mod( lhfi(1)+lli(1), 2) ) ! initial state parity

      CASE(0)
         check_initial: do ic = 1, ncsmxi
            equiv_orbitals_initial:IF( (lhfi(ic).EQ.lli(ic)).AND.(nhfi(ic).EQ.nmini(ic)) ) THEN 

               DO nr = 1, nhmx
                  dmxz(nr,it) = dmxz(nr,it) / dsqrt(2.0d0)
               END DO

            ENDIF equiv_orbitals_initial

            it = it + ndi(ic) 
         END DO check_initial

      CASE(1)

!final states
         check_final:DO ir = 1, ncsmxf

            SELECT CASE(idcs(ir))
            CASE(0)

               equiv_orbitals_final:IF( (lhff(ir).EQ.llf(ir)).AND.(nhff(ir).EQ.nminf(ir)) ) THEN
                  
                  DO nc = 1, nhmx
                     dmxz(it,nc) = dmxz(it,nc)/dsqrt(2.0d0)
                  END DO
                  
               ENDIF equiv_orbitals_final
               
            END SELECT
            it = it + ndf(ir)
            
         END DO check_final

      END SELECT


      WRITE(*,*) '# dmx2ebf::  normalization for equivalent orbitals done.'


      !allocations

      ALLOCATE(   eni(ntoti) )      
      ALLOCATE(  ciin(ntoti) )
      ALLOCATE(  cifin(ntotf, ncsmxf) )
      ALLOCATE(  tt(ntotf)   )

                                !  read final (continuum) states config. coef 

      CALL hfile(nh2e,"dat","h2e","bin",lof)

      READ(NH2E) ne
      CLOSE(NH2E) 


      IF(NMAXF.GT.NE) THEN
         NMAXF = NE
         WRITE(NOUT,*) ' MAXIMUM NUMBER STATE FOR FINAL SYMMETRY HAS CHANGED TO NE =', ne
      ENDIF
      


      ALLOCATE(enf(ne))
      ALLOCATE(ndef(ne))

      CALL hfile(nwf2e,"dat","wf2e","bin",linitial)              !initial states


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


      CALL dmxfile(nhdmx2e, "dat", "dmx2ebf","bin", gauge, loi, lof) ! initial-symmetry h-file

      WRITE(NOUT,*) '# dmx2ebf::    gauge =  ', mode      !(1) vel


      WRITE(UNIT=NHDMX2E) MODE
      WRITE(UNIT=NHDMX2E) NMININ, NMAXI
      WRITE(UNIT=NHDMX2E) NMINFIN, NMAXF

      dipole_2e_initial_states:DO ic = nminin, nmaxi

         READ(NWF2E)  eni(ic)
         READ(NWF2E)  (ciin(j), j = 1, ntoti)

         WRITE(*,*)  '# dmx2ebf                           e_i = ', ic, eni(ic)
         WRITE(*,*)  '# dmx2ebf                           nf = ', ne


         WRITE( UNIT=NHDMX2E) eni(ic)

         WRITE( NOUT, '(e14.6,2x,i5)') eni(ic)


         CALL hfile(nh2e,"dat","h2e","bin", lof)     !final-symmetry  h-file


            READ(NH2E) ne
            
            dipole_2e_final_states: DO ir = 1, nmaxf

               READ(NH2E) enf(ir), ndef(ir)
               READ(NH2E) ((cifin(n,j), n = 1, ntotf), j = 1, ndef(ir))

               ALLOCATE( dpz( ndef(ir) ) )

               alpha = 1.0D+00 
               beta  = 0.0D+00
               incx  = 1
               incy  = 1

               WRITE( NOUT, '(e14.6,2x,i5)') enf(ir), ndef(ir)               
               WRITE( UNIT=NHDMX2E)          enf(ir), ndef(ir)


               dpz = 0.0D+00
               channels_of_final_state: DO n = 1, ndef(ir)

                  !  calculate the final dipole

                  CALL dgemv('n', ntotf, ntoti, alpha, dmxz, nhmx, ciin, incx, beta, tt, incy)

                  dpz(n) = 0.0D+00
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
      
    END PROGRAM dmx2ebf

!#######################################################################

    SUBROUTINE trnsmx(hmx,nhmx,xm,nbsp,ihr,ihc,ir,ic)
      !
      USE PRECISION, ONLY: dpk
      !
      IMPLICIT NONE
      !
      REAL(dpk), DIMENSION(nhmx,nhmx) ::  hmx     ! T_FI =  < F | T | I >
      INTEGER                         ::  nhmx    !  dim of T_FI
      REAL(dpk), DIMENSION(nbsp,nbsp) ::  xm      ! t_fi = <a|t|b> 
      INTEGER                         ::  nbsp    !  dim of t_fi 
      INTEGER                         ::  ihr     !
      INTEGER                         ::  ihc     !
      INTEGER                         ::  ir      ! channel index of final state
      INTEGER                         ::  ic      ! channel index of initial state
      !
      INTEGER                         :: irf, icf, kr, k, kc, kk



!!!..............................
      
      !  allocate(xm(nsx,nsx))
      !  allocate(hmx(nhmx,nhmx))

      
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
            
!!!SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,BETA, Y, INCY )
!!!      DOUBLE PRECISION   ALPHA, BETA
!!!      INTEGER            INCX, INCY, LDA, M, N
!!!      CHARACTER*1        TRANS
!!!      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!!!*  Purpose
!!!*
!!!*  DGEMV  performs one of the matrix-vector operations
!!!*
!!!*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!!!*
!!!*  where alpha and beta are scalars, x and y are vectors and A is an
!!!*  m by n matrix.
!!!* 
         
!#######################################################################
! EOF
