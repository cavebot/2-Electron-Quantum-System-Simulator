!----- program modified to save memory  94.5.2 ------------------------------
!----- rewritten by j. zhang to accomodate mixed configuration
!-------------------------------------------------------------------
!!!        
!!!         laan/iesl/2001
!!!
!!!         06.09.2001: 
!!!                    Reconstructing the code. Mainly I/O facilities 
!!!                    have been added for a more automatic running
!!!                    of the whole package   
!!!
!!!           note  (0) old name was   dpbf.f90 (dpbf_mod.f90,dpbf_sub77.f)  
!!!
!!!
!!!           note  (1) The file opened by the Routine HWF2EFILE   
!!!                        cannot be read at once. The file closes  
!!!                        during the execution of the programm and
!!!                        and message error comes that do not exist
!!!                        The correspondent programm dmx2e.f for the  
!!!                        fxd code (compiled with pgf77) not such   
!!!                        problem appears.
!!!                        Temporarily, i close and open again this
!!!                        file.
n!!!          
!!!           note  (2)  NH2E has to be 1 as has been implemented now
!!!                           for single ionization files.
!!!                           and  2 for double ionization files.   
!!!  
!!!
!!!  Calculates un-normalized 2-e DME from a bound state | b > to 
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
!!!   nsx  : number of orbitals <=(?)    nsx  < = number of B-splines
!!!............................


PROGRAM oscmain
!  use jl
  USE wf_files
  USE dz_value
  USE configuration
  IMPLICIT NONE
  INTEGER, PARAMETER :: dpk = KIND(1.d0)
  INTEGER j1,j2, ls
  COMMON/ang/j1,j2,ls
  CHARACTER(len=16)  config, cdiag, outfile, gauge
  CHARACTER(len=16), ALLOCATABLE, DIMENSION(:) :: inidat
  INTEGER lof, lsf, ncsmxf, loi,lsi, ncsmxi, nconfig
  INTEGER nhmx, nsx, ncs, ncsi, ncsf, nsi, nsf, ntotf, ntoti, ne
  INTEGER kf, ki, ir, ic, nr, nc, i, j, it, n, ltotalMax, iout, nout
  INTEGER LINITIAL, NMININ, NMAXI, LFINAL, NMINFIN, NMAXF, k
  INTEGER MODE, NOUT, NWF2E, NH2E, NCFG, NHDMX2E
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nd
  INTEGER, DIMENSION(4) :: id
  REAL(dpk), DIMENSION(4) :: ang
  REAL(dpk) fli, flf, phase
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: dmxz, cifin
  REAL(dpk), ALLOCATABLE, DIMENSION(:) :: ciin, eni, enf, dpz, osz, tt
  EXTERNAL dgemv
  
! ------------------------------------------------------------------
!   nhmx : nhmx > total number of configurations ( how much ?)
!   ncs  : ncs >  number of conf. series , ( how much? )
!   nsx  : number of orbitals <=(?)    nsx  < = number of B-splines
!   nwf  : may be is obsolete
! ------------------------------------------------------------------
!!!...............................................

!!! NH2E = 1, SINGLE IONIZATION CONTINUUM
!!! NH2E = 2, DOUBLE IONIZATION CONTINUUM

      NWF2E   = 3         
      NHDMX2E = 4
      NOUT    = 16
      NCFG    = 18


!!!..........................

      OPEN(UNIT=9,file='INP/HDMX2E.INP',status='old')

      READ(UNIT=9,*) LINITIAL, NMININ,  NMAXI
      READ(UNIT=9,*) LFINAL,   NMINFIN, NMAXF
      READ(UNIT=9,*) GAUGE
      READ(UNIT=9,*) NH2E
      READ(UNIT=9,*) ltotalMax
      READ(UNIT=9,*) nhmx, ncsi, ncsf, nsx

      CLOSE(9)

!!!...........................


      OPEN(UNIT=NOUT,FILE='OUT/HDMX2E.LOG')


      IF(gauge.EQ.'V') THEN

         WRITE(*,*) ' 2e-DME IN VELOCITY GAUGE '

         MODE = 0     
      ELSE
            
         WRITE(*,*) ' 2e-DME IN LENGTH GAUGE '

         MODE = 1
      ENDIF

!!!..........................


  OPEN(9,file='inp/mchdmx2e.inp',status='old')

!  READ(9,*) ltotalMax
!  READ(9,*) nhmx, ncs, nsx 
!  READ(9,*) config
!!! Initial Bound States, ( E, C(E) )
!  READ(9,577) ( inidat(i),  i = 1,  2 * (ltotalMax+1)) 
!!! 1e - dmx  < a | D | b >  
!  READ(9,577) ( aadp(i),    i = 1,  ltotalMax )        
!!! 1e - dmx  < a | D | B_i(r) >
!  READ(9,577) ( badp(i),    i = 1,  2 * ltotalMax  )   
!!! overlap   < B_i || B_j >
!  READ(9,577) ( ba_over(i), i = 1,  ltotalMax + 1)     

  READ(9,*) nout

!  allocate(inidat(2*(ltotalMax + 1)))    
!  allocate(aadp(ltotalMax))         
!  allocate(badp(2*ltotalMax))       
!  allocate(ba_over(ltotalMax+1))      


  allocate(cdiag(nout), outfile(nout)) 

      do iout = 1, nout
  
         read(9,577)  cdiag(iout), outfile(iout)

      end do


 577  format(5a16)
! ------------------------------------------------------------------
    1 format(/2x,'the initial state of the transition --')
    2 format(/2x,'the final state of the transition --')
    3 format(2x,'total l and total s=',2i3,4x,'# of config. is ',i4)
    4 format(/2x,'energy eigenvalues of initial & final states')
    5 format(2x,1p8e15.7)
    6 format(/2x,'os values: length - top & velocity - bottom.   {lf is&
     & included. for absorp. lf=2*lof+1. for emiss.("-" sign) lf=2*loi+1&
     & }' /2x,'row - final states & column - initial states')
    7 format(2x)
    8 format(/2x,'lhff,llf,lhfi,lli,id(i),ang(i) --'/)
    9 format(2x,4i2,3x,4i2,1p4e12.4)
   11 format(/2x,'excitation energy in ryd.')
! -----------------------------------------------------------------
  call config_space(ncs)


  DO iout = 1, nout      

     WRITE(*,1)

     CALL openfl(1,loi,lsi,ncsmxi,nhfi,lhfi,lli,nmini,nmxi,ndi,&
          & ntoti,nsi,inidat)

     fli = dfloat(2*loi+1)
!  write(*,2)

     OPEN(16,file=outfile(iout), form='unformatted', access='sequential')
!!! ------------------------------------------------------------------
     CALL cxfin(config, nhff,lhff, llf, nminf, nmxf, nolf,&
          & isf, lof, lsf, ncsmxf,ntotf, ndf,idcs )
     nconfig=0
     DO i=1,ncsmxf
        DO j = 1, ndf(i)   
           nconfig=nconfig+1
        END DO
     END DO

  write(*,*) 'nconfig=', nconfig

  flf = dfloat( 2*lof + 1)
!  call openfl1(17,loi,lof,lsf)

  j1 = lof
  j2 = loi

  
  print*,'nhmx=', nhmx

  allocate(dmxz(nhmx, nhmx))
  dmxz = 0.0d0

  if(lsi .ne. lsf) stop
  ls = lsi

  if(ls .eq. 0) phase =   1.0d0
  if(ls .eq. 1) phase = - 1.0d0

  print*,'phase=', phase
!-------------------------------------------------------------------
  kf=1
  allocate(dz(nsx, nsx))   
  call read_dp(nsx,ltotalMax) 
  print*,'finished reading dp'

  do 100 ir=1,ncsmxf
    ki=1
    do 101 ic=1,ncsmxi

      call angfc(lhff(ir),llf(ir),lhfi(ic),lli(ic),lof,loi,id,ang)
!-------------------------------------------------------------------
      select case(idcs(ir))
      case(0)

         dz = 0.0d0
         call cal_nlnl(ir, ic, id, ang, phase)

      case(1)

         dz = 0.0d0
         call cal_nlsp(ir, ic, id, ang, phase)

      end select

!!!      write(*,*) 'dz()=', (dz(i,1),i=1,ndf(ir))
!!!      write(*,*) 'dz()=', (dz(i,2),i=1,ndf(ir))
!!!      select case(ir.ne.1 .or. ic.ne.1) 
!!!        case(.false.)
!!!        call mxprnt(dz,ndf(ir),ndi(ic))
!!!      end select

      call trnsmx(dmxz,nhmx,dz,nsx,kf,ki,ndf(ir),ndi(ic))

      ki=ki+ndi(ic)

  101 continue

      kf=kf+ndf(ir)

  100 continue
!-------------------------------------------------------------------

  deallocate(dz)   
  deallocate(vz)
  deallocate(vzba)
  deallocate(sp_a)
                                ! Taking care of sqrt(2)
      it=1
      select case( mod(lhfi(1)+lli(1),2) ) ! initial state parity
      case(0)
         do ic=1, ncsmxi
            select case(lhfi(ic).eq.lli(ic).and.nhfi(ic).eq.nmini(ic))
            case(.true.)
               do nr=1, nhmx
                  dmxz(nr,it)=dmxz(nr,it)/dsqrt(2.0d0)
               end do
            end select
            it = it + ndi(ic)
         end do
      case(1)
         do ir=1,ncsmxf
            select case(idcs(ir))
            case(0)
               select case((lhff(ir).eq.llf(ir)).and.(nhff(ir).eq.nminf(ir)))
               case(.true.)
                  do nc=1, nhmx
                     dmxz(it,nc) =dmxz(it,nc)/dsqrt(2.0d0)
                  end do
!                  write(*,*) 'it=', it
               end select
            end select 
            it = it + ndf(ir)
         end do
      end select

      print*,'sqrt2 finished'

  allocate(eni(ntoti))
  allocate(ciin(ntoti))
  allocate(cifin(nconfig,ncsmxf))
  allocate(tt(nconfig))

      open(2, file = cdiag(iout), form='unformatted', access='sequential')
      read(2) ne
      close(2)

  write(16) nsi,ne

  do ic = 1, nsi
     
!!!  read final (continuum) states config. coef 

      open(2, file = cdiag(iout), form='unformatted', access='sequential')
      read(2) ne
      allocate(enf(ne))
      allocate(nd(ne))
                                !  read initial states config. coef
      print*,' begin to read'
      read(1) eni(ic)
      read(1) (ciin(j), j=1, ntoti)

!      write(*,*)  'eni=',eni(ic)


      write(16)  eni(ic)

      do ir = 1, ne
         read(2) enf(ir), nd(ir)
         read(2) ((cifin(n,j), n=1,nconfig), j=1,nd(ir))
!    print*,' read cifin'
!
         allocate(dpz(nd(ir)))
         dpz=0.0d0
         write(16) enf(ir), nd(ir)
         do n = 1, nd(ir)

!!!  calculate the final dipole

            call dgemv('n', nconfig, ntoti, 1.0d0, dmxz, nhmx, ciin, 1, &
                  & 0.0d0, tt, 1)

!!!        dpz(n) = dot_product(cifin(:,n), tt)

                  dpz(n)=0.0d0

                  do j = 1, nconfig
                     dpz(n) = dpz(n) + cifin(j,n) * tt(j)
                  end do

               end do

!      write(*,*)  enf(ir)-eni(ic), (dpz(j), j=1,nd(ir))

               write(16) (dpz(j), j=1,nd(ir))

               deallocate(dpz)

!      deallocate(osz)
            end do
            close(2)

            deallocate(enf)
            deallocate(nd)

         end do 

         deallocate(eni)
         deallocate(ciin)
         deallocate(cifin)
         deallocate(tt)
         close(1)
     deallocate(dmxz)
     close(16)
  end do
end program 
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
subroutine trnsmx(hmx,nhmx,xm,nsx,ihr,ihc,ir,ic)
  implicit none
  integer nhmx, nsx
  integer ihr, ihc, ir,ic, irf, icf, kr, k, kc, kk
  real(8), dimension(nsx,nsx) ::  xm
  real(8), dimension(nhmx,nhmx) ::  hmx
  
!  allocate(xm(nsx,nsx))
!  allocate(hmx(nhmx,nhmx))

  irf=ihr+ir-1
  icf=ihc+ic-1
  do 10 k=ihr,irf
    do 10 kk=ihc,icf
      kr=k-ihr+1
      kc=kk-ihc+1
      hmx(k,kk)=xm(kr,kc)
10   continue
end subroutine trnsmx
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!EOF

