! Purpose: 
!         Calculate the 1-photon single + double ejection cross 
!         section for two-electron atoms
!         
!        cs_1( E ) 
!  
!   author:
!   l.aa.n/iesl/Jan 2002
!################################################################
MODULE configuration
  IMPLICIT NONE

  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhi, lhi, li, nmini, nmxi, ndi,  isi,  noli, idcsi
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhf, lhf, lf, nminf, nmxf, ndf,  isf,  nolf, idcsf


CONTAINS
  SUBROUTINE config_space(ncsi, ncsf)  
    IMPLICIT NONE
    INTEGER ncsi, ncsf

    ALLOCATE(nhi(ncsi))       ! initial states
    ALLOCATE(lhi(ncsi))
    ALLOCATE(li(ncsi))
    ALLOCATE(nmini(ncsi))
    ALLOCATE(nmxi(ncsi))
    ALLOCATE(ndi(ncsi))
    ALLOCATE(noli(ncsi))
    ALLOCATE(isi(ncsi))
    ALLOCATE(idcsi(ncsi))
    ALLOCATE(nhf(ncsf))      ! final states
    ALLOCATE(lhf(ncsf))
    ALLOCATE(lf(ncsf))
    ALLOCATE(nminf(ncsf))
    ALLOCATE(nmxf(ncsf))
    ALLOCATE(ndf(ncsf))
    ALLOCATE(nolf(ncsf))
    ALLOCATE(isf(ncsf))
    ALLOCATE(idcsf(ncsf))

  END SUBROUTINE config_space
END MODULE configuration

!!!###################################################
PROGRAM cs1ph_f

  USE PRECISION, ONLY: dpk
  USE units, ONLY:enau,alpha
  USE configuration
  USE ioroutines

  IMPLICIT NONE

  REAL(DPK) EI, EF
  INTEGER i, j, m, ndmax
  INTEGER nSymmetries, ncsm
  INTEGER LINITIAL, LFINAL, LOI,LOF
  INTEGER nei, ii, lsi, ncsmxi, ntoti, ncsi, nconfigi
  INTEGER nef, ir, lsf, ncsmxf, ntotf, ncsf
  INTEGER MODE, K
  INTEGER NOUTFILE, NBF, NCS1PH, NDCS1PH, NNDMX2E
  INTEGER N1E, NCFG
  INTEGER IC, IF, ID
  REAL(DPK) CS1PH_O
  !
  INTEGER llmax, llmin, lp, ie, nvalence, nmx                 ! 1e info
  INTEGER nl_core, ne_core 
  REAL(DPK),    ALLOCATABLE, DIMENSION(:,:)   :: ehf
  REAL(DPK)                                   :: e_core 
  REAL(DPK),    ALLOCATABLE, DIMENSION(:)     :: eni, enf     ! 2E energies, channels
  INTEGER,      ALLOCATABLE, DIMENSION(:)     :: nof
  INTEGER,      ALLOCATABLE, DIMENSION(:,:)   :: n1,l1,l2
  INTEGER,      ALLOCATABLE, DIMENSION(:,:)   :: nch_cs       ! cross sections
  REAL(DPK),    ALLOCATABLE, DIMENSION(:,:)   :: cs1ph_single
  REAL(DPK),    ALLOCATABLE, DIMENSION(:,:)   :: cs1ph_double, ratio_ds
  REAL(DPK),    ALLOCATABLE, DIMENSION(:,:)   :: Wph
  REAL(DPK),    ALLOCATABLE, DIMENSION(:,:,:) :: cs1ph_ch
  REAL(DPK),    ALLOCATABLE, DIMENSION(:,:,:) :: cs1ph_double_l
  COMPLEX(DPK), ALLOCATABLE, DIMENSION(:,:,:) :: dmxbf
  REAL(DPK)                                   :: threshold_2

  !
  CHARACTER( LEN = 6   )   sgauge
  CHARACTER( LEN = 100 )   argv 

!!!...................................


  CALL GETARG(1, ARGV)
  READ(ARGV,*) linitial
  CALL GETARG(2, ARGV)
  READ(ARGV,*) lfinal
  CALL GETARG(3, ARGV)
  READ(ARGV,*) sgauge


!...



  NBF        = 3        !  < g          |  D  |  E_I, I (C)  >
  NCFG       = 8        ! CONFIGURATION FILES L = 0, 1, 2, ...
  N1E        = 9        ! OUTPOUT FILE  1
  NCS1PH     = 10       ! OUTPOUT FILE  2
  NDCS1PH    = 19       ! OUTPOUT FILE  3
  NOUTFILE   = 16       ! LOG     FILE 


  threshold_2 = 0.0D+00 

  MODE = 1
  IF(sgauge.EQ.'V')          MODE = 0     


!...
  

  OPEN(UNIT=NOUTFILE,FILE='out/cs1ph.out')

  WRITE(*,*) '# cs1ph:: calculation gauge = ', sgauge




  CALL cfile(ncfg,"inp","cfg",linitial)      ! read initial cfg file
  READ(NCFG, 3) nSymmetries
  READ(NCFG, 3) loi, lsi
  READ(NCFG, 3) ncsi
  CLOSE(NCFG) 


  CALL cfile(ncfg,"inp","CFG",lfinal)      !   read final cfg file
  READ(ncfg, 3) nSymmetries
  READ(ncfg, 3) lof, lsf
  READ(ncfg, 3) ncsf
  CLOSE(NCFG) 
  

  !check
  IF(LOI.NE.LINITIAL) THEN
     WRITE(*,*) "# csh1ph::  LOI NOT EQUAL TO LINITIAL "
     WRITE(*,*) "# csh1ph::                      LOI = ",LOI
     WRITE(*,*) "# csh1ph::                 LINITIAL = ", LINITIAL
     STOP
  ENDIF
  !
  IF(LOF.NE.LFINAL) THEN
     WRITE(*,*) "#csh1e::    LOF NOT EQUAL TO LFINAL "
     WRITE(*,*) "#csh1e::                      LOF = ",LOF
     WRITE(*,*) "#csh1e::                   LFINAL = ", LFINAL
     STOP
  ENDIF

  !

  CALL config_space(ncsi, ncsf)      ! set dim for initial and final symmetries

  CALL cxfin(ncfg, nhf,lhf,lf,nminf,nmxf,nolf,isf,lof,lsf,ncsmxf,ntotf,ndf,idcsf)


  !..  total number of configurations for final states


!           Read 1-e core energies

      WRITE(*,*) '# cs1ph::  reading 1-e energies. use of en1e.dat.'
 
      OPEN(N1E, FILE='dat/en1e.dat')
      READ(N1E, *) LLMIN, LLMAX


      nl_core = ABS(llmax - llmin) + 1   !note: can be read from cfg-L.inp file
      ne_core = MAXVAL(nmxf) - 2 

      WRITE(*,*) '# cs1ph::  nof partial waves,         nl_core = ', nl_core
      WRITE(*,*) '# cs1ph::  max nof states in p. waves ne_core = ', ne_core

      ALLOCATE( ehf(nl_core, ne_core) )

      DO  LP = LLMIN, LLMAX
         READ(N1E, 4) NMX, NVALENCE
         DO  IE = 1, NMX
            READ(N1E,'(E25.14)') EHF(LP,IE)
         ENDDO
      ENDDO

      !check
      IF(nmx.NE.ne_core) THEN
         WRITE(*,*) ' ERROR IN NUMBER OF CORE STATES AND B-SPLINES'
         WRITE(*,*) ' dat/en1e.dat :                 nmx = ', nmx
         WRITE(*,*) ' INP/CFG.INP  :             ne_core = ', ne_core
         STOP
         ENDIF

      CLOSE(N1E)

      WRITE(*,*) '# cs1ph::                            lfn = ', lfinal 
      WRITE(*,*) '# cs1ph:                              Sf = ',  lsf
      WRITE(*,*) '# cs1ph::                           ncfg = ', ntotf
      WRITE(*,*) '# cs1ph::    max nof channels      ndmax = ', ncsf
      WRITE(*,*) '# cs1ph::       opening dmx2ebf file. ' 
      WRITE(*,*) '# cs1ph::                            Lin = ', linitial 
      WRITE(*,*) '# cs1ph:                             Sin = ',  lsi


      CALL dmxfile(nndmx2e, "dat","ndmx2ebf","bin",sgauge, linitial, lfinal)


!!! NMAXI : number of initial states
!!! NMAXF : number of final   states


      READ( UNIT = NNDMX2E) MODE
      READ( UNIT = NNDMX2E) nei, nef


      WRITE(*,*) '# cs1ph:          nof initial states, nei = ', nei
      WRITE(*,*) '# cs1ph:          nof finall  states, nef = ', nef

      ALLOCATE( eni(nei) )
      ALLOCATE( enf(nef) )
      ALLOCATE( nof(nef) )
      ALLOCATE( dmxbf( nei, nef, ncsf) )
      ALLOCATE( n1(nef, ncsf) )
      ALLOCATE( l1(nef, ncsf) )
      ALLOCATE( l2(nef, ncsf) )

      DO ic = 1,  nei

         READ(UNIT=NNDMX2E)   ii, eni(ic)

         WRITE(*,*) '# cs1ph::     initial state energy,  ein =', eni(ic), ic 
            WRITE(*,'(a30,1X,a12,I3,a5,1X,E15.7)') &
                 &'# cs1ph:: total initial energy', ' E(',ic,') = ', eni(ic)

         DO ir = 1,  nef

            READ(NNDMX2E)  if, enf(ir), nof(ir)

            WRITE(*,'(a28,1X,a14,I3,a5,1X,E15.7,a10,1X,i4)') &
                 &'# cs1ph:: total final energy', ' E(',ir,') = ', enf(ir),' nd = ', nof(ir) 

            READ(NNDMX2E) ( n1(ir, j),            j = 1, nof(ir) )
            READ(NNDMX2E) ( l1(ir, j),            j = 1, nof(ir) )
            READ(NNDMX2E) ( l2(ir, j),            j = 1, nof(ir) )
            READ(NNDMX2E) ( dmxbf(ic, ir, j),     j = 1, nof(ir) )


         ENDDO

      ENDDO

!
!        1           2                    j             nof(ir+1)     
! ir+1 ------     -------     .....    --------  ....   --------   |enf(ir+1)
!
!  
!        1           2                     j            nof(ir)
! ir   ------     -------      ....     -------- ....   --------   |enf(ir)
!        ^                                ^
!        |                                !
!        | dmx(ic,ir,1)                   ! dmx(ic,ir,j)
!        |                                !
!        |                                !
!ic  ----|------------------------------------------------------   | eni(ic)
!
!



! notes:
!
!  eni(ic) : Energy (Ryd) of the 2e initial state  ic = 1,2,..,nei(ic)
!  nei     : Total number of initial states

!  eni(ir) : Energy (Ryd) of the 2e final state  ir = 1, 2,.., nef(ir)
!  nef     : Total number of final states

!  ncsf    : maximum number of channels in the calculation
!  nof(ir) : number of open channels at energy  E == enf(ir) 
!           nof(ir) < ncsf      for all ir


      CALL dmxfile(ncs1ph,"dat", "cs1ph", "ascii", sgauge, linitial, lfinal)
      CALL dmxfile(ndcs1ph,"dat","dcs1ph","ascii", sgauge, linitial, lfinal)

      ALLOCATE( cs1ph_double_l( nei, nef, nl_core) )
      ALLOCATE(       cs1ph_ch( nei, nef,    ncsf) )
      ALLOCATE(   cs1ph_single( nei, nef) )
      ALLOCATE(   cs1ph_double( nei, nef) )
      ALLOCATE(         nch_cs( nei, nef) )
      ALLOCATE(       ratio_ds( nei, nef) )
      ALLOCATE(            Wph(nei, nef) )

      nch_cs = 0

      DO ic = 1,  nei

         cs1ph_single = 0.0D+00 ;
         cs1ph_double = 0.0D+00 ;
         cs1ph_double_l = 0.0D+00 ;

         DO ir = 1,  nef

            Wph(ic, ir) = ABS( enf(ir) - eni(ic) )/ 2  ! enf, eni in Rydbers

            WRITE(*,'(a24,1X,a14,I3,a5,1X,E15.7,a10,1X,i4)') &
                 &'# cs1ph:: photon energy', ' W(',ir,') = ', Wph(ir,ic),' nd = ', nof(ir) 

            cs1ph_single(ic, ir) = 0.0D+00 ;
            cs1ph_double(ic, ir) = 0.0D+00 ;
            
            cs1ph_o = 5.337819D+00 * 1.5D+00 * Wph(ic, ir)
            IF(mode.EQ.0) cs1ph_o = 5.337819D+00 * 1.5D+00 / Wph(ic, ir)


            cs1ph_ch(ic,ir,:) = cs1ph_o * ABS( dmxbf(ic,ir,:) )**2 


!            Wph(ic, ir) = ABS( enf(ir) - eni(ic) )/ 2  ! enf, eni in Rydbers

!            WRITE(NDCS1PH,'(32e15.7)') enau * wph(ic,ir), cs1ph_ch(ic,ir,:)
!            WRITE(NDCS1PH,'(32e15.7)') enf(ir), cs1ph_ch(ic,ir,:)

            e_core = ehf( l1(ir,1) + 1, n1(ir,1))

            WRITE(NDCS1PH,'(32e15.7)') enf(ir)-e_core, cs1ph_ch(ic,ir,:)



            DO id = 1, nof(ir)

!
!      separate in SI/DI according the target core energy e1=e_core
!
!           E = e_core + e    --->   e = ( E - e_core )
!
!           E       : total energy  
!           e_core  : energy of the inner electron '1'   e_1
!           e       : energy of the outer electron '2'   e_2
!
!
!
                  
        !  target energy for the channel id :    | nhf, lhf, lf >
                  
!                  e_core = ehf( lhf(id) + 1, nhf(id)+ lhf(id) )  

                  e_core = ehf( l1(ir,id) + 1, n1(ir,id))

                  si_or_di_ionization:IF(e_core < threshold_2) THEN

                     cs1ph_single(ic, ir) = cs1ph_single(ic, ir ) + cs1ph_ch(ic,ir,id)

                     WRITE(*,'(4I4,1X,2E15.7,1X,a20)') &
                          &id, n1(ir,id), l1(ir,id), l2(ir,id), e_core, cs1ph_single(ic, ir),' SI'
                  ELSE 

                     cs1ph_double(ic, ir) = cs1ph_double(ic, ir ) + cs1ph_ch(ic,ir,id) 

                     cs1ph_double_l(ic,ir,l1(ir,id)+1 ) = cs1ph_double_l(ic,ir,l1(ir,id)+1) + cs1ph_ch(ic,ir,id)
                     
                     WRITE(*,'(4I4,1X,3E15.7,1X,a20)')  id, n1(ir,id), l1(ir,id), l2(ir,id), e_core& 
                          &,cs1ph_double_l(ic, ir,l1(ir,id)+1), cs1ph_double(ic, ir),'DI'

                  ENDIF si_or_di_ionization
                     
               ENDDO
               
           
!            IF(IC.EQ.1) THEN
!               IF(cs1ph_single(ic,ir).NE.0.D+00 ) THEN
!                  ratio_ds(ic,ir) = cs1ph_double(ic,ir) / cs1ph_single(ic,ir)
!               ELSE
!                  ratio_ds(ic,ir) = 1.0D+00
!               ENDIF

               WRITE(NCS1PH,'(8e15.7)') EnAU * Wph(ic,ir), &
                                    & cs1ph_single(ic,ir), &
                                    & cs1ph_double(ic,ir),  &
                                    & cs1ph_double_l(ic,ir,:)

!, ratio_ds(ic,ir) 


!            ENDIF


            ENDDO
         ENDDO


      CLOSE(NCS1PH)
      CLOSE(NDCS1PH)



! clean up


   DEALLOCATE(ehf)
   DEALLOCATE(enf)      !  final states
   DEALLOCATE(nhf,lhf,lf)
   DEALLOCATE(n1,l1,l2)
   DEALLOCATE(nminf)
   DEALLOCATE(nmxf)
   DEALLOCATE(ndf)
   DEALLOCATE(nolf)
   DEALLOCATE(isf)
   DEALLOCATE(idcsf)
   DEALLOCATE(eni)      !initial states
   DEALLOCATE(nhi)
   DEALLOCATE(lhi)
   DEALLOCATE(li)
   DEALLOCATE(nmini)
   DEALLOCATE(nmxi)
   DEALLOCATE(ndi)
   DEALLOCATE(noli)
   DEALLOCATE(isi)
   DEALLOCATE(idcsi)
   DEALLOCATE(dmxbf)    ! dme, cross section, ...
   DEALLOCATE(cs1ph_ch)
   DEALLOCATE(cs1ph_single)
   DEALLOCATE(cs1ph_double)
   DEALLOCATE(cs1ph_double_l)
   DEALLOCATE(ratio_ds)
   DEALLOCATE(Wph)

!!!!!!!!!!!

3     FORMAT(6I5)
4     FORMAT(6i5)

!!!!!!!!!!


 END PROGRAM cs1ph_f
!!!#####################################################
!!!EOF

