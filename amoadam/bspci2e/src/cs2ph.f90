!!! Purpose: Normalizes the 2-e unormalized dme obtained 
!!!          from hdmx2eff.f90 program, by multipliing with 
!!!          data obtained from the kmtx.f90 program.
!!!         
!!!        < a | d | b >  ----> Norm (*) < a| d| b > 
!!!  
!!!   l.aa.n/iesl/2001
!!! 
!!!                     
!!!                  *   no input file as 'dpnorm.inp'
!!!                  *   modified to compiled by pgf90/pgf77 compilers 
!!!                  
!!!
!!!################################################################
MODULE configuration
  IMPLICIT NONE

  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhf, lhf, ll, nmin, nmx, nd,  is,  nol,  idcs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhfi,lhfi,lli,nmini,nmxi,ndi, isi, noli, idcsi
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhff,lhff,llf,nminf,nmxf, ndf,isf, nolf, idcsf


CONTAINS
  SUBROUTINE config_space(ncs, ncsi, ncsf)  
    IMPLICIT NONE
    INTEGER ncs, ncsi, ncsf

!!!...............  INITIAL STATE

    ALLOCATE(nhf(ncs))
    ALLOCATE(lhf(ncs))
    ALLOCATE(ll(ncs))
    ALLOCATE(nmin(ncs))
    ALLOCATE(nmx(ncs))
    ALLOCATE(nd(ncs))
    ALLOCATE(nol(ncs))
    ALLOCATE(is(ncs))
    ALLOCATE(idcs(ncs))

!!!...............  INTERMEDIATE STATES

    ALLOCATE(nhfi(ncsi))
    ALLOCATE(lhfi(ncsi))
    ALLOCATE(lli(ncsi))
    ALLOCATE(nmini(ncsi))
    ALLOCATE(nmxi(ncsi))
    ALLOCATE(ndi(ncsi))
    ALLOCATE(noli(ncsi))
    ALLOCATE(isi(ncsi))
    ALLOCATE(idcsi(ncsi))

!!!...............  FINAL STATE

    ALLOCATE(nhff(ncsf))
    ALLOCATE(lhff(ncsf))
    ALLOCATE(llf(ncsf))
    ALLOCATE(nminf(ncsf))
    ALLOCATE(nmxf(ncsf))
    ALLOCATE(ndf(ncsf))
    ALLOCATE(nolf(ncsf))
    ALLOCATE(isf(ncsf))
    ALLOCATE(idcsf(ncsf))

  END SUBROUTINE config_space
END MODULE configuration

!!!###################################################
PROGRAM cs2ph_f

  USE units
  USE configuration
  USE ioroutines

  IMPLICIT NONE

  INTEGER, PARAMETER::NPHOTONS = 2 
  CHARACTER(LEN=16) GAUGE
!  INTEGER EnAU
  REAL(8) CS_2, CS_2_SI 

  INTEGER i, j, m, ndmax
  INTEGER nSymmetries, l, ncsm
  INTEGER ii, ne,  n,  ls,  ncsmx,  ntot,  lo,  ncs,  nconfig
  INTEGER ic, ici, nei, ni, lsi, ncsmxi, ntoti, loi, ncsi, nconfigi
  INTEGER icb, icib, neib, nib, lsib, ncsmxib, ntotib, loib, ncsib, nconfigib
  INTEGER ir, irf, nef, nf, lsf, ncsmxf, ntotf, lof, ncsf, nconfigf
  INTEGER irb, nefb 
  INTEGER nEn, modebf, nchannels, modebf1

  INTEGER MODE, K
  INTEGER NDIPOLE,  NOUTFILE, NCFG, NFF, NBF, NBF1, NBB, NCS2PH
  INTEGER LINITIAL, LINTERMEDIATE, LFINAL

  REAL(8) echeck, energy 
  REAL(8) ENERGY_INITIALB, ENERGY_INITIAL
  REAL(8) cs2ph_E,  dmx_2_channel 
  REAL(8) cs2ph_Eb, dmx_2_channelb 
  REAL(8) dmx_2, DE, WEIGHT, DE_STEP
  REAL(8) cs2phbb, dwph, eph
  REAL(8) cs_n_au, cs_n_si

  INTEGER, ALLOCATABLE, DIMENSION(:)   :: no, noi, nof
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: nob, noib, nofb
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: umx
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: en, eni, enf, W_PHOTON
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: enibb, enib, enfb


  
  REAL(8), ALLOCATABLE, DIMENSION(:)        :: dmxbb
  COMPLEX(8), ALLOCATABLE, DIMENSION(:)     :: dmxbf
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:,:) :: dmxbfb
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)   :: dmxff
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:)   :: cs2ph, cs2phb,dmx2phb
  COMPLEX(8)  uim
  CHARACTER(LEN=100) ARGV

!!!...................................

  CS_2    = 2 * M_PI * 8 * M_PI * M_PI * M_PI * ALPHA**2
  CS_2_SI = 1.894D-53

!!!...................................


  NBB      = 1        !  < g          |  D  |  E_I (B)     >
  NBF1     = 2        !  <(B) E_I     |  D  |  E_F, J      > 
  NBF      = 3        !  < g          |  D  |  E_I, I (C)  >
  NFF      = 4        !  < E_I, I,(C) |  D  |  E_F, J, (C) > 

  NCFG     = 8        ! CONFIGURATION FILES L = 0, 1, 2
  NDIPOLE  = 9        ! OUTPOUT FILE  1
  NCS2PH   = 10       ! OUTPOUT FILE  2
  NOUTFILE = 16       ! LOG     FILE 

!!!..........................
!!!..........................

  CALL GETARG(1, ARGV)
  READ(ARGV,*) LINITIAL         
  CALL GETARG(2, ARGV)
  READ(ARGV,*) LFINAL       
  CALL GETARG(3, ARGV)
  READ(ARGV,*) GAUGE
  
  
  CS_N_AU  = 2 * M_PI * ( 2 * M_PI * ALPHA )**NPHOTONS
  CS_N_SI  = CS_N_AU * T0**(NPHOTONS-1) * A0**(2*NPHOTONS) 
           
!....................................

  WRITE(*,*) '# cs2ph::'
  WRITE(*,*) '# cs2ph::                          nphotons = ', NPHOTONS
  WRITE(*,*) '# cs2ph::                          cs_n_au  = ', CS_N_AU
  WRITE(*,*) '# cs2ph::                          cs_n_si  = ', CS_N_SI


  OPEN(UNIT=NOUTFILE,FILE='OUT/CS2PH.LOG')

!!!..........................

  MODE = 1
  IF(gauge.EQ.'V') MODE = 0
      
  WRITE(*,*) '# cs2ph::                       dmx in ', gauge, ' gauge'
  WRITE(*,*) '# cs2ph::                initial symmetry   li = ', LINITIAL
  WRITE(*,*) '# cs2ph::                final   symmetry   lf = ', LFINAL
  

!!!.............................................................

      ! configuration files 

3     FORMAT(6I5)

!!!.. initial       ( bound states ) 
      

      CALL CFGFILE(NCFG, LINITIAL)
      
      READ(NCFG, 3) nSymmetries
      READ(NCFG, 3) l, ls
      READ(NCFG, 3) ncsm
      
      CLOSE(NCFG) 

      ncs = ncsm

!!!..... intermediate   

      LINTERMEDIATE = LINITIAL + 1

!!!...   ( continuum states )

      CALL HCFGFILE(NCFG, LINTERMEDIATE)
      
      READ(ncfg, 3) nSymmetries
      READ(ncfg, 3) l, ls
      READ(ncfg, 3) ncsm
      
      CLOSE(NCFG) 

      ncsi = ncsm


!!!..... final ( continuum states) 

      CALL HCFGFILE(NCFG, LFINAL)
      
      READ(ncfg, 3) nSymmetries
      READ(ncfg, 3) l, ls
      READ(ncfg, 3) ncsm
      
      CLOSE(NCFG) 

      ncsf = ncsm


!!!...................................


      CALL config_space(ncs, ncsi, ncsf)

      

      WRITE(*,*) '# csh2ph::   read cfgs for intermediate symmetry   lt = ', LINTERMEDIATE


!!! intermediate continuum states

      loi = LINTERMEDIATE

      CALL cxfin(NCFG,nhfi,lhfi,lli,nmini,nmxi,noli,isi,loi,lsi,ncsmxi,ntoti,ndi,idcsi)


      !..  total number of configurations for intermediate states


      nconfigi = 0
      DO i = 1, ncsmxi
         DO j = 1, ndi(i)   
            
            nconfigi = nconfigi + 1

         END DO
      END DO

      WRITE(*,*) '# csh2ph::                         initial symmetry   ncsi = ', nconfigi



!!!.......... final states

      WRITE(*,*) '# csh2ph::   read cfgs for final  symmetry.'




      lof = LFINAL

      CALL cxfin(NCFG,nhff,lhff,llf,nminf,nmxf,nolf,isf,lof,lsf,ncsmxf,ntotf,ndf,idcsf )

      !..  total number of configurations for final states

      nconfigf = 0
      DO i = 1, ncsmxf
         DO j = 1, ndf(i)   
               
            nconfigf = nconfigf + 1
            
         END DO
      END DO

      WRITE(*,*) '# csh2ph::                        final symmetry   ncsf = ', nconfigf

!!!.....................................  checks


      WRITE(*,*) '# cs2ph::                            SI = ',  LSI
      WRITE(*,*) '# cs2ph::                            SF = ',  LSF
 
      IF(lsf.NE.lsi) THEN
   
         WRITE(*,*) '# cs2ph::  dipole forbidden transition'
         STOP

      ELSE 

         WRITE(*,*) '# cs2ph::              singlet symmetry'

         IF(lsf.EQ.3) WRITE(*,*) '# cs2ph::              triplet symmetry'

      ENDIF


!!!..........................

      
      IF(ncsf.GT.ncsi) THEN

         ndmax = ncsf
         
      ELSE IF(ncsi.GT.ncsf) THEN

         ndmax = ncsi

      ELSE

         ndmax = ncsi
         
      ENDIF

      WRITE(*,*) '# cs2ph::                 max number of channels  ndmax = ', ndmax

!!!...................................





      ALLOCATE(umx(ndmax,ndmax))

!      ALLOCATE(dmxbfb(ndmax))
      ALLOCATE(dmxbf(ndmax))
      ALLOCATE(dmxff(ndmax,ndmax))


!!!......................................

      CALL  BFDMX2EFILE(NBF, LINITIAL, LINTERMEDIATE, MODE)  ! 0 --> 1
      CALL  BBDMX2EFILE(NBB, LINITIAL, MODE)                 ! 0 --> 1
      CALL  FFDMX2EFILE(NFF,  LINTERMEDIATE, LFINAL, MODE)   ! 1 -- > 0 (2) 
      CALL  BFDMX2EFILE(NBF1, LINTERMEDIATE, LFINAL, MODE)   ! 1 -- > 0 (2) 

      WRITE(*,*) '# cs2ph::   read for  first photon done.', linitial,lintermediate
      WRITE(*,*) '# cs2ph::   read for second photon done.', lintermediate, lfinal

      OPEN(NDIPOLE, file="DAT/CS2PH.DAT")

! 1st photon 


      ! bound-bound part of cs2   ( 0 --> 1 )

!..................................


      READ(NBB)  MODE 
      READ(NBB)  lo, loi, n, ni

      ALLOCATE(  enibb( ni )  ) 
      ALLOCATE(  dmxbb( ni )  )

      READ(NBB) ENERGY_INITIALB
      READ(NBB) (enibb(ni), icib = 1, ni)

!      DO icb = 1,  ni    !   | i > = | icb, eni(icb)> 
         READ(NBB) (dmxbb(icb), icb = 1, ni)
!      ENDDO

      WRITE(*,*) dmxbb
!..................................

      modebf1 = mode


      WRITE(*,*) '# cs2ph::  bb :               0 -> 1      lo, n = ', lo, n 
      WRITE(*,*) '# cs2ph::  bb :               0 -> 1    loi, ni = ', loi, ni 
      WRITE(*,*) '# cs2ph::   bound-bound  data              eni  = ', ENERGY_INITIALB


      ! bound-free part of cs2 
      

!!!       | G, L, 1 > ---------> | IC, L + 1, ( ic = 1 - nei ) >


         READ(NBF) MODE
         READ(NBF) ne, nei
         READ(NBF) ii, energy_initial

         nEn    = nei
         modebf = MODE


         WRITE(*,*) '# cs2ph::  bf             0 -> 1  :  ne, nei = ', ne, nei
         WRITE(*,*) '# cs2ph:: bound-free data                eni = ', energy_initial


!!!..............

 
!!!     READ FOR SECOND PHOTON        FF
!!!     | IC, L + 1, ( ic = 1 - nei ) > -----> | IR, L +- 1, ( ir = 1 - nef) >


!!..........

         ! BOUND - FREE  PART  

         READ(NBF1) MODE
         READ(NBF1) neib, nefb

         ALLOCATE( enib( neib ))  ! b data
         ALLOCATE( noib( neib ))
         ALLOCATE( enfb( nefb ) )         
         ALLOCATE( nofb( nefb ) )
         ALLOCATE( dmxbfb( nefb, ndmax, neib ) )  


         DO icb = 1,  neib    !   | i > = | icb, eni(icb)> 

            READ(NBF1) icib,  enib(icb)  

            DO irb = 1,  nefb       !   | f > = | irb, enf(irb), nof(irb) >

               READ(NBF1) irf,  enfb(irb), nofb(irb)
               READ(NBF1) (dmxbfb(irb,j,icb), j = 1, nofb(irb))
            END DO
         END DO


         WRITE(*,*) '# cs2ph::  bf         1 -> 0  :  neib, nefb = ', ne, nei

         
         IF(mode.NE.modebf1) THEN 

            WRITE(*,*) '# cs2ph::  bb and bf dmx not in the same mode. stop.'
            WRITE(*,*) '# cs2ph::                                bb  mode = ', modebf1 
            WRITE(*,*) '# cs2ph                                  bf  mode = ', mode  
            STOP
           
         ENDIF



!!!..........

         ! free-free part of cs2


         READ(NFF) MODE
         READ(NFF) nei, nef


!!.....

         WRITE(*,*) '# cs2ph::  bf transition   to ', nEn, ' intermediate states'

         IF(nei.NE.nEn) THEN 
            WRITE(*,*) '# cs2ph::  bf and ff dmx not consistent. stop.'
            WRITE(*,*) '# cs2ph::  bf transition   to ', nEn, ' intermediate states'
            WRITE(*,*) '# cs2ph::  ff transition from ', nei, ' intermediate states'
            STOP
         ENDIF

!!!...........

         WRITE(*,*) '# cs2ph::     nof final states  nefb = ', nefb 

         IF(nefb.ne.nef) THEN 
            WRITE(*,*) '# cs2ph::  bf and ff dmx not consistent. stop.'
            WRITE(*,*) '# cs2ph::  bf transition   to ', nefb, ' intermediate states'
            WRITE(*,*) '# cs2ph::  ff transition from ', nef, ' intermediate states'
            STOP
         ENDIF

!!!..............

         IF(modebf.NE.mode) THEN 
            WRITE(*,*) '# cs2ph::  bb and bf dmx not in the same mode. stop.'
            WRITE(*,*) '# cs2ph::                                bb  mode = ', modebf 
            WRITE(*,*) '# cs2ph                                  ff  mode = ', mode  
            STOP
         ENDIF

         NE = 1         
!         WRITE(*,*) MODE, NE,NEI,NEF
   


!!!........................................................................
         

         WRITE(NDIPOLE,*) MODE
         WRITE(NDIPOLE,*) NEF

!!!.........................................


         ALLOCATE(en(ne))
         ALLOCATE(no(ne))
         ! intermediate
         ALLOCATE( eni(  nei  ))    ! c data
         ALLOCATE( noi(  nei  ))
         WRITE(*,*) '# cs2ph:: allocation for intermediate states done.  nei = ', nei
         ! final

         ALLOCATE( enf(  nef  ) )          ! c data
         ALLOCATE( nof(  nef  ) )

         ALLOCATE(W_PHOTON( nef ) )
         ALLOCATE( dmx2phb( nef, ndmax ) )  
         ALLOCATE( cs2phb( nef, ndmax ) )  
         ALLOCATE( cs2ph(  nef, ndmax ) )

         WRITE(*,*) '# cs2ph:: allocation for final states done.  nef = ', nef




!!!...............................................
! set units for complex number and real matrix


         uim = (0.0D+00, 1.0D+00)
         umx = 0.0D+00

         DO i = 1, ndmax

            umx(i,i) = 1.0D+00
         END DO
!!!


!!!!          BOUND PART OF CROSS SECTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



         WRITE(*,*) '# cs2ph:: bound part of cross section:'

         ! 1st - photon      < g | D | E, icb >
         


         dmx2phb = 0.0D+00
             
         DO irb = 1,  nefb       !    | f > = | irb, enf(irb), nof(irb) >

            w_photon(irb) = 0.25D+00 * ( enfb(irb) - energy_initialb )
            
            DO i = 1, nofb(irb)
                  DO icb = 2,  2     !   | i > = | icb, eni(icb)> 

                     de = 0.5D+00*energy_initial + w_photon(irb) - 0.5D+00*enib(icb)     !detuning
                     
                     dmx2phb(irb,i) = dmx2phb(irb,i) +  dmxbfb(irb,i,icb)
!                     dmx2phb(irb,i) = dmx2phb(irb,i) +  dmxbb(icb) * dmxbfb(irb,i,icb) / de
                  ENDDO ! neib
         ENDDO !nofb
      END DO  !nef 

!      WRITE(*,*) dmxbb

         !  note: j = initial states, i = final states

         
   CLOSE(NBB)
   CLOSE(NBF1)



   DO irb = 1,  nefb       ! final states looping (irb, nefb)

      eph = w_photon(irb)

      IF(MODE.EQ.0) THEN 
         DWPH = EPH**(-NPHOTONS)
      ELSE IF(MODE.EQ.1) THEN
         DWPH = EPH**(NPHOTONS)
      ELSE IF(MODE.EQ.2) THEN
         DWPH = EPH**(-NPHOTONS)
      ENDIF

      cs2phbb = 0.0D+00
      DO i = 1, nofb(irb)
         cs2phbb = cs2phbb + ( REAL(dmx2phb(irb,i))**2 + AIMAG(dmx2phb(irb,i))**2)*cs_n_si*dwph
      ENDDO ! nof

      !      cs2phb(irb) = cs_n_si*dwph*SUM( dmx2phb(irb,:)**2, dim=2) 

      WRITE(30,'(1x,e15.4,2x,e20.5)') w_photon(irb)*enau, cs2phbb

ENDDO


   WRITE(*,*) '# cs2ph:: bound part of cross section has been calculated.'


!!!!           END OF BOUND PART OF CROSS SECTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!          CONTINUUM PART OF CROSS SECTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




   WRITE(*,*) '# cs2ph:: continuum part of cross section:'


!!!
!!!     ***********     temporal implementation    ***********
!!! Here it is assumed that the energy step of the intermediate states are the 
!!! same as the final, though is not necessary.
!!!


   DE_STEP = 0.014D+00/2.0D+00 !!ABS( enib(21) - enib(20) )



   WRITE(*,*)  '# cs2ph::                          de_step  = ', DE_STEP

!!!.............................................
!!!xx  STAR LOOPING OVER THE INTERMEDIATE STATES IN TERMS OF ENERGY (ic,nei)

         cs2ph = 0.0D+00

!         DO ic = 1,  nei
         DO ic = 1,  0



            !   STATE  | I > = | ic, eni(ic), noi(ic) >
            
!!!..........................................


            READ(NBF) ii,  eni(ic), noi(ic)
            READ(NBF) (dmxbf(j), j = 1, noi(ic))


!!....


            energy    = eni(ic) 
            nchannels = noi(ic)  


!!......................


            READ(NFF) ici, eni(ic), noi(ic)

!!!......................................................................
         
            IF(ABS(eni(ic) - energy).GT.1.0D-08) THEN

            WRITE(*,*) ' NOT CONSISTENT FILES '
            WRITE(*,*) '        II = ', ii
            WRITE(*,*) ' BF ENERGY = ', energy 
            WRITE(*,*) '        IC = ', ic
            WRITE(*,*) ' FF ENERGY = ', eni(ici)
            STOP

         ENDIF

!!!...............

         IF(noi(ic).NE.nchannels) THEN

            WRITE(*,*) ' NOT CONSISTENT FILES '
            WRITE(*,*) '          II  = ', ii
            WRITE(*,*) ' BF CHANNELS  = ', nchannels 
            WRITE(*,*) '          IC  = ', ic
            WRITE(*,*) ' FF CHANNELS  = ', noi(ic)
            STOP
         ENDIF

!!!................


!!         WRITE(NDIPOLE, '(i5,2x,e14.6,2x,i5)')  ic, eni(ic), noi(ic)


!!         WRITE(*,*),'        IC = ', ICI
!!         WRITE(*,*),'    ENERGY = ', eni(ic)
!!         WRITE(*,*),' NCHANNELS = ', noi(ic)

!!!.........................................


!!!xx  STAR LOOPING OVER THE FINAL STATES IN TERMS OF ENERGY (ir,nef)
!!..................................................................


         DO ir = 1,  nef

            !   STATE  | F > = | ir, enf(ir), nof(ir) >

!!!..............................................



         READ(NFF)  irf, enf(ir), nof(ir)
         READ(NFF) ( (dmxff(j,i), j = 1, nof(ir)), i = 1, noi(ic) )


!!!..............................................
!!!         DE = 0.5D+00 * (ENERGY_INITIAL + enf(ir) ) - eni(ic) 
!!!         WRITE(NDIPOLE, *)   ir, enf(ir), nof(ir)
!!!         WRITE(*,*),'        IR = ', IR
!!!         WRITE(*,*),'    ENERGY = ', enf(ir)
!!!         WRITE(*,*),' NCHANNELS = ', nof(ir)


!!!             ir --> Final Energy
!!!             i  --> i-th channel at Final Enrgy ir  cs(ir,i)

!!!...................................................................


         W_PHOTON(ir) = 0.5D+00 * ( enf(ir) - ENERGY_INITIAL )

         DE = ENERGY_INITIAL + W_PHOTON(ir) - eni(ic)     !detuning
      

         DO i = 1, nof(ir)

            DO j = 1, noi(ic)


               IF( ABS(DE).LE.ABS(20*DE_STEP/2.0D+00)) THEN


                  IF(DE.LT.0.D+00) THEN

                  cs2ph(ir,i) = cs2ph(ir,i) !!- M_PI * uim *  dmxbf(j) * dmxff(j,i) * 2.0D+00/DE_STEP

                  ELSE

                  cs2ph(ir,i) = cs2ph(ir,i) !!+ M_PI * uim *  dmxbf(j) * dmxff(j,i) * 2.0D+00/DE_STEP

                  ENDIF


                  WRITE(*,*) ' STEP           IR = : ', IR  
                  WRITE(*,*) ' INITIAL STATE E_I = : ', ENERGY_INITIAL
                  WRITE(*,*) ' POLE    STATE E_P = : ', ENI(IC)
                  WRITE(*,*) ' PHOTON        W   = : ', W_PHOTON(IR)
                  WRITE(*,*) '                DE = : ', DE
                  WRITE(*,*) ' DE_STEP/2         = : ', DE_STEP/2.0D+00
                  WRITE(*,*) ' ####################################################'



               ELSE

                  cs2ph(ir,i) = cs2ph(ir,i) + dmxbf(j) * dmxff(j,i)/DE
               
               ENDIF

!!               WRITE(*,*) 'cs2ph = ', j, cs2ph(ir,i), dmxbf(j), dmxff(j,i) 

            END DO
         ENDDO

         !  note: j = initial states, i = final states

!!!nef 
      END DO

!!! nei
   END DO



!!   weight = 0.5D+00 * ( eni(3) - eni(2) ) 
   weight = DE_STEP/2.0D+00

   cs2ph  = weight  * cs2ph
   cs2phb = 2.0D+00 * cs2phb

!!!............................

         
   CLOSE(NFF)
   CLOSE(NBF)

   
!!!..............

   OPEN(NCS2PH, FILE='DAT/CS2PH_TOTAL.DAT')


   DO ir = 1, nef
      
      WRITE(NDIPOLE,*) ir, enf(ir), nof(ir)

      cs2ph_E  = 0.0D+00
      cs2ph_Eb = 0.0D+00

      DO i = 1, nof(ir)


         dmx_2_channel  =   REAL( cs2ph(ir, i) )**2  + AIMAG( cs2ph(ir,i)  )**2
         dmx_2_channelb =   REAL( cs2phb(ir, i) )**2 + AIMAG( cs2phb(ir,i) )**2
         dmx_2          = dmx_2_channelb + dmx_2_channel


         IF(MODE.EQ.0) THEN 

            cs2ph_Eb = cs2ph_Eb + cs_2 * dmx_2_channelb / (0.5D+00 * w_photon(ir) )**2 
            cs2ph_E  = cs2ph_E  + cs_2 * dmx_2_channel  / (0.5+00  * w_photon(ir) )**2 


            WRITE(NDIPOLE, *) cs_2 * dmx_2 /w_photon(ir)**2 , cs2phb(ir,i) + cs2ph(ir,i)
         
         ELSE

            cs2ph_Eb = cs2ph_Eb + cs_2 * dmx_2_channelb * (0.5D+00 * w_photon(ir) )**2
            cs2ph_E  = cs2ph_E  + cs_2 * dmx_2_channel  * (0.5D+00 * w_photon(ir) )**2


            WRITE(NDIPOLE, *) cs_2 * dmx_2 * w_photon(ir)**2 , cs2phb(ir,i) + cs2ph(ir,i)
            
         ENDIF


      ENDDO


      WRITE(NCS2PH,'(4e14.6)') 0.5D+00 * EnAU * w_photon(ir)   &
           &      , cs_2_si * ( cs2ph_Eb + cs2ph_E )  & 
           &      , cs_2_si * cs2ph_Eb                & 
           &      , cs_2_si * cs2ph_E  

      
   ENDDO


!!!............................


   CLOSE(NDIPOLE)
   CLOSE(NCS2PH)

!!!............................
   
   DEALLOCATE(cs2ph)
   DEALLOCATE(dmxbf)
   DEALLOCATE(dmxff)
   DEALLOCATE(umx)

   DEALLOCATE(en)
   DEALLOCATE(no)
   DEALLOCATE(nhf)
   DEALLOCATE(lhf)
   DEALLOCATE(ll)
   DEALLOCATE(nmin)
   DEALLOCATE(nmx)
   DEALLOCATE(nd)
   DEALLOCATE(nol)
   DEALLOCATE(is)
   DEALLOCATE(idcs)

   DEALLOCATE(eni)
   DEALLOCATE(noi)
   DEALLOCATE(nhfi)
   DEALLOCATE(lhfi)
   DEALLOCATE(lli)
   DEALLOCATE(nmini)
   DEALLOCATE(nmxi)
   DEALLOCATE(ndi)
   DEALLOCATE(noli)
   DEALLOCATE(isi)
   DEALLOCATE(idcsi)


   DEALLOCATE(enf)
   DEALLOCATE(nof)
   DEALLOCATE(nhff)
   DEALLOCATE(lhff)
   DEALLOCATE(llf)
   DEALLOCATE(nminf)
   DEALLOCATE(nmxf)
   DEALLOCATE(ndf)
   DEALLOCATE(nolf)
   DEALLOCATE(isf)
   DEALLOCATE(idcsf)



END PROGRAM cs2ph_f
!!!#####################################################
!!!EOF

