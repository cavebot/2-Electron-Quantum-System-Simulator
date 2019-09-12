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

  USE configuration
  USE ioroutines

  IMPLICIT NONE

  CHARACTER(LEN=16) GAUGE
  INTEGER EnAU
  REAL(8) PI, ALPHA, CS_2, CS_2_SI 

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
  REAL(8) dmx_2, DE, WEIGHT


  INTEGER, ALLOCATABLE, DIMENSION(:)   :: no, noi, nof
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: nob, noib, nofb
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: umx
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: en, eni, enf, W_PHOTON
  REAL(8), ALLOCATABLE, DIMENSION(:)   :: enibb, enib, enfb


  REAL(8) dmxbb
  COMPLEX(8), ALLOCATABLE, DIMENSION(:)   :: dmxbf, dmxbfb
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: dmxff
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: cs2ph
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: cs2phb
  COMPLEX(8)  uim

!!!...................................


  EnAU  = 27.211396181D+00
  PI    = 2.0D+00 * ASIN(1.0D+00)
  ALPHA = 1.0D+00 / 137.0359D+00
  CS_2    = 2 * PI * 8 * PI * PI * PI * ALPHA**2
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


      OPEN(UNIT=9,file='INP/CS2PH.INP',status='old')

      READ(UNIT=9,*) LINITIAL, LFINAL
      READ(UNIT=9,*) GAUGE

      CLOSE(9)


      OPEN(UNIT=NOUTFILE,FILE='OUT/CS2PH.LOG')

!!!..........................

      IF(gauge.EQ.'V') THEN

         WRITE(*,*) ' TWO-PHOTON MULTICHANNEL GENERALIZED CROSS SECTONS IN VELOCITY GAUGE '

         MODE = 0     
      ELSE
            
         WRITE(*,*) ' TWO-PHOTON MULTICHANNEL GENERALIZED CROSS SECTONS IN LENGTH GAUGE '

         MODE = 1
      ENDIF


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


!!$      CALL CFGFILE(NCFG, LINTERMEDIATE)
!!$      
!!$      READ(NCFG, 3) nSymmetries
!!$      READ(NCFG, 3) l, ls
!!$      READ(NCFG, 3) ncsm
!!$      
!!$      CLOSE(NCFG) 
!!$      ncs = ncsm



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

      

      WRITE(*,*) ' ######################################################################'
      WRITE(*,*) ' READ DATA FOR INTERMEDIATE  2e - CONTINUUM STATES L = ', LINTERMEDIATE
         

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


      WRITE(NOUTFILE, *) ' NUMBER OF CONFIGURATIONS FOR INITIAL STATE = ', nconfigi


!!!.......... final states

      WRITE(*,*) ' ######################################################################'
      WRITE(*,*) ' READ DATA FOR FINAL 2e - CONTINUUM STATES L = ', LFINAL



      lof = LFINAL

      CALL cxfin(NCFG,nhff,lhff,llf,nminf,nmxf,nolf,isf,lof,lsf,ncsmxf,ntotf,ndf,idcsf )

      !..  total number of configurations for final states

      nconfigf = 0
      DO i = 1, ncsmxf
         DO j = 1, ndf(i)   
               
            nconfigf = nconfigf + 1
            
         END DO
      END DO

      WRITE(NOUTFILE, *) ' NUMBER OF CONFIGURATIONS FOR FINAL STATE', nconfigf
      WRITE(*,*) ' ######################################################################'


!!!.....................................  checks

 
      IF(lsf.NE.lsi) THEN
   
         WRITE(*,*) 'NOT DIPOLE ALLOWED TRANSITION. STATES HAVE DIFFERENT PARITY'
         WRITE(*,*) 'INITIAL SPIN MOMENTUM ::  LSI = ',  LSI
         WRITE(*,*) 'FINAL  SPIN  MOMENTUM ::  LSF = ',  LSF
         STOP

      ELSE 

         IF(lsf.EQ.1) THEN
   
            WRITE(NOUTFILE,*) ' SINGLET SYMMETRY'
         ELSE IF(lsf.EQ.3)   THEN

            WRITE(NOUTFILE,*) ' TRIPLET SYMMETRY'
         ELSE
            
            WRITE(NOUTFILE,*) ' ALLOWED VALUES OF ISL IN r12.inp ARE 1 (Singlet), OR 3 (triplet) '
         ENDIF

      ENDIF

      WRITE(*,*) ' ######################################################################'
      WRITE(*,*) 'INITIAL SPIN MOMENTUM ::  LSI = ',  LSI
      WRITE(*,*) 'FINAL  SPIN  MOMENTUM ::  LSF = ',  LSF

!!!..........................

      
      IF(ncsf.GT.ncsi) THEN

         ndmax = ncsf
         
      ELSE IF(ncsi.GT.ncsf) THEN

         ndmax = ncsi

      ELSE

         ndmax = ncsi
         
      ENDIF

      WRITE(*,*), ' NDMAX = ', ndmax

!!!...................................


      WRITE(*,*) '######################################################################'


      ALLOCATE(umx(ndmax,ndmax))


      ALLOCATE(dmxbfb(ndmax))
      ALLOCATE(dmxbf(ndmax))
      ALLOCATE(dmxff(ndmax,ndmax))


!!!......................................


      CALL  BFDMX2EFILE(NBF, LINITIAL, LINTERMEDIATE, MODE)  ! 0 --> 1
      CALL  BBDMX2EFILE(NBB, LINITIAL, MODE)                 ! 0 --> 1
      
      CALL  FFDMX2EFILE(NFF,  LINTERMEDIATE, LFINAL, MODE)   ! 1 -- > 0 (2) 
      CALL  BFDMX2EFILE(NBF1, LINTERMEDIATE, LFINAL, MODE)   ! 1 -- > 0 (2) 



      WRITE(*,*) ' READ FOR  FIRST  PHOTON'
      WRITE(*,*) ' OPENING FILE ', LFINAL, ' TO ', LINTERMEDIATE


      WRITE(*,*) ' ######################################################################'


      WRITE(*,*) ' READ FOR 2-ND PHOTON'
      WRITE(*,*) ' OPENING FILE ', LFINAL, ' TO ', LINTERMEDIATE


      WRITE(*,*) ' ######################################################################'



      OPEN(NDIPOLE, file="DAT/CS2PH.DAT")

!!!..............

!!!       READ FOR FIRST PHOTON                


!!!......... 


      !    BOUND 


      READ(NBB)  MODE 
      READ(NBB)  lo, loi, n, ni



      WRITE(*,*) ' BB : 0 -----> 1   lo, loi, n, ni = ', lo, loi, n, ni 

      modebf1 = mode

      ALLOCATE(enibb(ni)) 

      READ(NBB) ENERGY_INITIALB
      READ(NBB) (enibb(ni), icib = 1, ni)

      WRITE(*,*) ' ENERGY_INITIALB = ', ENERGY_INITIALB


      !
      !      do nsr = 1, nsi
      !      READ(NBB) ( dmx(nsc, nsr, 1), nsc = 1, nsm)
      !   ENDDO
      !      CLOSE(NBB)

      
!!!...............

      
      !   BOUND - FREE 
      

!!!       | G, L, 1 > ---------> | IC, L + 1, ( ic = 1 - nei ) >


         READ(NBF) MODE
         READ(NBF) ne, nei
         READ(NBF) ii, ENERGY_INITIAL

         nEn    = nei
         modebf = MODE



         WRITE(*,*) ' BF 0 ---- > 1  :  ne, nei = ', ne, nei
         WRITE(*,*) ' ENERGY_INITIAL = ', ENERGY_INITIAL


!!!..............

 
!!!     READ FOR SECOND PHOTON        FF
!!!     | IC, L + 1, ( ic = 1 - nei ) > -----> | IR, L +- 1, ( ir = 1 - nef) >


!!..........

         ! BOUND - FREE  PART  

         READ(NBF1) MODE
         READ(NBF1) neib, nefb

         WRITE(*,*) ' BF 1 ---- > 0  :  neib, nefb = ', neib, nefb


!!$         IF(ni.NE.neib) THEN 
!!$
!!$            WRITE(*,*) ' NOT CONSISTENT FILES. STOP '
!!$            WRITE(*,*) ' BB TRANSITION    TO ', ni, ' INTERMEDIATE ENERGIES'
!!$            WRITE(*,*) ' BF TRANSITION  FROM ', neib, ' INTERMEDIATE ENERGIES'
!!$            STOP
!!$            
!!$         ELSE
!!$
!!$            WRITE(*,*) ' NUMBER OF BOUND INTERMEDIATE ENERGIES  neib = ', neib
!!$
!!$         ENDIF

         
         IF(mode.NE.modebf1) THEN 

            WRITE(*,*) ' NOT CONSISTENT FILES. STOP '
            WRITE(*,*) ' BB TRANSITION  MODE = ', modebf1 
            WRITE(*,*) ' BF TRANSITION  MODE = ', mode  
            STOP
            
         ELSE

            IF(MODE.EQ.0) THEN

               WRITE(*,*) ' 2-PHOTON CROSS SECTIONS IN VELOCITY GAUGE '

            ELSE

               WRITE(*,*) ' 2-PHOTON CROSS SECTIONS IN VELOCITY GAUGE '
            ENDIF


         ENDIF



!!!..........

         ! FREE - FREE  PART


         READ(NFF) MODE
         READ(NFF) nei, nef


!!.....

         IF(nei.NE.nEn) THEN 

            WRITE(*,*) ' NOT CONSISTENT FILES. STOP '
            WRITE(*,*) ' BF TRANSITION    TO ', nEn, ' INTERMEDIATE ENERGIES'
            WRITE(*,*) ' FF TRANSITION  FROM ', nei, ' INTERMEDIATE ENERGIES'
            STOP
            
         ELSE

            WRITE(*,*) ' NUMBER OF INTERMEDIATE ENERGIES  nEn = ', nEn

         ENDIF

!!!...........

         IF(nefb.NE.nef) THEN 

            WRITE(*,*) ' NOT CONSISTENT FILES. STOP '
            WRITE(*,*) ' BF TRANSITION    TO ', nefb, ' INTERMEDIATE ENERGIES'
            WRITE(*,*) ' FF TRANSITION  FROM ', nef, ' INTERMEDIATE ENERGIES'
            STOP
            
         ELSE

            WRITE(*,*) ' NUMBER OF FINAL  ENERGIES  n = ', nefb

         ENDIF

!!!..............

         IF(modebf.NE.mode) THEN 

            WRITE(*,*) ' NOT CONSISTENT FILES. STOP '
            WRITE(*,*) ' BF TRANSITION  MODE = ', modebf 
            WRITE(*,*) ' FF TRANSITION  MODE = ', mode  
            STOP
            
         ELSE

            IF(MODE.EQ.0) THEN

               WRITE(*,*) ' 2-PHOTON CROSS SECTIONS IN VELOCITY GAUGE '

            ELSE

               WRITE(*,*) ' 2-PHOTON CROSS SECTIONS IN VELOCITY GAUGE '
            ENDIF


         ENDIF

         NE = 1         
         WRITE(*,*) MODE, NE,NEI,NEF
   
         WRITE(*,*) ' ######################################################################'

!!!........................................................................
         

         WRITE(NDIPOLE,*) MODE
         WRITE(NDIPOLE,*) NEF

!!!.........................................


         ALLOCATE(en(ne))
         ALLOCATE(no(ne))

         
         WRITE(*,*) '==> ALLOCATIONS FOR INTERMEDIATE ENERGIES,  NEI X NEI = ', nei

         ! intermediate


         ! (B)

         ALLOCATE(enib(neib))
         ALLOCATE(noib(neib))

         ! (C)

         ALLOCATE(eni(nei))
         ALLOCATE(noi(nei))


         ! final

         WRITE(*,*) '==> ALLOCATIONS FOR FINAL ENERGIES,                NEF X NEF   = ', nef
         WRITE(*,*) ' ######################################################################'

         ALLOCATE(enfb(nefb))
         ALLOCATE(nofb(nefb))
         ALLOCATE(cs2phb(nef, ndmax))

         ALLOCATE(enf(nef))
         ALLOCATE(nof(nef))
         ALLOCATE(cs2ph(nef, ndmax))

         ALLOCATE(W_PHOTON(nef))


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



         WRITE(*,*) ' ==> BOUND PART OF CROSS SECTIONS :: '



!!!.............................................
!!!xx  START LOOPING OVER THE INTERMEDIATE STATES IN 
!!!    TERMS OF ENERGY (icb,neib)



         cs2phb = 0.0D+00

         DO icb = 1,  neib

      
            ! 1st - photon      < g | D | E, icb >
            !!..........................



            READ(NBB) dmxbb 


            
            !!............................!
            !! STATE  | I > = | ic, eni(ic), noi(ic) > 



            READ(NBF1) icib,  enib(icb)   !, noib(icb)





!!!xx  STAR LOOPING OVER THE FINAL STATES IN TERMS OF ENERGY (irb,nefbf)
!!..................................................................


         DO irb = 1,  nefb


!   STATE  | F > = | irb, enf(irb), nof(irb) >

            READ(NBF1) irb,  enfb(irb), nofb(irb)
            READ(NBF1) (dmxbfb(j), j = 1, nofb(irb))


!!!             ir --> Final Energy
!!!             i  --> i-th channel at Final Enrgy ir  cs(ir,i)

!!!...................................................................



            W_PHOTON(irb) = 0.5D+00 * ( enfb(irb) - ENERGY_INITIALB )

            DE = ENERGY_INITIAL + W_PHOTON(irb) - enib(icb)     !detuning


         DO i = 1, nofb(irb)
            DO j = 1, noib(icb)

               cs2phb(irb,i) = cs2ph(irb,i) + dmxbb * dmxbfb(j) / DE

               
!!               WRITE(*,*) 'cs2ph = ', j, cs2ph(ir,i), dmxbf(j), dmxff(j,i) 
            END DO
         ENDDO

         !  note: j = initial states, i = final states

!!!nef 
      END DO

!!! nei
   END DO


         
   CLOSE(NBB)
   CLOSE(NBF1)





!!!!           END  OF  BOUND  PART  OF   CROSS SECTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!          CONTINUUM  PART  OF  CROSS  SECTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



         WRITE(*,*) ' ==> BOUND PART OF CROSS SECTIONS :: '



!!!.............................................
!!!xx  STAR LOOPING OVER THE INTERMEDIATE STATES IN TERMS OF ENERGY (ic,nei)

         cs2ph = 0.0D+00

         DO ic = 1,  nei



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

               cs2ph(ir,i) = cs2ph(ir,i) + dmxbf(j) * dmxff(j,i)/DE
               

!!               WRITE(*,*) 'cs2ph = ', j, cs2ph(ir,i), dmxbf(j), dmxff(j,i) 
            END DO
         ENDDO

         !  note: j = initial states, i = final states

!!!nef 
      END DO

!!! nei
   END DO



         weight = 0.5D+00 * ( eni(3) - eni(2) ) 

         cs2ph  = sqrt(weight) * cs2ph
         

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


         dmx_2_channel  =   REAL( cs2ph(ir, i) )**2 + IMAG( cs2ph(ir,i) )**2
         dmx_2_channelb =   REAL( cs2phb(ir, i) )**2 + IMAG( cs2phb(ir,i) )**2
         dmx_2          = dmx_2_channelb + dmx_2_channel


         IF(MODE.EQ.0) THEN 

            cs2ph_Eb = cs2ph_Eb + cs_2 * dmx_2_channelb / w_photon(ir)**2 
            cs2ph_E  = cs2ph_E  + cs_2 * dmx_2_channel  / w_photon(ir)**2 


            WRITE(NDIPOLE, *) cs_2 * dmx_2 /w_photon(ir)**2 , cs2phb(ir,i) + cs2ph(ir,i)
         
         ELSE

            cs2ph_Eb = cs2ph_Eb + cs_2 * dmx_2_channelb * w_photon(ir)**2
            cs2ph_E  = cs2ph_E  + cs_2 * dmx_2_channel  * w_photon(ir)**2


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

