!!xxx
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
!!!          
!!!           note  (2)  NH2E has to be 1 as has been implemented now
!!!                           for single ionization files.
!!!                           and  2 for double ionization files.   
!!!  
!!!
!!!  Calculates un-normalized 2-e DME from a continuum state | c_i > to 
!!!  a multichannel continuum (single or double) | c_f >
!!!  
!!!  | c_i > ----> | c_f >       
!!!
!!!  ! c_i >  2-e continuum state calculated on the mixed basis states (free)
!!!  | c_f >  2-e continuum state calculated on the mixed basis states (free) 
!!!
!!!
!!!   nhmx : nhmx > total number of configurations ( how much ?)
!!!   ncsi : # configuration series for the initial state
!!!   ncsf : # configuration series for the final   states
!!!   nsx  : number of orbitals <=(?)    nsx  < = number of B-splines
!!!............................

program hdmx2eff

!  use jl

  USE wf_files
  USE dz_value
  USE configuration
  USE ioroutines

  IMPLICIT NONE
  INTEGER, PARAMETER :: dpk = KIND(1.d0)
  INTEGER j1, j2, ls
  COMMON/ang/j1,j2,ls
  CHARACTER(len=16)  configi, configf, outfile, gauge
  integer lof, lsf, ncsmxf, loi,lsi, ncsmxi, nconfigi, nconfigf
  INTEGER nhmx, nbsp, ncs, ncsi, ncsf, nsi, nsf, ntotf, ntoti, nei, nef, ltotalMax
  INTEGER kf, ki, ir, ic, nr, nc, i, j, it, n, m, nout, iout
  INTEGER LINITIAL, NMININ, NMAXI, LFINAL, NMINFIN, NMAXF, k
  INTEGER MODE, NOUTFILE, NWF2E, NH2EI, NH2EF, NCFG, NHDMX2E
  INTEGER incx, incy
  REAL(dpk) alpha, beta
  INTEGER, ALLOCATABLE, DIMENSION(:) :: noi, nof
  INTEGER, DIMENSION(4) :: id
  REAL(dpk), DIMENSION(4) :: ang
  REAL(dpk) fli, flf, phase
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: dmxz, cifin, ciin, dpz
  REAL(dpk), ALLOCATABLE, DIMENSION(:) :: eni, enf, osz, tt
  CHARACTER(LEN=100) ARGV 
  EXTERNAL dgemv
  
!!!..............................................

!!! NH2EI = 1, SINGLE IONIZATION CONTINUUM
!!! NH2EF = 2, SINGLE IONIZATION CONTINUUM

!!! NH2EI > 2, DOUBLE IONIZATION CONTINUUM
!!! NH2EF > 2, DOUBLE IONIZATION CONTINUUM

      NWF2E    = 3         
      NHDMX2E  = 4
      NOUTFILE = 16
      NCFG     = 18


      CALL GETARG(1, ARGV)
      READ(ARGV,*) LINITIAL         
      CALL GETARG(2, ARGV)
      READ(ARGV,*) LFINAL       
      CALL GETARG(3, ARGV)
      READ(ARGV,*) GAUGE


!!!..........................

      OPEN(9,file='INP/HDMX2EFF.INP',status='old')
!      READ(UNIT=9,*) NMININ,  NMAXI
!      READ(UNIT=9,*) NMINFIN, NMAXF
      READ(9,*) NH2EI,   NH2EF
      CLOSE(9)

      LTOTALMAX = 4 ! It is connected with nm=5 (h1e)
!.............................................

      MODE = 1
      IF(gauge.EQ.'V')          MODE = 0     

      WRITE(*,*) '#'
      WRITE(*,*) '# dmx2eff::                        dmx in ', gauge, ' gauge'
      WRITE(*,*) '# dmx2eff::       initial symmetry   li = ', LINITIAL
      WRITE(*,*) '# ndmx2eff::               nmini, nmaxi = ', nminin,   nmaxi
      WRITE(*,*) '# dmx2eff::       final   symmetry   lf = ', LFINAL
      WRITE(*,*) '# ndmx2eff::               nminf, nmaxf = ', nminfin, nmaxf
      WRITE(*,*) '#'

!!!..........................

      nout = 1
      DO iout = 1, nout


!!!... open outpout file
      
         OPEN(NOUTFILE, file="OUT/HDMX2EFF.LOG")

!!!.......... READ CONFIGURATION FILE FOR INITIAL FREE STATE



         ! NCSI # NUMBER OF CONFIGURATION SERIES FOR INITIAL STATE
         ! NCSF # NUMBER OF CONFIGURATION SERIES FOR FINAL   STATE


         CALL GETNUMBEROFCHANNELS(NCFG, LINITIAL, NCSI, NBSP)

         CALL GETNUMBEROFCHANNELS(NCFG, LFINAL, NCSF, NBSP)

         CALL CONFIG_SPACE(ncsi, ncsf)


         loi = LINITIAL

         CALL cxfin(NCFG, nhfi,lhfi, lli, nmini, nmxi, noli,&
              &        isi, loi, lsi, ncsmxi,ntoti, ndi, idcsi )

         
         nconfigi = 0

!!!.... CALCULATE TOTAL NUMBER OF CONFIGURATIONS (FREE+BOUND CHANNELS)

         DO i = 1, ncsmxi
            DO j = 1, ndi(i)   

               nconfigi = nconfigi + 1

            end do
         end do


         WRITE(*, *) '# dmx2eff::                  ncsi = ', nconfigi
         WRITE(*, *) '#'



         FLI = dfloat( 2*loi + 1 )

!!!.......... READ CONFIGURATION FILE FOR FINAL FREE STATE




         lof = LFINAL

         CALL cxfin(NCFG, nhff,lhff, llf, nminf, nmxf, nolf,&
              & isf, lof, lsf, ncsmxf,ntotf, ndf,idcsf )


         nconfigf = 0
         DO i = 1, ncsmxf
            DO j = 1, ndf(i)   
               
            nconfigf = nconfigf + 1
            
         END DO
      END DO

      WRITE(*, *) '# dmx2eff::                     ncsf = ', nconfigf
      WRITE(*, *) '#'


      FLF = dfloat( 2 * lof + 1 )

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

     
      WRITE(*,*) '# dmx2eff::      initial spin     si = ',  LSI
      WRITE(*,*) '# dmx2eff::      final   spin     sf = ',  LSI
      WRITE(*,*) '#'


!!!..........................


!!!
!!!     DON'T REMOVE. USED BY ROUTINE ANGLS FOR CALCULATION 
!!!
 
      j1 = lof
      j2 = loi
      ls = (lsi - 1)/2

!!!      IF(lsi.NE.lsf) STOP

      IF(ls.EQ.0)  phase =   1.0D+00
      IF(ls.EQ.1)  phase = - 1.0D+00

!      WRITE(*,*) ,'# dmxfileff::              phase= ', phase

!!!.........................

      nhmx = MAX(nconfigi, nconfigf)

      WRITE(*,*) '# ndmx2eff::   dimension of dmx    nhmx = ', nhmx

! or ?       ALLOCATE(dmxz(nconfigi, nconfigf))

      ALLOCATE(dmxz(nhmx, nhmx))

      dmxz = 0.0D+000


!!!.................   START READING 1-E DME



      ALLOCATE(dz(nbsp, nbsp))   

!!!
!!!  read    <  P_nl   |  d  | P_n'l'   >      0 <= l, l' < ltotalMax
!!!          <  P_nl   |  d  | B_n'(l') >      
!!!          <  B_n(l) |  d  | P_n'l'   >
!!!          <  B_n    |  d  | B_n      >
!!!


      CALL read_dp(nbsp, ltotalMax, mode) 



      WRITE(*,*) '# dmx2eff:: read dmx1e done.'

      
!!!..................................
!!!
!!!
!!!     | a, E,L  >   = S_{m_a,n'_a} c(m_a,n'_a) * F(m_a, n'_a) + S_{i_a}   c_i * X(i_a)  
!!!                         ----  free part ----------       ---- bound part ---      
!!!
!!!                   =   F(a) + B(a) 
!!!
!!!     | b, E',L' > = S_{m_b,n'_b} c(m_b,n'_b) * F(m_b, n'_b) + S_{i_b}  c_i * X(i_b)  
!!!                         ----  free part -------           ---- bound part --- 
!!!
!!!                   =   F(b) + B(b)
!!!
!!!
!!!     <a,E,L | D | b, E,L > =    < F(a) | D | F(b) >  +  < F(a) | D | B(b) > 
!!!        
!!!                                       FF                        FB     
!!!               
!!!                              + < B(a) | D | F(b) >  +  < B(a) | D | B(b) > 
!!!   
!!!                                       BF                        BB
!!!
!!!
!!!    Each of those terms includes four subterms because of antisymmetrization of 
!!!                        the 2-e wavefunction
!!!..............................



      kf = 1
      DO ir = 1, ncsmxf

         ki = 1

         DO  ic = 1, ncsmxi

!!!..................... CALCULATE THE ANGULAR FACTOR FOR THE CASES 1,2,3,4 


!!!
!!!     NOTE THAT THE ANGULAR FACTOR IS THE SAME FOR ALL CASES
!!!

!!$            WRITE(*,*)  '             IR, IC = ', IR , IC 
!!$            WRITE(*, *) '          lhff, llf = ', lhff(ir), llf(ir) 
!!$            WRITE(*, *) '          lhfi, lli = ', lhfi(ic), lli(ic) 
!!$            WRITE(*, *) ' L_FINAL, L_INITIAL = ', lof, loi


            CALL angfc( lhff(ir), llf(ir), lhfi(ic), lli(ic), lof, loi, id, ang)


!!$            WRITE(*,*) ' ANG(1), ID(1) = ', ang(1), id(1)  
!!$            WRITE(*,*) ' ANG(2), ID(2) = ', ang(2), id(2) 
!!$            WRITE(*,*) ' ANG(3), ID(3) = ', ang(3), id(3) 
!!$            WRITE(*,*) ' ANG(4), ID(4) = ', ang(4), id(4) 



!!!.....................


            SELECT CASE( 2 * idcsf(ir) + idcsi(ic) )

            CASE(0)      ! bb 2e-dmx


!               WRITE(*,*) ' MULTICHANNEL DME , BOUND -BOUND  ' , 2 * idcsf(ir) + idcsi(ic)
               dz = 0.0D+00

               CALL cal_nlnl(ir, ic, id, ang, phase, mode)

            CASE(1)    ! fb 2e-dmx 



!               WRITE(*,*) ' MULTICHANNEL DME , FREE - BOUND  ' , 2 * idcsf(ir) + idcsi(ic)
               dz = 0.0D+00

               CALL cal_spnl(ir, ic, id, ang, phase, mode)

            CASE(2)   ! bf 2e-dmx 

!               WRITE(*,*) '  MULTICHANNEL DME , BOUND - FREE   ' , 2 * idcsf(ir) + idcsi(ic)

               dz = 0.0D+00

               CALL cal_nlsp(ir, ic, id, ang, phase, mode)

            CASE(3)  ! ff 2e-dmx

!               WRITE(*,*) '  MULTICHANNEL DME , FREE - FREE   ' , 2 * idcsf(ir) + idcsi(ic)

               dz = 0.0D+00

               CALL cal_spsp(ir, ic, id, ang, phase, mode)

            END SELECT

!!!           WRITE(*,*) ' dz() = ', ( dz(i,1),i = 1, ndf(ir) )
!!!           WRITE(*,*) ' dz() = ', ( dz(i,2),i  =1, ndf(ir) )

            !      select case(ir.ne.1.or.ic.ne.1) 
            !        case(.false.)
            !        call mxprnt(dz,ndf(ir),ndi(ic))
            !      end select

!!!
!!!         dz ---> dmxz  (?)
!!!
!!!

            CALL trnsmx( dmxz, nhmx, dz, nbsp, kf, ki, ndf(ir), ndi(ic))

            ki = ki + ndi(ic)

         ENDDO

         kf = kf + ndf(ir)
            
      ENDDO

                WRITE(*,*) '# dmx2eff::       DMX2e(config,1),DMX2e(config,2) :'
                WRITE(*,'(5e14.6)')  (dmxz(j,1), j=1,5)
                WRITE(*,'(5e14.6)')  (dmxz(j,2), j=1,5)



         WRITE(*,*) '# dmx2eff: read dmx2e done.'


            DEALLOCATE(dz)   
            DEALLOCATE(vz)
            DEALLOCATE(vzba)
            DEALLOCATE(vzbb)
            DEALLOCATE(sp_a)
            DEALLOCATE(over)

!!!...............................


!!!... NORMALIZE CORRECTLY FOR EQUIVALENT ORBITALS 1/SQRT(2)

            it = 1
            SELECT CASE( MOD(lhfi(1)+lli(1),2) ) !!! initial state parity

            CASE(0)

!               PRINT*,'q2 in initial'

!!!...
               DO ic = 1, ncsmxi

                  SELECT CASE(idcsi(ic))

                  CASE(0)

                  IF( (lhfi(ic).EQ.lli(ic)).AND.(nhfi(ic).EQ.nmini(ic)) ) THEN 

!!!                     SELECT CASE(lhfi(ic).EQ.lli(ic).AND.nhfi(ic).EQ.nmini(ic))
!!!                     CASE(.TRUE.)

                        DO nr=1, nhmx

                           dmxz(nr,it) = dmxz(nr,it) / dsqrt(2.0D+00)

                        END DO

!!!                     END SELECT
                     ENDIF

                  END SELECT

                  it = it + ndi(ic)

               END DO
!!!...

            CASE(1)

!               PRINT*,'q2 in final'

               DO ir = 1, ncsmxf

                  SELECT CASE(idcsf(ir))

                  CASE(0)
                     IF( (lhff(ir).EQ.llf(ir)).AND.(nhff(ir).EQ.nminf(ir)) ) THEN

!!!                     SELECT CASE(lhff(ir).EQ.llf(ir).AND.nhff(ir).EQ.nminf(ir))
!!!                     CASE(.TRUE.)
   
                        DO nc = 1, nhmx

                           dmxz(it,nc) = dmxz(it,nc) / dsqrt(2.0D+00)

                        END DO

!!!END SELECT
                  ENDIF

               END SELECT

               it = it + ndf(ir)
               END DO

            END SELECT

            WRITE(*,*) '# dmx2eff:: equivalent orbitals normalization done.'

!!!............................



!!!# allocate data for initial-final states

            ALLOCATE(ciin( nconfigi, ncsmxi))
            ALLOCATE(cifin(nconfigf, ncsmxf))
            ALLOCATE(tt(nconfigf) )

!!!.................................
            
!!!#  read initial (continuum) states config. coef             


            CALL H2EFILE(NH2EI, LOI) 

            READ(NH2EI) nei
            CLOSE(NH2EI)            

            IF(NMAXI.GT.NEI)            NMAXI = NEI

            WRITE(*,*) '# dmx2eff::         initial states  nmaxi =', nmaxi


            ALLOCATE( eni( nei ) )
            ALLOCATE( noi( nei ) )

!!!........................................

!!!..........#  read final (continuum) states config. coef 


            CALL H2EFILE(NH2EF, LOF)
            READ(NH2EF) nef
            CLOSE(NH2EF) 
            
            IF(NMAXF.GT.NEF)      NMAXF = NEF
            WRITE(*,*) '# dmx2eff::         final states  nmaxf =', nmaxf            

            
            ALLOCATE(enf(nef))
            ALLOCATE(nof(nef))
            
            WRITE(*,*) '# dmx2eff::                       nei,nef =', nei,nef   

!!!.................................

      CALL  FFHDMX2EFILE(NHDMX2E, LOI, LOF, MODE) 


      IF(MODE.EQ.0) THEN
         WRITE(NOUTFILE,*) ' VELOCITY GAUGE '         
      ELSE
         WRITE(NOUTFILE,*) ' VELOCITY GAUGE '
      ENDIF


      WRITE(UNIT=NHDMX2E) MODE
      WRITE(UNIT=NHDMX2E) NEI, NEF


!      WRITE(UNIT=NHDMX2E) NMAXI, NMAXF 

!!!.............

      ALLOCATE( dpz(ncsmxf, ncsmxi))

!!!...........

!!!.............................................  initial states loop

!            ic : state index
!        en(ic) : energy for 'ic' state
!       noi(ic) : nof channels for 'ic' state
!
!       state  | i > = | ic, eni(ic), noi(ic) >


      
      CALL H2EFILE(NH2EI, LOI) 



      READ(NH2EI) nei

!!!      DO ic = 1, nei
         WRITE(*,*)'#'
         WRITE(*,*)'# dmx2eff:: index,      energy,    nof channels    i = '

         DO IC = 1, nei          ! initial states looping

         
            READ(NH2EI)  eni(ic), noi(ic)
            READ(NH2EI) ((ciin(n,j), n=1,nconfigi), j=1,noi(ic))



!!!               WRITE(*,*) ((ciin(n,j), n=1,nconfigi), j=1,noi(ic))
!!!................................               


!               IF(IC.GE.NMININ) THEN 

                  WRITE(*,'(14x,I4,1x,G20.8,1x,i4)') ic,eni(ic),noi(ic)

                  WRITE(UNIT=NHDMX2E) eni(ic), noi(ic) 
            !                  WRITE( NOUTFILE, '(e14.6,2x,i5)') eni(ic),noi(ic) 



!...................................... final states loop
!            ir : state index
!        en(ir) : energy for 'ir' state
!       noi(ir) : nof channels for 'ir' state
!
!       state  | F > = | ir, eni(ir), nof(ir) >

                  WRITE(*,*)'#'
                  WRITE(*,*)'# dmx2eff:: index,      energy,    nof channels    f = '

                  CALL H2EFILE(NH2EF, LOF) 
                  READ(NH2EF) nef

!                  DO ir = NMINFIN, NMAXF   ! final states loop

                  DO ir = 1, nef   ! final states loop

                     READ(NH2EF) enf(ir), nof(ir)
                     READ(NH2EF) ((cifin(n,j), n=1,nconfigf), j=1,nof(ir))

!!!                     WRITE(*,*) ((cifin(n,j), n=1,nconfigf), j=1,nof(ir))

                     WRITE(*,'(14x,I4,1x,G20.8,1x,i4)')ir,enf(ir),nof(ir)

                     WRITE(UNIT=NHDMX2E) enf(ir), nof(ir)

!                     WRITE(NOUTFILE, '(e15.7,2x,i5)') enf(ir), nof(ir)

!!! CALCULATE FINAL DIPOLE MOMENT < ic, noi | D | ir, nof >

                  alpha = 1.0D+00 
                  beta  = 0.0D+00
                  incx  = 1
                  incy  = 1

                  dpz = 0.0D+00
                  DO m = 1, noi(ic)
                     DO n = 1, nof(ir)

                        !  calculate the final dipole

                        CALL dgemv('n', nconfigf,nconfigi,alpha,dmxz,nhmx,ciin(:,m),incx,beta,tt,incy)

    
!!    WRITE(*,*) ' tt --------------->' 
!!    WRITE(*,*)  ( tt(j), j=1,5)         
!!    WRITE(*,*) ' cifin ------------> '     
!!    WRITE(*,*) ( cifin(j,n), j=1,5)             


!!!                        dpz(n,m) = 0.0D+00
!!!                        dpz(n,m) = DOT_PRODUCT(cifin(:,n), tt)


                        DO j = 1, nconfigf
                           
                           dpz(n,m) = dpz(n,m) + cifin(j,n) * tt(j)
                        END DO


                     END DO
                  END DO
            !  output

                  WRITE(NOUTFILE,'(4e15.7)') ((dpz(j,i), j=1,nof(ir)),i=1,noi(ic))

                  WRITE(UNIT=NHDMX2E) ((dpz(j,i), j = 1, nof(ir)), i = 1, noi(ic))

!                WRITE(*,*) '<',ic,n, '| D | ',ir,m,'> = ', dpz(n,m)

!            write(*,'(4e15.7)') ((dpz(j,i), j=1,nof(ir)),i=1,noi(ic))

               END DO
               
               CLOSE(NH2EF)               

!            ENDIF


         END DO

    CLOSE(NH2EI)


    DEALLOCATE(dpz)
    DEALLOCATE(eni)
    DEALLOCATE(enf)
    DEALLOCATE(ciin)
    DEALLOCATE(cifin)
    DEALLOCATE(noi)
    DEALLOCATE(nof)
    DEALLOCATE(tt)
    DEALLOCATE(dmxz)
    DEALLOCATE(nhfi)
    DEALLOCATE(nhff)
    DEALLOCATE(lhfi)
    DEALLOCATE(lhff)
    DEALLOCATE(lli)
    DEALLOCATE(llf)
    DEALLOCATE(nmini)
    DEALLOCATE(nminf)
    DEALLOCATE(nmxi)
    DEALLOCATE(nmxf)
    DEALLOCATE(ndi)
    DEALLOCATE(ndf)
    DEALLOCATE(noli)
    DEALLOCATE(nolf)
    DEALLOCATE(isi)
    DEALLOCATE(isf)
    DEALLOCATE(idcsi)
    DEALLOCATE(idcsf)
    
 END DO

!!!..........................................

!!! 577  format(5a16)

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

!!!..........................................


       END PROGRAM 
       
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

       SUBROUTINE trnsmx(hmx,nhmx,xm,nbsp,ihr,ihc,ir,ic)
         IMPLICIT NONE
         INTEGER nhmx, nbsp
         INTEGER ihr, ihc, ir,ic, irf, icf, kr, k, kc, kk
         REAL(8), DIMENSION(nbsp,nbsp) ::  xm
         REAL(8), DIMENSION(nhmx,nhmx) ::  hmx
         
         !  allocate(xm(nsx,nsx))
         !  allocate(hmx(nhmx,nhmx))
         
         irf=ihr+ir-1
         icf=ihc+ic-1
         DO 10 k=ihr,irf
            DO 10 kk=ihc,icf
               
               kr=k-ihr+1
               kc=kk-ihc+1
               hmx(k,kk)=xm(kr,kc)
               
10             CONTINUE
               
             END SUBROUTINE trnsmx

!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
       
       
