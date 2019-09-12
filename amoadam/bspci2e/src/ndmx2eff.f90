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

MODULE def 

  IMPLICIT NONE
  TYPE szeta
     REAL(8), POINTER, DIMENSION(:) :: zeta, phi
     COMPLEX(8), POINTER, DIMENSION(:,:) :: sx
  END TYPE szeta
END MODULE def

!!!################################################################
MODULE configuration
  IMPLICIT NONE
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhfi,lhfi,lli,nmini,nmxi,&
       & ndi,nhff,lhff,llf,nminf,nmxf, ndf,isi,isf, noli, nolf, idcsf, idcsi
CONTAINS
  SUBROUTINE config_space(ncsi, ncsf)  
    IMPLICIT NONE
    INTEGER ncsi, ncsf

!!!...............  INITIAL STATE

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
PROGRAM dpnormff

  USE def
  USE configuration
  USE ioroutines

  IMPLICIT NONE

  CHARACTER(LEN=16) GAUGE
  INTEGER EnAU
  REAL(8) PI 

  INTEGER i, j, m, ndmax
  INTEGER ncheck,  iout, nth, info
  INTEGER nSymmetries, l, ls, ncs
  INTEGER ic, nei, ni, lsi, ncsmxi, ntoti, loi, ncsi, nconfigi
  INTEGER ir, nef, nf, lsf, ncsmxf, ntotf, lof, ncsf, nconfigf
  INTEGER MODE, K
  INTEGER NDIPOLE,  NOUTFILE, NCFG, NHDMX2E, NNDMX2E
  INTEGER LINITIAL, NMININ,  NMAXI
  INTEGER LFINAL,   NMINFIN, NMAXF
  INTEGER NAFILEI,  NKFILEI, NZFILEI
  INTEGER NAFILEF,  NKFILEF, NZFILEF

  REAL(8) echeck, rmax, Irmax, Irmax_core

  REAL(8) ena,enk,enz
  INTEGER na,nk,nz

  INTEGER, ALLOCATABLE,DIMENSION(:) :: IPIV
  INTEGER, ALLOCATABLE, DIMENSION(:) :: noi, nof
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: selectwf
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: core

  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: dp, dpst, umx
  REAL(8), ALLOCATABLE, DIMENSION(:) :: eni, enf
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: iamxi, kmxi
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: iamxf, kmxf
  REAL(8), ALLOCATABLE, TARGET, DIMENSION(:) :: zetai, phii
  REAL(8), ALLOCATABLE, TARGET, DIMENSION(:) :: zetaf, phif

  COMPLEX(8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: sxi, sxf
  COMPLEX(8) uim
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: dpout
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: cikmxi, cikmxf
  COMPLEX(8), ALLOCATABLE, DIMENSION(:) :: zwork
  TYPE(szeta) phasei, phasef
  CHARACTER(LEN=100) ARGV

!!!...................................


  EnAU = 27.211396181D+00
  PI   = 2.0D+00 * ASIN(1.0D+00)


!!!...................................

!!! NH2EI = 1, SINGLE IONIZATION CONTINUUM
!!! NH2EF = 2, SINGLE IONIZATION CONTINUUM

!!! NH2EI > 2, DOUBLE IONIZATION CONTINUUM
!!! NH2EF > 2, DOUBLE IONIZATION CONTINUUM


      NNDMX2E = 2
      NHDMX2E = 3
      NCFG    = 4
      NDIPOLE = 8 

      NAFILEI  = 10
      NKFILEI  = 20 
      NZFILEI  = 30 

      NAFILEF  = 11
      NKFILEF  = 21
      NZFILEF  = 31

      NOUTFILE = 16

!!!..........................

      CALL GETARG(1, ARGV)
      READ(ARGV,*) LINITIAL         
      CALL GETARG(2, ARGV)
      READ(ARGV,*) LFINAL       
      CALL GETARG(3, ARGV)
      READ(ARGV,*) GAUGE


!      READ(UNIT=9,*) RMAX

!!!..........................

      MODE = 1
      
      IF(gauge.EQ.'V')      MODE = 0     

      WRITE(*,*) '#'
      WRITE(*,*) '# ndmx2eff::                       dmx in ', gauge, ' gauge'
      WRITE(*,*) '# ndmx2eff::      initial symmetry   li = ', LINITIAL
      WRITE(*,*) '# ndmx2eff::      final   symmetry   lf = ', LFINAL
      WRITE(*,*) '#'

!!!..........................

      OPEN(UNIT=NOUTFILE,FILE='OUT/DPNORMFF.LOG')

!!! READ CFG FILES FOR INITIAL-FINAL SYMMETRIES AND GET VALUES FOR
!!!
!!!    configi, configf, ncsi, ncsf, ndmax,id


!!!.......... READ CONFIGURATION FILE FOR INITIAL FREE STATE


!!! NCSI # NUMBER OF CONFIGURATION SERIES FOR INITIAL STATE
!!! NCSF # NUMBER OF CONFIGURATION SERIES FOR FINAL   STATE
!!!         CALL config_space(ncsi, ncsf)




!!!.....................................

 3    FORMAT(6i5)
      
      CALL HCFGFILE(NCFG, LINITIAL)
      
      READ(ncfg, 3) nSymmetries
      READ(ncfg, 3) l, ls
      READ(ncfg, 3) ncs
      
      CLOSE(NCFG) 

      ncsi = ncs

!!!... final 

      CALL HCFGFILE(NCFG, LFINAL)
      
      READ(ncfg, 3) nSymmetries
      READ(ncfg, 3) l, ls
      READ(ncfg, 3) ncs
      
      CLOSE(NCFG) 

      ncsf = ncs

!!!...................................


      CALL config_space(ncsi, ncsf)


      loi = LINITIAL


      CALL cxfin(NCFG, nhfi,lhfi, lli, nmini, nmxi, noli,&
           &        isi, loi, lsi, ncsmxi, ntoti, ndi, idcsi )


!!!.... CALCULATE TOTAL NUMBER OF CONFIGURATIONS (FREE+BOUND CHANNELS)

         nconfigi = 0
         DO i = 1, ncsmxi
            DO j = 1, ndi(i)   

               nconfigi = nconfigi + 1

            END DO
         END DO


         WRITE(*, *) '# nhdmx2eff::                  ncsi = ', nconfigi
         WRITE(*, *) '#'

!!!         FLI = dfloat( 2*loi + 1 )

!!!.......... READ CONFIGURATION FILE FOR FINAL FREE STATE

         


         lof = LFINAL

         CALL cxfin(NCFG, nhff,lhff, llf, nminf, nmxf, nolf,&
              & isf, lof, lsf, ncsmxf, ntotf, ndf,idcsf )


         nconfigf = 0
         DO i = 1, ncsmxf
            DO j = 1, ndf(i)   
               
            nconfigf = nconfigf + 1
            
         END DO
      END DO


      WRITE(*, *) '# nhdmx2eff::                     ncsf = ', nconfigf
      WRITE(*, *) '#'


!!!      FLF = dfloat( 2 * lof + 1 )

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

      WRITE(*,*) '#'
      WRITE(*,*) '# nhdmx2eff::      initial spin     si = ',  LSI
      WRITE(*,*) '# nhdmx2eff::      final   spin     sf = ',  LSI
      WRITE(*,*) '#'


!!!..........................


      ncsi = ncsmxi
      ncsf = ncsmxf

      
      IF(ncsf.GT.ncsi) THEN

         ndmax = ncsf
         
      ELSE IF(ncsi.GT.ncsf) THEN

         ndmax = ncsi

      ELSE

         ndmax = ncsi
         
      ENDIF


!!!...................................
      WRITE(*,*) '# ndmx2eff::    max nof channels  ndmax = ', ndmax 


      ALLOCATE( phii(  ndmax )  )                    !....... initial states
      ALLOCATE( zetai( ndmax ) )
      ALLOCATE( iamxi( ndmax, ndmax ) )
      ALLOCATE( kmxi(  ndmax, ndmax )  )
      ALLOCATE( sxi(   ndmax,  ndmax )   )
      ALLOCATE( phasei%zeta( ndmax ) )
      ALLOCATE( phasei%sx( ndmax, ndmax) )
      ALLOCATE( phasei%phi( ndmax )      )
      
      ALLOCATE( phif(  ndmax ) )                    !....... final states
      ALLOCATE( zetaf( ndmax ) )
      ALLOCATE( iamxf( ndmax, ndmax) )
      ALLOCATE( kmxf(  ndmax, ndmax) )
      ALLOCATE( sxf(   ndmax, ndmax) )
      ALLOCATE( phasef%phi( ndmax )  )
      ALLOCATE( phasef%zeta( ndmax ) )
      ALLOCATE( phasef%sx( ndmax, ndmax ) )
      
      ALLOCATE( umx( ndmax, ndmax ) )               !........
      ALLOCATE( cikmxi( ndmax, ndmax ) )
      ALLOCATE( cikmxf( ndmax, ndmax ) )
      ALLOCATE( dp( ndmax, ndmax ) )
      ALLOCATE( dpst( ndmax, ndmax ) )
      ALLOCATE( dpout( ndmax, ndmax ) )
      ALLOCATE( selectwf( ndmax, ndmax ) )
      ALLOCATE( core( ndmax, ndmax ) )
      
      
      phasef%phi  = 0.0D+00 
      phasef%zeta = 0.0D+00  
      phasef%sx   = 0.0D+00  
      
      phasei%phi  = 0.0D+00
      phasei%zeta = 0.0D+00
      phasei%sx   = 0.0D+00
      
      
!!!......................................

      CALL  FFHDMX2EFILE(NHDMX2E, LINITIAL, LFINAL, MODE) 
      CALL  FFDMX2EFILE(NNDMX2E, LINITIAL, LFINAL, MODE) 

      OPEN(NDIPOLE, file="DAT/DIPOLE.DAT")


!!!..............

         READ(NHDMX2E) MODE
         READ( UNIT = NHDMX2E) NEI,NEF


         
!         nei = nmaxi
!         nef = nmaxf

!         READ(NHDMX2E) nei, nef

         WRITE(*,*) '# ndmx2eff::                    mode = ', mode
         WRITE(*,*) '# ndmx2eff::                     nei = ', nei
         WRITE(*,*) '# ndmx2eff::                     nef = ', nef
         WRITE(*,*) '#'

         
!!!..............

         WRITE(NNDMX2E) MODE
         WRITE(NNDMX2E) nei, nef

         WRITE(NDIPOLE,*) MODE
         WRITE(NDIPOLE,*) nei, nef

!!!.........................................



         ALLOCATE( eni( nei ) )      ! initial states
         ALLOCATE( noi( nei ) )
         ALLOCATE( enf( nef ) )      ! final states
         ALLOCATE( nof( nef ) )


!!!...............................................


         ! set units for complex number and real matrix

         uim = (0.0D+00, 1.0D+00)
         umx =  0.0D+00

         DO i = 1, ndmax

            umx(i,i) = 1.0D+00
         END DO

!!!.............................................

!            ic : state index
!        en(ic) : energy for 'ic' state
!       noi(ic) : nof channels for 'ic' state
!
!       state  | i > = | ic, eni(ic), noi(ic) >



         CALL AMXFILE(NAFILEI, LINITIAL)
         CALL KMXFILE(NKFILEI, LINITIAL)
         CALL ZMXFILE(NZFILEI, LINITIAL)


!!!.............................................
         WRITE(*,*)'#'
         WRITE(*,*)'# ndmx2eff:: index,      energy,    nof channels    i = '


         DO ic = 1, nei ! initial states loop (ic,nei)


            READ(NHDMX2E) eni(ic), noi(ic)

!!!.........
         
            WRITE(NNDMX2E) ic, eni(ic), noi(ic)
!         WRITE(NDIPOLE, '(i5,2x,e14.6,2x,i5)')  ic, eni(ic), noi(ic)
            WRITE(*,'(14x,I4,1x,G20.8,1x,i4)') ic,eni(ic),noi(ic)



! check consistency with A,K,Z stored matrices

         READ(NAFILEI, '(1e15.7,2x,i5)') ena, na
         READ(NKFILEI, '(1e15.7,2x,i5)') enk, nk
         READ(NZFILEI, '(1e15.7,2x,i5)') enz, nz 

!!!.......      A-MATRIX        
         IF(na.NE.noi(ic).OR.ABS(ena-eni(ic)).GT.1e-8) THEN 

            WRITE(*,*) '# ndmx2eff::        INCOSISTENT A-MATRIX FILE:'
            WRITE(*,*) '# ndmx2eff::        AMX    FILE  ena = ', ena
            WRITE(*,*) '# ndmx2eff::        AMX    FILE   na = ', na
            WRITE(*,*) '# ndmx2eff::        EXIT'
            STOP
         ENDIF
!K-matrix
         IF(nk.NE.noi(ic).OR.ABS(enk-eni(ic)).GT.1e-8) THEN 
            
            WRITE(*,*) '# ndmx2eff::        INCOSISTENT K-MATRIX FILE:'
            WRITE(*,*) '# ndmx2eff::        KMX    FILE enk = ', enk
            WRITE(*,*) '# ndmx2eff::        KMX    FILE  nk = ', nk
            WRITE(*,*) '# ndmx2eff::        EXIT'
            STOP

         ENDIF
!!!.......     Z - MATRIX
         IF(nz.NE.noi(ic).OR.ABS(enz-eni(ic)).GT.1e-8) THEN
            WRITE(*,*) '# ndmx2eff::        INCOSISTENT Z-MATRIX FILE:'
            WRITE(*,*) '# ndmx2eff::        ZMX    FILE E_I = ', enz
            WRITE(*,*) '# ndmx2eff::        ZMX    FILE E_I = ', nz
            WRITE(*,*) '# ndmx2eff::        EXIT'
            STOP
         ENDIF

!!!.......

         iamxi = 0.0D+00  
         kmxi  = 0.0D+00
         zetai = 0.0D+00

         DO j = 1, noi(ic)

            READ(NAFILEI,'(5e15.7)') ( iamxi(m,j), m = 1, noi(ic) )
            READ(NKFILEI,'(5e15.7)') (  kmxi(m,j),  m = 1, noi(ic) )

         END DO


         READ(NZFILEI,'(5e15.7)') ( zetai(j), j = 1, noi(ic) )
         READ(NZFILEI,'(5e15.7)') ( phii(j), j = 1, noi(ic) )


         phasei%zeta => zetai
         phasei%phi => phii


!!!.....................................................


!            ir : state index
!        en(ir) : energy for 'ir' state
!       noi(ir) : nof channels for 'ir' state
!
!       state  | F > = | ir, eni(ir), nof(ir) >

!!! open files for A,K,Z data of final states
         
         CALL AMXFILE(NAFILEF, LFINAL)
         CALL KMXFILE(NKFILEF, LFINAL)
         CALL ZMXFILE(NZFILEF, LFINAL)

         WRITE(*,*)'#'
         WRITE(*,*)'# ndmx2eff:: index,      energy,    nof channels    f = '

!         DO ir = nminfin,  nmaxf   !  final states loop (ir,nef)
         DO ir = 1,  nef   !  final states loop (ir,nef)

         READ(NHDMX2E)      enf(ir), nof(ir)

         

         WRITE(NNDMX2E)  ir, enf(ir), nof(ir)
!         WRITE(NDIPOLE, *)   ir, enf(ir), nof(ir)

         WRITE(*,'(14x,I4,1x,G20.8,1x,i4)')ir,enf(ir),nof(ir)


!!!.......................................................
!!! OPEN DATA OF A, K, Z - MATRIX FOR THE INITIAL STATE IR


!!!...   A - MATRIX

         READ(NAFILEF, '(1e15.7,2x,i5)') echeck, ncheck
         

         IF(ncheck.NE.nof(ir).OR.ABS(echeck-enf(ir)).GT.1e-8) THEN

            WRITE(*,*) ' INCOSISTENT A-MATRIX FILE', echeck, ncheck
            WRITE(*,*) ' HDMX2E FILE E_F = ', enf(ir)
            WRITE(*,*) ' AMX    FILE E_F = ', echeck
            WRITE(*,*) ' HDMX2E FILE NDF = ', nof(ir)
            WRITE(*,*) ' AMX    FILE E_F = ', ncheck
            
            STOP
            
         ENDIF

!!!...   K - MATRIX


         READ(NKFILEF, '(1e15.7,2x,i5)') echeck, ncheck
         
            IF(ncheck.NE.nof(ir).OR.ABS(echeck-enf(ir)).GT.1e-8) THEN

               WRITE(*,*)' INCOSISTENT A-MATRIX FILE', echeck, ncheck
               WRITE(*,*)' HDMX2E FILE E_F = ', enf(ir)
               WRITE(*,*)' KMX    FILE E_F = ', echeck
               WRITE(*,*)' HDMX2E FILE NDF = ', nof(ir)
               WRITE(*,*)' KMX    FILE E_F = ', ncheck
               
               STOP
            
            ENDIF

!!!...   Z - MATRIX

         READ(NZFILEF, '(1e15.7,2x,i5)') echeck, ncheck

         
         IF(ncheck.NE.nof(ir).OR.ABS(echeck-enf(ir)).GT.1e-8) THEN

            WRITE(*,*) ' INCOSISTENT A-MATRIX FILE', echeck, ncheck
            WRITE(*,*) ' HDMX2E FILE E_F = ', enf(ir)
            WRITE(*,*) ' KMX    FILE E_F = ', echeck
            WRITE(*,*) ' HDMX2E FILE NDF = ', nof(ir)
            WRITE(*,*) ' KMX    FILE E_F = ', ncheck
         
            STOP
            
         ENDIF

!!!.......................................


         iamxf = 0.0D+00
         kmxf  = 0.0D+00
         zetaf = 0.0D+00

         DO j = 1, nof(ir)

            READ(NAFILEF,'(5e15.7)' ) ( iamxf(m,j), m = 1, nof(ir) )
            READ(NKFILEF,'(5e15.7)' ) ( kmxf(m,j),  m = 1, nof(ir) )
         END DO

!!!.........

         READ(NZFILEF,'(5e15.7)') ( zetaf(j), j = 1, nof(ir) )
         READ(NZFILEF,'(5e15.7)') ( phif(j), j = 1, nof(ir) )

         phasef%zeta => zetaf
         phasef%phi => phif

!!!.............................................................

!!! START RE-NORMALIZING 2-E DIPOLES 

         dp = 0.0D+00


         READ(NHDMX2E) ( (dp(j,i), j = 1, nof(ir)), i = 1, noi(ic) )



!!!...............
!!! calculate the dipole for standing wave b.c.
!!!        R_{sw} = A_INITIAL * R * A^(T)_FINAL 


            dp  = transpose(dp)
            dpst = matmul(iamxi, dp)
            dpst = matmul(dpst, transpose(iamxf)) 


!!!            write(9,'(4e14.6)') ((dpst(j,i), j=1,noi(ic)), i=1,nof(ir))


!!!!................
!!!! change to the S-matrix normal for ini. state


            cikmxi = umx - uim * kmxi


!!! imsl
!!!         call dlincg(noi(ic), cikmx, ndmax, cikmxinv, ndmax)
!!!         sxi = matmul(cikmxinv, (umx + uim*kmxi))
!!! imsl


!!       complex matrix inversion      s1x ---> s1x^{-1}

            ALLOCATE( IPIV(ndmax))
            ALLOCATE( ZWORK(64*ndmax) )

!!! Make the LU factorization

            CALL zgetrf(noi(ic), noi(ic), cikmxi, ndmax, ipiv, info) 

         IF(INFO == 0) THEN

            CALL zgetri(noi(ic), cikmxi, ndmax, ipiv, zwork, 64*ndmax,info) 

            ELSE

               WRITE(*,*) ' MATRIX IS SINGULAR.  K^(-1) MATRIX DOES NOT EXIST'
               STOP

            ENDIF



         DEALLOCATE( ZWORK )
         DEALLOCATE( IPIV )

!!!..........
!!!.....................   final states

!!       complex matrix inversion      s1x ---> s1x^{-1}


         cikmxf = umx + uim * kmxf


            ALLOCATE( IPIV(ndmax))
            ALLOCATE( ZWORK(64*ndmax) )

!!! Make the LU factorization

            CALL zgetrf(nof(ir), nof(ir), cikmxf, ndmax, ipiv, info) 

         IF(INFO == 0) THEN

            CALL zgetri(nof(ir), cikmxf, ndmax, ipiv, zwork, 64*ndmax,info) 

         ELSE

               WRITE(*,*) ' MATRIX IS SINGULAR.  K^(-1) MATRIX DOES NOT EXIST'
               STOP

            ENDIF



         DEALLOCATE( ZWORK )
         DEALLOCATE( IPIV )

!!!..........

!!! imsl
!!!            call dlincg(nof(ir), cikmx, ndmax, cikmxinv, ndmax)
!!!            sxf = matmul((umx - uim*kmxf),cikmxinv)
!!! imsl

!!!...........................................................


         sxi = MATMUL(cikmxi, (umx + uim * kmxi))
         sxf = MATMUL( (umx - uim * kmxf), cikmxf)


!!! .... S - matrix normalization

         dpout = MATMUL(cikmxi, dpst)
         dpout = MATMUL(dpout, TRANSPOSE(cikmxf))


!            write(9,'(4e14.6)') ((real(dpout(j,i))**2+imag(dpout(j,i))**2,&
!            & j=1,noi(ic)), i=1,nof(ir))


!!!.................................................................
!!!         
!!!
!!!         ADD   CORRECTIONS   TERMS  TO FINAL DIPOLE 
!!!
!!!.................................................................


         phasei%sx => sxi
         phasef%sx => sxf


         !  note: j = initial states, i = final states

         WRITE(NNDMX2E) ( ( dpout(j,i), j = 1, noi(ic)), i = 1, nof(ir) )
!         WRITE(NDIPOLE,'(3e15.7)')   ( ( dpout(j,i), j = 1, noi(ic)), i = 1, nof(ir) )


         WRITE(NDIPOLE,'(4e15.7)') enf(ir)-eni(ic),&
              & ( ( dpout(j,i), j = 1, noi(ic)), i = 1, nof(ir) ),&
              & ((REAL(dpout(j,i))**2+imag(dpout(j,i))**2,j=1,noi(ic)), i=1,nof(ir))

!         WRITE(NDIPOLE, *)   ir, enf(ir), nof(ir)
!         WRITE(NDIPOLE, '(i5,2x,e14.6,2x,i5)')  ic, eni(ic), noi(ic)

      END DO

      WRITE(*,*)'########################################################'

      CLOSE(NAFILEF) 
      CLOSE(NKFILEF) 
      CLOSE(NZFILEF)


   END DO


   CLOSE(NAFILEI) 
   CLOSE(NKFILEI) 
   CLOSE(NZFILEI) 

   CLOSE(NHDMX2E)
   CLOSE(NNDMX2E)
   CLOSE(NDIPOLE)
   
!!!............................
   
   DEALLOCATE(dp)
   DEALLOCATE(dpst)
   DEALLOCATE(dpout)
   DEALLOCATE(umx)


   DEALLOCATE(eni)
   DEALLOCATE(noi)
   DEALLOCATE(phii)
   DEALLOCATE(zetai)
   DEALLOCATE(iamxi)
   DEALLOCATE(kmxi)
   DEALLOCATE(sxi)
   DEALLOCATE(cikmxi)   
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
   DEALLOCATE(phif)
   DEALLOCATE(zetaf)
   DEALLOCATE(iamxf)
   DEALLOCATE(kmxf)
   DEALLOCATE(sxf)
   DEALLOCATE(cikmxf)
   DEALLOCATE(nhff)
   DEALLOCATE(lhff)
   DEALLOCATE(llf)
   DEALLOCATE(nminf)
   DEALLOCATE(nmxf)
   DEALLOCATE(ndf)
   DEALLOCATE(nolf)
   DEALLOCATE(isf)
   DEALLOCATE(idcsf)



DEALLOCATE(selectwf)
DEALLOCATE(core)

    END PROGRAM dpnormff
!!!#####################################################
!!!EOF

