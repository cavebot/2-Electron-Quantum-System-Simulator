!!! program to calculate the final dipole moment
!!!
!!!
!!!   l.aa.n/iesl/2001
!!! 
!!!         07/SEP/2001:
!!!                     
!!!                  *   no input file as 'dpnorm.inp'
!!!                  *   modified to compiled by pgf90/pgf77 compilers 
!!!                  *   data taken by the outpout of kmtx.f90 program
!!!                  *   information is entered trhough the CFG-L.INP
!!!                      files
!!!  
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
PROGRAM dpnorm

  USE ioroutines

  IMPLICIT NONE
  CHARACTER(LEN=16) GAUGE
  INTEGER EnAU
  INTEGER i, ii, j, m, ne, nd, ndmax, ncheck, info,ic
  INTEGER l_ini, s_ini, no_ini_states, l_fin
  INTEGER alloc_error, dealloc_error
  INTEGER, ALLOCATABLE,DIMENSION(:) :: IPIV
  INTEGER MODE, LINITIAL, NMININ, NMAXI, LFINAL, NMINFIN, NMAXF, K
  INTEGER NOUT, NCFG, NAFILE, NKFILE, NHDMX2E, NNDMX2E 
  INTEGER NCS1PH, NDIPOLE
  REAL(8) energy, energy_ini, echeck, total, factor
  REAL(8) de_au, cs1ph_o, alpha, pi
  REAL(8), ALLOCATABLE, DIMENSION(:) :: dp, dpst, dp_part, cs1
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: iamx, kmx , umx
  COMPLEX(8) uim
  COMPLEX(8), ALLOCATABLE, DIMENSION(:) :: dpout
  COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: cikmx,  sx
  COMPLEX(8), ALLOCATABLE, DIMENSION(:) :: zwork

!!!...................................


  EnAU = 27.211396181D+00

!!!...................................

!!! NH2E = 1, SINGLE IONIZATION CONTINUUM
!!! NH2E = 2, DOUBLE IONIZATION CONTINUUM


      NNDMX2E = 2
      NHDMX2E = 3
      NCFG    = 4
      NDIPOLE = 8 
      NCS1PH  = 9
      NAFILE  = 10
      NKFILE  = 10 + NAFILE
      NOUT    = 16

!!!..........................

      OPEN(UNIT=9,file='INP/HDMX2E.INP',status='old')

      READ(UNIT=9,*) LINITIAL, NMININ, NMAXI
      READ(UNIT=9,*) LFINAL,   NMINFIN, NMAXF
      READ(UNIT=9,*) GAUGE 


      CLOSE(9)



      OPEN(UNIT=NOUT,FILE='OUT/HDMX2E.LOG')

!!!..........................

      IF(gauge.EQ.'V') THEN

         WRITE(*,*) ' NORMALIZE 2e-DME IN VELOCITY GAUGE '

         MODE = 0     
      ELSE
            
         WRITE(*,*) ' NORMALIZE 2e-DME IN LENGTH GAUGE '

         MODE = 1
      ENDIF
!!!..........................
!!!        ndmax         : maximum number of open channels   
!!!        energy_ini    : energy of initial bound state
!!!             l_ini    : angular momentum of initial (bound) state
!!!             l_fin    : angular momentum of the final (continuum) state
!!!
  
!!!      CALL FindNoChannels(NCFG, LFINAL, ndmax)
!!!      WRITE(*,*) ' MAXIMUM NUMBER OF CHANNELS NDMAX = ', ndmax

!!$      ALLOCATE(dp_part(ndmax))
!!$      ALLOCATE(dp(ndmax))
!!$      ALLOCATE(dpst(ndmax))
!!$      ALLOCATE(dpout(ndmax))
!!$      ALLOCATE(cs1(ndmax))
!!$      ALLOCATE(iamx(ndmax,ndmax))
!!$      ALLOCATE(kmx(ndmax,ndmax))
!!$      ALLOCATE(umx(ndmax,ndmax))
!!$      ALLOCATE(sx(ndmax,ndmax))
!!$      ALLOCATE(cikmx(ndmax,ndmax))

!!!.......................................

      l_ini = LINITIAL
      s_ini = 1 
      no_ini_states = NMAXI - NMININ + 1

      
!!! set units for complex number and real matrix
!!!.............................................
      
      uim = (0.0D+00, 1.0D+00)
      
!!$      umx = 0.0D+00
!!$      DO i = 1, ndmax
!!$
!!$         umx(i,i) = 1.0D+00
!!$
!!$      END DO

!!!............................................

      CALL  HDMX2EFILE(NHDMX2E, LINITIAL, LFINAL, MODE) 
      CALL  NDMX2EFILE(NNDMX2E, LINITIAL, LFINAL, MODE) 
      CALL  CS1PH2EFILE(NCS1PH, LINITIAL, LFINAL, MODE) 
      
!!      OPEN(UNIT = NCS1PH,FILE ='DAT/CS1PH.DAT')

      OPEN(UNIT = NDIPOLE,FILE ='DAT/DIPOLE.DAT')


      READ( UNIT = NHDMX2E) MODE

      WRITE(UNIT = NNDMX2E) MODE
      WRITE(UNIT = NDIPOLE,*) MODE


      DO IC = NMININ, NMAXI


!!!.................................  Energy of initial state


         READ( UNIT=NHDMX2E)  energy_ini

         WRITE(NDIPOLE,'(e15.7)')  energy_ini*0.5D+00
         WRITE(UNIT=NNDMX2E)       energy_ini*0.5D+00


         WRITE(*,*)  'INITIAL STATE ENERGY     E_in   =', energy_ini, ic 

      
!!!.................................


      DO i = 1, NMAXF



!!! READ FINAL ENERGY AND   D_IF(J) = < I | D | F(J = 1, ND) >  

!!!........................................ 


         READ(NHDMX2E) energy, nd

!!!........................................

         ALLOCATE(dp_part(nd))
         ALLOCATE(dp(nd))
         ALLOCATE(dpst(nd))
         ALLOCATE(dpout(nd))
         ALLOCATE(cs1(nd))

!!!.......

!!real
         ALLOCATE(iamx(nd,nd))
         ALLOCATE(kmx(nd,nd))
         ALLOCATE(umx(nd,nd))
!! complex
         ALLOCATE(sx(nd,nd))
         ALLOCATE(cikmx(nd,nd))

!!!........................................ 



         READ(NHDMX2E) (dp(j), j = 1, nd)


!!!.......................................
         
         WRITE(*,*)  'FINAL STATE ENERGY     E_f   =', energy 
         WRITE(*,*)  'NUMBER OF CHANNELS  ND(E_F)  =', nd




         umx = 0.0D+00
         DO ii = 1, nd

            umx(ii,ii) = 1.0D+00

         END DO


!!!.........................................


!!! OPEN DATA OF A,K-MATRIX FOR NORMALIZING THE D_IF(J) DME

         ! check consistency in files

!!!.................
         


         CALL AMXFILE(NAFILE, LFINAL)


         READ(NAFILE, '(1e15.7,2x,i5)') echeck, ncheck


         IF(ncheck.NE.nd.OR.ABS(echeck - energy).GT.1.0D-08)      THEN
            
            WRITE(*,*) ,' INCOSISTENT A-MATRIX FILE', energy, ncheck
            WRITE(*,*) ,' HDMX2E FILE E_F = ', energy
            WRITE(*,*) ,' AMX    FILE E_F = ', echeck
            WRITE(*,*) ,' HDMX2E FILE ND  = ', nd
            WRITE(*,*) ,' AMX    FILE E_F = ', ncheck

            STOP

         ENDIF


!!!........................



         CALL KMXFILE(NKFILE, LFINAL)

         READ(NKFILE, '(1e15.7,2x,i5)') energy, nd

         IF(ncheck.NE.nd.OR.ABS(echeck - energy).GT.1.0D-08)      THEN
            
            WRITE(*,*) ,' INCOSISTENT K-MATRIX FILE', energy, ncheck
            WRITE(*,*) ,' HDMX2E FILE E_F = ', energy
            WRITE(*,*) ,' AMX    FILE E_F = ', echeck
            WRITE(*,*) ,' HDMX2E FILE ND  = ', nd
            WRITE(*,*) ,' AMX    FILE E_F = ', ncheck
            STOP

         ENDIF
         


!!!........................

                  
         DO j = 1, nd

            READ(NAFILE, '(5e15.7)') ( iamx(m, j), m = 1, nd)         ! A^-1 matrix elements
            READ(NKFILE, '(5e15.7)') ( kmx(m, j),  m = 1, nd)         ! pi * K matrix elements
            
         END DO
         

         WRITE(*,*) 'READ  A^-1 and  pi * K'
         
         !!! calculate the dipole for standing wave (sw) b.c. 
         !!!
         !!!        R_{sw} = A^{-1} * R 
         
!!!............................. 

         DO m = 1, nd
            
            dpst(m) = 0.0D+00
            
            DO j = 1, nd
               
               dpst(m) = dpst(m) + iamx(m, j) * dp(j)
               
            END DO
         END DO




!!!............................

         
         !  dpst = matmul(iamx, dp)
         !         write(9,'(5f13.8)') energy, 
         !         & (dpst(j), j = 1, nd)


         !  calculate  (1 + i *pi *K)
!!!............................



         
         cikmx = umx + uim * kmx



!!!....................................


         !  invert   (1 + i *pi *K)


!!!imsl
!!!         call dlincg(nd, cikmx, ndmax, cikmxinv, ndmax)
!!!         sx = matmul((umx - uim*kmx), cikmxinv)         
!!!         calculate the dipole for outgoing wave b.c.
!!!         R^{-} = (1 - i *pi *K)^{-1} * R_{sw}
!!!         cikmx --> cikmx^-1
!!!imsl


!!
!!       complex matrix inversion      s1x ---> s1x^{-1}
!!

         ALLOCATE(ipiv(nd))
         ALLOCATE( zwork(64*nd) )
 

!!! Make the LU factorization

         CALL zgetrf(nd, nd, cikmx, nd, ipiv, info) 



         IF(INFO == 0 ) THEN

            CALL zgetri(nd, cikmx, nd, ipiv, zwork, 64*nd, info) 


         ELSE

            WRITE(*,*) ' MATRIX IS SINGULAR.  K^(-1) MATRIX DOES NOT EXIST'
            STOP

         ENDIF

!!!   WRITE(*,*) ' S * S^(-1)  ',  MATMUL( ( ux + uim * kx ), s1x)



         DEALLOCATE( ZWORK )
         DEALLOCATE(IPIV)


!!!         CALL zgetri(nd, cikmx, ndmax, ipiv, cikmxinv, ndmax,info) 


         sx = MATMUL((umx - uim*kmx ), cikmx)

!!!........................................


!!!         WRITE(*,*) 'K^-1 Matrix ::'
!!!         WRITE(*,'(10f11.5)') ((cikmx(m,j),m=1,nd),j=1,nd) 
!!!         dpout = matmul(cikmxinv, dpst)


         do m = 1, nd

            dpout(m) = (0.0D+00, 0.0D+00)
            do j = 1, nd

               dpout(m) = dpout(m) + cikmx(m,j) * dpst(j)

            end do
         end do

!!!...........................................

!         

!         nd : number of open channels at energy E
!!
!!   prepare proportional factors/ results in  SI system
!!


      PI      = 3.141592654D+00
      alpha   = 1.0D+00/137.0359D+00


      de_au  = ABS( energy - energy_ini )/2


!!      cs1ph_o = 28.002 * ( 2 * pi * pi) * alpha


      IF(mode.EQ.0) THEN

         cs1ph_o = 5.337819D+00 * 3 / (2*de_au)

      ELSE

         cs1ph_o = 5.337819D+00 * (3.0D+00/2.0D+00) * de_au

      ENDIF


         DO j = 1, nd

            
            cs1(j) = cs1ph_o * ABS(dpout(j))**2

!!!    cs1(j) = cs1ph_o * ( 2.0D0/3.0D0) * float(2*LFINAL + 1) * ABS(dpout(j))**2/de_au

         ENDDO


         WRITE(NCS1PH,'(6e15.7)')   EnAU * ( energy - energy_ini ) / 2.0D+00, ( cs1(j), j = 1, nd)


         WRITE(NDIPOLE,'(e15.7,i5)')  energy/2.0D+00, nd 
         WRITE(NDIPOLE,'(6e15.7)')   (dpout(j), j = 1, nd)


         WRITE(NNDMX2E) 0.5 * energy,nd
         WRITE(NNDMX2E) (ABS(dpout(j)),j = 1,nd)


!!!.................................


         DEALLOCATE(kmx, STAT=dealloc_error)

         IF(dealloc_error /= 0) THEN 

            PRINT*," DEALLOCATION ERROR 1 " 
            STOP
         ENDIF


         DEALLOCATE(iamx,STAT=dealloc_error)

         IF(dealloc_error /= 0) THEN 
            
            PRINT*," DEALLOCATION ERROR 2 " 
            STOP
         ENDIF


         DEALLOCATE(cikmx, STAT=dealloc_error)

         IF(dealloc_error /= 0) THEN 

            PRINT*," DEALLOCATION ERROR 3 " 
            STOP
         ENDIF


         DEALLOCATE(sx, STAT=dealloc_error)

         IF(dealloc_error /= 0) THEN 

            PRINT*," DEALLOCATION ERROR 4 " 
            STOP
         ENDIF
         
         DEALLOCATE(umx, STAT=dealloc_error)

         IF(dealloc_error /= 0) THEN 

            PRINT*," DEALLOCATION ERROR 5 " 
            STOP
         ENDIF

      ENDDO

         DEALLOCATE(dp_part)
         DEALLOCATE(dp)
         DEALLOCATE(dpst)
         DEALLOCATE(dpout)
         DEALLOCATE(cs1)


   ENDDO

      CLOSE(NHDMX2E)
      CLOSE(NNDMX2E)
      CLOSE(NAFILE)
      CLOSE(NKFILE)
      CLOSE(NCS1PH)
      CLOSE(NDIPOLE)


!!!!!  Output
!#######################################################################
! Original statement
!         write(9,'(9e13.5)') (energy + 1.10635525)/2.0,(dpout(j), j=1,nd)
!
!
! Statement for output of the energy in a.u. with the same zero of energy 
!  as in the fixed-boundary codes.
!
!         write(9,'(9e15.7)') energy/2.0d0,(dpout(j), j=1,nd)
!
! New   output - total Dipole Matrix Element
!         total=abs(dpout(1))**4 + abs(dpout(2))**4     !FALSCH 
!         write(9,'(7e13.5)') (energy+1.10635525)/2.0,total
!
!
!      Very new output: 
!      
!       photon energy in eV  vs.  
!       one-photon crosss section in Mb :
! 
!sigma_{1-photon} (in Mb) = 8.006728 * (2/3)*(2*L_{final}+1) * (|D|^2/DE) 

! D is the dipole transition matrix element in velocity gauge and in a.u. 
!(Hartree) (for a continuum(!) final state normalised per unit energy), 
! 
!               DE = energy difference in Ry
!                       
!     and thus:
! sigma_{1-photon} (in Mb) = 5.337819 * (2*L_{final}+1) * (|D|^2/DE)  
!D in velocity gauge and in a.u. (Hartree) , DE = energy difference in Ry
!
!                1Mb = 10^{-18} cm^2 
!                a_o^2 = 28.00286 MB             ,    a_o  Bohr radius
!    sigma_{1-photon} (in a.u) = sigma_{1-photon} (in MB) / 28.00286
!#######################################################################
!
!       write(9,'(25e13.5)') (energy - energy_ini) * 13.6055,         &
!       ( ( 5.337819D+00 * float( 2*l_fin + 1) * abs(dpout(j))**2     &
!       / (energy - energy_ini) ), j=1,nd )
!!!        factor = 5.337819D+00 * float(2*l_fin+1) / (energy - energy_ini)
!
! 
!       do j = 1, nd
!            dp_part(j) = factor * abs(dpout(j))**2 
!       
!         end do
!
!       if (nd .eq. 1) then
!          
!          write(99,'(10e13.5)') (energy - energy_ini)*13.6055, dp_part(1)
!
!      endif
!      if (nd .EQ. 4) then                       
!
!   write(99,'(10e13.5)') (energy - energy_ini)*13.6055, dp_part(1),     &
! ( ( dp_part (2) + dp_part(3) + dp_part(4) ) / dp_part(1) ) * 100.D+0, &
! ( ( dp_part(3) + dp_part(4) ) / dp_part(2) )                  
!         
!      endif
!      if (nd .EQ. 9) then
!  write(99,'(10e13.5)') (energy - energy_ini)*13.6055, dp_part(1),  &
! ( ( dp_part (2) + dp_part(3) + dp_part(4) ) / dp_part(1) ) * 100.D+0,&
!( ( dp_part(5) + dp_part(6) + dp_part(7) + dp_part(8) + dp_part(9) ) & 
!   / dp_part(1) ) * 100.D+0,                            &
!   ( ( dp_part(3) + dp_part(4) ) / dp_part(2) ),        &  
!   ( ( dp_part(6) + dp_part(7) ) / dp_part(5) ),        &
!   ( ( dp_part(8) + dp_part(9) ) / dp_part(5) )
!   endif                        
!

    END PROGRAM dpnorm

