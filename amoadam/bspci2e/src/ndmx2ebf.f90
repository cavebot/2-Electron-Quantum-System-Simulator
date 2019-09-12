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
!!!              nd : open channels ata energy 'energy'
!!!            ndmax: maximum number of open channels   
!!!
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
PROGRAM ndmx2ebf

  USE PRECISION, ONLY: dpk
  USE units,     ONLY: enau, alpha
  USE ioroutines

  IMPLICIT NONE
  CHARACTER(LEN=16) GAUGE
  INTEGER i, ii, j, jj, m, ne, nd, ndmax, ncheck, info,ic
  INTEGER l_ini, s_ini, no_ini_states, l_fin
  INTEGER alloc_error, dealloc_error
  INTEGER, ALLOCATABLE,DIMENSION(:) :: IPIV
  INTEGER MODE, LINITIAL, NMININ, NMAXI, LFINAL, NMINFIN, NMAXF, K
  INTEGER NOUT, NAFILE, NKFILE, NHDMX2E, NNDMX2E 
  INTEGER ndipole
  REAL(dpk) energy, energy_ini, echeck, total, factor
  INTEGER, ALLOCATABLE,DIMENSION(:)    :: n1,l1,l2 
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: dp, dpst
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: iamx, kmx , umx
  COMPLEX(dpk) uim
  COMPLEX(dpk), ALLOCATABLE, DIMENSION(:) :: dpout
  COMPLEX(dpk), ALLOCATABLE, DIMENSION(:,:) :: cikmx,  sx
  COMPLEX(dpk), ALLOCATABLE, DIMENSION(:) :: zwork
  CHARACTER(LEN=100) ARGV

!............

  CALL GETARG(1, ARGV)
  READ(ARGV,*) LINITIAL         
  CALL GETARG(2, ARGV)
  READ(ARGV,*) LFINAL       
  CALL GETARG(3, ARGV)
  READ(ARGV,*) GAUGE

!

  uim = (0.0D+00, 1.0D+00)

  NNDMX2E    = 2
  NHDMX2E    = 3
  NDIPOLE    = 8 
  NAFILE     = 10
  NKFILE     = 10 + NAFILE
  NOUT       = 16

!!!..........................


  OPEN(UNIT=NOUT,FILE='out/ndmx2ebf.out')

  WRITE(*,*) '# ndmx2ebf::  normalization of dmx2e gauge = ',gauge
  
  
  MODE = 1
  IF(gauge.EQ.'V')          MODE = 0     
  
  l_ini = linitial
  s_ini = 1 
  no_ini_states = NMAXI - NMININ + 1
      
      !
  CALL dmxfile(nhdmx2e, "dat", "dmx2ebf"  ,"bin"  , gauge, linitial, lfinal)
  CALL dmxfile(nndmx2e, "dat", "ndmx2ebf" ,"bin"  , gauge, linitial, lfinal)
  CALL dmxfile(ndipole, "dat", "dipole"   ,"ascii", gauge, linitial, lfinal)
  

  READ( UNIT = NHDMX2E) MODE
  READ( UNIT = NHDMX2E) NMININ, NMAXI
  READ( UNIT = NHDMX2E) NMINFIN, NMAXF
  !
  WRITE(UNIT = NNDMX2E) MODE
  WRITE(UNIT = NNDMX2E) ABS(NMAXI - NMININ + 1), NMAXF
  !
  WRITE(NDIPOLE, *) MODE
  WRITE(NDIPOLE, *) ABS(NMAXI - NMININ + 1), NMAXF

!!........


      initial_states: DO IC = NMININ, NMAXI


         READ( UNIT=NHDMX2E )  energy_ini

         WRITE(NDIPOLE,'(i5,2x,e15.7)')  ic, energy_ini

         WRITE(UNIT=NNDMX2E)       ic, energy_ini
         WRITE(*,*)  '# ndmx2ebf                     e_in   =', energy_ini, ic 


         CALL hfile(nkfile,"dat","k2e","ascii",lfinal)
         CALL hfile(nafile,"dat","a2e","ascii",lfinal)

!


         final_states: DO i = 1, nmaxf


            READ(NHDMX2E) energy, nd

!            WRITE(*,*)  '# ndmx2ebf    e_f, nd(e_f),i = ', energy, nd,i 

            ALLOCATE(       n1( nd )     )
            ALLOCATE(       l1( nd )     )
            ALLOCATE(       l2( nd )     )
            ALLOCATE(       dp( nd )     )
            ALLOCATE(     dpst( nd )     )
            ALLOCATE(    dpout( nd )     )
            ALLOCATE(     iamx( nd, nd ) ) !real
            ALLOCATE(      kmx( nd, nd ) )
            ALLOCATE(      umx( nd, nd ) )
            ALLOCATE(       sx( nd, nd ) ) ! complex
            ALLOCATE(    cikmx( nd, nd ) )


! final energy e_i:                   < g | d | E_i, j >

            READ(NHDMX2E) ( dp(j), j = 1, nd)      !d_if(j) = < i | d | wf(j = 1, nd) >  
         
            umx = 0.0D+00 
            DO ii = 1, nd
               umx(ii,ii) = 1.0_dpk
            END DO


!!! open  A,K-matrix data files for the  d_if(j) normalization


         READ(NAFILE, '(1e15.7,2x,i5)') echeck, ncheck

         ! check for inconsistent files

         IF(ncheck.NE.nd.OR.ABS(echeck - energy).GT.1.0D-03)  THEN            
            WRITE(*,*) ' INCOSISTENT A-MATRIX FILE', echeck, ncheck
            WRITE(*,*) ' HDMX2E FILE        E_F = ', energy
            WRITE(*,*) ' AMX    FILE        E_F = ', echeck
            WRITE(*,*) ' HDMX2E FILE        ND  = ', nd
            WRITE(*,*) ' AMX    FILE        E_F = ', ncheck
            STOP
         ENDIF


         READ(NKFILE, '(1e15.7,2x,i5)') echeck, ncheck
         !
         IF(ncheck.NE.nd.OR.ABS(echeck - energy).GT.1.0D-03)      THEN 
            WRITE(*,*) ' INCOSISTENT K-MATRIX FILE', echeck, ncheck
            WRITE(*,*) ' HDMX2E FILE        E_F = ', energy
            WRITE(*,*) ' AMX    FILE        E_F = ', echeck
            WRITE(*,*) ' HDMX2E FILE        ND  = ', nd
            WRITE(*,*) ' AMX    FILE       E_F = ', ncheck
            STOP
         ENDIF
         

         DO j = 1, nd
            READ(NAFILE, '(4i5)') jj, n1(j),l1(j),l2(j)              ! (n,l,lp)
         ENDDO
         
         DO j = 1, nd
            READ(NAFILE, '(5e15.7)') ( iamx(m, j), m = 1, nd)         ! A^-1 matrix elements
            READ(NKFILE, '(5e15.7)') (  kmx(m, j), m = 1, nd)         ! pi * K matrix elements
         END DO

         ! transform box-normalized dipoles to outgoing waves normalized 

         dpout = 0.0_dpk
         dpst = 0.0_dpk


          dpst = MATMUL(iamx, dp)     ! standing-wave b.c.: R_{sw} = A^{-1} * R_box 
          cikmx = umx + uim * kmx     !                          C = 1 + i *pi *K

         !  invert C :  C --> C^{-1}


         ALLOCATE(ipiv(nd))
         ALLOCATE( zwork(64*nd) )
         CALL zgetrf(nd, nd, cikmx, nd, ipiv, info)  !lu factorization
         IF(INFO == 0 ) THEN
            CALL zgetri(nd, cikmx, nd, ipiv, zwork, 64*nd, info) 
         ELSE
            WRITE(*,*) ' matrix is singular. K^(-1) does not exist'
            STOP
         ENDIF
         DEALLOCATE( ZWORK )
         DEALLOCATE( IPIV  )


         sx = MATMUL((umx - uim*kmx ), cikmx)       ! S-matrix

                                       ! outgoing normalization
         dpout = MATMUL(cikmx, dpst)   ! R^{-} = C^{-1}* R_{sw} 
                                       !       = C^{-1}* A^{-1} * R_box

         ! save data

      
         WRITE(ndipole,'(i5,e15.7,i5)')  i, energy, nd 
         DO j = 1, nd
            WRITE(ndipole, '(4i5)') j, n1(j),l1(j),l2(j)              ! (n,l,lp)
         ENDDO
         WRITE(ndipole,'(86e15.7)')   (dpout(j), j = 1, nd)
      
         
         WRITE(NNDMX2E)   i, energy, nd
         WRITE(NNDMX2E) ( n1(j),    j = 1, nd)
         WRITE(NNDMX2E) ( l1(j),    j = 1, nd)
         WRITE(NNDMX2E) ( l2(j),    j = 1, nd)
         WRITE(NNDMX2E) ( dpout(j), j = 1, nd)

      ! deallocate


      DEALLOCATE(n1,l1,l2)
      DEALLOCATE(dp)
      DEALLOCATE(dpst)
      DEALLOCATE(dpout)
      DEALLOCATE(iamx)
      DEALLOCATE(kmx, STAT=dealloc_error)
      IF(dealloc_error /= 0) THEN 
         PRINT*," DEALLOCATION ERROR 1 " 
         STOP
      ENDIF
      DEALLOCATE(umx)
      DEALLOCATE(sx)
      DEALLOCATE(cikmx)
      

   ENDDO final_states
   
   
   CLOSE(NAFILE)
   CLOSE(NKFILE)
      
ENDDO initial_states

      CLOSE(NHDMX2E)
      CLOSE(NNDMX2E)
      CLOSE(NDIPOLE)

    END PROGRAM ndmx2ebf

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



