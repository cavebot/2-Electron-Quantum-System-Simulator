!/media/usbdisk_1/home/nlambros/andriana_01_10_2003/ATOM/ATOMLASER/tdse/mch/run/he 
!                                 written by          Xian  Tang,        Jan 1990
!                                 changes made by     Henri Rudholf,     Feb 1990
!                                 revised to f90 by   Jian  Zhang,       Feb 1996
!                                 changes made by    L.A.A. Nikolopoulos Jun 2000
!                                                 
!   changes made 3/31/90 by hr and xtang:
!   changed to remove the rapid oscillating part of the matrix
!   pulse shape changed to gaussian
!   pulse duration changed from -0.5* tau to 0.5*tau
!   pulse supported sin^2
!   pulse duration satrts from tau= 0. to --
!
!#   this program uses the NAG routine d02baf to integrate the diff. equation
!#   and it can be used for restarting an old job with a unit 33, with time and
!#   vector information.

!
!   this program uses the IMSL routine divprk to integrate the diff. equation
!   and it can be used for restarting an old job with a unit 33, with time and
!   vector information.
!   
!   
!
!   changes made by J. Zhang to multichannel TDSE,   feb  1996   
!   changes made by L. A.A. Nikolopoulos             June 2000
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      PROGRAM read_dmx2e_mch

        use precision
        USE units
        USE PARAMETER
        USE deriv
        USE io

        IMPLICIT NONE
!
        INTEGER, PARAMETER:: nv=10000
        REAL(dpk), DIMENSION(nv)::  y, yp

!

        REAL(dpk) de1, de2, thresh_1, thresh_2, thresh_3, ecut
!pulse
        REAL(dpk) ff, topt, tfac, t, t1, tim, tol
!
        REAL(dpk) pop_g, pop, sion, dion
        REAL(dpk) p_0, p_1n, p_n
        REAL(dpk) p_n_single, p_n_double
        REAL(dpk) ati1, ati2, ati3, ion1, ion2, ion3, slqp
!
        INTEGER ifail, count_n
        INTEGER irestart, irespec, nl, ido, nemax, index
        INTEGER i, j, p, q, k, k1, ni, nf, nb, lb,lz, iu, lstart, kt,lk
        INTEGER nout1, nout2, nt, li, l, ikt, ik, ii
!
        INTEGER, ALLOCATABLE, DIMENSION(:,:) ::noch
        INTEGER, ALLOCATABLE, DIMENSION(:) ::ne
!
        INTEGER, DIMENSION(8) :: n2, n3, nsgl
        INTEGER, DIMENSION(8) :: nbs
!
        REAL(dpk), ALLOCATABLE, DIMENSION(:,:) ::energy
        REAL(dpk), ALLOCATABLE, DIMENSION(:) :: dpr, dpi, time, dipole
        REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: w         !nag
!
        CHARACTER*25 coef_t,ion_t,pop_t, dipole_t, infile
        LOGICAL exists, hohg
        EXTERNAL fcn, d02baf

        !      interface
        !         real function ddot(nd,x,nx,y,ny)
        !            integer nd,nx,ny
        !            real(8), dimension(nd) :: x,y
        !         end function ddot
        !      end interface
        !
        !#######################################################################
!      write(*,'(a30, $)')' name of output file = '
!      call link("unit18=ntest4.in//")
!

!    read from input file (infile)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!        CALL getarg(1, argv)              !get arguments
!        READ(argv,*)   ff                 ! in SI 
!        CALL getarg(2, argv)              !get arguments
!        READ(argv,*)   omeg               ! in SI
!        CALL getarg(3, argv)              !get arguments
!        READ(argv,*)   tfac
        CALL getarg(1, argv)              !get arguments
        READ(argv,*)   nl

!        ff = ff/i0_td             
!        omeg = omeg/enau

!        READ(*,*) infile
        OPEN(18,file="tinp/tdse_bs_mch.inp")

        READ(18,'(a20)')   coef_t       ! c_i(t)    : coefficients  
        READ(18,'(a20)')   ion_t        ! Y(t)      : Y(+), Y(++)
        READ(18,'(a20)')   pop_t        ! pop(t)    : populations
        READ(18,'(a20)')   dipole_t     ! d(t)      : dipole
        READ(18,*)         irestart     !
        READ(18,*)         irespec      !
        READ(18,*)         nmax         ! Max no of States   in L = 0,1,2,. nl-1
        READ(18,*)         nemax        ! Max no of Energies in L = 0,1,2,. nl-1
        CLOSE(18)


!        READ(18,*)         ff           ! Io(a.u.)  : peak Intensity
!        READ(18,*)         omeg         ! w(a.u.)   : photon energy
!        READ(18,*)         tfac         ! # of cycles
!        READ(18,*)         nl           ! Total no of angul. included

!        hohg = .false.

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


!        OPEN(10,file=ion_t)

!        fi   = SQRT(ff)
!        topt = 2.0d+00 * m_pi / omeg
!        tau  = tfac * topt
        
!        WRITE(10,*)'# tdse_bs_mch::     peak intensity     Io = ',ff   , ' a.u.'
!        WRITE(10,*)'# tdse_bs_mch::  photon frequency     Wph = ',omeg , ' a.u.'
!        WRITE(10,*)'# tdse_bs_mch::     pulse duration     Tp = ',tau  , ' a.u.'


! some initial checks


!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxx  Memory allocations :  y, yp, dpr, dpi, dipole, ne, energy, noch
!
!            y(nmax * nl * 2) :                wf coefficients
!           yp(nmax * nl * 2) :                derivative of y
!            time (8192)        :
!            dipole(time) = dpr(nmax) + I * dpi(nmax) = < y(t) | d | y(t) >
!            ne :
!       energy  :  number of energies
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      CALL make_space(nmax,nl)



!!!      allocate(y(nmax*nl*2),yp(nmax*nl*2))


      ALLOCATE(            dpr(nmax) )
      ALLOCATE(            dpi(nmax) )
      ALLOCATE(             ne(nl) )
      ALLOCATE( energy( nemax, nl) )
      ALLOCATE(   noch( nemax, nl) )



! 17   format(2x,1p8e15.7)



!         Read data from dipole input files.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     lb         :  angular momentum of initial state
!     lz         :  angular momentum of final   state
!     n(i)       :  #  channel states in symmetry  L = i
!     nbs(i)     :  # of bound states in symmetry  L = i
!     ne(i)      :  # of energy states in symmetry L = i
!  energy(i,j)   :  energy of i-th state of symmetry L = j      
!                :      i = 1, ne , j = 1, nl + 1.     
!     noch(i,j)  :  # of channels at i-th energy state of symmetry L = j 
!     dzr(i,j,k) :  real part of dme from channel state i of symmetry L = k 
!                                    to channel state j of symmetry L = k + 1
!    dzi(i,j,k)  : imaginary part of dme from channel state i of symmetry L = k 
!                                    to channel state j of symmetry L = k + 1
! 
!               < i, k | d | j, k+1 >  
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      if(nl.le.1) stop

      do i = 1, nl-1

         CALL openfl(17,i-1,i)
                                ! Note: should check consistency lb+1 = i, 
                                !lz+1=i+1
         READ(17,*) lb,lz,n(i),n(i+1), nbs(i), nbs(i+1), ne(i), ne(i+1)
 
         IF(n(i).GT.nmax)   WRITE(*,*)  '# tdse_bs_mch::     ni > nmax,   ni = ', n(i)
         IF(n(i+1).GT.nmax) WRITE(*,*)  '# tdse_bs_mch:: n(i+1) > nmax, ni+1 = ', n(i+1)

         READ(17,*) ( energy(k,lb+1), k =          1,   ne(i) )        ! in Ryd
         READ(17,*) ( energy(k,lz+1), k =          1, ne(i+1) )        ! in Ryd
         READ(17,*) ( noch(k,i),      k = nbs(i) + 1,   ne(i) ) 
         READ(17,*) ( noch(k,i+1),    k = nbs(i+1)+1, ne(i+1) )

         DO  ni = 1, n(i)
            READ(17,*) ( dzr( nf, ni, i ), nf = 1, n(i+1) )
            READ(17,*) ( dzi( nf, ni, i ), nf = 1, n(i+1) )
         END DO

!         write(*,*) lb, lz, n(i), n(i+1), nbs(i), nbs(i+1), ne(i), ne(i+1)
         WRITE(*,*) "# tdse_bs_mch::                 lb = ", lb
         WRITE(*,*) "# tdse_bs_mch::                 lz = ", lz
         WRITE(*,*) "# tdse_bs_mch::               n(i) = ", n(i)
         WRITE(*,*) "# tdse_bs_mch::             n(i+1) = ", n(i+1)
         WRITE(*,*) "# tdse_bs_mch::             nbs(i) = ", nbs(i)
         WRITE(*,*) "# tdse_bs_mch::           nbs(i+1) = ", nbs(i+1)
         WRITE(*,*) "# tdse_bs_mch::              ne(i) = ", ne(i)
         WRITE(*,*) "# tdse_bs_mch::            ne(i+1) = ", ne(i+1)


         close(17)

      end do

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Expand the energy array

!!!   en(i,j) :: energy of channel state i= 1,2, ... index of symmetry j.
!!!                Must be        index == n(i) 

      DO i = 1, nl
         
         DO j = 1, nbs(i)

            noch(j,i) = 1
            
         END DO

         index = 0
         do k = 1, ne(i)
            do j = 1, noch(k,i)

               index = index + 1                   
               en(index,i) = energy(k,i)       ! in Ryd

            end do
         end do
      end do

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!     Basis is discretized. So dipole matrix elements have to be normalized
!     further. The algorithm is: 
!       (Cowan, 18-8, pg 535, Pseudo-discrete treatment of continuum problems)
!  
!     D(i,j) ---> D(i,j) * sqrt(de)
!                                        
!                 de = e_j - e_i
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      de1 = 0.008               ! Energy separation below E++ = 0. (in  a.u. )
      de2 = 0.0175              ! Energy separation above E++ = 0. (in  a.u. )


      normalize_dmx: DO i = 1, nl-1

!!!        BC transitions   d(n,m)   n : bound state      n of symmetry i 
!!!                                  m : continuum state  m of symmetry i+1

         normalize_dmx_bc: DO ni = 1, nbs(i)

            DO nf = nbs(i+1)+1, n(i+1)

               IF( en(nf,i+1).LE.0.0) THEN

                  dzr(nf,ni,i) = dzr(nf,ni,i) * SQRT(de1)
                  dzi(nf,ni,i) = dzi(nf,ni,i) * SQRT(de1)
                  
               ELSE

                  dzr(nf,ni,i) = dzr(nf,ni,i) * SQRT(de2)
                  dzi(nf,ni,i) = dzi(nf,ni,i) * SQRT(de2)
                  dzr(nf,1,i)  = dzr(nf,1,i)  * SQRT(0.0)
                  dzi(nf,1,i)  = dzi(nf,1,i)  * SQRT(0.0)

               END IF

            END DO
            
         END DO normalize_dmx_bc
         
!!!        CB transitions   d(n,m)   n : continuum state  n of symmetry i 
!!!                                  m : bound     state  m of symmetry i+1

         normalize_dmx_cb: DO ni = nbs(i) + 1, n(i)
            DO nf = 1, nbs(i+1)


               IF(en(ni,i).LE.0.0) THEN

                  dzr(nf,ni,i) = dzr(nf,ni,i) * SQRT(de1)
                  dzi(nf,ni,i) = dzi(nf,ni,i) * SQRT(de1)
                  
               ELSE

                  dzr(nf,ni,i) = dzr(nf,ni,i) * SQRT(de2)
                  dzi(nf,ni,i) = dzi(nf,ni,i) * SQRT(de2)
                  
               END IF

            END DO
         END DO normalize_dmx_cb

!!!        CC transitions   d(n,m)   n : continuum state  n of symmetry i 
!!!                                  m : continuum state  m of symmetry i+1

         normalize_dmx_cc: DO ni = nbs(i)+1, n(i)

            DO nf = nbs(i+1)+1, n(i+1)

               IF(en(ni,i).GT.0.0.AND.en(nf,i+1).GT.0.0) THEN

                  dzr(nf,ni,i) = dzr(nf,ni,i) * de2 
                  dzi(nf,ni,i) = dzi(nf,ni,i) * de2

               ELSE
                  
                  IF(en(ni,i).LE.0.0.AND.en(nf,i+1).LE.0.0) THEN

                     dzr(nf,ni,i) = dzr(nf,ni,i) * de1
                     dzi(nf,ni,i) = dzi(nf,ni,i) * de1

                  ELSE

                  dzr(nf,ni,i) = dzr(nf,ni,i) * sqrt(de1 * de2)
                  dzi(nf,ni,i) = dzi(nf,ni,i) * sqrt(de1 * de2)

               END IF
            END IF
            END DO
         END DO normalize_dmx_cc

      ENDDO normalize_dmx

!!     end of re-normalization
!!!xxxx


! write the dipole moments
      write_dmx: DO i= 2, 2

            index = 0

            DO k = 1, ne(i)

!               WRITE(*,*) energy(k,i)

               DO j = 1, noch(k,i)
                  
                  index = index + 1

                  IF(j.EQ.1) WRITE(*,'(I4,1X,4e14.5)') index &
                   &, dzr(2,index,i)**2+dzi(2,index,i)**2, dzr(2,index,i), dzi(2,index,i)


!         if(j.eq.2) write(*,'(2e13.5)') dzr(index,88,i-1),dzi(index,88,i-1)


               END DO
            END DO

         END DO write_dmx
         


       END PROGRAM read_dmx2e_mch
    
!!!.....................................................................


      subroutine openfl(n,l1,l2)
       IF(l1.EQ.0.AND.l2.EQ.1)   OPEN(n,file='dat/os01.dat')
       IF(l1.EQ.1.AND.l2.EQ.2)   OPEN(n,file='dat/os12.dat')
       IF(l1.EQ.2.AND.l2.EQ.3)   OPEN(n,file='dat/os23.dat')
       IF(l1.EQ.3.AND.l2.EQ.4)   OPEN(n,file='dat/os34.dat')
       IF(l1.EQ.4.AND.l2.EQ.5)   OPEN(n,file='dat/os45.dat')
       IF(l1.EQ.5.AND.l2.EQ.6)   OPEN(n,file='dat/os56.dat')
       IF(l1.EQ.6.AND.l2.EQ.7)   OPEN(n,file='dat/os67.dat')
       IF(l1.EQ.7.AND.l2.EQ.8)   OPEN(n,file='dat/os78.dat')
       IF(l1.EQ.8.AND.l2.EQ.9)   OPEN(n,file='dat/os89.dat')
       IF(l1.EQ.9.AND.l2.EQ.10)  OPEN(n,file='dat/os910.dat')
       IF(l1.EQ.10.AND.l2.EQ.11) OPEN(n,file='dat/os1011.dat')
       IF(l1.EQ.11.AND.l2.EQ.12) OPEN(n,file='dat/os1112.dat')
       IF(l1.EQ.12.AND.l2.EQ.13) OPEN(n,file='dat/os1213.dat')
       IF(l1.EQ.13.AND.l2.EQ.14) OPEN(n,file='dat/os1314.dat')
       IF(l1.EQ.14.AND.l2.EQ.15) OPEN(n,file='dat/os1415.dat')
       IF(l1.EQ.15.AND.l2.EQ.16) OPEN(n,file='dat/os1516.dat')
       IF(l1.EQ.16.AND.l2.EQ.17) OPEN(n,file='dat/os1617.dat')
       IF(l1.EQ.17.AND.l2.EQ.18) OPEN(n,file='dat/os1718.dat')
       IF(l1.EQ.18.AND.l2.EQ.19) OPEN(n,file='dat/os1819.dat')
      end subroutine openfl


!!!eof
