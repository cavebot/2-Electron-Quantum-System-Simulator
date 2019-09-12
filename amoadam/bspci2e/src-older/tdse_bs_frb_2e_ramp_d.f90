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


      PROGRAM tdse_bs_mch

        use precision
        USE units
        USE PARAMETER
        USE deriv
        USE io

        IMPLICIT NONE
!
        INTEGER,       PARAMETER::  nv = 7000
        REAL(dpk), DIMENSION(nv)::  y, yp

!

        REAL(dpk) de1, de2, thresh_1, thresh_2, thresh_3, ecut
!pulse
        REAL(dpk) ff, topt, tfac, t, t1, tim, tol, fi
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

        CALL getarg(1, argv)              !get arguments
        READ(argv,*)   ff                 ! in SI 
        CALL getarg(2, argv)              !get arguments
        READ(argv,*)   omeg               ! in SI
        CALL getarg(3, argv)              !get arguments
        READ(argv,*)   tfac
        CALL getarg(4, argv)              !get arguments
        READ(argv,*)   ltot

        ff = ff/i0_td             
        omeg = omeg/enau

!        READ(*,*) infile
        OPEN(18,file="tinp/tdse_bs_mch.inp")

        READ(18,'(a20)')   coef_t       ! c_i(t)    : coefficients  
        READ(18,'(a20)')   ion_t        ! Y(t)      : Y(+), Y(++)
        READ(18,'(a20)')   pop_t        ! pop(t)    : populations
        READ(18,'(a20)')   dipole_t     ! d(t)      : dipole
        READ(18,*)         irestart     !
        READ(18,*)         irespec      !
        READ(18,*)         nmax         ! Max no of States   in L = 0,1,2,. ltot-1
        READ(18,*)         nemax        ! Max no of Energies in L = 0,1,2,. ltot-1  
        CLOSE(18)


!        READ(18,*)         ff           ! Io(a.u.)  : peak Intensity
!        READ(18,*)         omeg         ! w(a.u.)   : photon energy
!        READ(18,*)         tfac         ! # of cycles
!        READ(18,*)         ltot         ! Total no of angul. included

        hohg = .false.

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


        OPEN(10,file=ion_t)

        fi   = SQRT(ff)
        topt = 2.0d+00 * m_pi / omeg
        tau  = tfac * topt
        
        WRITE(10,*)'# tdse_bs_mch::     peak intensity     Io = ',ff   , ' a.u.'
        WRITE(10,*)'# tdse_bs_mch::  photon frequency     Wph = ',omeg , ' a.u.'
        WRITE(10,*)'# tdse_bs_mch::     pulse duration     Tp = ',tau  , ' a.u.'


! some initial checks

!xxx  Is it a restarted run or not?

      if(irestart.eq.1.and.irespec.ne.1) then
         write(10,*) ' *** this is a restarted run *** '
      endif

!xxx   Restarted run: 
!xxx   check to see if the output file exists

      inquire(file=coef_t,exist=exists)

      if((irestart.eq.1).and..not.exists) then
 200     write(10,*) 'Error: input claims this to be a restart, but your'
         write(10,*) 'unit33 does not exist, the program will now abort'
         stop
      end if

!xxx   First run:
!xxx   check to make sure that a output file does not exist when not a restart

      if((irestart.ne.1).and.exists) then
 201     write(10,*) 'Error: input claims this not to be a restart, but '
         write(10,*) 'your unit33 is not empty, the program will now abort'
         stop
      end if

      open(33,file=coef_t,form='unformatted',access='sequential')
      open(34,file=pop_t)

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxx  Memory allocations :  y, yp, dpr, dpi, dipole, ne, energy, noch
!
!            y(nmax * ltot * 2) :                wf coefficients
!           yp(nmax * ltot * 2) :                derivative of y
!            time (8192)        :
!            dipole(time) = dpr(nmax) + I * dpi(nmax) = < y(t) | d | y(t) >
!            ne :
!       energy  :  number of energies
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      CALL make_space(nmax,ltot)

!!!      allocate(y(nmax*ltot*2),yp(nmax*ltot*2))

      ALLOCATE(  time(8192) )
      ALLOCATE(dipole(8192) )

      ALLOCATE(            dpr(nmax) )
      ALLOCATE(            dpi(nmax) )
      ALLOCATE(             ne(ltot) )
      ALLOCATE( energy( nemax, ltot) )
      ALLOCATE(   noch( nemax, ltot) )

!      ALLOCATE(ne(ltot), energy(nemax,ltot), noch(nemax,ltot))

! 17   format(2x,1p8e15.7)



!         Read data from dipole input files.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     lb         :  angular momentum of initial state
!     lz         :  angular momentum of final   state
!     n(i)       :  #  channel states in symmetry  L = i
!     nbs(i)     :  # of bound states in symmetry  L = i
!     ne(i)      :  # of energy states in symmetry L = i
!  energy(i,j)   :  energy of i-th state of symmetry L = j      
!                :      i = 1, ne , j = 1, ltot + 1.     
!     noch(i,j)  :  # of channels at i-th energy state of symmetry L = j 
!     dzr(i,j,k) :  real part of dme from channel state i of symmetry L = k 
!                                    to channel state j of symmetry L = k + 1
!    dzi(i,j,k)  : imaginary part of dme from channel state i of symmetry L = k 
!                                    to channel state j of symmetry L = k + 1
! 
!               < i, k | d | j, k+1 >  
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      nl = ltot

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

                  dzr(nf,ni,i) = dzr(nf,ni,i) * sqrt(de2)
                  dzi(nf,ni,i) = dzi(nf,ni,i) * sqrt(de2)
                  dzr(nf,1,i)  = dzr(nf,1,i)  * sqrt(0.0)
                  dzi(nf,1,i)  = dzi(nf,1,i)  * sqrt(0.0)

               END IF

            END DO
            
         END DO normalize_dmx_bc
         
!!!        CB transitions   d(n,m)   n : continuum state  n of symmetry i 
!!!                                  m : bound     state  m of symmetry i+1

         normalize_dmx_cb: DO ni = nbs(i) + 1, n(i)
            DO nf = 1, nbs(i+1)


               IF(en(ni,i).LE.0.0) THEN

                  dzr(nf,ni,i) = dzr(nf,ni,i) * sqrt(de1)
                  dzi(nf,ni,i) = dzi(nf,ni,i) * sqrt(de1)

               ELSE

                  dzr(nf,ni,i) = dzr(nf,ni,i) * sqrt(de2)
                  dzi(nf,ni,i) = dzi(nf,ni,i) * sqrt(de2)

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
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

       if(irestart.eq.1) then
          read(33) ff, omeg, tau, ltot, ecut
       else
          write(33) ff, omeg, tau, ltot, ecut
       endif

!      ecut=30.0d0
!      ecut=-thresh_1


!!!     Transform en to photoelectron energies         en --> k^2/2 
!!!     and ignore states with photoelectron energy    en > ecut.              
!!!      Energies in the calculation :                 en < ecut

!!! he+ thresholds

      thresh_1  = -4.0D+00             ! He+(1s)                 (in Ryd)
      thresh_2  = -1.0D+00             ! He+(2s,2p)
      thresh_3  = -0.4444444444D+00    ! He+(3s,3p,3d) 
       ecut     = 5.0d0               ! in a.u. 


!       ntot = 0
!       shift_energy_axis: DO  k = 1, nl   !relative to first ionization threshold
!          iu = 0         
!          DO  k1 = 1, n(k)
!             en(k1,k) = (en(k1,k) - thresh_1) / 2.d0
!             if( en(k1,k).gt.ecut) exit
!             iu = iu + 1
!          ENDDO



       ntot = 0
       shift_energy_axis: DO  k = 1, nl   !relative to first ionization threshold

          iu = 0
         

          DO 22 k1 = 1, n(k)

             en(k1,k) = (en(k1,k) - thresh_1) / 2.d0

             if(en(k1,k).gt.ecut) go to 22
             iu = iu + 1
22           CONTINUE



! Now number of channel states in the calculation is new one according to:
!
!       en < ecut
!

          n(k) = iu

!
          IF(irestart.EQ.1) THEN
             READ(33) n(k), ( en(k1,k), k1 = 1, n(k))
          ELSE
             WRITE(33) n(k),( en(k1,k), k1 = 1, n(k))
          ENDIF


!!! Total number of channel states included in the calculation

          ntot = ntot + n(k)
          
          WRITE(10,*)'# tdse_bs_mch::     angular momenta   l = ', k - 1
          WRITE(10,*)'# tdse_bs_mch::     nof channels   n(l) = ', n(k)

       ENDDO shift_energy_axis

       IF(irestart.EQ.1) THEN
          READ(33) ntot
       ELSE
          WRITE(33) ntot
       ENDIF


       WRITE(10,*) '# tdse_bs_mch::   nof channels propagated   n_tot =  ', ntot
       WRITE(*, *) '# tdse_bs_mch::   nof channels propagated   n_tot =  ', ntot
!--------------------------------------------------------------------


!!!  Make 2-d energy array en(i,j) 1-d array dag(nb) 
!!!
!!!     en(k1,k) ----> dag(nb)       nb = k * n(k) + k1

       count_n = 0
       nb = 0
       DO k = 1, nl
          DO k1 = 1, n(k)

             nb      = nb + 1
             dag(nb) = en(k1,k)

             count_n = count_n + 1
          ENDDO
       ENDDO

       WRITE(*,*)  '# tdse_bs_mch::                e_g = ', en(1,1) 
       WRITE(*,*)  '# tdse_bs_mch ::           dagm(1) = ', dag(1)
       WRITE(*,*)  '# tdse_bs_mch ::        nof_states = ', count_n
       WRITE(*,*)  '#  t(au),             p_g,             norm'


!!! if it is re-Run

       IF(irestart.EQ.1) THEN

!.0.5tau          t=-0.5d0*tau

!!!xxxxxx         read yp(t1) and assign to y(t1)
          
215       READ(33, END=210) t1, (yp(ii),ii = 1, 2*ntot)

          t = t1

          DO  ik = 1, 2 * ntot

             y(ik) = yp(ik)

          ENDDO

          GOTO 215
          
 210      lstart = int( 0.5 * (t + 0.5*tau) / topt ) + 1
                                !!!  lstart != 1 it is restart
                                !!!.0.5tau          t=t+0.5*tau
       ELSE

          lstart = 1

       ENDIF

       !xxx      Initial conditions at t = 0. y(1) : System in its ground state

       IF(lstart.EQ.1) THEN

          y    = 0.0d0
          y(1) = 1.d0
          t    = 0.0d0

       ENDIF


       IF(irespec.EQ.1) GOTO 249

!!
       IF(hohg) THEN 
          nsum(1) = 1
          DO i = 2, ltot
             nsum(i) = nsum(i-1) + n(i-1)
             
          END DO
          WRITE(*,*) '# tdse_bs_mch::  ntot = nsum(ltot) ' , ntot, ' = ', nsum(ltot) 
       ENDIF
    !!

       kt = 0  

       IF(hohg) THEN
          OPEN(41,file=dipole_t)
       ENDIF 

       ALLOCATE(yderiv(2*ntot))
       ALLOCATE( w( nmax * ltot * 2, 7))         !!! work array for nag

!!!  Propagate     from lstart ----> 8*tfac  

       propagate_wf: DO  l = lstart, INT(tfac*8)


!             WRITE(*,*) "# tdse_bs_mch::          lstart = ",lstart

          nout1 = kt + 1

          propagate_wf_inner: DO  lk = 1, 4 !break time-step for dipole(t)

             kt = kt + 1 
             tim = t + topt * 0.03125

!             WRITE(*,*) "# tdse_bs_mch::          lk = ",lk


!nag t-> t + dt

             neval = 0
               tol = 1.0d-8
              ifail = 0 

             CALL d02baf(t, tim, 2*ntot, y, tol, fcn, w, ifail)

             time(kt) = tim

             harmonics: IF(hohg) THEN

                dipole(kt) = 0.0D+00
                calculate_dipole_t: DO li = 1, nl-1              ! dipole(t)

                                !     < f | p | df/dt >

                   CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,&
                        &yderiv(nsum(li+1)),1,0.0d0,dpr,1)
                   CALL dgemv('t',n(li+1),n(li),-1.0d0,dzi(1,1,li),nmax,&
                        &yderiv(ntot+nsum(li+1)),1,1.0d0,dpr,1)
                   CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,&
                        &yderiv(ntot+nsum(li+1)),1,0.0d0,dpi,1)
                   CALL dgemv('t',n(li+1),n(li),1.0d0,dzi(1,1,li),nmax,&
                        &yderiv(nsum(li+1)),1,1.0d0,dpi,1)

                   dipole(kt)=dipole(kt)&
                        &+DOT_PRODUCT(y(nsum(li):nsum(li)+n(li)-1),dpr)&
                        &+DOT_PRODUCT(y(ntot+nsum(li):ntot+nsum(li)+n(li)-1),dpi)

                                !     < df/dt | p |f >

                   CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,&
                        &y(nsum(li+1)),1,0.0d0,dpr,1)
                   CALL dgemv('t',n(li+1),n(li),-1.0d0,dzi(1,1,li),nmax,&
                        &y(ntot+nsum(li+1)),1,1.0d0,dpr,1)
                   CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,&
                        &y(ntot+nsum(li+1)),1,0.0d0,dpi,1)                   
                   CALL dgemv('t',n(li+1),n(li),1.0d0,dzi(1,1,li),nmax,&
                        &y(nsum(li+1)),1,1.0d0,dpi,1)

                   dipole(kt) = dipole(kt)&
                        &+DOT_PRODUCT(yderiv(nsum(li):nsum(li)+n(li)-1),dpr)&
                        &+DOT_PRODUCT(yderiv(ntot+nsum(li):ntot+nsum(li)+n(li)-1),dpi) 


                   
                !     <f|dp/dt|f>

                !   t1=-0.5d0*tau+tim
                !   slqp=-fi*(2*dcos(t1/tau*pi)*pi/tau * dsin(t1/tau*pi)&
                !     &*dcos(omeg*t1)/omeg&
                !     & + dcos(t1/tau*pi)**2 * dsin(omeg*t1))

                !     dipole(kt)=dipole(kt) + slqp*(&
                !     &dot_product(y(nsum(li):nsum(li)+n(li)-1),dpr)&
                !     &+dot_product(y(ntot+nsum(li):ntot+nsum(li)+n(li)-1),dpi))

                END DO calculate_dipole_t

             ENDIF harmonics


          ENDDO propagate_wf_inner
!99           CONTINUE


!             deallocate(yderiv)
!          write(10,*)'number of fun. evaluations used=', neval





          IF(hohg) THEN       ! harmonics

             WRITE(*,*) '# tdse_bs_mch::  readout dipole(t) '

             nout2 = kt
          
             DO ikt = nout1, nout2
                
                t1 = -0.5d0 * tau + time(ikt)
                
                WRITE(41,'(2x,3e14.6)') time(ikt), dipole(ikt)
                
!0.5tau           write(41,'(2x,3e14.6)') time(ikt)-0.5e0*tau, dipole(ikt)
                
             END DO
          ENDIF                 !harmonics




!save data

          IF(irestart.EQ.1) THEN

             REWIND (33)
             WRITE(33) ff,omeg,tau,ltot,ecut

             DO  k = 1, nl

                WRITE(33) n(k), (en(k1,k),k1=1,n(k))  !save energies
             ENDDO

             WRITE(33) ntot

             WRITE(*,*) 'wrote ntot', ntot

             irestart = 0

          ENDIF


!0.5*tau          write(33) tim-0.5e0*tau,(y(ii),ii=1,2*ntot)
!0.5*tau          write(34,'(e15.6)') tim-0.5*tau


          WRITE(33) tim, ( y(ii), ii = 1, 2 * ntot)  ! save wf



          WRITE(34,'(e15.6)') tim

          pop = 0.0d+00

!             WRITE(*,*) '# tdse_bs_mch::  readout population '

! monitor population  
             nt = 0
             DO  k = 1, nl

!                WRITE(*,*) '# tdse_bs_mch::   l = ',k-1

             IF(k.LE.MIN(10, nl)) THEN

                WRITE(34,'(i3)') k - 1


                bound_states: DO  k1 = 1, 5              ! bound states population for l=0-4

                   WRITE(34,'(3e15.6)') en(k1,k),&
                        & ( y(k1 + nt)**2 + y( ntot + nt + k1)**2 )

                ENDDO bound_states
                

                norm: DO  k1 = 1, n(k)          ! the whole population (norm)

                   pop  = pop + y(k1 + nt)**2 + y( ntot + nt + k1 )**2
                END DO norm

             ENDIF

             nt = nt + n(k)
          ENDDO

!0.5tau          write(*,*)'t=', tim-0.5*tau, 'norm=', pop

          WRITE(*,'(3E15.8)') tim,  y(1)**2 + y( ntot + 1)**2, pop
          write(34,'(/)')


!   nsg(i) :  # of single-ionization channels in symmetry L = i,  i.e. he+

          nsgl(1) = 9    ! i.e. : 1ss, (2pp, 2ss), (3dd, 3pp, 3ss), (7dd, 7pp, 7ss)
          nsgl(2) = 9    ! i.e. : 1sp, (2ps, 2pd, 2sp),(3dp, 3df, 3ps, 3pd, 3sp) 
          nsgl(3) = 13   ! i.e. : 1sd, (2pp, 2pf, 2sd, 2dd), (3dd, 3dg, 3pp, 3pf, 3sd)
                         !                                   (7dd  7pp  7sd)
                   !          write(33) n(k),(en(k1,k),k1=1,n(k))

          pop  = 0.0d+00   !he
          sion = 0.0d+00   !he+
          dion = 0.0d+00   !he++

           p_0 = 0.0D+00   ! e_k < threshold_1
          p_1n = 0.0D+00   ! threshold_1 < e_k < threshold_n (n = infty)
          p_n  = 0.0D+00   ! e_k > threshold_n

          nt = 0
          populations: DO k = 1, nl          ! angular momenta
             DO  k1 = 1, n(k)   ! energies
                

                ! (e_k < 0.0 a.u)  : bound states , (ground + excited He) 
               

                IF(en(k1,k).LE.0.0) THEN 

                   p_0 = p_0 + ( y(k1+nt)**2 + y(ntot+nt+k1)**2 )
                ENDIF

                ! (0.0 < en < 2.0 a.u)   strictly single ionization channels 
                
                IF( ( en(k1,k).GT.0.0d+00).AND.(en(k1,k).LT.2.0d+00) ) THEN

                   p_1n = p_1n + ( y( k1 + nt )**2 + y( ntot + nt + k1 )**2 )
                ENDIF
 
                ! (en > 2.0 a.u)   single and double ionization channels
                
                IF(en(k1,k).GE.2.0) THEN
                   
                   p_n = p_n + ( y(k1+nt)**2 + y(ntot+nt+k1)**2 )
                ENDIF
                
                IF(en(k1,k).GE.2.0.AND.en(k1,k).EQ.en( k1 - nsgl(k), k)) THEN
                   dion = dion + ( y(k1+nt)**2 + y(ntot+nt+k1)**2 )
                ENDIF

                                !                endif
             ENDDO

             nt = nt + n(k)

          ENDDO populations

          pop_g = y(1)**2 + y(ntot+1)**2
          pop   = p_0 
          sion  = p_1n + p_n_single
          dion  = p_n  - p_n_single


!#######################################################################
!    pop + sion + dion = 1,
!    pop + ion = 1,
!    popg + pope + sion + dion = 1
!
!#######################################################################
!          y(k1 + nt + j)   : coefficient of j-th channel of angular 
!                             momentum L dependent on the value of nt 
!                             at energy en(k1, L)   
!                             
!                           k1 = 1, 2, .... n(L)      L = 0,1,2 .. Lmax
!                           j = 0, 1, 2,... nchannels,  at energy en(k1,L) 
!                           nt(i)  +=  n(L) 
!
!          nt = 0,           L = 0,        y(k1 + nt + j) 
!                               population of j-th channel with L = 0          
!                               at energy en en(k1,0)
!#######################################################################

!
!    popg :  Population in ground state       -     popg = y(1)^2 + y(ntot+1)^2 
!    pope :  Population in bound excited states -
!    pop  :  Ground + excited  states       -       pop = popg  + pope
!   sion  :  Single Ionization Channels     - 
!   dion  :  double ionization channels     -
!    ion  :  Ionization channels            -       ion = sion + dion


          WRITE(10,908) tim/topt, pop_g, pop, p_0, p_1n, p_n, dion  
!          WRITE(10,908) tim/topt, pop_g, 1.0 - pop - dion, dion

          pop = 0.0
          t = tim
          nt = 0

       ENDDO propagate_wf

       close(33)
       close(34)


!!!  Skip PES = PES(t)

       goto 999


!!   Monitor PES = PES(t)
!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 249   nt = 0

       OPEN(40,file='tdat/pes.dat')

                                !0.5tau       write(40,*) ' spectrum for time = ', t-0.5e0*tau
       write(40,*) ' spectrum for time = ', t/topt

       do  k = 1, nl

          write(40,*) ' Spectrum for l = ', k - 1

          do  k1 = 1, n(k)
             
             write(40,'(4e14.6)') en(k1,k),&
             &(y(k1+nt)**2+y(ntot+nt+k1)**2),y(k1+nt),y(ntot+nt+k1)
          
          enddo

          nt = nt + n(k)

       enddo

       close(40)
!xxxxxxxxxxxxxxxxxxxxxxxx
 
908    FORMAT( 1x, 7e14.6)
988    FORMAT( 2x, i4, i7, 5e14.6)
900    FORMAT( 2x, i7, 6(e10.3, 2x) )
400    FORMAT( 1x, '# tdse_bs_mch::  pulse duration (a.u) =', e15.8 )
300    FORMAT( 1x, '# tdse_bs_mch::  peak intensity (a.u) =', e15.8 )
       
!xxxxxxxxxxxxxxxxxxxxxxx


999    CONTINUE
       
       IF(hohg)   CLOSE(41)
       
      
!          nsgl(1) = 9     i.e. : 1ss, (2) (2pp, 2ss),      (3) (3dd, 3pp, 3ss), (7)(7ss,7pp,7dd)
!          nsgl(2) = 9     i.e. : 1sp, (2) (2ps, 2pd, 2sp), (3) (3dp, 3df, 3ps, 3pd, 3sp) 
!          nsgl(3) = 13    i.e. : 1sd, (2) (2pp, 2pf, 2sd, ), 
                                !      (3) (3ds, 3dd, 3dg, 3pp, 3pf, 3sd) 
                                !       (7dd  7pp  7sd)


      n2(1) = 2                 ! (2pp, 2ss),                         L = 0
      n3(1) = 3                 ! (3dd, 3pp, 3ss), 

      n2(2) = 3                 ! (2ps, 2pd, 2sp),                    L = 1
      n3(2) = 5                 ! (3dp, 3df, 3ps, 3pd, 3sp) 

      n2(3) = 3                 ! (2pp, 2pf, 2sd, 2dd),               L = 2
      n3(3) = 6                 ! (3dd, 3dg, 3pp, 3pf, 3sd)


!!! transform he+ thresholds to ionization-core thresholds in a.u.

      thresh_2 =  ( thresh_2 - thresh_1) * 0.5D+00
      thresh_3 =  ( thresh_3 - thresh_1) * 0.5D+00 

      WRITE(*, *) '# tdse_bs_mch::          E( He+(2s,2p) )   = ', thresh_2  
      WRITE(*, *) '# tdse_bs_mch::          E( He+(3s,3p,3d)) = ', thresh_3

!#######################################################################
!          y(k1 + nt + j)   : coefficient of j-th channel of angular 
!                             momentum L dependent on the value of nt 
!                             at energy en(k1, L)   
!                             
!                           k1 = 1, 2, .... n(L)      L = 0,1,2 .. Lmax
!                           j = 0, 1, 2,... nchannels,  at energy en(k1,L) 
!                           nt(i)  +=  n(L) 
!
!          nt = 0,           L = 0,        y(k1 + nt + j) 
!                               population of j-th channel with L = 0          
!                               at energy en en(k1,0)
!#######################################################################
 
       OPEN(40,file='tdat/ati.dat')
       OPEN(44,file='tdat/yield_n.dat')

       do  k = 1, nl

          ion1 = 0.0d0
          ion2 = 0.0d0
          ion3 = 0.0d0

          do  k1 = 1, n(k)             

                                ! Increase by energy step (k1 !=  k1 -1 ) 
                                ! and don't take bound states

             if(en(k1,k).ne.en(k1-1,k).and.k1.gt.nbs(k)) then


                ne   = ne + 1

                ati1 = y(k1 + nt)**2 + y(ntot + nt + k1)**2

                ion1 = ion1 + y(k1 + nt)**2 + y(ntot + nt + k1)**2

                                !xxxxxxxx   N = 2 threshold
                ati2 = 0.0d0
                if(en(k1,k).gt.thresh_2) then

                                ! Sum over channels at N = 2 
                   ati2 = 0.0d0
                   do j = 1, n2(k)

                     ati2 = ati2 + y(k1+nt+j)**2 + y(ntot+nt+k1+j)**2
                     ion2 = ion2 + y(k1+nt+j)**2 + y(ntot+nt+k1+j)**2

                   end do

                end if
                                !xxxxxxxx   N = 3 threshold
                ati3 = 0.0d0

                if(en(k1,k).gt.thresh_3) then

                  ati3 = 0.0d0
                                ! Sum over channels at N = 3 
                  do j = 1, n3(k)

                    ati3 = ati3 + &
                    &y(k1+nt+n2(k)+j)**2 + y(ntot+nt+k1+n2(k)+j)**2

                    ion3 = ion3 + &
                    &y(k1+nt+n2(k)+j)**2 + y(ntot+nt+k1+n2(k)+j)**2

                   end do

                end if

                write(40,'(5e13.5)') en(k1,k), ati1+ati2+ati3, ati1, ati2, ati3

             end if

          end do

          WRITE(44,'(I3,1X,5e13.5)') nl , ion1+ion2+ion3, ion1, ion2, ion3

          nt = nt + n(k)

       end do

       close(40)
       close(44)

    END 
    
!!!.....................................................................


!!!eof
