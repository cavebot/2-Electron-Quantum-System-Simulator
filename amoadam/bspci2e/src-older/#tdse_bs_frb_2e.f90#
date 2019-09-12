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
!   pulse duration starts from tau= 0. to --
!
!#   this program uses the NAG routine d02baf to integrate the diff. equation
!#   and it can be used for restarting an old job with a unit ncoe, with time and
!#   vector information.
!
!   this program uses the IMSL routine divprk to integrate the diff. equation
!   and it can be used for restarting an old job with a unit ncoe, with time and
!   vector information.
!   
!   
!
!   changes made by J. Zhang to multichannel TDSE,   feb  1996   
!   changes made by L. A.A. Nikolopoulos             June 2000
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM tdse_bs_frb_2e

        USE PRECISION
        USE UNITS
        USE PARAMETER
        USE deriv
        USE io
        
        IMPLICIT NONE
!
        REAL(dpk) de1, de2, thresh_1, thresh_2, thresh_3, ecut
!pulse
        REAL(dpk) i0, topt, tfac, t, t1, tim, tol
!
        REAL(dpk) pop_g, pop, pop_n, sion, dion, norm_wf
        REAL(dpk) p_0, p_1n, p_n
        REAL(dpk) p_n_single, p_n_double
        REAL(dpk) ati1, ati2, ati3, ion1, ion2, ion3, slqp
!
        INTEGER ifail, count_n
        INTEGER nl, ido, nemax
        INTEGER index, index_channel
        INTEGER il, i, j, p, q, k, k1, ni, nf, nb, lb,lz, iu, lstart, kt,lk
        INTEGER nout1, nout2, nt, li, l, ikt, ik, ii
        INTEGER it, ie
!
        INTEGER, ALLOCATABLE, DIMENSION(:,:)   ::noch
        INTEGER, ALLOCATABLE, DIMENSION(:)     ::ne
!
        INTEGER, DIMENSION(dpk)                :: n2, n3, nsgl
        INTEGER, DIMENSION(dpk)                :: nbs
!
        REAL(dpk), ALLOCATABLE, DIMENSION(:)   ::  y, yp
        REAL(dpk), ALLOCATABLE, DIMENSION(:,:) ::  energy
        REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: dpr, dpi, time, dipole
        REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: w         !nag
!
        INTEGER ndata
!        character(len=6) dmx_file
!
        CHARACTER*25 coe_t,ion_t,pop_t, dip_t, ati_t, pes_t, yld_t
        INTEGER  nion, npop, ncoe, ndip, nati, nyld, npes
        LOGICAL exists, hohg, pes
        EXTERNAL fcn, d02baf

        !      interface
        !         real function ddot(nd,x,nx,y,ny)
        !            integer nd,nx,ny
        !            real(8), dimension(nd) :: x,y
        !         end function ddot
        !      end interface
        !
        !#######################################################################

!    read from input file (infile)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        CALL getarg(1, argv)              !get arguments
        READ(argv,*)   i0                 ! in SI 
        CALL getarg(2, argv)              !get arguments
        READ(argv,*)   omeg               ! in SI
        CALL getarg(3, argv)              !get arguments
        READ(argv,*)   tfac
        CALL getarg(4, argv)              !get arguments
        READ(argv,*)   ltot
        CALL getarg(5, argv)              !get arguments
        READ(argv,*)   ndata


!        WRITE(dmx_file,'(i6)') n_data  ! convert integer 'run' to string 'run'

!

!        i0 = i0/i0_td             
        omeg = omeg/enau

!        READ(*,*) infile
        nion  = 10
        npop  = 11
        ncoe  = 12
        ndip  = 14
        nati  = 15
        npes  = 16
        nyld  = 17
!
        ion_t='tdat/ion.dat'
        pop_t='tdat/pop.dat'
        coe_t='tdat/coe.dat'
        dip_t='tdat/dip.dat'
        ati_t='tdat/ati.dat'
        pes_t='tdat/pes.dat'
        yld_t='tdat/yld.dat'
        !
        OPEN(nion,file=ion_t)
        IF(hohg)    OPEN(ndip,file=dip_t) 
        !
        OPEN(18,file="inp/tdse_bs_frb_2e.inp")
        READ(18,*)   nmax         ! Max no of States   in L = 0,1,2,. ltot-1
        READ(18,*)   nemax        ! Max no of Energies in L = 0,1,2,. ltot-1
        READ(18,*)   spulse       ! Max no of Energies in L = 0,1,2,. ltot-1
        CLOSE(18)

         pes = .false.
        hohg = .false.
        

!        READ(18,*)   coe_t        ! c_i(t)    : coefficients  
!        READ(18,*)   ion_t        ! Y(t)      : Y(+), Y(++)
!        READ(18,*)   pop_t        ! pop(t)    : populations
!        READ(18,*)   dipole_t     ! d(t)      : dipole
!        READ(18,*)   i0           ! Io(a.u.)  : peak Intensity
!        READ(18,*)   omeg         ! w(a.u.)   : photon energy
!        READ(18,*)   tfac         ! # of cycles
!        READ(18,*)   ltot         ! Total no of angul. included



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        e0   = SQRT(i0/i0_td)
        topt = 2.0d+00 * m_pi / omeg
        tau  = tfac * topt
        
        WRITE(nion,*)'# tdse_bs_mch::     peak intensity     Io = ',i0   , ' W/cm^2' 
        WRITE(nion,*)'# tdse_bs_mch::  photon frequency     Wph = ',omeg , ' a.u.'
        WRITE(nion,*)'# tdse_bs_mch::     pulse duration     Tp = ',tau  , ' a.u.'


! some initial checks


        OPEN(npop,file=pop_t)
        
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

      ALLOCATE(             ne(ltot) )
      ALLOCATE( energy( nemax, ltot) )
      ALLOCATE(   noch( nemax, ltot) )

!hohg
      IF(hohg) THEN
         ALLOCATE(            dpr(nmax) )
         ALLOCATE(            dpi(nmax) )
         ALLOCATE(  time(8192) )
         ALLOCATE(dipole(8192) )
      ENDIF
!hohg

!      ALLOCATE(ne(ltot), energy(nemax,ltot), noch(nemax,ltot))



!         Read data from dipole input files.
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!     lb         :  angular momentum of initial state
!     lz         :  angular momentum of final   state
!     n(l)       :  #  channel states in symmetry  L 
!     nbs(l)     :  # of bound states in symmetry  L 
!     ne(l)      :  # of energy states in symmetry L 
!  energy(i,l)   :  energy of i-th state of symmetry L 
!                :      i = 1,..., ne(l) , l = 1, ltot + 1.     
!     noch(i,l)  :  # of channels at i-th energy state of symmetry L 
!     dzr(i,j,l) :  real ( <i L | d | j L + 1> )
!     dzr(i,j,l) :  imag ( <i L | d | j L + 1> )
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      nl = ltot

      if(nl.le.1) stop

      DO i = 1, nl-1

         IF(ndata==0) THEN

            CALL openfl0(17,i-1,i)

            de1 = 0.008               ! Energy separation below E++ = 0. (in  a.u. )
            de2 = 0.0175              ! Energy separation above E++ = 0. (in  a.u. )
            
         ELSE IF(ndata==1) then

            CALL openfl1(17,i-1,i)
            de1 = 0.0175               ! Energy separation below E++ = 0. (in  a.u. )
            de2 = 0.0175               ! Energy separation above E++ = 0. (in  a.u. )

         ENDIF
                                ! Note: should check consistency lb+1 = i, 
                                !lz+1=i+1
         READ(17,*) lb,lz, n(i),n(i+1), nbs(i), nbs(i+1), ne(i), ne(i+1)
 
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
         WRITE(*,*) "# tdse_bs_mch::          dzr(i+1) = ",  dzr(100,1,i)
         WRITE(*,*) "# tdse_bs_mch::          dzi(i+1) = ",  dzi(100,1,i)

         CLOSE(17)

      END DO

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! Expand the energy array

!!!   en(i,j) :: energy of channel state i= 1,2, ... index of symmetry j.
!!!                Must be        index == n(i) 

      DO l = 1, nl

         noch(1:nbs(l), l) = 1

         index = 0
         DO k = 1, ne(l)
            DO j = 1, noch(k,l)

               index = index + 1
               
               en(index, l)    = energy(k,l)       ! in Ryd
              noch_en(index,l) =   noch(k,l)
          ch_index_en(index,l) =        j
     
       END DO
         END DO

      END DO

!         WRITE(*,*) "# tdse_bs_mch::                            neq = ", index 

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




      normalize_dmx: DO l = 1, nl-1


         WRITE(*,*) "# tdse_bs_mch::       normalize for dmx starting from l = ", i

!         BC transitions   d(n,m)   n : bound state      n of symmetry i 
!                                   m : continuum state  m of symmetry i+1

         normalize_dmx_bc: DO ni = 1, nbs(l)

            DO nf = nbs(l+1)+1, n(l+1)

               IF( en(nf,l+1).LE.0.0) THEN

                  dzr(nf,ni,l) = dzr(nf,ni,l) * SQRT(de1)
                  dzi(nf,ni,l) = dzi(nf,ni,l) * SQRT(de1)                  
               ELSE
                  dzr(nf,ni,l) = dzr(nf,ni,l) * sqrt(de2)
                  dzi(nf,ni,l) = dzi(nf,ni,l) * sqrt(de2)

!                  dzr(nf,1,l)  = dzr(nf,1,l)  * sqrt(0.0)
!                  dzi(nf,1,l)  = dzi(nf,1,l)  * sqrt(0.0)
               END IF

            END DO
         END DO normalize_dmx_bc
         
!         CB transitions   d(n,m)   n : continuum state  n of symmetry i 
!                                   m : bound     state  m of symmetry i+1

         normalize_dmx_cb: DO ni = nbs(l) + 1, n(l)
 
           DO nf = 1, nbs(l+1)

               IF(en(ni,l).LE.0.0) THEN

                  dzr(nf,ni,l) = dzr(nf,ni,l) * sqrt(de1)
                  dzi(nf,ni,l) = dzi(nf,ni,l) * sqrt(de1)
               ELSE

                  dzr(nf,ni,l) = dzr(nf,ni,l) * sqrt(de2)
                  dzi(nf,ni,l) = dzi(nf,ni,l) * sqrt(de2)
               END IF

            END DO
         END DO normalize_dmx_cb

!         CC transitions   d(n,m)   n : continuum state  n of symmetry i 
!                                   m : continuum state  m of symmetry i+1

         normalize_dmx_cc: DO ni = nbs(l)+1, n(l)

            DO nf = nbs(l+1)+1, n(l+1)

               IF(en(ni,l).GT.0.0.AND.en(nf,l+1).GT.0.0) THEN

                  dzr(nf,ni,l) = dzr(nf,ni,l) * de2 
                  dzi(nf,ni,l) = dzi(nf,ni,l) * de2

               ELSE
                  
                  IF(en(ni,l).LE.0.0.AND.en(nf,l+1).LE.0.0) THEN

                     dzr(nf,ni,l) = dzr(nf,ni,l) * de1
                     dzi(nf,ni,l) = dzi(nf,ni,l) * de1
                  ELSE

                  dzr(nf,ni,l) = dzr(nf,ni,l) * sqrt(de1 * de2)
                  dzi(nf,ni,l) = dzi(nf,ni,l) * sqrt(de1 * de2)
               END IF
            END IF
            
         END DO

      END DO normalize_dmx_cc

         WRITE(*,*) "# tdse_bs_mch::          dzr(i+1) = ",  dzr(100,1,l), l
         WRITE(*,*) "# tdse_bs_mch::          dzi(i+1) = ",  dzi(100,1,l), l
      ENDDO normalize_dmx

      ! end of re-normalization !

      ecut      = 5 !10.0d0               ! in a.u. 


!!!     Transform en to photoelectron energies         en --> k^2/2 
!!!     and ignore states with photoelectron energy    en > ecut.              
!!!      Energies in the calculation :                 en < ecut

!!! he+ thresholds

      thresh_1  = -4.0D+00             ! He+(1s)                 (in Ryd)
      thresh_2  = -1.0D+00             ! He+(2s,2p)
      thresh_3  = -0.4444444444D+00    ! He+(3s,3p,3d) 
      
      ntot = 0
      define_energy_origin: DO  l = 1, nl   !relative to first ionization threshold

          iu = 0
          DO  k1 = 1, n(l)

             en(k1, l) = ( en(k1, l) - thresh_1) / 2.0D+00

             IF(en(k1,l).GT.ecut) cycle
             iu = iu + 1
          ENDDO

          !       en < ecut
          n(l) = iu
          ntot = ntot + n(l)        ! nof states included 
          
          WRITE(nion,*)'# tdse_bs_mch::     angular momenta   l = ', k - 1
          WRITE(nion,*)'# tdse_bs_mch::     nof channels   n(l) = ', n(l)
          
       ENDDO define_energy_origin
       
       neq = 2*ntot
       ALLOCATE( y(neq), yp(neq))
       
       WRITE(*,*)    '# tdse_bs_mch::                              Emax = ', ecut
       WRITE(*, *)   '# tdse_bs_mch::   nof channels propagated   n_tot =  ', ntot


       !  Make 2-d energy array en(i,j) 1-d array dag(nb) 
       !
       !     en(k1,k) ----> dag(nb)       nb = k * n(k) + k1

       count_n = 0
       nb = 0
       DO l = 1, nl
          DO k1 = 1, n(l)
             
             nb      = nb + 1
             dag(nb) = en(k1, l)

             count_n = count_n + 1
          ENDDO
       ENDDO
       
       WRITE(*,*)  '# tdse_bs_mch::                e_g = ', en(1,1) 
       WRITE(*,*)  '# tdse_bs_mch ::           dagm(1) = ', dag(1)
       WRITE(*,*)  '# tdse_bs_mch ::        nof_states = ', count_n
       WRITE(*,*)  '#  t(au),             p_g,             norm'

       !xxx      Initial conditions at t = 0. y(1) : System in its ground state


       y    = 0.0D+00
       y(1) = 1.0D+00
       t    = 0.0D+00
       !
       nsum(1) = 1
       DO i = 2, ltot
          nsum(i) = nsum(i-1) + n(i-1)
       END DO
       WRITE(*,*) '# tdse_bs_mch::  ntot = nsum(ltot) ' , ntot, ' = ', nsum(ltot)
       
       kt = 0  
       
       IF(hohg)           ALLOCATE( yderiv(neq))

       ALLOCATE(      w( neq, 7))                               !!! work array for nag

!!!  Propagate     from 1 ----> 8*tfac  


       propagate_wf: DO  it = 1, INT(tfac*8)
          nout1 = kt + 1
          propagate_wf_inner: DO  lk = 1, 4    !break time-step for dipole(t)

             kt = kt + 1 
             tim = t + topt * 0.03125

!nag t-> t + dt

             neval = 0
             tol = 1.0d-8
             ifail = 0 

             CALL d02baf(t, tim, neq, y, tol, fcn, w, ifail)

             harmonics: IF(hohg) THEN
                time(kt) = tim
                dipole(kt) = 0.0d0
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
                !   slqp=-e0*(2*dcos(t1/tau*pi)*pi/tau*dsin(t1/tau*pi)&
                !     &*dcos(omeg*t1)/omeg&
                !     &+dcos(t1/tau*pi)**2*dsin(omeg*t1))
                !     dipole(kt)=dipole(kt) + slqp*(&
                !     &dot_product(y(nsum(li):nsum(li)+n(li)-1),dpr)&
                !     &+dot_product(y(ntot+nsum(li):ntot+nsum(li)+n(li)-1),dpi))
                END DO calculate_dipole_t
             ENDIF harmonics

          ENDDO propagate_wf_inner

!             deallocate(yderiv)
!          write(nion,*)'number of fun. evaluations used=', neval

          print_harmonic:IF(hohg) THEN       ! harmonics
             WRITE(*,*) '# tdse_bs_mch::  readout dipole(t) '
             nout2 = kt
             DO ikt = nout1, nout2
                t1 = -0.5d0 * tau + time(ikt)
                WRITE(ndip,'(2x,3e14.6)') time(ikt), dipole(ikt)
!0.5tau           write(ndip,'(2x,3e14.6)') time(ikt)-0.5e0*tau, dipole(ikt)
             END DO
          ENDIF print_harmonic                 !harmonics
          
!save coefficients data 

          WRITE(npop,'(e15.6)') tim

          !   nsg(i) :  # of single-ionization channels in symmetry L = i,  i.e. he+

          nsgl(1) = 9    ! i.e. : 1ss, (2pp, 2ss), (3dd, 3pp, 3ss), (7dd, 7pp, 7ss)
          nsgl(2) = 9    ! i.e. : 1sp, (2ps, 2pd, 2sp),(3dp, 3df, 3ps, 3pd, 3sp) 
          nsgl(3) = 13   ! i.e. : 1sd, (2pp, 2pf, 2sd, 2dd), (3dd, 3dg, 3pp, 3pf, 3sd)
                         !                                   (7dd  7pp  7sd)

          pop_g = y(1)**2 + y( ntot + 1 )**2
          norm_wf = 0.0d+00


          nt = 0
          population_monitor:DO  k = 1, nl            ! monitor population   
             IF(k.LE.MIN(10, nl)) THEN
                WRITE(npop,'(i3)') k - 1
                pop_n = y( k1 + nt )**2 + y( ntot + nt + k1 )**2 
                bound_states: DO  k1 = 1, 5        ! bound states population for nbs=1-4
                   WRITE(npop,'(3e15.6)') en(k1, k), pop_n
                ENDDO bound_states
                norm: DO  k1 = 1, n(k)          ! the whole population (norm)
                   norm_wf  = norm_wf + pop_n
                END DO norm
             ENDIF
             nt = nt + n(k)
          ENDDO population_monitor


!0.5tau          write(*,*)'t=', tim-0.5*tau, 'norm=', pop

          WRITE(*,'(3E15.8)') tim, pop_g, norm_wf
          write(npop,'(/)')

          !          write(ncoe) n(k),(en(k1,k),k1=1,n(k))

          pop  = 0.0d+00   !he
          sion = 0.0d+00   !he+
          dion = 0.0d+00   !he++

           p_0 = 0.0D+00        ! e_k < threshold_1
          p_1n = 0.0D+00        ! threshold_1 < e_k < threshold_n (n = infty)
          p_n  = 0.0D+00        ! e_k > threshold_n
          p_n_single = 0.0D+00  ! e_k > threshold_n (he+)
          p_n_double = 0.0D+00  ! e_k > threshold_n (he++)
          nt = 0

          index_channel = 1
          populations: DO il = 1, nl          ! loop on angular momenta
             DO  k1 = 1, n(il)                ! loop on energies
                
                pop_n = y( k1 + nt )**2 + y( ntot + nt + k1 )**2 


                ! (e_k < 0.0 a.u)       : bound states , (ground + excited He) 
                ! (0.0 < en < 2.0 a.u)  : strictly single ionization channels 
                ! (en > 2.0 a.u)        : single and double ionization channels

                IF(   en(k1,il).LE.0.0)                                     p_0 =  p_0  + pop_n                  
                IF( ( en(k1,il).GT.0.0D+00).AND.(en(k1,il).LT.2.0d+00) )    p_1n = p_1n + pop_n
                IF(   en(k1,il).GE.2.0D+00)                                 p_n =  p_n  + pop_n


                IF((en(k1,il).GE.2.0D+00) ) THEN 


                   IF (en(k1,il).EQ.en(k1-1,il))  THEN 
                      index_channel = index_channel + 1
                   ELSE
                         index_channel = 1
                   ENDIF

!                   WRITE(*,*) il,k1,en(k1,il), index_channel!

                   IF(index_channel.LT.nsgl(il)) THEN 

                      p_n_single = p_n_single + pop_n
                   ELSE
                      p_n_double = p_n_double + pop_n
                   ENDIF
                ENDIF

                                !                endif
             ENDDO

             nt = nt + n(il)

          ENDDO populations

!          stop

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


      WRITE(nion,908) tim, pop_g, pop, p_0, p_1n, p_n, sion, dion, p_n_single, p_n_double, p_n_single+p_n_double 
!          WRITE(nion,908) tim/topt, pop_g, 1.0 - pop - dion, dion

      pop = 0.0
      t = tim
      nt = 0

   ENDDO propagate_wf

   CLOSE(npop)
   CLOSE(nion)
   IF(hohg)   CLOSE(ndip)

!save coefficient data at final time.

   OPEN(nout,file="dat/yt.dat") 

   WRITE(nout,'(5e14.6)') i0,  omeg, tau
   WRITE(nout,'(2e14.6,1X,2I5)') tim, ecut, ltot, ntot
   nt = 0
   DO L = 1, LTOT
      WRITE(nout,'(2I8)')  N(L), L-1
      DO ie = 1, n(L)
         WRITE(nout,'(4I8,2X,3e14.6)') ie, L-1, noch_en(ie,L),ch_index_en(ie,L) &
         &,en(ie,L),y(ie+nt),y(ntot+nt+ie)!, y(ie+nt)**2+y(ntot+nt+ie)**2
      ENDDO
      nt = nt + n(L)
   ENDDO
   CLOSE(nout)

!save coefficient data at final time.

   OPEN(ncoe,file=coe_t,form='unformatted',access='sequential')
   WRITE(ncoe) i0, omeg, tau, ltot, ecut, tim
   WRITE(ncoe) (     n(l), l = 1, nl), ntot 
   DO l = 1, nl
      WRITE(ncoe) ( en(k1,l), k1 = 1, n(l))        !relative to he+(1s)
   END DO
   WRITE(ncoe) (    y(ii), ii = 1, neq)         ! save 
   CLOSE(ncoe)
   
! spectra
   spectra: IF(pes) THEN

      nt = 0
      OPEN(npes,file='tdat/pes.dat')
      
      WRITE(npes,*) ' spectrum for time = ', t/topt
      DO  k = 1, nl
         WRITE(npes,*) k - 1, nt, ntot
         DO  k1 = 1, n(k)            
            WRITE(npes,'(4e14.6)') en(k1,k),&
                 &(y(k1+nt)**2+y(ntot+nt+k1)**2)
         ENDDO
         nt = nt + n(k)
      ENDDO
      CLOSE(npes)
   END IF spectra
!xxxxxxxxxxxxxxxxxxxxxxxx
   
908 FORMAT( 1x, 12e14.6)
988 FORMAT( 2x, i4, i7, 5e14.6)
900 FORMAT( 2x, i7, 6(e10.3, 2x) )
400 FORMAT( 1x, '# tdse_bs_mch::  pulse duration (a.u) =', e15.8 )
300 FORMAT( 1x, '# tdse_bs_mch::  peak intensity (a.u) =', e15.8 )
   
   !xxxxxxxxxxxxxxxxxxxxxxx
   
   
999 CONTINUE
   
   
   
      
!          nsgl(1) = 9     i.e. : 1ss, (2) (2pp, 2ss),      (3) (3dd, 3pp, 3ss), (7)(7ss,7pp,7dd)
!          nsgl(2) = 9     i.e. : 1sp, (2) (2ps, 2pd, 2sp), (3) (3dp, 3df, 3ps, 3pd, 3sp) 
!          nsgl(3) = 13    i.e. : 1sd, (2) (2pp, 2pf, 2sd, ), 
                                !      (3) (3ds, 3dd, 3dg, 3pp, 3pf, 3sd)
                                !      (7dd  7pp  7sd)
   !
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
 
      OPEN(nati,file=ati_t)
      OPEN(nyld,file=yld_t)
      OPEN(1,   file="tdat/coe.out")

       DO  k = 1, nl

          ion1 = 0.0d0
          ion2 = 0.0d0
          ion3 = 0.0d0

          DO  k1 = 1, n(k)             ! all energies (n(k)) of symmetry k

             ! Increase by energy step (k1 !=  k1 - 1 ) 
             ! and don't take into account bound states

             IF(en(k1,k).NE.en(k1-1,k).AND.k1.GT.nbs(k)) THEN    !N = 1 (no degenerate spectrum)

                ne   = ne + 1
                
                ati1 =        y(k1 + nt)**2 + y(ntot + nt + k1)**2
                ion1 = ion1 + y(k1 + nt)**2 + y(ntot + nt + k1)**2


                                !xxxxxxxx   N = 2 threshold
                ati2 = 0.0d0
                IF(en(k1,k).GT.thresh_2) THEN

                                ! Sum over channels at N = 2 
                   ati2 = 0.0d0
                   DO j = 1, n2(k)

                     ati2 = ati2 + y( k1 + nt + j )**2 + y( ntot + nt + k1 + j )**2
                     ion2 = ion2 + y( k1 + nt + j )**2 + y( ntot + nt + k1 + j )**2

                  END DO

               END IF
                                !xxxxxxxx   N = 3 threshold
               ati3 = 0.0d0

               IF(en(k1,k).GT.thresh_3) THEN

                  ati3 = 0.0d0
                                ! Sum over channels at N = 3 
                  DO j = 1, n3(k)

                     ati3 = ati3 + y( k1 + nt + n2(k) + j )**2 + y( ntot + nt + k1 + n2(k) + j )**2
                     ion3 = ion3 + y( k1 + nt + n2(k) + j )**2 + y( ntot + nt + k1 + n2(k) + j )**2

                   end do

                end if

                write(nati,'(5e13.5)') en(k1,k), ati1 + ati2 + ati3, ati1, ati2, ati3
             END IF

             WRITE(1,'(4E13.5,1X,6I6)') en(k1,k),&
                  &  y( k1 + nt )**2 + y( ntot + nt + k1 )**2, &
                  &  y( k1 + nt ), y( ntot + nt + k1 ), &
                  &  nt+k1, ntot +nt+k1, k, k1, ntot 

                  END DO

          WRITE(nyld,'(I3,1X,5e13.5)') nl , ion1 + ion2 + ion3, ion1, ion2, ion3

          nt = nt + n(k)
       END DO

       CLOSE(1)
       CLOSE(nati)
       CLOSE(nyld)

     END PROGRAM tdse_bs_frb_2e
    
!!!.....................................................................


!!!eof
