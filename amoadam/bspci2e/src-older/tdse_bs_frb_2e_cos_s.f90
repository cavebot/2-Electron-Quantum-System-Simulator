!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx     
!                                         written by xian tang, jan 1990
!                                          changes made by  hr, feb 1990
!   changes made 3/31/90 by hr and xtang:
!   changed to remove the rapid oscillating part of the matrix
!   pulse shape changed to gaussian
!   pulse duration changed from -0.5* tau to 0.5*tau
!
!   this program uses the nag routine d02baf to integrate the diff. equation
!   and it can be used for restarting an old job with a unit 33, with time and
!   vector information
!
!     revised to FORTRAN 90 by Jian Zhang, Feb 1996 
!!!   changes made by  Lambros Nikolopoulos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program prop

      use parameter
      use deriv
      implicit none
!nag     
         INTEGER, PARAMETER    :: nv=7000
         REAL(8), DIMENSION(nv)::  y, yp
!nag

!imsl
!      integer, parameter :: nv=20000, nbb=20*2000
!      real(8) rwksp(nbb)
!      common/worksp/rwksp
!imsl
      real(8), dimension(50) :: parp
      REAL(8) ff, topt, tfac, pi, er, ecut, t, t1, tim, pop, tol, de1, fi 
      INTEGER ifail, tSteps                              !nag
      integer, dimension(8) :: nbs
      integer irestart, irespec, nl, ido
      INTEGER i, j, p, q, k, k1, ni, nf, nb, lb,lz, iu, lstart, kt,lk
      integer nout1, nout2, nt, li, l, ikt, ik, ii

!imsl      real(8), allocatable, dimension(:) :: y, yp, dpr, dpi, time, dipole

!!!      REAL(8), ALLOCATABLE, DIMENSION(:) :: dpr, dpi, time, dipole
      REAL(8), ALLOCATABLE, DIMENSION(:) :: dpr, dpi 
      REAL(8), ALLOCATABLE, DIMENSION(:,:) :: w         !nag
      character*25 coef_t,ion_t,pop_t, dipole_t
      logical exists
!nag
      external fcn, d02baf 
!nag

!imsl
!      external fcn, divprk
!     call iwkin(nbb)
!imsl

!      interface
!         real function ddot(nd,x,nx,y,ny)
!            integer nd,nx,ny
!            real(8), dimension(nd) :: x,y
!         end function ddot
!      end interface
!
!---------------------------------------------------------------------
!      write(*,'(a30, $)')' name of output file = '
!      call link("unit18=ntest4.in//")


!!!  input file
!!!...................................................................

!!!      read(*,*) infile
      OPEN(18,file='inp/tdse.inp')

      read(18,'(a20)')  coef_t
      read(18,'(a20)')  ion_t
      read(18,'(a20)')  pop_t
      read(18,'(a20)')  dipole_t
      read(18,*) ff
      read(18,*) omeg
      read(18,*) tfac
      read(18,*) ltot
      read(18,*) irestart
      read(18,*) irespec
      read(18,*) nmax

      close(18)

!!!...................................................................

      pi = 4.0e0 * ATAN(1.0e0)                ! pi = 3.141...
      fi   = sqrt(ff)                         ! E_o
      topt = 2.d0 * pi / omeg                 ! T_w  
      tau  = tfac * topt                      ! T_L


      open(10,file=ion_t)

      write(10,*)' Io  = ', ff,    ' a.u.'
      write(10,*)' Wph = ', omeg , ' a.u.'
      write(10,*)' T   = ', tau,   ' a.u.' 

!!!
!!!   restart or not ?
!!!
      if(irestart.eq.1.and.irespec.ne.1) then

         WRITE(10,*) '########  RESTARTED RUN ###########'

      endif



      inquire(file=coef_t,exist=exists)

!!!
!!!   check to see if the output file exists, if this is a restart
!!!


      if((irestart.eq.1).and..not.exists) then

 200     write(10,*)' warning: input claims this to be a restart, but your'
         write(10,*) ' unit33 does not exist, the program will now abort'
         stop

      end if

!!!
!!!   check to make sure that a output file does not exist when not a restart
!!!


      if((irestart.ne.1).and.exists) then
 201     write(10,*)' warning: input claims this not to be a restart, but '
         write(10,*)' your unit33 is not empty, the program will now abort'
         stop
      end if


!!!
!!!  Oen file with the time-dependent coefficients
!!!

      open(33,file=coef_t,form='unformatted',access='sequential')
      open(34,file=pop_t)

!!!
!!!  allocate space for matrices
!!!


      call make_space(nmax,ltot)


!!!imsl
      !
      !      allocate(y(nmax*ltot*2),yp(nmax*ltot*2))
!!!imsl      

!!!      ALLOCATE(dpr(nmax), dpi(nmax), time(8192), dipole(8192))

      ALLOCATE(dpr(nmax), dpi(nmax))

      
      nl = ltot

      IF(nl.LE.1) STOP

! 17   format(2x,1p8e15.7)


!!!
!!!   Number of bound states for each of the angular momenta 
!!!
!!!   L = 0 ---> 6
!!!   L = 1 ---> 5
!!!   ............
!!!   ............
!!!

      nbs(1) = 6 
      nbs(2) = 5 
      nbs(3) = 5  
      nbs(4) = 4  
      nbs(5) = 3  
      nbs(6) = 2  
      nbs(7) = 1


!!!
!!!  open files with energies and DME   E, < a | D | b > 
!!!

!      DO 899 i = 1, nl - 1
      DO  i = 1, nl - 1


         CALL openfl(17,i-1,i)

!!!         READ(17) lb,lz,n(i),n(i+1)

         READ(17,*) lb, lz, n(i), n(i+1)

         IF(n(i).GT.nmax)    WRITE(*,*) ' ni     > nmax', n(i)
         IF(n(i+1).GT.nmax)  WRITE(*,*) ' n(i+1) > nmax', n(i+1)


!!!         READ(17) ( en(k,lb+1), k = 1, n(i)  )
!!!         READ(17) ( en(k,lz+1), k = 1, n(i + 1) )

         READ(17,*) ( en(k,lb+1), k = 1, n(i)  )
         READ(17,*) ( en(k,lz+1), k = 1, n(i + 1) )

            WRITE(*,*) 'READ   DME(',1, '-',n(i),n(i+1), ')'
         DO  ni = 1, n(i)

!!!            READ(17) (dzr(nf,ni,i), nf=1,n(i+1))
!!!            READ(17) (dzi(nf,ni,i), nf=1,n(i+1))

            READ(17,*) (dzr(nf,ni,i), nf=1,n(i+1))
            READ(17,*) (dzi(nf,ni,i), nf=1,n(i+1))

!!!            WRITE(*,*) 'READ   DME(', i, ',', ni,')'

         END DO

         CLOSE(17)

         WRITE(*,*) 'READ   DME(', i, ',', i + 1,') '

!!!
!!!     Renormalize DME 
!!!

         de1 = 0.003

!!!
!!!        bound-free:    < n_b l | d | n_k l+1 > ---> < b | d | n_k l+1 > * norm( n_k, l + 1)
!!!

         DO ni = 1, nbs(i)
            DO nf = nbs(i + 1) + 1, n(i + 1)

               dzr(nf,ni,i) =  dzr(nf,ni,i) * SQRT(de1)
               dzi(nf,ni,i) =  dzi(nf,ni,i) * SQRT(de1)

            END DO
         END DO


!!!
!!!        free-bound:    < n_k l | d | n_b l+1 > --->  norm( n_k, l ) * < b | d | n_b l+1 > 
!!!

         DO ni = nbs(i) + 1, n(i)
            DO nf = 1, nbs(i + 1)

               dzr(nf,ni,i) = dzr(nf,ni,i) * SQRT(de1)
               dzi(nf,ni,i) = dzi(nf,ni,i) * SQRT(de1)

            END DO
         END DO


!!!
!!!        free-free:    < n_k l | d | n_k l+1 > ---> 
!!!
!!!                    ---->  norm( n_k, l ) * < b | d | n_b l+1 > * norm( n_k, l+1 )
!!!

         DO ni = nbs(i)+1, n(i)

            DO nf = nbs(i+1) + 1, n(i+1)

               dzr(nf,ni,i) = dzr(nf,ni,i) * de1
               dzi(nf,ni,i) = dzi(nf,ni,i) * de1

            END DO
         END DO

      ENDDO
!!899      CONTINUE

!      ecut=30.0d0
!      ecut=-er

!!!
!!!  Ignore states with total energy larger than   ECUT
!!!

         ecut = 2.5D+00

!!!       n(1) = 115;  n(2) = 114 ;  n(3) = 114

!!!
!!!   if restart get stored information from the previous run
!!!

         IF(irestart.EQ.1) THEN

            READ(33) ff,omeg,tau,ltot,ecut
         ELSE

            WRITE(33) ff,omeg,tau,ltot,ecut
         ENDIF


!!!
!!!   Set as Zero energy the grouns state of the atom.
!!!
!!!     
!!!    en(k1, l) > 0   k1-continuum state of the atom with angular momenta l
!!!    en(k1, l) < 0   k1-bound state of the atom with angular moment l   
!!!


       ntot = 0
       er = -1.1050719D+00    ! ground state of Magnesium


       DO 21 k = 1, nl

          iu = 0
          
          DO 22 k1 = 1, n(k)

             en(k1,k) = (en(k1,k) - er) / 2.d0

!!!   exclude   E > ecut

             IF(en(k1,k).GT.ecut) go to 22
             iu = iu + 1

22           CONTINUE

          n(k) = iu

!           en(1,1)=-0.5619705e0/2.e0

!!!
!!!  if restart get information from stored values
!!!

          IF(irestart.EQ.1) THEN

             read(33) n(k), ( en(k1,k), k1 = 1, n(k))
          ELSE

             write(33) n(k), ( en(k1,k), k1 = 1, n(k))
          ENDIF

!!!
!!!    total number of states (size of the time-dependent wavefunction)
!!!


          ntot = ntot + n(k)

          WRITE(10,*)' L = ', k - 1,' n(L) = ', n(k)

21        CONTINUE



!!!
!!!  if restart get information from stored values
!!!

          IF(irestart.EQ.1) THEN
          READ(33) ntot
       ELSE
          WRITE(33) ntot
       ENDIF

       WRITE(10,*)'Size of TD wavefunction = ', ntot
       WRITE(*,*) 'Size of TD wavefunction = ', ntot

!!!
!!!  Make the diagonial part of matrix  E(nl,nl)
!!!

       nb = 0
       DO  k = 1, nl
          DO  k1 = 1, n(k)

             nb = nb + 1
             dag(nb) = en(k1,k)

          ENDDO
       ENDDO

       WRITE(*,*) ' DAGM(1) = ', dag(1)

!!!
!!!   if restart get information from stored values
!!!

       IF(irestart.EQ.1) THEN

          t = -0.5D+00 * tau

215       READ(33,END=210) t1, ( yp(ii), ii = 1, 2 * ntot)

          t = t1

          DO 216 ik = 1, 2 * ntot

             y(ik) = yp(ik)

216          CONTINUE

             GOTO 215

210          lstart = INT(0.5 * ( t + 0.5 * tau ) /topt ) + 1

             t = t + 0.5 * tau

          ELSE

             lstart = 1

          ENDIF

!!!
!!!   Initial conditions if it is not a restart : system in its ground state  
!!!

          IF(lstart.EQ.1) THEN

             y    = 0.0D+00
             
             y(1) = 1.0D+00

             t    = 0.0D+00

          ENDIF



!!!imsl          ido = 1

          IF(irespec.EQ.1) GOTO 249

!!!imsl          parp(4) = 1000000
          nsum(1) = 1

          DO i = 2, ltot

             nsum(i) = nsum(i-1) + n(i-1)

          END DO
      
          kt = 0

!!!          OPEN(41,file=dipole_t)

          ALLOCATE( yderiv( 2 * ntot) )     !!! dy/dt
          ALLOCATE( w( nv, 7))           !!! work array for nag

          DO 100 l = lstart, INT( tfac * 0.5 )


!!!****** break more steps for calculating the time-dependent dipole ***

             tSteps = 1  !!!178
             nout1 = kt + 1

             DO  lk = 1, tSteps

                kt = kt + 1 

                WRITE(*,*) "lstart = ",lstart

                tim = t + topt * (2.0D+00/tSteps)

             neval = 0

!begin{nag}
!!!
!!!   tolerance
!!!
             tol = 1.0D-10
             ifail = 0 

             CALL d02baf(t, tim, 2*ntot, y, tol, fcn, w, ifail)

!end{nag}

!imsl             call divprk(ido,2*ntot,fcn,t,tim,tol,parp,y)


!!!
!!!         d(t) = < psi(t) | d(t) | psi(t) >
!!!

!!$             time(kt)   = tim
!!$             dipole(kt) = 0.0D+00
!!$
!!$             DO li = 1, nl - 1
!!$                
!!$                !!!     <f|p|df/dt>
!!$
!!$                CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,&
!!$                &yderiv(nsum(li+1)),1,0.0d0,dpr,1)
!!$
!!$                CALL dgemv('t',n(li+1),n(li),-1.0d0,dzi(1,1,li),nmax,&
!!$                &yderiv(ntot+nsum(li+1)),1,1.0d0,dpr,1)
!!$
!!$                CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,&
!!$                &yderiv(ntot+nsum(li+1)),1,0.0d0,dpi,1)
!!$
!!$                CALL dgemv('t',n(li+1),n(li),1.0d0,dzi(1,1,li),nmax,&
!!$                &yderiv(nsum(li+1)),1,1.0d0,dpi,1)
!!$
!!$                dipole(kt)=dipole(kt)&
!!$                &+dot_product(y(nsum(li):nsum(li)+n(li)-1),dpr)&
!!$                &+dot_product(y(ntot+nsum(li):ntot+nsum(li)+n(li)-1),dpi)
!!$
!!$                !!!     <df/dt|p|f>
!!$
!!$                CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,&
!!$                &y(nsum(li+1)),1,0.0d0,dpr,1)
!!$
!!$                CALL dgemv('t',n(li+1),n(li),-1.0d0,dzi(1,1,li),nmax,&
!!$                &y(ntot+nsum(li+1)),1,1.0d0,dpr,1)
!!$
!!$                CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,&
!!$                &y(ntot+nsum(li+1)),1,0.0d0,dpi,1)
!!$
!!$                CALL dgemv('t',n(li+1),n(li),1.0d0,dzi(1,1,li),nmax,&
!!$                &y(nsum(li+1)),1,1.0d0,dpi,1)
!!$
!!$                dipole(kt)=dipole(kt)&
!!$                &+dot_product(yderiv(nsum(li):nsum(li)+n(li)-1),dpr)&
!!$                &+dot_product(yderiv(ntot+nsum(li):ntot+nsum(li)+n(li)-1),dpi)
!!$
!!$             END DO

!!99           CONTINUE
          ENDDO
!!!             deallocate(yderiv)

!!!             WRITE(10,*)' # of fun. evaluations used = ', neval

!!$
!!$             nout2 = kt
!!$
!!$             DO ikt = nout1, nout2
!!$
!!$                WRITE(41,'(2x,3e14.6)') time(ikt) - 0.5e0*tau, dipole(ikt)
!!$
!!$             END DO

!!!
!!!   save for later restart
!!!
             IF(irestart.EQ.1) THEN

                REWIND (33)

                WRITE(33) ff,omeg,tau,ltot,ecut

                DO 260 k = 1, nl

                   WRITE(33) n(k), ( en(k1,k), k1 = 1, n(k) )

260                CONTINUE

                   WRITE(33) ntot

                   PRINT*, 'wrote ntot', ntot
                   
                   irestart = 0
                ENDIF

                WRITE(33) tim - 0.5e0 * tau, ( y(ii), ii = 1, 2*ntot )


!!!
!!!  write population for partial waves
!!!

                WRITE(34,'(e15.6)') tim-0.5*tau

                pop = 0.0D+00

                DO 265 k = 1, nl

                   IF(k.LE.MIN(10,nl)) THEN

!!!
!!!    angular momenta of the partial wave
!!!
                      WRITE(34,'(i3)') k - 1


                      DO 266 k1 = 1, 5

!!!
!!!                  | Y(nL) |^2  , n = 1-5, L = 0, Lmax
!!!

                         WRITE(34,'(3e15.6)') en(k1,k), (y(k1+nt)**2+y(ntot+nt+k1)**2)

266                      CONTINUE

!!!
!!!                  n(L)_max
!!!                  --
!!!       P(L) =     \   | Y(nL) |^2       L = 0, Lmax
!!!                  /__
!!!                   1  
!!!

                         DO  k1 = 1, n(k)

                            pop = pop + y(k1 + nt)**2 + y( ntot + nt + k1 )**2

                         END DO


                      ENDIF

                      nt = nt + n(k)

265                   CONTINUE

                      WRITE(*,*)' Norm( t = ', tim - 0.5 * tau,') = ', pop

                      WRITE(34,'(/)')

                      nt = 0

                      pop = 0.0D+00



!comment          write(33) n(k),(en(k1,k),k1=1,n(k))

!!!
!!!          bound state population
!!!
!!!                 Lmax   e(nL) < 0
!!!                 ---     --
!!!       Pb =      \     \   | Y(nL) |^2       L = 0, Lmax
!!!                 /__   /__
!!!                  0     1  
!!!


              DO  k = 1, nl
                 DO  k1 = 1, n(k)

                    IF( en(k1,k).LE.0.0) THEN
                               
                       pop = pop + ( y(k1 + nt)**2 + y( ntot+nt+k1)**2 )

                    ENDIF

                 ENDDO

                 nt = nt + n(k)

              ENDDO


!!!
!!!         t, P_g, 1 - pop
!!!

              WRITE(10,908) tim - 0.5D+00 * tau, ( y(1)**2+y(ntot+1)**2), 1.0D+00 - pop

              pop = 0.0D+00
              t   = tim
              nt  = 0

100           CONTINUE


              CLOSE(33)
              CLOSE(34)

!!!              GOTO 999

!!!
!!!   skip PES = PES(t)
!!!


 249   nt = 0



              WRITE(*,*) ' PES AT TIME t = ', t - 0.5D+00 * tau

              OPEN(40,file='dat/pes.dat')

              DO  k = 1, nl

                 WRITE(40,*) ' PES FOR  L = ', k - 1

                 DO  k1 = 1, n(k)

                    WRITE(40,'(4e14.6)') en(k1,k),&
                         & ( y(k1+nt)**2 + y(ntot+nt+k1)**2), y(k1+nt), y(ntot+nt+k1)
                 ENDDO


                 nt = nt + n(k)

              ENDDO

              CLOSE(40)

!!! format statements
!!!.............................................

 908  format(1x,6e14.6)
 988  format(2x,i4,i7,5e14.6)
 900  format(2x,i7,6(e10.3,2x))
 400  format(1x,'PULSE DURATION (a.u) = ',e15.8)
 300  format(1x,'PEAK INTENSITY (a.u) = ',e15.8)

!!!...........................................................

!!! 999  
!!!       CLOSE(41)

           END 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!EOF
