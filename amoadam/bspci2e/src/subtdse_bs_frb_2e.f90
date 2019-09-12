!!!imsl      subroutine fcn(neq,t,y,ypr)

      SUBROUTINE fcn(t,y,ypr)              !nag
        !
        USE PRECISION
        USE units
        USE PARAMETER 
        USE deriv
        !
        IMPLICIT NONE
        !
        REAL(dpk),         DIMENSION(neq) :: y
        REAL(dpk), TARGET, DIMENSION(neq) :: ypr
        !LOC
        INTEGER                           :: i, j, k, l, m                !nag
        REAL(dpk)                         :: t, t1, slq
        !EXE!
      !

        neval = neval + 1


      ! pulse


      !  velocity.gauge : use


      t1 = -0.5d0 * tau + t


      !  length.gauge   : 

!      IF(gauge=='v') THEN
         slq = e0 * dcos( omeg*t1) / omeg
!      ELSE IF
!         slq= e0 * dsin(omeg*t1)
!      ENDIF

      IF(spulse=='Gauss') THEN 
         slq =  dexp( -t1*t1 / (2.0d0*tau*tau) ) * slq 
      ELSE if(spulse=='SinSqr') then
         slq = dcos(t1/tau*M_PI)**2* slq
      ENDIF

      OPEN(20, file='tdat/pulse.dat', position='append')
      WRITE(20,*) t1, slq, slq**2
      CLOSE(20) 

      ! free-field propagation
      
      do i = 1, ntot
         ypr(i)      =  dag(i) * y(i+ntot)
         ypr(i+ntot) = -dag(i) * y(i)
      end do

      !l_in  --> l_in + 1       (0-->1) 
      call dgemv ('t', n(2), n(1), -slq, dzr(1,1,1), nmax, y(ntot+nsum(2)), 1, 1.0d0, ypr(1), 1)
      call dgemv ('t', n(2), n(1),  slq, dzi(1,1,1), nmax, y(nsum(2)     ), 1, 1.0d0, ypr(1), 1)

      call dgemv ('t', n(2), n(1), slq, dzr(1,1,1), nmax, y(nsum(2)), 1, 1.0d0, ypr(ntot+1), 1)
      call dgemv ('t', n(2), n(1), slq, dzi(1,1,1), nmax, y(ntot+nsum(2)), 1, 1.0d0, ypr(ntot+1), 1)



      DO k = 2, ltot - 1             !l  --> l + 1

         call dgemv ('t', n(1+k), n(k), -slq, dzr(1,1,k), nmax, &
         & y(ntot+nsum(k+1)), 1, 1.0d0, ypr(nsum(k)), 1)
         call dgemv ('t', n(1+k), n(k), slq, dzi(1,1,k), nmax, &
         & y(nsum(k+1)), 1, 1.0d0, ypr(nsum(k)), 1)

         call dgemv ('t', n(1+k), n(k), slq, dzr(1,1,k), nmax,&
         & y(nsum(k+1)), 1, 1.0d0, ypr(ntot+nsum(k)), 1)
         call dgemv ('t', n(1+k), n(k), slq, dzi(1,1,k), nmax,&
         & y(ntot+nsum(k+1)), 1, 1.0d0, ypr(ntot+nsum(k)), 1)

         call dgemv ('n', n(k), n(k-1), -slq, dzr(1,1,k-1), nmax, &
         &  y(ntot+nsum(k-1)), 1, 1.0d0, ypr(nsum(k)), 1)
         call dgemv ('n', n(k), n(k-1), -slq, dzi(1,1,k-1), nmax, &
         &  y(nsum(k-1)), 1, 1.0d0, ypr(nsum(k)), 1)

         call dgemv ('n', n(k), n(k-1), slq, dzr(1,1,k-1), nmax,&
         & y(nsum(k-1)), 1, 1.0d0, ypr(ntot+nsum(k)), 1)
         call dgemv ('n', n(k), n(k-1), -slq, dzi(1,1,k-1), nmax,&
         & y(ntot+nsum(k-1)), 1, 1.0d0, ypr(ntot+nsum(k)), 1)

      ENDDO

      ! l_max - 1 --> lmax
      CALL dgemv ('n', n(ltot), n(ltot-1), -slq, dzr(1,1,ltot-1), nmax,&
           y(ntot+nsum(ltot-1)), 1, 1.0d0, ypr(nsum(ltot)), 1)
      CALL dgemv ('n', n(ltot), n(ltot-1), -slq, dzi(1,1,ltot-1), nmax,&
           &y(nsum(ltot-1)), 1, 1.0d0, ypr(nsum(ltot)), 1)

      call dgemv ('n', n(ltot), n(ltot-1), slq, dzr(1,1,ltot-1), nmax,&
      &  y(nsum(ltot-1)), 1, 1.0d0, ypr(ntot+nsum(ltot)), 1)
      call dgemv ('n', n(ltot), n(ltot-1),-slq, dzi(1,1,ltot-1), nmax,&
      &  y(ntot+nsum(ltot-1)), 1, 1.0d0, ypr(ntot+nsum(ltot)), 1)

      yderiv => ypr

 990  format(2x,6e13.6)

    END SUBROUTINE fcn

    SUBROUTINE openfl0(n,l1,l2)
      IF(l1.EQ.0.AND.l2.EQ.1)   OPEN(n,file='dat/os01.v0.dat')
      IF(l1.EQ.1.AND.l2.EQ.2)   OPEN(n,file='dat/os12.v0.dat')
      IF(l1.EQ.2.AND.l2.EQ.3)   OPEN(n,file='dat/os23.v0.dat')
      IF(l1.EQ.3.AND.l2.EQ.4)   OPEN(n,file='dat/os34.v0.dat')
      IF(l1.EQ.4.AND.l2.EQ.5)   OPEN(n,file='dat/os45.v0.dat')
    END SUBROUTINE openfl0


    SUBROUTINE openfl1(n,l1,l2)
      IF(l1.EQ.0.AND.l2.EQ.1)   OPEN(n,file='dat/os01.v1.dat')
      IF(l1.EQ.1.AND.l2.EQ.2)   OPEN(n,file='dat/os12.v1.dat')
      IF(l1.EQ.2.AND.l2.EQ.3)   OPEN(n,file='dat/os23.v1.dat')
      IF(l1.EQ.3.AND.l2.EQ.4)   OPEN(n,file='dat/os34.v1.dat')
      IF(l1.EQ.4.AND.l2.EQ.5)   OPEN(n,file='dat/os45.v1.dat')
    END SUBROUTINE openfl1




