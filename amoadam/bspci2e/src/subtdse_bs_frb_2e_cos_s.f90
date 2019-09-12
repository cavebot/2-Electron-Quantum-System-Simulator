!!!imsl      subroutine fcn(neq,t,y,ypr)

      SUBROUTINE fcn(t,y,ypr)             !nag
      use parameter 
      use deriv
      implicit none
      INTEGER, PARAMETER:: nv=7000       !nag
      integer i, j, k, l, m
      REAL(8) t, t1, pi, slq, fi
      real(8), dimension(neq) :: y
      real(8), target, dimension(nv) :: ypr

!!!.......................................................................
      neval = neval + 1
      t1    = -0.5D+00 * tau + t
      pi    = 2.0 * dasin(1.0D+00)


      !  length.gauge : slq=sl*fi*dsin(omeg*t1)
      !  velocity.gauge : use




!!! gauss
!!!      slq=fi*dexp(-t1*t1/(2.0d0*tau*tau))*dcos(omeg*t1)/omeg  
!!! gauss

!!!  cos^2 

      slq = fi * dcos( t1/tau*pi )**2 * dcos( omeg*t1 ) / omeg

!!!  cos^2 

      DO i = 1, ntot

         ypr(i)      =   dag(i) * y(i + ntot)
         ypr(i+ntot) = - dag(i) * y(i)

      END DO


!!!    L = 0

      CALL dgemv ('t', n(2), n(1), -slq, dzr(1,1,1), nmax,& 
      &y(ntot+nsum(2)), 1, 1.0d0, ypr(1), 1)

      CALL dgemv ('t', n(2), n(1), slq, dzi(1,1,1), nmax,& 
      &y(nsum(2)), 1, 1.0d0, ypr(1), 1)

!!!.....

      CALL dgemv ('t', n(2), n(1), slq, dzr(1,1,1), nmax,&

      & y(nsum(2)), 1, 1.0d0, ypr(ntot+1), 1)
      CALL dgemv ('t', n(2), n(1), slq, dzi(1,1,1), nmax,&
           & y(ntot+nsum(2)), 1, 1.0d0, ypr(ntot+1), 1)

      DO k = 2, ltot - 1

         CALL dgemv ('t', n(1+k), n(k), -slq, dzr(1,1,k), nmax, &
              & y(ntot+nsum(k+1)), 1, 1.0d0, ypr(nsum(k)), 1)

         CALL dgemv ('t', n(1+k), n(k), slq, dzi(1,1,k), nmax, &
              & y(nsum(k+1)), 1, 1.0d0, ypr(nsum(k)), 1)

!!!.....
         CALL dgemv ('t', n(1+k), n(k), slq, dzr(1,1,k), nmax,&
              & y(nsum(k+1)), 1, 1.0d0, ypr(ntot+nsum(k)), 1)

         CALL dgemv ('t', n(1+k), n(k), slq, dzi(1,1,k), nmax,&
              & y(ntot+nsum(k+1)), 1, 1.0d0, ypr(ntot+nsum(k)), 1)
!!!.....
         CALL dgemv ('n', n(k), n(k-1), -slq, dzr(1,1,k-1), nmax, &
              &  y(ntot+nsum(k-1)), 1, 1.0d0, ypr(nsum(k)), 1)

         CALL dgemv ('n', n(k), n(k-1), -slq, dzi(1,1,k-1), nmax, &
              &  y(nsum(k-1)), 1, 1.0d0, ypr(nsum(k)), 1)
!!!.....
         CALL dgemv ('n', n(k), n(k-1), slq, dzr(1,1,k-1), nmax,&
              & y(nsum(k-1)), 1, 1.0d0, ypr(ntot+nsum(k)), 1)

         CALL dgemv ('n', n(k), n(k-1), -slq, dzi(1,1,k-1), nmax,&
              & y(ntot+nsum(k-1)), 1, 1.0d0, ypr(ntot+nsum(k)), 1)

      ENDDO


!!!    L = Lmax

      CALL dgemv ('n', n(ltot), n(ltot-1), -slq, dzr(1,1,ltot-1), nmax, &
      &  y(ntot+nsum(ltot-1)), 1, 1.0d0, ypr(nsum(ltot)), 1)

      CALL dgemv ('n', n(ltot), n(ltot-1), -slq, dzi(1,1,ltot-1), nmax, &
      &  y(nsum(ltot-1)), 1, 1.0d0, ypr(nsum(ltot)), 1)
!!!........
      CALL dgemv ('n', n(ltot), n(ltot-1), slq, dzr(1,1,ltot-1), nmax,&
      &  y(nsum(ltot-1)), 1, 1.0d0, ypr(ntot+nsum(ltot)), 1)

      CALL dgemv ('n', n(ltot), n(ltot-1),-slq, dzi(1,1,ltot-1), nmax,&
      &  y(ntot+nsum(ltot-1)), 1, 1.0d0, ypr(ntot+nsum(ltot)), 1)

      yderiv => ypr

990   FORMAT(2x,6e13.6)
      
    END SUBROUTINE fcn
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
!       if(l1.eq.6.and.l2.eq.7) open(n,file='/afs/ipp/u/jnz/mg/de60/os67.dat',&
 
      SUBROUTINE openfl(n,l1,l2)
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

     END SUBROUTINE openfl
!!!
!!! unformatted data
!!!
      subroutine openflbin(n,l1,l2)
       if(l1.eq.0.and.l2.eq.1) open(n,file='dat/os01.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.1.and.l2.eq.2) open(n,file='dat/os12.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.2.and.l2.eq.3) open(n,file='dat/os23.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.3.and.l2.eq.4) open(n,file='dat/os34.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.4.and.l2.eq.5) open(n,file='dat/os45.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.5.and.l2.eq.6) open(n,file='dat/os56.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.6.and.l2.eq.7) open(n,file='dat/os67.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.7.and.l2.eq.8) open(n,file='os78.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.8.and.l2.eq.9) open(n,file='os89.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.9.and.l2.eq.10) open(n,file='os910.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.10.and.l2.eq.11) open(n,file='os1011.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.11.and.l2.eq.12) open(n,file='os1112.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.12.and.l2.eq.13) open(n,file='os1213.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.13.and.l2.eq.14) open(n,file='os1314.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.14.and.l2.eq.15) open(n,file='os1415.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.15.and.l2.eq.16) open(n,file='os1516.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.16.and.l2.eq.17) open(n,file='os1617.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.17.and.l2.eq.18) open(n,file='os1718.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.18.and.l2.eq.19) open(n,file='os1819.dat',&    
      & form='unformatted',access='sequential')

      end subroutine openflbin
!!!EOF
