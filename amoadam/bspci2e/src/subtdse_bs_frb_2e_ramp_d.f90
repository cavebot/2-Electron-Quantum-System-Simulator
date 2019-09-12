!!!imsl      subroutine fcn(neq,t,y,ypr)
      subroutine fcn(t,y,ypr)          !!!nag
      use parameter 
      use deriv
      implicit none
      INTEGER,       PARAMETER::  nv = 7000
      integer i, j, k, l, m                !nag
      REAL(8) t, t1, pi, slq,rc, fi
      real(8), dimension(neq) :: y
      real(8), target, dimension(nv) :: ypr
!!!imsl      integer i, j, k, l, m, neq 

      neval=neval+1
      t1=-0.5d0*tau+t
      pi=2.0*dasin(1.0d0)
                                !  length.gauge : slq=sl*fi*dsin(omeg*t1)
                                !  velocity.gauge : use
!.gau      slq=fi*dexp(-t1*t1/(2.0d0*tau*tau))*dcos(omeg*t1)/omeg
!.cos      slq=fi*dcos(t1/tau*pi)**2*dcos(omeg*t1)/omeg
!      slq=fi*dcos(omeg*t1)/omeg
!      if(t.le.0.2*tau) slq=fi*dcos(omeg*t1)/omeg*t*5.0/tau

      rc = 3.0 
      slq=fi*dsin(omeg*t1)/omeg

      if(t.le.(tau/rc))   slq = fi * dsin(omeg*t1)/ omeg * t * rc/tau

      if(t.gt.((rc-1)*tau/rc)) slq = fi * dsin(omeg*t1)/ omeg * (tau - t)* rc/tau

	open(30,file='pulse.dat')
	write(30,*) t, slq
 
      do i=1,ntot
         ypr(i)=dag(i)*y(i+ntot)
         ypr(i+ntot)=-dag(i)*y(i)
      end do
      call dgemv ('t', n(2), n(1), -slq, dzr(1,1,1), nmax,& 
      &y(ntot+nsum(2)), 1, 1.0d0, ypr(1), 1)
      call dgemv ('t', n(2), n(1), slq, dzi(1,1,1), nmax,& 
      &y(nsum(2)), 1, 1.0d0, ypr(1), 1)

      call dgemv ('t', n(2), n(1), slq, dzr(1,1,1), nmax,&
      & y(nsum(2)), 1, 1.0d0, ypr(ntot+1), 1)
      call dgemv ('t', n(2), n(1), slq, dzi(1,1,1), nmax,&
      & y(ntot+nsum(2)), 1, 1.0d0, ypr(ntot+1), 1)
      do 2 k=2,ltot-1
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
 2    continue
      call dgemv ('n', n(ltot), n(ltot-1), -slq, dzr(1,1,ltot-1), nmax, &
      &  y(ntot+nsum(ltot-1)), 1, 1.0d0, ypr(nsum(ltot)), 1)
      call dgemv ('n', n(ltot), n(ltot-1), -slq, dzi(1,1,ltot-1), nmax, &
      &  y(nsum(ltot-1)), 1, 1.0d0, ypr(nsum(ltot)), 1)

      call dgemv ('n', n(ltot), n(ltot-1), slq, dzr(1,1,ltot-1), nmax,&
      &  y(nsum(ltot-1)), 1, 1.0d0, ypr(ntot+nsum(ltot)), 1)
      call dgemv ('n', n(ltot), n(ltot-1),-slq, dzi(1,1,ltot-1), nmax,&
      &  y(ntot+nsum(ltot-1)), 1, 1.0d0, ypr(ntot+nsum(ltot)), 1)

      yderiv => ypr
 990  format(2x,6e13.6)

      end subroutine fcn
! 
      subroutine openfl(n,l1,l2)
       if(l1.eq.0.and.l2.eq.1) open(n,file='../dat.L2/os01.dat',&
      & form='unformatted',access='sequential')
      print*,'file opened'
       if(l1.eq.1.and.l2.eq.2) open(n,file='../dat.L2/os12.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.2.and.l2.eq.3) open(n,file='os23.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.3.and.l2.eq.4) open(n,file='os34.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.4.and.l2.eq.5) open(n,file='os45.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.5.and.l2.eq.6) open(n,file='os56.dat',&
      & form='unformatted',access='sequential')
       if(l1.eq.6.and.l2.eq.7) open(n,file='os67.dat',&
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
      end subroutine openfl

