! functions taken from mcdh code
!
! *             
! *  LEGENDRE and SHERICAL-HARMONICS FUNCTIONS (legendre.f)            
! *                                                                    
! * Library module containing the following functions and subroutines: 
! *                                                                    
! *   plgndr:  computes the associated Legendre polynomial P_l^m (x).  
! *            Here m and l are integers satisfying 0<=m<=l, while x   
! *            lies in the range -1<=x<=1   (AJ, 8/96)                 
! *   dlegendre: computes P_n(cos(theta)) where P_n(x) is the Legendre 
! *              polynomial of order n.  (MB, 8/96)                    
! *   sphharm: Computes the spherical harmonics  Y_l^m (theta,phi).    
! *            Here m and l are integers satisfying 0<=l and -l<=m<=+l 
! *            general spherical harmonics function with               
! *            Condon-Shortley phase. (MCH 97)                         
! *   algndr:  computes the associated Legendre polynomial P_l^m (x).  
! *            Similar to plgndr but computes an array containing all   
! *            P_l^m (x) for m<=l<=lmax. This subroutine computes      
! *            the L^2-mormalized functions, if switch=.true. is set.  
! *                                                                    
!
module spherical_harmonics
  !
  use precision
  
!  public plgndr, dlegendre, sphharm, algndr

contains


!-----------------------------------------------------------------------
!                FUNCTION   PLGNDR
!
! PURPOSE:
!   Computes the associated Legendre polynomial P_l^m (x). Here m and l
!   are integers satisfying 0<=m<=l, while x lies in the range
!   -1<=x<=1
! INPUT PARAMETERS:
!   l: integer (see above)
!   m: integer (see above)
!   x: real*8  (see above)
! FUNCTION VALUE:
!     value of the "plgndr" function as defined above
!
! PHASE CONVENTION:
!   The phase convention differs by a factor (-1)**m from most textbooks
!   as e.g. Messiah (Quantum Mechanics) or Zare (Angular Momentum).
!-----------------------------------------------------------------------

  function plgndr(l,m,x)

    implicit none 
    
    integer l,m
    real*8  plgndr,x
    
    integer i,ll
    real*8  fact,pll,pmm,pmmp1,somx2
    
    if (m.lt.0.or.m.gt.l.or.abs(x).gt.1) then
       write(6,*) '# plgndr::   error: bad arguments in plgndr ', m, l, x
       write(2,*) '# plgndr::   error: bad arguments in plgndr ', m, l, x
       stop
    endif
    
    ! --- Compute P_m^m ---
    
    pmm=1.
    if (m.gt.0) then
       somx2=sqrt((1.d0-x)*(1.d0+x))
       fact=1.
       do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
       enddo
    endif
    if (l.eq.m) then
       plgndr=pmm
    else



       pmmp1 = x*(2*m+1)*pmm
       
       IF (l.EQ.m+1) THEN ! --- Compute P_m+1^m ---
          
          plgndr=pmmp1
          
       ELSE ! --- Compute P_l^m, l>m+1 ---


          DO ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          
          plgndr=pll
       endif
    endif
    
    return
  end function plgndr

! ----------------------------------------------------------------------
!                            FUNCTION DLEGENDRE
!
! PURPOSE:
!     Computes the value P_n(cos(theta)) of the Legendre polynomial 
!                 P_n(x) := 1/(2^n n!) d^n/dx^n (x^2-1)^n
!     using the recursion relation
!               nP_n(x) = (2n-1)xP_(n-1)(x) - (n-1)P_(n-2)(x),
!                        P_1(x) = x,  P_0(x) = 1.
!     The definition of the Legendre polynomials implies that
!           Integral_(-1)^(+1) P_n(x) P_m(x) dx = 2/(2n+1) delta_mn.
!     All calculations are performed with double precision. Note that
!     for |x| > 1 there is no guaranty that the above recurrence
!     relation is stable; the special form of the argument, cos(theta),
!     ensures that the routine works correctly for any parameter theta.
! INPUT PARAMETERS:
!     n,theta: parameters as stated above
! FUNCTION VALUE:
!     value of the "dlegendre" function as defined above
! ----------------------------------------------------------------------

  real*8 function dlegendre (n,theta)
    
    implicit none
    
    real*8  theta, x, old, veryold
    integer n,i
    
! --- COMPUTE FUNCTION VALUE ---

    x = cos(theta)
    if (n .eq. 0) then
       dlegendre = 1.0d0
    elseif (n .eq. 1) then
       dlegendre = x
    elseif (n .ge. 2) then
       old = 1.0d0
       dlegendre = x
       do i = 2,n
          veryold = old
          old = dlegendre
          dlegendre = (dble(2*i-1)*x*old-dble(i-1)*veryold)/dble(i)
       enddo
    endif
    
    return
  end function dlegendre
  
  !-----------------------------------------------------------------------
!                FUNCTION  SPHHARM
!
  ! PURPOSE:
  !   Computes the spherical harmonics  Y_l^m (theta,phi). Here m and l
  !   are integers satisfying 0<=l and -l<=m<=+l
  ! - general spherical harmonics function with Condon-Shortley phase -
  ! INPUT PARAMETERS:
  !   l: integer (see above)
  !   m: integer (see above)
  !   theta,phi : real*8 
  ! FUNCTION VALUE:
  !     value of the "sphharm" function as defined above
  ! CALL TO:
  ! function plgndr (l,m,x=cos(theta)) 
  !-----------------------------------------------------------------------
  complex*16 function sphharm(l,m, theta, phi)

    !
    implicit none
    !
    real*8     theta,phi,pi,ass, fac,x
    integer    i,im,m,l
    !
    
    pi  = 4.d0*atan(1.d0)
    im  = iabs(m)
    x = cos(theta)
    
    fac = 1.d0
    do i = l-im+1, l+im
       fac=fac*i
    enddo

    ass = plgndr(l,im,x)

    if (m.lt.0 .and. mod(im,2).eq.1 ) ass=-ass
    
    sphharm = sqrt((2*l+1)/ (fac*4.d0*pi)) * ass * cdexp(dcmplx(0.d0,m*phi))
    
    return

  end function sphharm


!-----------------------------------------------------------------------
!                SUBROUTINE   ALGNDR
!
!    PURPOSE:
!   Similar to PLGNDR, but stores all P_l^m (x) values in array plm.
!   l runs from m to lmax. 
!   Computes the associated Legendre polynomial P_l^m (x). Here m and l
!   are integers satisfying 0<=m<=l, while x lies in the range
!   -1<=x<=1
! INPUT PARAMETERS:
!   lmax: integer (see above)
!   m: integer (see above),  0<=m<=lmax
!   x: real*8  (see above),  x = cos(theta)
!   switch : logical , if .true. then compute L^2-mormalized functions.
! OUTPUT PARAMETERS:
!     Array plm(k), k=1,..,lmax+1-m.  l = m-1+k.
!
! PHASE !ONVENTION:
!   The phase differs by a factor (-1)**m from most textbooks as
!   e.g. Messiah (Quantum Mechanics) or Zare (Angular Momentum).
!-----------------------------------------------------------------------

  subroutine algndr(lmax,m,x,plm,switch)
    
    implicit none 
    
    logical switch
    integer lmax,m,i,ll
    real*8  x,fact, pll, pmm, pmmp1, somx2
    real*8  plm(*)
    
    if (m.lt.0 .or. m.gt.lmax .or. abs(x).gt.1) then
       write(6,*) 'bad arguments in algndr  - stop -'
       write(2,*) 'bad arguments in algndr  - stop -'
       stop
    end if

! --- Compute P_m^m ---

    pmm=1.
    if (m.gt.0) then
       somx2=sqrt((1.-x)*(1.+x))
       fact=1.
       do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
       enddo
    endif
    plm(1)=pmm
    if (switch)  then
       fact = 1.d0
       do i = 2, 2*m
          fact = fact*i
       end do
       plm(1)= plm(1)*sqrt((m+0.5d0)/fact)
    end if
    if (lmax.eq.m) return
    
    ! --- Compute P_m+1^m ---
    
    pmmp1=x*(2*m+1)*pmm
    plm(2)=pmmp1
    if (switch) then
       fact = fact*(2*m+1)
       plm(2)= plm(2)*sqrt((m+1.5d0)/fact)
    end if
    if (lmax.eq.m+1) return
    
    ! --- Compute P_l^m, l>m+1 ---
    
    do ll=m+2,lmax
       pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
       pmm=pmmp1
       pmmp1=pll
       plm(ll-m+1)=pll
    enddo
    
    if (switch) then
       do ll = m+2,lmax
          fact = fact*(ll+m)/dble(ll-m)
          plm(ll-m+1)=plm(ll-m+1)*sqrt((ll+0.5d0)/fact)
       end do
    end if
    return
  end subroutine algndr
  

end module spherical_harmonics
