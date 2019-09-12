! Library of angular momentum coupling coefficient routines in fortran 90
! Paul Stevenson, Oxford University/Oak Ridge National Laboratory.
! spaul@mail.phy.ornl.gov
!
!  contains: c3j,c6j,c9j, factorial(n),binomial(n,m)
!
module ang

  INTEGER, PARAMETER :: dpk = KIND(1.0D+00) 
!  INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(p=15)

contains


  
!!!c3j!   
!!!
!!! based on Messiah relation, appendix:vector addition coefficients and rotation matrices
!!!
!!! (j1, j2, j3 ; m1, m2, -m) = (-)^(j1-j2+m) * < j1 m1 j2 m2 | j m > / (2j+1)^(1/2)
!!!
  FUNCTION c3j(j1,m1,j2,m2,j,m)
    IMPLICIT NONE
    INTEGER    :: j1,m1,j2,m2,j,m 
    REAL(dpk)  :: c3j
    REAL(dpk)  :: phase

    phase=(-1)**(j1-j2+m)
       
    c3j = phase * cleb(j1,m1,j2,m2,j,m)/SQRT(DBLE(j)+1)
 
  END FUNCTION c3j

!!! calculate a clebsch-gordan coefficient < j1/2 m1/2 j2/2 m2/2 | j/2 m/2 >
!!! arguments are integer and twice the true value. 
!!!clebsch-Gordan!
  FUNCTION cleb(j1,m1,j2,m2,j,m)
    implicit none
    real(dpk)    :: cleb,factor,sum
    integer :: j1,m1,j2,m2,j,m,par,z,zmin,zmax
    !
    ! some checks for validity (let's just return zero for bogus arguments)
    !
    if (      2*(j1/2)-int(2*(j1/2.0)) /= 2*(abs(m1)/2)-int(2*(abs(m1)/2.0))     &
         .or. 2*(j2/2)-int(2*(j2/2.0)) /= 2*(abs(m2)/2)-int(2*(abs(m2)/2.0))     &
         .or. 2*( j/2)-int(2*( j/2.0)) /= 2*(abs( m)/2)-int(2*( abs(m)/2.0))     &
         .or.      j1 < 0  .or.         j2 < 0  .or.      j < 0                  &
         .or. abs(m1) > j1 .or.    abs(m2) > j2 .or. abs(m) > j                  &
         .or.  j1+j2  < j  .or. abs(j1-j2) > j  .or. m1+m2 /= m )                &
         then

       cleb= 0.0
    !
    else
    !
       factor = 0.0
       factor = binom(j1,(j1+j2-j)/2) / binom((j1+j2+j+2)/2,(j1+j2-j)/2)
       factor = factor * binom(j2,(j1+j2-j)/2) / binom(j1,(j1-m1)/2)
       factor = factor / binom(j2,(j2-m2)/2) / binom(j,(j-m)/2)
       factor = sqrt(factor)
       
       zmin = max(0,j2+(j1-m1)/2-(j1+j2+j)/2,j1+(j2+m2)/2-(j1+j2+j)/2)
       zmax = min((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2)
       
       sum=0.0
       do z = zmin, zmax
          par = 1

          if( 2*(z/2) - int(2*(z/2.0) ) /= 0) par=-1

          sum = sum + par * binom( ( j1 + j2 - j)/2,        z    )    &
                          * binom( ( j1 - j2 + j)/2, (j1-m1)/2-z )    &
                          * binom( (-j1 + j2 + j)/2, (j2+m2)/2-z )
       end do
       
       cleb = factor*sum/sqrt(dble(j)+1)
 
    end if

  END FUNCTION cleb
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  c6j
    ! calculates a Wigner 6-j symbol. Argument a-f are integer and are
    ! twice the true value of the 6-j's arguments, in the form
    ! { a b c }
    ! { d e f }
    ! Calculated using binomial coefficients to allow for (reasonably) high
    ! arguments.
!c6j!
  function c6j(a,b,c,d,e,f)
    implicit none
    integer, intent(in) :: a,b,c,d,e,f
    real(dpk) :: c6j
    integer :: phase, nlo, nhi, n
    real(dpk) :: outfactors, sum, sumterm

    ! First check for consistency of arguments:
    c6j=0.0
    if(mod(a+b,2)/=mod(c,2)) return
    if(mod(c+d,2)/=mod(e,2)) return
    if(mod(a+e,2)/=mod(f,2)) return
    if(mod(b+d,2)/=mod(f,2)) return
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(c-d)>e .or. c+d<e) return
    if(abs(a-e)>f .or. a+e<f) return
    if(abs(b-d)>f .or. b+d<f) return

    phase=(-1)**((a+c+d+f)/2)

    outfactors = angdelta(a,e,f)/angdelta(a,b,c)
    outfactors = outfactors * angdelta(b,d,f)*angdelta(c,d,e)

!    write(6,*) outfactors

    nlo = max( (a+b+c)/2, (c+d+e)/2, (b+d+f)/2, (a+e+f)/2 )
    nhi = min( (a+b+d+e)/2, (b+c+e+f)/2, (a+c+d+f)/2)

    sum=0.0
    do n=nlo,nhi
       sumterm = (-1)**n
       sumterm = sumterm * binom(n+1,n-(a+b+c)/2)
       sumterm = sumterm * binom((a+b-c)/2,n-(c+d+e)/2)
       sumterm = sumterm * binom((a-b+c)/2,n-(b+d+f)/2)
       sumterm = sumterm * binom((b-a+c)/2,n-(a+e+f)/2)
!       write(6,*) ',sumterm: ',sumterm
       sum=sum+sumterm
    end do

    c6j = phase * sum * outfactors

  end function c6j
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx c9j
    ! calculate a 9-j symbol. The arguments are given as integers twice the
    ! value of the true arguments in the form
    ! { a b c }
    ! { d e f }
    ! { g h i }
!c9j!
  function c9j(a,b,c,d,e,f,g,h,i)
    implicit none
    integer   :: a,b,c,d,e,f,g,h,i
    real(dpk) :: c9j, sum
    integer   :: xlo, xhi
    integer   :: x

    c9j=0.0

    ! first check for bogus arguments (and return zero if so)
    if(abs(a-b)>c .or. a+b<c) return
    if(abs(d-e)>f .or. d+e<f) return
    if(abs(g-h)>i .or. g+h<i) return
    if(abs(a-d)>g .or. a+d<g) return
    if(abs(b-e)>h .or. b+e<h) return
    if(abs(c-f)>i .or. c+f<i) return
    
    xlo = max(abs(b-f),abs(a-i),abs(h-d))
    xhi = min(b+f,a+i,h+d)
    
    sum=0.0
    do x=xlo,xhi,2
       sum=sum+(-1)**x*(x+1)*c6j(a,b,c,f,i,x)*c6j(d,e,f,b,x,h)*&
            c6j(g,h,i,x,a,d)
    end do
    c9j=sum

  end function c9j
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    G(n) =  n!
!G(n)!
  recursive function factorial(n) result(res)
    implicit none
    integer :: n
    real(dpk) :: res

    if (n==0 .or. n==1) then
       res=1.0
    else
       res=n*factorial(n-1)
    end if
  end function factorial
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    b(n m) 
!b(n,m)!
  recursive function binom(n,r) result(res)
    implicit none
    integer :: n,r
    real(dpk) :: res

    if(n==r .or. r==0) then
       res = 1.0
    else if (r==1) then
       res = real(n,dpk)
    else
       res = real(n,dpk)/real(n-r,dpk)*binom(n-1,r)
    end if
  end function binom
  ! calculate the function delta as defined in varshalovich et al. for
  ! use in 6-j symbol:
  !
  function angdelta(a,b,c)
    implicit none
    INTEGER  :: a,b,c
    real(dpk) :: angdelta, scr1
    scr1= factorial((a+b-c)/2)
    scr1=scr1/factorial((a+b+c)/2+1)
    scr1=scr1*factorial((a-b+c)/2)
    scr1=scr1*factorial((-a+b+c)/2)
    angdelta=sqrt(scr1)
  end function angdelta
end module ang

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!eof

