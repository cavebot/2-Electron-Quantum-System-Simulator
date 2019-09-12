!f2py double precision intent(in):: j1,m1,j2,m2,j3,m3
!f2py double precision intent(out):: clebsch

!Function to calculate Clebsch-Gordan coefficients,
!Original code by  David G. Simpson,
!NASA Goddard Space Flight Center


MODULE ang_utils

  USE PRECISION,ONLY: dpk
  !  integer, parameter :: rk = selected_real_kind(p=15)

CONTAINS

subroutine clebschgordan(j1,m1,j2,m2,j3,m3,clebsch)
implicit none
double precision:: j1,j2,j3,m1,m2,m3
double precision:: c, clebsch, fact4, fact5, fact6
double precision:: term, sumk, fact1, fact2, fact3
logical::          frac1,frac2,frac3,frac4,frac5,frac6
integer:: k

  call isfrac(j1+j2+j3,frac1)
  call isfrac(j1+m1,frac2)
  call isfrac(j2+m2,frac3)
  call isfrac(j3+m3,frac4)
  call isfrac(-j1+j3-m2,frac5)
  call isfrac(-j2+j3+m1,frac6)
  if (.not. frac1 .or.  .not.frac2 .or..not. frac3 .or..not.   &
      frac4  .or. .not. frac5 .or. .not. frac6) then
    c=0.0d0
  end if

!  Check for conditions that give C = 0.

!write(*,*)j1,j2,j3
  if ( (j3 .lt. abs(j1-j2)) .or.  &
       (j3 .gt. abs(j1+j2))    .or.  &
       (abs(m1) .gt. j1)    .or.  &
       (abs(m2) .gt. j2)    .or.  &
       (abs(m3) .gt. j3) .or. (m1+m2) .ne. m3) then
     c = 0.0d0
  else
!Compute clebsch-gordan coefficient.
     call factorial(j1+j2+j3+1,fact1)
     c = sqrt((j3+j3+1)/fact1)
     call factorial(j1+j2-j3,fact1)
     call factorial(j2+j3-j1,fact2)
     call factorial(j3+j1-j2,fact3)
     c = c * sqrt(fact1*fact2*fact3)
     call factorial(j1+m1,fact1)
     call factorial(j1-m1,fact2)
     call factorial(j2+m2,fact3)
     call factorial(j2-m2,fact4)
     call factorial(j3+m3,fact5)
     call factorial(j3-m3,fact6)
     c = c * sqrt(fact1*fact2*fact3*fact4*fact5*fact6)
     sumk = 0.0d0
     do k = 0, 99
        if (j1+j2-j3-dble(k) .lt. 0d0) cycle
        if (j3-j1-m2+dble(k) .lt. 0d0) cycle
        if (j3-j2+m1+dble(k) .lt. 0d0) cycle
        if (j1-m1-dble(k)    .lt. 0d0) cycle
        if (j2+m2-dble(k)    .lt. 0d0) cycle
        call factorial(j1+j2-j3-dble(k),fact1)
        call factorial(j3-j1-m2+dble(k),fact2)
        call factorial(j3-j2+m1+dble(k),fact3)
        call factorial(j1-m1-dble(k),fact4)
        call factorial(j2+m2-dble(k),fact5)
        call factorial(dble(k),fact6)
        term = fact1*fact2*fact3*fact4*fact5*fact6
        if (mod(k,2) .eq. 1) term = -term
        sumk = sumk + 1.0d0/term
     end do
     c = c * sumk
  end if
  clebsch= c
  return
end subroutine
! !-----------------------------------------------------------------------------
!f2py integer intent(in):: a,b,c,d,e,f
!f2py double precision intent(out) :: threej
subroutine threejsymbol(a,b,c,d,e,f, threej)
implicit none
double precision:: a, b, c, d, e, f
double precision:: threej, clebsch
logical::          cond

call clebschgordan(a,d,b,e,c,-f,clebsch)
threej= clebsch*((-1d0)**(b-a+f))/sqrt(2d0*c + 1d0)
call condition(a,b,c,cond)
if (cond) then
   threej=0d0
end if
return
end subroutine threejsymbol
!-----------------------------------------------------------------------------
subroutine factorial(x,fact)
!f2py double precision intent(in) :: x
!f2py double precision intent(out) :: fact
implicit none
double precision:: x, fact, y
logical::          frac

y=x-1
fact=x
do
   if (x .lt. 0) then
      fact=0d0
      exit
   else if (x .eq. 0) then
      fact=1d0
      y=-1
   end if
  if (y .le. 0) then
    exit
  end if
  fact=fact*y
  y=y-1
end do
call isfrac(x,frac)
if (.not. frac) then
  write(*,*)'not a fraction'
  fact=0d0
end if
return
end subroutine
!-----------------------------------------------------------------------------
subroutine triangle(x,y,z, triangleanswer)
!f2py double precision intent(in) :: x,y,z
!f2py double precision intent(out) :: triangleanswer
implicit none
double precision:: x,y,z, triangleanswer
double precision:: fact1,fact2,fact3,fact4
!write(*,*)x,y,z
call factorial(x+y-z,fact1)
call factorial(x-y+z,fact2)
call factorial(y+z-x,fact3)
call factorial(x+y+z+1,fact4)
triangleanswer= fact1*fact2*fact3/fact4
return
end subroutine
!-----------------------------------------------------------------------------
subroutine isfrac(x,answer)
!f2py double precision intent(in):: x
!f2py logical intent(out):: answer
implicit none
double precision:: x
logical:: answer

answer= .false.
if ((abs(2d0*x)-abs(int(2d0*x))) .lt. 1d-8) then
  answer = .true.
end if
return
end subroutine
!-----------------------------------------------------------------------------
subroutine sixjsymbol(j1,j2,j3,m1,m2,m3,sixj)
!f2py double precision intent(in):: j1,j2,j3,m1,m2,m3
!f2py double precision intent(out):: sixj
implicit none
double precision:: j1,j2,j3,m1,m2,m3
double precision:: a,b,c,d,e,f,g,x
double precision:: sumx, sixj,sqrtcoeff
double precision:: fact1,fact2,fact3
double precision:: tri1,tri2
integer:: t
logical:: cond1, cond2, cond3, cond4
sumx= 0d0

!write(*,*)j1,j2,j3,m1,m2,m3
cond1= .false.
cond2=.false.
call condition(j1,j2,j3,cond1)
call condition(m1,j2,m3,cond2)
!!$call condition(j1,m2,m3,cond3)
!!$call condition(m1,m2,j3,cond4)
if (cond1 .or. cond2) then
   sixj=0d0
   !write(*,*)"triangle"
else
   do t= -500,500,1
      a= t-j1-j2-j3
      b= t-j1-m2-m3
      c= t-m1-j2-m3
      d= t-m1-m2-j3
      e= j1+j2+m1+m2-t
      f= j2+j3+m2+m3-t
      g= j3+j1+m3+m1-t
      if (a .lt. 0d0 .or. b .lt. 0d0 .or. c .lt. 0d0 .or. d .lt. 0d0 &
           .or. e .lt. 0d0 .or. f .lt. 0d0 .or. g .lt. 0d0) cycle
      call factorial(a,fact1)
      call factorial(b,fact2)
      call factorial(c,fact3)
      x= fact1*fact2*fact3
      call factorial(d,fact1)
      call factorial(e,fact2)
      x= x*fact1*fact2
      call factorial(f,fact1)
      call factorial(g,fact2)
      x= x*fact1*fact2
      call factorial(dble(t)+1d0,fact3)
      x= ((-1)**t)*fact3/x
      sumx= sumx+x
   end do

   call triangle(j1,j2,j3,tri1)
   call triangle(j1,m2,m3,tri2)
   sqrtcoeff= tri1*tri2
   call triangle(m1,j2,m3,tri1)
   call triangle(m1,m2,j3,tri2)
   sqrtcoeff= sqrtcoeff*tri1*tri2
   sixj= sqrt(sqrtcoeff)*sumx
   !write(*,*)sumx,sqrtcoeff
end if
return
end subroutine
!-----------------------------------------------------------------------------
subroutine binomial(n,k,bin)
!f2py integer intent(in):: n,k
!f2py double precision intent(out):: bin
implicit none

integer::          n,k
double precision:: bin, fact1, fact2, fact3

call factorial(dble(n), fact1)
call factorial(dble(k), fact2)
call factorial(dble(n-k), fact3)
bin= fact1/(fact2*fact3)

return
end subroutine
!-----------------------------------------------------------------------------
subroutine condition(x,y,z,cond)
!f2py double precision intent(in):: x,y,z
!f2py logical intent(out):: cond
implicit none

double precision:: x,y,z 
logical:: cond
!if cond is true, then triangular inequalities not satisfied, should skip state
cond=.false.
if ((x+y) .lt. z .or. abs(x-y) .gt. z) then
   cond= .true.
   !write(*,*)x,y,z,1
end if

if ((x+z) .lt. y .or. abs(x-z) .gt. y) then
   cond= .true.
   !write(*,*)x,y,z,2
end if

if ((z+y) .lt. x .or. abs(z-y) .gt. x) then
   cond= .true.
   !write(*,*)x,y,z,3
end if
return
end subroutine

END MODULE ang_utils
