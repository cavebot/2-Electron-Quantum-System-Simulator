MODULE pad_utils

  USE PRECISION
  IMPLICIT NONE
  PUBLIC  pscoul, sigma_kl
  
CONTAINS
!
!  sigma_l(k) = argG(l+1 - i q/k)
!
  FUNCTION pscoul(l, qz, qk)
    
    USE units, ONLY: M_PI  
    IMPLICIT COMPLEX*16(a-h,r-z), REAL*8(o-q)
    INTEGER L
    
    q1 = l + 1
    ak = CMPLX(0.0D+00, qz/ qk )
    
    ga = q1 - ak
    
    cc=LOG(ga) + LOG(1.d0+ga)+ LOG(2.d0+ga)+LOG(3.d0+ga)
    cc=cc+LOG(4.d0+ga)+ LOG(5.d0+ga)+LOG(6.d0+ga)+LOG(7.d0+ga)
    
    z = ga + 8.0D+00
    
    gam = -cc + (z - 0.5d0)* LOG(z)-z+0.5d0* LOG(2.d0*M_PI)+1.d0/(12.d0*z)
    gam = gam-1.d0/(360.d0*z**3)+1.d0/(1260.d0*z**5)-1.d0/(1680.d0*z**7)
    pscoul= AIMAG(gam)
    
    RETURN
  END FUNCTION pscoul
  !
  ! f90 version
  !
  REAL(dpk) FUNCTION sigma_kl(l, q, k)
    !
    USE PRECISION, ONLY: DPK
    USE units, ONLY: M_PI  
    IMPLICIT NONE
    !
    INTEGER,   INTENT(in) :: l
    REAL(dpk), INTENT(in) :: q
    REAL(dpk), INTENT(in) :: k
    !locals
    DOUBLE COMPLEX :: gam, ak, ga, z, cc
    INTEGER i
    INTEGER q1
    
  ! executable statements
    
    !q1 = DBLE(l + 1)
    !ak = CMPLX(0.0D+00, q/k )
    !ga = q1 - ak
    
    ga = CMPLX( DBLE(l+1), -q/k )


    
    DO i = 0, 7
       cc = cc + LOG( DBLE(i)+ ga )
    ENDDO
    
    z = ga + 8.0D+00
  
    gam = - cc + (z - 0.5D+00) * LOG(z) - z + LOG( 2.0D+00 *M_PI )/2.0D+00 + 1.0D+00/( 12 * z )
    
    gam = gam - 1.d0/(360*z**3) + 1.d0/ (1260*z**5) - 1.d0/(1680* z**7)
    
    sigma_kl = AIMAG( gam )
  
!    RETURN sigma_kl
  END FUNCTION sigma_kl

! function gammln. Calculates the log(gamma(x))
! It is a conversion to complex argument version of 
! the Numerical Recepies by William. H. Press et al.
!
  FUNCTION GAMMLN(XX)
    !
    DOUBLE PRECISION COF(6),STP,HALF,FPF
    INTEGER J
    DOUBLE COMPLEX x,xx,gammln,tmp,ser,one
    DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,-1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/ 
    DATA HALF,ONE,FPF/0.5D0,(1.0D0,0.d0),5.5D0/
    !executables

    !
    X = XX - ONE
    TMP = X + FPF
    TMP = ( X + HALF) * LOG(TMP) - TMP
    SER = ONE
    !
    DO J = 1, 6
       X = X + ONE
       SER = SER + COF(J) / X
    ENDDO
    
    GAMMLN = TMP + LOG(STP*SER)
    RETURN

  END FUNCTION GAMMLN

  END MODULE pad_utils

!      ion = ion + aimag( sum( conjg(f(1,:)) * df(:) )  ) * dt
!
!  this function caculate arg of gamm function
!


!eof!
