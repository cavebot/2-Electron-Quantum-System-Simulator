!
! modules:param, one_e_matrix
!
!

MODULE param
  !
  USE PRECISION, ONLY:DPK
  !
  IMPLICIT NONE
  PUBLIC
  !
  INTEGER :: ns, nk, nl,np
  INCLUDE "parameter.1e.inc"
  INTEGER ntc, nsp, id
  INTEGER no,nos,nop,nod,isl,nw,lmin,lcore,lmax,itrmx,k
  INTEGER, DIMENSION(nl) :: n2
  REAL(DPK) znuc, rmin, rmax, redmass, crt, di, df
  REAL(DPK), DIMENSION(nl) :: alp1, r01, zk, rs
  REAL(DPK), DIMENSION(nl,np) :: rpw
CONTAINS
  
  SUBROUTINE input
    
    IMPLICIT NONE
    INTEGER i
    
!!!...............
    OPEN(14,file ='inp/d1e.inp',status='old')
    READ(14,*) znuc,rmin,rmax
    READ(14,*) no, redmass
    READ(14,*) nos,nop,nod
    READ(14,*) isl, nw
    READ(14,*) Lmin,Lcore,Lmax,itrmx
    READ(14,*) (alp1(i),i = 1, nw)
    READ(14,*) (r01(i), i = 1, nw)
    READ(14,*) (zk(i),i = 1, nos + nop + nod)
    READ(14,*) (n2(i),i = 1,nw),k
    READ(14,*) (rs(i),i = 1,nw)
    READ(14,*) crt,di,df,id
    CLOSE(14)
    
    nsp = nos + nop
    ntc = nsp + nod
    
    WRITE(*,*) '#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    WRITE(*,   '(a36,1x,I4)' )  ' # input:: Number of          nBsp = ', n2(1)
    WRITE(*,   '(a36,1x,I4)' )  ' # input:: order Bspl.           k = ', k
    WRITE(*,'(a36,1x,G10.3)' )  ' # input:: Box  radius           R = ', rmax
    WRITE(*,'(a36,1x,G10.3)' )  ' # input:: First point         fkn = ', rs(1)    
    WRITE(*,  '(a36,1x,3I4)' )  ' # input:: nof core s,p,d orbitals = ', nos, nop, nod
    IF((isl.NE.1).AND.(itrmx.NE.1)) THEN
       WRITE(*,'(a36,1x,I4)' )  '# input:: Hartree-Fock      lcore = ', lcore
    ELSE 
       WRITE(*,*)'# input:: Hyd - like functions or existing core coeff used'
    ENDIF
    
    IF(redmass.EQ.1) THEN
       WRITE(*,*)'# input:: Hyd - like functions are caclulated'
    ELSE IF(redmass.EQ.2) THEN
       WRITE(*,*)'# input:: Ps - like functions are calculated'
    ELSE 
       WRITE(*,*)'# input:: Not acceptable value for param redmass. Exiting...'
       STOP
    ENDIF
    
    IF((alp1(lmin).EQ.0.D+00)) THEN
       WRITE(*,*)'# input:: No core -polarization '
    ELSE
       WRITE(*,*)'# input:: Core - polarization is included'
    ENDIF
    WRITE(*,*) '#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    
566 FORMAT(5a16)
    
  END SUBROUTINE input
  
END MODULE param


MODULE set_grid

        USE PRECISION, only:dpk
        USE param
        IMPLICIT NONE
        PUBLIC
        INTEGER idbsp, nrp
        REAL(dpk) h, xi
        REAL(dpk), DIMENSION(np) :: r, dr
      
        INTERFACE r_grid
           MODULE PROCEDURE rin
        END INTERFACE
      
        PRIVATE rin, rsin, rexp
      CONTAINS
!!!!###############################################################
        SUBROUTINE rin
          IMPLICIT NONE
!..
          OPEN(10,file='INP/GRID.INP')

          READ(10,*) idbsp  
          SELECT CASE(idbsp)

          CASE(0)
             READ(10,*) xi
             CALL rsin

!             WRITE(*,3) r(2),h,rmax,no
          CASE(1)
             READ(10,*) xi,nrp
             CALL rexp
!             WRITE(*,12) xi,r(1),rmax,h,no,nrp
          END SELECT
          CLOSE(10)

!.................
3         FORMAT(2x,'sine scale:  r(2),h,rmax & no=',1p3e14.6,i5/)
12        FORMAT(2x,'xi,r(1),rmax & h=',1p4e13.5,', no & nrp=',2i5/)
!................
        END SUBROUTINE rin
        !#################################################
        SUBROUTINE rsin
          IMPLICIT NONE
          INTEGER j
          REAL(dpk)  pi, yy, y, drmax, fh, aug
          
!!! r = rmax sin[(pi/2) (x/rmax)**y]
!!! dr = rdx/dr
!!! dr = [2 rmax/(pi*y)]*
!!! [rmax/x]**(y-1) tan[(pi/2) (x/rmax)**y]
!!!
          h      = rmax/dfloat(no-1)
          pi     = 3.14159265358979324D+00
          y      = -dlog( 2.0d+00 * dasin(rmin/rmax) / pi) / dlog(dfloat(no-1) )
          drmax  = 2.0D+00*rmax**y/(pi*y)
          yy     = y - 1.0D+00
          r(1)   = 0.0D+00
          dr(1)  = 0.0D+00
          
          DO  j  = 2, no

             fh    = dfloat(j-1) * h
             aug   = (pi/2.0d0) * (fh/rmax)**y
             r(j)  = rmax  *  dsin(aug)
             dr(j) = drmax * dtan(aug) / fh**yy
             
          END DO
        END SUBROUTINE rsin
      !#######################################
        SUBROUTINE rexp
          USE PRECISION,ONLY:dpk
          IMPLICIT NONE
          INTEGER j
          REAL(dpk)  xmax,rho, x, ri
          !  r = ri*dexp( (x+xi)**(1/nrp))
          ! dr = dx  r/[nrp * (ln(r/ri))**(nrp-1)]
          
          rho = 1.0d0 / dfloat(nrp)
          ri  = rmin
          
          IF(xi.NE.0.d0)  ri = r(1) / dexp((xi)**rho)
          
          xmax = (dlog(rmax/ri))**nrp - xi
          
          h    = xmax/dfloat(no-1)
          
          DO j = 1, no
             
             x    = (j-1) * h
             r(j) = ri * dexp( (x+xi)**rho )
             
          END DO
          !        do j=1,no
          !          dr(j)=1.0d0
          !        end do
          !      if(nrp.eq.1) return
          DO j = 1, no
             
             dr(j) = dfloat(nrp) * (dlog(r(j) / ri))**( nrp - 1 )
             
          END DO
          
        END SUBROUTINE rexp
        !#######################################################################  
        SUBROUTINE cal_rpw
          IMPLICIT NONE
          INTEGER kk, j
!          INCLUDE "PARAMETER.INC"
        INTEGER, PARAMETER :: nl = 15
          !improve

          DO kk = 1, nl
             
             DO j = 1, no
                rpw(kk,j) = r(j)**kk
             END DO
             
          END DO
      
        END SUBROUTINE cal_rpw

      END MODULE set_grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      MODULE one_e_matrix

        IMPLICIT NONE
        PUBLIC
      CONTAINS 
!*******************************************************
!
!  input :
!         n = number of collocation points (no. B-splines)
!         k = order of spline > 2
!         a,b = 1st and last points on interval
!  output :
!         t(1) = t(2) = ... t(k) = a
!         t(n+1) = t(n+2) = ... t(n+k)
!         t(i) = a + (i-k)!h   i = k+1,n
!                                 with  h = (b-a)/(n-k+1)
!
!**********************************************************

        
        SUBROUTINE mkgrid(n,k,rmax,dr0,t)
          USE PRECISION,ONLY:DPK
          IMPLICIT NONE
          INTEGER i, n, k
          REAL(DPK) pi, y, dr0, hh, fh, aug, rmax
          REAL(DPK), INTENT(out), DIMENSION(:) :: t
          
          !........
          
          
          DO i = 1,k             
             t(i) = 0.0_dpk
          END DO
          
          hh = rmax/dfloat(n-k+1)
          pi = 3.14159265358979324d00
          y  = -dlog(2.0d0*dasin(dr0/rmax)/pi)/dlog(dfloat(n-k+1))
!!!  write(*,1000) ' t0 =',dr0,'  h =',hh
1000      FORMAT(/a,1pd14.6,a,1pd14.6/)
          
          !         write(*,*) 'Sine-knots is used for B-slpines'
          DO i = k+1,n             
             fh   =  dfloat(i-k)*hh
             aug  =  (pi/2.0d0)*(fh/rmax)**y
             t(i) =  rmax*dsin(aug)             
          END DO

         DO i = n+1,n+k
             
            t(i) = rmax

         END DO

       END SUBROUTINE mkgrid
     END MODULE one_e_matrix
     !eof
