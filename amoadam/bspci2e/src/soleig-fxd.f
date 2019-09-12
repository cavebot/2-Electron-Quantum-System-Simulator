
***********************************************************

*
*     call eispack to solve general eigenvalue problem
*
*     A c = eps B c
*
***********************************************************
      subroutine soleig(nm,n,a,b,c,er,fv1,fv2,work,OGP,iflag)
      implicit doubleprecision(a-h,o-z)
      dimension a(nm,1),b(nm,1),c(nm,1),er(nm),work(nm,nm),
     * fv1(1),fv2(1),OGP(1)
      data matz /1/

C    1 FORMAT(2X,'EN(',I2,') & NORMALIZATION =',1P2E19.10)
C1001        format(/'Eigenvalue # ',i2,' =  ',1pd14.6/
C    C       20x,'Eigenvector : '/(1p5d14.6))

*  preserve b-matrix
      do 100 j = 1,n
         do 100 i = 1,n
           work(i,j) = b(i,j)
 100  continue

      call rsg(nm,n,a,work,er,matz,c,fv1,fv2,ierr)

c     call rgg(nm,n,a,work,er,ei,beta,matz,c,ierr)

      DO 10 NB=1,n
      OGP(NB)=0.D0

      DO 10 IR=1,n
         DO 10 IC=1,n
 10         OGP(NB) = OGP(NB) + c(IR,NB) * b(IR,IC) * c(IC,NB)

            if(ierr.ne.0) go to 999

            do  iew = 1, n
cjian             er(iew) = - er(iew)
               WRITE(*,*) iew, er(iew), OGP(iew)
            enddo

            return

 901        write(6,*) ' complex eigenvalue encountered '
            iflag = 999
            return

 999        write(6,*) ' matrix eigenvalue routine fails'
            iflag = ierr
            return

            end
C################################################################
