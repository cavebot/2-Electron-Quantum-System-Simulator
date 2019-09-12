Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      subroutine wfhfin(nwf1e,p,l,no,nhf)
      include "parameter.hr12.inc"
      implicit real*8(a-h,o-z)
      dimension p(np)

      call wf1efile(nwf1e, L) 
      do k = 1, nhf

         read(nwf1e) ( p(j), j = 1, no)
      enddo
      close(nwf1e)

      return
      end
      
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine wfnlin(nwf1e,pp,l,no,nmin,noll)
      implicit real*8(a-h,o-z)
      include "parameter.hr12.inc"
C      parameter(np=2000,ns=102)
      dimension pp(ns,np),p(np)

      call wf1efile(nwf1e, L) 

c      write(*,*)' data file is ',indat(ll)

      kmin=nmin-1

      if(kmin.eq.0) go to 20

      do 15 k=1,kmin
   15 read(nwf1e) (p(j), j=1,no)

   20 do 10 k=1,noll
      read(nwf1e) (p(j), j=1,no)

      do 10 j=1,no
   10 pp(k,j)=p(j)

      close(nwf1e)

      return
      end

Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine wfspin(pp,l,no,nmin,noll,nini,nfin)

      implicit real*8(a-h,o-z)
      include "parameter.hr12.inc"
C      parameter(np=2000,ns=102,kx=15,nl=15)
      parameter(kx=15,nl=15)
      dimension pp(ns,np),p(np)
      dimension t(ns+kx), r(np), dr(np), coe(ns),nini(ns),nfin(ns)
      common/bsp/rs(nl),kb,nb,rmax
      external mkgrid,rin
C........................................

      ll = l + 1


      call mkgrid(nb,kb,rmax,rs(ll),t)

      call rin(r,dr,h,no,idr)

      do i = 1, nb

        coe(i) = 0.0D+00
      end do

      do i = 1, nb
         do j = 1, no

            pp(i,j) = 0.0D+00
         end do
      end do

      do k = 1, noll

        coe(k + nmin - 1) = 1.0D+00

        do j = 1, no
           p(j) = 0.0D+00
        end do

        do j = 1, no

           if( r(j).lt.t( max0(kb, k + nmin - 1 ) ) )  nini(k) = j

           if( r(j).le.t( k + nmin + kb - 1 ) )        nfin(k) = j + 1

        end do


        if(nfin(k).gt.no)  nfin(k) = no

C        write(*,*)'l=',l,'n=',k+nmin-1, 'nini=', nini(k), 'nfin=', nfin(k)

C#
C#      j takes the value zero (j = 0) here, a value that is not allowed for
C#      the fortran array p(j) given that the bounds are (1 - np=2000)
C#

        do j = nini(k), nfin(k) 

           p(j) = bvalue( t, coe, nb, kb, r(j), 0 )

        end do

        do j = nini(k), nfin(k)

           pp(k,j) = p(j)

        end do

        coe( k + nmin - 1 ) = 0.0D+00

      end do
      return
      end
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      subroutine trnsmx(hmx,xm,ihr,ihc,ir,ic)
      implicit real*8(a-h,o-z)
      include "parameter.hr12.inc"
c      parameter(mx=3000,ns=102)
      dimension hmx(ns,mx),xm(ns,ns)
      
C...........................................      
      irf=ihr+ir-1
      icf=ihc+ic-1

      do    k = ihr, irf
         do   kk = ihc, icf

            kr = k - ihr + 1
            kc = kk - ihc + 1

            hmx(k,kk) = xm(kr,kc)

         enddo
      enddo

      return
      end
C##########################################################
C#EOF
