C***
C***  program to calculate the h-matrix for ca-like atoms
C***  bsp hf wf are used for all configurations (normal parity)

      program r12b
      implicit real*8(a-h,o-z)
      include "parameter.hr12.inc"
      dimension pr1(np),pc1(np)
      dimension vxs(ns,ns), ehf(10,ns), nhf(ncs),lhf(ncs),ll(ncs)
      dimension nllmin(ncs),nllmax(ncs),noll(ncs),is(ncs),pr(ns,np)
      dimension ndi(ncs), idcs(ncs),pc(ns,np),yk(np)
      dimension n2i(ns),n2f(ns),n4i(ns),n4f(ns)
      common/hx/hmxs(ns,mx)
      common/rarr/rpw(12,np)
      common/bsp/rs(15),kb,nb,rmax
      common/rro/rho(6),alpha, beta
      common/rh/h,r(np),dr(np)
      CHARACTER*100 ARGV
C...................................



      NR12 = 3
      NWF1E = 2
      NCFG = 15

C.....


      CALL GETARG(1, ARGV)
      READ(ARGV,*) INPL         

      WRITE(*,*) 'hr12::             partial wave   L =', INPL


C~
C#    alpha :: dielectronic polarization potential: Vd
C#    beta  :: mass polarization potential        : Vm
C#

      OPEN(9,FILE='inp/r12.inp')
      READ(9,*) alpha, beta
      READ(9,*) ( rho(i), i = 1, lf)
      READ(9,*)  kb, nb, rmax
      READ(9,*) ( rs(i), i = 1, lf)
      CLOSE(9)


      OPEN(16,FILE='out/r12b.out')


C# read 1e target states


      OPEN(15,FILE='dat/en1e.dat')

      READ(15, 3) LLMIN,LLMAX

      DO  LP = LLMIN, LLMAX

         DO  IE = 1, NS
            EHF(LP,IE) = 0.0D+00
         ENDDO 

         READ(15,3) NMX, NVALENCE

         DO  IE = 1, NMX

            READ(15,'(E20.14)') EHF(LP,IE)

         ENDDO      
      ENDDO

      CLOSE(15)

C# endread 1e-target states


        WRITE(16,*) ' RUNNING        HR12.F         PROGRAM :  (2) '
        WRITE(16,*) ' '
        WRITE(16,*) ' PURPOSE: '
        WRITE(16,*) ' STORING 2-e HAMILTONIAN MATRIX ON MIXED BASIS 
     1  AT DAT/HR12-DAT FILES '
        WRITE(16,*) ' '
        WRITE(16,*) ' H(i,j) = < B(i) | H(1) + H(2) 
     1  + Vr(1,2) + Vd(1,2) | B(j) > '
        WRITE(16,*) ' '
        write(16,*) ' Order of  B - splines,  k  = ',  kb 
        write(16,*) ' Number of B - splines,  n  = ',  nb
        write(16,*) ' Box radius in a. units  r  = ',  rmax
        write(16,*) ' First knot point for l: rs = ', ( rs(i), i=1,lf)


C.... make grid for evaluation of integrals
 

      call rin(r, dr, h, no, idr)


C.. construct 1/r^l matrix


      do  kk=1,12
         do  j=1,no
            rpw(kk,j) = r(j)**kk
         enddo
      enddo

      
C... read configuration files

      lo = inpl

         call cxfin(ncfg, nhf,lhf,ll,nllmin,nllmax,noll,
     1        is,lo,ls, ncsmx, ncmx, ndi, idcs)


C.... ls=1(singlet), 3(triplet)


      if(inpl.ne.lo) then
         write(16,*) ' CFG-L.INP FILE : INCOSISTENT L. Exiting'
         stop
      endif


      CALL hfile(nr12,"dat","r12b","bin",lo)

C      call hr12file(NR12, LO)


 709  kr = 1

      do 100 ir = 1, ncsmx


         call wfhfin(nwf1e, pr1,lhf(ir),no,nhf(ir))

         if( idcs(ir).eq.0) then     ! free-bc 

            call wfnlin(nwf1e, pr,ll(ir),no,nllmin(ir),noll(ir))

            do i = 1, noll(ir)
               n2i(i) = 1
               n2f(i) = no
            end do

         else                        ! fixed-bc (correlation functions)

            call wfspin(pr,ll(ir),no,nllmin(ir),noll(ir),n2i,n2f)
         end if

C..

         kc = kr

!nullify
         do  i = 1, ns
           do  k = 1, mx
             hmxs(i,k) = 0.d0
           enddo
         enddo

C...  

               do 101 ic = ir, ncsmx

                 call wfhfin(nwf1e,pc1,lhf(ic),no,nhf(ic))
                 
                 if(idcs(ic).eq.0) then
                    
                    call wfnlin(nwf1e,pc,ll(ic),no,nllmin(ic),noll(ic))
                    
                    do i=1,noll(ic)
                      n4i(i)=1
                        n4f(i)=no
                     end do

                  else

                  call wfspin(pc,ll(ic),no,nllmin(ic),noll(ic),n4i,n4f)
               end if

C..     calculate the <pr|1/r_12|pc>

               
               call submx(pr1,pr,pc1,pc,r,dr,vxs,nhf,nllmin,noll,lhf,
     1                 ll,ehf,ir,ic,lo,ls,no,h,idr,n2i,n2f,n4i,n4f)


C# equivalent orbitals

      if(idcs(ic).eq.0.and.idcs(ir).eq.0) then 

         do nr=1,noll(ir)
            do nc=1,noll(ic)
               call q2tst1(lhf,ll,nhf,nllmin,ir,ic,nr,nc,vxs(nr,nc))
            end do
         end do

      end if
C
      if(idcs(ir).eq.1.and.idcs(ic).eq.0) then

         do nr=1,noll(ir)
            do nc=1,noll(ic)

               call q2tst0(lhf,ll,nhf,nllmin,ir,ic,nr,nc,vxs(nr,nc))

            end do
         end do

      end if

C................


C     call mxprnt(vxs,noll(ir),noll(ic))

      
      call trnsmx(hmxs,vxs,1,kc,noll(ir),noll(ic))


      kc = kc + noll(ic)

      write(*,*) '# r12b: < ir|V_11|ic> = ',  vxs(1, 1), ir,ic

c      write(*,*)ir, 'vxs(1,1)=',vxs(1,1)

 101  continue

      kr = kr - 1


! store V_12 configuration interaction matrix

      do  nr = 1, noll(ir)
        kr = kr + 1
        write(3) ( hmxs(nr,nc), nc = kr, ncmx)
      enddo


      kr = kr + 1

 100  continue


      close(3)
      close(4)
      close(16)
 
C........

    1 format(4d15.7)
    2 format(2x,1p9d14.6)
    3 format(6i5)
    4 format(1pd16.8)
    5 format(/2x,'row -- nl & ll=',3i3,',  column -- nl & ll=',3i3)
    6 format(2x,'l & ncmx=',2i4/)
    7 format(5a12)
C........

      end
c*******************************************************************
      subroutine submx(p1,pr,p3,pc,r,dr,vs,nhf,nllmin,noll,lhf,ll,ehf,
     1 ir,ic,lo,ls,no,h,idr,n2i,n2f,n4i,n4f)
      implicit real*8(a-h,o-z)
      include "parameter.hr12.inc"
C      parameter(np=2000,ns=102,mx=4000,ncs=50)
      dimension p1(np),pr(ns,np),p3(np),pc(ns,np),r(np),dr(np),
     1    vs(ns,ns), p2(np),p4(np), ehf(10,ns),
     2    nhf(ncs),noll(ncs),lhf(ncs),ll(ncs),nllmin(ncs),
     3    angd(10),kad(10),ange(10),kae(10),yd(np),ye(ns,np),yk(np)
 
      dimension n2i(ns),n2f(ns),n4i(ns),n4f(ns)
      common/rarr/rpw(12,np)
      common/rro/rho(6),alpha, beta
C..................................



      call angfk( lhf(ir),ll(ir),lhf(ic),ll(ic),lo,angd,kad,knod)




      call angfk( lhf(ir),ll(ir),ll(ic),lhf(ic),lo,ange,kae,knoe)


C.................................



      pha = (-1)**(iabs( lo + lhf(ic) - ll(ic) ) )



C.................................

      do  nr = 1, ns
         do  nc = 1, ns

            vs(nr,nc) = 0.0D+00

         end do
      end do 


      do j=1,no

         yd(j) = 0.0D+00
      end do


      do  k = 1, knod

         call set_yk(p1,p3,yk,r,dr,kad(k),no,h,idr,1,no)


         do  j = 1, no
        
            yd(j) = yd(j) + angd(k) * yk(j)
         enddo

      enddo


      do n = 1, noll(ic)

         do  j = 1, no

            ye(n,j) = 0.0D+00
            p4(j)   = pc(n,j)

         enddo


         do  k = 1, knoe

            call set_yk(p1,p4,yk,r,dr,kae(k),no,h,idr,n4i(n),n4f(n))

            do  j = 1, no

               ye(n,j) = ye(n,j) + ange(k) * yk(j)


               enddo
            enddo

            
       end do
       


      do 40 nr = 1, noll(ir)



         do  j = n2i(nr), n2f(nr)

            p2(j) = pr(nr,j)

         enddo


         do 40 nc = 1, noll(ic)


C            write(*,*) ' 1 === (nr,nc) = ', nr,nc
 
           do  j = n4i(nc), n4f(nc)

               p4(j) = pc(nc,j)
            enddo


            yk(1) = 0.0D+00

            if(dr(1).ne.0.d0) yk(1) = yd(1)*p2(1)*p4(1)/dr(1)

            nomin = max0( n4i(nc), n2i(nr),  2)

            nomax = min0( n4f(nc), n2f(nr), no)


            kj = 1
            if(nomin.gt.2) yk(kj) = 0.0D+00

            do j = nomin, nomax

               kj = kj + 1
               yk(kj) = yd(j) * p2(j) * p4(j) / dr(j)
            enddo


            if(nomin.le.nomax) then

               noint = kj
               vd = rint(yk, 1, noint, 14, h)
               vd = 2.D+00 * vd
            else

               vd = 0.0D+00

            end if

            yk(1) = 0.0D+00
            if(dr(1).ne.0.d0) yk(1) = ye(nc,1) * p2(1) * p3(1) / dr(1)

            kj = 1

            if(max0(n2i(nr),2).gt.2)  yk(kj) = 0.0D+00

            if(max0(n2i(nr),2).le.n2f(nr)) then

               do 26 j = max0(n2i(nr),2),n2f(nr)
                  kj = kj + 1
                  yk(kj) = ye(nc,j) * p2(j) * p3(j)/dr(j)
 26            continue


               no2 = kj

c     no2=n2f(nr)-max0(n2i(nr),2)+1

               ve = rint(yk, 1, no2, 14, h)
               ve = pha * 2.d0 * ve

            else

               ve = 0.0D+00
            end if

C            write(*,*) ' 2 === (nr,nc) = ', nr,nc
C###### here calculate the dielec. core polar. 
c	write(*,*) ' alph = ', alph

      vdd = 0.0D+00
      vde = 0.0D+00

C##  rho(0) do not exist

         rho1 = rho(lhf(ir))
         rho2 = rho(ll(ir))

      do 55 k = 1, knod

         if(kad(k) .ne. 1) goto 55

         vdd = vdd - 2.0D+00 * alpha * angd(k)
     1        * vdr( p1, lhf(ir), p3, lhf(ic),rho1, no)
     1        * vdr( p2, ll(ir),  p4, ll(ic), rho2, no)

 55   continue


      do 56 k = 1, knoe

         if(kae(k) .ne. 1) goto 56

         vde = vde - 2.0 * alpha * pha * ange(k)
     1             * vdr( p1, lhf(ir), p4, ll(ic),  rho1, no)
     1             * vdr( p2, ll(ir),  p3, lhf(ic), rho2, no)
c      print*,'vde=',vde
 56   continue

C#          write(*,*) ' 3 === (nr,nc) = ', nr,nc
C#          write(*,*) ' ===================================='
C########   MASS POLARIZATION TERM (POSITTONIUM LIKE SYSTEMS)

      beta = 0.0D+00

      VMD = 0.0D+00
      VME = 0.0D+00


c$$$      LHF1 = LHF(IR) + 1
c$$$      LHF3 = LHF(IC) + 1
c$$$      LL2  = LL(IR)  + 1
c$$$      LL4  = LL(IC)  + 1
c$$$
c$$$C*   'direct term'
c$$$
c$$$      DE13 = 0.5D+00 * ( EHF( LHF1, NHF(IR) ) - EHF( LHF3, NHF(IC) ))
c$$$      DE24 = 0.5D+00 * ( EHF( LL2, NR ) -  EHF(LL4, NC))
c$$$
c$$$      do 57 k = 1, knod
c$$$
c$$$         if(kad(k).ne.1) goto 57
c$$$
c$$$         VMD = VMD + 0.3333D+00 * beta * angd(k)
c$$$     1                       * DE13 * VMR( P1, P3, NO )
c$$$     1                       * DE24 * VMR( P2, P4, NO )
c$$$
c$$$c         write(*,*) 'vmd = ', vmd, 'xxxxxxxxx'
c$$$ 57   continue
c$$$
c$$$C*   'exchange term'
c$$$
c$$$      DE14 = 0.5D+00 * ( EHF( LHF1, NHF(IR) ) - EHF( LL4, NC ) )
c$$$      DE23 = 0.5D+00 * ( EHF( LL2, NR ) -  EHF(LHF3, NHF(IC)) )
c$$$
c$$$      do 58 k = 1, knoe
c$$$
c$$$         if(kae(k).ne.1) goto 58
c$$$
c$$$         VME  =  VME + 0.333D+00 * beta * PHA  * ange(k)      
c$$$     1                         * DE14 * VMR( P1, P4, NO)
c$$$     1                         * DE23 * VMR( P2, P3, NO)
c$$$
c$$$c          write(*,*) 'vmd = ', vme, 'xxxxxxxxx'
c$$$ 58   continue



           if(ls.eq.1) then
               
               vs(nr,nc) = vd + ve + vdd + vde + vmd + vme

            else if(ls.eq.3) then

               vs(nr,nc) = vd - ve + vdd - vde + vmd - vme

            endif

C#
C#            vs(nr,nc)=vd+ve+vdd+vde
C#            vt(nr,nc)=vd-ve+vdd-vde
C#

 40      continue

      return
      end

C###########################################################
      subroutine set_yk(p1,p2,yk,r,dr,k,no,h,idr,nmin,nmax)

      implicit real*8(a-h,o-z)
      include "parameter.hr12.inc"
C      parameter(np=2000)
      common/rarr/rpw(12,np)
      dimension zk(np),wy(np),wz(np),drkr(np),rkk(np)
      dimension p1(1),p2(1),yk(1),r(1),dr(1)
C.................................

      n  = 1
      kr = k + n - 1
      kk = k + 1

      if(kr.eq.0) go to 30

      drkr(1)=0.0d0

      if(r(1) .ne. 0.0d0) drkr(1) = 1.0D+00/rpw(kr,1)

      do 20 j = 2, no
   20 drkr(j) = 1.0D+00/rpw(kr,j)

   30 continue

      do 21 j=1,no
   21 rkk(j)=rpw(kk,j)

      call ykfctb(p1,p2,yk,zk,wy,wz,r,dr,drkr,rkk,h,k,no,idr,nmin,nmax)

      return
      end
C#################################################
      function vmr(PR, PC, NO)
      IMPLICIT REAL*8(a-h,o-z)
      INCLUDE "parameter.hr12.inc"
C      PARAMETER(NP = 2000 )
      DIMENSION PR(NP),PC(NP),F(NP)
      common/RH/H, R(NP), DR(NP)
C.................................................

      F(1) = 0.0D+00

      DO j = 2, NO

         F(J) = PR(J) * PC(J) * R(J)**2 / DR(j)

      END DO

      VMR = RINT(F, 1, NO, 14, H)

      RETURN
      END
C################################################
      function vdr(pr,lr,pc,lc,rho,no)
      implicit real*8(a-h,o-z)
      include "parameter.hr12.inc"
C      parameter(np=2000,ns=102,mx=4000,ncs=50)
      dimension pr(np),pc(np),vp(np),f(np)
      common/rh/h,r(np),dr(np)

C................................


      call vpol(r,vp,rho,h,no)
      
      f(1) = 0.0d0

      do j=2,no


         if(dr(j).eq.0.0D+00) then

            write(*,*) ' vdr::Division with zero dr(',j,' ) = ', dr(j) 
            stop
         endif


         f(j) = pr(j) * pc(j) * vp(j) / dr(j)

      end do


      vdr = rint(f,1,no,14,h)

c      print*,'vdr=',vdr
c      call dptnq(r,f,no,vdr,eest)

      return
      end
C#######################################################################
      subroutine vpol(r,vp,rho,h,no)
      implicit real*8(a-h,o-z)
      dimension r(1),vp(1)

      do  j = 2, no


         if(r(j).eq.0.0D+00) then

            write(*,*) ' vpol::Division with zero r(',j,' ) = ', r(j) 
            stop
         endif


         vp(j) = 1.0D+00

         if(rho.ne.0D+00) then 

            a     = (r(j)/rho)**6
         
         else

            a = 0.0D+00 

         endif


         if(a.gt.1.0D+02) then

         
            vp(j) = dsqrt(1.d0-dexp(-a))/r(j)

         else

            vp(j) = vp(j) / r(j)
         
         endif

      enddo

c      print*,'r(1)=',r(1),'r(2)=',r(2)
c      print*,'vp(1)=',vp(1),'vp(2)=',vp(2)
c      stop


      return
      end
C####################
      subroutine q2tst1(lhf,ll,nhf,nllmin,ir,ic,nr,nc,vxs)
      implicit real*8(a-h,o-z)
      dimension lhf(1),ll(1),nhf(1),nllmin(1)

      if(nc.gt.1 .and. nr.gt.1) return

      if(nhf(ir).eq.nllmin(ir) .and. lhf(ir).eq.ll(ir)) then
        if(nr .gt. 1) goto 10

	vxs=vxs/dsqrt(2.d0)

  10    endif

      if(nhf(ic).eq.nllmin(ic) .and. lhf(ic).eq.ll(ic)) then

        if(nc.gt.1) goto 11

	vxs=vxs/dsqrt(2.d0)
  11    endif

      return
      end
c-----------------------------------------------------------------------------
      subroutine q2tst0(lhf,ll,nhf,nllmin,ir,ic,nr,nc,vxs)
      implicit real*8(a-h,o-z)
      dimension lhf(1),ll(1),nhf(1),nllmin(1)

      if(nc.gt.1 ) return

      if(nhf(ic).eq.nllmin(ic) .and. lhf(ic).eq.ll(ic)) then

	vxs = vxs / dsqrt(2.d0)

      endif

      return
      end
c-----------------------------------------------------------------------------