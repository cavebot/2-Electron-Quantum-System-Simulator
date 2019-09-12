      program d1e

      implicit real*8(a-h,o-z)
      INCLUDE "parameter.1e.inc"
      dimension t(ns+nk)
      dimension amat(ns,ns),bmat(ns,ns),work(ns,ns)
      dimension cmat(ns,ns),er(ns)
      dimension fv1(ns),fv2(ns),ogp(ns)
      dimension r(np),p(nl,np),vd(np),ff(np),dr(np),ad(ns,ns),
     1 ae(ns,ns),rs(nl),yk(np),zk(np),wy(np),wz(np),
     2 br(np),bc(np),pod(nl,np),en(nl,ns),dtm(nl),n2(nl)
      common/params/alp1(nl),r01(nl),redmass,znuc,lang
      common/bspf/bf(ns-2,np)
      common/rarr/rpw(nl+3,np)


C...................................................................
      open(14,file='inp/d1e.inp',status='old')
      read(14,*) znuc, rmin, rmax
      read(14,*) no, redmass
      read(14,*) nos, nop, nod
      read(14,*) isl, nw
      read(14,*) lmin, lcore, lmax, itrmx
      read(14,*) ( alp1(i + lmin - 1), i = 1, nw)
      read(14,*) ( r01(i  + lmin - 1), i = 1, nw)
      read(14,*) ( zk(i),  i = 1, nos + nop + nod)
      read(14,*) ( n2(i),  i = 1, nw), k
      read(14,*) ( rs(i),  i = 1, nw)
      read(14,*)  crt, di, df, id
      close(14)

C...................................................................
 
      open(16,file='out/d1e.out')

      write(*,*)  '***********************************************'
      write(*,*)  'Number of          nBsp = ', n2(1)
      write(*,*)  'order Bspl.           k = ', k
      write(*,*)  'Box  radius           R = ', rmax
      write(*,*)  'First point         fkn = ', rs(1)    
      write(*,*)  ' # of s,p,d orbitals    = ', nos,nop,nod
      write(*,*)  '                    isl = ', lcore
      if((isl.ne.1).and.(itrmx.ne.1)) then
         write(*,*)'(isl.ne.1) Hartree-Fock     lcore = ', lcore
      else 
         write(*,*)'(isl.eq.1) Hydrogen or existing core coeff used'
      endif


      if(redmass.eq.1) then
         write(*,*)'Hyd - like functions are caclulated'
      else if(redmass.eq.2) then
         write(*,*)'Ps - like functions are calculated'
      else 
         write(*,*)'Not acceptable value for param redmass. Exiting...'
         stop
      endif

      if((alp1(lmin).eq.0.D0)) then
         write(*,*)' No core -polarization '
      else
         write(*,*)' Core - polarization is included'
      endif
      write(*,*) '**************************************************'

      nsp = nos + nop
      ntc = nsp + nod

C#     some checks first
C#      

      do  i = lmin, lmax - lmin + 1

         if(n2(i).gt.ns.or.n2(i).lt.(k+2)) then

            write(*,*)' nk+2 < n2  < k+2'
            write(*,*)'  n2 = ', n2(i)
            write(*,*)'  ns = ',  ns
            write(*,*)' k+2 = ',  k+2


            goto 999
         endif
      enddo

      if(k.gt.nk .or. k.lt.2) then
         write(*,*)' problem with k !! k=',k
         goto 999
      endif

      if(di.gt.1.0D+00 .or. df.gt.1.0D+00)then
         write(*,*)'   di  or   df greater than 1  '
          goto 999
      endif

C#########################
C#
C*  Generate knot points. Data provided from 'rin.din'
C*  Any waveFunction will be calculated on this knot sequence.
C*
C***********************************************************


      call rin(r, dr, h, no, idr)

 
c*** ....... write knots in 'knotFile'

      open(26, file ='dat/knot.dat')

      do ii = 1, no

         write(26,*) ii,'  ', r(ii)      

      enddo
 
      close(26)
C..................................
C#       prepare H-like core wf's

      do 28 kk = 1, nl+3
         do 28 j = 1, no

            rpw(kk,j) = r(j)**kk

 28   continue

      if (isl.eq.1) then

         do 40 i = 1, nos

           if (r(1).ne.0.0d+00) p(i,1) = hwf(i, 0, zk(i), r(1), er(i) )

           do 40 j = 2, no
 40           p(i,j) = hwf(i,0,zk(i),r(j),er(i))

         if(nos.eq.ntc) go to 45

         do 41 i = nos + 1, nsp

           if (r(1).ne.0.d0) p(i,1) = hwf( i, 1, zk(i), r(1), er(i))

           do 41 j = 2, no
 41           p(i,j) = hwf(i, 1, zk(i), r(j), er(i))

         if((nsp).eq.ntc)go to 45

         do 42 i = nsp + 1, ntc

           if (r(1).ne.0.0D+00) p(i,1)=hwf(i,2,zk(i),r(1),er(i))

           do 42 j = 2, no
 42           p(i,j) = hwf(i, 2, zk(i), r(j), er(i))

 45   endif

C#
C#  read existing coef. & generate core wf's.
C#

      if(isl.ne.1) then

         ni = 1
         nf = nos

         do 15 l = 1, lcore

            lang = l - 1

C........ open core wf files

            call corefile(1, lang) 

            read(1) rs1, rmx, n1, k1, lin

            call mkgrid(n1,k1,rmx,rs1,t)


            do 32 npp = 1, n1 - 2

               read(1) en(l, npp), (cmat(j,npp),j = 1, n1 - 2)

 32         continue


            write(16,*)' l = ',lang,'  Initially Input Energies =:  '

            write(16,8) ( 0.5 * en(l,i),i = n1 - 2,1,-1)


C#    subroutine corewf() takes cmat, r, dr, n, k, ni, nf, ns,no, ntc & h,
C#    returns core wf's p(i,j)(( i = ni, nf ),(j = 1, no) ).


            call corewf(cmat,p,ff,r,dr,ogp,t,n1,k1,ni,nf,ns,no,ntc,h)


            ni = nf + 1

            if(l.eq.1) nf = nf + nop

            if(l.eq.2) nf = nf + nod

            close(1)

 15      continue
      endif

C#       this loop improves coef. for core orbits for itrmx times.
C#       d1 = di + dfloat(itry) / dfloat(id) * (df - di)

       nm = ns

C#
C#     l > lcore   ---> do not calculate HF orbitals
C#       

      if (lmin.gt.lcore) goto 201


      do 200 itry = 1, itrmx

         d1  = di
         if(itry.gt.id) d1 = df

         write(*,*)'  #####################################'
         write(*,*)'  Iteration Number : ', itry

         write(16,*)'  #####################################'
         write(16,*)'  Iteration Number : ', itry

         ni = 1
         nf = nos
         dtma = 0.0d+00

         do 100  ll = 1, lcore

         write(*,*)'....................................... L = ', ll-1

            dr0 = rs(ll)
            n   = n2(ll)

            lang = ll - 1

            call d1efile(1, lang)

            write(16,3) znuc,h,rmin,dr0,rmax,n,k,no,idr,lang,d1
            
c*      geteig(...) takes lang,n,k,dr0,rmax,no,nos,nop,
c*                        ntc,n,m,r,dr,p,idr,
c*      returns coef. cmat & energies er. 
c*      other parameters only provide storage space to be 
c)      used inside it.

            call geteig(lang,n,k,dr0,rmax,no,nos,nop,ntc,nm,h,t,r,dr,
     1              p,vd,ff,yk,zk,wy,wz,br,bc,amat,bmat,cmat,ad,ae,er,
     2              redmass, fv1, fv2, work,ogp,idr)
               

               do  i = ni, nf
                  do  j = 1, no

                     pod(i,j) = p(i,j)

                  enddo
               enddo

               call corewf(cmat,p,ff,r,dr,ogp,t,n,k,ni,nf,ns,no,ntc,h)
               
C*          orn = |< wf_n | wf_n > | ^2   ( should be == 1)

               do 22 i = ni, nf

                  ff(1) = 0.0D+00

                  do j = 2, no

                     p(i,j) = d1 * p(i,j) + (1.D+00 - d1) * pod(i,j)

                     ff(j)  = p(i,j) * p(i,j) * r(j) / dr(j)

                  enddo

                     orn = rint(ff, 1, no, 14, h )
                     sq  = SQRT(orn)
                     
                     if(itrmx.ne.1) then

                        IF(LL.EQ.1.AND.I.EQ.1) THEN 

C                           WRITE(*,*) ' ORBITAL      P_1s  : ', i 
                           WRITE(*,*) ' <P_1s | P_1s >     = ',orn

                        ELSE IF (LL.EQ.1.AND.I.EQ.2) THEN 

C                           WRITE(*,*) ' ORBITAL      P_2s  : ', i
                           WRITE(*,*) ' <P_2s | P_2s >     = ',orn

                        ELSE IF (LL.EQ.1.AND.I.EQ.3) THEN 

C                           WRITE(*,*) ' ORBITAL      P_3s  : ', i
                           WRITE(*,*) ' <P_3s | P_3s >     = ',orn

                        ELSE IF (LL.EQ.2.AND.I.EQ.3) THEN 

C                           WRITE(*,*) ' ORBITAL      P_2p  : ', i
                           WRITE(*,*) ' <P_2p | P_2p >     = ', orn

                        ELSE IF (LL.EQ.2.AND.I.EQ.4) THEN 

C                           WRITE(*,*) ' ORBITAL      P_2p  : ', i
                           WRITE(*,*) ' <P_2p | P_2p >     = ', orn

                        ELSE IF (LL.EQ.2.AND.I.EQ.5) THEN 

C                           WRITE(*,*) ' ORBITAL      P_3p  : ', i
                           WRITE(*,*) ' <P_3p | P_3p >     = ', orn

                        ELSE

                           WRITE(*,*) ' ORBITAL  : ', i
                           WRITE(*,*) ' NORM     = ', orn

                        ENDIF

                     endif

                     do 22 j = 1, no

 22                     p(i,j) = p(i,j) / sq 
                        

                      
c***  from here to next comment line are auxilory. 

                      if((itry.eq.1).and.(n1.ne.n)) goto 52
         
                      dltem = 0.0D+00

                      do 51 ie = 1, n - 2
                              
                         dlte = dabs(er(ie) - en(ll,ie))

                         if(dltem.lt.dlte) then

                            kie   = (n-2) - ie + 1
                            delt  = er(ie) - en(ll,ie)
                            dltem = dlte

                         endif

 51                   continue

 52                   do 53 ie = 1, n - 2

 53                      en(ll,ie) = er(ie)
                              

                         dtm(ll) = 0.0D+00
                         
                         do i = ni, nf
                            do j = 1, no

                               dlt = dabs(p(i,j) - pod(i,j))

                               if(dtm(ll).lt.dlt) dtm(ll) = dlt

                            enddo
                         enddo

C**   ...........log

                         write(16,17) kie,delt,dtm(ll)
                         write(16,*)' ENERGIES :  '
                         write(16,8) (0.5 * er(j),j = n-2, 1, -1)

                         if(dtma.lt.dtm(ll))  dtma = dtm(ll)

                         if((dtm(ll).lt.crt).or.(itry.eq.itrmx)) then
  
                            
                            write(1) dr0, rmax, n, k, lang

c***  write in binaries in increasing eigenvalue order

                            do i = 1, n - 2

                               write(1) er(i), (cmat(j,i),j=1,n-2)              
                            enddo
                                    
                         endif

                         if(dtma.lt.crt.and.ll.eq.lcore) goto 201
                         ni = nf + 1
                        
                         if(lang.eq.0) nf = nf + nop
                         if(lang.eq.1) nf = nf + nod

                         close(1)

 100               continue
 200    continue


C#
C#       if  lcore < lmax     -->   leave
C#

 201    if(lcore.ge.lmax) goto 999

        if(lmin.gt.lcore) lmn = lmin

        if(lmin.le.lcore) lmn = lcore + 1
                           
c**   this loop uses coef for core to calculate coef. for non-core orbits.
        do 205 ll = lmn, lmax

           lang = ll - 1
           dr0  = rs( ll - lmn + 1 )
           n    = n2( ll - lmn + 1 )

           call d1efile( 1, lang)  

           write(*,*) ' ll ', ll
           write(*,*) ' lmin = ', lmn
           write(*,*) ' lmax = ', lmax
           write(*,*) ' r(',lang,') = ', dr0
           write(*,*) ' n(',lang,') = ', n
           write(*,*) '**********************'
          
           write(16,3) znuc,h,rmin,dr0,rmax,n,k,no,idr,lang

             call geteig(lang,n,k,dr0,rmax,no,nos,nop,ntc,nm,h,t,r,
     1          dr,p,vd,ff,yk,zk,wy,wz,br,bc,amat,bmat,cmat,ad,ae,er,
     2          redmass, fv1, fv2, work, ogp, idr)

             write(16,*)' En = '
             write(16,8) ( 0.5 * er(j),j = n - 2, 1, -1 )

C     #  write in binaries in increasing eigenvalue order

           write(1) dr0, rmax, n, k, lang
           
           do i = 1, n - 2
              write(1) er(i),( cmat(j,i), j = 1, n - 2 )
           enddo

C           write(*,*) er(40), ( cmat(j, 40), j = 1, n - 2 )

           close(1)
        
 205    continue

C     --- call tst(t,znuc,rmax,xi,r1,n,k,lang,no,nrp)

 999    stop

C...................................................................
    1 format(2x,'rs(ll) for b-spline, ll=llmin,llmax'/)
    2 format(2x,'n & k=?'/)
    3 format(/2x,'z,h,r(1),rs,rmx,n,k,no,idr & l= ',1p4e11.3,1pe14.6,
     1 5i5/2x,'fraction of new wf used: d1=',e9.2)
    5 format(10i5)
    7 format(5e13.4)
    8 format(1p7d18.10)
    9 format(2x,'llmin & llmax=?'/)
 16   format(' rmaxs in input data do not agree !'/
     * ' they are',f6.1,' in coe.dx and ',f6.1,' in bspeg.din ')
     & 
 17   format('   the',i3,'th dlt-e (max in magnitude) & max of dlt-p =',
     1 1p2e16.7)
      write(*,*) 'i got this far'
c**   read input data
 566  format(5a16)
C...................................................................
        end

C#######################################################################
C#
C#    hwf :: 
C#
C#
C#######################################################################

c      function hwf(n, l, z, r, e)
c      double precision hwf,z,r,e,fn,rk,rm,rll,fk,fm,fll,a,p,x
C----------------------------------------------------------------------
c      m   = n + l
c      k   = n - l - 1
c      ll  = 2*l + 1
c      fn  = n
c      rk  = k
c      rm  = m
c      rll = ll
c      fk  = 1.d0
c      fm  = 1.d0
c      fll = 1.d0
c      p   = 1.d0
c      a   = 1.d0
c      x   = -2.d0*z*r / fn
c
c      if (dabs(x) .gt. 160.d0) go to 20
c      
c     do  i = 1, m
c
c         fm = fm * rm
c         rm = rm - 1.d0
c
c      enddo
c
c      rm = m
c      do i = 1, ll
c
c         fll = fll * rll
c         rll = rll - 1.d0
c
c      enddo

c      if (k) 1,2,3

c    3 do 4 i = 1, k

c      p   = 1.d0+a/rk*p/rm*x
c     fk  = fk*rk
c      a   = a  + 1.d0
c      rk  = rk - 1.d0
c    4 rm  = rm - 1.d0
c    2 hwf = dsqrt(z*fm/fk)*p*dexp(x/2.d0)*(-x)**(l+1)/fll/fn
c      e  = - z**2 / fn**2

c      return

c   20 hwf = 0.d0
c      return

c    1 write(16,10) n, l, z, r
c      write(*,10)  n, l, z, r

c   10 format(52h  forbidden combination of n and l in hwf subprogram/
c     1     7h     n=,i3,5h   l=,i3,5h   z=,f6.3)
c      stop
c      end
C##################################################################


c****************************************************************************
c** description of the data in the input file:
c**  flm --  output file name;
c**  fln(1-5) -- unformatted in/out files (containing the in/out coefficients)
c**  znuc -- nuclear number fo the atom
c**  rmin, rmax -- minmum & maximum of radius to be used
c**  no -- total # of points; idr=0: sine-scale; idr=1: log-scale
c**  nos, nop & nod -- # of s, p & d orbits in the frozen core
C#  isl =  1 : use hydrogenic wavefunctions; 
C#  isl != 1 : use existing bsp-coe input 
C#
C#  redmass = 1 : Hydrogenic  - like
C#            2 : positronium - like
C#
c**  nw -- # of desired symmetries, e.g., if s to g wf's are needed, nw=5;
c**                                       if s to f needed, nw=4
c**  lcore -- largest l(+1) in core;
c**  lmin & lmax -- min. & max. l(+1) in this calculation: 1 <= lmin <= lmax
c**  itrmx --  max. # of iterations desired to do
c**  alp(i) static polarizabilities for each l of core-polarization potentials.
c**  r01(i) cut-off parameters for the  l- core polarization potentials 
c**  zk(i) here are the effective z's  of the core orbits.
c**  n2(i) --  number of collocation points for each l
c**  k --  order of spline
c**  rs(l) --  r value of first node point for corresponding l
c**  di, df initial & final fraction of new wf; crt-- criterion of calculation.
c**  flnm(1-lcore) --  useful only if isl!=1: unformatted coef files as input.)
c******************************************************************************
C#  EOF 












