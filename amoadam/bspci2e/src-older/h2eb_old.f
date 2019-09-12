C############################################################
      program diag2e
      implicit doubleprecision(b-h,o-z)
      INCLUDE "parameter.2e.inc" 
C      PARAMETER(nhx=2500, ncs=30, lf = 10)
      dimension h(nhx,nhx), u(nhx,nhx), ih(nhx)
      dimension  lhf(ncs),     ll(ncs) 
      dimension  nhf(ncs), nllmin(ncs), nllmax(ncs)
      dimension noll(ncs),     is(ncs),     nd(ncs), ndi(ncs)
      dimension n(lf)
      CHARACTER*100 ARGV 
C..........................


C#    Read symmetry from inp/wf2e.inp input file

      CALL GETARG(1, ARGV)
      READ(ARGV,*) INPL
      WRITE(*,*) ' CALCULATE FOR L =',INPL 

C#   read configurations from inp/cfg-L.inp file 

      CALL CFGFILE(15, INPL) 
      
      read(15, *) nsymx

C      do 200 nsy = 1,  nsymx         

      call cfin(nhf,lhf,ll,nllmin,nllmax,noll,is,lo,ls,ncsmx,nc)

      if(nc.gt.nhx) then
         write(*,*)' Error:: nc exceeds nhx!. Exiting.'
         stop
      endif

      do  i = 1, ncs
        n(i) = nhf(i) + lhf(i)
      enddo
  
      read(15,3) itry

      if(ncsmx.gt.NCS)  then
         write(*,*)' Error:: ncsmx or itry exceeds maximum! Exiting...'
         stop
      endif

C      do 17 it = 1, itry
      read(15,5) ( nd(k), k = 1, ncsmx), iouti, ioutf

C      do  k = 1, ncsmx
C        nd(k) = nllmax(k) - nllmin(k) + 1
C       enddo
 
      do 15 k = 1, ncsmx
         if(nd(k).lt.nllmin(k)) go to 15
         nd(k) = nd(k) - nllmin(k) + 1
 15      continue

C     # 17    enddo

C============================

      if(inpl.ne.lo) then
         write(*,*) ' Input for L from cfg file differ. Exiting'
         stop
      endif

C============================

C*                output on "diagxxal.dat" :
C------------------------------------------------------------------
C*    Lo and Ls   ::  Orbital and Spin Angular Momenta
C*    NCSMX       ::  Number of Configuration Series
C*    NHF(K), LHF(K), LL(K), NLLMIN(K), NLLMAX(K),   K = 1, ncsmx
C*    ITRY        ::  Number of data sets in the file
C*    ND(K),      ::  Number of conf. included from each series
C*                     K=1,NCSMX 
C*    NHMX        ::  total number of config. in the eigenvectors
C*    NDTOTAL     ::  Number of states with output in the file, 
C*                    repeating energy eigenvalue and
C*                    energy eigenvector  (1 to nhmx)


C......................................

      call d2efile(9, lo)

      write(9)    lo,ls
      write(9)    ncsmx
      write(9) ( nhf(k),    k = 1, ncsmx )
      write(9) ( lhf(k),    k = 1, ncsmx )
      write(9) ( ll(k),     k = 1, ncsmx )
      write(9) ( nllmin(k), k = 1, ncsmx )
      write(9) ( nllmax(k), k = 1, ncsmx )
      write(9)   itry

C.....................................


C      do 100 it = 1, itry


C#
C#     read   1/r_12 data
C#

         call r12file(3, lo)

C.........

         write(*,*) nc
         do 10 ir = 1, nc
            read(3) ( u(ir,ic), ic = ir, nc)

C#  fill matrix since symmetric

            do 10 ic = ir, nc

 10            u(ic,ir) = u(ir,ic)

               close(3)

C..........

               do  j = 1, ncsmx

                  ndi(j) = nd(j)
               enddo

               write(9) ( ndi(j), j = 1, ncsmx)

               do  i = 1, nc

                  ih(i) = 0
               enddo

               do 71 k = 1, ncsmx

                  if(ndi(k).eq.0) go to 71
        
                  ihd = is(k) - 1
            
                  do 72 id = 1, ndi(k)

                     ihd = ihd + 1
 72                  ih(ihd) = 1

 71               continue


C# include or not configurations
                  k = 0
                  do 20 ir = 1, nc
               
                     if(ih(ir).eq.0) go to 20

                     k = k + 1 
                     kk = 0
                     do 25 ic = 1, nc

                        if(ih(ic).eq.0) go to 25

                        kk = kk + 1
                        h(k, kk) = u(ir, ic)
                  
 25                  continue
 20               continue
            
            

            write( *, 9 ) ( n(i), lhf(i), ll(i), i = 1, 6 )


C*    diagonalize 2-e hamiltonian            
C...............

            call diag( h, u, k,iouti, ioutf, ndi)

C...............
            
C 100     continue

         close(16)
         close(9)
      
C# 200  continue

C..............................

C#
C#                             FORMAT STATEMENTS
C#

 1    format(/2x,'# of conf. series=',i4/)
 3    format(15i5)
 4    format(/2x,'data set # =',i3)
 5    format(10i4)
 6    format(/2x,'l and ls=',2i3/)
 8    format('  #',5x,'energies',5x,6(1x,3i2),'tot per')
 9    format('  #',5x,'En(Ryd)  ',6x,3i2,4x,5(1x,3i2),'   Sum_P(10)')
 7    format(3x,'total # of config. included: n=',i4/)
 76   format(10a14)
 77   format(2hok,i1)
 566  format(5a16)

C..............................

      end

C#############

      subroutine diag(h,u,nhmx,ndi,ndf,nd)
C     !     
      INCLUDE "parameter.2e.inc"
C      PARAMETER(NCS= 30)
      integer nhmx, nhx, ndtotal, ndi, ndf
      integer i, j, k, k1 
      integer lwork, il, iu, info
      double precision abstol, vl, vu, enorm, sumProb10
      integer ifail(nhmx), iwork(5*nhmx), nd(NCS)
      double precision en(nhmx), work(8*nhmx), p(ndf)
      double precision h(nhx,nhx), u(nhx,nhx) 
Cxxxxxxxx  format statements
    1 format(4i5)
    2 format(2x,1p9d14.6)
    3 format(1x,8f9.5)
    4 format(/2x,'eigen-vector --'/)
    6 format(2x)
    8 format(i3,1x,1pe15.8,e12.4,1x,5(f7.3),1x,e12.4)
   77 format(2hok,i1)
Cxxxxxxxxxxxxxxxxxxxx

c*    TOTAL NUMBER OF EIGENSTATES PRINTED OUT

      ndtotal = ndf - ndi + 1

c**  lapack routine

      lwork = 8 * nhmx
      vl = 0.0D0
      vu = 0.0D0
      il = 1
      iu = ndtotal 
      abstol = 0.0d0

      call dsyevx('V','I','U',nhmx, h,nhx,vl,vu,il,iu,abstol,M,
     1     en,u,nhx,work,lwork,iwork,ifail,info)


Cimsl      call devesf(nhmx,ndtotal,h,nhx,small,en,u,nhx)



C....              store 

      write(9) nhmx
      write(9) ndtotal

      do 30 k = 1, ndtotal 


C*    normalize eigenvectors
Cimsl      enorm = 0.d0
Cimsl      do j = 1, nhmx
Cimsl         enorm = enorm + u(j,k) * u(j,k)
Cimsl      enddo
Cimsl      do  j = 1, nhmx
Cimsl         u(j,k) = u(j,k) / dsqrt(enorm)
Cimsl      enddo


      write(9)  en(k) 
      write(9)  ( u(j,k), j = 1, nhmx )

C---------------------------------------------------------------
C*  CALCULATE % CONTRIBUTION (PROBABILITY DENSITY) OF 
C*  CONFIGURATION   |n1 l1 l2 > (SUM TO ALL N2) TO EIGENSTATES  
C---------------------------------------------------------------
C*  EXAMPLE :
C*
C*
C*  In the CI procedure the eigenstates of atom with quantum 
C*  number of LS is a sum over the 2e configurations  (n1 l1,n2 l2)
C*  (allowed by the Clebsch-Gordan coefficients).
C*  therefore:
C*
C*                 ____ 
C*                 \ 
C*   |He(gs):LS > = \     C(n1 l1,n2 l2)  * | LS : (n1 l1,n2 l2) > 
C*                  /
C*                 /___  
C*          (to all n1 l1,n2 l2)  
C*   
C*     
C*   We can calculate the contribution of the 1ss as :
C*
C*   1ss:      n1 = 1, l1 = s , l2 =s
C*
C*            ___  
C*   p(1ss) = \   C(1s n2 s) * C(1s n2 s)
C*            /__ 
C*           (n2)
C*          
C*  Also
C*  We can calculate the contribution of the 2ss as :
C*
C*   2ss:      n1 = 2, l1 = s , l2 = s
C*            ___  
C*
C*   p(2ss) = \   C(2s n2 s) * C(2s n2 s)
C*            /__ 
C*           (n2)
C*          
C*
C* 
C*              Energy      p(1ss)          p(2ss)   ..... and so on 
C*
C*   HE(1S^2) :: xxx       p(1) * 100      p(2)*100   ..... p(i)*100 
C*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C*     calculate for the first 10 configuration series
C*  
C*  for the  calculation of each of the | u(nln'l')|^2             
C*            uncomment the 'c*1' comments 
C* 
C*

      k1 = 0
      sumProb10 = 0.0d+00

      do 27 i = 1, 10
       
         p(i) = 0.0D+00
         
         do 26 j = 1, nd(i)
            
            k1 = k1 + 1

            if(k1.gt.nhmx)  goto 26

c*1            tmp  = u(k1,k) * u(k1,k) 

            p(i) = p(i) +  u(k1, k) * u(k1, k)

c*1            if(k1.lt.5) write(*,*) j, tmp

 26      continue

c*1         write(*,*) i, p(i)

         sumProb10  = sumProb10  + p(i)
        
 27   continue

      write(*,8) k, en(k)*0.5D+00,( p(i), i = 1, 6), sumProb10

   30 continue
      return
      end

C#######################################################################

      subroutine cfin(nhf,lhf,ll,nmin,nmx,nol,is,l,ls,ncsmx,ncmx)
      dimension nhf(1),lhf(1),ll(1),nmin(1),nmx(1),nol(1),is(1)
    3 format(5i5)

      read(15,3) l,ls
      read(15,3) ncsmx

      ncmx = 0
      do  k = 1, ncsmx

      read(15,3) nhf(k), lhf(k), ll(k), nmin(k), nmx(k)

      nol(k) = nmx(k) - nmin(k) + 1
      is(k)  = ncmx + 1
      ncmx   = ncmx + nol(k)

      end do

      write(*,*) ' Total angular momentum      L = ', l
      write(*,*) ' Total Spin                  S = ', ls
      write(*,*) ' No of channel series    ncsmx = ', ncsmx
      write(*,*) ' No of configurations     ncmx = ', ncmx
      write(*,*) 

      return
      end
C#######################################################################
C* EOF








