!###################################################################
      subroutine getncs(ncfg, nhf,lhf,ll,nmin,nmx,nol,is,
     &     l,ls,ncsmx,ncmx,ndi,idcs)
      parameter(nsx=100) 
      dimension nhf(nsx),lhf(nsx),ll(nsx),nmin(nsx),nmx(nsx),nol(nsx)
      dimension is(nsx),idcs(nsx),ndi(nsx)
!..................................................................

      
      CALL HCFGFILE(NCFG, l)
      
      read(ncfg, 3) nSymmetries
      read(ncfg, 3) l, ls
      read(ncfg, 3) ncsmx
      WRITE(*,*) NSYMMETRIES, L, LS, NCSMX
      ncmx = 0
      do   k = 1, ncsmx

         read(ncfg, 3) nhf(k),lhf(k),ll(k),nmin(k),nmx(k),idcs(k)

C         write(*,*) nhf(k), lhf(k), ll(k), nmin(k), nmx(k), idcs(k)

C#
C#     nol :  total number of B-splines
C#    is   :  global position of the vector of basis channels
C#    ncmx :  Number of B - splines for channel   
C#                                              | nhf, lhf, ll >
C#

         nol(k) = nmx(k) - nmin(k) + 1
         is(k)  = ncmx + 1 
         ncmx   = ncmx + nol(k)
         
      enddo

!      READ(NCFG, 3) itry
!      READ(NCFG, 7) ( ndi(k), k = 1, ncsmx)

      do k = 1, ncsmx

         ndi(k) = nol(k) 

         if(ndi(k).ne.nol(k)) then 
            write(*,*) ' error in CFG file '
            write(*,*) '   K = ', k 
            write(*,*) ' NDI = ', ndi(k)
            write(*,*) ' NOL = ', nol(k)

         endif

         write(*,*) k, ndi(k), nol(k) 

      enddo

      write(*,*)  ' TOTAL ANGULAR MOMENTUM           L  = ', l
      write(*,*)  ' TOTAL SPIN MOMENTUM              S  = ', ls
      write(*,*)  ' NUMBER OF CHANNELS              NCS = ', ncsmx
      write(*,*)  ' TOTAL NUMBER OF BASIS CHANNELS  NCF = ', ncmx

c      read(ncfg, *) (nd(k), k = 1, ncsmx)

      if(ls.eq.0.or.ls.eq.1) then 

         ls = 2*ls + 1 
      else

         ls = ls
      endif

      close(NCFG)
C...................

 3    FORMAT(6i5)
 7    FORMAT(10i4)
C..................

      return
      end
!!!###################################################################
      subroutine GETNUMBERCHANELLS(NCFG, L, NDMAX)
      IMPLICIT NONE
      INTEGER NCFG,NSYMMETRIES,L,LS,NDMAX

!..................................................................

      
      CALL HCFGFILE(NCFG, L)
      
      read(ncfg, 3) NSYMMETRIES
      read(ncfg, 3) L,LS
      read(ncfg, 3) NDMAX

      close(NCFG)
C...................

 3    FORMAT(6i5)
C..................

      return
      end
!!!###################################################################
!!!#  EOF
