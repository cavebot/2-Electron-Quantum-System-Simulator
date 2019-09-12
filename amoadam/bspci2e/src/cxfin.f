!###################################################################
      subroutine cxfin(ncfg,nhf,lhf,ll,nmin,nmx,nol,is,
     &     l,ls,ncsmx,ncmx,nd,idcs)
      dimension nhf(1),lhf(1),ll(1),nmin(1),nmx(1),nol(1),is(1)
      dimension nd(1), idcs(1)
!..................................................................
    3 format(6i5)



      
      read(ncfg, 3) nSymmetries
      read(ncfg, 3) l, ls
      read(ncfg, 3) ncsmx

      ncmx = 0
      do   k = 1, ncsmx

         read(ncfg, 3) nhf(k), lhf(k), ll(k), nmin(k), nmx(k), idcs(k)

         nol(k) = nmx(k) - nmin(k) + 1
         is(k)  = ncmx + 1
         ncmx   = ncmx + nol(k)
         
      enddo

C      read(ncfg,*) (nd(k), k = 1, ncsmx)


      write(*,*)  ' Angular Quantum Number L  = ', l
      write(*,*)  ' Spin    Quantum Number S  = ', ls
      write(*,*)  ' Number of channels    NCS = ', ncsmx
      write(*,*)  ' Numbers of configs    NCF = ', ncmx
C      close(ncfg)

      return
      end
!###################################################################
!  EOF
