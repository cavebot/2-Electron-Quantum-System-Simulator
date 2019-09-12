module jl
  implicit none
  integer j1, j2, ls
end module jl

module wf_files
  implicit none
  character(len=16), allocatable, dimension(:) :: aadp, badp, ba_over
  real(8), allocatable, dimension(:,:,:) :: vz, vzba, sp_a
  contains
    subroutine read_dp(nsx, ltotalMax)
      implicit none
      integer i, j, m, n, li, lf, nwfi, nwff, nsx, ir, ic, ltotalMax

      allocate(vz(ltotalMax,nsx,nsx))

      write(*,*) 'ltotalMax = ' , ltotalMax

      do i = 1, ltotalMax 

        open(4,file=aadp(i),form='unformatted',access='sequential')
        write(*,*) 'opening', i 


        read(4)  li, lf, nwfi, nwff
        write(*,*) 'li=', li, 'lf=',lf, 'nwfi=', nwfi, 'nwff=', nwff
        read(4) ((vz(i,ir,ic), ir=1,nwff), ic=1,nwfi)

        close(4)
      end do

      write(*,*) 'vz read'

      allocate( vzba(2*ltotalMax ,nsx,nsx) )

      do i = 1, 2 * ltotalMax

        open(4,file=badp(i),form='unformatted', access='sequential') 

        read(4) li,lf,nwfi,nwff
        read(4) ((vzba(i,ir,ic), ir=1,nwff), ic=1,nwfi)

        close(4)

      end do
      write(*,*) 'vzba read'

      allocate(sp_a(ltotalMax+1,nsx,nsx))

      do i=1, ltotalMax + 1

        open(4,file=ba_over(i),form='unformatted', access='sequential') 

        read(4) li,nwfi,nwff
        read(4) ((sp_a(li+1,ir,ic), ir=1,nwff), ic=1,nwfi) 

        close(4)

      end do

      write(*,*) 'sp_a read'

  end subroutine read_dp

end module wf_files
  

module configuration
  implicit none
  integer, allocatable, dimension(:) :: nhfi,lhfi,lli,nmini,nmxi,&
   &ndi,nhff,lhff,llf,nminf,nmxf, ndf, isf, nolf, idcs
  contains
    subroutine config_space(ncs)  
      implicit none
      integer ncs
      allocate(nhfi(ncs))
      allocate(nhff(ncs))
      allocate(lhfi(ncs))
      allocate(lhff(ncs))
      allocate(lli(ncs))
      allocate(llf(ncs))
      allocate(nmini(ncs))
      allocate(nminf(ncs))
      allocate(nmxi(ncs))
      allocate(nmxf(ncs))
      allocate(ndi(ncs))
      allocate(ndf(ncs))
      allocate(nolf(ncs))
      allocate(isf(ncs))
      allocate(idcs(ncs))
    end subroutine config_space
end module configuration  

module dz_value
  use wf_files
  implicit none
  real(8), allocatable, dimension(:,:) :: dz
  contains
    subroutine cal_nlsp(ir, ic, id, ang, phase) 
      use configuration
      implicit none
      integer i, nr, nc, mr, mc, ir, ic, idn
      integer, dimension(4) :: id


      real(8), dimension(4) :: ang
      real(8) phase

      do 20 i=1,4
      select case(id(i).eq.0)
        case(.false.)
        select case(i)
          case(1)
          do nr=1,ndf(ir)
            mr=nminf(ir)+nr-1
            do nc=1,ndi(ic)
              mc=nmini(ic)+nc-1
              if(mod(id(1),2).eq.0) then
                idn=id(1)/2
                dz(nr,nc)=dz(nr,nc)+ang(1)*& 
                  &(-vz(idn,nhfi(ic),nhff(ir)))*sp_a(lli(ic)+1,mr,mc)
              else
                idn=(id(1)+1)/2
                dz(nr,nc)=dz(nr,nc)+ang(1)*& 
                   &vz(idn,nhff(ir),nhfi(ic))*sp_a(lli(ic)+1,mr,mc)
              end if
            end do
          end do

! ---------------------------------------------------------------
          case(2)
          do nr=1, ndf(ir)
            mr=nminf(ir)+nr-1
            do nc=1,ndi(ic)
              mc=nmini(ic)+nc-1
              if(mod(id(2),2).eq.0) then
                idn=id(2)/2
                dz(nr,nc) = dz(nr,nc) + phase*ang(2)&
                 &*(-vz(idn,mc,nhff(ir)))*sp_a(lhfi(ic)+1,mr,nhfi(ic))
              else
                idn=(id(2)+1)/2
                dz(nr,nc) = dz(nr,nc) + phase*ang(2)*&
                   &vz(idn,nhff(ir),mc)*sp_a(lhfi(ic)+1,mr,nhfi(ic))
              end if
            end do
          end do
! ---------------------------------------------------------------
          case(3)
          do 33 nc=1,ndi(ic)
            mc=nmini(ic)+nc-1
            select case(mc .ne. nhff(ir)) 
              case(.false.)
              do 34 nr=1,ndf(ir)
                mr=nminf(ir)+nr-1
                dz(nr,nc)=dz(nr,nc)+phase*ang(3)*&
                   &vzba(id(3),mr,nhfi(ic))
   34         continue
            end select
   33     continue
!-------------------------------------------------------------------
          case(4)  
          select case(nhff(ir) .ne. nhfi(ic)) 
            case(.false.)
            do 35 nr=1,ndf(ir)
              mr=nminf(ir)+nr-1
              do 35 nc=1,ndi(ic)
                mc=nmini(ic)+nc-1
                dz(nr,nc)=dz(nr,nc)+ang(4)*&
                   &vzba(id(4),mr,mc)
   35       continue
          end select
        end select
      end select
 20 continue
!-------------------------------------------------------------------
    end subroutine cal_nlsp
!
    subroutine cal_nlnl(ir, ic, id, ang, phase) 
      use configuration
      implicit none
      integer i, nr, nc, mr, mc, ir, ic, idn
      integer, dimension(4) :: id
      real(8), dimension(4) :: ang
      real(8) phase
      do 20 i=1,4
      select case(id(i).eq.0)
        case(.false.)
        select case(i)
          case(1)
          do nr=1,ndf(ir)
            mr=nminf(ir)+nr-1
            do nc=1,ndi(ic)
              mc=nmini(ic)+nc-1
              select case(mc.eq.mr)
              case(.true.)
                 if(mod(id(1),2).eq.0) then
                    idn=id(1)/2
                    dz(nr,nc)=dz(nr,nc)+ang(1)*& 
                    &(-vz(idn,nhfi(ic),nhff(ir)))
                 else
                    idn=(id(1)+1)/2
                    dz(nr,nc)=dz(nr,nc)+ang(1)*& 
                    &vz(idn,nhff(ir),nhfi(ic))
              end if
              end select
            end do
          end do

! ---------------------------------------------------------------
          case(2)
          do nr=1, ndf(ir)
            mr=nminf(ir)+nr-1
            select case(mr.eq.nhfi(ic))
            case(.true.)
               do nc=1,ndi(ic)
                  mc=nmini(ic)+nc-1
                  if(mod(id(2),2).eq.0) then
                     idn=id(2)/2
                     dz(nr,nc) = dz(nr,nc) + phase*ang(2)&
                     &*(-vz(idn,mc,nhff(ir)))
                  else
                     idn=(id(2)+1)/2
                     dz(nr,nc) = dz(nr,nc) + phase*ang(2)*&
                     &vz(idn,nhff(ir),mc)
                  end if
               end do
            end select
          end do
! ---------------------------------------------------------------
          case(3)
          do 33 nc=1,ndi(ic)
            mc=nmini(ic)+nc-1
            select case(mc .eq. nhff(ir)) 
              case(.true.)
              do 34 nr=1,ndf(ir)
                mr=nminf(ir)+nr-1
                if(mod(id(3),2).eq.0) then
                   idn=id(3)/2
                   dz(nr,nc)=dz(nr,nc)+phase*ang(3)*&
                   &(-vz(idn,nhfi(ic),mr))
                else
                   idn=(id(3)+1)/2
                   dz(nr,nc) = dz(nr,nc) + phase*ang(3)*&
                   &vz(idn,mr,nhfi(ic))
                end if
   34         continue
            end select
   33     continue
!-------------------------------------------------------------------
          case(4)  
          select case(nhff(ir) .eq. nhfi(ic)) 
            case(.true.)
            do 35 nr=1,ndf(ir)
              mr=nminf(ir)+nr-1
              do 35 nc=1,ndi(ic)
                mc=nmini(ic)+nc-1
                if(mod(id(i),2).eq.0) then
                   idn=id(i)/2
                   dz(nr,nc)=dz(nr,nc)+ang(i)*(-vz(idn,mc,mr))
                else
                   idn=(id(i)+1)/2
                   dz(nr,nc)=dz(nr,nc)+ang(i)*vz(idn,mr,mc)
                end if
   35       continue
          end select
        end select
      end select
 20 continue
!-------------------------------------------------------------------
    end subroutine cal_nlnl
!
end module dz_value          
