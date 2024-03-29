MODULE jl
  IMPLICIT NONE
  INTEGER j1, j2, ls
END MODULE jl

MODULE wf_files

  USE ioroutines

  IMPLICIT NONE
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: vz, vzba, vzbb
  REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: sp_a, over

CONTAINS
!!!###################################################################

  SUBROUTINE read_dp(nsx, ltotalMax, mode)

    IMPLICIT NONE
    INTEGER NOUT, NDMX1E, NRWFB, NRBB, NRWFBB, NBMX
    INTEGER i, j, m, n, li, lf, nwfi, nwff, nsx, ir, ic, ltotalMax
    INTEGER mode

!!!......................

    NRWFB  = 1
    NRBB   = 2
    NRWFB  = 3
    NDMX1E = 4
    NBMX   = 8
    NOUT   = 16
    

!!!......................

    WRITE(NOUT, *) ' MAXIMUM ANGULAR MOMENTUM FOR 1-e DIPOLES ' , ltotalMax
    WRITE(NOUT, *) ' READ 1-e DME FOR CORE STATES (FXD-CODE)  '
!


      

! < p | d | p > 

    ALLOCATE(vz(ltotalMax, nsx, nsx))

    DO i = 1, ltotalMax  

       CALL HDMX1EFILE(NDMX1E, i-1, MODE)

       READ(NDMX1E) MODE
       READ(NDMX1E) li, lf, nwfi, nwff
       READ(NDMX1E) ((vz(i,ir,ic), ir=1,nwff), ic=1,nwfi)
         
       WRITE(NOUT,*) li,lf,nwfi,nwff
       
       CLOSE(NDMX1E)
         
    END DO

    WRITE(NOUT,*) '<p|d|p> read'
      

! < B | d | B > 

    ALLOCATE(vzbb(ltotalMax, nsx, nsx))

    DO i = 1, ltotalMax  

       CALL RWFBBFILE(NRWFBB, i-1, MODE)

       READ(NRWFBB) MODE
       READ(NRWFBB) li, lf, nwfi, nwff
       READ(NRWFBB) ((vzbb(i, ir, ic), ir = 1, nwff), ic = 1, nwfi)

       WRITE(NOUT,*) li,lf,nwfi,nwff
       
       CLOSE(NRWFBB)
         
    END DO

    WRITE(NOUT,*) '<B|d|B> read'



! < B | d | p(r) > 
! < p | d | B > 

    ALLOCATE(vzba(2*ltotalMax,nsx,nsx))

    j = 0
    DO i = 1, ltotalMax   


! < B | d | p(r) > 
       CALL RWFBFILE(NRWFB, i-1, i, MODE)

       READ(NRWFB) MODE
       READ(NRWFB) li,lf,nwfi,nwff
       READ(NRWFB) ((vzba(i + j, ir, ic), ir = 1, nwff), ic = 1, nwfi)

       CLOSE(NRWFB)

       WRITE(NOUT,*) li,lf,nwfi,nwff


! < p | d | B > 

       CALL RWFBFILE(NRWFB, i, i-1, MODE)

       READ(NRWFB) MODE
       READ(NRWFB) li,lf,nwfi,nwff
       READ(NRWFB) ((vzba(i+j+1, ir, ic), ir = 1, nwff), ic = 1, nwfi)

       j = j + 1

       CLOSE(NRWFB)


       WRITE(NOUT,*) li,lf,nwfi,nwff

    END DO

    WRITE(NOUT,*) '<B|d|p> read'
    WRITE(NOUT,*) '<p|d|B> read'

    ALLOCATE(sp_a(ltotalMax, nsx, nsx))

!< B | d | B > 

    DO i = 1, ltotalMax

       CALL RBBFILE(NRBB, i-1)

       READ(NRBB) li, nwfi, nwff
       READ(NRBB) ((sp_a(li+1,ir,ic), ir=1,nwff), ic=1,nwfi)

       CLOSE(NRBB)

       WRITE(NOUT,*) li,lf,nwfi,nwff
    END DO

    WRITE(NOUT,*) '<B|d|B> read'



!< B_i | B_j > 


    ALLOCATE(over(ltotalMax + 1, nsx, nsx))

    DO i = 1, ltotalMax + 1

       CALL BMXFILE(NBMX, i-1) 

       READ(NBMX) nwfi
       DO ic = 2, nwfi + 1
          READ(NBMX) (over(i,ir,ic), ir = ic, nwfi + 1)
       END DO

       CLOSE(NBMX)

       WRITE(NOUT,*) '<B|B> read'
        
!!! FILL THE MATRIX (SYMMETRIC) 
       
       DO ic = 2, nwfi + 1
          DO ir =  ic, nwfi + 1             
             over(i, ic,ir) = over(i, ir, ic)
          END DO
       END DO

    END DO

  END SUBROUTINE read_dp
END MODULE wf_files
  
!!!########################################################

MODULE configuration
  IMPLICIT NONE
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhfi,lhfi,lli,nmini,nmxi,&
       & ndi, nhff, lhff, llf, nminf, nmxf, ndf, isi, isf, noli, nolf, idcsf, idcsi

CONTAINS
  SUBROUTINE config_space(ncsi, ncsf)  
    IMPLICIT NONE
    INTEGER ncsi, ncsf

!!! initial state 

    ALLOCATE(nhfi(ncsi))
    ALLOCATE(lhfi(ncsi))
    ALLOCATE(lli(ncsi))
    ALLOCATE(nmini(ncsi))
    ALLOCATE(nmxi(ncsi))
    ALLOCATE(ndi(ncsi))
    ALLOCATE(noli(ncsi))
    ALLOCATE(isi(ncsi))
    ALLOCATE(idcsi(ncsi))

!!! final state
    
    ALLOCATE(nhff(ncsf))
    ALLOCATE(lhff(ncsf))
    ALLOCATE(llf(ncsf))
    ALLOCATE(nminf(ncsf))
    ALLOCATE(nmxf(ncsf))
    ALLOCATE(ndf(ncsf))
    ALLOCATE(nolf(ncsf))
    ALLOCATE(isf(ncsf))
    ALLOCATE(idcsf(ncsf))
    
  END SUBROUTINE config_space
END MODULE configuration

!!!########################################################

MODULE dz_value

  USE wf_files
  IMPLICIT NONE
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: dz
CONTAINS

!!!#############################################

  SUBROUTINE cal_nlsp(ir, ic, id, ang, phase, mode) 

    USE configuration

    IMPLICIT NONE
    INTEGER MODE
    INTEGER i, nr, nc, mr, mc, ir, ic, idn
    INTEGER, DIMENSION(4) :: id
    REAL(8), DIMENSION(4) :: ang
    REAL(8) phase


!!!............................


    DO 20 i = 1, 4

       IF(id(i).NE.0) THEN

          SELECT CASE(i)
          CASE(1)
             DO nr = 1, ndf(ir)

                mr = nminf(ir) + nr - 1

                DO nc = 1, ndi(ic)

                   mc = nmini(ic) + nc - 1

                   IF(MOD(id(1),2).EQ.0) THEN

                      idn = id(1)/2
                      IF(MODE.EQ.0) THEN

                         dz(nr,nc) = dz(nr,nc)&
                              &+ang(1)* (-vz(idn,nhfi(ic),nhff(ir)))*sp_a(lli(ic)+1,mr,mc)

                      ELSE

                         dz(nr,nc) = dz(nr,nc)&
                              &+ang(1)* vz(idn,nhfi(ic),nhff(ir)) * sp_a(lli(ic)+1,mr,mc)
                      ENDIF

                   ELSE

                      idn=(id(1)+1)/2
                      dz(nr,nc) = dz(nr,nc) &
                           &+ang(1) * vz(idn,nhff(ir),nhfi(ic)) * sp_a(lli(ic)+1,mr,mc)

                   END IF

                END DO
             END DO

!!!........................................

          CASE(2)
             DO nr = 1, ndf(ir)

                mr=nminf(ir)+nr-1

                DO nc = 1, ndi(ic)

                   mc=nmini(ic)+nc-1

                   IF(MOD(id(2),2).EQ.0) THEN

                      idn=id(2)/2
                      IF(MODE.EQ.0) THEN


                      dz(nr,nc) = dz(nr,nc)&
                           & + phase*ang(2)*(-vz(idn,mc,nhff(ir)))*sp_a(lhfi(ic)+1,mr,nhfi(ic))
                      ELSE
                      dz(nr,nc)=dz(nr,nc)&
                           &+phase*ang(2)*vz(idn,mc,nhff(ir))*sp_a(lhfi(ic)+1,mr,nhfi(ic))
                      ENDIF
  
                   ELSE
                      idn=(id(2)+1)/2
                      dz(nr,nc) = dz(nr,nc)&
                           &+phase*ang(2)*vz(idn,nhff(ir),mc)*sp_a(lhfi(ic)+1,mr,nhfi(ic))
                   END IF
                END DO
             END DO

!!!.................................

          CASE(3)
             DO  nc = 1, ndi(ic)

                mc = nmini(ic) + nc - 1

                IF(mc.EQ.nhff(ir)) THEN

                   DO  nr=1,ndf(ir)

                      mr = nminf(ir)+nr-1

                      dz(nr,nc) = dz(nr,nc)+ phase * ang(3) * vzba(id(3),mr,nhfi(ic))

                   ENDDO
                   
                ENDIF

                ENDDO

!!!..................................
                CASE(4)  
                   IF(nhff(ir).EQ.nhfi(ic)) THEN

                      DO  nr=1,ndf(ir)
                         
                         mr=nminf(ir)+nr-1
                         DO  nc=1,ndi(ic)
                            
                            mc=nmini(ic)+nc-1

                            dz(nr,nc) = dz(nr,nc)+ang(4)* vzba(id(4),mr,mc)
                            
                         ENDDO
                      ENDDO

                   ENDIF

                END SELECT

             ENDIF

20           CONTINUE


           END SUBROUTINE cal_nlsp
!!!#######################################################################
!!$    subroutine cal_nlsp(ir, ic, id, ang, phase) 
!!$      use configuration
!!$      implicit none
!!$      integer i, nr, nc, mr, mc, ir, ic, idn
!!$      integer, dimension(4) :: id
!!$      real(8), dimension(4) :: ang
!!$      real(8) phase
!!$      do 20 i=1,4
!!$      select case(id(i).eq.0)
!!$        case(.false.)
!!$        select case(i)
!!$          case(1)
!!$          do nr=1,ndf(ir)
!!$            mr=nminf(ir)+nr-1
!!$            do nc=1,ndi(ic)
!!$              mc=nmini(ic)+nc-1
!!$              if(mod(id(1),2).eq.0) then
!!$                idn=id(1)/2
!!$                dz(nr,nc)=dz(nr,nc)+ang(1)*& 
!!$                  &(-vz(idn,nhfi(ic),nhff(ir)))*sp_a(lli(ic)+1,mr,mc)
!!$              else
!!$                idn=(id(1)+1)/2
!!$                dz(nr,nc)=dz(nr,nc)+ang(1)*& 
!!$                   &vz(idn,nhff(ir),nhfi(ic))*sp_a(lli(ic)+1,mr,mc)
!!$              end if
!!$            end do
!!$          end do
!!$
!!$! ---------------------------------------------------------------
!!$          case(2)
!!$          do nr=1, ndf(ir)
!!$            mr=nminf(ir)+nr-1
!!$            do nc=1,ndi(ic)
!!$              mc=nmini(ic)+nc-1
!!$              if(mod(id(2),2).eq.0) then
!!$                idn=id(2)/2
!!$                dz(nr,nc) = dz(nr,nc) + phase*ang(2)&
!!$                 &*(-vz(idn,mc,nhff(ir)))*sp_a(lhfi(ic)+1,mr,nhfi(ic))
!!$              else
!!$                idn=(id(2)+1)/2
!!$                dz(nr,nc) = dz(nr,nc) + phase*ang(2)*&
!!$                   &vz(idn,nhff(ir),mc)*sp_a(lhfi(ic)+1,mr,nhfi(ic))
!!$              end if
!!$            end do
!!$          end do
!!$! ---------------------------------------------------------------
!!$          case(3)
!!$          do 33 nc=1,ndi(ic)
!!$            mc=nmini(ic)+nc-1
!!$            select case(mc .ne. nhff(ir)) 
!!$              case(.false.)
!!$              do 34 nr=1,ndf(ir)
!!$                mr=nminf(ir)+nr-1
!!$                dz(nr,nc)=dz(nr,nc)+phase*ang(3)*&
!!$                   &vzba(id(3),mr,nhfi(ic))
!!$   34         continue
!!$            end select
!!$   33     continue
!!$!-------------------------------------------------------------------
!!$          case(4)  
!!$          select case(nhff(ir) .ne. nhfi(ic)) 
!!$            case(.false.)
!!$            do 35 nr=1,ndf(ir)
!!$              mr=nminf(ir)+nr-1
!!$              do 35 nc=1,ndi(ic)
!!$                mc=nmini(ic)+nc-1
!!$                dz(nr,nc)=dz(nr,nc)+ang(4)*&
!!$                   &vzba(id(4),mr,mc)
!!$   35       continue
!!$          end select
!!$        end select
!!$      end select
!!$ 20 continue
!!$    end subroutine cal_nlsp
!!$!----------------------------------------------------------------------


!!!#############################################

           SUBROUTINE cal_spnl(ir, ic, id, ang, phase, mode) 

             USE configuration
             IMPLICIT NONE
             INTEGER i, nr, nc, mr, mc, ir, ic, idn, mode
             INTEGER, DIMENSION(4) :: id
             REAL(8), DIMENSION(4) :: ang
             REAL(8) phase

!!!.................................................

             DO 20 i = 1, 4


!!         SELECT CASE(id(i).EQ.0)
!!         CASE(.FALSE.)


                IF(id(i).NE.0)  THEN

!!                   WRITE(*,*) ' CAL_SPNL, ---> MODE, i, id(i) ', MODE,i, id(i) 

                   SELECT CASE(i)

                   CASE(1)
                      DO nr = 1, ndf(ir)

                         mr = nminf(ir) + nr - 1
                         
                         DO nc = 1, ndi(ic)

                            mc = nmini(ic) + nc - 1

                            IF(MOD(id(1),2).EQ.0) THEN

                               idn = id(1)/2
                               IF(MODE.EQ.0) THEN

                                  dz(nr,nc) = dz(nr,nc) + ang(1)*& 
                                       &(-vz(idn,nhfi(ic),nhff(ir)))*sp_a(llf(ir)+1,mc,mr)

                               ELSE
                                  dz(nr,nc) = dz(nr,nc) + ang(1)*& 
                                       &vz(idn,nhfi(ic),nhff(ir)) * sp_a(llf(ir)+1,mc,mr)
                               ENDIF

                            ELSE

                               idn = (id(1)+1)/2
                               dz(nr,nc) = dz(nr,nc) + ang(1) * & 
                                    &vz(idn,nhff(ir),nhfi(ic)) * sp_a(llf(ir)+1,mc,mr)

                            END IF

                         END DO
                      END DO

!!!................................................
                   CASE(2)

                      DO nr = 1, ndf(ir)
                         mr = nminf(ir) + nr - 1
                         
                         DO nc=1,ndi(ic)

                            mc = nmini(ic) + nc - 1

!!!                     SELECT CASE(mr .EQ. nhfi(ic)) 
!!!                     CASE(.TRUE.)

                            IF(mr.EQ.nhfi(ic) ) THEN

                               IF(MOD(id(2),2).EQ.0) THEN

                                  idn = id(2) - 1
                               ELSE
                                  idn = id(2) + 1
                               END IF

                               IF(MODE.EQ.0) THEN
!!!                                  WRITE(*,*) ' dz, ang(1), vz, sp_a',dz(nr,nc), ang(2),vzba(idn,mc,nhff(ir))
                                  dz(nr,nc) = dz(nr,nc) + phase*ang(2) *(-vzba(idn,mc,nhff(ir)))
                               ELSE

                                  dz(nr,nc) = dz(nr,nc) + phase*ang(2) * vzba(idn,mc,nhff(ir))

                               ENDIF
                            ENDIF

!!!                  END SELECT

                         END DO
                      END DO
!!!................................................
                   CASE(3)

                      DO  nc = 1, ndi(ic)
                         mc = nmini(ic) + nc - 1

                         DO  nr = 1, ndf(ir)
                            mr = nminf(ir) + nr - 1

                            IF(MOD(id(3),2).EQ.0) THEN

                               idn = id(3)/2
                               IF(MODE.EQ.0) THEN

                                  dz(nr,nc)=dz(nr,nc)+phase*ang(3)*&
                                       &(-vz(idn,nhfi(ic),mr))*sp_a(lhff(ir)+1,mc,nhff(ir))
                               ELSE
                                  dz(nr,nc)=dz(nr,nc)+phase*ang(3)*&
                                       &vz(idn,nhfi(ic),mr) * sp_a(lhff(ir)+1,mc,nhff(ir))

                               ENDIF

                            ELSE
                               idn = ( id(3) + 1 ) /2
                               dz(nr,nc) = dz(nr,nc) + phase*ang(3)*&
                                    &vz(idn,mr,nhfi(ic))*sp_a(lhff(ir)+1,mc,nhff(ir))
                            END IF

                         ENDDO
                      ENDDO
!!!......................................................

                   CASE(4)
  
!!!                     SELECT CASE(nhff(ir) .EQ. nhfi(ic)) 
!!!                     CASE(.TRUE.)

                      IF(nhff(ir) .EQ. nhfi(ic))  THEN

                         DO  nr = 1, ndf(ir)

                            mr = nminf(ir) + nr - 1
                            DO  nc = 1, ndi(ic)
                               mc = nmini(ic) + nc - 1

                               IF(MOD(id(4),2).EQ.0) THEN

                                  idn = id(4) - 1
                               ELSE
                                 
                                  idn = id(4) + 1
                              END IF
                              
                              IF(MODE.EQ.0) THEN

                                 dz(nr,nc)  = dz(nr,nc) + ang(4) * ( -vzba(idn,mc,mr) )

                              ELSE

                                 dz(nr,nc)  = dz(nr,nc) + ang(4) * vzba(idn,mc,mr)
                              ENDIF

                           ENDDO
                        ENDDO

!!!                              END SELECT
                     ENDIF

                  END SELECT
!!!                        END SELECT
!!!               WRITE(*,*) ' dz, ang(4), vzba,',dz(1,1), ang(4), vzba(1,1,1) 
               ENDIF

20             CONTINUE

             END SUBROUTINE cal_spnl

!!!########################################################

    SUBROUTINE cal_spsp(ir, ic, id, ang, phase, mode) 

      USE configuration
      IMPLICIT NONE
      INTEGER i, nr, nc, mr, mc, ir, ic, idn, mode
      INTEGER, DIMENSION(4) :: id
      REAL(8), DIMENSION(4) :: ang
      REAL(8) phase
!!!...........................................
!!!    WRITE(*,*) ' CAL_SPSP '
      DO 20 i = 1, 4

!!!         SELECT CASE(id(i).EQ.0)
!!!         CASE(.FALSE.)
         IF(id(i).NE.0) THEN

            SELECT CASE(i)

            CASE(1)
               DO nr = 1, ndf(ir)

                  mr = nminf(ir) + nr - 1

                  DO nc = 1, ndi(ic)

                     mc = nmini(ic) + nc - 1

                     IF(MOD(id(1),2).EQ.0) THEN

                        idn = id(1)/2
                        IF(MODE.EQ.0) THEN

                        dz(nr,nc) = dz(nr,nc) + ang(1)*& 
                             &(-vz(idn,nhfi(ic),nhff(ir)))*over(lli(ic)+1,mr,mc)
                        ELSE

                        dz(nr,nc) = dz(nr,nc)+ang(1)*& 
                             &vz(idn,nhfi(ic),nhff(ir))*over(lli(ic)+1,mr,mc)

                        ENDIF

                     ELSE
                        idn = ( id(1)+1 ) / 2

                        dz(nr,nc) = dz(nr,nc) + ang(1) * & 
                             &vz(idn,nhff(ir),nhfi(ic))*over(lli(ic)+1,mr,mc)

                     END IF

                  END DO
               END DO
               
!!!...........................................

            CASE(2)
               DO nr = 1,  ndf(ir)

                  mr = nminf(ir) + nr - 1

                  DO nc = 1, ndi(ic)

                     mc = nmini(ic) + nc - 1

                     IF(MOD(id(2),2).EQ.0) THEN
                        idn = id(2) - 1
                     ELSE
                        idn = id(2) + 1
                     END IF
                     
                     IF(MODE.EQ.0) THEN

                        dz(nr,nc) = dz(nr,nc) + phase*ang(2)&
                             &*(-vzba(idn,mc,nhff(ir)))*sp_a(lhfi(ic)+1,mr,nhfi(ic))
                     ELSE
                        dz(nr,nc) = dz(nr,nc) + phase*ang(2)&
                             & * vzba(idn,mc,nhff(ir)) * sp_a(lhfi(ic)+1,mr,nhfi(ic))
                      ENDIF

                   END DO
               END DO

!!!..............................................

            CASE(3)

               DO  nc = 1, ndi(ic)

                  mc = nmini(ic) + nc - 1

                  DO  nr = 1, ndf(ir)

                     mr = nminf(ir) + nr - 1

                     dz(nr,nc) = dz(nr,nc) + phase*ang(3)*&
                          &vzba(id(3),mr,nhfi(ic))*sp_a(lhff(ir)+1,mc,nhff(ir))

                  ENDDO
               ENDDO
!!!..............................................

            CASE(4)  
               
!!!               SELECT CASE(nhff(ir) .NE. nhfi(ic)) 
!!!               CASE(.FALSE.)
               IF( nhff(ir).EQ.nhfi(ic) ) THEN

                  DO  nr = 1, ndf(ir)

                     mr = nminf(ir) + nr - 1

                     DO  nc = 1, ndi(ic)

                        mc = nmini(ic) + nc - 1

                        IF(MOD(id(4),2).EQ.0) THEN

                           idn = id(4)/2
                           IF(MODE.EQ.0) THEN

                              dz(nr,nc) = dz(nr,nc) + ang(4) * (-vzbb(idn,mc,mr))
                           ELSE

                              dz(nr,nc) = dz(nr,nc) + ang(4) * vzbb(idn,mc,mr)
                           ENDIF

                        ELSE
                           idn = ( id(4) + 1 ) / 2
                           dz(nr,nc) = dz(nr,nc) + ang(4) * vzbb(idn,mr,mc)
                        END IF

                     ENDDO
                  ENDDO

                  END IF
                  
               END SELECT

               ENDIF

20             CONTINUE


             END SUBROUTINE cal_spsp
!!!!
!!!#######################################################################


!                 SUBROUTINE cal_nlnl(ir, ic, id, ang, phase)
 
             SUBROUTINE cal_nlnl(ir, ic, id, ang, phase, mode) 
               USE configuration
               IMPLICIT NONE
                                    
               INTEGER  mode
               INTEGER i, nr, nc, mr, mc, ir, ic, idn
               INTEGER, DIMENSION(4) :: id
               REAL(8), DIMENSION(4) :: ang
               REAL(8) phase
                   
!!!..............................
!!!               WRITE(*,*) ' CAL_NLNL '                                    
               DO 20 i = 1, 4

                  IF(id(i).NE.0) THEN
                     
                     SELECT CASE(i)
                            
                     CASE(1)

                        DO nr=1,ndf(ir)
                                                
                           mr=nminf(ir)+nr-1

                           DO nc=1,ndi(ic)

                              mc=nmini(ic)+nc-1

!!!..........
                              IF(mc.EQ.mr) THEN

                                 IF(MOD(id(1),2).EQ.0) THEN

                                    idn=id(1)/2
                                        
                                    IF(MODE.EQ.0) THEN
                                       dz(nr,nc)= dz(nr,nc)+ang(1)*(-vz(idn,nhfi(ic),nhff(ir)))
                                       dz(nr,nc) = 1.0D+00
                                    ELSE
                                       dz(nr,nc)= dz(nr,nc)+ang(1) * vz(idn,nhfi(ic),nhff(ir))
                                    ENDIF
                                    
                                 ELSE
                                       
                                    idn=(id(1)+1)/2
                                    dz(nr,nc)= dz(nr,nc)+ang(1) * vz(idn,nhff(ir),nhfi(ic))

                                 END IF

                              ENDIF
!!....................

                           END DO
                        END DO
                            
                            
                     CASE(2)

                        DO nr=1, ndf(ir)

                           mr=nminf(ir)+nr-1

                           IF(mr.EQ.nhfi(ic)) THEN
                                                   
                              DO nc=1,ndi(ic)
                                 mc=nmini(ic)+nc-1

                                 IF(MOD(id(2),2).EQ.0) THEN

                                    idn=id(2)/2
                                    IF(MODE.EQ.0) THEN
                                      
                                       dz(nr,nc) = dz(nr,nc) + phase * ang(2) *(-vz(idn,mc,nhff(ir)))
                                       dz(nr,nc) = 1.0D+00
                                    ELSE
                                       dz(nr,nc) = dz(nr,nc) + phase * ang(2) * vz(idn,mc,nhff(ir))
                                    ENDIF

                                 ELSE
                                    idn = ( id(2) + 1 ) / 2
                                    dz(nr,nc) = dz(nr,nc) + phase * ang(2) * vz(idn,nhff(ir),mc)
                                 END IF
                                                      

                              END DO

                           ENDIF
                                                
                        END DO
!!!..............................................                            
                     CASE(3)
                        DO 33 nc=1,ndi(ic)
                           mc=nmini(ic)+nc-1
                           
                           IF(mc.EQ.nhff(ir)) THEN

                              DO 34 nr=1,ndf(ir)

                                 mr=nminf(ir)+nr-1

                                 IF(MOD(id(3),2).EQ.0) THEN

                                    idn=id(3)/2
                                    IF(MODE.EQ.0) THEN  
                                       dz(nr,nc) = dz(nr,nc) + phase*ang(3)*(-vz(idn,nhfi(ic),mr))
                                       dz(nr,nc) = 1.0D+00
                                    ELSE
                                       dz(nr,nc) = dz(nr,nc) + phase*ang(3)*vz(idn,nhfi(ic),mr)
                                    ENDIF
                                    
                                 ELSE

                                    idn = (id(3)+1)/2     
                                    dz(nr,nc) = dz(nr,nc) + phase*ang(3)*vz(idn,mr,nhfi(ic))
                                  
                                 END IF
34                               CONTINUE

                              ENDIF
                              
33                            CONTINUE
!!!.....................................................                                  
                            CASE(4)
                               
                               IF(nhff(ir).EQ.nhfi(ic)) THEN
!                           
                                  DO  nr = 1, ndf(ir)

                                     mr = nminf(ir) + nr - 1

                                     DO  nc = 1, ndi(ic)

                                        mc = nmini(ic) + nc - 1

                                        IF(MOD(id(i),2).EQ.0) THEN

                                           IF(MODE.EQ.0) THEN  
                                              
                                              dz(nr,nc) = dz(nr,nc) + ang(i) * (-vz(idn,mc,mr))
                                              dz(nr,nc) = 1.0D+00
                                           ELSE
                                              dz(nr,nc) = dz(nr,nc) + ang(i)*vz(idn,mc,mr)
                                           ENDIF

                                           idn=id(i)/2
                                        ELSE

                                           idn=(id(i)+1)/2
                                           dz(nr,nc) = dz(nr,nc) + ang(i) * vz(idn,mr,mc)

                                        END IF

                                     ENDDO
                                  ENDDO

                               ENDIF
                            END SELECT
                         ENDIF

20                       CONTINUE

                       END SUBROUTINE cal_nlnl
                     END MODULE dz_value          
!!!#######################################################################
!!$    subroutine cal_nlnl(ir, ic, id, ang, phase) 
!!$      use configuration
!!$      implicit none
!!$      integer i, nr, nc, mr, mc, ir, ic, idn
!!$      integer, dimension(4) :: id
!!$      real(8), dimension(4) :: ang
!!$      real(8) phase
!!$      do 20 i=1,4
!!$      select case(id(i).eq.0)
!!$        case(.false.)
!!$        select case(i)
!!$          case(1)
!!$          do nr=1,ndf(ir)
!!$            mr=nminf(ir)+nr-1
!!$            do nc=1,ndi(ic)
!!$              mc=nmini(ic)+nc-1
!!$              select case(mc.eq.mr)
!!$              case(.true.)
!!$                 if(mod(id(1),2).eq.0) then
!!$                    idn=id(1)/2
!!$                    dz(nr,nc)=dz(nr,nc)+ang(1)*& 
!!$                    &(-vz(idn,nhfi(ic),nhff(ir)))
!!$                 else
!!$                    idn=(id(1)+1)/2
!!$                    dz(nr,nc)=dz(nr,nc)+ang(1)*& 
!!$                    &vz(idn,nhff(ir),nhfi(ic))
!!$              end if
!!$              end select
!!$            end do
!!$          end do
!!$
!!$! ---------------------------------------------------------------
!!$          case(2)
!!$          do nr=1, ndf(ir)
!!$            mr=nminf(ir)+nr-1
!!$            select case(mr.eq.nhfi(ic))
!!$            case(.true.)
!!$               do nc=1,ndi(ic)
!!$                  mc=nmini(ic)+nc-1
!!$                  if(mod(id(2),2).eq.0) then
!!$                     idn=id(2)/2
!!$                     dz(nr,nc) = dz(nr,nc) + phase*ang(2)&
!!$                     &*(-vz(idn,mc,nhff(ir)))
!!$                  else
!!$                     idn=(id(2)+1)/2
!!$                     dz(nr,nc) = dz(nr,nc) + phase*ang(2)*&
!!$                     &vz(idn,nhff(ir),mc)
!!$                  end if
!!$               end do
!!$            end select
!!$          end do
!!$! ---------------------------------------------------------------
!!$          case(3)
!!$          do 33 nc=1,ndi(ic)
!!$            mc=nmini(ic)+nc-1
!!$            select case(mc .eq. nhff(ir)) 
!!$              case(.true.)
!!$              do 34 nr=1,ndf(ir)
!!$                mr=nminf(ir)+nr-1
!!$                if(mod(id(3),2).eq.0) then
!!$                   idn=id(3)/2
!!$                   dz(nr,nc)=dz(nr,nc)+phase*ang(3)*&
!!$                   &(-vz(idn,nhfi(ic),mr))
!!$                else
!!$                   idn=(id(3)+1)/2
!!$                   dz(nr,nc) = dz(nr,nc) + phase*ang(3)*&
!!$                   &vz(idn,mr,nhfi(ic))
!!$                end if
!!$   34         continue
!!$            end select
!!$   33     continue
!!$!-------------------------------------------------------------------
!!$          case(4)  
!!$          select case(nhff(ir) .eq. nhfi(ic)) 
!!$            case(.true.)
!!$            do 35 nr=1,ndf(ir)
!!$              mr=nminf(ir)+nr-1
!!$              do 35 nc=1,ndi(ic)
!!$                mc=nmini(ic)+nc-1
!!$                if(mod(id(i),2).eq.0) then
!!$                   idn=id(i)/2
!!$                   dz(nr,nc)=dz(nr,nc)+ang(i)*(-vz(idn,mc,mr))
!!$                else
!!$                   idn=(id(i)+1)/2
!!$                   dz(nr,nc)=dz(nr,nc)+ang(i)*vz(idn,mr,mc)
!!$                end if
!!$
!!$   35       continue
!!$
!!$          end select
!!$        end select
!!$      end select
!!$
!!$ 20 continue
!!$!-------------------------------------------------------------------
!!$    end subroutine cal_nlnl
!!$!

