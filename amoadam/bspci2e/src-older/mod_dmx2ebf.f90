MODULE jl
  IMPLICIT NONE
  INTEGER j1, j2, ls
END MODULE jl

MODULE wf_files
  USE PRECISION, ONLY: dpk
  USE ioroutines

  IMPLICIT NONE
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:) :: vz, vzba, sp_a

CONTAINS
!!!#############################################################

  SUBROUTINE read_dp(nsx, ltotalMax,mode)
    !
    USE utils, ONLY: dfile, datafile
    !
    IMPLICIT NONE
    !arg!
    integer          :: nsx
    integer          :: ltotalmax
    integer          :: mode
    !loc!
    CHARACTER(len=6) :: gauge
    !
    INTEGER          :: m, n, li, lf, nwfi, nwff, ir, ic
    !io!    
    INTEGER         :: nout, nbin
    INTEGER         :: i,j

!!!......................
    WRITE(*, '(a60)') 'subroutine::read_dp in.'
    WRITE(*, '(a60)') 'dipole matrix elements for target states.'

    !
    nbin   = 2
    nout   = 16


    IF(mode==0) THEN
       gauge = 'v'
    ELSE IF(mode==1) THEN
       gauge = 'l'
    ELSE IF(mode==2) THEN
       gauge = 'a'
    ENDIF
       
       !
    WRITE(*,'(a60,i10)') ' lmax = ' , ltotalMax
    WRITE(*,'(a60,i10)') ' nof basis functions nb = ', nsx
    WRITE(*,'(a60,a10)') '                  gauge = ', gauge


    ALLOCATE(vz(ltotalMax, nsx, nsx))             ! <P_nl|T|P_ml+1>, l=0,...,lmax-1

    DO i = 1, ltotalMax  

!       CALL HDMX1EFILE(NDMX1E, i-1, mode)

       CALL dfile(nbin,i-1,i,'d1e-',gauge)
       READ(nbin) mode
       READ(nbin) li, lf, nwfi, nwff
       READ(nbin) ( ( vz(i,ir,ic), ir = 1, nwff), ic = 1, nwfi)       
       CLOSE(nbin)
         
       WRITE(*,'(a60,5i5)') '(mode,la,lb,na,nb) = ', mode, li,lf,nwfi,nwff
    END DO

    WRITE(*,'(a60)') 'dipoles  <P|T|P> read.'

    ALLOCATE( vzba(2*ltotalMax, nsx, nsx) )

    j = 0
    get_overlap_dipoles: DO i = 1, ltotalMax   

!       CALL RWFBFILE(nbin, i-1, i, mode)
       CALL dfile(nbin,i-1,i,'overlap-d-',gauge)

       READ(nbin) mode
       READ(nbin) li,lf,nwfi,nwff
       READ(nbin) ((vzba(i + j, ir, ic), ir = 1, nwff), ic = 1, nwfi)

       CLOSE(nbin)


!       CALL RWFBFILE(NRWFB, i, i-1, MODE)
       CALL dfile(nbin,i,i-1,'overlap-d-',gauge)

       READ(nbin) mode
       READ(nbin) li,lf,nwfi,nwff
       READ(nbin) ((vzba(i+j+1, ir, ic), ir = 1, nwff), ic = 1, nwfi)

       j = j + 1

       CLOSE(nbin)
       
       WRITE(*,'(a60,5i5)') '(mode,la,lb,na,nb) = ', mode, li,lf,nwfi,nwff
    END DO get_overlap_dipoles

    WRITE(*,'(a60)') 'overlap dipoles  <P|T|B> read.'



    ALLOCATE(sp_a(ltotalMax, nsx, nsx))          ! <P_nl|B_j> overlaps

    get_overlaps:DO i = 1, ltotalMax

!       CALL RBBFILE(nbin, i-1)
       CALL datafile(nbin,i-1,'overlap-1-') 

       READ(nbin) li,nwfi,nwff
       READ(nbin) ((sp_a(li+1,ir,ic), ir=1,nwff), ic=1,nwfi)
       CLOSE(nbin)

       WRITE(*,'(a60,3i5)') '(la,na,j) = ', li,nwfi,nwff
    END DO get_overlaps

    WRITE(*,'(a60)') 'overlap <P|B> read.'


    WRITE(*, '(a60)') 'subroutine::read_dp out.'
  END SUBROUTINE read_dp

END MODULE wf_files
      

MODULE configuration
  IMPLICIT NONE
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhfi,lhfi,lli,nmini,nmxi, ndi
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nhff,lhff,llf,nminf,nmxf, ndf
  INTEGER, ALLOCATABLE, DIMENSION(:) :: isf, nolf, idcs
CONTAINS
!sub
  SUBROUTINE config_space(ncsi, ncsf)  
    IMPLICIT NONE
    INTEGER ncsi, ncsf

    ALLOCATE( nhfi (ncsi) )
    ALLOCATE( lhfi (ncsi) )
    ALLOCATE( lli  (ncsi) )
    ALLOCATE( nmini(ncsi) )
    ALLOCATE( nmxi (ncsi) )
    ALLOCATE( ndi  (ncsi) )

    ALLOCATE( nhff (ncsf) )
    ALLOCATE( lhff (ncsf) )
    ALLOCATE( llf  (ncsf) )
    ALLOCATE( nminf(ncsf) )
    ALLOCATE( nmxf (ncsf) )
    ALLOCATE( ndf  (ncsf) )


    ALLOCATE( nolf(ncsf) )
    ALLOCATE( isf (ncsf) )
    ALLOCATE( idcs(ncsf) )

    
  END SUBROUTINE config_space
END MODULE configuration

MODULE dz_value
  USE wf_files
  IMPLICIT NONE
  REAL(8), ALLOCATABLE, DIMENSION(:,:) :: dz
CONTAINS
  !S
  !S
  !S
  SUBROUTINE cal_nlsp(ir, ic, id, ang, phase, mode)
    !
    USE PRECISION, only: dpk
    !
    USE configuration
    !
    IMPLICIT NONE
    INTEGER MODE
    INTEGER i, nr, nc, mr, mc, ir, ic, idn
    INTEGER, DIMENSION(4) :: id          
    REAL(dpk), DIMENSION(4) :: ang
    REAL(dpk) phase
    !

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

                      dz(nr,nc) = dz(nr,nc)+ang(1)* (-vz(idn,nhfi(ic),nhff(ir)))*sp_a(lli(ic)+1,mr,mc)

                      ELSE

                      dz(nr,nc) = dz(nr,nc)+ang(1)* vz(idn,nhfi(ic),nhff(ir)) * sp_a(lli(ic)+1,mr,mc)

                      ENDIF

                   ELSE

                      idn=(id(1)+1)/2
                      dz(nr,nc)= dz(nr,nc)+ang(1) * vz(idn,nhff(ir),nhfi(ic)) * sp_a(lli(ic)+1,mr,mc)

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
                      dz(nr,nc)=dz(nr,nc)+phase*ang(2)*(-vz(idn,mc,nhff(ir)))*sp_a(lhfi(ic)+1,mr,nhfi(ic))
                      ELSE
                      dz(nr,nc)=dz(nr,nc)+phase*ang(2)*vz(idn,mc,nhff(ir))*sp_a(lhfi(ic)+1,mr,nhfi(ic))
                      ENDIF
  
                   ELSE
                      idn=(id(2)+1)/2
                      dz(nr,nc) = dz(nr,nc)+phase*ang(2)*vz(idn,nhff(ir),mc)*sp_a(lhfi(ic)+1,mr,nhfi(ic))
                   END IF
                END DO
             END DO

!!!.................................

          CASE(3)
             DO 33 nc=1,ndi(ic)

                mc=nmini(ic)+nc-1

                IF(mc.EQ.nhff(ir)) THEN

                   DO 34 nr=1,ndf(ir)

                      mr=nminf(ir)+nr-1

                      dz(nr,nc)= dz(nr,nc)+ phase * ang(3) * vzba(id(3),mr,nhfi(ic))

34                    CONTINUE

                   ENDIF

33                 CONTINUE

!!!..................................
                CASE(4)  
                   IF(nhff(ir).EQ.nhfi(ic)) THEN

                      DO 35 nr=1,ndf(ir)
                         
                         mr=nminf(ir)+nr-1
                         DO 35 nc=1,ndi(ic)
                            
                            mc=nmini(ic)+nc-1

                            dz(nr,nc) = dz(nr,nc)+ang(4)* vzba(id(4),mr,mc)

35                          CONTINUE


                         ENDIF
                      END SELECT
                   ENDIF

20                 CONTINUE

                 END SUBROUTINE cal_nlsp
!!!#######################################################################
                 
                 SUBROUTINE cal_nlnl(ir, ic, id, ang, phase, mode) 
                   !mod!
                   USE PRECISION, only: dpk
                   USE configuration
                   !
                   IMPLICIT NONE
                   !arg!
                   INTEGER                    :: ir
                   INTEGER                    :: ic
                   INTEGER,      DIMENSION(4) :: id
                   REAL(dpk),    DIMENSION(4) :: ang
                   REAL(dpk)                  :: phase
                   INTEGER                    :: mode
                   INTEGER i, nr, nc, mr, mc, idn
                   !exe!

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
                                        ELSE
                                           dz(nr,nc) = dz(nr,nc) + phase * ang(2) * vz(idn,mc,nhff(ir))
                                        ENDIF

                                     ELSE
                                        idn=(id(2)+1)/2
                                        dz(nr,nc) = dz(nr,nc) + phase * ang(2) * vz(idn,nhff(ir),mc)
                                     END IF


                                  END DO

                               ENDIF

                            END DO
                            
                         CASE(3)
                            DO 33 nc=1,ndi(ic)
                               mc=nmini(ic)+nc-1

                               IF(mc.EQ.nhff(ir)) THEN

                                  DO 34 nr=1,ndf(ir)

                                     mr=nminf(ir)+nr-1

                                     IF(MOD(id(3),2).EQ.0) THEN

                                        idn=id(3)/2
                                        IF(MODE.EQ.0) THEN  
                                        dz(nr,nc)=dz(nr,nc)+phase*ang(3)*(-vz(idn,nhfi(ic),mr))
                                        ELSE
                                        dz(nr,nc)=dz(nr,nc)+phase*ang(3)*vz(idn,nhfi(ic),mr)
                                        ENDIF

                                     ELSE
                                        idn=(id(3)+1)/2
                                        dz(nr,nc) = dz(nr,nc) + phase*ang(3)*vz(idn,mr,nhfi(ic))
                                     END IF
34                                   CONTINUE

                                  ENDIF

33                                CONTINUE
                                  
                            CASE(4)
                               
                               IF(nhff(ir).EQ.nhfi(ic)) THEN
!                           
                                  DO 35 nr = 1, ndf(ir)

                                     mr=nminf(ir)+nr-1

                                     DO  35 nc = 1, ndi(ic)

                                        mc=nmini(ic)+nc-1

                                        IF(MOD(id(i),2).EQ.0) THEN
                                           idn=id(i)/2
                                           IF(MODE.EQ.0) THEN  
                                              dz(nr,nc) = dz(nr,nc) + ang(i) * (-vz(idn,mc,mr))
                                           ELSE
                                              dz(nr,nc) = dz(nr,nc) + ang(i)*vz(idn,mc,mr)
                                           ENDIF


                                        ELSE

                                           idn=(id(i)+1)/2
                                           dz(nr,nc) = dz(nr,nc) + ang(i) * vz(idn,mr,mc)

                                        END IF

35                                      CONTINUE

                               ENDIF
                            END SELECT
                         ENDIF

20                       CONTINUE
                                  
                                  
                       END SUBROUTINE cal_nlnl

                       !
                     END MODULE dz_value          
                     !eof         
                     
