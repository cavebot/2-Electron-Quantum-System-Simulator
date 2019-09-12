
!
!  CALCULATES       < P(nb,lb) | d | B(lz) >
!    
!                      d =  r   d = d/dr
!

MODULE dmx_1e_matrix
  
  USE PRECISION,ONLY:dpk
  IMPLICIT NONE
  PUBLIC
  
CONTAINS

  SUBROUTINE dmxbb(NRWFBB, lb, lz, ndim, mode)
    !
    USE PRECISION, ONLY:dpk
    USE param
    USE one_e_matrix
    !
    IMPLICIT NONE
    !
    INTEGER lb, lz, i,j, lt, l, nderiv, m, jj, ndim, mode
    REAL(dpk), DIMENSION(ndim)             :: coef
    REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: dvx
    REAL(dpk)                              ::  t1(ns+nk), t(ns+nk)
    REAL(dpk)                              :: seh, xm, wfx
    REAL(dpk), DIMENSION(nk,nk)            :: db1, db
    REAL(dpk), DIMENSION(nk)               :: xg,wg
    INTEGER                                :: NRWFBB, NOUT
    DATA nderiv/2/
    REAL(DPK)  bvalue
!!!.......................
    

    !
    NOUT = 16
    
1   FORMAT(/2x,' <n1,l1 | D |n2,l2> = ',4i3)
2   FORMAT(2x,1p8e16.8)
3   FORMAT(/2x,'LENGTH'/)
4   FORMAT(/2x,'VELOCITY'/)
    
    
    
    WRITE(*,'(a60)') ' dmxbb in'


!!!........................

  db  = 0.0D+00
  db1 = 0.0D+00
  xg  = 0.0D+00
  wg  = 0.0D+00


  CALL gauss(k,xg,wg)

  ! 
  !  Generate B-spline points for the two partial waves
  !

  CALL mkgrid( n2(lz + 1), k, rmax, rs(lz + 1), t1)
  CALL mkgrid( n2(lb +1),  k, rmax, rs(lb + 1),  t)


  WRITE(*,'(a60)'),' Calculation of T(i, la ; j,lb) = < B_i(r),la| T_G | B_j, lb >  '
  WRITE(*,'(a60)') '  i = 1,.., nb'
  WRITE(*,'(a60,i10)') ' la = ', lb      
  WRITE(*,'(a60)') '  j = 1,.., nb'
  WRITE(*,'(a60,i10)') ' lb = ', lz
  WRITE(*,'(a60,i10)'),' n_b(la) = ', n2(lb+1)
  WRITE(*,'(a60,i10)'),' n_b(lb) = ', n2(lz+1)
  WRITE(*,'(a60)') ' Matrix <B_i, la| T_G | B_j, lb> is stored as  T_G(j,i).'


  ALLOCATE( dvx( n2(lz + 1), n2(lb + 1) ) )

  dvx = 0.0D+00

  DO i = 1, n2(lb + 1)
     DO j = 1, n2(lz + 1)

        coef    = 0.0D+00
        coef(j) = 1.0D+00

        DO lt = max0(k,i), min0( i + k - 1, n2(lb + 1))
                                !        print*, 'lt=', lt
           seh = 0.0d+00

           DO m = 1, k

              xm = ( t(lt + 1 ) - t(lt)) * xg(m) + t(lt)
              
              CALL bsplvd(t, k, xm, lt, db, nderiv)

              wfx = bvalue(t1, coef, n2( lz + 1), k, xm, 0)

              IF(MODE.EQ.0) THEN                 !!! velocity gauge

!!! original
!!!                 SELECT CASE( lb - lz)
                 SELECT CASE( lz - lb)
                 
                 CASE(1)

                    seh = seh + wg(m) * wfx * &
                         &( db( i - lt + k, 2) + db( i - lt + k,1) * dble(lb) /xm)

                 CASE(-1)

                    seh = seh + wg(m) * wfx * &
                         &(db( i - lt + k, 2) - db( i - lt + k, 1)*dble(lb + 1) /xm)
                 END SELECT

              ELSE                                   !!! length gauge

                 seh = seh + wg(m) * wfx * db( i - lt + k, 1 ) * xm

              ENDIF

!!!               write(*,*) xm, wfx, db(j-lt+k,1), db(j-lt+k,2)

           END DO

              dvx(j,i) = dvx(j,i) - ( t(lt+1) - t(lt) ) * seh     

        END DO
     END DO
  END DO


  WRITE(NRWFBB) mode
  WRITE(NRWFBB) lb, lz, n2( lb + 1), n2( lz + 1)
  WRITE(NRWFBB) ( (dvx(j,i), j = 1, n2( lz + 1)), i = 1, n2( lb + 1) )


  DO  i = n2(lb+1), n2(lb+1) 
     WRITE(*,'(3E20.10)') dvx(i,i), dvx(i,i-1), dvx(i-1,i)
  ENDDO

  DEALLOCATE(dvx)

  WRITE(*,'(a60)') ' dmxbb out.'
END SUBROUTINE dmxbb


!
!
!
!  Calculate the  T_(nj) = < B_j |T_G| P_{nl} > matrix
!
!          n = 1,2,...,nb - 2
!          j = 1,2,.. ,nb
!
!          G = length, velocity
!

  
  SUBROUTINE dmx(nrwfb, lb, lz, coef, ndim, mode)
    !mod!    
    USE PRECISION,ONLY: dpk
    USE param
    USE one_e_matrix
    USE ioroutines
    !
    IMPLICIT NONE
    !arg!
    INTEGER                           :: nrwfb
    INTEGER                           :: lb
    INTEGER                           :: lz
    REAL(dpk),   DIMENSION(ndim,ndim) :: coef
    INTEGER                           :: ndim
    INTEGER                           :: mode
    !loc!
    INTEGER                                :: nb, kb
    INTEGER                                :: nderiv, m, jj
    REAL(dpk),             DIMENSION(ndim) :: bcoef
    REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: dvx
    REAL(dpk)                              ::  t1(ns+nk), t(ns+nk)
    REAL(dpk)                              ::  seh, xm, wfx
    REAL(dpk),            DIMENSION(nk,nk) :: db
    REAL(dpk),               DIMENSION(nk) :: xg,wg
    REAL(DPK)                              :: bvalue
    INTEGER                                :: i,j, lt, l
    INTEGER                                :: nout, ifile
    integer                                :: nb_plot
    DATA nderiv/2/
!!!.........................

    WRITE(*,'(a60)') 'dmx in '

    ifile = 10

    nb = n2(lz+1)
    kb = k

    
    NOUT = 16
         
    db = 0.0D+00
    xg = 0.0D+00
    wg = 0.0D+00

    CALL gauss(k,xg,wg)
    CALL mkgrid( n2(lz+1), k, rmax, rs(lz+1), t )  !make grid for lz
    CALL mkgrid( n2(lb+1), k, rmax, rs(lb+1), t1 ) !make grid for lb
        
    WRITE(*,'(a60)') ' Calculation of <P_{n_al_a}| T_G | B_j, lb> '
    WRITE(*,'(a60)') '  j = 1,.., nb'
    WRITE(*,'(a60,i10)') ' lb = ', lz
    WRITE(*,'(a60)') ' n_a = 1,.., nb-2'
    WRITE(*,'(a60,i10)') ' la = ', lb      
    WRITE(*,'(a60)') ' Matrix <P_{n_al_a}| T_G | B_j, lb> is stored as  T_G(j, n_a).'

    ALLOCATE( dvx( n2(lz+1), n2(lb+1) - 2 ) )
    dvx = 0.0D+00
    loop_over_fxd_states: DO i = 1, n2(lb+1) - 2
              
       DO jj = 1, n2(lb+1) - 1
          bcoef(jj) = coef(jj,i)
       END DO
       !              1,2,......nb-1,nb
       ! bcoef(nb) = (0,x,.....,   x,0)
       !
       loop_over_b_splines: DO j = 1, n2(lz+1)
          loop_over_segments: DO lt = max0(k,j), min0(j+k-1,n2(lz+1))

             !   wfx = P_a (r) 
             !
             seh = 0.0D+00
             DO m = 1, k
                
                xm = ( t(lt+1) - t(lt) ) * xg(m) + t(lt)

                CALL bsplvd(t, k, xm, lt, db, nderiv)
         

                wfx = bvalue(t1, bcoef, n2(lb+1), k, xm, 0 )


                IF(MODE.EQ.0) THEN  !velocity gauge

                   SELECT CASE(lz-lb)
                   CASE(1)

!                      seh = seh + wg(m) * wfx * db(j-lt+k,2) 
                      seh = seh + wg(m) * wfx * &
                           &(db(j-lt+k,2) + db(j-lt+k,1) *DBLE(lz) / xm )

!                      seh = seh + wg(m) * wfx * &
!                           &(db(j-lt+k,1) *DBLE(lz) / xm )

                   CASE(-1)
                      seh = seh + wg(m) * wfx * &
                           &(db(j-lt+k,2) - db(j-lt+k,1) * DBLE(lz + 1) / xm )
                   END SELECT

                ELSE               !length gauge
                   
                   seh = seh + wg(m) * wfx * db(j-lt+k, 1) * xm
                   
                ENDIF
                
             END DO


             dvx(j,i) = dvx(j,i) - ( t(lt+1) - t(lt) ) * seh 

!                 IF(j==nb) THEN
!                     WRITE(*,*) j, (kb-1) * rMax * bcoef(nb-2) / ( 2.0d+00* ( rMax - t(nb) ) ), dvx(j,i)
!                     dvx(j,i) = dvx(j,i) + (kb-1) * rMax*bcoef(nb-2) / ( 2.0d+00* ( rMax - t(nb) ) )
!                     WRITE(*,*) j, (kb-1) * rMax * bcoef(nb-2) / ( 2.0d+00* ( rMax - t(nb) ) ), dvx(j,i)
!                     WRITE(*,*) '--------'
!                  ENDIF

          END DO loop_over_segments
       END DO loop_over_b_splines
    END DO loop_over_fxd_states
    
    WRITE(NRWFB) MODE
    WRITE(NRWFB) lb, lz, n2(lb+1) - 2, n2(lz+1)
    WRITE(NRWFB) ( (dvx(j,i), j = 1, n2(lz+1)), i = 1, n2(lb+1)-2 )
    
!!!.............................

    !l.gg      write(3) ((drx(j,i), j=1,nwfz), i=1,nwfb)
         
    WRITE(NOUT,1) lb, lz, n2(lb+1)-2, n2(lz+1)
    DO  j = nb, nb
       WRITE(nout, 2) (dvx(j,i), i=1,n2(lz+1)-2)
    ENDDO
    WRITE(*,'(a60,2i10)')  ' T(n,i) = <P_nl|T|B_i>: '
    WRITE(*,'(a5,1x,6a14)') " i\n","1","2","3","4","5","6"
    DO  i = 1, 6        
       WRITE(*,'(2X,I3,1X,1P8E16.8)') i, (dvx( i,j ), j = 1, 6)
    ENDDO
    WRITE(*,*)'.................................................'
    DO  i = nb-5, nb
       WRITE(*,'(2X,I3,1X,1P8E16.8)') i, (dvx( i,j ), j = 1, 6)
    ENDDO

     IF(lb.EQ.0) THEN
        nb_plot = 4
        OPEN(ifile,file='out/rwfb-bp-mpq.out')
        DO j = 1, ndim
           WRITE(ifile,'(i10,2E25.14,2i5)') j, ABS(dvx(j,nb_plot))**2, dvx(j,nb_plot), nb_plot, lb
        ENDDO
        CLOSE(ifile)
     ENDIF

!            DO  i = 1, 2
!               WRITE(NOUT, 2) (dvx(j,i), j=1,n2(lz+1))
!            ENDDO

    DEALLOCATE(dvx)

    WRITE(*,'(a60)') 'dmx out. '



1   FORMAT(/2x,' < L1 N1 | R | L2, N2 >  = ',4i3)
2   FORMAT(2x,1p8e16.8)
3   FORMAT(/2x,'LENGTH ::::'/) 
4   FORMAT(/2x,'VELOCITY ::::'/)
    

  END SUBROUTINE dmx

  !
  !
  !
  !  Calculate the  O_(nj) = < B_j |P_{nl} > matrix
  !
  !          n = 1,2,...,nb - 2
  !          j = 1,2,.. ,nb
  !

SUBROUTINE overlap( nrbb, lb, coef, ndim)
  !
  USE PRECISION,only:dpk
  USE param
  USE one_e_matrix
  USE ioroutines
  !
  IMPLICIT NONE
  !
  INTEGER lb, i,j, lt, l, nderiv, m, jj, ndim
  INTEGER NOUT, nrbb
  REAL(dpk),        DIMENSION(ndim,ndim) :: coef
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) ::  over
  REAL(dpk)  t(617)
  REAL(dpk)  seg, xm, wfx
  REAL(dpk), DIMENSION(15,15) :: db
  REAL(dpk), DIMENSION(15) :: xg,wg
  DATA nderiv/1/
  REAL(DPK)  bvalue


!!!.....................

1 FORMAT(/2x,'<B_I|R|B_J> = ',4i3)
2 FORMAT(2x,1p8e16.8)
3 FORMAT(/2x,'LENGTH   ::::'/)
4 FORMAT(/2x,'VELOCITY ::::'/)

!!!....................

  NOUT = 16 

  WRITE(*,'(a60)') 'overlap in.'

  !
  db = 0.0D+00
  xg = 0.0D+00
  wg = 0.0D+00


  CALL gauss(k, xg, wg)
  CALL mkgrid(n2(lb+1), k, rmax, rs(lb+1),t)  !make grid for lb
  !

  WRITE(*,'(a60)') ' Calculation of <B_j|P_{nl}> '
  WRITE(*,'(a60)') ' j = 1,.., nb'
  WRITE(*,'(a60)') ' n = 1,.., nb-2'
  WRITE(*,'(a60,i10)') 'lb = ', lb   
  WRITE(*,'(a60,i10)') 'nb = ', n2(lb+1)


  !
  ALLOCATE(over( n2(lb + 1), n2(lb + 1) - 2 ) )

  over  = 0.0D+00

  loop_over_b_splines_1: DO i = 1, n2(lb+1) - 2
     loop_over_b_splines_2: DO j = 1, n2(lb+1)
        loop_over_segments_t: DO lt = max0(k,j), min0(j+k-1,n2(lb+1))

           seg = 0.0D+00

           DO m = 1, k

              xm = ( t(lt+1) - t(lt)) * xg(m) + t(lt)

              CALL bsplvd(t,k,xm,lt,db,nderiv)

              wfx = bvalue(t, coef(:,i), n2(lb+1), k, xm, 0 )

              seg = seg + wg(m) * wfx * db( j - lt + k, 1 )

           END DO

           over(j,i) = over(j,i) + ( t(lt + 1) - t(lt) ) * seg

        END DO loop_over_segments_t
     END DO loop_over_b_splines_2
  END DO loop_over_b_splines_1


  WRITE(NRBB)    lb, n2(lb+1) - 2,  n2(lb+1)
  WRITE(NRBB) ( (over(j,i), j = 1, n2(lb+1) ), i = 1, n2(lb+1) - 2 )

  WRITE(*,'(a60,2i10)')  ' O(n,i) = <P_nl|B_i>: '
  WRITE(*,'(a5,1x,6a14)')" i\n","1","2","3","4","5","6"
  DO  j = 1, 6
     WRITE(*,'(2X,I3,1X,1P8E16.8)') j, ( over( j,i ), i = 1, 6)
  ENDDO
  WRITE(*,*)'.................................................'
  DO  j = n2(lb+1)-5, n2(lb+1)
     WRITE(*,'(2X,I3,1X,1P8E16.8)') j, ( over( j,i ), i = 1, 6)
  ENDDO


 DO i = 1, 2
    WRITE(NOUT, 2) ( over(j, i), j = 1, n2(lb+1) )
 ENDDO

 WRITE(*,'(a60)') 'overlap out.'

 DEALLOCATE(over)
END SUBROUTINE overlap


END MODULE dmx_1e_matrix
!eof
