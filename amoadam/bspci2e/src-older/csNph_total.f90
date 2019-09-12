!
!  add the various partial waves and provides the total cross sections
!  
!  sN = sum_j sN_j
!
!
!
! (note: the various partial waves have different final energies
!        for this reason an interpolation is required for the summation
!
!
!
PROGRAM CSNPH_TOTAL
  !
  use precision, only:dpk
  use csNph_utils, only:linear_interpolation,beta_k,pad
  use units
  IMPLICIT NONE

  !..................... input section 
  INTEGER,          PARAMETER :: nstates = 8000
  CHARACTER(LEN=*), PARAMETER ::  scsnph = "csNph-l"
  CHARACTER(LEN=*), PARAMETER ::    sdat = "dat"
  CHARACTER(LEN=*), PARAMETER ::    spad = "pad"
  CHARACTER(LEN=*), PARAMETER ::     scs = "cs"
  CHARACTER(LEN=*), PARAMETER ::   sbeta = "beta"
  !
  INTEGER                                :: nout,i,j
  INTEGER                                :: ncsfile, nbetafile, npadfile
  INTEGER                                :: nphotons
  INTEGER                                :: nfinal_l, l_initial 
  CHARACTER( LEN = 6  )                  :: sl, sph, sn_pad
  CHARACTER( LEN = 30 )                  :: csfile, betafile, padfile 
  CHARACTER( LEN = 30 ), ALLOCATABLE, DIMENSION(:) :: CS_L_FILE
  INTEGER, ALLOCATABLE,     DIMENSION(:) :: L, NSTATES_L
  REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: WPH, S_AU, S_SI
  REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: D_R,D_I,PHASE,AN
  REAL(DPK), ALLOCATABLE, DIMENSION(:)   :: x,y,yr,yi,wn
  REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: csNph_L,dmx_l_r,dmx_l_i, delta_l
  REAL(DPK), ALLOCATABLE, DIMENSION(:)   :: csNph
  REAL(DPK), ALLOCATABLE, DIMENSION(:)   :: min_xl, max_xl
  REAL(DPK)                              :: w_in, w_fin
  REAL(DPK)                              :: unit_n,cs_n_si
  INTEGER                                :: lp
! pad
  REAL(DPK), ALLOCATABLE, DIMENSION(:)   :: bp
  INTEGER                                :: kp, n_pad
  REAL(dpk)                              :: m_n_2 
  INTEGER                                :: check
  INTEGER                                :: nxn
  CHARACTER(LEN=1)                       ::   sxm  
  CHARACTER(LEN = 100)                   :: argv
!!!!
!.............................  input data section

  ncsfile   = 21
  nbetafile = 22
  npadfile  = 23

!.............................

  CALL GETARG(1, ARGV)
  READ(ARGV,*) nphotons
  CALL GETARG(2, ARGV)
  READ(ARGV,*) n_pad

  NOUT = 16
                                  ! convert number to string
  WRITE(sph,   '(i6)') nphotons   !nof photons 
  WRITE(sn_pad,'(i6)') n_pad      !

!..............................

  UNIT_N = 1.0_dpk

  IF(NPHOTONS.EQ.1) THEN         
     UNIT_N = 1.0e+18_dpk
  ENDIF


  CS_N_SI  =  2.0_dpk * M_PI * ( 2.0_dpk * M_PI * ALPHA )**NPHOTONS & ! cs_au
                &  * T0**(NPHOTONS - 1 ) * A0**(2*NPHOTONS) 

      WRITE(*,*) "#  READ_CS::         cs_n_si = ",cs_n_si

  ! determine number of final angular symmetries

      NFINAL_L = 0
      DO I = NPHOTONS, 0, -2     
         NFINAL_L = NFINAL_L + 1 
      ENDDO
      
      L_INITIAL = I + 2

      ALLOCATE(           L(NFINAL_L) )
      ALLOCATE(   NSTATES_L(NFINAL_L) )
      ALLOCATE( CS_L_FILE(NFINAL_L) ) 
      ALLOCATE(   bp(0:NPHOTONS) )
      
      ! initialize final angular symmetry
      
      DO I = 1, NFINAL_L
         L(I) = L_INITIAL + 2*(I-1)
      ENDDO
      

      WRITE(*,*) "#  READ_CS::         INPUT DATA:              "
      WRITE(*,*) "#"
      WRITE(*,*) "#  READ_CS::                       NPHOTONS = ", NPHOTONS
      WRITE(*,*) "#  READ_CS::     CALCULATED DATA:             "
      WRITE(*,*) "#"
      WRITE(*,*) "#  READ_CS:: NOF FINAL SYMMETRIES  NFINAL_L = ", NFINAL_L
      WRITE(*,*) "#  READ_CS:: INITIAL SYMMETRY     L_INITIAL = ", L_INITIAL
      WRITE(*,*) "#  READ_CS::               FINAL SYMMETRIES   "
      WRITE(*,*) "#  READ_CS:: L = ", L
      WRITE(*,*) "#"
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !
      !  get the names of the files with the partial waves
      !

      DO I = 1, NFINAL_L
         WRITE(SL,'(I6)') L(I)   !convert number to string
         CS_L_FILE(I) = sdat//"/"//SCSNPH//TRIM(ADJUSTL(SL))//"."//sdat
      ENDDO
      
      WRITE(*,'(10a20)') (CS_L_FILE(I), I = 1, NFINAL_L)
      
!!!.......  read   ( energy, cs_si, cs_au, dr, di, phase )
      
      ALLOCATE(  WPH(NFINAL_L, NSTATES) )
      ALLOCATE( S_SI(NFINAL_L, NSTATES) )
      ALLOCATE( S_AU(NFINAL_L, NSTATES) )
      ALLOCATE(  D_I(NFINAL_L, NSTATES) )
      ALLOCATE(  D_R(NFINAL_L, NSTATES) )
      ALLOCATE(  AN(NFINAL_L, NSTATES) )
      ALLOCATE(  PHASE(NFINAL_L, NSTATES) )
      
      DO I = 1, NFINAL_L
         
         OPEN(UNIT=I+10, FILE = CS_L_FILE(I))         
         READ(I+10,*)    sxm, nstates_l(i)
         
         DO j =  1, nstates_l(i)
            
            READ(I+10,'(6(1X,1PE15.8))') &
                 &   WPH(I,J)                 &
                 &  ,S_SI(I,J)                &
                 &  ,D_R(I,J)                 &
                 &  ,D_I(I,J)                 &
                 &  ,AN(I,J)                  &
                 &  ,PHASE(I,J)                   
         ENDDO
         
         CLOSE(I+10)
         
      ENDDO
      
!!!
      !  do j = 1, nstates
      !     write(*,*)(phase(i,j), i =1,nfinal_l)
      !  enddo
      !  stop
!!!........ Make the linear interpolation between the partial waves
!!!  the input values to be interpolated are those of symmetry lp
      
      lp = 1  !?

      ALLOCATE(       wn( nstates_l(lp) ) )
      ALLOCATE( csnph_L( nfinal_l, nstates_l(lp) ) )
      ALLOCATE(   dmx_l_r( nfinal_l, nstates_l(lp) ) )
      ALLOCATE(   dmx_l_i( nfinal_l,nstates_l(lp) ) )
      ALLOCATE(   delta_l( nfinal_l,nstates_l(lp) ) )
      ALLOCATE(     csNph( nstates_l(lp) ) )
      
      wn = wph(lp,1:nstates_l(lp))
      nxn = SIZE( wn )                   ! nof_states_of_symmetry l
      
      WRITE(*,*) "#  csNph_int::           nof sampled points = ", nxn
      
      
!!!
      !xxxxxxxxxxx                 interpolate the cross sections first
!!!  
      
      dmx_l_r =  0.0_dpk
      dmx_l_i =  0.0_dpk
      delta_l =  0.0_dpk
      csNph   =  0.0_dpk
      csNph_l =  0.0_dpk
  

      interpolate_partial_wave_l:DO i = 1, nfinal_l 

         ALLOCATE(  x( nstates_l(i) ) )
         ALLOCATE(  y( nstates_l(i) ) ) 
         ALLOCATE(  yr(nstates_l(i) ) ) 
         ALLOCATE(  yi(nstates_l(i) ) ) 
         

         WRITE(*,*) "#  csNph_int::   partial wave interpolation:  L = ", L(i)
         
         x =        wph(i, 1:nstates_l(i) )     
         y =      s_si( i, 1:nstates_l(i) )
         yi =      d_r( i, 1:nstates_l(i) )
         yr =    phase( i, 1:nstates_l(i) )
         
         WRITE(*,*) "# csNph_int::    size_x = ", SIZE(x) 
         
         IF(i.EQ.lp) THEN 
            
            WRITE(*,*) "# csNph_int      lp =  ",i
            
            csnph_l(i, 1:nstates_l(i)) =  y(1:nstates_l(i)) 
            dmx_l_r(i,1:nstates_l(i)) = yi(1:nstates_l(i))  
            delta_l(i,1:nstates_l(i)) = yr(1:nstates_l(i))  
            
            !           phase(i,1:nstates_l(i)) = yr(1:nstates_l(i))  
            


!!!........  input values
            
            OPEN(nout,file='out/xy.out') 
            DO j = 1, SIZE(x)
               WRITE(nout,*) x(j),y(j),yi(j) !,yi(j) 
            ENDDO
            CLOSE(nout)

         ELSE 
            
            check = 1 ! no monitor output (1 = monitor) 
            !csN
            !phase shift
            !dmxN_r
            !dmxN_i
            CALL linear_interpolation( x, y,  SIZE(x), wn, csnph_l(i,:), nxn, 0 )
            CALL linear_interpolation( x, yr, SIZE(x), wn, delta_l(i,:), nxn, 0 )
            CALL linear_interpolation( x, yi, SIZE(x), wn, dmx_l_r(i,:), nxn, 0 )
            !       call linear_interpolation(x,yi,size(x), wn, dmx_l_i(i,:),nxn,0)
            
         ENDIF
         
         DEALLOCATE(x)
         DEALLOCATE(y) 
         DEALLOCATE(yr) 
         DEALLOCATE(yi) 
         
         
      ENDDO interpolate_partial_wave_l
      
      !
      !...  get bounds of interpolation
      !
      !

      ALLOCATE(min_xl(nfinal_l))
      ALLOCATE(max_xl(nfinal_l))

      DO  i = 1, nfinal_l         
         min_xl(i) = MINVAL( wph(i,:) )
         max_xl(i) = MAXVAL( wph(i,:) )         
      ENDDO      
      w_in  = MAXVAL( min_xl )
      w_fin = MINVAL( max_xl )
 

      DEALLOCATE(min_xl,max_xl)

!
!        cs_n = Sum_over_L( s_n_l )
!
!
      csfile  = scs//"/"//scs//TRIM(ADJUSTL(sph))//"ph."//sdat      
      OPEN(ncsfile , file=csfile)
      

      csnph = SUM(csnph_l,dim=1)
     
      save_total_cross_section:DO j = 1, nxn
         
         IF(wn(j) > w_in.AND.wn(j) < w_fin) THEN 
            !       WRITE(ncsfile,'(5(1X,1PE15.8),1X,I4)') &
            WRITE(ncsfile,'(5(1X,1PE15.8))' ) &
                 & wn(j), csnph(j), ( csnph_l(i,j), i = 1, nfinal_l )
            !,j
            !write(21,'(1X,I4,5(1X,1PE15.8))') j, wn(j),csnph(j), ( csnph_l(i,j), i = 1, nfinal_l )
            
         ENDIF
         
      ENDDO save_total_cross_section
      
      CLOSE(ncsfile)
      
      !  
      !  beta_k parameters
      !
      
      betafile  = scs//"/"//sbeta//TRIM(ADJUSTL(sph))//"ph."//sdat
      OPEN(nbetafile,file=betafile)
      
      IF(n_pad.EQ.0) THEN 
         padfile  = sdat//"/"//spad//TRIM(ADJUSTL(sph))//"ph."//TRIM(ADJUSTL(sn_pad))//"."//sdat     
         OPEN(npadfile, file=padfile)
      ENDIF

      !      padfile  = scs//"/"//spad//TRIM(ADJUSTL(sph))//"ph."//TRIM(ADJUSTL(sn_pad))//"."//sdat

      
      !for each photon energy = electron kinetic energy, ek(j)


      calculate_pads:DO j = 1, nxn          
         

         !total N-photon transition amplitude
         m_n_2 = SUM( dmx_l_r(:,j)**2 ) + SUM( dmx_l_i(:,j)**2 )  
         
         !     w_fin = 4.32
         IF(wn(j) > w_in.AND.wn(j) < w_fin) THEN 
            
            ! b_k assymetry parameter

            bp = 0.0_dpk
            DO kp = 0, nphotons
               bp(kp) = beta_k(2*kp, m_n_2, l, dmx_l_r(:,j), delta_l(:,j), nfinal_l)           
            ENDDO
            WRITE(nbetafile,'(6(1X,1PE15.8))') wn(j), (bp(i), i = nphotons,0,-1)

            
            !        WRITE(nbetafile,'(6(1X,1PE15.8))') wn(j), (bp(i)/bp(1),i = nphotons,0,-1)


        
            !        bp = bp
            !/m_n_2
            !

            IF(j==n_pad) THEN 
           
               padfile  = scs//"/"//spad//TRIM(ADJUSTL(sph))//"ph."//TRIM(ADJUSTL(sn_pad))//"."//sdat
               PRINT*, " open padfile = ", j, n_pad, wn(j), nphotons, padfile, npadfile 
               OPEN(npadfile, file=padfile)
               CALL pad(n_pad, wn(j), bp, nphotons, npadfile)
               CLOSE(npadfile)
               
            ELSE  IF(n_pad==0) THEN               
               CALL pad(n_pad,wn(j),bp, nphotons,npadfile)           
            ENDIF
!
!!!!!!!                    endif
!
         ENDIF
      ENDDO calculate_pads
      
      WRITE(*,*) "# pad::            pad stored in  ", padfile
  

      CLOSE(nbetafile)
      CLOSE(npadfile)
      
      
      !  return
    END PROGRAM CSNPH_TOTAL
!#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
