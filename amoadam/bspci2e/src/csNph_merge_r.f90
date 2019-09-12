!
!25082004LAANKYOTO:
!
!
! 
!   csNph_read : takes the output partial cs values of csNph.f for different 
!                boxes joins (sort?) them and saves one file for each 
!                partial wave. 
!
!   input   txt files: inp/csNph_read.inp
!   input   dat files: dat/csNph-L.rXXX.dat  
!                                        L = 0,2,4,...  nphotons = even
!                                     or L = 1, 3,...   nphotons = odd
!   output  txt files: out/csNph_read.out
!   output  dat files: dat/csNph-L.dat   
!                                        L = 0,2,4,...  nphotons = even
!                                     or L = 1, 3,...   nphotons = odd
!
!   notes: 
!             (1) input for partial-waves should go to a module
!             (2) final values are stored unsorted. Should add 
!                 a sort routine. An initial try with dpsort did 
!                 succeed   
!
!29082004:  
!
!

PROGRAM CSNPH_merge_r
  use precision, only:dpk
  IMPLICIT NONE
  INTEGER,          PARAMETER :: nstates = 8000
  CHARACTER(LEN=*), PARAMETER :: scsnph = "csNph-l"
  INTEGER NINP,NOUT ,NBIN
  INTEGER NPHOTONS, NBOXES
  INTEGER  NFINAL_L, L_INITIAL
  INTEGER I,J,IE,IT
  INTEGER NTOTAL_STATES, NDIM 
  REAL(DPK) W1,W2,W3,W4,W5,W6,W7,W8,SCALE
  !
  CHARACTER( LEN = 6  ) SL
  CHARACTER( LEN = 6  ), ALLOCATABLE, DIMENSION(:)   :: SR_BOXES
  CHARACTER( LEN = 30 ), ALLOCATABLE, DIMENSION(:)   :: CS_L_FILE
  CHARACTER( LEN = 30 ), ALLOCATABLE, DIMENSION(:,:) :: CS_LR_FILE
  INTEGER, ALLOCATABLE,     DIMENSION(:) :: L,SGN_R
  INTEGER, ALLOCATABLE,     DIMENSION(:) :: IPERM      !DPSORT
  INTEGER IFAIL, IP                                    !DPSORT
  REAL(DPK), ALLOCATABLE, DIMENSION(:) :: WPH, S_SI, S_AU
  REAL(DPK), ALLOCATABLE, DIMENSION(:) :: D_2, D_R,  D_I 
  REAL(DPK), ALLOCATABLE, DIMENSION(:) :: PHASE,AN
!
  CHARACTER(LEN = 100) ARGV

!
  
  WRITE(*,*) "#  READ_CS:: EXECUTING AT "
  
!....

  NBIN = 1
  NINP = 5
  NOUT = 16
  
  !..   read input data

  CALL GETARG(1, ARGV)
  READ(ARGV,*) NPHOTONS

        
  OPEN(UNIT=NINP,FILE="inp/csNph_merge_r.inp")
  READ(NINP,*) NBOXES
  ALLOCATE( SR_BOXES(NBOXES) ) 
  ALLOCATE( SGN_R(NBOXES) ) 
  
  READ(NINP,*) (SR_BOXES(I),I=1,NBOXES)
  READ(NINP,*) (SGN_R(I),I=1,NBOXES)
  
  CLOSE(NINP)

! make a check of input sign for boxes

  DO I = 1, NBOXES
     IF( ABS(SGN_R(I)).NE.1) THEN 
        WRITE(*,*)"# csNph_r::   error in input file "
        WRITE(*,*)"# csNph_r::           i,sgn(i)) = ", i,sgn_r(i)
        stop
     ENDIF
     
!     write(*,*) i,sr_boxes(i),sgn_r(i)
  ENDDO

!xxx scale the 4-photon cs  below 10-{-100} (for xmgr reasons)

  SCALE = 1.0_dpk
  IF(NPHOTONS.EQ.4) THEN 
     SCALE = 1.0e+113_dpk
     WRITE(*,*) "# csNph_merge_r:: cross sections have re-scaled by a factor equal to ", scale
  ENDIF

!xxx



  ! determine number of final angular symmetries

  NFINAL_L = 0
  DO I = NPHOTONS, 0, -2     
     NFINAL_L = NFINAL_L + 1 
  ENDDO

  L_INITIAL = I + 2

  ALLOCATE(         L(NFINAL_L) )
  ALLOCATE( CS_L_FILE(NFINAL_L) ) 
  ALLOCATE( CS_LR_FILE(NFINAL_L,NBOXES) ) 

  ! initialize final angular symmetry

  DO I = 1, NFINAL_L
     L(I) = L_INITIAL + 2*(I-1)
  ENDDO

  WRITE(*,*) "#  READ_CS::         INPUT DATA:              "
  WRITE(*,*) "#"
  WRITE(*,*) "#  READ_CS::                       NPHOTONS = ", NPHOTONS
  WRITE(*,*) "#  READ_CS::                         NBOXES = ", NBOXES
  WRITE(*,*) "#  READ_CS::                   RADIUS BOX ARE "  
  WRITE(*,*) "#  R = ", SR_BOXES
  WRITE(*,*) "#"
  WRITE(*,*) "#"
  WRITE(*,*) "#  READ_CS::     CALCULATED DATA:             "
  WRITE(*,*) "#"
  WRITE(*,*) "#  READ_CS:: NOF FINAL SYMMETRIES  NFINAL_L = ", NFINAL_L
  WRITE(*,*) "#  READ_CS:: INITIAL SYMMETRY     L_INITIAL = ", L_INITIAL
  WRITE(*,*) "#  READ_CS::               FINAL SYMMETRIES   "
  WRITE(*,*) "#  READ_CS::                              L = ", L
  WRITE(*,*) "#"


!.................................

  DO I = 1, NFINAL_L

     WRITE(SL,'(I6)') L(I)   !convert number to string

     DO J = 1, NBOXES        
        CS_LR_FILE(I,J) = "dat/"//SCSNPH//TRIM(ADJUSTL(SL))//".r"//TRIM(SR_BOXES(J))//".dat"
     ENDDO

     CS_L_FILE(I) = "dat/"//SCSNPH//TRIM(ADJUSTL(SL))//".dat"
  ENDDO

  WRITE(*,*) "# csNph_r::    CS DATA FILES:             "
  WRITE(*,*) "#"

  DO J = 1, NBOXES
     WRITE(*,'(1x,a30)') ( CS_LR_FILE(I,J), I =1,NFINAL_L)
  ENDDO
  
  WRITE(*,*) (CS_L_FILE(I), I = 1, NFINAL_L)


  NDIM = NBOXES * NSTATES

  ALLOCATE(  WPH( ndim ) )
  ALLOCATE( S_SI( ndim ) )
  ALLOCATE( S_AU( ndim ) )
  ALLOCATE(  D_2( ndim) )
  ALLOCATE(  D_R( ndim) )
  ALLOCATE(  D_I( ndim) )
  ALLOCATE(  PHASE( ndim) )
  ALLOCATE(  AN( ndim) )



  DO I = 1, NFINAL_L

     WRITE(*,*) "# csNph_r::       partial wave l = ", L(I) 

     IT = 1
     DO J = 1, NBOXES
       
        WRITE(*,*) "# csNph_r::                   R = ", SR_BOXES(J)  
        WRITE(*,*) "# csNph_r::            data file :", CS_LR_FILE(I,J) 

        OPEN(UNIT=J+10, FILE=CS_LR_FILE(I,J))

        IT = IT - 1
        DO IE = 1, NSTATES
           
           READ (J+10, FMT='(8(1X,1PE15.8))', END = 201) &
                &              W1,W2,W3,W4,W5,W6,W7,W8 
           
           IT = IT + 1


           WPH(IT)  = W1
           S_SI(IT) = W2
           S_AU(IT) = W3
            D_2(IT) = W4 
            D_R(IT) = W5 * SGN_R(J)
            D_I(IT) = W6 * SGN_R(J)
          PHASE(IT) = W7
             AN(IT) = W8
        ENDDO

        STOP

201     CLOSE(J+10)

        WRITE(*,*) "# csNph_r::           nof states_l = ", IE 
        
     ENDDO

     NTOTAL_STATES  = IT

     WRITE(*,*) "#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
     WRITE(*,*) "#"
     WRITE(*,*) "# csNph_merge_r::          nof total states  = ", NTOTAL_STATES 
     WRITE(*,*) "#"

     ALLOCATE(IPERM(NTOTAL_STATES))     

     ifail = 0 

!slatec
     CALL dpsort(WPH, NTOTAL_STATES, IPERM, -1,IFAIL)

     if(ifail.ne.0) then 

        write(*,*) "#  error in dpsort routine "
        write(*,*) "#                  ifail = ",ifail
       
     endif


     OPEN(UNIT=I+10, FILE = CS_L_FILE(I))

     WRITE(I+10,*)'& wph (a.u), csN_if, dmxN_if_real, dmxN_if_im, A_n, phase'
     WRITE(I+10,*)'&', NTOTAL_STATES

     DO IT =  NTOTAL_STATES , 1, -1

        IP = IPERM(IT)

        S_SI(IP) = SCALE*S_SI(IP)            

        WRITE(I+10,FMT='(6(1X,1PE15.8))')  &
             &   WPH(IP)                   &
             &  ,S_SI(IP)                  &
             &  ,D_R(IP)                   &
             &  ,D_I(IP)                   &
             &  ,AN(IP)                    &
             &  ,PHASE(IP)
     ENDDO
     
     CLOSE(I+10)

     DEALLOCATE(IPERM)          
  ENDDO
  !
END PROGRAM CSNPH_MERGE_R
!#EOF
