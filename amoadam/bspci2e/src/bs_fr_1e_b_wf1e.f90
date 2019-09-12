!###################################################################
! laan jan/feb 2003 iesl days : bsp1ef
!
!
!   Author         :  Lambros Nikolopoulos
!   I/0 files      : inp/h1e.inp inp/hwf1e.inp out/hwf1e.out
!   I/O data files :   DAT/HD1E-L.LIN.DAT,  L = 0, 1, ..., NW
!                  :   DAT/HD1E-L.INV.DAT,  L = 0, 1, ..., NW
!                  :   DAT/HWF1E-L.INV.DAT, L = 0, 1, ..., NW
!                  :   DAT/HWF1E-L.INV.DAT, L = 0, 1, ..., NW
!                  :   DAT/HWF1E-L.INV.DAT, L = 0, 1, ..., NW
!                  :   dat/hwf1e.inv.dat, dat/hwf1e.lin.dat
!
!    JAN/FEB 2003 :
!
!    hwf1e.f90 :
!               1. reads energies + coefficients for l = 0, 1, 2,..,nw
!                  calulated from h1e.f90 program and stores the 
!                  wavefunctions P_nl(r) in interval [0,R].
!               2. Plots the radial wavefunction P(n, l) and 
!                  stores it in dat/hwf1e.lin.dat.
!                  nplot, lplot are read from 'inp/hwf1e.inp' input file.
!
!
!   I/0 files      : inp/h1e.inp inp/wf1e.inp out/hwf1e.out
!   I/O data files :   HD1E-L.LIN.DAT, L = 0, 1, ...NW
!                  :   HD1E-L.INV.DAT, L = 0, 1, ...NW
!                  :   dat/hwf1e.inv.dat, dat/hwf1e.lin.dat
!
!  2003.02.17      :routine reaf_target_wf added for code reliability      
!  2006.02.10      :routine reaf_target_wf removed
!
!

PROGRAM hwf1e

  USE param
  USE DATA, ONLY: read_v, read_v_mx, write_wf1e, write_fxd_wf1e, write_target_energies
  
  implicit none

  INTEGER I, J 
  INTEGER L, NS, NVALENCE
  INTEGER LMN,LMX
  INTEGER NOUT
  REAL(DPK), DIMENSION(:),   POINTER   :: EN
  REAL(DPK), DIMENSION(:,:), POINTER   :: CE
  
  CHARACTER*100 ARGV
!.............................. get arguments!

  CALL GETARG(1, ARGV)
  READ(ARGV,*)   L


! to be modified 

  NVALENCE    = 1

  IF(L.LT.0) THEN

     LMN = 0
     LMX = -L

     CALL WRITE_TARGET_ENERGIES(LMN,LMX,"dat/en1e.dat")
     STOP
  ENDIF


!note a second argument is taken from the write_wf1e sbroutine

!...............

  WRITE(*, '(a2,1X,a50,i3)') "#"," hwf1e.f90:   Rwf1e_f  -l = ", L
  WRITE(*, '(a2,1X,a50,i3)') "#"," Input file:   inp/hwf1e.inp"


  call input

!..............

  NOUT       = 16
  OPEN(NOUT, FILE='out/wf1e.log')

!  ns = nbs + ncs
     
!...............


        
  IF(code=='bs_fx')     THEN 

     WRITE(*,*) '# main:: open d1e-files for the bs_fx code'

     CALL READ_V_MX( L, EN, CE, "d1e-")     ! get e_{nl} and ce_{nl}
     CALL WRITE_FXD_WF1E(L, EN, CE)
                                            ! store radial wfs p_nl(r),p'_nl(r)
                                            ! for use with bs_fx_1e/bs_fx_2e
  ELSE
     WRITE(*,*) '# main:: open h1e-files for the bs_fr code'

     CALL READ_V_MX( L, EN, CE, "h1e-")     ! get e_{nl} and ce_{nl}
     CALL WRITE_WF1E(L, EN, CE)             ! store radial wfs p_nl(r),p'_nl(r)

  ENDIF


  DEALLOCATE( EN )
  DEALLOCATE( CE )

!..............

END PROGRAM HWF1E
!#################################################################
!EOF!
             
