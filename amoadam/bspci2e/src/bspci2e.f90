!
!    bsci2e.f90 : inactive, not finished

!


PROGRAM bsci2e
  !
  USE PRECISION,only:dpk
  USE param
  USE M_kracken

  !..................
  IMPLICIT NONE
  
  CHARACTER(len=255)  ::  argument
  INTEGER             ::  len
  INTEGER             ::  ier
  REAL                ::  arg_value
  
  !
  REAL(dpk)           :: v_0
  REAL(dpk)           :: z_n
  REAL(dpk)           :: r_b
  REAL(dpk)           :: n_b
  !  INTEGER         :: n_b
  
  ! I/O
  INTEGER NOUT 
  !EXE!


  !define command line options  
  CALL kracken(   'Rbspci2e',' -v 4.96 -zn 2.0 -rb 100 -nb 100 ')
  !dump out command
  CALL  retrev('Rbspci2e_oo', argument,len,ier)
  WRITE(*,*)'#',argument(:len)
  !atomic charge
  CALL retrev('Rbspci2e_zn',argument,len,ier)
  CALL string_to_real(argument,arg_value,ier)
  z_n = arg_value
  !
  CALL retrev('Rbspci2e_v',argument,len,ier)
  CALL string_to_real(argument,arg_value,ier)
  v_0 = arg_value
  !box radius
  CALL retrev('Rbspci2e_rb',argument,len,ier)
  CALL string_to_real(argument,arg_value,ier)
  r_b = arg_value
  !bsplines size
  CALL retrev('Rbspci2e_nb',argument,len,ier)
  CALL string_to_real(argument,arg_value,ier)
  n_b = INT(arg_value)

  WRITE(*,*)'# bspci2e:      value for  -v =',v_0
  WRITE(*,*)'# bspci2e:      value for -zn =',z_n
  WRITE(*,*)'# bspci2e:      value for -rb =',r_b
  WRITE(*,*)'# bspci2e:      value for -nb =',n_b
  
  
  CALL input_read_write(v_0, z_n, r_b, INT(n_b))
  
  nout  = 16  
  OPEN(nout, file='log/bsci2e.log')
  WRITE(*,*) '@ main program bspci2e end.'

END PROGRAM bsci2e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EOF


