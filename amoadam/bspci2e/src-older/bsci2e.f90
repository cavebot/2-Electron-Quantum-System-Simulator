!
!    bsci2e.f90 : inactive, not finished

!


PROGRAM bsci2e
  !
!  USE PRECISION,only:dpk
  USE param

  !..................

  IMPLICIT NONE
  INTEGER  i,j
  INTEGER  l, ns   ! nof energy states at symmetry 'l' (nbs + ncs)
  ! I/O
  INTEGER NOUT                            
  !
  CHARACTER*100 ARGV
  CHARACTER*80 LINE, EXE, WHAT
  CHARACTER*40  CMD
  INTEGER       NARG,IARG, NXT_ARG
  INTEGER  LENGTH, ISTATUS
  INTEGER  COMMAND_ARGUMENT_COUNT
  !ARG!
  REAL(dpk)                              :: z_n
  REAL(dpk)                              :: r_b
  INTEGER                                :: n_b
  CHARACTER(len=15)                      :: action 
!!
  INTERFACE
     SUBROUTINE get_command_line_arg(action, z_n, r_b, n_b)
       !
       USE                PRECISION, ONLY: DPK
       IMPLICIT NONE
       !
       REAL(dpk),         INTENT(inout) :: z_n
       REAL(dpk),         INTENT(inout) :: r_b
       INTEGER,           INTENT(inout) :: n_b
       CHARACTER(len=15), INTENT(inout) :: action 
       !
     END SUBROUTINE get_command_line_arg
  END INTERFACE
  !EXE!


  nout  = 16  
  ! default values
  z_n    = 2.0
  n_b    = 100
  r_b    = 100.0_dpk
  action='write'


  CALL get_command_line_arg(action,z_n,r_b,n_b)
  CALL input_read_write(z_n, r_b, n_b)

  OPEN(nout, file='log/bsci2e.log')

  WRITE(*,*) '@ main program bspci2e end.'

END PROGRAM bsci2e
!
!
!
SUBROUTINE get_command_line_arg(action,z_n,r_b,n_b)
  !
  USE PRECISION, ONLY: DPK
  !
  IMPLICIT NONE
  !ARG!
  REAL(dpk),         INTENT(inout) :: z_n
  REAL(dpk),         INTENT(inout) :: r_b
  INTEGER,           INTENT(inout) :: n_b
  CHARACTER(len=15), INTENT(inout) :: action 
!  CHARACTER(len=6),  INTENT(inout) :: gauge
  !LOC!
  CHARACTER*100 ARGV
  CHARACTER*180 LINE, EXE, WHAT
  CHARACTER*40  CMD 
  INTEGER       NARG,IARG, NXT_ARG
  INTEGER       LENGTH, ISTATUS
  INTEGER       COMMAND_ARGUMENT_COUNT
  CHARACTER(len=10)       :: date,time,zone
  INTEGER, DIMENSION(dpk) :: values 
  INTEGER nascii
  !EXE!
 
  nascii = 1
  
  CALL DATE_AND_TIME(date, time, zone, values)
  
  NARG = COMMAND_ARGUMENT_COUNT()    
  CALL GET_COMMAND( LINE, LENGTH, ISTATUS)     
  CALL GET_COMMAND_ARGUMENT(0,EXE,LENGTH,ISTATUS) 


  WRITE(*,*)'#'
  WRITE(*,'(a45,1X,a40)')'# w1e::            executable command exe = ', line
  WRITE(*,'(a45,1X,i2)') '# wf1e::        nof arguments        narg = ', narg  
  IF(narg.LE.1) THEN
     WRITE(*,'(a45,1X,G8.3)') '# bsci2e::                                zn = ', z_n
     WRITE(*,'(a45,1X,G8.3)') '# bsci2e::                                rb = ', r_b
     WRITE(*,'(a45,1X,i3)')   '# bsci2e::                                nb = ', n_b
     WRITE(*,'(a45,1X,a40)')  '# bsci2e::                            action = ', action
  ENDIF

  get_arguments: DO iarg = 1, narg
     !
     CALL GET_COMMAND_ARGUMENT(IARG,CMD,LENGTH,ISTATUS)
     !
     IF(MOD(iarg,2)==1) THEN    !iarg = 1, 3, 5, ...

        READ(cmd,*) what

        IF(what=='-o') THEN 
           WRITE(*,*) '# bsci2e::'
           WRITE(*,*) '# bsci2e:: Available options:'
           WRITE(*,*) '# bsci2e:: 0. -o information '
           WRITE(*,*) '# bsci2e:: 1. zn= nuclei charge,       [2.0]'
           WRITE(*,*) '# bsci2e:: 2. rb= box radius,        [100.0]'
           WRITE(*,*) '# bsci2e:: 3. nb= nof splines,         [100]'
           WRITE(*,*) '# bsci2e:: 4.  a= action           (character,l=40),  [save]'
     
      WRITE(*,*) '# bsci2e::        (write)   : writes inp/h1e.inp file'
           WRITE(*,*) '# bsci2e:: happy end.'  
           STOP
        ENDIF

        nxt_arg = iarg + 1

        WRITE(*,'(a45,1X,a40)') '# bsci2e::                             what = ', what
        IF(     (what.NE.'zn='        ).AND.     &
                (what.NE.'rb='        ).AND.     &
             &  (what.NE.'nb='        )) THEN  
           WRITE(*,'(a60)') 'available options for :'
           WRITE(*,'(a60)') ' (zn= nuclei charge)'
           WRITE(*,'(a60)') ' (rb= box radius)'
           WRITE(*,'(a60)') ' (nb= nof of b-splines)'
           WRITE(*,'(a60)') ' (a= (write) )'
           STOP
        ENDIF
        CYCLE
     ENDIF

!....................................................................................!
!                                                                                    !

     IF(iarg == nxt_arg ) THEN
        IF(what=='zn=') THEN 
           READ(cmd,'(e8.3)')  z_n    
           WRITE(*,'(a45,1X,G10.5)')  '# bsci2e::                                zn = ', z_n
        ELSE IF(what=='rb=') THEN 
           READ(cmd,'(e8.3)')  r_b    
           WRITE(*,'(a45,1X,G10.5)')  '# bsci2e::                                rb = ', r_b
        ELSE IF(what=='nb=') THEN
           READ(cmd,*)  n_b 
           WRITE(*,'(a45,1X,i3 )')     '# bsci2e::                                nb = ', n_b
        ELSE IF(what=='a=') THEN
           action = cmd
           IF(  (action.NE.'write'       ) ) THEN
              WRITE(*,'(a60)') ' available options for a='
              WRITE(*,'(a60)') ' (write)'
           ENDIF
           WRITE(*,'(a45,1X,a40)') '# bsci2e::                           action = ', action
        ELSE
           CYCLE
        ENDIF
        
     ENDIF
  ENDDO get_arguments
  !
  !............... save history
  !
  OPEN(nascii, file='log/bspci_history.log', position='append')
  WRITE(nascii,'(a10,1x,a10,2x,a60)') date,time, line    
  CLOSE(nascii)
  !
  RETURN
  !
END SUBROUTINE get_command_line_arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EOF


