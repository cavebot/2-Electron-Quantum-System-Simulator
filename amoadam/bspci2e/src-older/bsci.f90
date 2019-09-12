!
!    bsci2e.f90 : inactive, not finished

!


PROGRAM bsci2e

  USE PRECISION,only:dpk
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
  INTEGER                                :: n_a, l_a
  CHARACTER(len=15)                      :: action 
  CHARACTER(len=6)                       :: gauge
!!
  INTERFACE
     SUBROUTINE get_command_line_arg(action, rmax, nb)
       IMPLICIT NONE
       REAL(dpk),         INTENT(inout) :: rb
       INTEGER,           INTENT(inout) :: nb
       CHARACTER(len=15), INTENT(inout) :: action 
!       CHARACTER(len=6),  INTENT(inout) :: gauge
     END SUBROUTINE get_command_line_arg
  END INTERFACE
  !EXE!


  nout  = 16  
  nb    = 100
  rb    = 100.0_dpk
  action='write'
!  gauge='no'

  CALL get_command_line_arg(action,rb,nb)
  CALL input
  
  OPEN(nout, file='out/bsci2e.log')

  WRITE(*,*) '@ main program bspci end.'

END PROGRAM bsci2e
!
!
!
SUBROUTINE get_command_line_arg(action,rb,nb)
  !
  USE PRECISION, ONLY: DPK
  !
  IMPLICIT NONE
  !ARG!
  REAL(dpk),         INTENT(inout) :: rb
  INTEGER,           INTENT(inout) :: nb
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
     WRITE(*,'(a45,1X,i3 )') '# bspci::                                rb = ', rb
     WRITE(*,'(a45,1X,i3)')  '# bspci::                                nb = ', nb
     WRITE(*,'(a45,1X,a40)') '# bspci::                            action = ', action
  ENDIF

  get_arguments: DO iarg = 1, narg

     CALL GET_COMMAND_ARGUMENT(IARG,CMD,LENGTH,ISTATUS)

     IF(MOD(iarg,2)==1) THEN    !iarg = 1, 3, 5, ...

        READ(cmd,*) what

        IF(what=='-o') THEN 
           WRITE(*,*) '# bspci::'
           WRITE(*,*) '# bspci:: Available options:'
           WRITE(*,*) '# bspci:: 0. -o information '
           WRITE(*,*) '# bspci:: 1. rb= box radius,        [100.0]'
           WRITE(*,*) '# bspci:: 2. nb= nof splines,       [100]'
           WRITE(*,*) '# bspci:: 3.  a= action           (character,l=40),  [save]'
     
      WRITE(*,*) '# bspci::        (write)   : writes inp/bsci2e.inp file'
           WRITE(*,*) '# bspci:: happy end.'  
           STOP
        ENDIF

        nxt_arg = iarg + 1

        WRITE(*,'(a45,1X,a40)') '# bspci::                             what = ', what
        IF(     (what.NE.'rb='        ).AND.     &
             &  (what.NE.'nb='        )) THEN  
           WRITE(*,'(a60)') 'available options for :'
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
        IF(what=='rb=') THEN 
           READ(cmd,'(G8)')  rb    
           WRITE(*,'(a45,1X,i3)')  '# bspci::                                rb = ', rb
        ELSE IF(what=='nb=') THEN
           READ(cmd,*)  nb 
           WRITE(*,'(a45,1X,i3 )') '# bspci::                                nb = ', nb
        ELSE IF(what=='a=') THEN
           action = cmd
           IF(  (action.NE.'write'       ) ) THEN
              WRITE(*,'(a60)') ' available options for a='
              WRITE(*,'(a60)') ' (write)'
           ENDIF
           WRITE(*,'(a45,1X,a40)') '# bspci::                           action = ', action
        ELSE
           CYCLE
        ENDIF
        
     ENDIF
  ENDDO get_arguments
  
  !............... save history
  OPEN(nascii, file='log/bspci_history.log', position='append')
  WRITE(nascii,'(a10,1x,a10,2x,a60)') date,time, line    
  CLOSE(nascii)


  RETURN
END SUBROUTINE get_command_line_arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EOF


