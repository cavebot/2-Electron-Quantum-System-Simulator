!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
MODULE  ioroutines
  IMPLICIT NONE
  PUBLIC
CONTAINS
 
!!!###############################################
   SUBROUTINE D1EFILE(N, L)
    IMPLICIT NONE
    INTEGER N,L
    
!    WRITE(*,*) ' OPENING FILE DATA  FOR L = ', L
          
    IF(L.EQ.0) OPEN(N,FILE='dat/d1e-0.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.1) OPEN(N,FILE='dat/d1e-1.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.2) OPEN(N,FILE='dat/d1e-2.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.3) OPEN(N,FILE='dat/d1e-3.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.4) OPEN(N,FILE='dat/d1e-4.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.5) OPEN(N,FILE='dat/d1e-5.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.6) OPEN(N,FILE='dat/d1e-6.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.7) OPEN(N,FILE='dat/d1e-7.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.8) OPEN(N,FILE='dat/d1e-8.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.9) OPEN(N,FILE='dat/d1e-9.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.10) OPEN(N,FILE='dat/d1e-10.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.11) OPEN(N,FILE='dat/d1e-11.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    RETURN

  END SUBROUTINE D1EFILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE HD1EFILE(N, L, SPECTRUM_TYPE)

    IMPLICIT NONE
    INTEGER N,L
    CHARACTER(LEN=2, KIND=1) SPECTRUM_TYPE
    
!    WRITE(*,*) '#             OPENING FILE DATA  FOR L = ', L

    IF(SPECTRUM_TYPE.EQ."FR") THEN
          
       if(l.eq.0)  open(n, file='dat/hd1e-0.dat',  form='unformatted',access='sequential')
       if(l.eq.1)  open(n, file='dat/hd1e-1.dat',  form='unformatted',access='sequential')
       if(l.eq.2)  open(n, file='dat/hd1e-2.dat',  form='unformatted',access='sequential')
       if(l.eq.3)  open(n, file='dat/hd1e-3.dat',  form='unformatted',access='sequential')
       if(l.eq.4)  open(n, file='dat/hd1e-4.dat',  form='unformatted',access='sequential')
       if(l.eq.5)  open(n, file='dat/hd1e-5.dat',  form='unformatted',access='sequential')
       if(l.eq.6)  open(n, file='dat/hd1e-6.dat',  form='unformatted',access='sequential')
       if(l.eq.7)  open(n, file='dat/hd1e-7.dat',  form='unformatted',access='sequential')
       if(l.eq.8)  open(n, file='dat/hd1e-8.dat',  form='unformatted',access='sequential')
       if(l.eq.9)  open(n, file='dat/hd1e-9.dat',  form='unformatted',access='sequential')
       if(l.eq.10) open(n, file='dat/hd1e-10.dat', form='unformatted',access='sequential')
       if(l.eq.11) open(n, file='dat/hd1e-11.dat', form='unformatted',access='sequential')
      
    else if(spectrum_type.eq."FX") then

       if(l.eq.0)  open(n, file='dat/d1e-0.dat',  form='unformatted',access='sequential')
       if(l.eq.1)  open(n, file='dat/d1e-1.dat',  form='unformatted',access='sequential')
       if(l.eq.2)  open(n, file='dat/d1e-2.dat',  form='unformatted',access='sequential')
       if(l.eq.3)  open(n, file='dat/d1e-3.dat',  form='unformatted',access='sequential')
       if(l.eq.4)  open(n, file='dat/d1e-4.dat',  form='unformatted',access='sequential')
       if(l.eq.5)  open(n, file='dat/d1e-5.dat',  form='unformatted',access='sequential')
       if(l.eq.6)  open(n, file='dat/d1e-6.dat',  form='unformatted',access='sequential')
       if(l.eq.7)  open(n, file='dat/d1e-7.dat',  form='unformatted',access='sequential')
       if(l.eq.8)  open(n, file='dat/d1e-8.dat',  form='unformatted',access='sequential')
       if(l.eq.9)  open(n, file='dat/d1e-9.dat',  form='unformatted',access='sequential')
       if(l.eq.10) open(n, file='dat/d1e-10.dat', form='unformatted',access='sequential')
       if(l.eq.11) open(n, file='dat/d1e-11.dat', form='unformatted',access='sequential')

    else if(spectrum_type.eq."MX") then

       if(l.eq.0)  open(n, file ='dat/md1e-0.dat',  form='unformatted',access='sequential')
       if(l.eq.1)  open(n, file ='dat/md1e-1.dat',  form='unformatted',access='sequential')
       if(l.eq.2)  open(n, file ='dat/md1e-2.dat',  form='unformatted',access='sequential')
       if(l.eq.3)  open(n, file ='dat/md1e-3.dat',  form='unformatted',access='sequential')
       if(l.eq.4)  open(n, file ='dat/md1e-4.dat',  form='unformatted',access='sequential')
       if(l.eq.5)  open(n, file ='dat/md1e-5.dat',  form='unformatted',access='sequential')
       if(l.eq.6)  open(n, file ='dat/md1e-6.dat',  form='unformatted',access='sequential')
       if(l.eq.7)  open(n, file ='dat/md1e-7.dat',  form='unformatted',access='sequential')
       if(l.eq.8)  open(n, file ='dat/md1e-8.dat',  form='unformatted',access='sequential')
       if(l.eq.9)  open(n, file ='dat/md1e-9.dat',  form='unformatted',access='sequential')
       if(l.eq.10) open(n, file ='dat/md1e-10.dat', form='unformatted',access='sequential')
       if(l.eq.11) open(n, file ='dat/md1e-11.dat', form='unformatted',access='sequential')

    ENDIF


    RETURN

  END SUBROUTINE HD1EFILE
!!!###############################################
!   SUBROUTINE HWF1EFILE(N, L, SPECTRUM_TYPE)
   SUBROUTINE HWF1EFILE(N, L)
    IMPLICIT NONE
!    CHARACTER(LEN=2, KIND=1) SPECTRUM_TYPE
    INTEGER N,L
    
!    WRITE(*,*) '#   HWF1E:    OPENING FILE DATA  FOR L = ', L
          
!    IF(SPECTRUM_TYPE.EQ."FR") THEN

       if(l.eq.0) open(n,file='dat/hwf1e-0.dat')
       if(l.eq.1) open(n,file='dat/hwf1e-1.dat')
       if(l.eq.2) open(n,file='dat/hwf1e-2.dat')
       if(l.eq.3) open(n,file='dat/hwf1e-3.dat')
       if(l.eq.4) open(n,file='dat/hwf1e-4.dat')
       if(l.eq.5) open(n,file='dat/hwf1e-5.dat')
       if(l.eq.6) open(n,file='dat/hwf1e-6.dat')
       if(l.eq.7) open(n,file='dat/hwf1e-7.dat')
       if(l.eq.8) open(n,file='dat/hwf1e-8.dat')
       if(l.eq.9) open(n,file='dat/hwf1e-9.dat')
       if(l.eq.10) open(n,file='dat/hwf1e-10.dat')
       if(l.eq.11) open(n,file='dat/hwf1e-11.dat')
       if(l.eq.12) open(n,file='dat/hwf1e-12.dat')
       if(l.eq.13) open(n,file='dat/hwf1e-13.dat')
       if(l.eq.14) open(n,file='dat/hwf1e-14.dat')
       if(l.eq.15) open(n,file='dat/hwf1e-15.dat')
       if(l.eq.20) open(n,file='dat/hwf1e-20.dat')
       if(l.eq.25) open(n,file='dat/hwf1e-25.dat')


    RETURN

  END SUBROUTINE HWF1EFILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE HNORMFILE(N, L )
    IMPLICIT NONE
    INTEGER N, L

    
!    IF(L.GT.14) THEN
       
!       WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', &
!            &       ' MODIFY SUBROUTINE NORMFILE. EXITING...'
!       STOP
!    ELSE

     WRITE(*,*) '#NORM1E ::                    DATA FOR L = ', L

!    ENDIF
         
       IF(L.EQ.0)  OPEN(N, FILE = 'dat/hnorm-0.dat',  FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.1)  OPEN(N, FILE = 'dat/hnorm-1.dat',  FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.2)  OPEN(N, FILE = 'dat/hnorm-2.dat',  FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.3)  OPEN(N, FILE = 'dat/hnorm-3.dat',  FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.4)  OPEN(N, FILE = 'dat/hnorm-4.dat',  FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.5)  OPEN(N, FILE = 'dat/hnorm-5.dat',  FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.6)  OPEN(N, FILE = 'dat/hnorm-6.dat',  FORM = 'UNFORMATTED', ACCESS = 'SEQUENTIAL')
       IF(L.EQ.7)  OPEN(N, FILE = 'dat/hnorm-7.dat',  FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.8)  OPEN(N, FILE = 'dat/hnorm-8.dat',  FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.9)  OPEN(N, FILE = 'dat/hnorm-9.dat',  FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.10) OPEN(N, FILE = 'dat/hnorm-10.dat', FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.11) OPEN(N, FILE = 'dat/hnorm-11.dat', FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.12) OPEN(N, FILE = 'dat/hnorm-12.dat', FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.13) OPEN(N, FILE = 'dat/hnorm-13.dat', FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.14) OPEN(N, FILE = 'dat/hnorm-14.dat', FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.15) OPEN(N, FILE = 'dat/hnorm-15.dat', FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.20) OPEN(N, FILE = 'dat/hnorm-20.dat', FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
       IF(L.EQ.25) OPEN(N, FILE = 'dat/hnorm-25.dat', FORM ='UNFORMATTED',  ACCESS = 'SEQUENTIAL')
    
       

    RETURN

    END SUBROUTINE HNORMFILE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE HNORMASCIIFILE(N, L)

      IMPLICIT NONE
      INTEGER N, L

!      IF(L.GT.14) THEN
         
!         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', &
!              &       ' MODIFY SUBROUTINE NORMFILE. EXITING...'
!         STOP
         
!      ELSE

         WRITE(*,*) '#NORM1E ::                     OUT FOR L = ', L

!      ENDIF

      IF(L.EQ.0)  OPEN( N, FILE = 'out/hnorm-0.out'  )
      IF(L.EQ.1)  OPEN( N, FILE = 'out/hnorm-1.out'  )
      IF(L.EQ.2)  OPEN( N, FILE = 'out/hnorm-2.out'  )
      IF(L.EQ.3)  OPEN( N, FILE = 'out/hnorm-3.out'  )
      IF(L.EQ.4)  OPEN( N, FILE = 'out/hnorm-4.out'  )
      IF(L.EQ.5)  OPEN( N, FILE = 'out/hnorm-5.out'  )
      IF(L.EQ.6)  OPEN( N, FILE = 'out/hnorm-6.out'  )
      IF(L.EQ.7)  OPEN( N, FILE = 'out/hnorm-7.out'  )
      IF(L.EQ.8)  OPEN( N, FILE = 'out/hnorm-8.out'  )
      IF(L.EQ.9)  OPEN( N, FILE = 'out/hnorm-9.out'  )
      IF(L.EQ.10) OPEN( N, FILE = 'out/hnorm-10.out' )
      IF(L.EQ.11) OPEN( N, FILE = 'out/hnorm-11.out' )
      IF(L.EQ.12) OPEN( N, FILE = 'out/hnorm-12.out' )
      IF(L.EQ.13) OPEN( N, FILE = 'out/hnorm-13.out' )
      IF(L.EQ.14) OPEN( N, FILE = 'out/hnorm-14.out' )
      IF(L.EQ.15) OPEN( N, FILE = 'out/hnorm-15.out' )
      IF(L.EQ.20) OPEN( N, FILE = 'out/hnorm-20.out' )
      IF(L.EQ.25) OPEN( N, FILE = 'out/hnorm-25.out' )

      RETURN

    END SUBROUTINE HNORMASCIIFILE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE HPHASEASCIIFILE(N, L)

      IMPLICIT NONE
      INTEGER N, L

!      IF(L.GT.14) THEN
         
!         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', &
!              &       ' MODIFY SUBROUTINE NORMFILE. EXITING...'
!         STOP
!      ELSE

         WRITE(*,*) '#PHASE_1E ::                   OUT FOR L = ', L
!      ENDIF

      IF(L.EQ.0)  OPEN( N, FILE = 'out/hphase-0.out'  )
      IF(L.EQ.1)  OPEN( N, FILE = 'out/hphase-1.out'  )
      IF(L.EQ.2)  OPEN( N, FILE = 'out/hphase-2.out'  )
      IF(L.EQ.3)  OPEN( N, FILE = 'out/hphase-3.out'  )
      IF(L.EQ.4)  OPEN( N, FILE = 'out/hphase-4.out'  )
      IF(L.EQ.5)  OPEN( N ,FILE = 'out/hphase-5.out'  )
      IF(L.EQ.6)  OPEN( N, FILE = 'out/hphase-6.out'  )
      IF(L.EQ.7)  OPEN( N, FILE = 'out/hphase-7.out'  )
      IF(L.EQ.8)  OPEN( N, FILE = 'out/hphase-8.out'  )
      IF(L.EQ.9)  OPEN( N, FILE = 'out/hphase-9.out'  )
      IF(L.EQ.10) OPEN( N, FILE = 'out/hphase-10.out' )
      IF(L.EQ.11) OPEN( N, FILE = 'out/hphase-11.out' )
      IF(L.EQ.12) OPEN( N, FILE = 'out/hphase-12.out' )
      IF(L.EQ.13) OPEN( N, FILE = 'out/hphase-13.out' )
      IF(L.EQ.14) OPEN( N, FILE = 'out/hphase-14.out' )
      IF(L.EQ.15) OPEN( N, FILE = 'out/hphase-15.out' )
      IF(L.EQ.20) OPEN( N, FILE = 'out/hphase-20.out' )
      IF(L.EQ.25) OPEN( N, FILE = 'out/hphase-25.out' )

      RETURN

    END SUBROUTINE HPHASEASCIIFILE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  SUBROUTINE COREFILE(N, L)
    IMPLICIT NONE
    INTEGER N,L

!    WRITE(*,*) '#             OPENING FILE OUTA  FOR L = ', L    

    IF(L.EQ.0) OPEN(N,FILE='dat/core-0.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.1) OPEN(N,FILE='dat/core-1.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.2) OPEN(N,FILE='dat/core-2.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    IF(L.EQ.3) OPEN(N,FILE='dat/core-3.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    
    RETURN
  END SUBROUTINE COREFILE
!!!#####################################################################
  subroutine hmxfile(n, l)
    implicit none
    integer n,l
    

!    write(*,*) '#             opening file data  for l = ', l    

    if(l.eq.0) open(n,file='dat/hmx-0.dat', form='unformatted',access='sequential')
    if(l.eq.1) open(n,file='dat/hmx-1.dat', form='unformatted',access='sequential')
    if(l.eq.2) open(n,file='dat/hmx-2.dat', form='unformatted',access='sequential')
    if(l.eq.3) open(n,file='dat/hmx-3.dat', form='unformatted',access='sequential')
    if(l.eq.4) open(n,file='dat/hmx-4.dat', form='unformatted',access='sequential')
    if(l.eq.5) open(n,file='dat/hmx-5.dat', form='unformatted',access='sequential')
    if(l.eq.6) open(n,file='dat/hmx-6.dat', form='unformatted',access='sequential')
    if(l.eq.7) open(n,file='dat/hmx-7.dat', form='unformatted',access='sequential')
    if(l.eq.8) open(n,file='dat/hmx-8.dat', form='unformatted',access='sequential')
    if(l.eq.9) open(n,file='dat/hmx-9.dat', form='unformatted',access='sequential')
    if(l.eq.10) open(n,file='dat/hmx-10.dat', form='unformatted',access='sequential')
    if(l.eq.11) open(n,file='dat/hmx-11.dat', form='unformatted',access='sequential')
    return
  end subroutine hmxfile
!!!#####################################################################
  subroutine bmxfile(n, l)
    implicit none
    integer n,l

!    write(*,*) '#             opening file data  for l = ', l        

    
    if(l.eq.0) open(n,file='dat/bmx-0.dat', form='unformatted',access='sequential')
    if(l.eq.1) open(n,file='dat/bmx-1.dat', form='unformatted',access='sequential')
    if(l.eq.2) open(n,file='dat/bmx-2.dat', form='unformatted',access='sequential')
    if(l.eq.3) open(n,file='dat/bmx-3.dat', form='unformatted',access='sequential')
    if(l.eq.4) open(n,file='dat/bmx-4.dat', form='unformatted',access='sequential')
    if(l.eq.5) open(n,file='dat/bmx-5.dat', form='unformatted',access='sequential')
    if(l.eq.6) open(n,file='dat/bmx-6.dat', form='unformatted',access='sequential')
    if(l.eq.7) open(n,file='dat/bmx-7.dat', form='unformatted',access='sequential')
    if(l.eq.8) open(n,file='dat/bmx-8.dat', form='unformatted',access='sequential')
    if(l.eq.9) open(n,file='dat/bmx-9.dat', form='unformatted',access='sequential')
    if(l.eq.10) open(n,file='dat/bmx-10.dat', form='unformatted',access='sequential')
    if(l.eq.11) open(n,file='dat/bmx-11.dat', form='unformatted',access='sequential')
    return
  end subroutine bmxfile
!!!#################################################################################
  SUBROUTINE HR12FILE(N, L)
    IMPLICIT NONE
    INTEGER N,L
    
    IF(L.GT.14) THEN
       
       WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14 MODIFY SUBROUTINE HR12FILE. EXITING...'
         STOP
      ELSE

         WRITE(*,*) '#             R12: OPENING FILE DATA  FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(N,FILE='DAT/HR12-0.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='DAT/HR12-1.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='DAT/HR12-2.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='DAT/HR12-3.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='DAT/HR12-4.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='DAT/HR12-5.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='DAT/HR12-6.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='DAT/HR12-7.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='DAT/HR12-8.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='DAT/HR12-9.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='DAT/HR12-10.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='DAT/HR12-11.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='DAT/HR12-12.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='DAT/HR12-13.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='DAT/HR12-14.DAT',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      RETURN
    
    END SUBROUTINE HR12FILE

!!!#######################################################################
  SUBROUTINE HCFGFILE(N, L)
    IMPLICIT NONE
    INTEGER N, L

!!!........................

    WRITE(*,*) ' OPENING FILE DATA  FOR L = ', L
    

    IF(L.EQ.0) OPEN(N,FILE='INP/CFG-0.INP')
    IF(L.EQ.1) OPEN(N,FILE='INP/CFG-1.INP')
    IF(L.EQ.2) OPEN(N,FILE='INP/CFG-2.INP')
    IF(L.EQ.3) OPEN(N,FILE='INP/CFG-3.INP')
    IF(L.EQ.4) OPEN(N,FILE='INP/CFG-4.INP')
    IF(L.EQ.5) OPEN(N,FILE='INP/CFG-5.INP')
    IF(L.EQ.6) OPEN(N,FILE='INP/CFG-6.INP')
    IF(L.EQ.7) OPEN(N,FILE='INP/CFG-7.INP')
    IF(L.EQ.8) OPEN(N,FILE='INP/CFG-8.INP')
    IF(L.EQ.9) OPEN(N,FILE='INP/CFG-9.INP')
    IF(L.EQ.10) OPEN(N,FILE='INP/CFG-10.INP')
    IF(L.EQ.11) OPEN(N,FILE='INP/CFG-11.INP')
    IF(L.EQ.12) OPEN(N,FILE='INP/CFG-12.INP')
    IF(L.EQ.13) OPEN(N,FILE='INP/CFG-13.INP')
    IF(L.EQ.14) OPEN(N,FILE='INP/CFG-14.INP')

    RETURN
  END SUBROUTINE HCFGFILE
!!!#######################################################################
  SUBROUTINE CFGFILE(N, L)
    IMPLICIT NONE
    INTEGER N, L
!!!........................

    WRITE(*,*) ' CFG:: OPENING FILE DATA  FOR L = ', L
    
    IF(L.EQ.0) OPEN(N,FILE='inp/cfg-0.inp')
    IF(L.EQ.1) OPEN(N,FILE='inp/cfg-1.inp')
    IF(L.EQ.2) OPEN(N,FILE='inp/cfg-2.inp')
    IF(L.EQ.3) OPEN(N,FILE='inp/cfg-3.inp')
    IF(L.EQ.4) OPEN(N,FILE='inp/cfg-4.inp')
    IF(L.EQ.5) OPEN(N,FILE='inp/cfg-5.inp')
    IF(L.EQ.6) OPEN(N,FILE='inp/cfg-6.inp')
    IF(L.EQ.7) OPEN(N,FILE='inp/cfg-7.inp')
    IF(L.EQ.8) OPEN(N,FILE='inp/cfg-8.inp')
    IF(L.EQ.9) OPEN(N,FILE='inp/cfg-9.inp')
    IF(L.EQ.10) OPEN(N,FILE='inp/cfg-10.inp')
    IF(L.EQ.11) OPEN(N,FILE='inp/cfg-11.inp')
    IF(L.EQ.12) OPEN(N,FILE='inp/cfg-12.inp')
    IF(L.EQ.13) OPEN(N,FILE='inp/cfg-13.inp')
    IF(L.EQ.14) OPEN(N,FILE='inp/cfg-14.inp')

    RETURN
  END SUBROUTINE CFGFILE
!!!#####################################################################
  SUBROUTINE H2EFILE(N, L)
    IMPLICIT NONE
    INTEGER N,L
    
    IF(N.EQ.1.OR.N.EQ.2) THEN 

       WRITE(*,*) ' OPENING FILE DATA  FOR L = ', L
       WRITE(*,*) ' WRITE DATA FOR SINGLE IONIZATION SPECTRUM OF ATOMIC SYSTEM '

       IF(L.EQ.0) OPEN(N,FILE='DAT/H2ES-0.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.1) OPEN(N,FILE='DAT/H2ES-1.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.2) OPEN(N,FILE='DAT/H2ES-2.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.3) OPEN(N,FILE='DAT/H2ES-3.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.4) OPEN(N,FILE='DAT/H2ES-4.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.5) OPEN(N,FILE='DAT/H2ES-5.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.6) OPEN(N,FILE='DAT/H2ES-6.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.7) OPEN(N,FILE='DAT/H2ES-7.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.8) OPEN(N,FILE='DAT/H2ES-8.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.9) OPEN(N,FILE='DAT/H2ES-9.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.10) OPEN(N,FILE='DAT/H2ES-10.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.11) OPEN(N,FILE='DAT/H2ES-11.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       
    ELSE

       WRITE(*,*) ' WRITE DATA FOR DOUBLE IONIZATION SPECTRUM OF ATOMIC SYSTEM '
       
       IF(L.EQ.0) OPEN(N,FILE='DAT/H2ED-0.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.1) OPEN(N,FILE='DAT/H2ED-1.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.2) OPEN(N,FILE='DAT/H2ED-2.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.3) OPEN(N,FILE='DAT/H2ED-3.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.4) OPEN(N,FILE='DAT/H2ED-4.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.5) OPEN(N,FILE='DAT/H2ED-5.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.6) OPEN(N,FILE='DAT/H2ED-6.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.7) OPEN(N,FILE='DAT/H2ED-7.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.8) OPEN(N,FILE='DAT/H2ED-8.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.9) OPEN(N,FILE='DAT/H2ED-9.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.10) OPEN(N,FILE='DAT/H2ED-10.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       IF(L.EQ.11) OPEN(N,FILE='DAT/H2ED-11.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
       
    ENDIF


    RETURN
  END SUBROUTINE H2EFILE
!!!#####################################################################
  SUBROUTINE AMXFILE(N, L)
    IMPLICIT NONE
    INTEGER N,L
        
    IF(N.EQ.10.OR.N.EQ.11) THEN 

       WRITE(*,*) ' AMXFILE:: OPENING FILE DATA  FOR L = ', L
       WRITE(*,*) ' READ DATA FOR SINGLE IONIZATION SPECTRUM OF ATOMIC SYSTEM '

       IF(L.EQ.0) OPEN(N, FILE='DAT/AMXS-0.DAT')
       IF(L.EQ.1) OPEN(N, FILE='DAT/AMXS-1.DAT')
       IF(L.EQ.2) OPEN(N, FILE='DAT/AMXS-2.DAT')
       IF(L.EQ.3) OPEN(N, FILE='DAT/AMXS-3.DAT')
       IF(L.EQ.4) OPEN(N, FILE='DAT/AMXS-4.DAT')
       IF(L.EQ.5) OPEN(N, FILE='DAT/AMXS-5.DAT')
       IF(L.EQ.6) OPEN(N, FILE='DAT/AMXS-6.DAT')
       IF(L.EQ.7) OPEN(N, FILE='DAT/AMXS-7.DAT')
       IF(L.EQ.8) OPEN(N, FILE='DAT/AMXS-8.DAT')
       IF(L.EQ.9) OPEN(N, FILE='DAT/AMXS-9.DAT')
       IF(L.EQ.10) OPEN(N,FILE='DAT/AMXS-10.DAT')
       IF(L.EQ.11) OPEN(N,FILE='DAT/AMXS-11.DAT')
       
    ELSE
       WRITE(*,*) ' AMXFILE:: OPENING FILE DATA  FOR L = ', L
       WRITE(*,*) ' READ DATA FOR DOUBLE IONIZATION SPECTRUM OF ATOMIC SYSTEM '

       IF(L.EQ.0) OPEN(N, FILE='DAT/AMXD-0.DAT')
       IF(L.EQ.1) OPEN(N, FILE='DAT/AMXD-1.DAT')
       IF(L.EQ.2) OPEN(N, FILE='DAT/AMXD-2.DAT')
       IF(L.EQ.3) OPEN(N, FILE='DAT/AMXD-3.DAT')
       IF(L.EQ.4) OPEN(N, FILE='DAT/AMXD-4.DAT')
       IF(L.EQ.5) OPEN(N, FILE='DAT/AMXD-5.DAT')
       IF(L.EQ.6) OPEN(N, FILE='DAT/AMXD-6.DAT')
       IF(L.EQ.7) OPEN(N, FILE='DAT/AMXD-7.DAT')
       IF(L.EQ.8) OPEN(N, FILE='DAT/AMXD-8.DAT')
       IF(L.EQ.9) OPEN(N, FILE='DAT/AMXD-9.DAT')
       IF(L.EQ.10) OPEN(N,FILE='DAT/AMXD-10.DAT')
       IF(L.EQ.11) OPEN(N,FILE='DAT/AMXD-11.DAT')
       


    ENDIF

    RETURN
  END SUBROUTINE AMXFILE

!!!#####################################################################
    SUBROUTINE KMXFILE(N, L)
      IMPLICIT NONE
      INTEGER N,L
      
      IF(N.EQ.20.OR.N.EQ.21) THEN 

         WRITE(*,*) ' KMXFILE:: OPENING FILE DATA  FOR L = ', L
         WRITE(*,*) ' READ DATA FOR SINGLE IONIZATION SPECTRUM OF ATOMIC SYSTEM '


       IF(L.EQ.0) OPEN(N, FILE='DAT/KMXS-0.DAT')
       IF(L.EQ.1) OPEN(N, FILE='DAT/KMXS-1.DAT')
       IF(L.EQ.2) OPEN(N, FILE='DAT/KMXS-2.DAT')
       IF(L.EQ.3) OPEN(N, FILE='DAT/KMXS-3.DAT')
       IF(L.EQ.4) OPEN(N, FILE='DAT/KMXS-4.DAT')
       IF(L.EQ.5) OPEN(N, FILE='DAT/KMXS-5.DAT')
       IF(L.EQ.6) OPEN(N, FILE='DAT/KMXS-6.DAT')
       IF(L.EQ.7) OPEN(N, FILE='DAT/KMXS-7.DAT')
       IF(L.EQ.8) OPEN(N, FILE='DAT/KMXS-8.DAT')
       IF(L.EQ.9) OPEN(N, FILE='DAT/KMXS-9.DAT')
       IF(L.EQ.10) OPEN(N,FILE='DAT/KMXS-10.DAT')
       IF(L.EQ.11) OPEN(N,FILE='DAT/KMXS-11.DAT')
         
      ELSE
         WRITE(*,*) ' KMXFILE:: OPENING FILE DATA  FOR L = ', L
         WRITE(*,*) ' WRITE DATA FOR DOUBLE IONIZATION SPECTRUM OF ATOMIC SYSTEM '

       IF(L.EQ.0) OPEN(N, FILE='DAT/KMXD-0.DAT')
       IF(L.EQ.1) OPEN(N, FILE='DAT/KMXD-1.DAT')
       IF(L.EQ.2) OPEN(N, FILE='DAT/KMXD-2.DAT')
       IF(L.EQ.3) OPEN(N, FILE='DAT/KMXD-3.DAT')
       IF(L.EQ.4) OPEN(N, FILE='DAT/KMXD-4.DAT')
       IF(L.EQ.5) OPEN(N, FILE='DAT/KMXD-5.DAT')
       IF(L.EQ.6) OPEN(N, FILE='DAT/KMXD-6.DAT')
       IF(L.EQ.7) OPEN(N, FILE='DAT/KMXD-7.DAT')
       IF(L.EQ.8) OPEN(N, FILE='DAT/KMXD-8.DAT')
       IF(L.EQ.9) OPEN(N, FILE='DAT/KMXD-9.DAT')
       IF(L.EQ.10) OPEN(N,FILE='DAT/KMXD-10.DAT')
       IF(L.EQ.11) OPEN(N,FILE='DAT/KMXD-11.DAT')
         
      ENDIF

      RETURN
    END SUBROUTINE KMXFILE

!!!#####################################################################
      SUBROUTINE ZMXFILE(N, L)
        IMPLICIT NONE
        INTEGER N,L
        
        IF(N.EQ.30.OR.N.EQ.31) THEN 
           
           WRITE(*,*) ' ZMXFILE:: OPENING FILE DATA  FOR L = ', L
           WRITE(*,*) ' READ DATA FOR SINGLE IONIZATION SPECTRUM OF ATOMIC SYSTEM '

       IF(L.EQ.0) OPEN(N, FILE='DAT/ZMXS-0.DAT')
       IF(L.EQ.1) OPEN(N, FILE='DAT/ZMXS-1.DAT')
       IF(L.EQ.2) OPEN(N, FILE='DAT/ZMXS-2.DAT')
       IF(L.EQ.3) OPEN(N, FILE='DAT/ZMXS-3.DAT')
       IF(L.EQ.4) OPEN(N, FILE='DAT/ZMXS-4.DAT')
       IF(L.EQ.5) OPEN(N, FILE='DAT/ZMXS-5.DAT')
       IF(L.EQ.6) OPEN(N, FILE='DAT/ZMXS-6.DAT')
       IF(L.EQ.7) OPEN(N, FILE='DAT/ZMXS-7.DAT')
       IF(L.EQ.8) OPEN(N, FILE='DAT/ZMXS-8.DAT')
       IF(L.EQ.9) OPEN(N, FILE='DAT/ZMXS-9.DAT')
       IF(L.EQ.10) OPEN(N,FILE='DAT/ZMXS-10.DAT')
       IF(L.EQ.11) OPEN(N,FILE='DAT/ZMXS-11.DAT')
           
           
        ELSE
           
           WRITE(*,*) ' AMXFILE:: OPENING FILE DATA  FOR L = ', L
           WRITE(*,*) ' WRITE DATA FOR DOUBLE IONIZATION SPECTRUM OF ATOMIC SYSTEM '

       IF(L.EQ.0) OPEN(N, FILE='DAT/ZMXD-0.DAT')
       IF(L.EQ.1) OPEN(N, FILE='DAT/ZMXD-1.DAT')
       IF(L.EQ.2) OPEN(N, FILE='DAT/ZMXD-2.DAT')
       IF(L.EQ.3) OPEN(N, FILE='DAT/ZMXD-3.DAT')
       IF(L.EQ.4) OPEN(N, FILE='DAT/ZMXD-4.DAT')
       IF(L.EQ.5) OPEN(N, FILE='DAT/ZMXD-5.DAT')
       IF(L.EQ.6) OPEN(N, FILE='DAT/ZMXD-6.DAT')
       IF(L.EQ.7) OPEN(N, FILE='DAT/ZMXD-7.DAT')
       IF(L.EQ.8) OPEN(N, FILE='DAT/ZMXD-8.DAT')
       IF(L.EQ.9) OPEN(N, FILE='DAT/ZMXD-9.DAT')
       IF(L.EQ.10) OPEN(N,FILE='DAT/ZMXD-10.DAT')
       IF(L.EQ.11) OPEN(N,FILE='DAT/ZMXD-11.DAT')
           
        ENDIF


        RETURN
      END SUBROUTINE ZMXFILE

!!!#####################################################################
      SUBROUTINE RBBFILE(N, L)
        IMPLICIT NONE
        INTEGER N,L
        
        WRITE(*,*) '#                       HDMX1E:: < B| R | P_nl > DATA  FOR L = ', L

        IF(L.EQ.0) OPEN(N,FILE='DAT/RBB-0.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.1) OPEN(N,FILE='DAT/RBB-1.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.2) OPEN(N,FILE='DAT/RBB-2.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.3) OPEN(N,FILE='DAT/RBB-3.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.4) OPEN(N,FILE='DAT/RBB-4.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.5) OPEN(N,FILE='DAT/RBB-5.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.6) OPEN(N,FILE='DAT/RBB-6.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.7) OPEN(N,FILE='DAT/RBB-7.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.8) OPEN(N,FILE='DAT/RBB-8.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.9) OPEN(N,FILE='DAT/RBB-9.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.10) OPEN(N,FILE='DAT/RBB-10.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.11) OPEN(N,FILE='DAT/RBB-11.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.12) OPEN(N,FILE='DAT/RBB-12.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.13) OPEN(N,FILE='DAT/RBB-13.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF(L.EQ.14) OPEN(N,FILE='DAT/RBB-14.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        
        RETURN
      END SUBROUTINE RBBFILE

!!!##########################################################
      subroutine rwfbfile(n, l1, l2, mode)
        implicit none
        integer n,l1,l2,mode

      if(l1.gt.14) then

         write(*,*) '#               l must be less or equal 14  . for l > 14'&
            &' modify subroutine cs2phfile. exiting...'
         stop

      else

         write(*,*) '#                hdmx1e :: < b_i l | r | b_j, l+1 > data for ', l1, ' -> ', l2

      endif

      if(mode.eq.0) then

      if(l1.eq.0.and.l2.eq.1) open(n,file='dat/rwfb-01.v.dat', form='unformatted',access='sequential')
      if(l1.eq.1.and.l2.eq.0) open(n,file='dat/rwfb-10.v.dat', form='unformatted',access='sequential')

      if(l1.eq.1.and.l2.eq.2) open(n,file='dat/rwfb-12.v.dat', form='unformatted',access='sequential')
      if(l1.eq.2.and.l2.eq.1) open(n,file='dat/rwfb-21.v.dat', form='unformatted',access='sequential')

      if(l1.eq.2.and.l2.eq.3) open(n,file='dat/rwfb-23.v.dat', form='unformatted',access='sequential')
      if(l1.eq.3.and.l2.eq.2) open(n,file='dat/rwfb-32.v.dat', form='unformatted',access='sequential')

      if(l1.eq.3.and.l2.eq.4)  open(n,file='dat/rwfb-34.v.dat', form='unformatted',access='sequential')
      if(l1.eq.4.and.l2.eq.3)  open(n,file='dat/rwfb-43.v.dat', form='unformatted',access='sequential')

      if(l1.eq.4.and.l2.eq.5)  open(n,file='dat/rwfb-45.v.dat', form='unformatted',access='sequential')
      if(l1.eq.5.and.l2.eq.4)  open(n,file='dat/rwfb-54.v.dat', form='unformatted',access='sequential')

      if(l1.eq.5.and.l2.eq.6)  open(n,file='dat/rwfb-56.v.dat', form='unformatted',access='sequential')
      if(l1.eq.6.and.l2.eq.5)  open(n,file='dat/rwfb-65.v.dat', form='unformatted',access='sequential')

      if(l1.eq.6.and.l2.eq.7)  open(n,file='dat/rwfb-67.v.dat', form='unformatted',access='sequential')
      if(l1.eq.7.and.l2.eq.6)  open(n,file='dat/rwfb-76.v.dat', form='unformatted',access='sequential')

      if(l1.eq.7.and.l2.eq.8)  open(n,file='dat/rwfb-78.v.dat', form='unformatted',access='sequential')
      if(l1.eq.8.and.l2.eq.7)  open(n,file='dat/rwfb-87.v.dat', form='unformatted',access='sequential')

      if(l1.eq.8.and.l2.eq.9)  open(n,file='dat/rwfb-89.v.dat', form='unformatted',access='sequential')
      if(l1.eq.9.and.l2.eq.8)  open(n,file='dat/rwfb-98.v.dat', form='unformatted',access='sequential')

      if(l1.eq.9.and.l2.eq.10) open(n,file='dat/rwfb-910.v.dat', form='unformatted',access='sequential')
      if(l1.eq.10.and.l2.eq.9) open(n,file='dat/rwfb-109.v.dat', form='unformatted',access='sequential')

      if(l1.eq.10.and.l2.eq.11) open(n,file='dat/rwfb-1011.v.dat', form='unformatted',access='sequential')
      if(l1.eq.11.and.l2.eq.10) open(n,file='dat/rwfb-1110.v.dat', form='unformatted',access='sequential')

      if(l1.eq.11.and.l2.eq.12) open(n,file='dat/rwfb-1112.v.dat', form='unformatted',access='sequential')
      if(l1.eq.12.and.l2.eq.11) open(n,file='dat/rwfb-1211.v.dat', form='unformatted',access='sequential')

      if(l1.eq.12.and.l2.eq.13) open(n,file='dat/rwfb-1213.v.dat', form='unformatted',access='sequential')
      if(l1.eq.13.and.l2.eq.12) open(n,file='dat/rwfb-1312.v.dat', form='unformatted',access='sequential')

      if(l1.eq.13.and.l2.eq.14) open(n,file='dat/rwfb-1314.v.dat', form='unformatted',access='sequential')
      if(l1.eq.4.and.l2.eq.13)  open(n,file='dat/rwfb-1415.v.dat', form='unformatted',access='sequential')

      else
         
      if(l1.eq.0.and.l2.eq.1) open(n,file='dat/rwfb-01.l.dat', form='unformatted',access='sequential')
      if(l1.eq.1.and.l2.eq.0) open(n,file='dat/rwfb-10.l.dat', form='unformatted',access='sequential')

      if(l1.eq.1.and.l2.eq.2) open(n,file='dat/rwfb-12.l.dat', form='unformatted',access='sequential')
      if(l1.eq.2.and.l2.eq.1) open(n,file='dat/rwfb-21.l.dat', form='unformatted',access='sequential')

      if(l1.eq.2.and.l2.eq.3) open(n,file='dat/rwfb-23.l.dat', form='unformatted',access='sequential')
      if(l1.eq.3.and.l2.eq.2) open(n,file='dat/rwfb-32.l.dat', form='unformatted',access='sequential')

      if(l1.eq.3.and.l2.eq.4)   open(n,file='dat/rwfb-34.l.dat', form='unformatted',access='sequential')
      if(l1.eq.4.and.l2.eq.3)   open(n,file='dat/rwfb-43.l.dat', form='unformatted',access='sequential')

      if(l1.eq.4.and.l2.eq.5)   open(n,file='dat/rwfb-45.l.dat', form='unformatted',access='sequential')
      if(l1.eq.5.and.l2.eq.4)   open(n,file='dat/rwfb-54.l.dat', form='unformatted',access='sequential')

      if(l1.eq.5.and.l2.eq.6)   open(n,file='dat/rwfb-56.l.dat', form='unformatted',access='sequential')
      if(l1.eq.6.and.l2.eq.5)   open(n,file='dat/rwfb-65.l.dat', form='unformatted',access='sequential')

      if(l1.eq.6.and.l2.eq.7)   open(n,file='dat/rwfb-67.l.dat', form='unformatted',access='sequential')
      if(l1.eq.7.and.l2.eq.6)   open(n,file='dat/rwfb-76.l.dat', form='unformatted',access='sequential')

      if(l1.eq.7.and.l2.eq.8)   open(n,file='dat/rwfb-78.l.dat', form='unformatted',access='sequential')
      if(l1.eq.8.and.l2.eq.7)   open(n,file='dat/rwfb-87.l.dat', form='unformatted',access='sequential')

      if(l1.eq.8.and.l2.eq.9)   open(n,file='dat/rwfb-89.l.dat', form='unformatted',access='sequential')
      if(l1.eq.9.and.l2.eq.8)   open(n,file='dat/rwfb-98.l.dat', form='unformatted',access='sequential')

      if(l1.eq.9.and.l2.eq.10)  open(n,file='dat/rwfb-910.l.dat', form='unformatted',access='sequential')
      if(l1.eq.10.and.l2.eq.9)  open(n,file='dat/rwfb-109.l.dat', form='unformatted',access='sequential')

      if(l1.eq.10.and.l2.eq.11) open(n,file='dat/rwfb-1011.l.dat', form='unformatted',access='sequential')
      if(l1.eq.11.and.l2.eq.10) open(n,file='dat/rwfb-1110.l.dat', form='unformatted',access='sequential')

      if(l1.eq.11.and.l2.eq.12) open(n,file='dat/rwfb-1112.l.dat', form='unformatted',access='sequential')
      if(l1.eq.12.and.l2.eq.11) open(n,file='dat/rwfb-1211.l.dat', form='unformatted',access='sequential')

      if(l1.eq.12.and.l2.eq.13) open(n,file='dat/rwfb-1213.l.dat', form='unformatted',access='sequential')
      if(l1.eq.13.and.l2.eq.12) open(n,file='dat/rwfb-1312.l.dat', form='unformatted',access='sequential')

      if(l1.eq.13.and.l2.eq.14) open(n,file='dat/rwfb-1314.l.dat', form='unformatted',access='sequential')
      if(l1.eq.14.and.l2.eq.13)  open(n,file='dat/rwfb-1413.l.dat', form='unformatted',access='sequential')

   endif

   return

 end subroutine rwfbfile

!!!#########################################

      subroutine rwfbbfile(n, l, mode)

        implicit none
        integer n,l,mode

      if(l.gt.14) then

         write(*,*) '#               l must be less or equal 14  . for l > 14'&
            &' modify subroutine rwfbbfile. exiting...'
         stop

      else

         write(*,*) '#                   rwfbb :: < wf | r | b > data for ', l, ' -> ', l+1

      endif

      if(mode.eq.0) then

      if(l.eq.0)  open(n,file='dat/rwfbb-01.v.dat', form='unformatted',access='sequential')
      if(l.eq.1)  open(n,file='dat/rwfbb-12.v.dat', form='unformatted',access='sequential')
      if(l.eq.2)  open(n,file='dat/rwfbb-23.v.dat', form='unformatted',access='sequential')
      if(l.eq.3)  open(n,file='dat/rwfbb-34.v.dat', form='unformatted',access='sequential')
      if(l.eq.4)  open(n,file='dat/rwfbb-45.v.dat', form='unformatted',access='sequential')
      if(l.eq.5)  open(n,file='dat/rwfbb-56.v.dat', form='unformatted',access='sequential')
      if(l.eq.6)  open(n,file='dat/rwfbb-67.v.dat', form='unformatted',access='sequential')
      if(l.eq.7)  open(n,file='dat/rwfbb-78.v.dat', form='unformatted',access='sequential')
      if(l.eq.8)  open(n,file='dat/rwfbb-89.v.dat', form='unformatted',access='sequential')
      if(l.eq.9)  open(n,file='dat/rwfbb-910.v.dat', form='unformatted',access='sequential')
      if(l.eq.10) open(n,file='dat/rwfbb-1011.v.dat', form='unformatted',access='sequential')
      if(l.eq.11) open(n,file='dat/rwfbb-1112.v.dat', form='unformatted',access='sequential')
      if(l.eq.12) open(n,file='dat/rwfbb-1213.v.dat', form='unformatted',access='sequential')
      if(l.eq.13) open(n,file='dat/rwfbb-1314.v.dat', form='unformatted',access='sequential')

      else

      if(l.eq.0)  open(n,file='dat/rwfbb-01.l.dat', form='unformatted',access='sequential')
      if(l.eq.1)  open(n,file='dat/rwfbb-12.l.dat', form='unformatted',access='sequential')
      if(l.eq.2)  open(n,file='dat/rwfbb-23.l.dat', form='unformatted',access='sequential')
      if(l.eq.3)  open(n,file='dat/rwfbb-34.l.dat', form='unformatted',access='sequential')
      if(l.eq.4)  open(n,file='dat/rwfbb-45.l.dat', form='unformatted',access='sequential')
      if(l.eq.5)  open(n,file='dat/rwfbb-56.l.dat', form='unformatted',access='sequential')
      if(l.eq.6)  open(n,file='dat/rwfbb-67.l.dat', form='unformatted',access='sequential')
      if(l.eq.7)  open(n,file='dat/rwfbb-78.l.dat', form='unformatted',access='sequential')
      if(l.eq.8)  open(n,file='dat/rwfbb-89.l.dat', form='unformatted',access='sequential')
      if(l.eq.9)  open(n,file='dat/rwfbb-910.l.dat', form='unformatted',access='sequential')
      if(l.eq.10) open(n,file='dat/rwfbb-1011.l.dat', form='unformatted',access='sequential')
      if(l.eq.11) open(n,file='dat/rwfbb-1112.l.dat', form='unformatted',access='sequential')
      if(l.eq.12) open(n,file='dat/rwfbb-1213.l.dat', form='unformatted',access='sequential')
      if(l.eq.13) open(n,file='dat/rwfbb-1314.l.dat', form='unformatted',access='sequential')

   endif

   return

 end subroutine rwfbbfile
!!!##############################################################
      subroutine rwfwffile(n, l, mode)

        implicit none
        integer n,l,mode

      if(l.gt.14) then

         write(*,*) '#               l must be less or equal 14  . for l > 14'&
            &' modify subroutine rwfbbfile. exiting...'
         stop

      else

         write(*,*) '#         rwfwf :: < p(na,la ;r) | r | p(nb,lb;r) >  : ', l, ' -> ', l+1

      endif

      if(mode.eq.0) then

         if(l.eq.0)  open(n,file='dat/rwfwf-01.v.dat', form='unformatted',access='sequential')
         if(l.eq.1)  open(n,file='dat/rwfwf-12.v.dat', form='unformatted',access='sequential')
         if(l.eq.2)  open(n,file='dat/rwfwf-23.v.dat', form='unformatted',access='sequential')
         if(l.eq.3)  open(n,file='dat/rwfwf-34.v.dat', form='unformatted',access='sequential')
         if(l.eq.4)  open(n,file='dat/rwfwf-45.v.dat', form='unformatted',access='sequential')
         if(l.eq.5)  open(n,file='dat/rwfwf-56.v.dat', form='unformatted',access='sequential')
         if(l.eq.6)  open(n,file='dat/rwfwf-67.v.dat', form='unformatted',access='sequential')
         if(l.eq.7)  open(n,file='dat/rwfwf-78.v.dat', form='unformatted',access='sequential')
         if(l.eq.8)  open(n,file='dat/rwfwf-89.v.dat', form='unformatted',access='sequential')
         if(l.eq.9)  open(n,file='dat/rwfwf-910.v.dat', form='unformatted',access='sequential')
         if(l.eq.10) open(n,file='dat/rwfwf-1011.v.dat', form='unformatted',access='sequential')
         if(l.eq.11) open(n,file='dat/rwfwf-1112.v.dat', form='unformatted',access='sequential')
         if(l.eq.12) open(n,file='dat/rwfwf-1213.v.dat', form='unformatted',access='sequential')
         if(l.eq.13) open(n,file='dat/rwfwf-1314.v.dat', form='unformatted',access='sequential')

      else if(mode.eq.1) then

      if(l.eq.0)  open(n,file='dat/rwfwf-01.l.dat', form='unformatted',access='sequential')
      if(l.eq.1)  open(n,file='dat/rwfwf-12.l.dat', form='unformatted',access='sequential')
      if(l.eq.2)  open(n,file='dat/rwfwf-23.l.dat', form='unformatted',access='sequential')
      if(l.eq.3)  open(n,file='dat/rwfwf-34.l.dat', form='unformatted',access='sequential')
      if(l.eq.4)  open(n,file='dat/rwfwf-45.l.dat', form='unformatted',access='sequential')
      if(l.eq.5)  open(n,file='dat/rwfwf-56.l.dat', form='unformatted',access='sequential')
      if(l.eq.6)  open(n,file='dat/rwfwf-67.l.dat', form='unformatted',access='sequential')
      if(l.eq.7)  open(n,file='dat/rwfwf-78.l.dat', form='unformatted',access='sequential')
      if(l.eq.8)  open(n,file='dat/rwfwf-89.l.dat', form='unformatted',access='sequential')
      if(l.eq.9)  open(n,file='dat/rwfwf-910.l.dat', form='unformatted',access='sequential')
      if(l.eq.10) open(n,file='dat/rwfwf-1011.l.dat', form='unformatted',access='sequential')
      if(l.eq.11) open(n,file='dat/rwfwf-1112.l.dat', form='unformatted',access='sequential')
      if(l.eq.12) open(n,file='dat/rwfwf-1213.l.dat', form='unformatted',access='sequential')
      if(l.eq.13) open(n,file='dat/rwfwf-1314.l.dat', form='unformatted',access='sequential')


      else if(mode.eq.2) then

         if(l.eq.0)  open(n,file='dat/rwfwf-01.a.dat', form='unformatted',access='sequential')
         if(l.eq.1)  open(n,file='dat/rwfwf-12.a.dat', form='unformatted',access='sequential')
         if(l.eq.2)  open(n,file='dat/rwfwf-23.a.dat', form='unformatted',access='sequential')
         if(l.eq.3)  open(n,file='dat/rwfwf-34.a.dat', form='unformatted',access='sequential')
         if(l.eq.4)  open(n,file='dat/rwfwf-45.a.dat', form='unformatted',access='sequential')
         if(l.eq.5)  open(n,file='dat/rwfwf-56.a.dat', form='unformatted',access='sequential')
         if(l.eq.6)  open(n,file='dat/rwfwf-67.a.dat', form='unformatted',access='sequential')
         if(l.eq.7)  open(n,file='dat/rwfwf-78.a.dat', form='unformatted',access='sequential')
         if(l.eq.8)  open(n,file='dat/rwfwf-89.a.dat', form='unformatted',access='sequential')
         if(l.eq.9)  open(n,file='dat/rwfwf-910.a.dat', form='unformatted',access='sequential')
         if(l.eq.10) open(n,file='dat/rwfwf-1011.a.dat', form='unformatted',access='sequential')
         if(l.eq.11) open(n,file='dat/rwfwf-1112.a.dat', form='unformatted',access='sequential')
         if(l.eq.12) open(n,file='dat/rwfwf-1213.a.dat', form='unformatted',access='sequential')
         if(l.eq.13) open(n,file='dat/rwfwf-1314.a.dat', form='unformatted',access='sequential')         
   endif

   return

 end subroutine rwfwffile
!!!##############################################################
      SUBROUTINE HWF2EFILE(N, L)
       IMPLICIT NONE
       INTEGER N,L

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14',&
         &' MODIFY SUBROUTINE WF2EFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'WF2E :: DATA FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(UNIT=N,FILE='dat/wf2e-0.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/wf2e-1.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/wf2e-2.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/wf2e-3.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/wf2e-4.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/wf2e-5.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/wf2e-6.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/wf2e-7.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/wf2e-8.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/wf2e-9.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/wf2e-10.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/wf2e-11.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/wf2e-12.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/wf2e-13.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/wf2e-14.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
    END SUBROUTINE HWF2EFILE

!!#######################################################################

      SUBROUTINE HDMX1EFILE(N, L, MODE)

        IMPLICIT NONE
        INTEGER N, L,MODE

        IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14',&
            &' MODIFY SUBROUTINE DMX1EFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'DMX1E :: DATA FOR L = ', L

      ENDIF

      IF(MODE.EQ.0) THEN 

         
      IF(L.EQ.0) OPEN(N,FILE='dat/dmx1e-01.v.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/dmx1e-12.v.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/dmx1e-23.v.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/dmx1e-34.v.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/dmx1e-45.v.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/dmx1e-56.v.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/dmx1e-67.v.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/dmx1e-78.v.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/dmx1e-89.v.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/dmx1e-910.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/dmx1e-1011.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/dmx1e-1112.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/dmx1e-1213.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/dmx1e-1314.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/dmx1e-1415.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ELSE

      IF(L.EQ.0) OPEN(N,FILE='dat/dmx1e-01.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/dmx1e-12.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/dmx1e-23.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/dmx1e-34.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/dmx1e-45.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/dmx1e-56.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/dmx1e-67.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/dmx1e-78.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/dmx1e-89.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/dmx1e-910.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/dmx1e-1011.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/dmx1e-1112.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/dmx1e-1213.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/dmx1e-1314.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/dmx1e-1415.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ENDIF

      RETURN
    END SUBROUTINE HDMX1EFILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE HCS1PH1EFILE(N, L1, L2, MODE)

      IMPLICIT NONE      
      INTEGER N,L1,L2,MODE

!.......

      IF(L1.GT.5) THEN

         WRITE(*,*) ' L MUST BE LESS OR EQUAL 5  . FOR L > 14'&
              &' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
         STOP

      ELSE

         IF(MODE.EQ.0) THEN

            WRITE(*,*) 'HCS1E :: D_IF(J) = < I, L | D_v | F_J, L+1 > DATA FOR ', L1, ' -> ', L2
         ELSE

            WRITE(*,*) 'HCS1E :: D_IF(j) = < I, L | D_l | F_J J, L+1 > DATA FOR ', L1, ' -> ', L2
         END IF
            
      ENDIF

!........

      IF(MODE.EQ.0) THEN

         IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='dat/hcs1ph-01.v.dat')
         IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='dat/hcs1ph-10.v.dat')
         IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs1ph-12.v.dat')
         IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='dat/hcs1ph-21.v.dat')
         IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs1ph-23.v.dat')
         IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs1ph-32.v.dat')
         IF(L1.EQ.3.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs1ph-34.v.dat')
         IF(L1.EQ.4.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs1ph-43.v.dat')
         IF(L1.EQ.4.AND.L2.EQ.5) OPEN(N,FILE='dat/hcs1ph-45.v.dat')
         IF(L1.EQ.5.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs1ph-54.v.dat')

      ELSE IF(MODE.EQ.1) THEN 
         IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='dat/hcs1ph-01.l.dat')
         IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='dat/hcs1ph-10.l.dat')
         IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs1ph-12.l.dat')
         IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='dat/hcs1ph-21.l.dat')
         IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs1ph-23.l.dat')
         IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs1ph-32.l.dat')
         IF(L1.EQ.3.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs1ph-34.l.dat')
         IF(L1.EQ.4.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs1ph-43.l.dat')
         IF(L1.EQ.4.AND.L2.EQ.5) OPEN(N,FILE='dat/hcs1ph-45.l.dat')
         IF(L1.EQ.5.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs1ph-54.l.dat')

      ELSE IF(MODE.EQ.2) THEN 

         IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='dat/hcs1ph-01.a.dat')
         IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='dat/hcs1ph-10.a.dat')
         IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs1ph-12.a.dat')
         IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='dat/hcs1ph-21.a.dat')
         IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs1ph-23.a.dat')
         IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs1ph-32.a.dat')
         IF(L1.EQ.3.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs1ph-34.a.dat')
         IF(L1.EQ.4.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs1ph-43.a.dat')
         IF(L1.EQ.4.AND.L2.EQ.5) OPEN(N,FILE='dat/hcs1ph-45.a.dat')
         IF(L1.EQ.5.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs1ph-54.a.dat')

   ENDIF

   RETURN

 END SUBROUTINE HCS1PH1EFILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 
 SUBROUTINE HCS2PH1EFILE(N,L1,L2, MODE)

   IMPLICIT NONE
   INTEGER N, L1, L2, MODE
   
   IF(L1.GT.8) THEN
         
      WRITE(*,*) ' L MUST BE LESS OR EQUAL 8  . FOR L > 8', &
&     ' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
      STOP

   ELSE
      WRITE(*,*) '#          HCS2PH :: DATA FOR ', L1, ' -> ', L2
   ENDIF

   IF(MODE.EQ.0) THEN

      IF(L1.EQ.0.AND.L2.EQ.0) OPEN(N,FILE='dat/hcs2ph-00.v.dat')
      IF(L1.EQ.0.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs2ph-02.v.dat')
      IF(L1.EQ.1.AND.L2.EQ.1) OPEN(N,FILE='dat/hcs2ph-11.v.dat')
      IF(L1.EQ.1.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs2ph-13.v.dat')
      IF(L1.EQ.2.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs2ph-22.v.dat')
      IF(L1.EQ.2.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs2ph-24.v.dat')
      IF(L1.EQ.3.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs2ph-33.v.dat')
      IF(L1.EQ.3.AND.L2.EQ.5) OPEN(N,FILE='dat/hcs2ph-35.v.dat')
      IF(L1.EQ.4.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs2ph-44.v.dat')
      IF(L1.EQ.4.AND.L2.EQ.6) OPEN(N,FILE='dat/hcs2ph-46.v.dat')
      IF(L1.EQ.5.AND.L2.EQ.5) OPEN(N,FILE='dat/hcs2ph-55.v.dat')
      IF(L1.EQ.5.AND.L2.EQ.7) OPEN(N,FILE='dat/hcs2ph-57.v.dat')
      IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hcs2ph-66.v.dat')
      IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hcs2ph-68.v.dat')
      IF(L1.EQ.7.AND.L2.EQ.7) OPEN(N,FILE='dat/hcs2ph-77.v.dat')
      IF(L1.EQ.7.AND.L2.EQ.9) OPEN(N,FILE='dat/hcs2ph-79.v.dat')
      IF(L1.EQ.8.AND.L2.EQ.8) OPEN(N,FILE='dat/hcs2ph-88.v.dat')
      IF(L1.EQ.8.AND.L2.EQ.10) OPEN(N,FILE='dat/hcs2ph-810.v.dat')

      ELSE IF (MODE.EQ.1) THEN
         
         IF(L1.EQ.0.AND.L2.EQ.0) OPEN(N,FILE='dat/hcs2ph-00.l.dat')
         IF(L1.EQ.0.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs2ph-02.l.dat')
         IF(L1.EQ.1.AND.L2.EQ.1) OPEN(N,FILE='dat/hcs2ph-11.l.dat')
         IF(L1.EQ.1.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs2ph-13.l.dat')
         IF(L1.EQ.2.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs2ph-22.l.dat')
         IF(L1.EQ.2.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs2ph-24.l.dat')
         IF(L1.EQ.3.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs2ph-33.l.dat')
         IF(L1.EQ.3.AND.L2.EQ.5) OPEN(N,FILE='dat/hcs2ph-35.l.dat')
         IF(L1.EQ.4.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs2ph-44.l.dat')
         IF(L1.EQ.4.AND.L2.EQ.6) OPEN(N,FILE='dat/hcs2ph-46.l.dat')
         IF(L1.EQ.5.AND.L2.EQ.5) OPEN(N,FILE='dat/hcs2ph-55.l.dat')
         IF(L1.EQ.5.AND.L2.EQ.7) OPEN(N,FILE='dat/hcs2ph-57.l.dat')
         IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hcs2ph-66.l.dat')
         IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hcs2ph-68.l.dat')
         IF(L1.EQ.7.AND.L2.EQ.7) OPEN(N,FILE='dat/hcs2ph-77.l.dat')
         IF(L1.EQ.7.AND.L2.EQ.9) OPEN(N,FILE='dat/hcs2ph-79.l.dat')
         IF(L1.EQ.8.AND.L2.EQ.8) OPEN(N,FILE='dat/hcs2ph-88.l.dat')
         IF(L1.EQ.8.AND.L2.EQ.10) OPEN(N,FILE='dat/hcs2ph-810.l.dat')
         
      ELSE IF (MODE.EQ.2) THEN
         
         IF(L1.EQ.0.AND.L2.EQ.0) OPEN(N,FILE='dat/hcs2ph-00.a.dat')
         IF(L1.EQ.0.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs2ph-02.a.dat')
         IF(L1.EQ.1.AND.L2.EQ.1) OPEN(N,FILE='dat/hcs2ph-11.a.dat')
         IF(L1.EQ.1.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs2ph-13.a.dat')
         IF(L1.EQ.2.AND.L2.EQ.2) OPEN(N,FILE='dat/hcs2ph-22.a.dat')
         IF(L1.EQ.2.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs2ph-24.a.dat')
         IF(L1.EQ.3.AND.L2.EQ.3) OPEN(N,FILE='dat/hcs2ph-33.a.dat')
         IF(L1.EQ.3.AND.L2.EQ.5) OPEN(N,FILE='dat/hcs2ph-35.a.dat')
         IF(L1.EQ.4.AND.L2.EQ.4) OPEN(N,FILE='dat/hcs2ph-44.a.dat')
         IF(L1.EQ.4.AND.L2.EQ.6) OPEN(N,FILE='dat/hcs2ph-46.a.dat')
         IF(L1.EQ.5.AND.L2.EQ.5) OPEN(N,FILE='dat/hcs2ph-55.a.dat')
         IF(L1.EQ.5.AND.L2.EQ.7) OPEN(N,FILE='dat/hcs2ph-57.a.dat')
         IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hcs2ph-66.a.dat')
         IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hcs2ph-68.a.dat')
         IF(L1.EQ.7.AND.L2.EQ.7) OPEN(N,FILE='dat/hcs2ph-77.a.dat')
         IF(L1.EQ.7.AND.L2.EQ.9) OPEN(N,FILE='dat/hcs2ph-79.a.dat')
         IF(L1.EQ.8.AND.L2.EQ.8) OPEN(N,FILE='dat/hcs2ph-88.a.dat')
         IF(L1.EQ.8.AND.L2.EQ.10) OPEN(N,FILE='dat/hcs2ph-810.a.dat')
         
      ENDIF

   RETURN
 END SUBROUTINE HCS2PH1EFILE
!!!##########################################################
 SUBROUTINE HDMX2PHFILE(N, L1,L2, MODE)

   IMPLICIT NONE
   INTEGER N, L1, L2, MODE
   
   IF(L1.GT.8) THEN
         
      WRITE(*,*) ' L MUST BE LESS OR EQUAL 8  . FOR L > 8', &
&     ' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
      STOP

   ELSE
      WRITE(*,*) '#          HDMX2PH :: DATA FOR ', L1, ' -> ', L2
   ENDIF

   IF(MODE.EQ.0) THEN

      IF(L1.EQ.0.AND.L2.EQ.0) OPEN(N,FILE='dat/hdmx2ph-00.v.dat')
      IF(L1.EQ.0.AND.L2.EQ.2) OPEN(N,FILE='dat/hdmx2ph-02.v.dat')
      IF(L1.EQ.1.AND.L2.EQ.1) OPEN(N,FILE='dat/hdmx2ph-11.v.dat')
      IF(L1.EQ.1.AND.L2.EQ.3) OPEN(N,FILE='dat/hdmx2ph-13.v.dat')
      IF(L1.EQ.2.AND.L2.EQ.2) OPEN(N,FILE='dat/hdmx2ph-22.v.dat')
      IF(L1.EQ.2.AND.L2.EQ.4) OPEN(N,FILE='dat/hdmx2ph-24.v.dat')
      IF(L1.EQ.3.AND.L2.EQ.3) OPEN(N,FILE='dat/hdmx2ph-33.v.dat')
      IF(L1.EQ.3.AND.L2.EQ.5) OPEN(N,FILE='dat/hdmx2ph-35.v.dat')
      IF(L1.EQ.4.AND.L2.EQ.4) OPEN(N,FILE='dat/hdmx2ph-44.v.dat')
      IF(L1.EQ.4.AND.L2.EQ.6) OPEN(N,FILE='dat/hdmx2ph-46.v.dat')
      IF(L1.EQ.5.AND.L2.EQ.5) OPEN(N,FILE='dat/hdmx2ph-55.v.dat')
      IF(L1.EQ.5.AND.L2.EQ.7) OPEN(N,FILE='dat/hdmx2ph-57.v.dat')
      IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hdmx2ph-66.v.dat')
      IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hdmx2ph-68.v.dat')
      IF(L1.EQ.7.AND.L2.EQ.7) OPEN(N,FILE='dat/hdmx2ph-77.v.dat')
      IF(L1.EQ.7.AND.L2.EQ.9) OPEN(N,FILE='dat/hdmx2ph-79.v.dat')
      IF(L1.EQ.8.AND.L2.EQ.8) OPEN(N,FILE='dat/hdmx2ph-88.v.dat')
      IF(L1.EQ.8.AND.L2.EQ.10) OPEN(N,FILE='dat/hdmx2ph-810.v.dat')

      ELSE IF (MODE.EQ.1) THEN
         
         IF(L1.EQ.0.AND.L2.EQ.0) OPEN(N,FILE='dat/hdmx2ph-00.l.dat')
         IF(L1.EQ.0.AND.L2.EQ.2) OPEN(N,FILE='dat/hdmx2ph-02.l.dat')
         IF(L1.EQ.1.AND.L2.EQ.1) OPEN(N,FILE='dat/hdmx2ph-11.l.dat')
         IF(L1.EQ.1.AND.L2.EQ.3) OPEN(N,FILE='dat/hdmx2ph-13.l.dat')
         IF(L1.EQ.2.AND.L2.EQ.2) OPEN(N,FILE='dat/hdmx2ph-22.l.dat')
         IF(L1.EQ.2.AND.L2.EQ.4) OPEN(N,FILE='dat/hdmx2ph-24.l.dat')
         IF(L1.EQ.3.AND.L2.EQ.3) OPEN(N,FILE='dat/hdmx2ph-33.l.dat')
         IF(L1.EQ.3.AND.L2.EQ.5) OPEN(N,FILE='dat/hdmx2ph-35.l.dat')
         IF(L1.EQ.4.AND.L2.EQ.4) OPEN(N,FILE='dat/hdmx2ph-44.l.dat')
         IF(L1.EQ.4.AND.L2.EQ.6) OPEN(N,FILE='dat/hdmx2ph-46.l.dat')
         IF(L1.EQ.5.AND.L2.EQ.5) OPEN(N,FILE='dat/hdmx2ph-55.l.dat')
         IF(L1.EQ.5.AND.L2.EQ.7) OPEN(N,FILE='dat/hdmx2ph-57.l.dat')
         IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hdmx2ph-66.l.dat')
         IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hdmx2ph-68.l.dat')
         IF(L1.EQ.7.AND.L2.EQ.7) OPEN(N,FILE='dat/hdmx2ph-77.l.dat')
         IF(L1.EQ.7.AND.L2.EQ.9) OPEN(N,FILE='dat/hdmx2ph-79.l.dat')
         IF(L1.EQ.8.AND.L2.EQ.8) OPEN(N,FILE='dat/hdmx2ph-88.l.dat')
         IF(L1.EQ.8.AND.L2.EQ.10) OPEN(N,FILE='dat/hdmx2ph-810.l.dat')
         
      ELSE IF (MODE.EQ.2) THEN
         
         IF(L1.EQ.0.AND.L2.EQ.0) OPEN(N,FILE='dat/hdmx2ph-00.a.dat')
         IF(L1.EQ.0.AND.L2.EQ.2) OPEN(N,FILE='dat/hdmx2ph-02.a.dat')
         IF(L1.EQ.1.AND.L2.EQ.1) OPEN(N,FILE='dat/hdmx2ph-11.a.dat')
         IF(L1.EQ.1.AND.L2.EQ.3) OPEN(N,FILE='dat/hdmx2ph-13.a.dat')
         IF(L1.EQ.2.AND.L2.EQ.2) OPEN(N,FILE='dat/hdmx2ph-22.a.dat')
         IF(L1.EQ.2.AND.L2.EQ.4) OPEN(N,FILE='dat/hdmx2ph-24.a.dat')
         IF(L1.EQ.3.AND.L2.EQ.3) OPEN(N,FILE='dat/hdmx2ph-33.a.dat')
         IF(L1.EQ.3.AND.L2.EQ.5) OPEN(N,FILE='dat/hdmx2ph-35.a.dat')
         IF(L1.EQ.4.AND.L2.EQ.4) OPEN(N,FILE='dat/hdmx2ph-44.a.dat')
         IF(L1.EQ.4.AND.L2.EQ.6) OPEN(N,FILE='dat/hdmx2ph-46.a.dat')
         IF(L1.EQ.5.AND.L2.EQ.5) OPEN(N,FILE='dat/hdmx2ph-55.a.dat')
         IF(L1.EQ.5.AND.L2.EQ.7) OPEN(N,FILE='dat/hdmx2ph-57.a.dat')
         IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hdmx2ph-66.a.dat')
         IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/hdmx2ph-68.a.dat')
         IF(L1.EQ.7.AND.L2.EQ.7) OPEN(N,FILE='dat/hdmx2ph-77.a.dat')
         IF(L1.EQ.7.AND.L2.EQ.9) OPEN(N,FILE='dat/hdmx2ph-79.a.dat')
         IF(L1.EQ.8.AND.L2.EQ.8) OPEN(N,FILE='dat/hdmx2ph-88.a.dat')
         IF(L1.EQ.8.AND.L2.EQ.10) OPEN(N,FILE='dat/hdmx2ph-810.a.dat')
         
      ENDIF

   RETURN
 END SUBROUTINE HDMX2PHFILE
!!!##########################################################
      SUBROUTINE BFHDMX2EFILE(N, L1, L2, MODE)
        IMPLICIT NONE
        INTEGER N,L1,L2,MODE

      IF(L1.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS OR EQUAL 14  . FOR L > 14'&
            &' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
         STOP

      ELSE

         IF(MODE.EQ.0) THEN

            WRITE(*,*) 'BFHDMX2E :: D_IF(J) = < I, L | D_v | F_J, L+1 > DATA FOR ', L1, ' -> ', L2
         ELSE

            WRITE(*,*) 'BFHDMX2E :: D_IF(j) = < I, L | D_l | F_J J, L+1 > DATA FOR ', L1, ' -> ', L2
         END IF
            
      ENDIF

      IF(MODE.EQ.0) THEN

      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/HDMX2EBF-01.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/HDMX2EBF-10.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/HDMX2EBF-12.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/HDMX2EBF-21.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/HDMX2EBF-23.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/HDMX2EBF-32.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.4)  OPEN(N,FILE='DAT/HDMX2EBF-34.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3)  OPEN(N,FILE='DAT/HDMX2EBF-43.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.5)  OPEN(N,FILE='DAT/HDMX2EBF-45.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4)  OPEN(N,FILE='DAT/HDMX2EBF-54.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.6)  OPEN(N,FILE='DAT/HDMX2EBF-56.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5)  OPEN(N,FILE='DAT/HDMX2EBF-65.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.7)  OPEN(N,FILE='DAT/HDMX2EBF-67.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.6)  OPEN(N,FILE='DAT/HDMX2EBF-76.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.8)  OPEN(N,FILE='DAT/HDMX2EBF-78.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.7)  OPEN(N,FILE='DAT/HDMX2EBF-87.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.9)  OPEN(N,FILE='DAT/HDMX2EBF-89.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.9.AND.L2.EQ.8)  OPEN(N,FILE='DAT/HDMX2EBF-98.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.9.AND.L2.EQ.10) OPEN(N,FILE='DAT/HDMX2EBF-910.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.10.AND.L2.EQ.9) OPEN(N,FILE='DAT/HDMX2EBF-109.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.10.AND.L2.EQ.11) OPEN(N,FILE='DAT/HDMX2EBF-1011.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.11.AND.L2.EQ.10) OPEN(N,FILE='DAT/HDMX2EBF-1110.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.11.AND.L2.EQ.12) OPEN(N,FILE='DAT/HDMX2EBF-1112.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.12.AND.L2.EQ.11) OPEN(N,FILE='DAT/HDMX2EBF-1211.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.12.AND.L2.EQ.13) OPEN(N,FILE='DAT/HDMX2EBF-1213.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.13.AND.L2.EQ.12) OPEN(N,FILE='DAT/HDMX2EBF-1312.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.13.AND.L2.EQ.14) OPEN(N,FILE='DAT/HDMX2EBF-1314.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.13)  OPEN(N,FILE='DAT/HDMX2EBF-1415.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ELSE
         
      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/HDMX2EBF-01.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/HDMX2EBF-10.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/HDMX2EBF-12.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/HDMX2EBF-21.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/HDMX2EBF-23.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/HDMX2EBF-32.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.4)   OPEN(N,FILE='DAT/HDMX2EBF-34.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3)   OPEN(N,FILE='DAT/HDMX2EBF-43.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.5)   OPEN(N,FILE='DAT/HDMX2EBF-45.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4)   OPEN(N,FILE='DAT/HDMX2EBF-54.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.6)   OPEN(N,FILE='DAT/HDMX2EBF-56.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5)   OPEN(N,FILE='DAT/HDMX2EBF-65.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.7)   OPEN(N,FILE='DAT/HDMX2EBF-67.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.6)   OPEN(N,FILE='DAT/HDMX2EBF-76.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.8)   OPEN(N,FILE='DAT/HDMX2EBF-78.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.7)   OPEN(N,FILE='DAT/HDMX2EBF-87.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.9)   OPEN(N,FILE='DAT/HDMX2EBF-89.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.9.AND.L2.EQ.8)   OPEN(N,FILE='DAT/HDMX2EBF-98.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.9.AND.L2.EQ.10)  OPEN(N,FILE='DAT/HDMX2EBF-910.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.10.AND.L2.EQ.9)  OPEN(N,FILE='DAT/HDMX2EBF-109.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.10.AND.L2.EQ.11) OPEN(N,FILE='DAT/HDMX2EBF-1011.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.11.AND.L2.EQ.10) OPEN(N,FILE='DAT/HDMX2EBF-1110.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.11.AND.L2.EQ.12) OPEN(N,FILE='DAT/HDMX2EBF-1112.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.12.AND.L2.EQ.11) OPEN(N,FILE='DAT/HDMX2EBF-1211.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.12.AND.L2.EQ.13) OPEN(N,FILE='DAT/HDMX2EBF-1213.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.13.AND.L2.EQ.12) OPEN(N,FILE='DAT/HDMX2EBF-1312.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.13.AND.L2.EQ.14) OPEN(N,FILE='DAT/HDMX2EBF-1314.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.13)  OPEN(N,FILE='DAT/HDMX2EBF-1415.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

   ENDIF

   RETURN

 END SUBROUTINE BFHDMX2EFILE
!!!##########################################################
 SUBROUTINE BFDMX2EFILE(N, L1, L2, MODE)
   IMPLICIT NONE
   INTEGER N,L1,L2,MODE
   
   IF(L1.GT.14) THEN

      WRITE(*,*) ' L MUST BE LESS OR EQUAL 14  . FOR L > 14'&
           &' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
      STOP
         
   ELSE

      IF(MODE.EQ.0) THEN
            
            WRITE(*,*) 'NDMX2E :: D_IF(J) = < I, L | D_v | F_J, L+1 > DATA FOR ', L1, ' -> ', L2
         ELSE

            WRITE(*,*) 'NDMX2E :: D_IF(j) = < I, L | D_l | F_J J, L+1 > DATA FOR ', L1, ' -> ', L2
         END IF

      ENDIF

      IF(MODE.EQ.0) THEN

      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/DMX2EBF-01.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/DMX2EBF-10.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/DMX2EBF-12.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/DMX2EBF-21.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/DMX2EBF-23.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/DMX2EBF-32.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.4)  OPEN(N,FILE='DAT/DMX2EBF-34.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3)  OPEN(N,FILE='DAT/DMX2EBF-43.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.5)  OPEN(N,FILE='DAT/DMX2EBF-45.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4)  OPEN(N,FILE='DAT/DMX2EBF-54.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.6)  OPEN(N,FILE='DAT/DMX2EBF-56.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5)  OPEN(N,FILE='DAT/DMX2EBF-65.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.7)  OPEN(N,FILE='DAT/DMX2EBF-67.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.6)  OPEN(N,FILE='DAT/DMX2EBF-76.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.8)  OPEN(N,FILE='DAT/DMX2EBF-78.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.7)  OPEN(N,FILE='DAT/DMX2EBF-87.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.9)  OPEN(N,FILE='DAT/DMX2EBF-89.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.9.AND.L2.EQ.8)  OPEN(N,FILE='DAT/DMX2EBF-98.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.9.AND.L2.EQ.10) OPEN(N,FILE='DAT/DMX2EBF-910.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.10.AND.L2.EQ.9) OPEN(N,FILE='DAT/DMX2EBF-109.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.10.AND.L2.EQ.11) OPEN(N,FILE='DAT/DMX2EBF-1011.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.11.AND.L2.EQ.10) OPEN(N,FILE='DAT/DMX2EBF-1110.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.11.AND.L2.EQ.12) OPEN(N,FILE='DAT/DMX2EBF-1112.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.12.AND.L2.EQ.11) OPEN(N,FILE='DAT/DMX2EBF-1211.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.12.AND.L2.EQ.13) OPEN(N,FILE='DAT/DMX2EBF-1213.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.13.AND.L2.EQ.12) OPEN(N,FILE='DAT/DMX2EBF-1312.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.13.AND.L2.EQ.14) OPEN(N,FILE='DAT/DMX2EBF-1314.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.13)  OPEN(N,FILE='DAT/DMX2EBF-1415.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ELSE
         
      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/DMX2EBF-01.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/DMX2EBF-10.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/DMX2EBF-12.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/DMX2EBF-21.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/DMX2EBF-23.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/DMX2EBF-32.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.4)   OPEN(N,FILE='DAT/DMX2EBF-34.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3)   OPEN(N,FILE='DAT/DMX2EBF-43.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.5)   OPEN(N,FILE='DAT/DMX2EBF-45.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4)   OPEN(N,FILE='DAT/DMX2EBF-54.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.6)   OPEN(N,FILE='DAT/DMX2EBF-56.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5)   OPEN(N,FILE='DAT/DMX2EBF-65.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.7)   OPEN(N,FILE='DAT/DMX2EBF-67.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.6)   OPEN(N,FILE='DAT/DMX2EBF-76.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.8)   OPEN(N,FILE='DAT/DMX2EBF-78.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.7)   OPEN(N,FILE='DAT/DMX2EBF-87.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.9)   OPEN(N,FILE='DAT/DMX2EBF-89.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.9.AND.L2.EQ.8)   OPEN(N,FILE='DAT/DMX2EBF-98.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.9.AND.L2.EQ.10)  OPEN(N,FILE='DAT/DMX2EBF-910.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.10.AND.L2.EQ.9)  OPEN(N,FILE='DAT/DMX2EBF-109.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.10.AND.L2.EQ.11) OPEN(N,FILE='DAT/DMX2EBF-1011.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.11.AND.L2.EQ.10) OPEN(N,FILE='DAT/DMX2EBF-1110.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.11.AND.L2.EQ.12) OPEN(N,FILE='DAT/DMX2EBF-1112.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.12.AND.L2.EQ.11) OPEN(N,FILE='DAT/DMX2EBF-1211.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.12.AND.L2.EQ.13) OPEN(N,FILE='DAT/DMX2EBF-1213.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.13.AND.L2.EQ.12) OPEN(N,FILE='DAT/DMX2EBF-1312.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.13.AND.L2.EQ.14) OPEN(N,FILE='DAT/DMX2EBF-1314.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.13)  OPEN(N,FILE='DAT/DMX2EBF-1415.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

   ENDIF

   RETURN

 END SUBROUTINE BFDMX2EFILE
!!!##########################################################
 SUBROUTINE FFHDMX2EFILE(N, L1, L2, MODE)

   IMPLICIT NONE
   INTEGER N,L1,L2,MODE

   IF(L1.GT.14) THEN

      WRITE(*,*) ' L MUST BE LESS OR EQUAL 14  . FOR L > 14'&
           &' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
      STOP

      ELSE

         IF(MODE.EQ.0) THEN

            WRITE(*,*) 'FFHDMX2E :: D_IF(J) = < I, L | D_v | F_J, L+1 > DATA FOR ', L1, ' -> ', L2
         ELSE

            WRITE(*,*) 'FFHDMX2E :: D_IF(j) = < I, L | D_l | F_J J, L+1 > DATA FOR ', L1, ' -> ', L2
         END IF
            
      ENDIF

      IF(MODE.EQ.0) THEN

      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/HDMX2EFF-01.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/HDMX2EFF-10.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/HDMX2EFF-12.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/HDMX2EFF-21.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/HDMX2EFF-23.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/HDMX2EFF-32.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.4)  OPEN(N,FILE='DAT/HDMX2EFF-34.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3)  OPEN(N,FILE='DAT/HDMX2EFF-43.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.5)  OPEN(N,FILE='DAT/HDMX2EFF-45.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4)  OPEN(N,FILE='DAT/HDMX2EFF-54.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.6)  OPEN(N,FILE='DAT/HDMX2EFF-56.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5)  OPEN(N,FILE='DAT/HDMX2EFF-65.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.7)  OPEN(N,FILE='DAT/HDMX2EFF-67.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.6)  OPEN(N,FILE='DAT/HDMX2EFF-76.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.8)  OPEN(N,FILE='DAT/HDMX2EFF-78.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.7)  OPEN(N,FILE='DAT/HDMX2EFF-87.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.9)  OPEN(N,FILE='DAT/HDMX2EFF-89.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.9.AND.L2.EQ.8)  OPEN(N,FILE='DAT/HDMX2EFF-98.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.9.AND.L2.EQ.10) OPEN(N,FILE='DAT/HDMX2EFF-910.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.10.AND.L2.EQ.9) OPEN(N,FILE='DAT/HDMX2EFF-109.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.10.AND.L2.EQ.11) OPEN(N,FILE='DAT/HDMX2EFF-1011.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.11.AND.L2.EQ.10) OPEN(N,FILE='DAT/HDMX2EFF-1110.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.11.AND.L2.EQ.12) OPEN(N,FILE='DAT/HDMX2EFF-1112.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.12.AND.L2.EQ.11) OPEN(N,FILE='DAT/HDMX2EFF-1211.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.12.AND.L2.EQ.13) OPEN(N,FILE='DAT/HDMX2EFF-1213.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.13.AND.L2.EQ.12) OPEN(N,FILE='DAT/HDMX2EFF-1312.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.13.AND.L2.EQ.14) OPEN(N,FILE='DAT/HDMX2EFF-1314.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.13)  OPEN(N,FILE='DAT/HDMX2EFF-1415.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ELSE
         
      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/HDMX2EFF-01.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/HDMX2EFF-10.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/HDMX2EFF-12.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/HDMX2EFF-21.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/HDMX2EFF-23.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/HDMX2EFF-32.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.4)   OPEN(N,FILE='DAT/HDMX2EFF-34.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3)   OPEN(N,FILE='DAT/HDMX2EFF-43.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.5)   OPEN(N,FILE='DAT/HDMX2EFF-45.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4)   OPEN(N,FILE='DAT/HDMX2EFF-54.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.6)   OPEN(N,FILE='DAT/HDMX2EFF-56.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5)   OPEN(N,FILE='DAT/HDMX2EFF-65.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.7)   OPEN(N,FILE='DAT/HDMX2EFF-67.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.6)   OPEN(N,FILE='DAT/HDMX2EFF-76.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.8)   OPEN(N,FILE='DAT/HDMX2EFF-78.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.7)   OPEN(N,FILE='DAT/HDMX2EFF-87.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.9)   OPEN(N,FILE='DAT/HDMX2EFF-89.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.9.AND.L2.EQ.8)   OPEN(N,FILE='DAT/HDMX2EFF-98.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.9.AND.L2.EQ.10)  OPEN(N,FILE='DAT/HDMX2EFF-910.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.10.AND.L2.EQ.9)  OPEN(N,FILE='DAT/HDMX2EFF-109.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.10.AND.L2.EQ.11) OPEN(N,FILE='DAT/HDMX2EFF-1011.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.11.AND.L2.EQ.10) OPEN(N,FILE='DAT/HDMX2EFF-1110.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.11.AND.L2.EQ.12) OPEN(N,FILE='DAT/HDMX2EFF-1112.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.12.AND.L2.EQ.11) OPEN(N,FILE='DAT/HDMX2EFF-1211.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.12.AND.L2.EQ.13) OPEN(N,FILE='DAT/HDMX2EFF-1213.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.13.AND.L2.EQ.12) OPEN(N,FILE='DAT/HDMX2EFF-1312.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.13.AND.L2.EQ.14) OPEN(N,FILE='DAT/HDMX2EFF-1314.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.13)  OPEN(N,FILE='DAT/HDMX2EFF-1415.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

   ENDIF

   RETURN

 END SUBROUTINE FFHDMX2EFILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
 SUBROUTINE FFDMX2EFILE(N, L1, L2, MODE)
   IMPLICIT NONE
   INTEGER N,L1,L2,MODE
   
   IF(L1.GT.14) THEN

      WRITE(*,*) ' L MUST BE LESS OR EQUAL 14  . FOR L > 14'&
           &' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
      STOP
         
   ELSE

      IF(MODE.EQ.0) THEN
            
            WRITE(*,*) 'NDMX2E :: D_IF(J) = < I, L | D_v | F_J, L+1 > DATA FOR ', L1, ' -> ', L2
         ELSE

            WRITE(*,*) 'NDMX2E :: D_IF(j) = < I, L | D_l | F_J J, L+1 > DATA FOR ', L1, ' -> ', L2
         END IF

      ENDIF

      IF(MODE.EQ.0) THEN

      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/DMX2EFF-01.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/DMX2EFF-10.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/DMX2EFF-12.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/DMX2EFF-21.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/DMX2EFF-23.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/DMX2EFF-32.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.4)  OPEN(N,FILE='DAT/DMX2EFF-34.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3)  OPEN(N,FILE='DAT/DMX2EFF-43.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.5)  OPEN(N,FILE='DAT/DMX2EFF-45.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4)  OPEN(N,FILE='DAT/DMX2EFF-54.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.6)  OPEN(N,FILE='DAT/DMX2EFF-56.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5)  OPEN(N,FILE='DAT/DMX2EFF-65.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.7)  OPEN(N,FILE='DAT/DMX2EFF-67.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.6)  OPEN(N,FILE='DAT/DMX2EFF-76.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.8)  OPEN(N,FILE='DAT/DMX2EFF-78.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.7)  OPEN(N,FILE='DAT/DMX2EFF-87.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.9)  OPEN(N,FILE='DAT/DMX2EFF-89.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.9.AND.L2.EQ.8)  OPEN(N,FILE='DAT/DMX2EFF-98.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.9.AND.L2.EQ.10) OPEN(N,FILE='DAT/DMX2EFF-910.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.10.AND.L2.EQ.9) OPEN(N,FILE='DAT/DMX2EFF-109.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.10.AND.L2.EQ.11) OPEN(N,FILE='DAT/DMX2EFF-1011.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.11.AND.L2.EQ.10) OPEN(N,FILE='DAT/DMX2EFF-1110.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.11.AND.L2.EQ.12) OPEN(N,FILE='DAT/DMX2EFF-1112.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.12.AND.L2.EQ.11) OPEN(N,FILE='DAT/DMX2EFF-1211.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.12.AND.L2.EQ.13) OPEN(N,FILE='DAT/DMX2EFF-1213.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.13.AND.L2.EQ.12) OPEN(N,FILE='DAT/DMX2EFF-1312.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.13.AND.L2.EQ.14) OPEN(N,FILE='DAT/DMX2EFF-1314.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.13)  OPEN(N,FILE='DAT/DMX2EFF-1415.V.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ELSE
         
      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/DMX2EFF-01.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/DMX2EFF-10.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/DMX2EFF-12.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/DMX2EFF-21.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/DMX2EFF-23.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/DMX2EFF-32.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.4)   OPEN(N,FILE='DAT/DMX2EFF-34.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3)   OPEN(N,FILE='DAT/DMX2EFF-43.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.5)   OPEN(N,FILE='DAT/DMX2EFF-45.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4)   OPEN(N,FILE='DAT/DMX2EFF-54.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.6)   OPEN(N,FILE='DAT/DMX2EFF-56.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5)   OPEN(N,FILE='DAT/DMX2EFF-65.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.7)   OPEN(N,FILE='DAT/DMX2EFF-67.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.6)   OPEN(N,FILE='DAT/DMX2EFF-76.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.8)   OPEN(N,FILE='DAT/DMX2EFF-78.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.7)   OPEN(N,FILE='DAT/DMX2EFF-87.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.9)   OPEN(N,FILE='DAT/DMX2EFF-89.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.9.AND.L2.EQ.8)   OPEN(N,FILE='DAT/DMX2EFF-98.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.9.AND.L2.EQ.10)  OPEN(N,FILE='DAT/DMX2EFF-910.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.10.AND.L2.EQ.9)  OPEN(N,FILE='DAT/DMX2EFF-109.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.10.AND.L2.EQ.11) OPEN(N,FILE='DAT/DMX2EFF-1011.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.11.AND.L2.EQ.10) OPEN(N,FILE='DAT/DMX2EFF-1110.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.11.AND.L2.EQ.12) OPEN(N,FILE='DAT/DMX2EFF-1112.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.12.AND.L2.EQ.11) OPEN(N,FILE='DAT/DMX2EFF-1211.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.12.AND.L2.EQ.13) OPEN(N,FILE='DAT/DMX2EFF-1213.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.13.AND.L2.EQ.12) OPEN(N,FILE='DAT/DMX2EFF-1312.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.13.AND.L2.EQ.14) OPEN(N,FILE='DAT/DMX2EFF-1314.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.13)  OPEN(N,FILE='DAT/DMX2EFF-1415.L.DAT', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

   ENDIF

   RETURN

 END SUBROUTINE FFDMX2EFILE
!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 SUBROUTINE BBDMX2EFILE(N, L, MODE)

   IMPLICIT NONE
   INTEGER MODE, N,L

   IF(L.GT.14) THEN

      WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', &
           &       ' MODIFY SUBROUTINE DMX2EFILE. EXITING...'
      STOP

   ELSE

      WRITE(*,*) 'DMX2E :: DATA FOR L = ', L

   ENDIF
         
   IF(MODE.EQ.0) THEN

      IF(L.EQ.0) OPEN(N,FILE='dat/dmx2e-01.v.dat',FORM ='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/dmx2e-12.v.dat',FORM ='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/dmx2e-23.v.dat',FORM ='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/dmx2e-34.v.dat',FORM ='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/dmx2e-45.v.dat',FORM ='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/dmx2e-56.v.dat',FORM ='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/dmx2e-67.v.dat',FORM ='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/dmx2e-78.v.dat',FORM ='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/dmx2e-89.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/dmx2e-910.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/dmx2e-1011.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/dmx2e-1112.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/dmx2e-1213.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/dmx2e-1314.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/dmx2e-1415.v.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      
      ELSE
         
         IF(L.EQ.0) OPEN(N,FILE='dat/dmx2e-01.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.1) OPEN(N,FILE='dat/dmx2e-12.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.2) OPEN(N,FILE='dat/dmx2e-23.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.3) OPEN(N,FILE='dat/dmx2e-34.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.4) OPEN(N,FILE='dat/dmx2e-45.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.5) OPEN(N,FILE='dat/dmx2e-56.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.6) OPEN(N,FILE='dat/dmx2e-67.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.7) OPEN(N,FILE='dat/dmx2e-78.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.8) OPEN(N,FILE='dat/dmx2e-89.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.9) OPEN(N,FILE='dat/dmx2e-910.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.10) OPEN(N,FILE='dat/dmx2e-1011.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.11) OPEN(N,FILE='dat/dmx2e-1112.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.12) OPEN(N,FILE='dat/dmx2e-1213.l.dat', FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.13) OPEN(N,FILE='dat/dmx2e-1314.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         IF(L.EQ.14) OPEN(N,FILE='dat/dmx2e-1415.l.dat',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         
      ENDIF
      
      RETURN
    END SUBROUTINE BBDMX2EFILE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CS1PH2EFILE(N, L1, L2, MODE)
        IMPLICIT NONE
        INTEGER N,L1,L2,MODE

      IF(L1.GT.5) THEN

         WRITE(*,*) ' L MUST BE LESS OR EQUAL 5  . FOR L > 14'&
            &' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
         STOP

      ELSE

         IF(MODE.EQ.0) THEN

  WRITE(*,*) 'HDMX2E :: D_IF(J) = < I, L | D_v | F_J, L+1 > DATA FOR ', L1, ' -> ', L2
         ELSE

  WRITE(*,*) 'HDMX2E :: D_IF(j) = < I, L | D_l | F_J J, L+1 > DATA FOR ', L1, ' -> ', L2
         END IF
            
      ENDIF

      IF(MODE.EQ.0) THEN

         IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/CS1PH-01.V.DAT')
         IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/CS1PH-10.V.DAT')
         IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/CS1PH-12.V.DAT')
         IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/CS1PH-21.V.DAT')
         IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/CS1PH-23.V.DAT')
         IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/CS1PH-32.V.DAT')
         IF(L1.EQ.3.AND.L2.EQ.4) OPEN(N,FILE='DAT/CS1PH-34.V.DAT')
         IF(L1.EQ.4.AND.L2.EQ.3) OPEN(N,FILE='DAT/CS1PH-43.V.DAT')
         IF(L1.EQ.4.AND.L2.EQ.5) OPEN(N,FILE='DAT/CS1PH-45.V.DAT')
         IF(L1.EQ.5.AND.L2.EQ.4) OPEN(N,FILE='DAT/CS1PH-54.V.DAT')


      ELSE

         IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='DAT/CS1PH-01.L.DAT')
         IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='DAT/CS1PH-10.L.DAT')
         IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='DAT/CS1PH-12.L.DAT')
         IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='DAT/CS1PH-21.L.DAT')
         IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='DAT/CS1PH-23.L.DAT')
         IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='DAT/CS1PH-32.L.DAT')
         IF(L1.EQ.3.AND.L2.EQ.4) OPEN(N,FILE='DAT/CS1PH-34.L.DAT')
         IF(L1.EQ.4.AND.L2.EQ.3) OPEN(N,FILE='DAT/CS1PH-43.L.DAT')
         IF(L1.EQ.4.AND.L2.EQ.5) OPEN(N,FILE='DAT/CS1PH-45.L.DAT')
         IF(L1.EQ.5.AND.L2.EQ.4) OPEN(N,FILE='DAT/CS1PH-54.L.DAT')


   ENDIF

   RETURN

 END SUBROUTINE CS1PH2EFILE

!!!#######################################################################

END MODULE ioroutines
!!!##############################################################
!!!EOF
