C#
C#   d1e.f I/O files
C#

C#
C#    Open binary file for storing Energies and coefficient data
C################################################################

      SUBROUTINE D1EFILE(N,L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE D1EFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'D1E :: OPEN COEFF FILE FOR L = ', L

      ENDIF
         
         
      IF(L.EQ.0) OPEN(N,FILE='dat/d1e-0.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/d1e-1.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/d1e-2.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/d1e-3.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/d1e-4.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/d1e-5.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/d1e-6.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/d1e-7.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/d1e-8.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/d1e-9.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/d1e-10.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/d1e-11.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/d1e-12.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/d1e-13.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/d1e-14.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
      END
C#######################################################################
      SUBROUTINE HMXFILE(N,L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE HMXFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'HMXFILE :: OPEN COEFF FILE FOR L = ', L

      ENDIF
         
         
      IF(L.EQ.0) OPEN(N,FILE='dat/Hmx-0.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/Hmx-1.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/Hmx-2.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/Hmx-3.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/Hmx-4.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/Hmx-5.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/Hmx-6.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/Hmx-7.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/Hmx-8.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/Hmx-9.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/Hmx-10.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/Hmx-11.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/Hmx-12.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/Hmx-13.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/Hmx-14.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
      END
C#######################################################################
      SUBROUTINE BMXFILE(N,L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE HMXFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'BMXFILE :: OPEN COEFF FILE FOR L = ', L

      ENDIF
         
         
      IF(L.EQ.0) OPEN(N,FILE='dat/Bmx-0.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/Bmx-1.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/Bmx-2.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/Bmx-3.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/Bmx-4.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/Bmx-5.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/Bmx-6.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/Bmx-7.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/Bmx-8.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/Bmx-9.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/Bmx-10.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/Bmx-11.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/Bmx-12.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/Bmx-13.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/Bmx-14.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
      END
C#######################################################################
C#
C#    Open binary files of core symmetries Energies and coefficient data
C#######################################################################
      SUBROUTINE COREFILE(N,L)


      IF(L.GT.4) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 4 . FOR L > 4 ',
     1    'MODIFY SUBROUTINE COREFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'D1E :: OPEN CORE FILE FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(N,FILE='dat/core-0.dat',STATUS='OLD',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/core-1.dat',STATUS='OLD',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/core-2.dat',STATUS='OLD',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/core-3.dat',STATUS='OLD',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/core-4.dat',STATUS='OLD',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      RETURN 
      END
C#######################################################################
C#
C#   wf1e.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################

      SUBROUTINE WF1EFILE(N,L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE WF1EFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'WF1E  :: DATA FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(N,FILE='dat/wf1e-0.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/wf1e-1.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/wf1e-2.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/wf1e-3.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/wf1e-4.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/wf1e-5.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/wf1e-6.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/wf1e-7.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/wf1e-8.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/wf1e-9.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/wf1e-10.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/wf1e-11.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/wf1e-12.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/wf1e-13.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/wf1e-14.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
      END
C#######################################################################
C#
C#   r12.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################
      SUBROUTINE R12FILE(N, L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE R12FILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'R12 :: DATA FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(N,FILE='dat/r12-0.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/r12-1.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/r12-2.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/r12-3.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/r12-4.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/r12-5.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/r12-6.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/r12-7.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/r12-8.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/r12-9.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/r12-10.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/r12-11.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/r12-12.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/r12-13.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/r12-14.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
      END
C#######################################################################
C#
C#   r12.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################
      SUBROUTINE CFGFILE(N, L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 15 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE CFGFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'CFG :: CONFIGURATIONS FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(N, FILE ='inp/cfg-0.inp')
      IF(L.EQ.1) OPEN(N, FILE ='inp/cfg-1.inp')
      IF(L.EQ.2) OPEN(N, FILE ='inp/cfg-2.inp')
      IF(L.EQ.3) OPEN(N, FILE ='inp/cfg-3.inp')
      IF(L.EQ.4) OPEN(N, FILE ='inp/cfg-4.inp')
      IF(L.EQ.5) OPEN(N, FILE ='inp/cfg-5.inp')
      IF(L.EQ.6) OPEN(N, FILE ='inp/cfg-6.inp')
      IF(L.EQ.7) OPEN(N, FILE ='inp/cfg-7.inp')
      IF(L.EQ.8) OPEN(N, FILE ='inp/cfg-8.inp')
      IF(L.EQ.9) OPEN(N, FILE ='inp/cfg-9.inp')
      IF(L.EQ.10) OPEN(N,FILE ='inp/cfg-10.inp')
      IF(L.EQ.11) OPEN(N,FILE ='inp/cfg-11.inp')
      IF(L.EQ.12) OPEN(N,FILE ='inp/cfg-12.inp')
      IF(L.EQ.13) OPEN(N,FILE ='inp/cfg-13.inp')
      IF(L.EQ.14) OPEN(N,FILE ='inp/cfg-14.inp')

      RETURN
      END
C#######################################################################
C#
C#   d2e.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################
      SUBROUTINE D2EFILE(N, L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE D2EFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'D2E :: DATA FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(N,FILE='dat/d2e-0.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/d2e-1.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/d2e-2.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/d2e-3.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/d2e-4.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/d2e-5.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/d2e-6.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/d2e-7.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/d2e-8.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/d2e-9.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/d2e-10.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/d2e-11.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/d2e-12.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/d2e-13.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/d2e-14.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
      END
C#######################################################################
C#
C#   wf2e.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################

      SUBROUTINE WF2EFILE(N, L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE WF2EFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'WF2E :: DATA FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(N,FILE='dat/wf2e-0.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/wf2e-1.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/wf2e-2.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/wf2e-3.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/wf2e-4.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/wf2e-5.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/wf2e-6.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/wf2e-7.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/wf2e-8.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/wf2e-9.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/wf2e-10.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/wf2e-11.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/wf2e-12.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/wf2e-13.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/wf2e-14.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      RETURN
      END
C#######################################################################
C#
C#   dmx1e.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################

      SUBROUTINE DMX1EFILE(N, L, MODE)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE DMX1EFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'DMX1E :: DATA FOR L = ', L

      ENDIF

      IF(MODE.EQ.0) THEN 

         
      IF(L.EQ.0) OPEN(N,FILE='dat/dmx1e-01.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/dmx1e-12.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/dmx1e-23.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/dmx1e-34.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/dmx1e-45.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/dmx1e-56.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/dmx1e-67.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/dmx1e-78.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/dmx1e-89.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/dmx1e-910.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/dmx1e-1011.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/dmx1e-1112.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/dmx1e-1213.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/dmx1e-1314.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/dmx1e-1415.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ELSE

      IF(L.EQ.0) OPEN(N,FILE='dat/dmx1e-01.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/dmx1e-12.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/dmx1e-23.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/dmx1e-34.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/dmx1e-45.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/dmx1e-56.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/dmx1e-67.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/dmx1e-78.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/dmx1e-89.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/dmx1e-910.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/dmx1e-1011.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/dmx1e-1112.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/dmx1e-1213.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/dmx1e-1314.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/dmx1e-1415.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ENDIF

      RETURN
      END
C#######################################################################
C#
C#   dmx2e.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################

      SUBROUTINE DMX2EFILE(N, L, MODE)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE DMX2EFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'DMX2E :: DATA FOR L = ', L

      ENDIF
         
      IF(MODE.EQ.0) THEN

      IF(L.EQ.0) OPEN(N,FILE='dat/dmx2e-01.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/dmx2e-12.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/dmx2e-23.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/dmx2e-34.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/dmx2e-45.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/dmx2e-56.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/dmx2e-67.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/dmx2e-78.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/dmx2e-89.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/dmx2e-910.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/dmx2e-1011.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/dmx2e-1112.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/dmx2e-1213.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/dmx2e-1314.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/dmx2e-1415.v.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ELSE

         IF(L.EQ.0) OPEN(N,FILE='dat/dmx2e-01.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/dmx2e-12.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/dmx2e-23.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/dmx2e-34.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/dmx2e-45.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/dmx2e-56.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/dmx2e-67.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/dmx2e-78.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/dmx2e-89.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='dat/dmx2e-910.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='dat/dmx2e-1011.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='dat/dmx2e-1112.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='dat/dmx2e-1213.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='dat/dmx2e-1314.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='dat/dmx2e-1415.l.dat',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ENDIF
         
      RETURN
      END
C#######################################################################
C#
C#   cs1ph.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################
      SUBROUTINE CS1PHFILE(N, L, MODE)

      IF(L.GT.8) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 15. FOR L > 14', 
     1       ' MODIFY SUBROUTINE CS1PHFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'CS1PH :: DATA FOR ', L,' -> ', L + 1

      ENDIF

      IF(MODE.EQ.0) THEN
         
      IF(L.EQ.0) OPEN(N,FILE='dat/cs1ph-01.v.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/cs1ph-12.v.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/cs1ph-23.v.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/cs1ph-34.v.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/cs1ph-45.v.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/cs1ph-56.v.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/cs1ph-67.v.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/cs1ph-78.v.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/cs1ph-89.v.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ELSE

      IF(L.EQ.0) OPEN(N,FILE='dat/cs1ph-01.l.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='dat/cs1ph-12.l.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='dat/cs1ph-23.l.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='dat/cs1ph-34.l.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='dat/cs1ph-45.l.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='dat/cs1ph-56.l.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='dat/cs1ph-67.l.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='dat/cs1ph-78.l.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='dat/cs1ph-89.l.dat')
c     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      
      ENDIF

      RETURN
      END
C#######################################################################
C#
C#   dmx2e.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################
 
      SUBROUTINE CS2PHFILE(N,L1,L2, MODE)

      IF(L1.GT.8) THEN

         WRITE(*,*) ' L MUST BE LESS OR EQUAL 8  . FOR L > 8', 
     1       ' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'CS2PH :: DATA FOR ', L1, ' -> ', L2

      ENDIF

      IF(MODE.EQ.0) THEN

      IF(L1.EQ.0.AND.L2.EQ.0) OPEN(N,FILE='dat/cs2ph-00.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.0.AND.L2.EQ.2) OPEN(N,FILE='dat/cs2ph-02.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.1) OPEN(N,FILE='dat/cs2ph-11.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.3) OPEN(N,FILE='dat/cs2ph-13.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.2) OPEN(N,FILE='dat/cs2ph-22.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.4) OPEN(N,FILE='dat/cs2ph-24.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.3) OPEN(N,FILE='dat/cs2ph-33.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.5) OPEN(N,FILE='dat/cs2ph-35.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.4) OPEN(N,FILE='dat/cs2ph-44.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.6) OPEN(N,FILE='dat/cs2ph-46.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.5) OPEN(N,FILE='dat/cs2ph-55.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.7) OPEN(N,FILE='dat/cs2ph-57.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/cs2ph-66.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/cs2ph-68.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.7) OPEN(N,FILE='dat/cs2ph-77.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.9) OPEN(N,FILE='dat/cs2ph-79.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.8) OPEN(N,FILE='dat/cs2ph-88.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.10) OPEN(N,FILE='dat/cs2ph-810.v.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ELSE
         
      IF(L1.EQ.0.AND.L2.EQ.0) OPEN(N,FILE='dat/cs2ph-00.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.0.AND.L2.EQ.2) OPEN(N,FILE='dat/cs2ph-02.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.1.AND.L2.EQ.1) OPEN(N,FILE='dat/cs2ph-11.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.3) OPEN(N,FILE='dat/cs2ph-13.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.2.AND.L2.EQ.2) OPEN(N,FILE='dat/cs2ph-22.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.4) OPEN(N,FILE='dat/cs2ph-24.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.3.AND.L2.EQ.3) OPEN(N,FILE='dat/cs2ph-33.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.5) OPEN(N,FILE='dat/cs2ph-35.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.4.AND.L2.EQ.4) OPEN(N,FILE='dat/cs2ph-44.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.6) OPEN(N,FILE='dat/cs2ph-46.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.5.AND.L2.EQ.5) OPEN(N,FILE='dat/cs2ph-55.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.7) OPEN(N,FILE='dat/cs2ph-57.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/cs2ph-66.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.8) OPEN(N,FILE='dat/cs2ph-68.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.7.AND.L2.EQ.7) OPEN(N,FILE='dat/cs2ph-77.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.7.AND.L2.EQ.9) OPEN(N,FILE='dat/cs2ph-79.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      IF(L1.EQ.8.AND.L2.EQ.8) OPEN(N,FILE='dat/cs2ph-88.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.8.AND.L2.EQ.10) OPEN(N,FILE='dat/cs2ph-810.l.dat')
c     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      ENDIF

      RETURN
      END
C########################################################################


C#
C#               MANY CHANNEL PROJECT ( FREE BOUNDARY CONDITIONS ) 
C#

C#######################################################################
C#
C#   r12.f I/O files
C# 
C#
C#    Open binary file for storing Energies and coefficient data
C#######################################################################
      SUBROUTINE HR12FILE(N, L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 14 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE HR12FILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'R12 :: DATA FOR L = ', L

      ENDIF
         
      IF(L.EQ.0) OPEN(N,FILE='DAT/HR12-0.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.1) OPEN(N,FILE='DAT/HR12-1.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.2) OPEN(N,FILE='DAT/HR12-2.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.3) OPEN(N,FILE='DAT/HR12-3.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.4) OPEN(N,FILE='DAT/HR12-4.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.5) OPEN(N,FILE='DAT/HR12-5.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.6) OPEN(N,FILE='DAT/HR12-6.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.7) OPEN(N,FILE='DAT/HR12-7.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.8) OPEN(N,FILE='DAT/HR12-8.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.9) OPEN(N,FILE='DAT/HR12-9.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.10) OPEN(N,FILE='DAT/HR12-10.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.11) OPEN(N,FILE='DAT/HR12-11.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.12) OPEN(N,FILE='DAT/HR12-12.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.13) OPEN(N,FILE='DAT/HR12-13.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L.EQ.14) OPEN(N,FILE='DAT/HR12-14.DAT',
     1        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')

      RETURN
      END
C#######################################################################
      SUBROUTINE HCFGFILE(N, L)

      IF(L.GT.14) THEN

         WRITE(*,*) ' L MUST BE LESS THAN 15 . FOR L > 14', 
     1       ' MODIFY SUBROUTINE HCFGFILE. EXITING...'
         STOP

      ELSE

         WRITE(*,*) 'HCFG :: CONFIGURATIONS FOR L = ', L

      ENDIF
         
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
      END
C#######################################################################

C#EOF
