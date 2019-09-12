! LOG NOTES
!
!   **  THIS PROGRAM WRITES THE INPUT FILES NEEDED FOR THE   **
!   **  BSPCI2EF PACKAGE IN A CONSISTENT WAY                 **
!
!
!      AUTHOR     :      L. A.A. NIKOLOPOULOS   
!      DATE       :      Sat Feb 16 20:04:48 EET 2002
!      INSTITUTE  :      IESL/FORTH 
!
!
!
!     16022002     :  REVISION O.O      LAANIESLFORTH
!                     
!            MODE  :  PREPARATION / NOT INCLUDED IN THE PACKAGE
!
!     INPUT  FILES :  WF2E.INP
!     OUTPUT FILES :  INP/H1E.INP       INP/HDMX1E.INP 
!                     INP/GRID.INP,
!                     INP/H2E.INP,      INP/KMTX.INP
!                     INP/HDMX2EBF.INP  INP/HDMX2EFF.INP
!                     INP/CS1PH.INP     INP/CS2PH.INP
!
!    NOTE          :    
!
!!!#####################################################################

      PROGRAM INPUT

      IMPLICIT NONE      
      INTEGER ZERO, ONE
      PARAMETER(ZERO=0, ONE=1)
      INTEGER I, J, NWRITE, NREAD

C.... GRID.INP

      INTEGER idr_grid, np_grid, nrp_grid
      REAL*8   xi_grid, rmin_grid, rmx_grid
      

C....................................

C..... h1e.inp

      INTEGER L1E_H1E, L2E_H2E      

C.....  r12.inp
c$$$      INTEGER kb, nb
c$$$      REAL*8   alpha, beta, rmax
c$$$      REAL*8  rho(lf), rs(lf)
c$$$C..... grid.inp


C...  H2E.INP

      DOUBLE PRECISION EINITIAL_H2E, DE_H2E
      INTEGER  NE_H2E, NCORE_H2E, NOCNSTR_H2E 

C... KMTX.INP

      INTEGER NAFILE_KMTX 
      DOUBLE PRECISION Z_EFF_KMTX, ALPHA_KMTX, VP_KMTX

C...  HDMX2EBF.INP

      INTEGER LINITIAL_HDMX, NMININ_HDMX, NMAXIN_HDMX
      INTEGER LINT_HDMX,   NMININT_HDMX,  NMAXINT_HDMX
      CHARACTER*16 GAUGE_HDMX
      INTEGER NH2E_HDMX, LTOTALMAX_HDMX

C...  HDMX2EFF.INP

      INTEGER LFINAL_HDMXf,   NMINFIN_HDMXF,  NMAXFIN_HDMXF
      INTEGER NH2EI_HDMXF, NH2EF_HDMXF

C.... CS2PH.INP

      INTEGER LINITIAL_CS2PH, LFINAL_CS2PH

C....    start the body


      NREAD     = 12
      NWRITE    = 16
      NCORE_H2E = 5

C#
C#   read from single-channel code (FXD)
C#

      OPEN(NREAD,file='inp/grid.inp')

      READ(NREAD, *) idr_grid
      READ(NREAD, *) xi_grid, rmin_grid, rmx_grid
      READ(NREAD, *) np_grid, nrp_grid

      CLOSE(NREAD) 

C######################################################

      OPEN(NREAD, FILE = 'INP/HE.INP') 

      READ(NREAD,*) L1E_H1E, L2E_H2E

C... h2e.inp      

      READ(NREAD,*) EINITIAL_H2E, DE_H2E, NE_H2E
      READ(NREAD,*) NOCNSTR_H2E  

C... kmtx.inp

      READ(NREAD,*) NAFILE_KMTX
      READ(NREAD,*) Z_EFF_KMTX, ALPHA_KMTX, VP_KMTX

C... hdmx2ebf.inp

      READ(NREAD, * ) GAUGE_HDMX
      READ(NREAD, * ) NH2E_HDMX
      READ(NREAD, * ) LINITIAL_HDMX, NMININ_HDMX,  NMAXIN_HDMX
      READ(NREAD, * ) LINT_HDMX, NMININT_HDMX, NMAXINT_HDMX
C... hdmx2eff.inp

      READ(NREAD, * ) LFINAL_HDMXF, NMINFIN_HDMXF,NMAXFIN_HDMXF
      READ(NREAD, * ) NH2EI_HDMXF,  NH2EF_HDMXF

C...  cs2ph.inp
C      READ(UNIT= NREAD,*) LINITIAL_CS2PH, LFINAL_CS2PH

      CLOSE(NREAD)
C...

C################################################################

C#
C#   re-write again for multichannel code  (FREE)
C#

C..............................................

      OPEN(NWRITE,file='INP/GRID.INP')
      WRITE(NWRITE, '(I4)') idr_grid
      WRITE(NWRITE, '(E10.2,1X,I5)'  ) xi_grid,   nrp_grid
      WRITE(NWRITE, '(2E10.2,1X,I5)' ) rmin_grid, rmx_grid, np_grid
      CLOSE(NWRITE)

C..... H2E.INP
C      OPEN(NWRITE, file='INP/WF2E.INP')
C      WRITE(NWRITE,  '(I3)' ) L2E_H2E
C      CLOSE(NWRITE)

C..... H1E.INP
C      OPEN(NWRITE, file='INP/H1E.INP')
C      WRITE(NWRITE,  '(I3)' ) L1E_H1E            
C      CLOSE(NWRITE)

C.... H2E.INP

      OPEN(NWRITE, file='INP/H2E.INP')
      
      WRITE(NWRITE,'(2E15.4,1X,I5)') EINITIAL_H2E, DE_H2E, NE_H2E
      WRITE(NWRITE, '(I4)')          NOCNSTR_H2E  
      WRITE(NWRITE, '(I4)')          NCORE_H2E
      CLOSE(NWRITE)

C.... KMTX.INP
C      OPEN(NWRITE, file='INP/K2E.INP')
C      WRITE(NWRITE, '(I4)') NAFILE_KMTX
C      WRITE(NWRITE, '(I4)') NP_GRID
C      WRITE(NWRITE,  '(3E10.2)') Z_EFF_KMTX, ALPHA_KMTX, VP_KMTX
C      CLOSE(8)

C.... HDMX1E.INP
C      OPEN(NWRITE, file='INP/HDMX1E.INP')
C      IF(GAUGE_HDMX.EQ.'V') THEN
C         WRITE(NWRITE, '(3X,I4)' ) ZERO
C      ELSE
C         WRITE(NWRITE, '(3X,I4)' ) ONE
C      ENDIF
C      CLOSE(NWRITE)

C..............


      IF(NMAXIN_HDMX.GT.NE_H2E)            NMAXIN_HDMX = NE_H2E
      IF(NMAXINT_HDMX.GT.NE_H2E)          NMAXINT_HDMX = NE_H2E
      IF(NMAXINT_HDMX.GT.NE_H2E)         NMAXFIN_HDMXF = NE_H2E




C............. HDMX2EBF.INP      | INITIAL > ---> | INTERMEDIATE  >


      OPEN(UNIT=NWRITE,file='INP/HDMX2EBF.INP')

      WRITE(NWRITE,'(3I4)')    NMININ_HDMX,  NMAXIN_HDMX
      WRITE(NWRITE,'(3I4)')   NMININT_HDMX
      CLOSE(NWRITE)
C      WRITE(NWRITE,'(3I4)')   LINITIAL_HDMX, NMININ_HDMX,  NMAXIN_HDMX
C      WRITE(NWRITE,'(3I4)')   LINT_HDMX,     NMININT_HDMX
C      WRITE(NWRITE,'(3I4)')   LINT_HDMX,     NMININT_HDMX, NMAXINT_HDMX
C      WRITE(NWRITE,'(3X,A4)') GAUGE_HDMX
C      WRITE(NWRITE,'(I4)' )   NH2E_HDMX
C      WRITE(NWRITE,'(I4)' )   L1E_H1E - 1


C............. HDMX2EFF.INP :  | INTERMEDIATE  > ---> | FINAL >


      OPEN(UNIT=NWRITE,file='INP/HDMX2EFF.INP')

C      WRITE(NWRITE,'(3I4)')   LINT_HDMX,    NMININT_HDMX,  NMAXINT_HDMX
      WRITE(NWRITE,'(3I4)')   LINITIAL_HDMX, NMININ_HDMX,  NMAXIN_HDMX
      WRITE(NWRITE,'(3I4)')   LFINAL_HDMXF, NMINFIN_HDMXF, NMAXFIN_HDMXF
      WRITE(NWRITE,'(3X,A4)') GAUGE_HDMX
      WRITE(NWRITE,'(I4)' )   NH2EI_HDMXF, NH2EF_HDMXF
      WRITE(NWRITE,'(I4)' )   L1E_H1E - 1

      CLOSE(NWRITE)


C............ CS1PH.INP 


      OPEN(UNIT=NWRITE,file='INP/CS1PH.INP')

      WRITE(NWRITE, '(I4)' )     LINITIAL_HDMX 
      WRITE(NWRITE,  '(3X,A4)' ) GAUGE_HDMX

      CLOSE(NWRITE)

C............

      OPEN(UNIT=NWRITE,file='INP/CS2PH.INP')

      WRITE(NWRITE, '(2I4)' )   LINITIAL_HDMX, LFINAL_HDMXF
      WRITE(NWRITE,  '(3X,A4)' ) GAUGE_HDMX
      
      CLOSE(NWRITE)

     
      END

C################################################################



