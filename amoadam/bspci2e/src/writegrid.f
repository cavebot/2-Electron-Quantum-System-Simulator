C#######################################################################
C#

      PROGRAM WRITEGRID

      IMPLICIT NONE
      

      INTEGER I, J, NWRITE, NREAD

C..... h1e.inp

      INTEGER L1E_H1E      

C.....  r12.inp
c$$$      INTEGER kb, nb
c$$$      REAL*8   alpha, beta, rmax
c$$$      REAL*8  rho(lf), rs(lf)
c$$$C..... grid.inp



C...  H2E.INP

      DOUBLE PRECISION E0_H2E, DE_H2E
      INTEGER  NE_H2E, NCORE_H2E, NOCNSTR_H2E, NCHL_H2E

C... KMTX.INP

      INTEGER NAFILE, NPP
      DOUBLE PRECISION Z_EFF_KMTX, ALPHA_KMTX, VP_KMTX

C...  HDMX2EBF.INP

      INTEGER LINITIAL_HDMX, NMININ_HDMX, NMAXI_HDMX
      INTEGER LFINAL_HDMX,   NMINFIN_HDMX,  NMAXF_HDMX
      CHARACTER*16 GAUGE
      INTEGER NH2E_HDMX, LTOTALMAX_HDMX
      INTEGER NHMX_HDMX, NCSI_HDMX, NCSF_HDMX, NSX_HDMX
C.... CS2PH.INP

      INTEGER LINITIAL_CS2PH, LFINAL_CS2PH

C....

      INTEGER idr, np, nrp
      REAL*8   xi, rmin, rmx
      

C....................................

c$$$      OPEN(9,FILE='inp/r12.inp') 
c$$$
c$$$      READ(9,*) alpha, beta
c$$$      READ(9,*) ( rho(i), i = 1, lf)
c$$$      READ(9,*)  kb, nb, rmax
c$$$      READ(9,*) ( rs(i), i = 1, lf)
c$$$
c$$$      CLOSE(9)
c$$$
c$$$      OPEN(9,FILE='INP/HR12.INP') 
c$$$
c$$$      WRITE(9, '(2E3.7)')      alpha, beta
c$$$      WRITE(9, '(10E3.7)') (     rho(i), i = 1, lf)
c$$$      WRITE(9, '(2I5,1X,E3,7)')  kb, nb, rmax
c$$$      WRITE(9,  '(10E3.7)' ) ( rs(i), i = 1, lf)
c$$$
c$$$      CLOSE(9)

C....................................


      NREAD  = 12
      NWRITE = 16


C....................................
C#
C#   read from single-channel code (FXD)
C#

      OPEN(NREAD,file='inp/grid.inp')

      READ(NREAD, *) idr
      READ(NREAD, *) xi, rmin, rmx
      READ(NREAD, *) np, nrp

      CLOSE(NREAD) 



C######################################################

      OPEN(NREAD, FILE = 'WF2E.INP') 

      READ(NREAD,*) L1E_H2E

C... h2e.inp      

      READ(NREAD,*) E0_H2E,    DE_H2E,      NE_H2E
      READ(NREAD,*) NCHL_H2E,  NOCNSTR_H2E  
      READ(NREAD,*) NCORE_H2E

C... kmtx.inp

      READ(NREAD,*) NAFILE_KMTX
      READ(NREAD,*) NPP_KMTX
      READ(NREAD,*) Z_EFF_KMTX, ALPHA_KMTX, VP_KMTX

C... hdmx2ebf.inp

      READ(UNIT = NREAD, * ) LINITIAL_HDMX, NMININ_HDMX,  NMAXI_HDMX
      READ(UNIT = NREAD, * ) LFINAL_HDMX,   NMINFIN_HDMX, NMAXF_HDMX
      READ(UNIT = NREAD, * ) GAUGE_HDMX
      READ(UNIT = NREAD, * ) NH2E_HDMX
      READ(UNIT = NREAD, * ) LTOTALMAX_HDMX
      READ(UNIT = NREAD, * ) NHMX_HDMX, NCSI_HDMX, NCSF_HDMX, NSX_HDMX

C... hdmx2eff.inp

      READ(UNIT = NREAD, * ) LINITIAL_HDMXF, NMININ_HDMXF,  NMAXI_HDMXF
      READ(UNIT = NREAD, * ) LFINAL_HDMXF,   NMINFIN_HDMXF, NMAXF_HDMXF
      READ(UNIT = NREAD, * ) NH2EI_HDMXF, NH2EF_HDMXF

C      READ(UNIT = NREAD, * ) LTOTALMAX_HDMX
C      READ(UNIT = NREAD, * ) NHMX_HDMX, NCSI_HDMX, NCSF_HDMX, NSX_HDMX



C...  cs2ph.inp

      READ(UNIT= NREAD,*) LINITIAL_CSH2PH, LFINAL_CS2PH

      CLOSE(NREAD)





C...

C################################################################

C#
C#   re-write again for multichannel code  (FREE)
C#

C..............................................

      OPEN(NWRITE,file='INP/GRID.INP')

      WRITE(NWRITE, '(I3)') idr
      WRITE(NWRITE, '(E10.2,1X,I5)'  ) xi,   nrp
      WRITE(NWRITE, '(2E10.2,1X,I5)' ) rmin, rmx, np

      CLOSE(NWRITE)
C.......................................

C..... H1E.INP

      OPEN(NWRITE, file='INP/H1E.INP')

      WRITE(NWRITE, * ) L1E_H2E
            
      CLOSE(NWRITE)

C.... HDMX1E.INP

      OPEN(NWRITE, file='INP/HDMX1E.INP')

      WRITE(NWRITE, * ) GAUGE_HDMX
            
      CLOSE(NWRITE)

C.... H2E.INP

      OPEN(NWRITE, file='INP/H2E.INP')
      
      WRITE(NWRITE,*) E0_H2E,    DE_H2E,      NE_H2E
      WRITE(NWRITE,*) NCHL_H2E,  NOCNSTR_H2E  
      WRITE(NWRITE,*) NCORE_H2E

      CLOSE(NWRITE)

!!!..............................


      OPEN(NWRITE, file='INP/KMTX.INP')

      WRITE(NWRITE,*) NAFILE_KMTX
      WRITE(NWRITE,*) NPP_KMTX
      WRITE(NWRITE,*) Z_EFF_KMTX, ALPHA_KMTX, VP_KMTX

      CLOSE(8)

C............. HDMX2EBF.INP


      OPEN(UNIT=NWRITE,file='INP/HDMX2EBF.INP')

      WRITE(UNIT = NWRITE, * ) LINITIAL_HDMX, NMININ_HDMX,  NMAXI_HDMX
      WRITE(UNIT = NWRITE, * ) LFINAL_HDMX,   NMINFIN_HDMX, NMAXF_HDMX
      WRITE(UNIT = NWRITE, * ) GAUGE_HDMX
      WRITE(UNIT = NWRITE, * ) NH2E_HDMX
      WRITE(UNIT = NWRITE, * ) LTOTALMAX_HDMX
      WRITE(UNIT = NWRITE, * ) NHMX_HDMX, NCSI_HDMX, NCSF_HDMX, NSX_HDMX

      CLOSE(NWRITE)

C............. HDMX2EFF.INP


      OPEN(UNIT=NWRITE,file='INP/HDMX2EBF.INP')

      WRITE(UNIT = NWRITE, * ) LINITIAL_HDMXF, NMININ_HDMXF, NMAXI_HDMXF
      WRITE(UNIT = NWRITE, * ) LFINAL_HDMXF,  NMINFIN_HDMXF, NMAXF_HDMXF
      WRITE(UNIT = NWRITE, * ) GAUGE_HDMX
      WRITE(UNIT = NWRITE, * ) NH2EI_HDMXF, NH2EF_HDMXF
      WRITE(UNIT = NWRITE, * ) LTOTALMAX_HDMX
      WRITE(UNIT = NWRITE, * ) NHMX_HDMX, NCSI_HDMX, NCSF_HDMX, NSX_HDMX

      CLOSE(NWRITE)


C............ CS1PH.INP 

      OPEN(UNIT=NWRITE,file='INP/CS1PH.INP',status='OLD')

      WRITE(UNIT=9,*) LINITIAL_HDMX
      WRITE(UNIT=9,*) GAUGE_HDMX
      
      END(NWRITE)


C............

      OPEN(UNIT=NWRITE,file='INP/CS1PH.INP',status='OLD')

      WRITE(UNIT=9,*) LINITIAL_NHDMX
      WRITE(UNIT=9,*) GAUGE_HDMX
      
      END(NWRITE)

C.............

      OPEN(UNIT=NWRITE,file='INP/CS2PH.INP',status='OLD')

      WRITE(UNIT=9,*) LINITIAL_CSH2PH, LFINAL_CS2PH
      WRITE(UNIT=9,*) GAUGE_HDMX
      
      END(NWRITE)

      return
      end

C################################################################



