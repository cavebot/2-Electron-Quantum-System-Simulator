C#####################################################################
C#
C#    written  by L. A.A. Nikolopoulos 2001
C#
C#        10/2001 : modified to accomodate for negative ion 
C#                  systems ( i.e. H-) 
C#                  modified to accomodate for 1/r^2 as well as
C#                  1/r^4 interactions
C#
C#    Calculates the phase and normalization constant for the 
C#    continuum states of a 2-e atom. 
C#    Normalization is done according the Burgess article (1963)  
C#    Outpout file used as input files in the calculation of 
C#    2-e DME.
C# 
C#    
C#    A. Burgess, PROC. PHYS. SOC. 1963, VOL. 81, pg. 442
C#   ' The determination of Phases and Amplitudes of Wave Functions'
C#
C#    inputs :
C#             L     :  2-e angular symmetry  
C#             Ethr  :  ionization threshold    
C#             kll   :  channel number
C#             zeff  :  charge seen in large distance
C#             apol :  short-range (CI interaction) potential ( a/r^2 )
C#             vpol :  short-range (polarization) potential   ( vp/r^4) 
C#
C#        l.aa.n   2001/iesl
C#
C######################################################################
      PROGRAM NORM2E
C     
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.2e.inc"
      INCLUDE "units.inc"
C...........

      integer  nwf2efile, ninp, n2e_bin, n2e_asc
      CHARACTER*100 ARGV

C..........

C     
      DOUBLE PRECISION E_ION
      DIMENSION CIIN(NHX,NWF2E),ENI(NWF2E)
      DIMENSION NHFI(NCS),LHFI(NCS),LLI(NCS)
      DIMENSION NMINI(NCS),NMXI(NCS)
      DIMENSION NDI(NCS),DR(NP),P(NP)
      DIMENSION R(NP), YK(NP)
      COMMON/PCA/PC(NS,NP)

C.........................................


      nwf2efile  = 1
      n2e_bin    = nwf2efile + 1
      n2e_asc    = nwf2efile + 2
      ninp       = nwf2efile + 3
C........................................

C#
C#  read program's start inputs 
C#


C GET THE PARAMETERS FROM THE COMMAND LINE ARGUMENT



      CALL GETARG(1, ARGV)
      READ(ARGV,*) L
 
      OPEN(NINP,file='inp/n2e.inp')
      READ(NINP,*)  KLL, KLL_MAX 
      READ(NINP,*)  ZEFF
      READ(NINP,*)  APOL, VPOL
      CLOSE(NINP)

      print*, "#      kll = ", kll, kll_max
      print*, "#     zeff = "

C#
C#  get the ionization threshold of the ion he+,mg+,ca+,li2+
C#


      CALL READ_ION_THRESHOLD(E_ION)   ! return e_ionic in a.u.

C  WRITE IN RYD

      E_ION = 2*E_ION

C      E_IONIC = -1.1025

C#
C#    read data from the 2e-angular symmetry.
C#

      LOI = L

      CALL hfile(nwf2efile,"dat","w2eb","bin",loi)              !initial states

      READ(unit=NWF2EFILE)    LOI, LSI
      READ(unit=NWF2EFILE)    NCMXI
      READ(unit=NWF2EFILE)  ( NHFI(K), K = 1, NCMXI)
      READ(unit=NWF2EFILE)  ( LHFI(K), K = 1, NCMXI)
      READ(unit=NWF2EFILE)  ( LLI(K),  K = 1, NCMXI)
      READ(unit=NWF2EFILE)  ( NMINI(K),K = 1, NCMXI)
      READ(unit=NWF2EFILE)  ( NMXI(K), K = 1, NCMXI)
      READ(unit=NWF2EFILE)  ( NDI(K),  K = 1, NCMXI)
      READ(unit=NWF2EFILE)    NHMXI
      READ(unit=NWF2EFILE)    NSI

      FLI = DBLE( 2*LOI + 1 )

      DO NSC = 1,  NSI

         READ(NWF2EFILE) ENI(NSC)
         READ(NWF2EFILE) (CIIN(K,NSC), K = 1, NHMXI)

      ENDDO

      CLOSE(NWF2EFILE)


C#  make the grid 

      CALL RIN(r, dr, h, no, idr)

C      print*, h, no, idr
C      print*,r

      print*, "#"
      print*, "# CONFIGURATION SERIES   K = ", KLL
      print*, "#"
      print*, "# n1,l1,l2, n2min,n2max = (" 
      print*, NHFI(KLL),LHFI(KLL) ,LLI(KLL), NMINI(KLL), NMXI(KLL), ")"
      print*, "#"

C# get data for the channel wfs | n1,l1,l2; n2 >, n2 =nmini(kll)-nmxi(kll)


      call wfnlin( yk, lli(kll), no, nmini(kll), nmxi(kll) )


C# total number of channel functions

      ncs_1e = 0
      DO  iu = 1, kll - 1
         DO  ik = 1, nmxi(iu) - nmini(iu) + 1

          ncs_1e = ncs_1e + 1

         ENDDO
      ENDDO

C.................

      CALL hfile(n2e_bin,  "dat","n2e" ,"bin"  , loi)        !initial states
      CALL hfile(n2e_asc,  "out","n2e" ,"ascii", loi)        !initial states
      CALL hfile(nwf2efile,"out","w1ek","ascii", kll_out)    !partial channel

C      CALL NORMFILE(NNORM, LOI) 


!      OPEN(NNORMASCII, FILE='out/norm.out') 


      WRITE(n2e_bin) NSI, LLI(KLL)


      IF(KLL_MAX.EQ.0)  KLL_MAX = NSI


C# loop over all states 

      DO  KL = 1, KLL_MAX

        print*, "#" 
        print*, "#             state       kl = ", kl, eni(kl)

         DO j = 1, no

            p(j) = 0.0D+00
            kzz2 = ncs_1e

            DO kzz = 1, nmxi(kll) - nmini(kll) + 1

               kzz2 = kzz2 + 1
               p(j) = p(j) + pc(kzz2, j) * ciin(kzz2, kl)
            ENDDO

            if(kl.eq.kl_out)  then 
               write(nwf2efile,*) r(j), p(j)
            endif
         ENDDO


C#     calculate normalization factor and phase shift for state KL


      CALL nom(n2e_bin,kl,eni(kl),lli(kll),no,r,p,zeff,e_ion,apol,vpol)
               
      ENDDO

      close(n2e_bin)
      close(n2e_asc) 
      close(nwf2efile)

      WRITE(*,*)'# norm2e::              NWF2E = ', NSI
      WRITE(*,*)'# norm2e::              L_2E  = ', L
      WRITE(*,*)'# norm2e::                l_2 = ', LLI(KLL)
C.........................................

      end
C#######################################################################
C# GET THE ENERGY OF THE FIRST IONIZATION THRESHOLD
C# THE FOLLOWING APPLIES ONLY TO NS^2 ATOMS (Hn,He,Mg,Ca)
C# 

      SUBROUTINE READ_ION_THRESHOLD(E_ION)
!
      IMPLICIT NONE
!
      DOUBLE PRECISION ENAU
      PARAMETER ( ENAU     = 27.211 396 181  D+00 )
      INTEGER NINP,I
      INTEGER LMIN_1E, LMAX_1E
      INTEGER NOF_STATES, NVALENCE
      CHARACTER*20 EN1EFILE
      DOUBLE PRECISION E_ION
C..............................
      NINP = 9 
C................
      EN1EFILE ='dat/en1e.dat'

      OPEN(NINP, FILE=EN1EFILE)

      READ(NINP,*) LMIN_1E, LMAX_1E

      IF(LMIN_1E.NE.1) THEN         
         WRITE(*,*) '# NORM1E::  ION ENERGIES NOT FOR L = 0. EXITING..'
         STOP
      ENDIF
 
      READ(NINP,*) NOF_STATES, NVALENCE

      DO  I = 1, NVALENCE
         READ(NINP,*)  E_ION
      ENDDO

      CLOSE(NINP)

      E_ION = 0.5D+00 * E_ION

      WRITE(*,*)'# read_ion_threshold:: IONIZATION THRESHOLD FOR ION'
      WRITE(*,*)'# read_ion_threshold::    E_ION =', E_ION,     " AU "
      WRITE(*,*)'# read_ion_threshold::    E_ION =', E_ION*ENAU," eV "

      RETURN
      END
C####################################################################
C*EOF









