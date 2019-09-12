C*
C***********************************************************************
C*
      PROGRAM tdse
C*
C***********************************************************************
C*
C*    This program performs the time propagation describing a
C*    two-electron atom in a time-dependent laser field expanding
C*    the time-dependent wavefunction in a basis of time-independent
C*    field-free atomic two-electron states.
C*
C*    At the moment, three different pulse shapes are implemeted, 
C*    a Gaussian and a cos^2 type pulse with respect to the vector
C*    potential A(t) (both pulse forms are e.g. given in eqs.(66) and 
C*    (67) of Phys.Rep. vol.305, 203 (1998)) and a cos^2 pulse with
C*    respect to the field E(t). 
C*
C*    In the case of the cos^2 pulse the integration is performed 
C*    from -0.5\tau to +0.5\tau. In the case of a Gaussian pulse the
C*    integration limit is taken from the input file (and should 
C*    usually be at least 3.0).
C*
C*    It is possible to specify (via the input file) either the 
C*    number of optical cycles or the pulse duration (the other 
C*    quantity is then obtained from the specified one). The other
C*    pulse specifications are the photon energy E_phot = \hbar*omega 
C*    and the peak intensity I. 
C*
C*    This program uses the NAG routine D02BAF (old version) or the
C*    IMSL routine DIVPRK (new version) to integrate the differential
C*    equation.
C*
C*    This program can be used also for restarts based on an existing
C*    data file (unit 33) containing time and vector information.
C*
C*    Required routines: GETSUB.F, NAG or IMSL library, and ESSL library
C*    ==================
C*
C*----------------------------------------------------------------------
C*
C*    Author: Xian Tang
C*    Date  : January 1990
C*    Modified by:  H. Rudolph (?), Jian Zhang, Paul Maragakis, and
C*         Alejandro Saenz (since February 1999)
C*    Last Change:  29.06.2000
C*
C*    29.06.2000: New determination of variable LSTART in the case of 
C*       a restart is implemented. (The original one seemed not to 
C*       work.) 
C*    10.09.1999: The output into file no.34 is changed. Now this file
C*       is unformatted and contains the complete final wave function
C*       at the end of the pulse. (Before, it was formatted and contained
C*       the populations of the five lowest states for the first ten
C*       angular momenta at the end of every interval.) 
C*    07.09.1999: Implement the choice in between three pulse shapes, 
C*       and the option to use either the number of optical cycles or
C*       the total pulse duration for specifying the pulse length.
C*    30.08.1999: Change of the pulse (in accordance to the code of 
C*       Lampros it has now cos^2 behaviour of E(t) and not of A(t)).  
C*       In addition it is made sure that the pulse starts and 
C*       ends with zero field. Also some constants have been defined
C*       with higher accuracy. The pulse is now calculated and 
C*       stored in a separate file. Finally, the array dimensioning 
C*       is now done via preprogram TDSDIM.F that generates the 
C*       parameter file 'param_tdse.f'.
C*    29.03.1999: A complete restructuring of the code. The input/output
C*       structure is changed. The code is now commented.
C*    31.03.1990: Changes made in order to remove the rapidly
C*       oscillating part of the matrix. Pulse shape changed to Gaussian
C*       form. Pulse duration changed from -1.0 times \tau to 1.0 times
C*       \tau.
C*
C*======================================================================
C*
      IMPLICIT REAL*8 (a-h,o-z)
C*
C*  
C*         Parameter:   
C*         ==========
C*
C*             lvlmax: max. number of L values.
C*             nvlmax: max. number of eigenstates per L value.
C*             ncfmax (=lvlmax*nvlmax): max. number of coefficients,
C*                     i.e. max. number of field-free states.
C*             nstmax: max. number of time steps (considering also 
C*                     sub-intervals).
C*
      INCLUDE 'param_tdse.f' 
C*
      PARAMETER ( ncfmax = lvlmax * nvlmax )
C*
C*IMSL  1
      PARAMETER ( nbbmax = 20 * ncfmax ) 
C*
C*
C*NAG   1
C      DIMENSION y(ncfmax), yp(nvmax), w(nvmax,7)
C*
C*IMSL  1
      DIMENSION y(2*ncfmax), yp(2*ncfmax), parp(50)
C*
C*
      DIMENSION dpr(nvlmax), dpi(nvlmax), time(nstmax)
      DIMENSION dipole(nstmax)
C*
      CHARACTER*3 ustrng,vstrng
      CHARACTER*10 rstrng
      CHARACTER*132 filnam
      CHARACTER*132 filinp(lvlmax)
C*
      LOGICAL exists,lcut
C*
      EXTERNAL fcn,fcn2
C*
      COMMON/f/am(nvlmax,nvlmax,lvlmax-1), en(nvlmax,lvlmax), n(lvlmax), 
     +     nsum(lvlmax), n1, ntot, fi, omag, ta, neval
      COMMON/f1/dag(ncfmax)
      COMMON/pival/pi
      COMMON/worksp/rwksp(nbbmax)
      COMMON/ivar/ipulse
      COMMON/dvar/dlimit
C*
C*----------------------------------------------------------------------
C*
C*IMSL  1
      CALL iwkin(nbbmax)
C***
C*
      pi = ACOS( - 1.0D+00 )
      nvl = nvlmax
C*
C*
C*         Open the log file:
C*
      CALL getenv('LOGFIL',filnam)
      OPEN ( UNIT=10, FILE=filnam, STATUS='UNKNOWN',
     +       FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C*
C*
C*         Write program header into the log file:
C*
      WRITE(10,1002)
      WRITE(10,1004)
      WRITE(10,1005)
      WRITE(10,1006)
      WRITE(10,1005)
      WRITE(10,1004)
      WRITE(10,1002)
      WRITE(*,1002)
      WRITE(*,1004)
      WRITE(*,1005)
      WRITE(*,1006)
      WRITE(*,1005)
      WRITE(*,1004)
      WRITE(*,1002)
C*
C*
C*         Open and read the input file (previously 'ntestCI.in'):
C*
      CALL getenv('INPFIL',filnam)
      OPEN ( UNIT=18, FILE=filnam, STATUS='OLD',
     +       FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C*
C*
C*         Read the peak intensity in W/cm^2, transform it into 
C*         atomic units, and calculate the corresponding value of 
C*         E_0:
C*
      CALL getd(18,ff)
      WRITE(10,1060) ff
      WRITE(*,1060) ff
      ff = ff / 3.509338D+16
      WRITE(10,1062) ff
      WRITE(*,1062) ff
      fi = SQRT(ff)
C*
C*
C*         Read the photon energy in eV, transform it into atomic 
C*         units where it is identical to the frequency omega:
C*
      CALL getd(18,omag)      
      WRITE(10,1070) omag
      WRITE(*,1070) omag
      omag = omag / 27.211396181D+00
      WRITE(10,1072) omag
      WRITE(*,1072) omag
C*
      topt = 2.0D+00 * pi / omag
C*
C*
C*         Read the pulse type:
C*
      CALL geti(18,ipulse)
      IF ((ipulse. LT. 1) .OR. (ipulse .GT.3)) THEN
        WRITE(*,1035)
        WRITE(*,1038) ipulse
        WRITE(10,1035)
        WRITE(10,1038) ipulse
        CLOSE(18)
        CLOSE(10)
        STOP
      ENDIF
C*
      IF (ipulse .EQ. 1) THEN
        WRITE(*,1260) 
        WRITE(10,1260) 
      ENDIF
C*
      IF (ipulse .EQ. 2) THEN
        WRITE(*,1270) 
        WRITE(10,1270) 
      ENDIF
C*
      IF (ipulse .EQ. 3) THEN
        WRITE(*,1280) 
        WRITE(10,1280) 
      ENDIF
C*
C*
C*         Read the number of cycles and pulse duration tau:
C*
      CALL get2d(18,cycle,tau)
      IF ((cycle .LE. 0.0D+00) .AND. (tau .GT. 0.0D+00)) THEN
        tau = tau * 41.34127637D+00
        cycle = tau / topt
      ELSE
        IF ((cycle .GT. 0.0D+00) .AND. (tau .LE. 0.0D+00)) THEN
          tau = cycle * topt
        ELSE
          WRITE(*,1035)
          WRITE(*,1300)
          WRITE(*,1310)
          WRITE(*,1320) cycles,tau
          WRITE(10,1035)
          WRITE(10,1300)
          WRITE(10,1310)
          WRITE(10,1320) cycles,tau
          CLOSE(18)
          CLOSE(10)
          STOP
        ENDIF
      ENDIF
C*
      WRITE(10,1080) tau / 41.34127637D+00
      WRITE(*,1080) tau / 41.34127637D+00
      WRITE(10,1081) tau
      WRITE(*,1081) tau
      WRITE(10,1105) cycle
      WRITE(*,1105) cycle
C*
      ta = 0.5D+00 * tau
      tfac = ta / topt
C*
C*
C*         Read the (lower or upper) integration limit:
C*
      CALL getd(18,dlimit)
      IF (ipulse .EQ. 1) dlimit = 0.5D+00
      IF (ipulse .EQ. 2) dlimit = 0.5D+00
C*
C*
C*         Read the number of intervals (INTERV) in which the propagation 
C*         should be split: 
C*
      CALL geti(18,interv)
      WRITE(10,1075) interv
      WRITE(*,1075) interv
C*
C*
C*         Read the number of subintervals (LKMAX) per interval for 
C*         which results should be saved for calculating the higher-
C*         harmonics spectra:
C*
      CALL geti(18,lkmax)
      WRITE(10,1077) lkmax
      WRITE(10,1078)
      WRITE(*,1077) lkmax
      WRITE(*,1078)
C*
C*
C*         Read the tolerance parameter (for the numerical solution
C*         of the differential equation):
C*
      CALL getd(18,tol)
      WRITE(10,1150) tol
      WRITE(*,1150) tol
C*
C*
C*         The integration range is 2 times the integration limit
C*         (-limit ... +limit). A second factor two comes from the
C*         fact that the integration range is expressed in units
C*         of 0.5*\tau (and not \tau): 
C*
      range = 4.0D+00 * dlimit
      step = range / DFLOAT( interv * lkmax )
C*
C* 
C*         Read the maximum value of L:
C*
      CALL geti(18,llmax)
      IF (llmax .LE. 0) THEN
        WRITE(10,1035)
        WRITE(10,1082) llmax
        WRITE(*,1035)
        WRITE(*,1082) llmax
        GOTO 990
      ENDIF
      WRITE(10,1085) llmax
      WRITE(*,1085) llmax
C*
C*
C*         Read the threshold energy in atomic units Hartree:
C*
      CALL getd(18,er)
      IF ( er .GT. 0.0D+00 ) er = - er
      WRITE(10,1090) er
      WRITE(*,1090) er
C*
C*
C*         Read the cut-off energy (i.e. the highest energy to be 
C*         included) in atomic units Hartree:
C*
      CALL getd(18,ecut)
      IF (ecut .GT. er) THEN
        WRITE(10,1095) ecut
        WRITE(*,1095) ecut
      ELSE
        WRITE(10,1092)
        WRITE(*,1092)
      ENDIF
C*
C*
C*         Read the parameter deciding whether the names of the input 
C*         files should be generated automatically (only possible if 
C*         the configurations used for all values of L stem from the 
C*         same xxx.cfs file) or if they should be read explicitly:
C*  
      CALL geti(18,ifltyp)
C*
      IF (ifltyp.EQ.1) THEN
C*
C*
C*            Read the filenames:
C*
        WRITE(10,1210)
        WRITE(*,1210)
        DO 40 i=0,llmax-1
          CALL getstr(18,filnam)
          filinp(i+1) = filnam
   40   CONTINUE
      ELSE
C*
C*
C*            Generate the file names:
C*
        WRITE(10,1220)
        WRITE(*,1220)
        ustrng = 'UUU' 
        vstrng = 'VVV' 
        CALL getenv('BASNAM',filnam)
        ilianf = INDEX(filnam,ustrng)
        ilfanf = INDEX(filnam,vstrng)
        iliend = ilianf + 2
        ilfend = ilfanf + 2
      ENDIF
C*
      CLOSE(18)
C*
      irntyp = 0
      CALL getenv('RUNTYP',rstrng)
      IF (rstrng .EQ. 'RESTART   ') THEN 
        irntyp = 1
        WRITE(10,1200)
        WRITE(*,1200)
      ELSE
        IF (rstrng .EQ. 'POPULATION') THEN 
          irntyp = -1
          WRITE(10,1205)
          WRITE(*,1205)
        ENDIF
      ENDIF
C*
      IF (rstrng .EQ. 'PULSE     ') THEN
        t = 0.0D+00
        GOTO 245
      ENDIF
C*
C*
C*        Check whether the data file exists, if this is a restart:
C*
      CALL getenv('DATFIL',filnam)
      INQUIRE ( FILE=filnam, EXIST=exists )
C*
      IF ( (IABS(irntyp) .EQ. 1) .AND. .NOT. exists ) THEN
        WRITE(10,1035)
        WRITE(10,1050)
        WRITE(10,1055)
        WRITE(*,1035)
        WRITE(*,1050)
        WRITE(*,1055)
        GOTO 990
      ENDIF
C*
C*
C*        Check to make sure that a data file does not exist,
C*        if this is not a restart:
C*
      IF ( (IABS(irntyp) .NE. 1) .AND. exists ) THEN
        WRITE(10,1035)
        WRITE(10,1040)
        WRITE(10,1045)
        WRITE(*,1035)
        WRITE(*,1040)
        WRITE(*,1045)
        GOTO 990
      ENDIF
C*
C*
C*        Open the unformatted data file (unit 33):
C*
      CALL getenv('DATFIL',filnam)
      OPEN ( UNIT=33, FILE=filnam, STATUS='UNKNOWN',
     +       FORM='UNFORMATTED', ACCESS='SEQUENTIAL' )
C*
C*
C*        Open and read the unformatted files containing the 
C*        energies and transition dipole matrix elements between
C*        two respective values of the two-electron angular-momentum
C*        quantum number L:
C*        ( The read information comprises:
C*              lb: L_{initial},  lz: L_{final}
C*              n(i)/n(i+1): number of initial/final states
C*              en(j,i): energy of state no.:j with L=i-1
C*              am(k,l,i): dipole transition matrix element 
C*                     between state no.:k (L=i) and state no.:l
C*                     (L=i-1). )
C*
      n1 = llmax + 1
      DO 140 i=1,n1-1
         IF (ifltyp .EQ. 1) THEN
           filnam = filinp(i)
         ELSE
           CALL getenv('BASNAM',filnam)
           WRITE(filnam(ilianf:iliend),'(I3.3)') i-1
           WRITE(filnam(ilfanf:ilfend),'(I3.3)') i
         ENDIF
         OPEN ( UNIT=17, FILE=filnam, STATUS='OLD',
     +          ACCESS='SEQUENTIAL', FORM='UNFORMATTED' )
         READ(17) lb,lz,n(i),n(i+1)
         READ(17) ( en(k,lb+1), k=1,n(i) ) 
         READ(17) ( en(k,lz+1), k=1,n(i+1) )
         DO 130 nr=1,n(i)
            READ(17) (am(nc,nr,i), nc=1,n(i+1))
  130    CONTINUE
         CLOSE(17)
  140 CONTINUE
C*
      IF (IABS(irntyp) .EQ. 1) THEN
         READ(33) ff,omag,ta,llmin,llmax,ecut
      ELSE
         WRITE(33) ff,omag,ta,llmin,llmax,ecut
      ENDIF
C*
C*
C*         Transform the state energies (in Ryd) into photoelectron
C*         energies (in E_h), and select only those states with energy
C*         below E_cut (=ecut):
C*
      ntot = 0
      IF (ecut .GT. er) THEN
        lcut = .TRUE.
      ELSE
        lcut = .FALSE.
      ENDIF
      DO 160 k=1,n1
         iu = 0
         DO 150 k1=1,n(k)
            en(k1,k) = ( en(k1,k) / 2.0D+00 ) - er
            IF ((lcut .EQV. .TRUE.).AND.(en(k1,k) .GT. ecut)) 
     +                                                   GOTO 150
            iu = iu + 1
  150    CONTINUE
         n(k) = iu
         IF (IABS(irntyp) .EQ. 1) THEN
            READ(33) n(k),(en(k1,k),k1=1,n(k))
         ELSE
            WRITE(33) n(k),(en(k1,k),k1=1,n(k))
         ENDIF
         ntot = ntot + n(k)
         WRITE(10,1230) k-1,n(k)
         WRITE(*,1230) k-1,n(k)
  160 CONTINUE
C*
      IF (IABS(irntyp) .EQ. 1) THEN
         READ(33) ntot
      ELSE
         WRITE(33) ntot
      ENDIF
      WRITE(10,1240) ntot
      WRITE(*,1240) ntot
C*
C*
C*         Restore the energies in the one-dimensional array DAG:
C* 
      nb1 = 0
      DO 180 k=1,n1
         DO 170 k1=1,n(k)
            nb1 = nb1 + 1
            dag(nb1) = en(k1,k)
  170    CONTINUE
  180 CONTINUE
C*
      lstart = 1
C*
C*
C*         Set up the coefficients and the time variable in the case
C*         of a restart:
C*
      IF (IABS(irntyp) .EQ. 1) THEN
         t = - dlimit * 2.0D+00 * ta
  190    CONTINUE
            READ(33,END=210) t1,(yp(ii),ii=1,2*ntot)
            lstart = lstart + 1
            t = t1
            t2 = t1
            DO 200 ik=1,2*ntot
               y(ik) = yp(ik)
  200       CONTINUE
         GOTO 190
  210    CONTINUE
         t = t + dlimit * 2.0D+00 * ta
         ktanf = ( lstart - 1 ) * lkmax
      ENDIF
C*
C*
C*        Set up the coefficients and time variable for a new 
C*        (i.e. a non-restarted) calculation by setting all
C*        coefficients except the first one to zero (and the first
C*        one to one) and setting t=0, as it corresponds to the 
C*        situation that the system starts at t=0 in its ground 
C*        state:
C*
      IF ( lstart .EQ. 1 ) THEN
         DO 230 i=1,nvmax
            y(i) = 0.0D+00
  230    CONTINUE
         y(1) = 1.0D+00
         t = 0.0D+00
         ktanf = 0
      ENDIF
C*
      IF ( irntyp .EQ. -1 ) GOTO 420
C*
      ido = 1
C***
C      parp(1) = 0.1D+00
C      parp(2) = 0.01D+00
C***
      parp(4) = 1000000
      nsum(1) = 1
      DO 240 i=2,llmax+1
         nsum(i) = nsum(i-1) + n(i-1)
  240 CONTINUE
      kt = ktanf
      prop = 0.0D+00
      yield = 0.0D+00
C*
C*
C*        Calculate and store the pulse:
C*
  245 CONTINUE
      tend = t + ( ta * range )
      CALL pulse(t,tend,fi,omag,ta,dlimit,ipulse)
      IF (rstrng .EQ. 'PULSE     ') GOTO 990
C*
C*
C*        Split the complete calculation into INTERV intervals 
C*        so that restarts become possible, since after each 
C*        interval the intermediate result, i.e. the wavefunction
C*        at a given time is stored in the unformatted data file 33:
C*
      DO 400 l=lstart,interv
C*
C*
C*          Break into more sub-intervals for calculating the time-
C*          dependent dipole (but without storing the intermediate 
C*          result to allow restarts):
C*
         DO 260 lk=1,lkmax
            kt = kt + 1
            tend = t + ( ta * step )
            neval = 0
C*
            WRITE(*,1160) kt,t,tend
            WRITE(10,1160) kt,t,tend
C*
C*
C*             Integrate the differential equation for given initial
C*             values from time T to time TEND:
C*
C*NAG   1
C            CALL D02BAF(t,tend,2*ntot,y,tol,fcn,w,ifail)
C***
C*
C*IMSL  1
            CALL DIVPRK(ido,2*ntot,fcn,t,tend,tol,parp,y)
C***
C*
C*             Calculate the time-dependent dipole moment:
C*
            time(kt) = tend
C*
C* li = 0, lmax

            DO 250 li=1,n1-1   

               CALL DGEMV('n',n(li+1),n(li),1.0D+00,am(1,1,li),nvl,
     +                    y(nsum(li)),1,0.0D+00,dpr,1)
               CALL DGEMV('n',n(li+1),n(li),1.0D+00,am(1,1,li),nvl,
     +                    y(ntot+nsum(li)),1,0.0D+00,dpi,1)
               dipole(kt)=
     +              - DDOT(n(li+1),y(ntot+nsum(li+1)),1,dpr,1)
     +              + DDOT(n(li+1),y(nsum(li+1)),1,dpi,1)
  250       CONTINUE
C*
C*       End of the loop over subintervals lk:
  260    CONTINUE
C*
         WRITE(10,1250) neval
         WRITE(*,1250) neval
C*
C*
C*          Store the intermediate (or final) result in the 
C*          unformatted data file 33:
C*
         IF ( IABS(irntyp) .EQ. 1 ) THEN
            REWIND (33)
            WRITE(33) ff,omag,ta,llmin,llmax,ecut
            DO 270 kp=1,n1
               WRITE(33) n(kp),(en(k1,kp),k1=1,n(kp))
  270       CONTINUE
            WRITE(33) ntot
            irntyp = 0
         ENDIF
         WRITE(33) tend - dlimit * 2.0D+00 * ta, (y(ii),ii=1,2*ntot)
C*
C*
C*          Store the actual time together with the population of 
C*          the ground state and the population in the electronic 
C*          bound and continuum states in the formatted log file:
C*
         nt = 0
         DO 310 k=1,n1
            DO 300 k1=1,n(k)
               IF ( en(k1,k) .LE. 0.0D+00 ) THEN
                  prop = prop + ( y(k1+nt)**2 + y(ntot+nt+k1)**2 )
               ELSE
                  yield = yield + ( y(k1+nt)**2 + y(ntot+nt+k1)**2 )
               ENDIF
  300       CONTINUE
            nt = nt + n(k)
  310    CONTINUE
         ppgrnd = (y(1)**2+y(ntot+1)**2)
         prop = prop - ppgrnd
         dum = ppgrnd + prop + yield
         WRITE(10,1125) ppgrnd, prop, yield, dum
         WRITE(*,1125) ppgrnd, prop, yield, dum
         prop = 0.0D+00
         yield = 0.0D+00
C*
         t = tend
         nt = 0
C*
C*    End of the loop over intervals l:
  400 CONTINUE
C*
      ktend = kt
C*
      CLOSE (33)
C*
C*
C*        Open the unformatted output file (unit 34):
C*
      CALL getenv('OUTFIL',filnam)
      OPEN ( UNIT=34, FILE=filnam, STATUS='UNKNOWN',
     +       FORM='UNFORMATTED', ACCESS='SEQUENTIAL' )
C*
      WRITE(34) n1
      DO 410 kp=1,n1
         WRITE(34) n(kp)
         WRITE(34) (en(k1,kp),k1=1,n(kp))
  410 CONTINUE
      WRITE(34) ntot
      WRITE(34) (y(ii),ii=1,2*ntot)
C*
      CLOSE (34)
C*
      GOTO 950
C*
  420 CONTINUE
      nt = 0
C*
C*
C*       Store the energies, populations, and coefficients (real 
C*       and imaginary part) of all field-free states:
C* 
      CALL getenv('POPFIL',filnam)
      OPEN ( 40, FILE=filnam, STATUS='UNKNOWN',
     +       FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C*
      WRITE(40,*) ' Populations for time = ', t - dlimit*2.0D+00*ta
      DO 440 k=1,n1
         WRITE(40,*) ' Populations for l = ',k-1
         DO 430 k1=1,n(k)
            WRITE(40,'(2X,4(1PE14.6))') en(k1,k),
     +           (y(k1+nt)**2+y(ntot+nt+k1)**2),y(k1+nt),y(ntot+nt+k1)
  430    CONTINUE
         nt = nt + n(k)
  440 CONTINUE
      CLOSE (40)
C*
      GOTO 990
C*
  950 CONTINUE
C*
C*
C*       Store the time-dependent dipole spectrum as a function of
C*       the time:
C*
      CALL getenv('DIPOUT',filnam)
      OPEN ( 41, FILE=filnam, STATUS='UNKNOWN',
     +       FORM='FORMATTED', ACCESS='SEQUENTIAL' )
C*
      DO 970 kt = ktanf, ktend
         WRITE(41,'(2X,3E14.6)') 
     +                        time(kt)-dlimit*2.0D+00*ta, dipole(kt)
  970 CONTINUE
C*
      CLOSE (41)
C*
  990 CONTINUE
C*
C*       Format statements:
C*
 1002 FORMAT(//)
 1004 FORMAT('                         *********************')
 1005 FORMAT('                         ***               ***')
 1006 FORMAT('                         ***      TDSE     ***')
 1035 FORMAT(//,'   *****   ERROR !!!   *****',/)
 1038 FORMAT('   Unknown pulse type: ',I3,///)
 1040 FORMAT('   Input claims this is not a restart, but unit 33')
 1045 FORMAT('   is not empty, the program will stop!',///)
 1050 FORMAT('   Input claims this to be a restart, but unit 33')
 1055 FORMAT('   does not exist, the program will stop!',///)
 1060 FORMAT(3X,'Max. intensity (Epsilon_0^2): ',1PD26.10,' W/cm^2')
 1062 FORMAT(3X,'                              ',F16.10,' a.u.',/)
 1070 FORMAT(3X,'Photon energy               : ',F16.10,' eV')
 1072 FORMAT(3X,'                              ',F16.10,' a.u.',/)
 1075 FORMAT(3X,'The whole calculation is split into ',I3,
     +   ' intervals.')
 1077 FORMAT(3X,'Every interval is split into ',I3,
     +   ' sub-intervals')
 1078 FORMAT(6X,'at whose end the dipole spectrum is evaluated.',/)
 1080 FORMAT(3X,'Pulse length tau            : ',F16.9,' fs')
 1081 FORMAT(3X,'                              ',F16.9,' a.u.',/)
 1082 FORMAT(3X,'L_{max.} = ',I5,' ?!?',///)
 1085 FORMAT(3X,'L_{max.}            : ',I3)
 1090 FORMAT(3X,'E_{thr.} (in E_h)   : ',F16.10)
 1092 FORMAT(/,3X,'   +++  All states will be included  +++',/)
 1095 FORMAT(3X,'E_{cut-off} (in E_h): ',F16.10,/)
 1100 FORMAT(4I5)
 1105 FORMAT(3X,'Corresponding to ',F5.1,' opticle cycles!',/)
 1110 FORMAT(2X,8E15.7)
 1125 FORMAT(5X,4(E16.9,X),/)
 1130 FORMAT(1X,' Pulse duration (in a.u) = ',E15.8)
 1140 FORMAT(1X,' Peak intensity (in a.u) = ',E15.8)
 1150 FORMAT(3X,'Tolerance (for propagation): ',1PD12.5,/)
 1160 FORMAT(/,2X,'Interval no.: ',I4,'  Time(in a.u.): ',F16.10,
     +    ' to ',F16.10)
 1200 FORMAT(/,3X,'+++++  This is a restarted run  +++++',/)
 1205 FORMAT(/,3X,'+++++  This run only rewrites existing data',
     +   '  +++++',/)
 1210 FORMAT(3X,'===  The names of the input files are read',
     +   ' ! ===',/) 
 1220 FORMAT(3X,'===  The names of the input files are generated',
     +   ' !  ===',/)
 1230 FORMAT(3X,'For L = ',I3,'  the number of states included is: ',
     +   I6)
 1240 FORMAT(/,3X,'   --->  Total number of states is: ',I7,//)
 1250 FORMAT(6X,'Number of function evaluations: ',I8)
 1260 FORMAT(3X,'***  A cos^2 pulse (for E(t)) is used !  ***',/)
 1270 FORMAT(3X,'***  A cos^2 pulse (for A(t)) is used !  ***',/)
 1280 FORMAT(3X,'***  A Gaussian pulse (for A(t)) is used !  ***',/)
 1300 FORMAT('   Non-unique definition of the pulse duration!')
 1310 FORMAT('   Either number of cycles or its duration should be ',
     +   '.LE. 0')
 1320 FORMAT('   Found: Cycles = ',D15.8,'  duration = ',D15.8,
     +   ///)
C*
      STOP
      END
C*
C*
C***********************************************************************
C*
C*IMSL  1
      SUBROUTINE fcn(neq,t,y,ypr)
C*
C*NAG   1
C      SUBROUTINE fcn(t,y,ypr)
C***
C*
C*======================================================================
C*
C*    This subroutine evaluates the values of the derivative YPR of Y 
C*    where Y is a function of the independent variable T (which is 
C*    here equivalent to the time t). 
C*
C*    Function y_n(t) is the time-dependent coefficient of the field-
C*    free wavefunction n forming the basis for the time-dependent
C*    wavefunction to be calculated by this program. More accurately, 
C*    for  (0 .LE. n .LE. ntot) y_n is the real part of the coefficient 
C*    of wavefunction n, while for ( ntot .LT. n .LE. 2*ntot) it is the 
C*    imaginary part of that coefficient. 
C*
C*    Due to the ortho-normality of the field-free states the time-
C*    dependent Schroedinger equation is transformed into (see e.g. 
C*    eq.(69) in Phys.Rep. vol.305, 203 (1998))
C*  
C*                   ___
C*         d y_n     \    /                               \
C*       i -----  =   >   | E_m \delta_{n,m} - V_{n,m}(t) |  y_m    
C*          d t      /    \                               /
C*                   ---
C*                    m
C* 
C*    what can be transformed into
C*                                         ___
C*       d y_n                             \                      
C*       -----  =  - i * E_n * y_n  +  i *  >   V_{n,m}(t) *  y_m    
C*        d t                              /     
C*                                         ---
C*                                          m
C*
C*    Here E_n is the (field-free) energy of state n, while V_{n,m} is
C*    the interaction matrix element connecting states n and m. Here it
C*    is the dipole interaction operator, either in length or velocity
C*    gauge. Thus V_{n,m}(t) is the dipole matrix element (in either 
C*    gauge) between wavefunctions n and m times the time-dependent 
C*    pulse amplitude A(t).
C* 
C*    In the actual calculation the set of equations (for all n 
C*    coefficients) is manipulated in algebraic form by rewriting the
C*    set of equations in matrix form. In addition it is used that due
C*    to the dipole selction rule V_{n,m} is non-zero only, if the 
C*    angular momenta corresponding to states n and m exactly differ by
C*    one. Finally, the manipulations are done separately for the real
C*    and imaginary parts of dy/dt.
C*
C*    The matrix/vector calculations are performed using ESSL routine
C*    DGEMV. 
C* 
C*----------------------------------------------------------------------
C*
C*    Author: Xian Tang
C*    Date  : January 1990
C*    Modified by:  H. Rudolph (?), Jian Zhang, Paul Maragakis, and
C*         Alejandro Saenz (since February 1999)
C*    Last Change:  30.08.1999
C*
C*----------------------------------------------------------------------
C*
      IMPLICIT REAL*8 (a-h,o-z)
C*
C*  
C*         Parameter:   
C*         ==========
C*
C*             lvlmax: max. number of L values.
C*             nvlmax: max. number of eigenstates per L value.
C*             nstmax: (not used in this routine).
C*             ncfmax (=lvlmax*nvlmax): max. number of coefficients,
C*                     i.e. max. number of field-free states.
C*
      INCLUDE 'param_tdse.f'  
C*
      PARAMETER ( ncfmax = lvlmax * nvlmax )
C*
      DIMENSION y(2*ncfmax),ypr(2*ncfmax)
C*
      COMMON/f/am(nvlmax,nvlmax,lvlmax-1), en(nvlmax,lvlmax), n(lvlmax), 
     +     nsum(lvlmax), n1, ntot, fi, omag, ta, neval
      COMMON/f1/dag(ncfmax)
      COMMON/pival/pi
      COMMON/ivar/ipulse
      COMMON/dvar/dlimit
C*
C*----------------------------------------------------------------------
C*
      zero = 0.0D+00
C*
      nvl = nvlmax
      neval = neval + 1
C*
C*
C*          Calculate the amplitude A(t):
C*
C*               1) for a cos^2 pulse (with respect to E):
C*
      IF (ipulse .EQ. 1) THEN
        t1 = t
        om = pi / ( 2.0D+00 * ta )
        om1 = omag + ( 2.0D+00 * om )
        om2 = omag - ( 2.0D+00 * om )
        half = 0.5D+00
        slq = - half * fi * ( 
     +          ( (1.0D+00 - DCOS(omag*t1) ) / omag )
     +          - half * ( ( ( 1.0D+00 - DCOS(om1*t1) )/ om1 ) 
     +                   + ( ( 1.0D+00 - DCOS(om2*t1) ) /om2 ) )
     +                  )
      ENDIF
C*
C*
C*               2) for a cos^2 pulse (with respect to A):
C*
      IF (ipulse .EQ. 2) THEN
        t1 = - ( dlimit * 2.0D+00 * ta ) + t
        sl = DCOS( ( t1 * pi ) / ( 2.0D+00 * ta ) )**2
        slq = sl * fi * DCOS( omag * t1 ) / omag
      ENDIF
C*
C*
C*               3) for a Gaussian pulse (with respect to A):
C*
      IF (ipulse .EQ. 3) THEN
        t1 = - ( dlimit * 2.0D+00 * ta ) + t
        t2 = t1 * t1 / ( 2.0D+00 * ta * ta )
        sl = DEXP( - t2 )
        slq = sl * fi * DCOS( omag * t1 ) / omag
      ENDIF
C*
C*
C*          Calculate ypr(n) = - i * E(n) * y(n), 
C*                  i.e. Re{ypr(n)} = + E(n) * Im{y(n)}
C*                       Im{ypr(n)} = - E(n) * Re{y(n)}:
C*
      DO  i=1,ntot
        ypr(i)      =   dag(i) * y(ntot + i)
        ypr(ntot+i) = - dag(i) * y(i)
       enddo
C*
C*
C*          Calculate ypr(n) = ypr(n) + i * \sum_m  V_{n,m} * y(m),
C*            i.e. Re{ypr(n)} = Re{ypr(n)} - \sum_m V_{n,m} * Im{y(m)}
C*                 Im{ypr(n)} = Im{ypr(n)} + \sum_m V_{n,m} * Re{y(m)}
C*          for all angular momenta:
C* 
      CALL DGEMV('t', n(2), n(1), -slq, am(1,1,1), nvl,
     +     y(ntot+1+n(1)), 1, 1.0D+00, ypr(1), 1)
      CALL DGEMV('t', n(2), n(1), slq, am(1,1,1), nvl,
     +     y(1+n(1)), 1, 1.0D+00, ypr(ntot+1), 1)
C*

C n1 = lmax + 1

      DO 80 k=2,n1-1
         CALL DGEMV('t', n(1+k), n(k), -slq, am(1,1,k), nvl,
     +        y(ntot+nsum(k+1)), 1, 1.0D+00, ypr(nsum(k)), 1)
         CALL DGEMV('t', n(1+k), n(k), slq, am(1,1,k), nvl,
     +        y(nsum(k+1)), 1, 1.0D+00, ypr(ntot+nsum(k)), 1)

         CALL DGEMV('n', n(k), n(k-1), -slq, am(1,1,k-1), nvl,
     +        y(ntot+nsum(k-1)), 1, 1.0D+00, ypr(nsum(k)), 1)
         CALL DGEMV('n', n(k), n(k-1), slq, am(1,1,k-1), nvl,
     +        y(nsum(k-1)), 1, 1.0D+00, ypr(ntot+nsum(k)), 1)
   80 CONTINUE
C*
      CALL DGEMV('n', n(n1), n(n1-1), -slq, am(1,1,n1-1), nvl,
     +     y(ntot+nsum(n1-1)), 1, 1.0D+00, ypr(nsum(n1)), 1)
      CALL DGEMV('n', n(n1), n(n1-1), slq, am(1,1,n1-1), nvl,
     +     y(nsum(n1-1)), 1, 1.0D+00, ypr(ntot+nsum(n1)), 1)
C*
      RETURN
      END
C*
C*
C***********************************************************************
C*
      SUBROUTINE pulse(t0,tend,fi,omag,ta,dlimit,ipulse)
C*
C*======================================================================
C*
C*    Calculate and store the pulse.
C* 
C*----------------------------------------------------------------------
C*
C*    Author: Alejandro Saenz
C*    Date  : 21.06.1999
C*    Last Change: 
C*
C*----------------------------------------------------------------------
C*
      IMPLICIT REAL*8 (a-h,o-z)
C*
      CHARACTER*256 filnam
C*
C*----------------------------------------------------------------------
C*
      pi = DACOS(-1.0D+00)
      zero = 0.0D+00
C*
      tstep = ( tend - t0 ) / DFLOAT(1001)
C*
      CALL getenv('PULSE',filnam)
      OPEN ( UNIT=99, FILE=filnam, STATUS='UNKNOWN',
     +       ACCESS='SEQUENTIAL', FORM='FORMATTED' )
C*
      DO 100 i=1,1000
        t = t0 + DFLOAT(i-1) * tstep
C*
C*
C*          Calculate the amplitude A(t):
C*
C*               1) for a cos^2 pulse (with respect to E):
C*
        IF (ipulse .EQ. 1) THEN
          t1 = t
          om = pi / ( 2.0D+00 * ta )
          om1 = omag + ( 2.0D+00 * om )
          om2 = omag - ( 2.0D+00 * om )
          half = 0.5D+00
          slq = - half * fi * ( 
     +          ( (1.0D+00 - DCOS(omag*t1) ) / omag )
     +          - half * ( ( ( 1.0D+00 - DCOS(om1*t1) )/ om1 ) 
     +                   + ( ( 1.0D+00 - DCOS(om2*t1) ) /om2 ) )
     +                  )
        ENDIF
C*
C*
C*               2) for a cos^2 pulse (with respect to A):
C*
        IF (ipulse .EQ. 2) THEN
          t1 = - ( dlimit * 2.0D+00 * ta ) + t
          sl = DCOS( ( t1 * pi ) / ( 2.0D+00 * ta ) )**2
          slq = sl * fi * DCOS( omag * t1 ) / omag
        ENDIF
C*
C*
C*               3) for a Gaussian pulse (with respect to A):
C*
        IF (ipulse .EQ. 3) THEN
          t1 = - ( dlimit * 2.0D+00 * ta ) + t
          t2 = t1 * t1 / ( 2.0D+00 * ta * ta )
          sl = DEXP( - t2 )
          slq = sl * fi * DCOS( omag * t1 ) / omag
        ENDIF
C*
        WRITE(99,*) t,slq
C*
  100 CONTINUE
C*
      CLOSE(99)
C*
      RETURN
      END
C*
C*
C***********************************************************************
C*
C*IMSL  1
      SUBROUTINE fcn2(neq,t,y,ypr)
C*
C*NAG   1
C      SUBROUTINE fcn2(t,y,ypr)
C***
C*
C*======================================================================
C*
C*    This subroutine evaluates the values of the derivative YPR of Y 
C*    where Y is a function of the independent variable T (which is 
C*    here equivalent to the time t). 
C*
C*    Function y_n(t) is the time-dependent coefficient of the field-
C*    free wavefunction n forming the basis for the time-dependent
C*    wavefunction to be calculated by this program. More accurately, 
C*    for  (0 .LE. n .LE. ntot) y_n is the real part of the coefficient 
C*    of wavefunction n, while for ( ntot .LT. n .LE. 2*ntot) it is the 
C*    imaginary part of that coefficient. 
C*
C*    Due to the ortho-normality of the field-free states the time-
C*    dependent Schroedinger equation is transformed into (see e.g. 
C*    eq.(69) in Phys.Rep. vol.305, 203 (1998))
C*  
C*                   ___
C*         d y_n     \    /                               \
C*       i -----  =   >   | E_m \delta_{n,m} - V_{n,m}(t) |  y_m    
C*          d t      /    \                               /
C*                   ---
C*                    m
C* 
C*    what can be transformed into
C*                                         ___
C*       d y_n                             \                      
C*       -----  =  - i * E_n * y_n  +  i *  >   V_{n,m}(t) *  y_m    
C*        d t                              /     
C*                                         ---
C*                                          m
C*
C*    Here E_n is the (field-free) energy of state n, while V_{n,m} is
C*    the interaction matrix element connecting states n and m. Here it
C*    is the dipole interaction operator, either in length or velocity
C*    gauge. Thus V_{n,m}(t) is the dipole matrix element (in either 
C*    gauge) between wavefunctions n and m times the time-dependent 
C*    pulse amplitude A(t).
C* 
C*    In the actual calculation the set of equations (for all n 
C*    coefficients) is manipulated in algebraic form by rewriting the
C*    set of equations in matrix form. In addition it is used that due
C*    to the dipole selction rule V_{n,m} is non-zero only, if the 
C*    angular momenta corresponding to states n and m exactly differ by
C*    one. Finally, the manipulations are done separately for the real
C*    and imaginary parts of dy/dt.
C*
C*    The matrix/vector calculations are performed using ESSL routine
C*    DGEMV. 
C* 
C*----------------------------------------------------------------------
C*
C*    Author: Xian Tang
C*    Date  : January 1990
C*    Modified by:  H. Rudolph (?), Jian Zhang, Paul Maragakis, and
C*         Alejandro Saenz (since February 1999)
C*    Last Change:  30.08.1999
C*
C*----------------------------------------------------------------------
C*
      IMPLICIT REAL*8 (a-h,o-z)
C*
C*  
C*         Parameter:   
C*         ==========
C*
C*             lvlmax: max. number of L values.
C*             nvlmax: max. number of eigenstates per L value.
C*             nstmax: (not used in this routine).
C*             ncfmax (=lvlmax*nvlmax): max. number of coefficients,
C*                     i.e. max. number of field-free states.
C*
      INCLUDE 'param_tdse.f' 
C*
      PARAMETER ( ncfmax = lvlmax * nvlmax )
C*
      DIMENSION y(2*ncfmax),ypr(2*ncfmax)
C*
      COMMON/f/am(nvlmax,nvlmax,lvlmax-1), en(nvlmax,lvlmax), n(lvlmax), 
     +     nsum(lvlmax), n1, ntot, fi, omag, ta, neval
      COMMON/f1/dag(ncfmax)
      COMMON/pival/pi
      COMMON/ivar/ipulse
      COMMON/dvar/dlimit
C*
C*----------------------------------------------------------------------
C*
      nvl = nvlmax
      neval = neval + 1
C*
C*
C*          Calculate the amplitude A(t):
C*
C*               1) for a cos^2 pulse (with respect to E):
C*
      IF (ipulse .EQ. 1) THEN
        t1 = t
        om = pi / ( 2.0D+00 * ta )
        om1 = omag + ( 2.0D+00 * om )
        om2 = omag - ( 2.0D+00 * om )
        half = 0.5D+00
        slq = - half * fi * ( 
     +          ( (1.0D+00 - DCOS(omag*t1) ) / omag )
     +          - half * ( ( ( 1.0D+00 - DCOS(om1*t1) )/ om1 ) 
     +                   + ( ( 1.0D+00 - DCOS(om2*t1) ) /om2 ) )
     +                  )
      ENDIF
C*
C*
C*               2) for a cos^2 pulse (with respect to A):
C*
      IF (ipulse .EQ. 2) THEN
        t1 = - ( dlimit * 2.0D+00 * ta ) + t
        sl = DCOS( ( t1 * pi ) / ( 2.0D+00 * ta ) )**2
        slq = sl * fi * DCOS( omag * t1 ) / omag
      ENDIF
C*
C*
C*               3) for a Gaussian pulse (with respect to A):
C*
      IF (ipulse .EQ. 3) THEN
        t1 = - ( dlimit * 2.0D+00 * ta ) + t
        t2 = t1 * t1 / ( 2.0D+00 * ta * ta )
        sl = DEXP( - t2 )
        slq = sl * fi * DCOS( omag * t1 ) / omag
      ENDIF
C*
C*
C*          Calculate ypr(n) = - i * E(n) * y(n), 
C*                  i.e. Re{ypr(n)} = + E(n) * Im{y(n)}
C*                       Im{ypr(n)} = - E(n) * Re{y(n)}:
C*
      DO  i=1,ntot
        ypr(i) = dag(i) * y(i+ntot)
        ypr(i+ntot) = - dag(i) * y(i)
      enddo
C     *
C*       i*ypr = e*yr + V*yr ==>    (V == real) 
C*
C*          Calculate ypr(n) = ypr(n) + i * \sum_m  V_{n,m} * y(m),
C*            i.e. Re{ypr(n)} = Re{ypr(n)} - \sum_m V_{n,m} * Im{y(m)}
C*                 Im{ypr(n)} = Im{ypr(n)} + \sum_m V_{n,m} * Re{y(m)}
C*          for all angular momenta:
C* 
      ind = 0
      DO  i=1,n1
        DO  j=1,n(i)

          ind = ind + 1

          IF (i.NE.1) THEN
            jnd = nsum(i-1) - 1
            DO  jj=1,n(i-1)

              jnd = jnd + 1 
           
      ypr(   0+ind) = ypr(ind)      - slq * am(j,jj,i-1) * y(ntot + jnd) !
      ypr(ntot+ind) = ypr(ntot+ind) + slq * am(j,jj,i-1) * y(  0  + jnd) !

      enddo

      ENDIF
      
      IF (i.NE.n1) THEN

         jnd = nsum(i+1) - 1
         DO  jj=1,n(i+1)

           jnd = jnd + 1
      ypr(ind)      = ypr(ind)      + slq * am(jj,j,i) * y(ntot+jnd) !ok
      ypr(ntot+ind) = ypr(ntot+ind) - slq * am(jj,j,i) * y(jnd)      !ok
         enddo

      ENDIF
      
      enddo
      enddo
C*
      RETURN
      END
C*EOF
