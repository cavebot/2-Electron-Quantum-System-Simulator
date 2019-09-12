C#######################################################################
C#
C#      PROGRAM TO CALCULATE THE B-SPLINE WAVE FUNCTIONS
C#      FUNCTION BVALUE(T,C,N,K,X,NDR) IS USED
C#
C#
C#######################################################################


C#
C#     l.aa.n/iesl/2001 :  
C# 
C#                     spring 2001 :
C#                 --- Large modifications on I/O code fragments
C#                 --- Removal of several obsolete variables 
C#                  
C#
C#     20/09/2001 :   working on PARAMETER inputs of the programms
C#
C#

C#
C#   FOR  NS,NK,NP parameters see d1e.f 
C#
C#      RS  : FIRST KNOT POINT
C#      RMX : RADIUS OF BOX 
C#      NB  : NUMBER OF B-SPLINES  
C#      K   : ORDER OF B-SPLINES
C#      LIN : ?  

      PROGRAM WF1E

      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.1e.inc"
      DIMENSION R(NP),T(NS+NK)
      DIMENSION COE(NS), COEIN(NS-2,NS)
      DIMENSION P(NP),DP(NP)             
      DIMENSION NVALENCE(NK)
      DIMENSION FF(NP),DR(NP),EN(NS)
      INTEGER NPLOT,LPLOT
      double precision de,rho
      double precision enau
C....................................................

      enau = 27.211396181D+00

C#   read the standard input file

      OPEN(10,FILE='inp/wf1e.inp')
      READ(10, *) llmin, llmx
      READ(10, *) (nvalence(i), i = 1, llmx - llmin + 1) 
      READ(10, *) nplot,lplot 
      CLOSE(10)

C======

      OPEN(15,FILE ='dat/en1e.dat')
      OPEN(16,FILE ='out/wf1e.out')
      OPEN(17,FILE ='out/en1e.out')
      OPEN(27,FILE ='dat/knot1e.dat')

      WRITE(15,11) LLMIN, LLMX
      WRITE(17,*) 'n, en(n) (eV), de (a.u.), rho (a.u.)' 

      CALL RIN( R, DR, H, NO, IDR)

      DO  IL = LLMIN, LLMX

         LANG = IL - 1

      WRITE(17,'(I3)') LANG

C* read from data files (en+coefficients) 

         call d1efile(1, lang) 

         READ(1) RS, RMX, NB, K, LIN

         IF(LIN.NE.LANG)               STOP
         IF(NB.GT.NS .OR. NB.LT.(K+2)) STOP
         IF(K.GT.NK .OR. K.LT.2)       STOP


C#      make grid points and store in:  t(nb + k)

         CALL MKGRID(NB, K, RMX, RS, T)
                  

C#        Write knot points for l = llmin in 'KnotFile'


         IF(IL.EQ.LLMIN) THEN

            DO I = 1, NB + K
               WRITE(27,*)  I,' ', T(I) 
            ENDDO

            CLOSE(27)

         ENDIF

C#          COEIN(n,m)   n (state index),  m ( b-spline index)
C#
C#      set    P(0)=P(R)=0 ==>  c_1 = c_n = 0
C#      and    read energy and B-splines coeff of state 'npo'

         NPMX = NB - 2
         NM   = NB - 1
         NJJ  = 0

         DO  NPP = 1, NPMX
            
            NPO  = NPMX - NPP + 1

            COEIN( NPO, 1)  = 0.0D+00
            COEIN( NPO, NB) = 0.0D+00
           
            READ(1) EN(NPO), ( COEIN(NPO, I), I = 2, NM)

            NJJ = NJJ + 1
         ENDDO


         WRITE(15, 11) NJJ, NVALENCE(IL)
C         WRITE(17, 11) NJJ, NVALENCE(IL)

         CLOSE(1)

C...

C# calculate radial wf  P(r) = Sum_i c_i * B_i(r)


         call wf1efile(3,lang) 

         DO NPP = 1, NJJ
            
            DO  I = 1, NB
               COE(I) = COEIN(NPP, I)
            ENDDO

            P(1) = 0.0D+00 
            DO  J = 2, NO 
               P(J) = BVALUE(T, COE, NB, K, R(J), 0)
            ENDDO

            DO  J = 1, NO 
               DP(J) = BVALUE(T, COE, NB, K, R(J), 1)
            ENDDO

            DP(no) = -(k-1) *coe(nb-1)/(t(nb+1) - t(nb))


c            FF(1) = 0.0D+00
c            DO  J = 2, NO
c               FF(J) = P(J) * P(J) * R(J) / DR(J)
c            ENDDO
c
c            OTH(NP) = RINT(FF,1,NO,9,H)
c            SQ2     = DSQRT(OTH(NP))
c            DO 46  J = 1, NO
c            P(J) = P(J) / SQ2
c            ENDDO
            
            IF(LANG.EQ.LPLOT.AND.NPP.EQ.NPLOT) THEN
               DO J = 1, NO
                 WRITE(16,*)  R(J), P(J), DP(J)
               ENDDO
            ENDIF

            WRITE(3) ( P(J), J = 1, NO )

            de = en(npp+1)-en(npp-1) 
            rho = 2.0d+00/abs(de)

            WRITE(15,'(E20.14)')       en(npp)
            if(en(npp).lt.8) then
            WRITE(17,'(I5,1X,3G15.5)')   npp,en(npp)*enau*0.5,DE,rho
            else
            WRITE(17,'(I5,1X,2P3E15.5)') npp,en(npp)*enau*0.5,DE,rho
            endif 

            IF(EN(NPP).LE.0.0D+00) THEN
               WRITE(*,*) NPP,'.........', EN(NPP)
            ENDIF

          ENDDO
         
          CLOSE(3)
      
       ENDDO
      
       close(15)
       close(16)
       close(17)

       write(*,*)'# wf1e:: P(r) in wf1e.out n,l = ',nplot,lplot
C=====================    format statements

    3 FORMAT(1X,'RS,RMX,H,NB,K,NO,L=',1P3E13.5,4I5/)
    4 FORMAT(1X,1PE11.3,1P6E17.9)
    5 FORMAT(/1X,'KNOTS: T'/(1X,1P5D14.6))
    6 FORMAT(/2X,'ORTHOGONALITY --'/)
    7 FORMAT(1X,1P10E13.4)
    8 FORMAT(/1X,'EG(n=',I3,') =', 1P2E19.10/)
    9 FORMAT(1X)
   10 FORMAT(1X,1P5E15.7)
   11 FORMAT(5I5)
   12 FORMAT(1X,'LLMIN & LLMAX=?'/)
 566  FORMAT(5A16)


      END

C###################################################################
C#
C#     OTH,SQ2 
C#      
C#     OTH  ==  | < wf_n | wf_n > |^2   ( == 1 )   
C#
C#     Discretized WaveFunctions  (bound and continuum)  have a norm 
C#     inside the interval [0,R] equal to 1.
C#
C#     Normalization of bound states is correct ( = 1 ) ,
C#   
C#     Bound States :
C#                               {  0   n != m
C#     <wf_n | wf_m > =   d_nm = {                   E_n , E_m < 0 
C#                               {  1   n  = m  
C#     however 
C#     continuum solutions have to be energy- normalized according to
C#
C#     Continuum States :
C#                                        {     0      n != m
C#     <wf_n | wf_n > = delta (En - Em) = {
C#                                        {  infinity  n  = m 
C#
C#      Therefore for calculations of Dipole Matrix Elements of 
C#      type
C#   
C#      BC         Bound      -  Continuum
C#      CC         Continuum  -  Continuum
C#
C#      an additional normalization factor is needed for the 
C#      contnuum state.
C#
C#      Here the "unnormalized" discrete - continuum states 
C#      are stored. 
C#
C###################################################################

C#EOF












