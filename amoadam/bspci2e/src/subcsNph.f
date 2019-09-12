C#
C#   csNph.f subroutines
C#
C#
C# GET THE ENERGY OF THE FIRST IONIZATION THRESHOLD
C# THE FOLLOWING APPLIES ONLY TO NS^2 ATOMS (Hn,He,Mg,Ca)
C# 
C####################################################################
      SUBROUTINE READ_ION_THRESHOLD(E_IONIC)

      IMPLICIT NONE
      DOUBLE PRECISION ENAU
      PARAMETER ( ENAU     = 27.211 396 181  D+00 )
      INTEGER NINP,I
      INTEGER LMIN_1E, LMAX_1E
      INTEGER NOF_STATES, NVALENCE
      CHARACTER*20 EN1EFILE
      DOUBLE PRECISION E_IONIC
C..............................

      NINP = 9 
C................
      EN1EFILE ='dat/en1e.dat'

      OPEN(NINP, FILE=EN1EFILE)

      READ(NINP,*) LMIN_1E, LMAX_1E

      IF((LMIN_1E - 1).NE.0) THEN
         
         WRITE(*,*) '# NORM1E::  ION ENERGIES NOT FOR L = 0. EXITING..'
         STOP
      ENDIF
 
      READ(NINP,*) NOF_STATES, NVALENCE

      DO  I = 1, NVALENCE
         READ(NINP,*)  E_IONIC
      ENDDO

      CLOSE(NINP)

      E_IONIC = 0.5D+00 * E_IONIC

      WRITE(*,*)'# read_ion_threshold::'
      WRITE(*,*)'# read_ion_threshold::    A+ =', E_IONIC," au = "
     & ,E_IONIC*ENAU," eV "


      RETURN
      END
C####################################################################
c      FUNCTION SIGMA(CHEMIN,W,DIP,NNINIT,NBPHOT,NNFIN,EPSILON)
      FUNCTION SIGMA(CHEMIN,W,DIP,NNINIT,NBPHOT,NNFIN,EPSILON,
     &     ENERG,RNORME,PHASE,NDIMENSION)
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               

      PARAMETER ( MAXPHOT = 4, MAXPATH = 10)
      PARAMETER ( MAXDIM = 2000, MAXL = 4 )
      COMPLEX*16 SIGMA, B, A1, A2, VF, VI, VTEMP, DNIF
      COMPLEX*16 ONE, ZERO, ZDOTU
      PARAMETER ( ONE=(1.D0,0.D0), ZERO=(0.D0,0.D0) ) 
      INTEGER   CHEMIN,ARBRELM                               
      DIMENSION VI(MAXDIM),VF(MAXDIM),B(MAXDIM,MAXDIM)                  
      DIMENSION VTEMP(MAXDIM)
      DIMENSION A1(MAXDIM,MAXDIM),A2(MAXDIM,MAXDIM)                     
      DIMENSION DIP(MAXDIM,MAXDIM,MAXL)
!      DIMENSION ENERG(MAXDIM,MAXL)
      DIMENSION EPSILON( MAXL ) 
      INTEGER   NDIMENSION
      DIMENSION ENERG(MAXDIM,0:MAXL),RNORME(MAXDIM,0:MAXL)
      DIMENSION PHASE(MAXDIM,0:MAXL)
      COMMON /PATHS/ ARBRELM(MAXPATH,MAXPHOT,2),ANGULAR(MAXPATH)
c      COMMON /STATE/ ENERG(MAXDIM,0:MAXL),RNORME(MAXDIM,0:MAXL), 
c     &               PHASE(MAXDIM,0:MAXL),NDIMENSION   
C......................................
CAN(MAXDIM,0:MAXL)


      LINIT = ARBRELM( CHEMIN, 1, 1 )                                       
       LFIN = ARBRELM( CHEMIN, NBPHOT+1, 1 )                                 
      EINIT =   ENERG( NNINIT, LINIT )


C#                                                                       
C# 1 PHOTON ABSORPTION CASE
C#

      IF (NBPHOT.EQ.1) THEN

	  L1 = ARBRELM(CHEMIN,1,1)
	  M1 = ARBRELM(CHEMIN,1,2)
	  L2 = ARBRELM(CHEMIN,2,1)
	  M2 = ARBRELM(CHEMIN,2,2)
           L = MAX0(L1,L2)
         NN1 = NNINIT
         NN2 = NNFIN

         IF (L1.LT.L2) THEN

	    DNIF = DCMPLX( DIP(NN1,NN2,MAX0(L1,L2)) *
     S           RNORME( NNFIN,ARBRELM(CHEMIN,NBPHOT+1,1)),0.D0 ) 
         ELSE

	    DNIF = DCMPLX( DIP(NN2,NN1,MAX0(L1,L2)) *
     S           RNORME( NNFIN,ARBRELM(CHEMIN,NBPHOT+1,1)),0.D0 )
         ENDIF
C

	 SIGMA = DNIF * DCMPLX( ANGULAR(CHEMIN), 0.D0 )




	 RETURN
      ELSE
C                                                                       
C CHARGEMENT DE VII1                                                    
C                                                                       
                                                                        
      L1 = ARBRELM(CHEMIN,1,1)                                          
      M1 = ARBRELM(CHEMIN,1,2)                                          
      L2 = ARBRELM(CHEMIN,2,1)                                          
      M2 = ARBRELM(CHEMIN,2,2)                                          
      NN1 = NNINIT 

      IF (L1.LT.L2) THEN                                                     
         DO NN2=1,NDIMENSION                                             
            VI(NN2) = DCMPLX( DIP(NN1,NN2,MAX0(L1,L2)), 0.D0 ) /
     S                DCMPLX( EINIT-ENERG(NN2,L2)+W, EPSILON(1) )
         ENDDO
      ELSE
         DO NN2=1,NDIMENSION                                             
            VI(NN2) = DCMPLX( DIP(NN2,NN1,MAX0(L1,L2)), 0.D0 )/
     S                DCMPLX( EINIT-ENERG(NN2,L2)+W, EPSILON(1) )
         ENDDO
      ENDIF
C                                                                       
C CHARGEMENT DE VFIN-1                                                  
C                                                                       
      L1 = ARBRELM(CHEMIN,NBPHOT,1)                                     
      M1 = ARBRELM(CHEMIN,NBPHOT,2)                                     
      L2 = ARBRELM(CHEMIN,NBPHOT+1,1)                                   
      M2 = ARBRELM(CHEMIN,NBPHOT+1,2)                                   
      NN2 = NNFIN                                                       
      IF (L1.LT.L2) THEN
         DO NN1=1,NDIMENSION                                             
            VF(NN1) = DCMPLX( DIP(NN1,NN2,MAX0(L1,L2)), 0.D0 )
         ENDDO
      ELSE
         DO NN1=1,NDIMENSION                                             
            VF(NN1) = DCMPLX( DIP(NN2,NN1,MAX0(L1,L2)), 0.D0 )
         ENDDO
      ENDIF
C                                                                       
C CHARGEMENT DE LA PREMI}RE MATRICE DE DROITE POUR LE CALCUL DES        
C PRODUITS.                                                             
C                                                                       
      L1 = ARBRELM(CHEMIN,NBPHOT-1,1)                                   
      M1 = ARBRELM(CHEMIN,NBPHOT-1,2)                                   
      L2 = ARBRELM(CHEMIN,NBPHOT,1)                                     
      M2 = ARBRELM(CHEMIN,NBPHOT,2)                                     
      IF (NBPHOT.LE.2) THEN                                     
         DO NN2=1,NDIMENSION                                             
            DO NN1=1,NDIMENSION     
               B(NN1,NN2) = 0.D0 
            ENDDO
         ENDDO
         DO NN1=1,NDIMENSION                                             
            B(NN1,NN1) = 1.D0 
         ENDDO
      ELSE
         IF (L1.LT.L2) THEN
            DO NN2=1,NDIMENSION   
               DO NN1=1,NDIMENSION                                          
		  B(NN1,NN2) = DCMPLX( DIP(NN1,NN2,MAX0(L1,L2)), 0.D0) /
     S                         DCMPLX( EINIT-ENERG(NN2,L2)+(NBPHOT-1)*W, 
     S                                 EPSILON( NBPHOT - 1 ) )
               ENDDO
            ENDDO
         ELSE
            DO NN1=1,NDIMENSION                                             
               DO NN2=1,NDIMENSION                                          
		  B(NN1,NN2) = DCMPLX( DIP(NN2,NN1,MAX0(L1,L2)),
     S                                 0.D0 ) /
     S                       DCMPLX( EINIT-ENERG(NN2,L2)+(NBPHOT-1)*W, 
     S                                EPSILON( NBPHOT - 1 ) )
               ENDDO
            ENDDO
         ENDIF
      ENDIF
C                                                                       
C CALCUL DES PRODUITS DE MATRICES INTERNES                              
C        
                                                               
      DO NIVEAU=NBPHOT-3,1,-1                                         
         L1 = ARBRELM(CHEMIN,NIVEAU+1,1)                                
         M1 = ARBRELM(CHEMIN,NIVEAU+1,2)                                
         L2 = ARBRELM(CHEMIN,NIVEAU+2,1)                                
         M2 = ARBRELM(CHEMIN,NIVEAU+2,2)                                
         IF (L1.LT.L2) THEN
            DO NN2=1,NDIMENSION                                     
               DO NN1=1,NDIMENSION                                
	          A1(NN1,NN2) = DCMPLX( DIP(NN1,NN2,MAX0(L1,L2)),0.D0 ) / 
     S              DCMPLX( EINIT-ENERG(NN2,L2)+(NIVEAU+1)*W, 
     S                      EPSILON( NIVEAU + 1 ) )
                  A2(NN1,NN2) = B(NN1,NN2)        
               ENDDO               
            ENDDO
         ELSE
            DO NN1=1,NDIMENSION                                     
               DO NN2=1,NDIMENSION      
	          A1(NN1,NN2) = DCMPLX( DIP(NN2,NN1,MAX0(L1,L2)),
     S                                  0.D0 ) /
     S              DCMPLX( EINIT-ENERG(NN2,L2)+(NIVEAU+1)*W,      
     S                      EPSILON( NIVEAU + 1 ) )
                  A2(NN1,NN2) = B(NN1,NN2)   
               ENDDO               
            ENDDO
         ENDIF

C     CALL PMAT(A1,A2,B,NDIMENSION,MAXDIM)

	 CALL ZGEMM('N','N',NDIMENSION,NDIMENSION,NDIMENSION,ONE,
     S              A1,MAXDIM,A2,MAXDIM,ZERO,B,MAXDIM)

      ENDDO                                                          

      CALL ZGEMV('N',NDIMENSION,NDIMENSION,ONE,B,MAXDIM,VF,1,
     S                                     ZERO,VTEMP,1)

      DNIF = ZDOTU(NDIMENSION,VI,1,VTEMP,1) *
     S       DCMPLX( RNORME(NNFIN,ARBRELM(CHEMIN,NBPHOT+1,1)), 0.D0 )
C                                                                       
      SIGMA = DNIF * ANGULAR(CHEMIN)

      ENDIF
C
      RETURN                                                            
      END                                                               
C                                                                       
C##############################################################################
      SUBROUTINE READDIPO(DIPFILE,ENFILE,D1,EN,RN,PHASE,AN,NDIM,
     S             LMIN,LMAX,MAXDIM,MAXL,ANGULPARTINCLUDED, MODE)

C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER MODE
      LOGICAL ANGULPARTINCLUDED
      CHARACTER*14 DIPFILE,ENFILE
      DIMENSION D1(MAXDIM,MAXDIM,MAXL)
      DIMENSION EN(MAXDIM,0:MAXL),RN(MAXDIM,0:MAXL)
      DIMENSION PHASE(MAXDIM,0:MAXL),AN(MAXDIM,0:MAXL)
      INTEGER NDIM
C==================================
 

      DIPFILE ='dat/dmx.dat'
      ENFILE = 'dat/en.dat' 

      WRITE(*,*) "# csNph_L::read_dmx   file opened:" 
      WRITE(*,*) "# csNph_L::      en.dat, dmx.dat  " 


C.............

      OPEN(2,FILE='dat/dmx.dat',FORM='UNFORMATTED')

      READ(2) MODE
      READ(2) LMIN, LMAX, NDIM, ANGULPARTINCLUDED

      DO L=LMIN,LMAX-1
	 LFIN = L + 1
	 DO N1 = 1, NDIM
	    DO N2 = 1, NDIM

	       READ(2) D1(N1,N2,LFIN)

	    ENDDO
	 ENDDO
      ENDDO

      CLOSE(2)

      WRITE(*,*) "# csNph_L::read_dmx       read dipoles(nl,ml+1)" 
C......................

      OPEN(2,FILE=ENFILE,FORM='UNFORMATTED')

      READ(2) LMIN1, LMAX1, NDIM1

      IF ((LMIN.NE.LMIN1).OR.(LMAX.NE.LMAX1).OR.(NDIM.NE.NDIM1)) THEN
	 WRITE(6,*)' CORRUPTED FILES: ',DIPFILE,' AND ',ENFILE,
     S             ' ARE NOT CONSISTENT'
	 WRITE(6,*)' LMIN,LMAX,NDIMENSION FOR BOTH FILES:'
	 WRITE(6,*)LMIN,LMAX,NDIM
	 WRITE(6,*)LMIN1,LMAX1,NDIM1
	 STOP

      ENDIF

      DO L = LMIN, LMAX
	 DO N1 = 1, NDIM

	    READ(2) EN(N1,L),RN(N1,L),phase(n1,L),an(n1,L)

	 ENDDO
      ENDDO

      CLOSE(2)

      WRITE(*,*) "# csNph_L::read_dmx   energies(n), phase(n), norm(n)" 
C.........................

      RETURN                                                            
      END                                                               
C
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE ANGULCOEFF(ANGULAR,POLAR,NBPHOT,NBPATH,ARBRELM,
     S                      MAXPATH,MAXPHOT,ANGULPARTINCLUDED)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ANGULPARTINCLUDED
      INTEGER POLAR,ARBRELM,CHEMIN
      DIMENSION ARBRELM(MAXPATH,MAXPHOT,2),ANGULAR(MAXPATH)
C
      DO CHEMIN=1,NBPATH
	 COEFF = 1.D0
	 DO NUM=1,NBPHOT
	    L1 = ARBRELM(CHEMIN,NUM,1)
	    M1 = ARBRELM(CHEMIN,NUM,2)
	    L2 = ARBRELM(CHEMIN,NUM+1,1)
	    M2 = ARBRELM(CHEMIN,NUM+1,2)
	    COEFF = COEFF * ANGUL(L1,M1,L2,M2,POLAR)
	 ENDDO
         IF ( ANGULPARTINCLUDED ) THEN
            ANGULAR(CHEMIN) = 1.D0
         ELSE   
	    ANGULAR(CHEMIN) = COEFF
         ENDIF
      ENDDO
C
      RETURN
      END
C                                                                       
C                                                                       
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE CHECK(LMIN,LMAX,NDIMENSION,NNINIT,NBPHOT,LINIT,MINIT)
C
      INTEGER ERROR
C
      ERROR = 0
      IF ((LINIT.LT.LMIN).OR.(LINIT.GT.LMAX)) ERROR = 1
      IF ((NNINIT.LT.1).OR.(NNINIT.GT.NDIMENSION)) ERROR = 2
      IF ((LINIT+NBPHOT.GT.LMAX).OR.(MAX0(LINIT-NBPHOT,0).LT.LMIN)) then
         
         write(*,*) '&  check::            LINIT = ', LINIT
         write(*,*) '&  check::           NBPHOT = ', NBPHOT
         write(*,*) '&  check::             LMIN = ', LMIN
         write(*,*) '&  check::             LMAX = ', LMAX

         ERROR = 3
      endif

      IF ((MINIT.LT.-LINIT).OR.(MINIT.GT.LINIT)) ERROR = 4
C
      IF (ERROR.NE.0) THEN
	 WRITE(6,*)'LACK OF STATES TO PERFORM MULTIPHOTON CALCULATION.'
	 WRITE(6,*)'ERROR #:',ERROR
	 STOP
      ENDIF
C
      RETURN
      END
C
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE FINDFINALSTATE(L,NN,E,ENER,NDIMENSION,MAXDIM,MAXL)     
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      LOGICAL FOUND                                                     
      DIMENSION ENER(MAXDIM,0:MAXL)                                     
C                                                                       
      NHI=NDIMENSION                                                    
      LO=1                                                              
      MID=NHI/2                                                         
      FOUND = .FALSE.                                                   
 1    IF (.NOT.FOUND) THEN                                              
         IF (E.GE.ENER(MID,L)) THEN                                     
              LO=MID                                                    
         ELSE                                                           
              NHI=MID                                                   
         ENDIF                                                          
         IF (NHI.EQ.LO+1) FOUND=.TRUE.                                  
         MID=LO+(NHI-LO)/2                                              
         GO TO 1                                                        
      ENDIF                                                             
      IF (DABS(ENER(LO,L)-E).GT.DABS(E-ENER(LO+1,L))) THEN              
         NN = LO + 1                                                    
      ELSE                                                              
         NN = LO                                                        
      ENDIF                    
                                         
      RETURN                                                            
      END                                                               
C
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C                                                                       
C CETTE ROUTINE CALCULE TOUS LES CHEMINS POSSIBLES ENTRE LES SITES      
C (L1,M1) ET (L2,M2) EN RESPECTANT LES RªGLES DE S{LECTIONS :           
C       L2 = L1 +- 1                                                    
C       M2 = M1 + MDELTA                                                
C EN PARTANT DU SITE INITIAL (LINIT,MINIT). IL Y A DANS CHAQUE CHEMIN   
C NBPHOT {TAPES. LE NOMBRE TOTAL DE CHEMINS TROUV{S EST STOCK{ DANS     
C NBPATH. POUR UN CHEMIN DONN{, APR{S UNE {TAPE DONN{E, LE SITE ACC{D{  
C EST DONN{ PAR :                                                       
C           L = ARBRELM(CHEMIN,NIVEAU,1)                                
C           M = ARBRELM(CHEMIN,NIVEAU,2)                                
C                                                                       
C                                                                       
      SUBROUTINE PATHLM(NBPHOT,LINIT,MINIT,MDELTA,NBPATH,ARBRELM)
      IMPLICIT INTEGER*4 (A-Z)                                          

      PARAMETER (MAXPATH =10,NIVMAX=10)
C             
      DIMENSION ARBRE(MAXPATH,NIVMAX)                                   
      DIMENSION ARBRELM(MAXPATH,NIVMAX,2)                               
C         
                                                              
      CALL PATHL(NBPHOT,LINIT,NBPATH,ARBRE)                             
C                                                                       
      CHEMIN = 0                                                        
      DO 1 J=1,NBPATH                                                   
         CHEMIN = CHEMIN + 1                                            
         DO 3 I=1,NBPHOT+1                                              
              L = ARBRE(J,I)                                            
              M = MINIT + MDELTA*(I-1)                                  
              IF ((M.GE.-L).AND.(M.LE.L)) THEN                          
                   ARBRELM(CHEMIN,I,1) = L                              
                   ARBRELM(CHEMIN,I,2) = M                              
                   GO TO 3                                              
              ELSE                                                      
                   CHEMIN = CHEMIN -1                                   
                   GO TO 1                                              
              ENDIF                                                     
 3       CONTINUE                                                       
 1    CONTINUE                                                          
      NBPATH = CHEMIN                                                   
      END                                                               
C                                                                       
C                                                                       
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE PATHL(NBPHOT,LINIT,NBPATH,ARBRE)
C                                                                       
C CETTE ROUTINE CALCULE TOUS LES CHEMINS POSSIBLES ENTRE LES SITES      
C (L1) ET (L2) EN RESPECTANT LES RªGLES DE S{LECTIONS :                 
C       L2 = L1 +- 1                                                    
C EN PARTANT DU SITE INITIAL (LINIT). IL Y A DANS CHAQUE CHEMIN         
C NBPHOT {TAPES. LE NOMBRE TOTAL DE CHEMINS TROUV{S EST STOCK{ DANS     
C NBPATH. POUR UN CHEMIN DONN{, APR{S UNE {TAPE DONN{E, LE SITE ACC{D{  
C EST DONN{ PAR :                                                       
C           L = ARBRE(CHEMIN,NIVEAU)                                    
C                                                                       
C                                                                       
C                                                                       
      IMPLICIT INTEGER*4 (A-Z)                                          
      PARAMETER (MAXPATH=10,NIVMAX=10,LMAX=2*NIVMAX)
      LOGICAL SITE(-1:LMAX,NIVMAX)                                      
      LOGICAL FINDPATH                                                  
      DIMENSION ARBRE(MAXPATH,NIVMAX)                                   
      COMMON /SENTIER/ NIV,SITE                                         
C                                                                       
      NIV = NBPHOT + 1                                                  
      CHEMIN = 0                                                        
C                                                                       
      CALL INITGRID(SITE,NIV,LINIT,NIVMAX,LMAX)                         
C                                                                       
 1    CHEMIN = CHEMIN + 1                                               
      IF (CHEMIN.GT.MAXPATH) THEN                                       
         WRITE(6,*)'TROP DE CHEMINS TROUV{S, AUGMENTER MAXPATH :',      
     S              CHEMIN                                              
      ENDIF                                                             
      IF (FINDPATH(ARBRE,CHEMIN,LINIT)) GO TO 1                         
C                                                                       
      NBPATH = CHEMIN                                                   
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE INITGRID(SITE,NIV,LINIT,NIVMAX,LMAX)                   
C                                                                       
C INITIALISATION DU TABLEAU SITE UTILIS{ POUR D{TERMINER LES CHEMINS    
C POSSIBLES.                                                            
C                                                                       
      IMPLICIT INTEGER*4 (A-Z)                                          
      LOGICAL SITE(-1:LMAX,NIVMAX)                                      
C                                                                       
      DO 1 N=1,NIV                                                      
         SITE(-1,N) = .FALSE.                                           
         DO 2 L=0,LINIT+NIV                                             
              SITE(L,N) = .TRUE.                                        
 2       CONTINUE                                                       
 1    CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      FUNCTION FINDPATH(ARBRE,CHEMIN,LINIT)                             
C                                                                       
C CHERCHE IT{RATIVEMENT UN CHEMIN COMPLET D'UN SITE INITIAL A UN SITE   
C FINAL, DIFF{RENT DE TOUS LES PR{CEDENTS. RENVOIE .TRUE. S'IL EN A     
C TROUV{ UN. LE CHEMIN AINSI TROUV{ EST AJOUT{ DANS LE TABLEAU : ARBRE. 
C                                                                       
C                                                                       
      IMPLICIT INTEGER*4 (A-Z)                                          
      PARAMETER (MAXPATH=10,NIVMAX=10,LMAX=2*NIVMAX)
      LOGICAL SITE(-1:LMAX,NIVMAX)                                      
      LOGICAL FINDPATH                                                  
      DIMENSION ARBRE(MAXPATH,NIVMAX)                                   
      COMMON /SENTIER/ NIV,SITE                                         
C                                                                       
C                                                                       
C TEST SI ON A ATTEINT LE DERNIER NIVEAU, AUQUEL CAS ON VIENT DE TROUVER
C UN CHEMIN.                                                            
C                                                                       
      L = LINIT                                                         
      N = 1                                                             
      ARBRE(CHEMIN,N) = L                                               
C                                                                       
 1    IF (N.EQ.NIV) THEN                                                
         ARBRE(CHEMIN,N) = L                                            
         SITE(L,N) = .FALSE.                                            
         FINDPATH = .TRUE.                                              
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      IF (SITE(L+1,N+1)) THEN                                           
         N = N+1                                                        
         L = L+1                                                        
         ARBRE(CHEMIN,N) = L                                            
         GO TO 1                                                        
      ENDIF                                                             
C                                                                       
      IF (SITE(L-1,N+1)) THEN                                           
         N = N+1                                                        
         L = L-1                                                        
         ARBRE(CHEMIN,N) = L                                            
         GO TO 1                                                        
      ENDIF                                                             
C LORSQUE L'ON ARRIVE ICI CECI SIGNIFIE QUE LE SITE (N,L) EST UN CUL-DE-
C SAC. IL FAUT DONC LIBERER LES SITES (N+1,L+1),(N+1,L-1) ET INTERDIRE  
C LE SITE (N,L).                                                        
C                                                                       
      SITE(L,N) = .FALSE.                                               
      SITE(L+1,N+1) = .TRUE.                                            
      IF (L.NE.0) SITE(L-1,N+1) = .TRUE.                                
C                                                                       
C TEST SI IL N'Y A PLUS DE CHEMIN, ALORS LE DERNIER CHEMIN ETAIT LE     
C PR{C{DENT.                                                            
C                                                                       
      IF (.NOT.SITE(LINIT,1)) THEN                                      
         FINDPATH = .FALSE.                                             
         CHEMIN = CHEMIN - 1                                            
         RETURN                                                         
      ENDIF                                                             
C                                                                       
C IL FAUT RECULER D'UN NIVEAU POUR POUVOIR REPARTIR                     
C                                                                       
      N = N-1                                                           
      L = ARBRE(CHEMIN,N)                                               
      GO TO 1                                                           
      END                                                               
C                                                                       
C                                                                       
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE SELECTPATH(NBPHOT,LFIN,MFIN,NBPATH,ARBRELM)
C                                                                       
C CETTE ROUTINE S{LECTIONNE DANS ARBRELM LES CHEMINS POUR LESQUELS L ET 
C M DU DERNIER NIVEAU (NIVEAU ATTEINT PAR LE DERNIER PHOTON) SONT TELS  
C QUE :                                                                 
C       L = LFIN                                                        
C       M = MFIN                                                        
C LES CHEMINS AINSI S{LECTIONN{S SONT RERANG{S DANS ARBRELM ET NBPATH   
C EST MODIFI{ (LES AUTRES CHEMINS SONT PERDUS).                         
C                                                                       
C                                                                       
      IMPLICIT INTEGER*4 (A-Z)                                          
      PARAMETER (MAXPATH=10,NIVMAX=10)
C                                                                       
      DIMENSION ARBRELM(MAXPATH,NIVMAX,2)                               
      DIMENSION ARBRE2(MAXPATH,NIVMAX,2)                                
C                                                                       
      NB=0                                                              
      DO I=1,NBPATH                                                     
         IF ((ARBRELM(I,NBPHOT+1,1).EQ.LFIN).AND.                       
     S       (ARBRELM(I,NBPHOT+1,2).EQ.MFIN)) THEN                      
              NB = NB+1                                                 
              DO J=1,NBPHOT+1                                           
                   ARBRE2(NB,J,1) = ARBRELM(I,J,1)                      
                   ARBRE2(NB,J,2) = ARBRELM(I,J,2)                      
              ENDDO                                                     
         ENDIF                                                          
      ENDDO                                                             
C                                                                       
      DO I=1,NB                                                         
         DO J=1,NBPHOT+1                                                
              ARBRELM(I,J,1) = ARBRE2(I,J,1)                            
              ARBRELM(I,J,2) = ARBRE2(I,J,2)                            
         ENDDO                                                          
      ENDDO                                                             
C                                                                       
      NBPATH = NB                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE COUNTFINALSTATE(NBPHOT,NBPATH,ARBRELM,TAB,NB,
     S                           MAXPATH,NIVMAX)
C                                                                       
C CETTE ROUTINE REL}VE DANS ARBRELM LES DIFF{RENTES VALEURS POSSIBLES   
C (L,M) DU DERNIER NIVEAU ET RENVOIE CES COUPLES DANS TAB. LE NOMBRE    
C DE CES COUPLES C EST NB.                                              
C                                                                       
C                                                                       
      IMPLICIT INTEGER*4 (A-Z)                                          
C                                                                       
      DIMENSION ARBRELM(MAXPATH,NIVMAX,2)                               
      DIMENSION TAB(NIVMAX,2)                                           
C                                                                       
      NB=1                                                              
      TAB(NB,1) = ARBRELM(1,NBPHOT+1,1)                                 
      TAB(NB,2) = ARBRELM(1,NBPHOT+1,2)                                 
      DO I=2,NBPATH                                                     
         L = ARBRELM(I,NBPHOT+1,1)                                      
         M = ARBRELM(I,NBPHOT+1,2)                                      
         DO J=1,NB                                                      
            IF ((L.EQ.TAB(J,1)).AND.(M.EQ.TAB(J,2))) THEN               
               GOTO 1                                                   
            ELSE                                                        
               IF(J.EQ.NB) THEN                                         
                  NB = NB+1                                             
                  TAB(NB,1) = L                                         
                  TAB(NB,2) = M                                         
                  GOTO 1                                                
               ENDIF                                                    
            ENDIF                                                       
         ENDDO                                                          
 1       CONTINUE                                                       
      ENDDO                                                             
C                                                                       
      RETURN                                                            
      END                                                               
C
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      FUNCTION ANGUL(L1,M1,L,M,POL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POL
      ANGULAIRE=0.D0
C
      IF (POL.EQ.0) THEN

	 A=DSQRT(DBLE((L+M+1)*(L-M+1))/DBLE((2*L+1)*(2*L+3)))
	 B=DSQRT(DBLE((L+M)*(L-M))/DBLE((2*L+1)*(2*L-1)))

	 IF ((M1.EQ.M).AND.(L1.EQ.L+1)) ANGULAIRE = A
	 IF ((M1.EQ.M).AND.(L1.EQ.L-1)) ANGULAIRE = B

      ENDIF
C
      IF (POL.EQ.-1) THEN

	 A=-DSQRT(DBLE((L+M+1)*(L+M+2))/DBLE((2*L+1)*(2*L+3)))
	 B=DSQRT(DBLE((L-M)*(L-M-1))/DBLE((2*L+1)*(2*L-1)))

	 IF ((M1.EQ.M+1).AND.(L1.EQ.L+1)) ANGULAIRE = A
	 IF ((M1.EQ.M+1).AND.(L1.EQ.L-1)) ANGULAIRE = B

      ENDIF
C
      IF (POL.EQ.1) THEN
	 A=DSQRT(DBLE((L-M+1)*(L-M+2))/DBLE((2*L+1)*(2*L+3)))
	 B=-DSQRT(DBLE((L+M)*(L+M-1))/DBLE((2*L+1)*(2*L-1)))
	 IF ((M1.EQ.M-1).AND.(L1.EQ.L+1)) ANGULAIRE = A
	 IF ((M1.EQ.M-1).AND.(L1.EQ.L-1)) ANGULAIRE = B
      ENDIF
C
      ANGUL = ANGULAIRE
C
      RETURN
      END
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE CSPHFILE(N, L, MODE)
      INTEGER L,MODE,N
      
      IF(L.GT.5) THEN

         WRITE(*,*) ' L MUST BE LESS OR EQUAL 8  . FOR L > 8', 
     1        ' MODIFY SUBROUTINE CS2PHFILE. EXITING...'
         STOP

      ELSE

C         WRITE(16,*) 'CS2PH :: DATA FOR ', L1, ' -> ', L2

      ENDIF

      IF(MODE.EQ.0) THEN

         IF(L.EQ.0) OPEN(N,FILE='dat/csNph-0.v.dat')
         IF(L.EQ.1) OPEN(N,FILE='dat/csNph-1.v.dat')
         IF(L.EQ.2) OPEN(N,FILE='dat/csNph-2.v.dat')
         IF(L.EQ.3) OPEN(N,FILE='dat/csNph-3.v.dat')
         IF(L.EQ.4) OPEN(N,FILE='dat/csNph-4.v.dat')
         IF(L.EQ.5) OPEN(N,FILE='dat/csNph-5.v.dat')
         
      ELSE 

         IF(L.EQ.0) OPEN(N,FILE='dat/csNph-0.l.dat')
         IF(L.EQ.1) OPEN(N,FILE='dat/csNph-1.l.dat')
         IF(L.EQ.2) OPEN(N,FILE='dat/csNph-2.l.dat')
         IF(L.EQ.3) OPEN(N,FILE='dat/csNph-3.l.dat')
         IF(L.EQ.4) OPEN(N,FILE='dat/csNph-4.l.dat')
         IF(L.EQ.5) OPEN(N,FILE='dat/csNph-5.l.dat')
         
      ENDIF

      RETURN
      END
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C  (C) Copr. 1986-92 Numerical Recipes Software #12-)L.
      SUBROUTINE ratint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n),TINY
      PARAMETER (NMAX=100,TINY=1.d-25)
      INTEGER i,m,ns
      DOUBLE PRECISION dd,h,hh,t,w,c(NMAX),d(NMAX)
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.d0)then
          y=ya(i)
          dy=0.0d0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.d0) then              
          write(*,*) 'failure in ratint dd = 0.'
             stop
          endif
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C#EOF
