!     
      PROGRAM  DMX1E
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.1e.inc"
C      PARAMETER(NP = 2000, NL = 15)
      DIMENSION R(NP),DR(NP),ni(NL),nf(NL)
      character*16 gauge
      CHARACTER*100 ARGV
!     
 577  format(5A14)
 566  format(5A16)
!

      CALL GETARG(1, ARGV)
      READ(ARGV,*) GAUGE

      open(9,file='inp/dmx1e.inp',status='old')
      read(9,*) lmin, lmax
      read(9,*) (ni(i), i = 1, lmax - lmin + 1 )
      read(9,*) (nf(i), i = 1, lmax - lmin + 1 )
      close(9)
C.....................................................................


      write(*,*)  lmin, lmax, ' dipoles in ', gauge, 'gauge'

C*
C*    Calculate   DME(njj-1,njj) up to DME(nj-1,nj)
C*    in total    |njj - nj|  DME


      OPEN(16,FILE ='out/dmx1e.out')


      WRITE(16,'(a10)') "# dmx1e:: "
      WRITE(16,'(a10)') "# dmx1e:: grid parameters:"
      CALL RIN(R,DR,H,NO,idr)

      IF(GAUGE.EQ.'v') THEN
         WRITE(16,'(a10)') " Velocity:"
      ELSE
         WRITE(16,'(a10)') " Length:"
      ENDIF


      if((lmin-1).gt.0) then
         lmin = lmin -1 
      endif
      
      do  l = lmin, lmax - 1

        mode = 1 
        if(gauge.eq.'v')         mode = 0 

      write(16,'(a4,1X,I2,a4,I2,a4)') 'd(n', l,' ->m', l+1, ')'


      call dmx1efile(3, l, mode)
      call dmx(gauge,r,dr,z,h, l,l+1, nf(l-lmin+1), nf(l-lmin+2),no)
      close(3)

      enddo

      stop
      end
      
C*********************************************************************
      SUBROUTINE DMX(gauge,r,dr,z,h,lb,lz,nwfb,nwfz,no)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.1e.inc"
C      PARAMETER(NS = 802, NP = 2000 )
      DIMENSION PZ(NS-2,NP),DVX(NS-2,NS-2),FF(NP)
      dimension PB(NP),PD(NP),DR(NP),R(NP)
      character*16 gauge
C...........
      call wf1efile(2, lz) 

      do  i = 1, nwfz
         do   j = 1, no
            pz(i,j) = 0.0d+00
         enddo
      enddo


      DO  NC = 1, NWFZ
         READ(2) (FF(J), J = 1, no)
         DO  J = 1, NO
           PZ(NC,J) = FF(J)         
        enddo
      enddo

         CLOSE(2)

         
         call wf1efile(2,lb) 

         DO 200 NR = 1, NWFB

           DO  J = 1, NO
             PB(J) = 0.0D+00
           ENDDO
          
           READ(2) ( PB(J), J = 1, NO )


           CALL DERIV(PB, PD, R, DR, Z, H, LB, LZ, NO)

C*********************************************************************
C*     dipole matrix elements : 
C*           velocity gauge (mode = 0) ,
C*           length   gauge (mode = 1)
C*********************************************************************       


      DO 200 NC = 1, NWFZ

         IF(gauge.EQ.'v') THEN
 
            MODE = 0 

            DO  J = 1, NO
               FF(J) = PZ(NC,J) * PD(J)
            ENDDO   
            
            DELR  = 0.0D+00

            IF(R(1) .NE. 0.0D0) DELR = 0.50D+00 * FF(1) * DR(1)

            DVX(NR,NC) = RINT(FF, 1, NO, 14, H ) + DELR
            
         ELSE

            MODE = 1

            FF(1) = 0.0D+00

            IF(DR(1).NE.0.0D0) FF(1) = PZ(NC, 1) * PB(1) * R(1)/DR(1)

            DELR = 0.0D+00

            IF(R(1).NE.0.0D0) DELR = 0.50D+00 * FF(1) * R(1)

            DO  J = 2, NO
               
               FF(J) = PZ(NC,J) * PB(J) * R(J)**2 / DR(J)

            ENDDO

            DVX(NR,NC) = RINT(FF, 1, NO, 14, H ) + DELR

         ENDIF


 200     CONTINUE


         CLOSE(2)
        

C***     write dme in oupout files.
C...................................


C
C     l-->l+1 (lb=l, lz=lb+1)
C
         WRITE(3) mode
         WRITE(3) LB,LZ,NWFB,NWFZ
         WRITE(3) ((DVX(NR,NC), NC = 1, NWFZ), NR = 1, NWFB)

C.....................................

         WRITE(16,'(a5,1x,6a14)') " n\m","1","2","3","4","5","6"
         DO  NR = 1, 6
            WRITE(16,2) NR, ( DVX( NR,NC ), NC = 1, 6)
        enddo

C-------------------
    1 FORMAT(/2X,'MATRIX <PZ,O,PB>; LB, NB(ROW) & LZ, NZ(COL.)=',4I3)
    2 FORMAT(2X,I3,1X,1P8E16.8)
 17   format(4i5)
 16   format(1p6d21.13)


C------------------
 
      RETURN
      END
C#####################################################################
      SUBROUTINE DERIV(PB,F1,R,DR,Z,H,LB,LZ,NO)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "parameter.1e.inc"
C      PARAMETER (NP = 2000)
      DIMENSION PB(3),F1(1),R(NP),DR(NP)

C**   PB  --  INPUT FUNCTION  ( ANGULAR MOMENTUM = LB) ?
C**   F1  --  OUTPUT FUNCTION ( ANGULAR MOMENTUM = LZ) ?
      
      Z = 1.0D+00 
      LL  = LB + 1
      LB3 = LB + 3

      ZL = Z / DFLOAT(LL)

      D = ((1.0D0-ZL*R(2))*R(3)**2-(1.0D0-ZL*R(3))*R(2)**2)

      A = ( PB(2) * R(3)**2/R(2)**LL - PB(3)* R(2)**2/R(3)**LL )/D

      B = (PB(3)*(1.0D0-ZL*R(2))/R(3)**LL-PB(2)*(1.0D0-ZL*R(3))
     1  /R(2)**LL)/D

      DO 10 J = 1, 3
   10 F1(J) = A * (DFLOAT(LL) - DFLOAT(LB+2)*ZL*R(J))*R(J)**LL+
     1        B * DFLOAT(LB3) * R(J)**LB3

      NN = NO - 3

      DO 12 J = 4, NN
   12 F1(J)=((PB(J+3)-PB(J-3))/60.D0-3.D0*(PB(J+2)-PB(J-2))/20.D0+
     1      3.D0*(PB(J+1)-PB(J-1))/4.D0)/H

      F1(NO-2) = ( 2.D+00 * ( PB(NO-1) - PB(NO-3) )/ 3.0D+00 
     1           - ( PB(NO) - PB(NO-4) ) / 12.0D+00 )/ H

      F1(NO-1) = 0.5D+00 * ( PB(NO) - PB(NO-2) ) / H
      F1(NO)   = (PB(NO)-PB(NO-1))/H
      F1(1)    = 0.0D+00

      IF (LB .GT. LZ) GO TO 4
      IF(DR(1).NE.0.0D0) F1(1)=F1(1)-DFLOAT(LL)*PB(1)/DR(1)

      DO 3 J=2,NO
    3 F1(J)=F1(J)-DFLOAT(LL)*PB(J)/DR(J)
      RETURN

    4 IF(DR(1).NE.0.0D0) F1(1)=F1(1)+DFLOAT(LB)*PB(1)/DR(1)

      DO 5 J = 2, NO
    5 F1(J)  = F1(J) + DFLOAT(LB) * PB(J) / DR(J)

      RETURN
      END
C###################################################################
C#EOF

C*******************************************************************
C*
C*    Purpose      ::  Calculates One-electron dme  with B-splines 
C*                     technique
C*    Authors      ::  X. Tang (1985-1990), Jian Zhang (1993-1996)
C*    Libs         ::  None
C*    Dependecies  ::  rinbat.f 
C*    Language     ::  f77, g77
C*    Compilers    ::  xlf IBM Risc 6000, g77 on Redhat linux
C*    Modified     ::  L.A.A. Nikolopoulos  (1999),
C*                     L.A.A. Nikolopoulos  (Feb/2001/cineca) 
C*                     to calculate in velocity or length dme       
C*    Compiled as  ::
C*                   IBM/RISC 6000 :
C*                   xlf -c O3      aedmxn.f rinbat.f 
C*                   xlf -o Raedmx  aedmxn.o rinbat.o    
C*
C*
C*   L. AA. N.  22/09/1999 iesl days
C*
C*
C*    Z         ==  Atomic Number 
C*    R(i)      ==  knot sequence
C*    ni        ==  Initial Symmetry  ni = 0, 1, ....., LF
C*    nf        ==  Finaly Symmetry   nf = 0, 1, ....., LF
C*    lmin,lmax ==  Calculate  dme
C*    From dme(lmin-1 -> lmin) up to  dme(lmax-1 -> lmax)
C*
C*    Notes::
C*
C*    -1-
C*    Z  set to 1 
C*    
C*    DME in length and velocity gauge are independent on Z
C*    therefore the value of Z should be set to 1
C*    This happens because in velocity gauge is needed 
C*    the derivative of the WF  dy/dx and the formula which 
C*    is used for the first 3 - points contains Z (see DERIV subroutine).
C*    However the final result for the velocity gauge (DME_v) doesn't depend on Z.
C*
C*   -2- 
C*    Generalize to include Length gauge
C*    ..made it soon ...
C*
C*   -3- 
C*    Accuracy is of eight (8) digits  (why is this happens?)
C*
C*   -4-
C*    Change it with my own C++ code 
C*
C*    READ  input file "dmx1e.inp" 
C*    indat
C*    outdat
C*    outdat
C*eof



