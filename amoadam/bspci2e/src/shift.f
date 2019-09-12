C******************************************************************************
C   This programm calculate AC-Stark shift of states bound or continuum.
C
C
C   This Program calculate AC Stark Shift of bound states
C
C                                               Jian Zhang , May 1991
C
C   modified and extended to account for continuum states
C
C                                              L. A.A. Nikolopoulos 2001
C 
C   note 1 :
C   For continuum states the interest is on the shif of the 
C   autoinizing states compared with the smooth continuum 
C   (ponderomotive shift).
C
C
C
C******************************************************************************

      PROGRAM ACSHIFT
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION AM(500,500,5),N(5),EN(1100,5),X1(500),NNZ(500)

      factor = (1.D0/( 4*350.9D0 )) * 27.211396D+00 


C#.............................................................................
C#
C#         INPUT L(initial), L(inter) N(initial), N(inter) 
C#            and the ENERGY of the photon in(a.u) 
C#
C#.............................................................................

      OPEN(18,FILE='inp/shift.inp')
      READ(18,*) la, na, lm, omega
      READ(18,*) Eth2
      READ(18,*) id

      CLOSE(18)

C#
C#      INPUT THE ENERGY LEVELS OF L = la and L = lm 
C#      and the dipole matrix in between
C#

c      CALL OPENFL1(17,la,lm,0)


c      OPEN(17,file='dmx2e-s01.dat')

      open(17,file ='dmx2e-s01.dat', 
     1     form ='unformatted', access='sequential')


      if(id.eq.0) then

      READ(17)   MODE
      READ(17) LB, LZ, N(la+1), N(lm+1)

      if(lb.ne.la.or.lz.ne.lm) then 
         write(*,*),' Wrong file! ' 
         stop
      endif 

      READ(17) ( EN(K, La + 1), K = 1, N( la + 1) )
      READ(17) ( EN(K, Lm + 1), K = 1, N( lm + 1) )

      DO  NR = 1, N( la + 1 )

         READ(17) (AM(nr,nc,lm + 1), NC = 1, N(lm+1))


      ENDDO

      ELSE

      READ(17)   MODE
      READ(17) LZ, LB, N(lm+1), N(la+1)

      write(*,*) 'lb,lz = ',  lb,lz
      write(*,*) 'la,lm = ' , la,lm
      if(lb.ne.la.or.lz.ne.lm) then 
         write(*,*),' Wrong file! ' 
         stop
      endif 

      READ(17) ( EN(K, Lm + 1), K = 1, N( lm + 1) )
      READ(17) ( EN(K, La + 1), K = 1, N( la + 1) )

      DO  nr = 1, N( lm + 1 )

         READ(17) ( AM(nc, nr, lm + 1), NC = 1, N(la + 1))


      ENDDO


      ENDIF
      CLOSE(17)


C.......................................................................

      WRITE(*,*)' W_PH    = ', omega,  ' a.u. '
      WRITE(*,*)' Eth2    = ', Eth2,   ' Ryd  '
      WRITE(*,*)' NA      = ', na 
      WRITE(*,*)' LA      = ', la
      WRITE(*,*)' Lm      = ', lm
      WRITE(*,*)' LB      = ', lb
      WRITE(*,*)' LZ      = ', lz
      WRITE(*,*)' N(la+1) = ', N(la+1)
      WRITE(*,*)' N(lm+1) = ', N(lm+1)

C#
C#  Renormalize the energy according to the appropriate threshold  ******
C#

      DO  k = la + 1, lm + 1
         DO  K1 = 1, N(K)

            EN(K1,K) = ( EN(K1,K) + Eth2) / 2.0D+00

         ENDDO
      ENDDO

C......................................................................


      WRITE(*,*)' E(g)  = ', EN(1,1), ' a.u. '

C      WRITE(*,*)' E(g)  = ', EN(1,1) * 109737.31d0 * 2.d0, ' cm(-1)'

C.......................................................................

C#
C#    If the state is in the continuum calculate the renormalization 
C#    factor from the density of states
C#  

      if(en(na, la+1).gt.0.0D0) then 

       rnorm = sqrt( 2.0D+00/dabs( en(na+1, la+1) - en(na-1,la+1) ) )

      else

         rnorm = 1.0D+00

      endif

       WRITE(*,*) ' EN(SHIFT) = ' , en(na,la+1), ' a.u. '
       WRITE(*,*) ' RNORM     = ' , rnorm       , ' a.u. ' 

C.......................................................................
C#
C#                   At this point the sum rule 
C#
C#                                  ----  
C#                  Sum( n = 0 ) =  \     w_an  * | < a | d | n > |^2 
C#                                  /___
C#
C#               is evaluated for chechking basis reliability
C#
C.......................................................................
          
          
          sum0 = 0.D+00
          sum0_pos = 0.0D+00
          sum0_neg = 0.0D+00

          do  k1 = 1, n(lm + 1)

             
             de_am = en(k1, lm + 1) - en(na, la + 1) 

             sum0 = de_am * am(na, k1, lm + 1) * am(na, k1, lm + 1)

C#
C#       add negative and positive values separately
C#
            if(sum0.gt.0.0D0)  then

               sum0_pos  = sum0_pos  + sum0
            else

               sum0_neg  = sum0_neg  + sum0
            endif


         enddo

         sum_total = sum0_pos + sum0_neg 
         sum_norm = (sum0_pos + sum0_neg ) * rnorm**2


         WRITE(*,*), 'Sum rule S(0) = ', sum_total,sum_norm


C..................................................................
C#
C#          SUMMATION OF ALL THE lm INTERMEDIATE STATES
C#               FOR THE | na, L = la >  state
C#                      DOWN ONE PHOTON
C#
C#
C#              ---------------------------- |na,la>
C#                                   \ 
C#                                    \   --------------------------
C#                                    _\| --------------------------
C#                                 ----------------------------------
C#




       omega0 = en(na,la+1) - omega

       WRITE(*,*),' E_down = ', omega0, ' a.u. '


       if(omega0.le.0.d0.or.omega0.lt.en(1,lm+1)) then

          
          sum1 = 0.D+00
          acmd1_pos = 0.0D+00
          acmd1_neg = 0.0D+00
          

          do  k1 = 1, n(lm + 1)

             
             de_pos = ( en(na, la + 1) - omega - en(k1, lm + 1) )


             sum1 = am(na, k1, lm + 1) * am(na, k1, lm + 1)/de_pos
           

C#
C#       add negative and positive values separately
C#

            if(sum1.gt.0.0D0)  then

               acmd1_pos = acmd1_pos + sum1 
            else

               acmd1_neg = acmd1_neg + sum1 
            endif


         enddo


C#
C#    down one-photon contribution, no pole
C#

         acmd1 =  acmd1_pos + acmd1_neg 


         WRITE(*,*), ' M(Lm, -, no pole) = ', acmd1 * factor




      else


C#
C#   FOR CONTINUUM INTERMEDIATE STATES L = lm 
C#

         do  k1 = 1, n(lm+1)

            x1(k1)  = dabs( omega0 - en(k1,lm+1) )
         enddo

C#         
C#    LOOKING FOR THE NEAREST STATE  
C#

         call indexx( n(lm+1), x1, nnz)

         nzd = nnz(1)
         
        if( dabs(omega0 - en(nzd, lm+1)).gt.9.0d-03 ) goto 410

         amd1     = 0.0D+00
         amd1_pos = 0.0D+00
         amd1_neg = 0.0D+00
         sum3     = 0.0D00

         omad = en(nzd,lm+1) - en(na, la+1)


         do  k1 = 1, nzd - 1

            sum3  =   am(na,k1,lm+1) * am(na,k1,lm+1)
     1            /  (en(na,la+1) + omad - en(k1,lm+1))


C# 
C#       add negative and positive values separately
C#

            if(sum3.gt.0.0D0)  then

            amd1_pos = amd1_pos + sum3
            else

            amd1_neg = amd1_neg + sum3
            endif

         enddo

C#
C#       omitt the pole  k1 = nz1
C#

         do k1 = nzd + 1, n( lm + 1)

           sum3 =   am(na,k1,lm+1) * am(na,k1,lm+1) 
     1            /  ( en(na,la+1) + omad - en(k1,lm+1))


            if(sum3.gt.0.0D0)  then

            amd1_pos = amd1_pos + sum3
            else

            amd1_neg = amd1_neg + sum3
            endif

         enddo

C#
C#    down one-photon contribution, pole
C#

         acmd1 =  amd1_pos + amd1_neg 


         WRITE(*,*),' M(Lm, -, pole)  = ', acmd1 * factor


      endif



C#
C#
C#         SUMMATION OF ALL THE la + 1 INTERMEDIATE STATES
C#                FOR THE | na, L = la > state
C#                      UP ONE PHOTON
C#
C#

      omega1 = en(na,la+1) + omega


      WRITE(*,*),' E_up = ', omega1 , ' a.u.'

C#
C#   w_1 < 0.0 , the intermediate states are bound 
C#


      if(omega1.le.0.d0.or.omega1.gt.en(n(lm+1),lm+1)) then


         acmu1 = 0.D+00
         sum2  = 0.D+00
         acmu1_pos = 0.0D+00
         acmu1_neg = 0.0D+00

         do  k1 = 1, n( lm + 1)

            de_neg = (en(na,la+1) + omega - en(k1, lm+1))

            sum2 = am(na,k1,lm+1) * am(na,k1,lm+1)/de_neg


C# 
C#       add negative and positive values separately
C#

            if(sum2.gt.0.0D0)  then

            acmu1_pos = acmu1_pos + sum2 
            else

            acmu1_neg = acmu1_neg + sum2 
            endif



         enddo


         acmu1 =  acmu1_pos + acmu1_neg 

         WRITE(*,*),' M(Lm, +, no pole)  = ', acmu1 * factor


      else


C#
C#   FOR CONTINUUM INTERMEDIATE STATES L = lm 
C#

         do  k1 = 1, n(lm+1)

            x1(k1)  = dabs( omega1 - en(k1,lm+1) )
         enddo

C#         
C#    LOOKING FOR THE NEAREST STATE  
C#

         call indexx( n(lm+1), x1, nnz)

         nzu = nnz(1)
         
        if( dabs(omega1 - en(nzu, lm+1)).gt.9.0d-03 ) goto 410

         amu1     = 0.0D+00
         amu1_pos = 0.0D+00
         amu1_neg = 0.0D+00
         sum3     = 0.0D00

         omau = en(nzu,lm+1) - en(na, la+1)


         do  k1 = 1, nzu - 1

            sum3  =   am(na,k1,lm+1) * am(na,k1,lm+1)
     1            /  (en(na,la+1) + omau - en(k1,lm+1))


C# 
C#       add negative and positive values separately
C#

            if(sum3.gt.0.0D0)  then

            amu1_pos = amu1_pos + sum3
            else

            amu1_neg = amu1_neg + sum3
            endif

         enddo

C#
C#       omitt the pole  k1 = nz1
C#

         do k1 = nzu + 1, n( lm + 1)

           sum3 =   am(na,k1,lm+1) * am(na,k1,lm+1) 
     1            /  ( en(na,la+1) + omau - en(k1,lm+1))


            if(sum3.gt.0.0D0)  then

            amu1_pos = amu1_pos + sum3
            else

            amu1_neg = amu1_neg + sum3
            endif

         enddo

         acmu1 =  amu1_pos + amu1_neg 


         WRITE(*,*),' M(Lm, +, pole)  = ', acmu1 * factor

c         WRITE(*,801), oma1, omega
c         WRITE(* 802), acmu1

      ENDIF

c1.10244781D0       amc = ( acmd1 + acmu1 ) * 0.0156235d0

       amc = ( acmd1 + acmu1 ) 
C* rnorm
**2 


c/0.80655
         WRITE(*,800)  amc * factor*8.065479D0
         WRITE(*, 805) amc * factor 

      STOP

 410     print*,'There is no state at the pole'
 801     format('  Omega 1,2,3(LA+1)=',3d15.6)
 802     format('  M    1,2,3(LA+1)=',3d15.6)
 803     format('  Omega 1,2,3(LA-1)=',3d15.6)
 804     format('  M    1,2,3(LA-1)=',3d15.6)
c 800  FORMAT(/'AC Stark shift:',D15.7,' cm-1 / ', 
 800  FORMAT(/'AC Stark shift:  ',F16.5,' cm-1 / ', 
     F                                       '10e10 W.cm-2'/)
 805  FORMAT(/' or             ',F16.5,' eV / ',
     2                                       '10e14 W.cm-2'/)    
 908  FORMAT(1x,6E14.6)
 988  FORMAT(2X,I4,I7,5E14.6)
 900  FORMAT(2X,I7,6(D10.3,2X))
 400  FORMAT(1X,'PULSE DURATION(A.U)=',D15.8)
 300  FORMAT(1X,'PEAK INTENSITY(A.U)',D15.8)
 999  STOP
      END
C###########################################################################
      SUBROUTINE OPENFL1(N,L1,L2,LS)
      IF(LS.EQ.0) THEN
      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='OSSP1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='OSPS1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='OSPD1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='OSDP1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='OSDF1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='OSFD1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.4) OPEN(N,FILE='OSFG1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3) OPEN(N,FILE='OSGF1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.5) OPEN(N,FILE='OSGH1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4) OPEN(N,FILE='OSHG1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.6) OPEN(N,FILE='OSHI1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5) OPEN(N,FILE='OSIH1.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      ENDIF
      IF(LS.EQ.1) THEN
      IF(L1.EQ.0.AND.L2.EQ.1) OPEN(N,FILE='OSSP3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.0) OPEN(N,FILE='OSPS3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.1.AND.L2.EQ.2) OPEN(N,FILE='OSPD3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.1) OPEN(N,FILE='OSDP3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.2.AND.L2.EQ.3) OPEN(N,FILE='OSDF3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.2) OPEN(N,FILE='OSFD3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.3.AND.L2.EQ.4) OPEN(N,FILE='OSFG3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.3) OPEN(N,FILE='OSGF3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.4.AND.L2.EQ.5) OPEN(N,FILE='OSGH3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.4) OPEN(N,FILE='OSHG3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.5.AND.L2.EQ.6) OPEN(N,FILE='OSGI3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IF(L1.EQ.6.AND.L2.EQ.5) OPEN(N,FILE='OSIH3.DAT',
     1    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      ENDIF
      RETURN
      END
C###########################################################################
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRIN(1),INDX(1)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
C###########################################################################
C#EOF
