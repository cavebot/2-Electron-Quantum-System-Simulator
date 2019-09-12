C#######################################################################
C# This program calculates the AC-STARK SHIFT for a two-electron state
C# in L ang V gauges, as it given by the formula (49) page 241 in
C# in P.Lambropoulos et al Physics Reports, 305, (1998), 203-293.
C#                       _____
C#                 E*E    \        2 * w_in | < n | D | i > |^2
C#    S_( W ) == - ----    \    _________________________________
C#                  4     /       
C#                       /____      W_in * W_in  -  W * W 
C#                         i   
C#
C#          W_in  == W_i - W_n 
C#
C#######################################################################
C#    :  infile:  contains energy levels and dipole moments            #
C#    : outfile:  photon energies and ionization cross-section         #
C#    :    NINI:  INDEX of the INITIAL STATE                           #
C#    :     NPT:  number of points for output                          #
C#######################################################################

      PROGRAM ACSTARK

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(NHMX=500)
      DIMENSION  DMX(NHMX,NHMX,2), ENI(NHMX),ENF(NHMX),ENM(NHMX)
      character*16 idmxfile1, idmxfile2, choice

C#
C#    open input file (cs2ph.inp) and read data
C#

      open(9,file='inp/stark.inp',status='old')
      read(9,*) idmxfile1, idmxfile2
      read(9,*) threshold, w1, w2
      read(9,*) nstep, nnf
      read(9,*) id
      read(9,*) id1

C#######################################################################
C#
C#                       open outfile
C#


      open(16,file ='out/stark.out') 

      write(*,*),'read input file', dmxfile1


C#######################################################################
C#
C#                 read for first photon
C#                            hw
C#      |INITIAL STATES > ---------> |INTERMEDIATE STATES >
C#
C#       read energies and 2e - dme      id = 0 ?
C#
C#


      open(17,file=idmxfile1, form='unformatted', access='sequential')


      if( id.eq.0 ) then

C#
C#           dmx_1( n L , m L + 1) 
C#       

         read(17)  mode  
         read(17)  loi, lom, nsi, nsm
         read(17) (eni(nsc), nsc = 1, nsi)
         read(17) (enm(nsr), nsr = 1, nsm)

         write(*,*) ' read energies'
         write(*,*) ' nsi = ', nsi
         write(*,*) ' nsm = ', nsm

         do nsr = 1, nsi

            read(17) ( dmx(nsc, nsr, 1), nsc = 1, nsm)
         enddo

         write(*,*),' Reading first transition:', loi, '---->', lom

      else

C#
C#               dmx_1( m L + 1 , n L ) 
C#       

         read(17)  mode  
         read(17)  lom, loi, nsm, nsi
         read(17) (enm(nsc), nsc = 1, nsm)
         read(17) (eni(nsr), nsr = 1, nsi)

         do nsc = 1, nsm
          
            read(17) ( dmx( nsc, nsr, 1), nsr = 1, nsi)
         enddo

      write(*,*),' Reading first transition:', loi, '---->', lom

      end if 



      close(17)



C#######################################################################
c***      Read second photon    DMX(m,f)
c***                                    hw
c***      | INTERMEDIATE STATES > ----------------------> | FINAL STATES >

      open(17,file = idmxfile2, form ='unformatted',access='sequential')

C#
C#      id1 = 0
C#
      if( id1.eq.0 ) then

C#
C#       dmx_2( m L+1, k L + 2 ) 
C#


         read(17)  mode
         read(17)  lom, lof, nsm, nsf
         read(17) (enm(nsc), nsc = 1, nsm)
         read(17) (enf(nsr), nsr = 1, nsf)

         do nsc = 1, nsm

            read(17) ( dmx(nsc, nsr, 2), nsr = 1, nsf)
         enddo

      write(*,*),' Reading second transition:', lom, '---->', lof


      else
C#
C#       dmx_2( k L + 2, m L + 1) 
C#

         read(17)  mode
         read(17)  lof, lom, nsf, nsm
         read(17) (enf(nsc), nsc = 1, nsf)
         read(17) (enm(nsr), nsr = 1, nsm)

         do nsr = 1, nsf

            read(17) ( dmx( nsc, nsr, 2), nsc = 1, nsm)

         enddo

      write(*,*),' Reading second transition:', lom, '---->', lof

      endif

      close(17)
C#######################################################################


      WRITE(*,*)' W_PH    = ', w1 ,  ' a.u. '
      WRITE(*,*)' Eth2    = ', threshold,   ' Ryd '
      WRITE(*,*)' NA      = ', nnf 


      en_ryd  = 13.605D+00

      write(*,*) ' Ionization Potential of A+  E_A(eV) = ', threshold 
     1 * en_ryd
      write(*,*) ' Range of photon energies = ', 
     1 w1 * 27.211396, w2 * 27.211396
      write(*,*) ' En(eV)  = ', eni(nnf) * en_ryd, eni(nnf) * en_ryd 
      write(*,*) ' NSM = ',  nsm

      write(*,*)  dmx(1,1,1), dmx(2,1,1)

      pi = 3.141592654D+00

 
C#
C#     Find the normalization factor of the shifted state enf(nnf)
C#
C#                 dens_j = sqr(2/DE(nnf)) 


         if(enf(nnf).lt.threshold)    then 


C#
C#      if shifted state belongs to bound spectrum no renormalization 
C#      factor is needed.
C#
C#                    dens_j  = 1 ;

            rho_j = 1.0D+00 ;

            write(*,*) ' Shifted state is bound state. En = ',enf(nnf)
         else

            rho_j = 1.0D0/sqrt(dabs( enf(nnf + 1) - enf(nnf - 1) )/4)

            write(*,*) ' Shifted state is continuum ', enf(nnf)

         endif


C#######################################################################
C#
C#   For photon frequencies in (w1, w2) calculate (-) dipole 
C#   polarizability  a_w(-)  
C#
C#              calculate   negative part ad_2 
C#
C#
C#                  _____
C#                  \      | D (i,n ) * rho(n) |^2 
C#    a_w(-) ==      \     ________________________
C#                   /       
C#                  /____     W_in   -  W 
C#                    i 
C#
C#


      step = abs(w2 - w1)/nstep 
 

c     do eph = w1, w2, step



      eph = w1

      write(*,*) ' w1 = ', w1

C#
C#      First look if there is pole and if any change accorrdingly 
C#      the photon energy to match exactly to the pole state.
C#

         de1 = eni(nnf)/2 + eph


         do im = 1, nsm
            

            deu   = ( enm(im + 1) + enm(im) ) / 4 
            deb = ( enm(im) + enm(im - 1) ) / 4 
            
            if(de1.gt.deb.and.de1.le.deu.and.enm(im).gt.threshold) then
               
               ewph = ( enm(im) - eni(nnf) )/2 


c               ewph = eph

c               write(*,*) ' w2 = ', ewh
               ipole = im 
               write(*,*) ' Pole : E(',im,' ) = ', enm(ipole)

               goto 90


            else


               ewph = eph

            endif


         enddo

 90      continue


C#
C#       Sum over the intermediate states 
C#
C#

        ad_1        = 0.D+00
        ad_1_pos    = 0.0D+00
        ad_1_neg    = 0.0D+00
        ad_1_im     = 0.D+00

        ad_2        = 0.D+00
        ad_2_pos    = 0.0D+00
        ad_2_neg    = 0.0D+00
        
        do im = 1, nsm


C#
C#   No normalization factor needed for the intermediate states
C#
C#              

C#      Find the detuning    DE(m) = (En - Ei)^2 - w^2
C#
C#     W         in a.u. (read by stark.inp) file
C#     En, Ei    in Ryd  (read by dmx2e-sLL+1.dat files)
C#


         de_in = ( enm(im) - eni(nnf)) / 2 





C#
C#     Let the off-pole part with the same photon   eph = w1
C#

         de_pos = de_in  + eph





C#
C#     calculate  off-resonant part first
C#

         ad_2  = ( dmx(im, nnf,1) * dmx(im, nnf, 2) ) /de_pos

c         ad_2  =  dmx(im,nnf,1)**2 /de_pos


C#     For numerical accuracy avoid subtraction between 
C#     small numbers.
C#    


         if(ad_2.gt.0D0) then
            
              ad_2_pos = ad_2_pos + ad_2

           else
            
              ad_2_neg = ad_2_neg + ad_2

           endif

C#
C#    Change the photon to the on-resonace part to match exactly to
C#    the pole
C#
C#           ewph = ...
C#


         de_neg = de_in  - ewph

           if(im.eq.ipole) then
              
              
              ad_1_im = - pi * dmx(im,nnf,1) * dmx(im,nnf,2)  
     1                  * rho_im**2

c              write(*,*) 'de_pos, de_neg', de_neg,de_pos,ewph
        else


           ad_1  = ( dmx(im, nnf,1) * dmx(im, nnf, 2) ) / de_neg

c           ad_1  = dmx(im,nnf,1)**2 / de_neg
           

C# 
C#       add negative and positive values separately
C#
           if(ad_1.gt.0D0) then

              ad_1_pos = ad_1_pos + ad_1

           else
            
              ad_1_neg = ad_1_neg + ad_1

           endif


        endif

c              write(*,*) 'de_neg = ', de_neg
      enddo


C#
C#     Renormalize result because DME's are discrete
C#

       ad_2 = rho_j * rho_j * ( ad_2_pos + ad_2_neg )

       ad_1 = rho_j * rho_j * ( ad_1_pos + ad_1_neg )

       ad_1_im = rho_j * rho_j * ad_1_im 


       write(*,*)  'M(Lm,-) = ', ad_2
       write(*,*)  'M(Lm,+) = ', ad_1

C#
C#
C#    stark = - S_(W) / (E^2/4)  == a_d(W) 
C#
C#
C#   Write to the standard oupout:    mode = 0  velocity
C#                                    mode = 1  length
C#

C#######################################################################

           if(mode.eq.1) then

             
              stark_re =   ad_1 + ad_2
           
              stark_im =  ad_1_im
           else
              
              stark_re =  (ad_1 + ad_2 )/ (ewph**2) 

              stark_im =  ad_1_im/ (ewph**2) 

           endif


           write(16,5)  ewph * 27.207696, stark_re, stark_im

 
           write(*,*)  'ad(',ewph * 27.207696,') = ' , stark_re


c     end do


 180  CONTINUE


      close(16) 
C#                    FORMAT STATEMENTS
C#####################################################################

 577  FORMAT(8A14)
    4 FORMAT(/2X,'ENERGY EIGENVALUES OF INITIAL & FINAL STATES')
    5 FORMAT(2X,1P8E15.7)
    6 FORMAT(/2X,'OS Values: Length - top & Velocity - bottom.   {LF is
     & included. For Absorp. LF=2*LOF+1. For Emiss.("-" sign) LF=2*LOI+1
     1 }' /2X,'ROW - FINAL STATES & COLUMN - INITIAL STATES')
    7 FORMAT(2X)
    8 FORMAT(/2X,'LHFF,LLF,LHFI,LLI,ID(I),ANG(I) --'/)
    9 FORMAT(2X,4I2,3X,4I2,1P4E12.4)
   11 FORMAT(/2X,'EXCITATION ENERGY IN RYD.')
 55   format(4i5)


      END
C#####################################################################
C#EOF
