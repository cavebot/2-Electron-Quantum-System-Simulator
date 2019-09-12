 
      PROGRAM tdse_pes_bs_fxd_2e
        !
        USE PRECISION
        USE units
        USE parameter_tdse_fxd, ONLY: atomic_name, lmax, nmax, en_ion_1, en_ground 
        USE parameter_tdse_fxd, ONLY: bin_type, en_cut, gauge, hohg, irestart, tol
        USE parameter_tdse_fxd, ONLY: input_tdse_fxd
!        USE PARAMETER
        USE io
        !
        IMPLICIT NONE
        !
        REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: en_nl, pop
        REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: yr, yi
        REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: den, pes
        INTEGER  , ALLOCATABLE, DIMENSION(:)   :: ns
        !
        INTEGER,   ALLOCATABLE, DIMENSION(:)   :: ndi_l                ! pop with  E(ndi_l,l) = E++
        REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pop_si_l             ! pop with  E(ndi_l,l) < E++ 
        REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pop_di_l             ! pop with  E(ndi_l,l) > E++   
        REAL(dpk)                              :: pop_si               ! pop with  E < E++ 
        REAL(dpk)                              :: pop_di               ! pop with  E > E++
        !
        INTEGER                                :: l, l_pes, nl
        INTEGER                                :: ne, nn, ne_tot, ne_tot_l
        INTEGER                                :: nnout, npes, npop, npop_l_t, npop_sdi
        REAL(dpk)                              :: frm
        CHARACTER(len=1)                       :: ctmp
        INTEGER                                :: ntmp
        !
        !pop_l_t
        REAL(dpk)                              :: time
        REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: p_l_t
        REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: p_el_t
        INTEGER                                :: nof_t_steps
        INTEGER                                :: l_max_state
        INTEGER                                :: nof_e_states
        !
        INTEGER                                :: ie, il, it, it_index
        INTEGER                                :: nt, time_steps
        !
        CHARACTER(len=15)                      :: action         
        INTERFACE
!           SUBROUTINE get_command_line_arg(action,na,la,nb,lb,gauge,interaction)
           SUBROUTINE get_command_line_arg(action)
             IMPLICIT NONE
             CHARACTER(len=15), INTENT(inout) :: action 
           END SUBROUTINE get_command_line_arg
        END INTERFACE
        
        !


        !EX!

        ! first read tdse input file (tinp/tdse_bs.inp)
        
        CALL input_tdse_fxd          ! tinp/tdse_bs_fxd_2e.inp
        
        ! default parameters (if not command line present)


        action  ='pes'
        CALL get_command_line_arg(action)
        



        
        calculate_pes:IF(action == 'pes') THEN

           WRITE(*,'(a1,a60)')'&',' pes::  Calculating photoelectron energy spectrum.'

           npes   = 10
           
        ! read input partial pes file
           
        OPEN(npes,FILE='tdat/pesl.dat')              !open spectra file
           
        READ(npes,'(a1,1x,2i6)') ctmp, nl, ne_tot    !'&', nof partial waves, nof states

        IF(lmax.NE.(nl-1)) THEN

           WRITE(*,*) "# tdse_pes:: variable (lmax) from (tinp/tdse_bs.inp) differ than variable (nl-1) in (tdat/pesl.dat) file."
           WRITE(*,*) "# tdse_pes:: exit"
           !           lmax = nl - 1

        END IF
           

        ne_tot_l = INT(ne_tot/nl) ! fails if nl = 1,3,5,...

        WRITE(*,*) '& pes           nl = ', nl
        WRITE(*,*) '& pes         lmax = ', lmax
        WRITE(*,*) '& pes       ne_tot = ', ne_tot
        WRITE(*,*) '& pes     ne_tot_l = ', ne_tot_l


        ALLOCATE(               ns(0:lmax) )
        ALLOCATE(  en_nl(ne_tot_l, 0:lmax) )
        ALLOCATE(    pop(ne_tot_l, 0:lmax) )
        ALLOCATE(     yr(ne_tot_l, 0:lmax) )
        ALLOCATE(     yi(ne_tot_l, 0:lmax) )

        !
        ! pop(ne,l) =  y( ne + ntot_l )**2 + y(ne + ntot_l + ntot)**2 
        !

        
        partial_wave_l:DO  l = 0, lmax                  ! read for each partial wave

           READ(npes,'(a1,1x,2i6)') ctmp, l_pes, ntmp  ! partial wave, nof states for 'l' partial wave
           
           WRITE(*,*) '& pes:     l = ', l_pes
           eigenstate_of_partial_wave_l: DO  ne = 1, ne_tot_l

              READ(npes,'(4e14.6)') en_nl(ne,l), pop(ne,l), yr(ne,l), yi(ne,l)
!              WRITE(*,*) ne,l,ne_tot_l, en_nl(ne,l)

           ENDDO eigenstate_of_partial_wave_l

           ns(l) = ne_tot_l
        ENDDO partial_wave_l
        
        CLOSE(npes)
        
        WRITE(*,*) "& pes: partial spectra file read."

       
        !# divide energy axis in DE pieces (according 'lmax' discretization)

        ALLOCATE(  den(ne_tot_l) )

        DO ne = 1, ns(lmax) - 1  
           den(ne) = ABS( en_nl(ne+1,lmax) - en_nl(ne,lmax))
        ENDDO

        WRITE(*,*) "& pes: discretization done."
        

!      write(*,*)  ' energies read :', nmax - 1,lmax
!      write(*,*)  en_nl(1,lmax), en_nl(ns(lmax)-1,lmax)

        ALLOCATE(  pes(ne_tot_l) )


        
        pes = 0.0_dpk
         sum_partial_wave: DO  l = 0, lmax

            nn = 1
            DO  ne = 2, ns(l) - 1
             
               nn = ne
 
               ! determine the states (n,l) with
               !   en_nl(n,l) 
               !        in 
               !  [en_nl(ne,lmax)-0.5*den(ne), en_nl(ne,lmax)+0.5*den(ne)]
               ! ....
!          
230            nn = nn - 1
               IF( en_nl(nn,l).GT.(en_nl(ne,lmax) + 0.5D+00 * den(ne)) )  GOTO 230 
!                  
240            nn = nn + 1
               IF( en_nl(nn,l).LT.(en_nl(ne,lmax) - 0.5D+00 * den(ne)) )  GOTO 240 

               !.. and add their contribution

             pes(ne) = pes(ne) + pop(nn,l) 
             !             WRITE(*,'(3i3,1x,2E15.6)') ne, l, nn, pes(ne), pop(nn,l)

          ENDDO
       ENDDO sum_partial_wave
            

       npop = 10
       OPEN(npop,file='tdat/pes.dat')

!            write(*,*) '$DATA=CURVE3D' 
!            write(*,*) '$DATA=VECTOR' 
!            write(*,*) 'linecolor=0' 
            
       frm = 1.0D+00
       store_pes:DO ne = 1, ns(lmax)-1               
          continuum_spectrum:IF(en_nl(ne-1,lmax).GT.0.0D+00) THEN
             
             frm = 2.0D+00/(den(ne)+den(ne-1)) 

             WRITE(npes,'(2e14.6)') en_nl(ne,lmax) * enau, pes(ne)*frm/enau   !  S = S(E)
             
          ENDIF continuum_spectrum
       ENDDO store_pes

       close(npes)

       WRITE(*,'(a1,a60)')'&',' pes calculated.'


       



    !
    !
    !  calculate continuum population below (pop_si) and above (pop_di) DI threshold (E++)
    !
    !
       !

    PRINT*, en_ion_1

    ALLOCATE(  pop_di_l(0:lmax) )
    ALLOCATE(  pop_si_l(0:lmax) )
    ALLOCATE(     ndi_l(0:lmax) ) 

    pop_di_l = 0.0_dpk
    pop_si_l = 0.0_dpk
    pop_si = 0.0_dpk
    pop_di = 0.0_dpk
    partial_waves: DO  l = 0, lmax

       
       find_ndi_l_threshold_for_di:DO  ne = 1, ns(l)
          
          ndi_l(l) = ne
          IF( (en_nl(ne,l)+en_ion_1).GE.0.0_dpk) EXIT
          
       ENDDO find_ndi_l_threshold_for_di
       
       WRITE(*,*) "tdse_pes::   ndi_l(",l,") = ", ndi_l(l),  " ----  E = ", en_nl(ndi_l(l),l) 

       
       pop_si_l(l) = SUM( pop(         1: ndi_l(l), l) )
       pop_di_l(l) = SUM( pop(  ndi_l(l): ns(l)-1,  l) )
       
    ENDDO partial_waves


    pop_si = SUM(pop_si_l)   ! all population  with (E+  < E < E++)
    pop_di = SUM(pop_di_l)   ! all population  with (E++ < E ) ! mixed with
    


    npop_sdi = 11
    OPEN(npop_sdi,file='tdat/pop_sdi.dat')
       
    WRITE(npop_sdi, '(20e14.6)') pop_si, pop_si_l(0:lmax)
    WRITE(npop_sdi, '(20e14.6)') pop_di, pop_di_l(0:lmax)

    CLOSE(npop_sdi)

    WRITE(*,'(a1,a60)')'&',' populations for E<E++ and E>E++ are calculated.'


    

 END IF calculate_pes


 !
 !
 !
 ! POPULATION(S)
 !
 !
 !

    
    calculate_pop_l_t:IF(action == 'pop') THEN


       
       WRITE(*,'(a1,a60)')'&',' pop:: calculate population in time'
       npop   = 11
       npop_l_t  = 12

       OPEN(npop, file='tdat/pop.dat')
       OPEN(npop_l_t,file='tdat/pop_l_t.dat')
       !       OPEN(npopt,file='tdat/pop_e_t.dat')
       
       READ(npop,*) nof_t_steps, l_max_state, nof_e_states

       ALLOCATE(     p_l_t( 1:nof_t_steps,  0:l_max_state) )
       ALLOCATE(    p_el_t( 1:nof_e_states, 0:l_max_state) )



       time_loop:DO it = 1, nof_t_steps

          READ(npop,'(i8, 1X,e15.6)') it_index, time

       read_pop_partial_waves: DO il = 0, l_max_state
          READ(npop,'(I6,1X,5e25.14)') l, p_l_t(it, il), ( p_el_t(ie,il), ie=1,nof_e_states)
          
       ENDDO read_pop_partial_waves
       READ(npop,'(/)')

!       WRITE(npop_l_t,'(1e15.6,2X,<l_max_state+1>e25.14)')                    & 
       WRITE(npop_l_t,'(1e15.6,2X,10e25.14)')                    & 
            &    time                                                         & 
            &,   (p_l_t(it, il), il = 0, l_max_state) 
       !    &,     (p_el_t, ie = 1, nof_e_states)

       
    ENDDO time_loop
    
    
    close(npop)
    CLOSE(npop_l_t)
    

    time_steps = 10
    nt = int(nof_t_steps/time_steps)
    it = 1
    OPEN(npop_l_t,file='tdat/pop_l.dat')
    partial_waves_1:DO il = 0, l_max_state
!       WRITE(npop_l_t,'(I6,1X,<nt>e25.14)') il, (p_l_t(it, il), it = 1, nof_t_steps, time_steps)
       WRITE(npop_l_t,'(I6,1X, 10e25.14)') il, (p_l_t(it, il), it = 1, nof_t_steps, time_steps)
    ENDDO partial_waves_1
    CLOSE(npop_l_t)

       WRITE(*,'(a1,a60)')'&',' pop:: population in time calculated'

    ENDIF calculate_pop_l_t
    


            
!            ifail = 0            
!... find the first 4 peaks in the pes!!
!            call m01caf(pes,1, ns(lmax),'D',ifail)
!            write(*,'(5e14.6)') i0, (pes(n), n = 1, 3)


!	   call qsort(pes, ns(lmax), 8, 1)
!           ntdout = 88
!           open(ntdout,file='pes/pes_i.dat', position='append')
!           read(ntdout,'(5e14.6)', end = 99) ps
! 99        write(ntdout,'(5e14.6)') i0, (pes(n), n = 1, 4)
!           write(ntdout,'(5e14.6)') i0, (pes(n), n = 1, 4)
!           close(ntdout)
!            write(*,*)'$END' 
!            close(npes)
!         enddo




          END PROGRAM tdse_pes_bs_fxd_2e
!comman_line_sub

SUBROUTINE get_command_line_arg(action)
  !
  USE PRECISION, ONLY: DPK
  !
  IMPLICIT NONE
  !ARG!
  CHARACTER(len=15), INTENT(inout) :: action 
  !LOC!
  CHARACTER*100 ARGV
  CHARACTER*180 LINE, EXE, WHAT
  CHARACTER*40  CMD 
  INTEGER       NARG,IARG, NXT_ARG
  INTEGER       LENGTH, ISTATUS
  INTEGER       COMMAND_ARGUMENT_COUNT
  CHARACTER(len=10)       :: date,time,zone
  INTEGER, DIMENSION(dpk) :: values 
  INTEGER nascii
  !EXE!
 
  nascii = 1
  
  CALL DATE_AND_TIME(date, time, zone, values)
  
  NARG = COMMAND_ARGUMENT_COUNT()    
  CALL GET_COMMAND( LINE, LENGTH, ISTATUS)     
  CALL GET_COMMAND_ARGUMENT(0,EXE,LENGTH,ISTATUS) 


  WRITE(*,*)'#'
  WRITE(*,'(a45,1X,a40)')'# w1e::            executable command exe = ', line
  WRITE(*,'(a45,1X,i2)') '# wf1e::        nof arguments        narg = ', narg  

  IF(narg.LE.1) THEN
!!%     WRITE(*,'(a45,1X,i3 )') '# wf1e::                                na = ', na
!!%     WRITE(*,'(a45,1X,i3)')  '# wf1e::                                la = ', la
!!%     WRITE(*,'(a45,1X,i3)')  '# wf1e::                                nb = ', nb
!!%     WRITE(*,'(a45,1X,i3 )') '# wf1e::                                lb = ', lb
     WRITE(*,'(a45,1X,a40)') '# wf1e::                            action = ', action
!!%     WRITE(*,'(a45,1X,a40)') '# wf1e::                             gauge = ', gauge
!!%     WRITE(*,'(a45,1X,a40)') '# wf1e::                       interaction = ', interaction
  ENDIF

  get_arguments: DO iarg = 1, narg

     CALL GET_COMMAND_ARGUMENT(IARG,CMD,LENGTH,ISTATUS)

     IF(MOD(iarg,2)==1) THEN    !iarg = 1, 3, 5, ...

        READ(cmd,*) what

        IF(what=='-o') THEN 
           WRITE(*,*) '# wf1e::'
           WRITE(*,*) '# wf1e:: Available options:'
!!%           WRITE(*,*) '# wf1e:: 0. -o information '
!!%           WRITE(*,*) '# wf1e:: 1. na= eigenstate index (integer,   > 0),        [1]'
!!%           WRITE(*,*) '# wf1e:: 2. la= angular symmetry (integer,  >= 0),        [0]'
!!%           WRITE(*,*) '# wf1e:: 3. nb= eigenstate index (integer,   > 0),        [1]'
!!%           WRITE(*,*) '# wf1e:: 4. lb= angular symmetry (integer,  >= 0),     [la+1]'
           WRITE(*,*) '# wf1e:: 1. a= action           (character,l=40),  [save]'
           WRITE(*,*) '# wf1e::            (pes):  photoelectron energy spectrum'
           WRITE(*,*) '# wf1e::            (pop):  population as a function of time'  
           WRITE(*,*) '# wf1e::            (hhg):  high-harmonic generation signal'  
           WRITE(*,*) '# wf1e::            (pad):  photoelectron angular distribution'  
           WRITE(*,*) '# wf1e::            (rho): calculates all <P_nl|B_j> integrals'  
!!%           WRITE(*,*) '# wf1e:: 6. g= gauge          (character,l=6),  [l],v,a,lva'
!!%           WRITE(*,*) '# wf1e:: 7. i= interaction    (character,l=6), save data for [tdse],lopt
           WRITE(*,*) '# wf1e:: happy end.'  
           STOP
        ENDIF

        nxt_arg = iarg + 1

        WRITE(*,'(a45,1X,a40)') '# wf1e::                             what = ', what

        IF(what.NE. 'a='        ) THEN     
           WRITE(*,'(a60)') 'available options for :'
           WRITE(*,'(a60)') ' (a= (pes), pop, hhg, pad, rho)'
           STOP
        ENDIF
         
        CYCLE
     
     ENDIF

!...................     

     IF(iarg == nxt_arg ) THEN

!!%        IF(what=='na=') THEN 
!!%           READ(cmd,'(I3)')  na    
!!%           WRITE(*,'(a45,1X,i3)') '# wf1e::                                na = ', na
!!%        ELSE IF(what=='la=') THEN
!!%           READ(cmd,*)  la 
!!%           WRITE(*,'(a45,1X,i3 )') '# wf1e::                               la = ', la
!!%        ELSE IF(what=='nb=') THEN
!!%           READ(cmd,*)  nb 
!!%           WRITE(*,'(a45,1X,i3 )') '# wf1e::                               nb = ', nb
!!%        ELSE IF(what=='lb=') THEN
!!%           READ(cmd,*)  lb 
!!%           WRITE(*,'(a45,1X,i3 )') '# wf1e::                               lb = ', lb

        IF(what=='a=') THEN

           action = cmd

           IF(   (action.NE.'pes' ).AND.&
      &          (action.NE.'pop' ).AND.&
      &          (action.NE.'hhg' ).AND.&
      &          (action.NE.'pad' ).AND.&
      &          (action.NE.'rho' ) ) then 
              WRITE(*,'(a60)') ' available options for a='
              WRITE(*,'(a60)') ' (pes)'
              WRITE(*,'(a60)') ' (pop)'
              WRITE(*,'(a60)') ' (hhg)'
              WRITE(*,'(a60)') ' (pad)'
              WRITE(*,'(a60)') ' (rho)'
           ENDIF

           WRITE(*,'(a45,1X,a40)') '# wf1e::                           action = ', action
!!%        ELSE IF(what=='g=') THEN
!!%          gauge = cmd
!!%          WRITE(*,'(a45,1X,a40)') '# wf1e::                           gauge = ', gauge
!!%        ELSE IF(what=='i=') THEN
!!%          interaction = cmd
!!%          WRITE(*,'(a45,1X,a40)') '# wf1e::                      interaction = ', interaction
        ELSE
           CYCLE
        ENDIF
        
     ENDIF
  ENDDO get_arguments
  
  !............... save history
  OPEN(nascii, file='log/wf1e_history.log', position='append')
  WRITE(nascii,'(a10,1x,a10,2x,a60)') date,time, line    
  CLOSE(nascii)


  RETURN
END SUBROUTINE get_command_line_arg
!eof



!#EOF
