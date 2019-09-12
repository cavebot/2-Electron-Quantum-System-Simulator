C#############################
      subroutine interv(xt,lxt,x,left,mflag)
      implicit real*8(a-h,o-z)
      dimension xt(lxt)
      data ilo /1/

      ihi = ilo + 1
      if(ihi .lt. lxt)                go to 20
      if(x .ge. xt(lxt))              go to 110
      if(lxt . le. 1)                 go to 90
      ilo = lxt - 1
      ihi = lxt

 20   if(x .ge. xt(ihi))              go to 40
      if(x .ge. xt(ilo))              go to 100

      istep = 1
 31   continue
         ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)              go to 35
         if (x .ge. xt(ilo))          go to 50
         istep = istep*2
      go to 31
 35   ilo = 1
      if(x .lt. xt(1))                go to 90
                                      go to 50

 40   istep = 1
 41   continue
         ilo = ihi
         ihi = ilo + istep
         if(ihi .ge. lxt)             go to 45
         if(x .lt. xt(ihi))           go to 50
         istep = istep*2
      go to 41

 45   if (x .ge. xt(lxt))             go to 110
      ihi = lxt
 50   continue
            middle = (ilo + ihi)/2
            if(middle .eq. ilo)       go to 100
            if(x .lt. xt(middle))     go to 53
            ilo = middle
          go to 50
 53       ihi = middle
      go to 50

 90   mflag = -1
      left = 1
      RETURN

 100  mflag = 0
      left = ilo
      RETURN

 110  mflag = 1
      left = lxt
      RETURN
      end
C#######################################################################
C#
C#   input:
C#
C#        t()    =  knot sequence from 'mkgrid' routine
C#        k      =  order of B-splines
C#        x      =  r (at chosen point for gaussian integration)
C#        left   =  index for the knot
C#        nderiv =  how many derivatives , e.g. two : 0th and first
C#
C#   output:
C#
C#        dbiatx() =  derivatives,e.g. two 0th=(:,1), 1st=(:,2)
C#   external function used = bsplvb()
C#
C#######################################################################

      subroutine bsplvd(t, k, x, left, dbiatx, nderiv)
      implicit real*8(a-h,o-z)
C#   changes from de Boor 
      parameter(KX=15)
      dimension a(KX,KX),dbiatx(KX,nderiv)
      dimension t(*)
C#
      mhigh = max0(min0(nderiv,k),1)
      kp1 = k + 1
      call bsplvb(t,kp1-mhigh,1,x,left,dbiatx)
      if(mhigh .eq. 1) go to 99

      ideriv = mhigh
      do 15 m = 2,mhigh
          jp1mid = 1
          do 11 j = ideriv,k
             dbiatx(j,ideriv) = dbiatx(jp1mid,1)
             jp1mid = jp1mid + 1
 11       continue
          ideriv = ideriv - 1
          call bsplvb(t,kp1-ideriv,2,x,left,dbiatx)
 15   continue

      jlow = 1
      do 20 i = 1,k
         do 19 j = jlow,k
            a(j,i) = 0.0D0
 19      continue
         jlow = i
         a(i,i) = 1.0D0
 20   continue

      do 40 m = 2,mhigh
         kp1mm = kp1 - m
         fkp1mm = DBLE(kp1mm)
         il = left
         i = k

         do 25 ldummy = 1,kp1mm
            factor = fkp1mm/(t(il+kp1mm) - t(il))
            do 24 j = 1,i
               a(i,j) = (a(i,j) -a(i-1,j))*factor
 24         continue
            il = il - 1
            i = i - 1
 25      continue

         do 36 i = 1,k
              sum = 0.0D0
              jlow = max0(i,m)
              do 35 j = jlow,k
                 sum = a(j,i)*dbiatx(j,m) + sum
 35           continue
              dbiatx(i,m) = sum
 36      continue
 40   continue

 99   RETURN
      end
********************************************************
      subroutine bsplvb( t, jhigh, index, x, left, biatx)
      implicit real*8(a-h,o-z)
      parameter(JMAX=20)
      dimension biatx(jhigh),t(*),deltal(jmax),deltar(jmax)
      data j/1/

      go to (10,20),index

 10   j = 1
      biatx(1) = 1.0D0
      if(j .ge. jhigh)                 go to 99

 20   continue
         jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.0D0
         do 26 i = 1,j
             term = biatx(i)/(deltar(i) + deltal(jp1-i))
             biatx(i) = saved + deltar(i)*term
             saved = deltal(jp1-i)*term
 26      continue
         biatx(jp1) = saved
         j = jp1
         if(j .lt. jhigh)
     *go to 20

 99   RETURN
      end

C------------------------------------------------------------------------
      function bvalue(t,bcoef,n,k,x,jderiv)
      implicit real*8(a-h,o-z)
      parameter(KMAX=15)
      dimension bcoef(n),t(*),aj(KMAX),dl(KMAX),dr(KMAX)

      bvalue = 0.0D0
      if(x.eq.t(n+k)) then
         bvalue=bcoef(n)
      end if 
      if(jderiv .ge. k)                go to 99
      call interv(t,n+k,x,i,mflag)
      if(mflag .ne. 0)                 go to 99

      km1 = k - 1
      if(km1 . gt. 0)                  go to 1
      bvalue = bcoef(i)
                                       go to 99

 1    jcmin = 1
      imk = i - k
      if(imk .ge. 0)                    go to 8
      jcmin = 1 - imk

      do 5 j = 1,i
         dl(j) = x - t(i+1-j)
 5    continue
      do 6 j = i,km1
         aj(k-j)  = 0.0D0
         dl(j) = dl(i)
 6    continue
                                      go to 10

 8    do 9 j = 1,km1
         dl(j) = x - t(i+1-j)
 9    continue

 10   jcmax = k
      nmi = n - i
      if(nmi .ge. 0)                  go to 18
      jcmax = k + nmi
      do 15 j = 1,jcmax
          dr(j) = t(i+j) - x
 15   continue
      do 16 j = jcmax,km1
          aj(j+1) = 0.0D0
          dr(j) = dr(jcmax)
 16   continue
                                     go to 20

 18   do 19 j = 1,km1
        dr(j) = t(i+j) - x
 19   continue

 20   do 21 jc = jcmin,jcmax
          aj(jc) = bcoef(imk + jc)
 21   continue

      if(jderiv .eq. 0)             go to 30
      do 23 j = 1,jderiv
         kmj = k - j
         fkmj = DBLE(kmj)
         ilo = kmj
         do 22 jj = 1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
            ilo = ilo - 1
 22      continue
 23   continue

 30   if(jderiv .eq. km1)           go to 39
      do 33 j = jderiv+1,km1
         kmj = k - j
         ilo = kmj
         do 32 jj = 1,kmj
           aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
           ilo = ilo - 1
 32      continue
 33   continue
 39   bvalue = aj(1)
 
 99   RETURN
      end
c------------------------------------------------------------------------

      subroutine gauss(k,x,w)
************************************************
*
*   x(i) = gaussian coordinates for [0,1]
*   w(i) = gaussian coordinates for [0,1]
*   1<= k <= 15   for k point case
*
*************************************************
      implicit double precision(a-h,o-z)
      dimension x(k),w(k)

      if(k.lt.1.or.k.gt.15) go to 901
      goto(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),k
 10   x(1) = 0.5D00
      w(1) = 1.0D0
      go to 900
 20   x(1) = .211324865405187d0
      x(2) = .788675134594813d0
      w(1) = 0.5D00 
      w(2) = 0.5D00
      go to 900
 30   x(1) = .112701665379258d0 
      x(2) = 0.5D00               
      x(3) = .887298334620742d0
      w(1) = .277777777777778d0 
      w(2) = .444444444444444d0
      w(3) = .277777777777778d0
      go to 900
 40   x(1) = .0694318442029737d0
      x(2) = .330009478207572d0
      x(3) = .669990521792428d0
      x(4) = .930568155797026d0
      w(1) = .173927422568727d0
      w(2) = .326072577431273d0
      w(3) = .326072577431273d0
      w(4) = .173927422568727d0
      go to 900
 50   x(1) = .046910077030668d0
      x(2) = .230765344947158d0
      x(3) = 0.5D00               
      x(4) = .769234655052842d0
      x(5) = .953089922969332d0
      w(1) = .118463442528095d0
      w(2) = .239314335249683d0
      w(3) = .284444444444444d0
      w(4) = .239314335249683d0
      w(5) = .118463442528095d0
      go to 900
 60   x(1) = .033765242898424d0
      x(2) = .169395306766868d0
      x(3) = .380690406958402d0
      x(4) = .619309593041598d0
      x(5) = .830604693233132d0
      x(6) = .966234757101576d0
      w(1) = .0856622461895852d0
      w(2) = .180380786524069d0
      w(3) = .233956967286346d0
      w(4) = .233956967286346d0
      w(5) = .180380786524069d0
      w(6) = .0856622461895852d0
      go to 900
 70   x(1) = .0254460438286207d0
      x(2) = .129234407200303d0
      x(3) = .297077424311301d0
      x(4) = 0.5D00
      x(5) = .702922575688699d0
      x(6) = .870765592799697d0
      x(7) = .974553956171379d0
      w(1) = .0647424830844348d0
      w(2) = .139852695744638d0
      w(3) = .19091502525256d0
      w(4) = .208979591836735d0
      w(5) = .19091502525256d0
      w(6) = .139852695744638d0
      w(7) = .0647424830844348d0
      go to 900
 80   x(1) = .0198550717512319d0
      x(2) = .101666761293187d0
      x(3) = .237233795041835d0
      x(4) = .408282678752175d0
      x(5) = .591717321247825d0
      x(6) = .762766204958164d0
      x(7) = .898333238706813d0
      x(8) = .980144928248768d0
      w(1) = .0506142681451881d0
      w(2) = .111190517226687d0
      w(3) = .156853322938944d0
      w(4) = .181341891689181d0
      w(5) = .181341891689181d0
      w(6) = .156853322938944d0
      w(7) = .111190517226687d0
      w(8) = .0506142681451881d0
      go to 900
 90   x(1) = .015919880246187d0
      x(2) = .0819844463366821d0
      x(3) = .193314283649705d0
      x(4) = .337873288298095d0
      x(5) = 0.5D00                
      x(6) = .662126711701904d0
      x(7) = .806685716350295d0
      x(8) = .918015553663318d0
      x(9) = .984080119753813d0
      w(1) = .0406371941807872d0
      w(2) = .0903240803474287d0
      w(3) = .130305348201468d0
      w(4) = .156173538520001d0
      w(5) = .16511967750063d0
      w(6) = .156173538520001d0
      w(7) = .130305348201468d0
      w(8) = .0903240803474287d0
      w(9) = .0406371941807872d0
      go to 900
 100  x(1) = .0130467357414141d0
      x(2) = .0674683166555077d0
      x(3) = .160295215850488d0
      x(4) = .283302302935376d0
      x(5) = .425562830509184d0
      x(6) = .574437169490816d0
      x(7) = .716697697064624d0
      x(8) = .839704784149512d0
      x(9) = .932531683344492d0
      x(10)= .986953264258586d0
      w(1) = .0333356721543441d0
      w(2) = .0747256745752903d0
      w(3) = .109543181257991d0
      w(4) = .134633359654998d0
      w(5) = .147762112357376d0
      w(6) = .147762112357376d0
      w(7) = .134633359654998d0
      w(8) = .109543181257991d0
      w(9) = .0747256745752903d0
      w(10)= .0333356721543441d0
      go to 900
 110  x(1) = .0108856709269715d0
      x(2) = .0564687001159523d0
      x(3) = .134923997212975d0
      x(4) = .240451935396594d0
      x(5) = .365228422023827d0
      x(6) = 0.5D00                
      x(7) = .634771577976172d0
      x(8) = .759548064603406d0
      x(9) = .865076002787025d0
      x(10)= .943531299884048d0
      x(11)= .989114329073028d0
      w(1) = .0278342835580868d0
      w(2) = .0627901847324523d0
      w(3) = .0931451054638672d0
      w(4) = .116596882295995d0
      w(5) = .131402272255123d0
      w(6) = .13646254338895d0
      w(7) = .131402272255123d0
      w(8) = .116596882295995d0
      w(9) = .0931451054638672d0
      w(10)= .0627901847324523d0
      w(11)= .0278342835580868d0
      go to 900
 120  x(1) = .00921968287664038d0
      x(2) = .0479413718147626d0
      x(3) = .115048662902848d0
      x(4) = .206341022856691d0
      x(5) = .31608425050091d0
      x(6) = .437383295744266d0
      x(7) = .562616704255734d0
      x(8) = .68391574949909d0
      x(9) = .793658977143309d0
      x(10)= .884951337097152d0
      x(11)= .952058628185237d0
      x(12)= .99078031712336d0
      w(1) = .0235876681932559d0
      w(2) = .0534696629976592d0
      w(3) = .0800391642716731d0
      w(4) = .101583713361533d0
      w(5) = .116746268269177d0
      w(6) = .124573522906701d0
      w(7) = .124573522906701d0
      w(8) = .116746268269177d0
      w(9) = .101583713361533d0
      w(10)= .0800391642716731d0
      w(11)= .0534696629976592d0
      w(12)= .0235876681932559d0
      go to 900
 130  x(1) = .00790847264070593d0
      x(2) = .041200800388511d0
      x(3) = .099210954633345d0
      x(4) = .17882533027983d0
      x(5) = .275753624481777d0
      x(6) = .384770842022433d0
      x(7) = 0.5D00                 
      x(8) = .615229157977567d0
      x(9) = .724246375518223d0
      x(10)= .82117466972017d0
      x(11)= .900789045366655d0
      x(12)= .958799199611489d0
      x(13)= .992091527359294d0
      w(1) = .0202420023826579d0
      w(2) = .0460607499188642d0
      w(3) = .0694367551098937d0
      w(4) = .0890729903809729d0
      w(5) = .103908023768444d0
      w(6) = .113141590131449d0
      w(7) = .116275776615437d0
      w(8) = .113141590131449d0
      w(9) = .103908023768444d0
      w(10)= .0890729903809729d0
      w(11)= .0694367551098937d0
      w(12)= .0460607499188642d0
      w(13)= .0202420023826579d0
      go to 900
 140  x(1) = .00685809565159383d0
      x(2) = .0357825581682132d0
      x(3) = .0863993424651175d0
      x(4) = .156353547594157d0
      x(5) = .242375681820923d0
      x(6) = .340443815536055d0
      x(7) = .445972525646328d0
      x(8) = .554027474353672d0
      x(9) = .659556184463945d0
      x(10)= .757624318179077d0
      x(11)= .843646452405843d0
      x(12)= .913600657534882d0
      x(13)= .964217441831787d0
      x(14)= .993141904348406d0
      w(1) = .0175597301658759d0
      w(2) = .0400790435798801d0
      w(3) = .0607592853439516d0
      w(4) = .0786015835790968d0
      w(5) = .092769198738969d0
      w(6) = .102599231860648d0
      w(7) = .107631926731579d0
      w(8) = .107631926731579d0
      w(9) = .102599231860648d0
      w(10)= .092769198738969d0
      w(11)= .0786015835790968d0
      w(12)= .0607592853439516d0
      w(13)= .0400790435798801d0
      w(14)= .0175597301658759d0
      go to 900
 150  x(1) = .00600374098975728d0
      x(2) = .031363303799647d0
      x(3) = .0758967082947864d0
      x(4) = .137791134319915d0
      x(5) = .214513913695731d0
      x(6) = .302924326461218d0
      x(7) = .399402953001283d0
      x(8) = 0.5D00                 
      x(9) = .600597046998717d0
      x(10)= .697075673538782d0
      x(11)= .785486086304269d0
      x(12)= .862208865680085d0
      x(13)= .924103291705214d0
      x(14)= .968636696200353d0
      x(15)= .993996259010243d0
      w(1) = .0153766209980586d0
      w(2) = .0351830237440541d0
      w(3) = .053579610233586d0
      w(4) = .0697853389630772d0
      w(5) = .083134602908497d0
      w(6) = .0930805000077812d0
      w(7) = .0992157426635559d0
      w(8) = .101289120962781d0
      w(9) = .0992157426635559d0
      w(10)= .0930805000077812d0
      w(11)= .083134602908497d0
      w(12)= .0697853389630772d0
      w(13)= .053579610233586d0
      w(14)= .0351830237440541d0
      w(15)= .0153766209980586d0

 900  return

 901  stop 'error in gauss'
      end
