!
!
! ``Bargaining over Taxes and Entitlements in the Era of Unequal Growth''
!                 Azzimonti, Karpuska, Mihalache (IER, 2022)
!
!
PROGRAM unequal
  USE, INTRINSIC :: iso_Fortran_env, ONLY: wp => real64, output_unit
  IMPLICIT NONE
  
  CHARACTER (LEN=*), PARAMETER :: outDir = "./results/"
  CHARACTER (LEN=1), PARAMETER :: tabChar = CHAR(9)
  REAL(wp), PARAMETER :: eulerMasc = 0.5772156649015328606065120_wp
  
  !
  ! Parameters
  !
  REAL(wp), PARAMETER :: q = 0.8_wp
  REAL(wp), PARAMETER :: beta = 0.97_wp
  REAL(wp), PARAMETER :: theta = 0.5_wp  ! weight on u(g) in objective
  
  REAL(wp), PARAMETER :: yP = 0.1_wp
  REAL(wp), PARAMETER :: yR = 1.1_wp
  REAL(wp), PARAMETER :: Y = yR + yP
  REAL(wp), PARAMETER :: DeltaY = yR - yP
  
  REAL(wp), PARAMETER :: xBar = 0.1_wp

  ! Numerical
  REAL(wp), PARAMETER :: errTol = 1.0D-12
  REAL(wp), PARAMETER :: dchoice = 5.0D-4
  REAL(wp), PARAMETER :: veryNegative = -1.0D+9
  INTEGER, PARAMETER :: maxIter = 2000
  INTEGER, PARAMETER :: n = 141, nn = n*((n-1)/2+1) ! Grid size

  REAL(wp), DIMENSION(n) :: tauGrid, eGrid
  
  !
  ! Rules
  !
  INTEGER, PARAMETER :: RULE_TAU = 1       ! tau only
  INTEGER, PARAMETER :: RULE_E = 2         ! e only
  INTEGER, PARAMETER :: RULE_BOTH = 3      ! tau and e
  INTEGER, PARAMETER :: RULE_G = 4         ! g
  INTEGER, PARAMETER :: RULE_ALL = 5       ! tau, e, and g
  INTEGER, PARAMETER :: RULE_NONE = 6      ! no status quo
  INTEGER, PARAMETER :: RULE = RULE_BOTH
 
  !
  ! Discretion
  !
  REAL(wp), PARAMETER :: cDiscretion = (Y - xBar) / (1 + theta)
  REAL(wp), PARAMETER :: gDiscretion = (Y - xBar) * theta / (1 + theta)
  REAL(wp), PARAMETER :: tauDiscretionR = yR - cDiscretion
  REAL(wp), PARAMETER :: eDiscretionR = xBar - yP
  REAL(wp), PARAMETER :: tauDiscretionP = yR - xBar
  REAL(wp), PARAMETER :: eDiscretionP = cDiscretion - yP
  REAL(wp) :: Wdisc, Vdisc, discTmp0(2), discTmp1(2), discU(2), discTran(2, 2)

  !
  ! Values and policies
  !
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: Vr0, Vp0, Wr0, Wp0
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: Vr1, Vp1, Wr1, Wp1
  REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: polR, polP
  
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: tau, e, g, cR, cP
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: uValR, uValP
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: uValOutR, uValOutP
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: outsideP, outsideR
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: gotoVal, gotoOther, probs
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: vvalid

  REAL(wp), ALLOCATABLE, DIMENSION(:, :) :: markovHuge, markovTiny
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: dist0, dist1, valTmp

  REAL(wp), ALLOCATABLE, DIMENSION(:) :: polRcR, polRcP, polRtau, polRe, polRg
  REAL(wp), ALLOCATABLE, DIMENSION(:) :: polPcR, polPcP, polPtau, polPe, polPg

  INTEGER :: ix, tauIx, eIx, iter, iunit
  REAL(wp) :: epsVal, theMax, theSum, theStat
  REAL(wp) :: errwp, errvp, errwr, errvr, evalR, evalP
  
  CALL tstamp()

  SELECT CASE (RULE)
    CASE (RULE_TAU)
      WRITE (*, *) "Rule: tau"
    CASE (RULE_E)
      WRITE (*, *) "Rule: e"
    CASE (RULE_BOTH)
      WRITE (*, *) "Rule: tau + e"
    CASE (RULE_G)
      WRITE (*, *) "Rule: g"
    CASE (RULE_ALL)
      WRITE (*, *) "Rule: tau + e + g"
    CASE DEFAULT ! No status quo
      WRITE (*, *) "Rule: no status quo"
  END SELECT

  !
  ! Allocate
  !
  WRITE (*, *) "Allocate..."
  ALLOCATE( Vr0(nn), Vp0(nn), Wr0(nn), Wp0(nn) )
  ALLOCATE( Vr1(nn), Vp1(nn), Wr1(nn), Wp1(nn) )
  ALLOCATE( polR(nn, nn), polP(nn, nn) )
  
  ALLOCATE( tau(nn), e(nn), g(nn), cR(nn), cP(nn) )
  ALLOCATE( uValR(nn), uValP(nn) )
  ALLOCATE( uValOutR(nn), uValOutP(nn) )
  ALLOCATE( outsideP(nn), outsideR(nn) )
  ALLOCATE( gotoVal(nn), gotoOther(nn), probs(nn), vvalid(nn) )
 
  !
  ! WrDisc = (1-beta) uOut + beta ( (1-q) VrDisc +     q WrDisc )
  ! VrDisc = (1-beta) uRin + beta ( q     VrDisc + (1-q) WrDisc )
  !
  discU = (1 - beta) * [ u(xBar, gDiscretion), u(cDiscretion, gDiscretion) ]
  discTran(1, :) = beta * [ 1-q, q ]
  discTran(2, :) = beta * [ q, 1-q ]

  discTmp0 = 0_wp
  discTmp1 = 0_wp
  iter = 1
  epsVal = 1.0_wp
  DO WHILE ( iter <= maxIter .AND. epsVal > errTol ) 
    discTmp1 = discU + MATMUL(discTran, discTmp0)
    epsVal = MAXVAL(ABS(discTmp1 - discTmp0))
    discTmp0 = discTmp1
    iter = iter + 1
  END DO

  Wdisc = discTmp1(1)
  Vdisc = discTmp1(2)
  WRITE (*, "(A,2F15.5)") "Discretion values: ", Wdisc, Vdisc

  Vp0 = u(xBar, xBar)
  Vp1 = Vp0
  Wp0 = Vp0
  Wp1 = Vp0
  Vr0 = Vp0
  Vr1 = Vp0
  Wr0 = Vp0
  Wr1 = Vp0
  
  !
  ! Setup grid
  !
  WRITE (*, *) "Setup..."
  WRITE (*, "(A,A20,I15)") tabChar, "nn: ", nn
  WRITE (*, "(A,A20,ES15.4)") tabChar, "taste shock: ", dchoice

  tauGrid = linspace(xBar + xBar - yP, yR - xBar, n)
  eGrid = linspace(xBar - yP, yR - xBar - xBar, n)

  WRITE (*, "(A,2F10.4)") "tau range: ", tauGrid(1), tauGrid(n)
  WRITE (*, "(A,2F10.4)") "  e range: ", eGrid(1), eGrid(n)
  
  ix = 0
  DO tauIx = 1,n
    DO eIx = 1,n
      IF (tauGrid(tauIx) - eGrid(eIx) >= xBar - 1.0D-8) THEN
        ! Setup triangular grid
        ix = ix + 1
        tau(ix) = tauGrid(tauIx)
        e(ix) = eGrid(eIx)
        g(ix) = tau(ix) - e(ix)
        cR(ix) = yR - tau(ix)
        cP(ix) = yP + e(ix)

        ! Flow values from consuming at location `ix'
        uValR(ix) = u(cR(ix), g(ix))
        uValP(ix) = u(cP(ix), g(ix))
    
        ! Flow values in outside option/status quo at location `ix'
        SELECT CASE (RULE)
          CASE (RULE_TAU)
            uValOutR(ix) = u(yR - tau(ix), xBar)           
            uValOutP(ix) = u(xBar, xBar)
          CASE (RULE_E)
            uValOutR(ix) = u(yR - (xBar + e(ix)), xBar)   
            uValOutP(ix) = u(yP + e(ix), xBar)
          CASE (RULE_BOTH)
            uValOutR(ix) = u(yR - tau(ix), xBar)           
            uValOutP(ix) = u(yP + e(ix), xBar)
           CASE (RULE_G)
            uValOutR(ix) = u(xbar, tau(ix) - e(ix))         
            uValOutP(ix) = u(xbar, tau(ix) - e(ix))
          CASE (RULE_ALL)
            uValOutR(ix) = u(yR - tau(ix), tau(ix) - e(ix)) 
            uValOutP(ix) = u(yP + e(ix), tau(ix) - e(ix))
          CASE DEFAULT ! No status quo
            uValOutR(ix) = u(xBar, gDiscretion)             
            uValOutP(ix) = u(xBar, gDiscretion)
        END SELECT
      END IF
    END DO
  END DO

  !
  ! Main Loop
  !
  WRITE (*, *) "Main loop..."
  iter = 1
  epsVal = 1.0_wp
  DO WHILE ( iter <= maxIter .AND. epsVal > errTol ) 
    !
    ! R in power
    !
    gotoVal = (1 - beta) * uValR + beta * ( q * Vr0 + (1 - q) * Wr0 )     ! Value to proposer R, for all possible choices
    gotoOther = (1 - beta) * uValP + beta * ( q * Wp0 + (1 - q) * Vp0 )   ! Corresponding value to other
    outsideP = (1 - beta) * uValOutP + beta * ( q * Wp0 + (1 - q) * Vp0 ) ! Outside option of other, given status quo
    !$OMP PARALLEL DO PRIVATE(ix,probs,theMax,theSum,vvalid)
    DO ix = 1,nn
      vvalid = (gotoOther >= outsideP(ix) - errTol)

      IF (.NOT. ANY(vvalid)) THEN
        WRITE (*, "(A,4ES10.3)") "R: no valid option ", MAXVAL(gotoOther), outsideP(ix), tau(ix), e(ix)
        STOP
      END IF

      theMax = MAXVAL(gotoVal, vvalid)
      WHERE (vvalid)
        probs = EXP( (gotoVal - theMax) / dchoice )
      ELSEWHERE
        probs = 0.0_wp
      END WHERE
      theSum = SUM(probs, vvalid)
      probs = probs / theSum
      polR(ix, :) = probs
      Vr1(ix) = theMax + dchoice * (LOG(theSum) - eulerMasc)
      Wp1(ix) = DOT_PRODUCT(probs, gotoOther)
    END DO
      
    !
    ! P in power
    !
    gotoVal = (1 - beta) * uValP + beta * ( q * Vp0 + (1 - q) * Wp0 )     ! Same as above, but P in power
    gotoOther = (1 - beta) * uValR + beta * ( q * Wr0 + (1 - q) * Vr0 )   ! R out of power
    outsideR = (1 - beta) * uValOutR + beta * ( q * Wr0 + (1 - q) * Vr0 ) ! R's outside option
    !$OMP PARALLEL DO PRIVATE(ix,probs,theMax,theSum,vvalid)
    DO ix = 1,nn
      vvalid = (gotoOther >= outsideR(ix) - errTol)

      IF (.NOT. ANY(vvalid)) THEN
        WRITE (*, "(A,4ES10.3)") "P: no valid option ", MAXVAL(gotoOther), outsideR(ix), tau(ix), e(ix)
        STOP
      END IF

      theMax = MAXVAL(gotoVal, vvalid)
      WHERE (vvalid)
        probs = EXP( (gotoVal - theMax) / dchoice )
      ELSEWHERE
        probs = 0.0_wp
      END WHERE
      theSum = SUM(probs, vvalid)
      probs = probs / theSum
      polP(ix, :) = probs
      Vp1(ix) = theMax + dchoice * (LOG(theSum) - eulerMasc)
      Wr1(ix) = DOT_PRODUCT(probs, gotoOther)
    END DO
    
    errwp = MAXVAL(ABS(Wp1 - Wp0))
    errvp = MAXVAL(ABS(Vp1 - Vp0))
    errwr = MAXVAL(ABS(Wr1 - Wr0))
    errvr = MAXVAL(ABS(Vr1 - Vr0))

    IF (MOD(iter, 100) == 0) WRITE (*, "(I5,4ES12.4)") iter, errwp, errvp, errwr, errvr
    FLUSH(output_unit)

    epsVal = MAXVAL((/ errwp, errvp, errwr, errvr /))
    iter = iter + 1
    Wp0 = Wp1
    Vp0 = Vp1
    Wr0 = Wr1
    Vr0 = Vr1
  END DO

  CALL tstamp()

  DEALLOCATE(uValR, uValP)
  DEALLOCATE(uValOutR, uValOutP)
  DEALLOCATE(outsideP, outsideR)
  DEALLOCATE(gotoVal, gotoOther, probs)

  !
  ! Stationary distribution
  !
  WRITE (output_unit, *) "Stationary distribution..."
  FLUSH(output_unit)
  ! Convention: top half of column R in power, bottom half P in power
  ALLOCATE(markovHuge(nn * 2, nn * 2))
  ALLOCATE(dist0(nn * 2), dist1(nn * 2))

  ! Initial condition: discretion of R (as close as possible on the grid)
  dist0 = 0.0_wp
  dist0(MINLOC( (tau - tauGrid(n/2))**2 + (e - eGrid(n/2))**2, 1)) = 1.0_wp

  !
  !
  !  R -> R  |  R -> P
  !  P -> R  |  P -> P
  !
  !
  markovHuge(1:nn, 1:nn)                   = q * polR
  markovHuge(1:nn, (nn+1):(2*nn))          = (1 - q) * polR
  markovHuge((nn+1):(2*nn), 1:nn)          = (1 - q) * polP
  markovHuge((nn+1):(2*nn), (nn+1):(2*nn)) = q * polP

  epsVal = 1.0_wp
  iter = 1
  DO WHILE ( iter <= maxIter .AND. epsVal > errTol )
    !$OMP PARALLEL DO PRIVATE(ix)
    DO ix = 1,(nn * 2)
      dist1(ix) = DOT_PRODUCT(dist0, markovHuge(:, ix))
    END DO

    epsVal = MAXVAL(ABS( dist1 - dist0 ))
    dist0 = dist1
    IF (MOD(iter, 100) == 0) WRITE (*, "(I5,E12.4)") iter, epsVal
    iter = iter + 1
  END DO

  DEALLOCATE(markovHuge, dist0)
  ALLOCATE(valTmp(nn * 2))
  ALLOCATE(polRcR(nn), polRcP(nn), polRtau(nn), polRe(nn), polRg(nn))
  ALLOCATE(polPcR(nn), polPcP(nn), polPtau(nn), polPe(nn), polPg(nn))
  polRcR = MATMUL(polR, cR)
  polRcP = MATMUL(polR, cP)
  polRtau = MATMUL(polR, tau)
  polRe = MATMUL(polR, e)
  polRg = MATMUL(polR, g)
  polPcR = MATMUL(polP, cR)
  polPcP = MATMUL(polP, cP)
  polPtau = MATMUL(polP, tau)
  polPe = MATMUL(polP, e)
  polPg = MATMUL(polP, g)

  WRITE (*, *) "Parameters:"
  WRITE (*, "(A,A20,F15.5)") tabChar, "Y: ", Y
  WRITE (*, "(A,A20,F15.5)") tabChar, "yR: ", yR
  WRITE (*, "(A,A20,F15.5)") tabChar, "yP: ", yP
  WRITE (*, "(A,A20,F15.5)") tabChar, "q: ", q
  WRITE (*, "(A,A20,F15.5)") tabChar, "beta: ", beta
  WRITE (*, "(A,A20,F15.5)") tabChar, "theta: ", theta
  WRITE (*, "(A,A20,F15.5)") tabChar, "xBar: ", xBar
  WRITE (*, *) tabChar
  WRITE (*, "(A,A20,F15.5)") tabChar, "cDiscretion: ", cDiscretion
  WRITE (*, "(A,A20,F15.5)") tabChar, "gDiscretion: ", gDiscretion
  WRITE (*, "(A,A20,F15.5)") tabChar, "tauDiscretionR: ", tauDiscretionR
  WRITE (*, "(A,A20,F15.5)") tabChar, "eDiscretionR: ", eDiscretionR
  WRITE (*, "(A,A20,F15.5)") tabChar, "tauDiscretionP: ", tauDiscretionP
  WRITE (*, "(A,A20,F15.5)") tabChar, "eDiscretionP: ", eDiscretionP
  WRITE (*, *) tabChar

  WRITE (*, *) "Ergodic:"
  WRITE (*, "(A,A20,F15.5)") tabChar, "Y: ", Y
  WRITE (*, "(A,A20,F15.5)") tabChar, "yR/Y: ", yR / Y
  WRITE (*, "(A,A20,F15.5)") tabChar, "yP/Y: ", yP / Y

  valTmp(1:nn) = polRcR
  valTmp((nn+1):(2*nn)) = polPcR
  theStat = DOT_PRODUCT( valTmp, dist1 )
  WRITE (*, "(A,A20,F15.5)") tabChar, "Mean cR/Y: ", theStat / Y
  WRITE (*, "(A,A20,F15.5)") tabChar, "Mean tau/Y: ", (yR - theStat) / Y
  valTmp(1:nn) = (polRcR - theStat)**2
  valTmp((nn+1):(2*nn)) = (polPcR - theStat)**2
  theStat = DOT_PRODUCT( valTmp, dist1 )
  WRITE (*, "(A,A20,F15.5)") tabChar, "Std cR: ", SQRT(theStat)

  valTmp(1:nn) = polRcP
  valTmp((nn+1):(2*nn)) = polPcP
  theStat = DOT_PRODUCT( valTmp, dist1 )
  WRITE (*, "(A,A20,F15.5)") tabChar, "Mean cP/Y: ", theStat / Y
  WRITE (*, "(A,A20,F15.5)") tabChar, "Mean e/Y: ", (theStat - yP) / Y
  valTmp(1:nn) = (polRcP - theStat)**2
  valTmp((nn+1):(2*nn)) = (polPcP - theStat)**2
  theStat = DOT_PRODUCT( valTmp, dist1 )
  WRITE (*, "(A,A20,F15.5)") tabChar, "Std cP: ", SQRT(theStat)

  valTmp(1:nn) = polRg
  valTmp((nn+1):(2*nn)) = polPg
  theStat = DOT_PRODUCT( valTmp, dist1 )
  WRITE (*, "(A,A20,F15.5)") tabChar, "Mean g/Y: ", theStat / Y
  valTmp(1:nn) = (polRg - theStat)**2
  valTmp((nn+1):(2*nn)) = (polPg - theStat)**2
  theStat = DOT_PRODUCT( valTmp, dist1 )
  WRITE (*, "(A,A20,F15.5)") tabChar, "Std g: ", SQRT(theStat)
  
  valTmp(1:nn) = polRe / (polRe + polRg)
  valTmp((nn+1):(2*nn)) = polPe / (polPe + polPg)
  theStat = DOT_PRODUCT( valTmp, dist1 )
  WRITE (*, "(A,A20,F15.5)") tabChar, "Mean e/(e+g): ", theStat

  WRITE (*, *) "Welfare:"
  valTmp(1:nn) = Vr0
  valTmp((nn+1):(2*nn)) = Wr0
  evalR = DOT_PRODUCT( valTmp, dist1 )
  WRITE (*, "(A,A20,F15.5)") tabChar, "Expected value R: ", evalR

  valTmp(1:nn) = Wp0
  valTmp((nn+1):(2*nn)) = Vp0
  evalP = DOT_PRODUCT( valTmp, dist1 ) 
  WRITE (*, "(A,A20,F15.5)") tabChar, "Expected value P: ", evalP 
 
  WRITE (*, *) "Saving..."
  OPEN(newunit=iunit, file=outDir // "parameters.tab")
  WRITE (iunit, "(I10)") nn
  WRITE (iunit, "(I10)") RULE
  WRITE (iunit, "(E20.10)") q 
  WRITE (iunit, "(E20.10)") beta
  WRITE (iunit, "(E20.10)") theta
  WRITE (iunit, "(E20.10)") Y
  WRITE (iunit, "(E20.10)") DeltaY
  WRITE (iunit, "(E20.10)") yR
  WRITE (iunit, "(E20.10)") yP
  WRITE (iunit, "(E20.10)") xBar 
  WRITE (iunit, "(E20.10)") errTol
  WRITE (iunit, "(E20.10)") dchoice
  WRITE (iunit, "(E20.10)") cDiscretion 
  WRITE (iunit, "(E20.10)") gDiscretion
  WRITE (iunit, "(E20.10)") tauDiscretionR
  WRITE (iunit, "(E20.10)") eDiscretionR
  WRITE (iunit, "(E20.10)") tauDiscretionP
  WRITE (iunit, "(E20.10)") eDiscretionP
  WRITE (iunit, "(E20.10)") evalR
  WRITE (iunit, "(E20.10)") evalP
  CLOSE(iunit)

!
!  Too big! Don't save/load choice probabilities unless you want to debug.
!
!  OPEN(newunit=iunit, file=outDir // "polR.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
!  WRITE (iunit) polR
!  CLOSE(iunit)
!
!  OPEN(newunit=iunit, file=outDir // "polP.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
!  WRITE (iunit) polP
!  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "tau.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) tau
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "e.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) e
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "g.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) g
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "cR.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) cR
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "cP.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) cP
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "Wp.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) Wp0
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "Vp.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) Vp0
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "Wr.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) Wr0
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "Vr.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) Vr0
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polRcR.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polRcR
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polRcP.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polRcP
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polRg.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polRg
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polRtau.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polRtau
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polRe.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polRe
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polPcR.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polPcR
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polPcP.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polPcP
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polPg.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polPg
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polPtau.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polPtau
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "polPe.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) polPe
  CLOSE(iunit)

  OPEN(newunit=iunit, file=outDir // "stationary.bin", FORM="unformatted", ACCESS="stream", STATUS="unknown")
  WRITE (iunit) dist1(1:nn)+dist1((nn+1):(2*nn))
  CLOSE(iunit)

  CALL tstamp()

!
!
!
CONTAINS
!
!
!

  PURE FUNCTION u(cc, gg)
    REAL(wp), INTENT(IN) :: cc, gg
    REAL(wp) :: u
  
    u = LOG(cc) + theta * LOG(gg)
  END FUNCTION u
 
  PURE FUNCTION linspace(startVal, endVal, noEl)
    REAL(wp), INTENT(IN) :: startVal, endVal
    INTEGER, INTENT(IN) :: noEl
    REAL(wp), DIMENSION(noEl) :: linspace
    INTEGER :: ix
    REAL(wp) :: dx

    IF (noEl == 1) THEN
      linspace = 0.5_wp * (endVal + startVal)
    ELSE
      dx = (endVal - startVal) / (noEl - 1)
      linspace = [ ( startVal + (ix - 1) * dx, ix = 1,noEl ) ]
    END IF
  END FUNCTION linspace
 
  SUBROUTINE tstamp()
    character ( len = 8 ) ampm
    integer ( kind = 4 ) d
    integer ( kind = 4 ) h
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mm
    character ( len = 9 ), parameter, dimension(12) :: month = (/ &
      'January  ', 'February ', 'March    ', 'April    ', &
      'May      ', 'June     ', 'July     ', 'August   ', &
      'September', 'October  ', 'November ', 'December ' /)
    integer ( kind = 4 ) n
    integer ( kind = 4 ) s
    integer ( kind = 4 ) values(8)
    integer ( kind = 4 ) y

    call date_and_time ( values = values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if ( h < 12 ) then
      ampm = 'AM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Noon'
      else
        ampm = 'PM'
      end if
    else
      h = h - 12
      if ( h < 12 ) then
        ampm = 'PM'
      else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
          ampm = 'Midnight'
        else
          ampm = 'AM'
        end if
      end if
    end if

    write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
      d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
    return
  END SUBROUTINE tstamp

END PROGRAM unequal
