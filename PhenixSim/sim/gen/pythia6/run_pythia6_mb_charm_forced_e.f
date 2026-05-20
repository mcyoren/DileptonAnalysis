C=======================================================================
C  run_pythia6_mb_charm_forced_e.f
C
C  PYTHIA6 minimum-bias / all-QCD Tune A
C  Accept events with open charm
C  Force weak open-charm hadrons to decay only through electron channels
C  Weight each event by product of inclusive BR(charm hadron -> e)
C  Write accepted charm events and accepted e+e- pairs
C=======================================================================

      PROGRAM RUN_PYTHIA6_MB_CHARM_FORCED_E

      IMPLICIT NONE

C     PYTHIA6 event record
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

C     PYTHIA6 process info
      INTEGER MSTI
      DOUBLE PRECISION PARI
      COMMON /PYINT1/ MSTI(200), PARI(200)

C     PYTHIA6 decay table
      INTEGER MDCY, MDME, KFDP
      DOUBLE PRECISION BRAT
      COMMON /PYDAT3/ MDCY(500,3), MDME(8000,2),
     &                BRAT(8000), KFDP(8000,5)

C     User settings
      INTEGER TARGET_ACC, MAX_GEN
      DOUBLE PRECISION SQRTS

C     Counters
      INTEGER IGEN, IACC, IPAIR
      INTEGER I, J
      INTEGER ISUB
      INTEGER NOPEN, NWEAK, NELE

C     Particle storage
      INTEGER ELE(4000)
      INTEGER PDG1, PDG2
      INTEGER MOM1, MOM2

C     Kinematics
      DOUBLE PRECISION PT, Y
      DOUBLE PRECISION PT1, PT2, Y1, Y2, ETA1, ETA2, PHI1, PHI2
      DOUBLE PRECISION PAIRPX, PAIRPY, PAIRPZ, PAIRE
      DOUBLE PRECISION PAIRPT, PAIRM, PAIRY, MASS2

C     Weights
      DOUBLE PRECISION EVENT_WEIGHT
      DOUBLE PRECISION BR

C     Function declarations
      INTEGER OPENCHARMMOTHER
      DOUBLE PRECISION GETPT, GETY, GETETA, GETPHI
      DOUBLE PRECISION FORCEELECTRONDECAYS, ELECTRONBR
      LOGICAL ISOPENCHARM, ISWEAKOPENCHARM

C=======================================================================
C     Settings
C=======================================================================

C     Number of accepted charm events to collect.
C     Accepted event = event containing at least one weak open-charm
C     hadron.
      TARGET_ACC = 100

C     Safety limit.
      MAX_GEN = 200000000

C     RHIC energy
      SQRTS = 200.0D0

C=======================================================================
C     PYTHIA6 setup
C=======================================================================

C     Minimum-bias / all QCD.
C     This allows charm from flavor creation, flavor excitation,
C     gluon splitting, parton showers, and MPI.
      CALL PYGIVE('MSEL=1')

C     Tune A-like parameters
      CALL PYGIVE('MSTP(51)=7')
      CALL PYGIVE('PARP(67)=4.0')
      CALL PYGIVE('PARP(82)=2.0')
      CALL PYGIVE('PARP(84)=0.4')
      CALL PYGIVE('PARP(85)=0.9')
      CALL PYGIVE('PARP(86)=0.95')
      CALL PYGIVE('PARP(89)=1800.0')
      CALL PYGIVE('PARP(90)=0.25')
      CALL PYGIVE('PARP(91)=1.5')

C     Minimum hard-process pT.
C     Keep 1.5 for consistency with your previous setup.
      CALL PYGIVE('CKIN(3)=1.5')

C     Less verbose PYTHIA output
      CALL PYGIVE('MSTU(11)=6')
      CALL PYGIVE('MSTU(12)=12345')

C     Initialize p+p at sqrt(s)
      CALL PYINIT('CMS','p','p',SQRTS)

C=======================================================================
C     Force weak open-charm hadrons to electron decay channels only
C=======================================================================

      WRITE(*,*) 'Forcing weak open-charm hadrons to e channels'
      WRITE(*,*) 'D+       BR_e = ', FORCEELECTRONDECAYS(411)
      WRITE(*,*) 'D0       BR_e = ', FORCEELECTRONDECAYS(421)
      WRITE(*,*) 'Ds+      BR_e = ', FORCEELECTRONDECAYS(431)
      WRITE(*,*) 'Lambda_c BR_e = ', FORCEELECTRONDECAYS(4122)
      WRITE(*,*) 'Xi_c0    BR_e = ', FORCEELECTRONDECAYS(4132)
      WRITE(*,*) 'Xi_c+    BR_e = ', FORCEELECTRONDECAYS(4232)
      WRITE(*,*) 'Omega_c0 BR_e = ', FORCEELECTRONDECAYS(4332)

      CALL PYSTAT(1)

C=======================================================================
C     Output files
C=======================================================================

C     Event-level accepted charm events
      OPEN(11, FILE='pythia6_mb_charm_events.dat',
     &     STATUS='UNKNOWN')

      WRITE(11,*) '# igen iacc isub nopen nweak nele weight'

C     Pair-level output
      OPEN(10, FILE='pythia6_mb_charm_pairs.dat',
     &     STATUS='UNKNOWN')

      WRITE(10,*) '# igen iacc ipair isub nopen nweak weight',
     &            ' pdg1 pdg2 mother1 mother2',
     &            ' pt1 pt2 y1 y2 eta1 eta2 phi1 phi2',
     &            ' pair_mass pair_pt pair_y'

C=======================================================================
C     Event loop
C=======================================================================

      IGEN = 0
      IACC = 0

 10   CONTINUE

      IF (IACC.GE.TARGET_ACC) GOTO 999
      IF (IGEN.GE.MAX_GEN) GOTO 999

      IGEN = IGEN + 1

      IF (MOD(IGEN,100000).EQ.0) THEN
        WRITE(*,*) 'Generated ', IGEN,
     &             ' accepted charm events ', IACC,
     &             ' target ', TARGET_ACC
      ENDIF

      CALL PYEVNT

      ISUB = MSTI(1)

      NOPEN = 0
      NWEAK = 0
      NELE = 0
      EVENT_WEIGHT = 1.0D0

C-----------------------------------------------------------------------
C     Scan event
C-----------------------------------------------------------------------

      DO I = 1, N

C       Count all open-charm hadrons for diagnostics
        IF (ISOPENCHARM(K(I,2))) THEN
          NOPEN = NOPEN + 1
        ENDIF

C       Count weak open-charm hadrons and build BR weight.
C       K(I,1) < 20 avoids most documentation lines.
        IF (K(I,1).LT.20 .AND. ISWEAKOPENCHARM(K(I,2))) THEN
          BR = ELECTRONBR(K(I,2))
          IF (BR.GT.0.0D0) THEN
            NWEAK = NWEAK + 1
            EVENT_WEIGHT = EVENT_WEIGHT * BR
          ENDIF
        ENDIF

C       Final-state e+ or e-
        IF (K(I,1).NE.1) GOTO 100
        IF (ABS(K(I,2)).NE.11) GOTO 100

        PT = GETPT(I)
        Y  = GETY(I)

C       PHENIX-like midrapidity electron cuts
        IF (PT.LT.0.2D0) GOTO 100
        IF (ABS(Y).GT.0.5D0) GOTO 100

C       Require open-charm ancestor
        MOM1 = OPENCHARMMOTHER(I)
        IF (MOM1.EQ.0) GOTO 100

        NELE = NELE + 1
        ELE(NELE) = I

 100    CONTINUE

      ENDDO

C-----------------------------------------------------------------------
C     Accept charm events
C-----------------------------------------------------------------------

      IF (NWEAK.LE.0) GOTO 10

      IACC = IACC + 1

      WRITE(11,9100) IGEN, IACC, ISUB, NOPEN, NWEAK, NELE,
     &                EVENT_WEIGHT

C-----------------------------------------------------------------------
C     Write unlike-sign e+e- pairs if present
C-----------------------------------------------------------------------

      IF (NELE.LT.2) GOTO 10

      IPAIR = 0

      DO I = 1, NELE
        DO J = I+1, NELE

          PDG1 = K(ELE(I),2)
          PDG2 = K(ELE(J),2)

C         Unlike-sign only
          IF (PDG1*PDG2.GT.0) GOTO 300

          IPAIR = IPAIR + 1

          MOM1 = OPENCHARMMOTHER(ELE(I))
          MOM2 = OPENCHARMMOTHER(ELE(J))

          PT1  = GETPT(ELE(I))
          PT2  = GETPT(ELE(J))
          Y1   = GETY(ELE(I))
          Y2   = GETY(ELE(J))
          ETA1 = GETETA(ELE(I))
          ETA2 = GETETA(ELE(J))
          PHI1 = GETPHI(ELE(I))
          PHI2 = GETPHI(ELE(J))

          PAIRPX = P(ELE(I),1) + P(ELE(J),1)
          PAIRPY = P(ELE(I),2) + P(ELE(J),2)
          PAIRPZ = P(ELE(I),3) + P(ELE(J),3)
          PAIRE  = P(ELE(I),4) + P(ELE(J),4)

          PAIRPT = SQRT(PAIRPX*PAIRPX + PAIRPY*PAIRPY)

          MASS2 = PAIRE*PAIRE
     &          - PAIRPX*PAIRPX
     &          - PAIRPY*PAIRPY
     &          - PAIRPZ*PAIRPZ

          IF (MASS2.GT.0.0D0) THEN
            PAIRM = SQRT(MASS2)
          ELSE
            PAIRM = 0.0D0
          ENDIF

          IF ((PAIRE-PAIRPZ).GT.0.0D0 .AND.
     &        (PAIRE+PAIRPZ).GT.0.0D0) THEN
            PAIRY = 0.5D0*LOG((PAIRE+PAIRPZ)/(PAIRE-PAIRPZ))
          ELSE
            PAIRY = 999.0D0
          ENDIF

          WRITE(10,9000) IGEN, IACC, IPAIR, ISUB,
     &      NOPEN, NWEAK, EVENT_WEIGHT,
     &      PDG1, PDG2, MOM1, MOM2,
     &      PT1, PT2, Y1, Y2, ETA1, ETA2, PHI1, PHI2,
     &      PAIRM, PAIRPT, PAIRY

 300      CONTINUE

        ENDDO
      ENDDO

      GOTO 10

C=======================================================================
C     End
C=======================================================================

 999  CONTINUE

 9000 FORMAT(I10,1X,I10,1X,I6,1X,I6,1X,
     &       I6,1X,I6,1X,E16.8,1X,
     &       I6,1X,I6,1X,I8,1X,I8,1X,
     &       11(E16.8,1X))

 9100 FORMAT(I10,1X,I10,1X,I6,1X,I6,1X,
     &       I6,1X,I6,1X,E16.8)

      CLOSE(10)
      CLOSE(11)

      WRITE(*,*) 'Done.'
      WRITE(*,*) 'Generated events: ', IGEN
      WRITE(*,*) 'Accepted charm events: ', IACC
      WRITE(*,*) 'Event output: pythia6_mb_charm_events.dat'
      WRITE(*,*) 'Pair output:  pythia6_mb_charm_pairs.dat'

      END


C=======================================================================
C     Find open-charm mother of particle IPART
C=======================================================================

      INTEGER FUNCTION OPENCHARMMOTHER(IPART)

      IMPLICIT NONE

      INTEGER IPART
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

      INTEGER CUR, MOM, DEPTH, PDG
      LOGICAL ISOPENCHARM

      CUR = IPART
      OPENCHARMMOTHER = 0

      DO DEPTH = 1, 50

        IF (CUR.LE.0 .OR. CUR.GT.N) RETURN

        PDG = K(CUR,2)

        IF (ISOPENCHARM(PDG)) THEN
          OPENCHARMMOTHER = PDG
          RETURN
        ENDIF

        MOM = K(CUR,3)

        IF (MOM.LE.0 .OR. MOM.GT.N) RETURN
        IF (MOM.EQ.CUR) RETURN

        CUR = MOM

      ENDDO

      RETURN
      END


C=======================================================================
C     Open charm hadrons, including excited charm for ancestry tracing
C=======================================================================

      LOGICAL FUNCTION ISOPENCHARM(PDG)

      IMPLICIT NONE

      INTEGER PDG, APDG

      APDG = ABS(PDG)

      ISOPENCHARM = .FALSE.

C     Weak D mesons
      IF (APDG.EQ.411)  ISOPENCHARM = .TRUE.
      IF (APDG.EQ.421)  ISOPENCHARM = .TRUE.
      IF (APDG.EQ.431)  ISOPENCHARM = .TRUE.

C     Excited D mesons
      IF (APDG.EQ.413)  ISOPENCHARM = .TRUE.
      IF (APDG.EQ.423)  ISOPENCHARM = .TRUE.
      IF (APDG.EQ.433)  ISOPENCHARM = .TRUE.

C     Weak charm baryons
      IF (APDG.EQ.4122) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4132) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4232) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4332) ISOPENCHARM = .TRUE.

C     Excited/light charm baryons often decay strongly to weak charm
      IF (APDG.EQ.4112) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4212) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4222) ISOPENCHARM = .TRUE.

      RETURN
      END


C=======================================================================
C     Weakly decaying open-charm hadrons
C     These are the ones we force to electron channels and weight by
C     BR_e.
C=======================================================================

      LOGICAL FUNCTION ISWEAKOPENCHARM(PDG)

      IMPLICIT NONE

      INTEGER PDG, APDG

      APDG = ABS(PDG)

      ISWEAKOPENCHARM = .FALSE.

      IF (APDG.EQ.411)  ISWEAKOPENCHARM = .TRUE.
      IF (APDG.EQ.421)  ISWEAKOPENCHARM = .TRUE.
      IF (APDG.EQ.431)  ISWEAKOPENCHARM = .TRUE.
      IF (APDG.EQ.4122) ISWEAKOPENCHARM = .TRUE.
      IF (APDG.EQ.4132) ISWEAKOPENCHARM = .TRUE.
      IF (APDG.EQ.4232) ISWEAKOPENCHARM = .TRUE.
      IF (APDG.EQ.4332) ISWEAKOPENCHARM = .TRUE.

      RETURN
      END


C=======================================================================
C     Force particle and antiparticle to electron decay channels
C=======================================================================

      DOUBLE PRECISION FUNCTION FORCEELECTRONDECAYS(KF)

      IMPLICIT NONE

      INTEGER KF
      DOUBLE PRECISION FORCEONEELECTRONDECAY
      DOUBLE PRECISION BR1, BR2

      BR1 = FORCEONEELECTRONDECAY(KF)
      BR2 = FORCEONEELECTRONDECAY(-KF)

      FORCEELECTRONDECAYS = BR1

      RETURN
      END


      DOUBLE PRECISION FUNCTION FORCEONEELECTRONDECAY(KF)

      IMPLICIT NONE

      INTEGER KF

      INTEGER MDCY, MDME, KFDP
      DOUBLE PRECISION BRAT
      COMMON /PYDAT3/ MDCY(500,3), MDME(8000,2),
     &                BRAT(8000), KFDP(8000,5)

      INTEGER KC, PYCOMP
      INTEGER IDC1, IDC2, IDC, J
      LOGICAL HASELECTRON

      KC = PYCOMP(KF)

      FORCEONEELECTRONDECAY = 0.0D0

      IF (KC.LE.0) RETURN
      IF (MDCY(KC,1).EQ.0) RETURN

      IDC1 = MDCY(KC,2)
      IDC2 = MDCY(KC,2) + MDCY(KC,3) - 1

      DO IDC = IDC1, IDC2

        HASELECTRON = .FALSE.

        DO J = 1, 5
          IF (ABS(KFDP(IDC,J)).EQ.11) HASELECTRON = .TRUE.
        ENDDO

        IF (HASELECTRON) THEN
          FORCEONEELECTRONDECAY = FORCEONEELECTRONDECAY
     &                            + BRAT(IDC)
          MDME(IDC,1) = 1
        ELSE
          MDME(IDC,1) = 0
        ENDIF

      ENDDO

      RETURN
      END


C=======================================================================
C     Return inclusive electron branching ratio from PYTHIA table
C=======================================================================

      DOUBLE PRECISION FUNCTION ELECTRONBR(KF)

      IMPLICIT NONE

      INTEGER KF

      INTEGER MDCY, MDME, KFDP
      DOUBLE PRECISION BRAT
      COMMON /PYDAT3/ MDCY(500,3), MDME(8000,2),
     &                BRAT(8000), KFDP(8000,5)

      INTEGER KC, PYCOMP
      INTEGER IDC1, IDC2, IDC, J
      LOGICAL HASELECTRON
      LOGICAL ISWEAKOPENCHARM

      ELECTRONBR = 1.0D0

      IF (.NOT.ISWEAKOPENCHARM(KF)) RETURN

      ELECTRONBR = 0.0D0

      KC = PYCOMP(KF)

      IF (KC.LE.0) RETURN
      IF (MDCY(KC,1).EQ.0) RETURN

      IDC1 = MDCY(KC,2)
      IDC2 = MDCY(KC,2) + MDCY(KC,3) - 1

      DO IDC = IDC1, IDC2

        HASELECTRON = .FALSE.

        DO J = 1, 5
          IF (ABS(KFDP(IDC,J)).EQ.11) HASELECTRON = .TRUE.
        ENDDO

        IF (HASELECTRON) THEN
          ELECTRONBR = ELECTRONBR + BRAT(IDC)
        ENDIF

      ENDDO

      RETURN
      END


C=======================================================================
C     Kinematic helper functions
C=======================================================================

      DOUBLE PRECISION FUNCTION GETPT(I)

      IMPLICIT NONE

      INTEGER I
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

      GETPT = SQRT(P(I,1)*P(I,1) + P(I,2)*P(I,2))

      RETURN
      END


      DOUBLE PRECISION FUNCTION GETY(I)

      IMPLICIT NONE

      INTEGER I
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

      IF ((P(I,4)-P(I,3)).GT.0.0D0 .AND.
     &    (P(I,4)+P(I,3)).GT.0.0D0) THEN
        GETY = 0.5D0*LOG((P(I,4)+P(I,3))/
     &                    (P(I,4)-P(I,3)))
      ELSE
        GETY = 999.0D0
      ENDIF

      RETURN
      END


      DOUBLE PRECISION FUNCTION GETETA(I)

      IMPLICIT NONE

      INTEGER I
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

      DOUBLE PRECISION PP

      PP = SQRT(P(I,1)*P(I,1) + P(I,2)*P(I,2)
     &        + P(I,3)*P(I,3))

      IF ((PP-P(I,3)).GT.0.0D0 .AND.
     &    (PP+P(I,3)).GT.0.0D0) THEN
        GETETA = 0.5D0*LOG((PP+P(I,3))/(PP-P(I,3)))
      ELSE
        GETETA = 999.0D0
      ENDIF

      RETURN
      END


      DOUBLE PRECISION FUNCTION GETPHI(I)

      IMPLICIT NONE

      INTEGER I
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

      GETPHI = ATAN2(P(I,2), P(I,1))

      RETURN
      END
