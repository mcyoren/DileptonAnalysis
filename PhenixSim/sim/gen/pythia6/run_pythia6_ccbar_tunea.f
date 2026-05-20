      PROGRAM RUN_PYTHIA6_CCBAR_TUNEA

      IMPLICIT NONE

C     PYTHIA6 event record
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

      INTEGER NEV, IEV, I, J
      INTEGER PDG1, PDG2, MOM1, MOM2
      INTEGER ELE(4000), NELE
      DOUBLE PRECISION SQRTS
      DOUBLE PRECISION PT, Y
      DOUBLE PRECISION PT1, PT2, Y1, Y2, ETA1, ETA2, PHI1, PHI2
      DOUBLE PRECISION PX, PY, PZ, EE
      DOUBLE PRECISION PAIRPX, PAIRPY, PAIRPZ, PAIRE
      DOUBLE PRECISION PAIRPT, PAIRM, PAIRY
      DOUBLE PRECISION MASS2

      INTEGER OPENCHARMMOTHER
      DOUBLE PRECISION GETPT, GETY, GETETA, GETPHI

      NEV = 10000
      SQRTS = 200.0D0

C     ------------------------------------------------------------
C     PYTHIA6 setup
C     ------------------------------------------------------------

C     Dedicated charm production.
C     MSEL=4 is ccbar production.
C     The table you showed has MSEL=1 for all QCD, but for ccbar use 4.
      CALL PYGIVE('MSEL=4')

C     Tune A parameters from your table
      CALL PYGIVE('MSTP(51)=7')
      CALL PYGIVE('PARP(67)=4.0')
      CALL PYGIVE('PARP(82)=2.0')
      CALL PYGIVE('PARP(84)=0.4')
      CALL PYGIVE('PARP(85)=0.9')
      CALL PYGIVE('PARP(86)=0.95')
      CALL PYGIVE('PARP(89)=1800.0')
      CALL PYGIVE('PARP(90)=0.25')
      CALL PYGIVE('PARP(91)=1.5')
      CALL PYGIVE('CKIN(3)=1.5')

C     Less verbose output
      CALL PYGIVE('MSTU(11)=6')
      CALL PYGIVE('MSTU(12)=12345')

C     Initialize pp at sqrt(s)
      CALL PYINIT('CMS','p','p',SQRTS)

C     Print PYTHIA settings summary
      CALL PYSTAT(1)

C     ------------------------------------------------------------
C     Output file
C     ------------------------------------------------------------

      OPEN(10, FILE='pythia6_ccbar_tunea_dielectrons.dat',
     &     STATUS='UNKNOWN')

      WRITE(10,*) '# event pdg1 pdg2 mother1 mother2',
     &            ' pt1 pt2 y1 y2 eta1 eta2 phi1 phi2',
     &            ' pair_mass pair_pt pair_y'

C     ------------------------------------------------------------
C     Event loop
C     ------------------------------------------------------------

      DO IEV = 1, NEV

        IF (MOD(IEV,1000).EQ.0) THEN
          WRITE(*,*) 'Event ', IEV, ' / ', NEV
        ENDIF

        CALL PYEVNT

        NELE = 0

C       Find final-state electrons/positrons from open charm
        DO I = 1, N

          IF (K(I,1).NE.1) GOTO 100
          IF (ABS(K(I,2)).NE.11) GOTO 100

          PT = GETPT(I)
          Y  = GETY(I)

          IF (PT.LT.0.2D0) GOTO 100
          IF (ABS(Y).GT.0.5D0) GOTO 100

          MOM1 = OPENCHARMMOTHER(I)
          IF (MOM1.EQ.0) GOTO 100

          NELE = NELE + 1
          ELE(NELE) = I

 100      CONTINUE
        ENDDO

C       Build unlike-sign e+e- pairs
        DO I = 1, NELE
          DO J = I+1, NELE

            PDG1 = K(ELE(I),2)
            PDG2 = K(ELE(J),2)

C           Unlike sign: e- has PDG 11, e+ has PDG -11
            IF (PDG1*PDG2.GT.0) GOTO 200

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
     &            - PAIRPX*PAIRPX
     &            - PAIRPY*PAIRPY
     &            - PAIRPZ*PAIRPZ

            IF (MASS2.GT.0.0D0) THEN
              PAIRM = SQRT(MASS2)
            ELSE
              PAIRM = 0.0D0
            ENDIF

            IF ((PAIRE-PAIRPZ).GT.0.0D0 .AND.
     &          (PAIRE+PAIRPZ).GT.0.0D0) THEN
              PAIRY = 0.5D0*LOG((PAIRE+PAIRPZ)/(PAIRE-PAIRPZ))
            ELSE
              PAIRY = 999.0D0
            ENDIF

            WRITE(10,9000) IEV, PDG1, PDG2, MOM1, MOM2,
     &        PT1, PT2, Y1, Y2, ETA1, ETA2, PHI1, PHI2,
     &        PAIRM, PAIRPT, PAIRY

 200        CONTINUE
          ENDDO
        ENDDO

      ENDDO

 9000 FORMAT(I8,1X,I6,1X,I6,1X,I8,1X,I8,1X,
     &       11(E16.8,1X))

      CLOSE(10)

      WRITE(*,*) 'Done.'
      WRITE(*,*) 'Output: pythia6_ccbar_tunea_dielectrons.dat'

      END


C ================================================================
C Return open charm mother PDG code for particle I.
C Returns 0 if no open charm mother is found.
C ================================================================

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


C ================================================================
C Open charm hadrons.
C This excludes charmonia like J/psi.
C ================================================================

      LOGICAL FUNCTION ISOPENCHARM(PDG)

      IMPLICIT NONE

      INTEGER PDG, APDG

      APDG = ABS(PDG)

      ISOPENCHARM = .FALSE.

C     D mesons
      IF (APDG.EQ.411)  ISOPENCHARM = .TRUE.
      IF (APDG.EQ.421)  ISOPENCHARM = .TRUE.
      IF (APDG.EQ.431)  ISOPENCHARM = .TRUE.

C     Charm baryons
      IF (APDG.EQ.4122) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4112) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4212) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4222) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4132) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4232) ISOPENCHARM = .TRUE.
      IF (APDG.EQ.4332) ISOPENCHARM = .TRUE.

      RETURN
      END


C ================================================================
C Kinematic helper functions
C ================================================================

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
        GETY = 0.5D0*LOG((P(I,4)+P(I,3))/(P(I,4)-P(I,3)))
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
