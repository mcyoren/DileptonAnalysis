C=======================================================================
C run_pythia6_mb_ccbar_dielectrons_forced_QA.f
C
C PYTHIA6 MB/all-QCD Tune A
C Force weak open-charm hadrons to electron decay channels
C
C Production output:
C   keep only PHENIX-perfect accepted c-cbar dielectron pairs:
C     OS e+e-, each from weak open charm, opposite-sign charm parents
C     |y_e| < 0.5, pT_e > 0.2 GeV
C   stop when TARGET_PAIRS production pairs are reached
C
C QA output:
C   write ALL open-HF OS dielectron pairs before acceptance cuts
C   write ALL weak open-charm hadrons before pair acceptance cuts
C
C Outputs:
C   pythia6_pairs.dat              production pair-level file
C   pythia6_tracks.dat             production track-level file
C   pythia6_summary.dat            normalization summary
C   pythia6_qa_pairs.dat           all-pair QA flat file
C   pythia6_qa_charmhadrons.dat    all weak open-charm hadron QA file
C=======================================================================

      PROGRAM RUN_PYTHIA6_MB_CCBAR_DIELECTRONS_QA

      IMPLICIT NONE

C----- PYTHIA event record
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

C----- PYTHIA process info
      INTEGER MSTI
      DOUBLE PRECISION PARI
      COMMON /PYINT1/ MSTI(200), PARI(200)

C----- PYTHIA cross sections
      INTEGER NGEN, NXSEC
      DOUBLE PRECISION XSEC
      COMMON /PYINT5/ NGEN(0:500,3), XSEC(0:500,3), NXSEC(0:500)

C----- PYTHIA decay table
      INTEGER MDCY, MDME, KFDP
      DOUBLE PRECISION BRAT
      COMMON /PYDAT3/ MDCY(500,3), MDME(8000,2),
     &                BRAT(8000), KFDP(8000,5)

C----- PYTHIA random number state
      INTEGER MRPY
      DOUBLE PRECISION RRPY
      COMMON /PYDATR/ MRPY(6), RRPY(100)

C----- Settings
      INTEGER TARGET_PAIRS, MAX_GEN
      DOUBLE PRECISION SQRTS

C----- Condor / command line
      INTEGER NARG, IARGC
      INTEGER JOBID
      CHARACTER*128 ARG

C----- Counters
      INTEGER IGEN, IACC_PAIR, IQA_PAIR
      INTEGER I, J
      INTEGER ISUB, SRCBIN
      INTEGER NELE
      INTEGER ELE(4000)

C----- Particle and parent variables
      INTEGER PDG1, PDG2
      INTEGER PARENT1, PARENT2
      INTEGER MOMIDX1, MOMIDX2
      INTEGER SPECIES
      DOUBLE PRECISION BR1, BR2, BRI
      DOUBLE PRECISION WPAIR, WREL
      DOUBLE PRECISION BRD0, BRD0SQ

C----- Electron and pair kinematics
      DOUBLE PRECISION PT1, PT2, Y1, Y2, ETA1, ETA2, PHI1, PHI2
      DOUBLE PRECISION PAIRPX, PAIRPY, PAIRPZ, PAIRE, PAIRP
      DOUBLE PRECISION PAIRPT, PAIRM, PAIRY, PAIRETA, MASS2
      DOUBLE PRECISION PTMOM1, PTMOM2, YMOM1, YMOM2

C----- Acceptance flags
      INTEGER PASS_PHENIX, PASS_STAR

C----- Weights and normalization sums
      DOUBLE PRECISION SUM_WPAIR, SUM_WREL
      DOUBLE PRECISION SUM_QA_WPAIR, SUM_QA_WREL

C----- Functions
      INTEGER OPENCHARMMOTHERINDEX, GETPROCESSBIN, GETPARENTBIN
      DOUBLE PRECISION GETPT, GETY, GETETA, GETPHI, GETMASS
      DOUBLE PRECISION FORCEELECTRONDECAYS, ELECTRONBR
      LOGICAL ISWEAKOPENCHARM

C=======================================================================
C User settings and command-line arguments
C=======================================================================

      TARGET_PAIRS = 100
      MAX_GEN      = 200000000
      SQRTS        = 200.0D0
      JOBID        = 0

C     Arguments:
C       arg1 = target number of PHENIX-perfect accepted dielectron pairs
C       arg2 = job id / seed / Condor process number
      NARG = IARGC()

      IF (NARG.GE.1) THEN
        CALL GETARG(1, ARG)
        READ(ARG,*) TARGET_PAIRS
      ENDIF

      IF (NARG.GE.2) THEN
        CALL GETARG(2, ARG)
        READ(ARG,*) JOBID
      ENDIF

      WRITE(*,*) 'TARGET_PAIRS = ', TARGET_PAIRS
      WRITE(*,*) 'JOBID        = ', JOBID

      SUM_WPAIR = 0.0D0
      SUM_WREL  = 0.0D0
      SUM_QA_WPAIR = 0.0D0
      SUM_QA_WREL  = 0.0D0

C=======================================================================
C PYTHIA6 setup: MB/all-QCD Tune A
C=======================================================================

      CALL PYGIVE('MSEL=1')

C----- Tune A parameters PARP(67) is 4 and PARP(91)=1.5  for phenix
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

C----- Less verbose
      CALL PYGIVE('MSTU(11)=6')
      CALL PYGIVE('MSTU(12)=12345')

C----- Random seed. PYTHIA6 uses MRPY(1).
      MRPY(1) = 100000 + 1000*JOBID
      MRPY(2) = 0

      CALL PYINIT('CMS','p','p',SQRTS)

C=======================================================================
C Force weak open-charm hadrons to electron channels
C=======================================================================

      WRITE(*,*) 'Forcing weak open-charm hadrons to e channels'
      WRITE(*,*) 'D+       BR_e = ', FORCEELECTRONDECAYS(411)
      WRITE(*,*) 'D0       BR_e = ', FORCEELECTRONDECAYS(421)
      WRITE(*,*) 'Ds+      BR_e = ', FORCEELECTRONDECAYS(431)
      WRITE(*,*) 'Lambda_c BR_e = ', FORCEELECTRONDECAYS(4122)
      WRITE(*,*) 'Xi_c0    BR_e = ', FORCEELECTRONDECAYS(4132)
      WRITE(*,*) 'Xi_c+    BR_e = ', FORCEELECTRONDECAYS(4232)
      WRITE(*,*) 'Omega_c0 BR_e = ', FORCEELECTRONDECAYS(4332)

      BRD0 = ELECTRONBR(421)
      BRD0SQ = BRD0 * BRD0

      WRITE(*,*) 'Reference BR(D0->e) = ', BRD0
      WRITE(*,*) 'Reference BR(D0)^2  = ', BRD0SQ

      IF (BRD0.LE.0.0D0) THEN
        WRITE(*,*) 'ERROR: BRD0 <= 0'
        STOP
      ENDIF

      CALL PYSTAT(1)

C=======================================================================
C Output files
C=======================================================================

      OPEN(10, FILE='pythia6_pairs.dat', STATUS='UNKNOWN')
      OPEN(11, FILE='pythia6_tracks.dat', STATUS='UNKNOWN')
      OPEN(12, FILE='pythia6_summary.dat', STATUS='UNKNOWN')
      OPEN(13, FILE='pythia6_qa_pairs.dat', STATUS='UNKNOWN')
      OPEN(14, FILE='pythia6_qa_charmhadrons.dat', STATUS='UNKNOWN')
      OPEN(15, FILE='pythia6_qa_electrons.dat', STATUS='UNKNOWN')

      WRITE(10,*) '# pair_id gen_event isub srcbin weight_br',
     &            ' weight_rel parent1 parent2 ptmom1 ptmom2',
     &            ' ymom1 ymom2 pdg1 pdg2',
     &            ' pt1 pt2 y1 y2 eta1 eta2 phi1 phi2',
     &            ' pair_mass pair_pt pair_y pair_eta'

      WRITE(11,*) '# pair_id itrack pid parent mass energy',
     &            ' px py pz vx vy vz'

      WRITE(13,*) '# qa_pair_id gen_event isub srcbin weight_br',
     &            ' weight_rel parent1 parent2 ptmom1 ptmom2',
     &            ' ymom1 ymom2 pdg1 pdg2',
     &            ' pt1 pt2 y1 y2 eta1 eta2 phi1 phi2',
     &            ' pair_mass pair_pt pair_y pair_eta',
     &            ' pass_phenix pass_star'

      WRITE(14,*) '# gen_event isub srcbin pdg species br_e',
     &            ' pt y eta phi px py pz e'
      WRITE(15,*) '# gen_event isub srcbin pid parent species br_e',
     &            ' pt y eta phi px py pz e'

C=======================================================================
C Event loop
C=======================================================================

      IGEN = 0
      IACC_PAIR = 0
      IQA_PAIR = 0

 1000 CONTINUE

      IF (IACC_PAIR.GE.TARGET_PAIRS) GOTO 9999
      IF (IGEN.GE.MAX_GEN) GOTO 9999

      IGEN = IGEN + 1

      IF (MOD(IGEN,10000).EQ.0) THEN
        WRITE(*,*) 'Generated ', IGEN,
     &             ' PHENIX pairs ', IACC_PAIR,
     &             ' QA pairs ', IQA_PAIR,
     &             ' target ', TARGET_PAIRS
      ENDIF

      CALL PYEVNT

      ISUB = MSTI(1)
      SRCBIN = GETPROCESSBIN(ISUB)

C----- Write weak open-charm hadrons for unbiased mother pT QA.
      DO I = 1, N
        IF (K(I,1).LT.20 .AND. ISWEAKOPENCHARM(K(I,2))) THEN
          BRI = ELECTRONBR(K(I,2))
          SPECIES = GETPARENTBIN(K(I,2))
          WRITE(14,9400) IGEN, ISUB, SRCBIN, K(I,2), SPECIES, BRI,
     &      GETPT(I), GETY(I), GETETA(I), GETPHI(I),
     &      P(I,1), P(I,2), P(I,3), P(I,4)
        ENDIF
      ENDDO

C----- Collect all final-state electrons from weak open-charm parents.
C     No PHENIX/STAR cuts here: those are QA categories.
      NELE = 0

      DO I = 1, N

        IF (K(I,1).NE.1) GOTO 110
        IF (ABS(K(I,2)).NE.11) GOTO 110

        MOMIDX1 = OPENCHARMMOTHERINDEX(I)
        IF (MOMIDX1.LE.0) GOTO 110

        PARENT1 = K(MOMIDX1,2)
        IF (.NOT.ISWEAKOPENCHARM(PARENT1)) GOTO 110

        BRI = ELECTRONBR(PARENT1)
        SPECIES = GETPARENTBIN(PARENT1)

        WRITE(15,9500) IGEN, ISUB, SRCBIN,
     &      K(I,2), PARENT1, SPECIES, BRI,
     &      GETPT(I), GETY(I), GETETA(I), GETPHI(I),
     &      P(I,1), P(I,2), P(I,3), P(I,4)

        NELE = NELE + 1
        ELE(NELE) = I

 110    CONTINUE
      ENDDO

      IF (NELE.LT.2) GOTO 1000

C----- Build all unlike-sign c-cbar dielectron pairs for QA.
      DO I = 1, NELE
        DO J = I+1, NELE

          PDG1 = K(ELE(I),2)
          PDG2 = K(ELE(J),2)

C--------- unlike-sign e+e-
          IF (PDG1*PDG2.GT.0) GOTO 220

          MOMIDX1 = OPENCHARMMOTHERINDEX(ELE(I))
          MOMIDX2 = OPENCHARMMOTHERINDEX(ELE(J))

          IF (MOMIDX1.LE.0 .OR. MOMIDX2.LE.0) GOTO 220

          PARENT1 = K(MOMIDX1,2)
          PARENT2 = K(MOMIDX2,2)

C--------- require one charm hadron and one anti-charm hadron
          IF (PARENT1*PARENT2.GT.0) GOTO 220

          BR1 = ELECTRONBR(PARENT1)
          BR2 = ELECTRONBR(PARENT2)

          IF (BR1.LE.0.0D0 .OR. BR2.LE.0.0D0) GOTO 220

          WPAIR = BR1 * BR2
          WREL  = WPAIR / BRD0SQ

          PT1  = GETPT(ELE(I))
          PT2  = GETPT(ELE(J))
          Y1   = GETY(ELE(I))
          Y2   = GETY(ELE(J))
          ETA1 = GETETA(ELE(I))
          ETA2 = GETETA(ELE(J))
          PHI1 = GETPHI(ELE(I))
          PHI2 = GETPHI(ELE(J))

          PTMOM1 = GETPT(MOMIDX1)
          PTMOM2 = GETPT(MOMIDX2)
          YMOM1  = GETY(MOMIDX1)
          YMOM2  = GETY(MOMIDX2)

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

          PAIRP = SQRT(PAIRPX*PAIRPX + PAIRPY*PAIRPY
     &          + PAIRPZ*PAIRPZ)
          IF ((PAIRP-PAIRPZ).GT.0.0D0 .AND.
     &        (PAIRP+PAIRPZ).GT.0.0D0) THEN
            PAIRETA = 0.5D0*LOG((PAIRP+PAIRPZ)/(PAIRP-PAIRPZ))
          ELSE
            PAIRETA = 999.0D0
          ENDIF

          PASS_PHENIX = 0
          IF (PT1.GT.0.2D0 .AND. PT2.GT.0.2D0 .AND.
     &        ABS(Y1).LT.0.5D0 .AND. ABS(Y2).LT.0.5D0) THEN
            PASS_PHENIX = 1
          ENDIF

          PASS_STAR = 0
          IF (PT1.GT.0.2D0 .AND. PT2.GT.0.2D0 .AND.
     &        ABS(ETA1).LT.1.0D0 .AND. ABS(ETA2).LT.1.0D0) THEN
            PASS_STAR = 1
          ENDIF

C--------- QA pair output before final production selection.
          IQA_PAIR = IQA_PAIR + 1
          SUM_QA_WPAIR = SUM_QA_WPAIR + WPAIR
          SUM_QA_WREL  = SUM_QA_WREL  + WREL

          WRITE(13,9300) IQA_PAIR, IGEN, ISUB, SRCBIN,
     &      WPAIR, WREL,
     &      PARENT1, PARENT2, PTMOM1, PTMOM2, YMOM1, YMOM2,
     &      PDG1, PDG2,
     &      PT1, PT2, Y1, Y2, ETA1, ETA2, PHI1, PHI2,
     &      PAIRM, PAIRPT, PAIRY, PAIRETA,
     &      PASS_PHENIX, PASS_STAR

C--------- Production output only for PHENIX-perfect pairs.
          IF (PASS_PHENIX.NE.1) GOTO 220

          IACC_PAIR = IACC_PAIR + 1
          SUM_WPAIR = SUM_WPAIR + WPAIR
          SUM_WREL  = SUM_WREL  + WREL

          WRITE(10,9000) IACC_PAIR, IGEN, ISUB, SRCBIN,
     &      WPAIR, WREL,
     &      PARENT1, PARENT2, PTMOM1, PTMOM2, YMOM1, YMOM2,
     &      PDG1, PDG2,
     &      PT1, PT2, Y1, Y2, ETA1, ETA2, PHI1, PHI2,
     &      PAIRM, PAIRPT, PAIRY, PAIRETA

C--------- track-level production output: track 1
          WRITE(11,9100) IACC_PAIR, 1,
     &      K(ELE(I),2), PARENT1,
     &      GETMASS(ELE(I)), P(ELE(I),4),
     &      P(ELE(I),1), P(ELE(I),2), P(ELE(I),3),
     &      V(ELE(I),1), V(ELE(I),2), V(ELE(I),3)

C--------- track-level production output: track 2
          WRITE(11,9100) IACC_PAIR, 2,
     &      K(ELE(J),2), PARENT2,
     &      GETMASS(ELE(J)), P(ELE(J),4),
     &      P(ELE(J),1), P(ELE(J),2), P(ELE(J),3),
     &      V(ELE(J),1), V(ELE(J),2), V(ELE(J),3)

          IF (IACC_PAIR.GE.TARGET_PAIRS) GOTO 9999

 220      CONTINUE
        ENDDO
      ENDDO

      GOTO 1000

C=======================================================================
C Finish
C=======================================================================

 9999 CONTINUE

 9000 FORMAT(I12,1X,I12,1X,I6,1X,I6,1X,
     &       E16.8,1X,E16.8,1X,
     &       I8,1X,I8,1X,
     &       E16.8,1X,E16.8,1X,E16.8,1X,E16.8,1X,
     &       I6,1X,I6,1X,
     &       12(E16.8,1X))

 9100 FORMAT(I12,1X,I4,1X,I8,1X,I8,1X,
     &       8(E16.8,1X))

 9300 FORMAT(I12,1X,I12,1X,I6,1X,I6,1X,
     &       E16.8,1X,E16.8,1X,
     &       I8,1X,I8,1X,
     &       E16.8,1X,E16.8,1X,E16.8,1X,E16.8,1X,
     &       I6,1X,I6,1X,
     &       12(E16.8,1X),I3,1X,I3)

 9400 FORMAT(I12,1X,I6,1X,I6,1X,I8,1X,I3,1X,
     &       9(E16.8,1X))
 9500 FORMAT(I12,1X,I6,1X,I6,1X,I8,1X,I8,1X,I3,1X,
     &       9(E16.8,1X))
      CALL PYSTAT(1)

C----- Summary / normalization info
      WRITE(12,*) '# PYTHIA6 MB ccbar dielectron forced-e summary'
      WRITE(12,*) 'sqrt_s_GeV ', SQRTS
      WRITE(12,*) 'job_id ', JOBID
      WRITE(12,*) 'target_phenix_pairs ', TARGET_PAIRS
      WRITE(12,*) 'generated_events ', IGEN
      WRITE(12,*) 'accepted_phenix_pairs ', IACC_PAIR
      WRITE(12,*) 'qa_all_pairs ', IQA_PAIR
      WRITE(12,*) 'sum_phenix_pair_BR_weights ', SUM_WPAIR
      WRITE(12,*) 'sum_phenix_relative_rep_weights ', SUM_WREL
      WRITE(12,*) 'sum_qa_pair_BR_weights ', SUM_QA_WPAIR
      WRITE(12,*) 'sum_qa_relative_rep_weights ', SUM_QA_WREL
      WRITE(12,*) 'BRD0_to_e ', BRD0
      WRITE(12,*) 'BRD0_to_e_squared ', BRD0SQ
      WRITE(12,*) 'pythia_XSEC_0_3_mb ', XSEC(0,3)
      WRITE(12,*) 'pythia_XSEC_96_3_mb ', XSEC(96,3)
      WRITE(12,*) 'estimated_phenix_pair_cross_section_mb ',
     &             XSEC(0,3) * SUM_WPAIR / DBLE(IGEN)
      WRITE(12,*) 'estimated_all_qa_pair_cross_section_mb ',
     &             XSEC(0,3) * SUM_QA_WPAIR / DBLE(IGEN)

      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)

      WRITE(*,*) 'Done.'
      WRITE(*,*) 'Generated events: ', IGEN
      WRITE(*,*) 'Accepted PHENIX pairs: ', IACC_PAIR
      WRITE(*,*) 'QA all pairs: ', IQA_PAIR
      WRITE(*,*) 'Sum PHENIX BR weights: ', SUM_WPAIR
      WRITE(*,*) 'Sum QA BR weights: ', SUM_QA_WPAIR
      WRITE(*,*) 'XSEC(0,3) mb: ', XSEC(0,3)
      WRITE(*,*) 'PHENIX pair xsec mb: ',
     &             XSEC(0,3)*SUM_WPAIR/DBLE(IGEN)
      WRITE(*,*) 'All QA pair xsec mb: ',
     &             XSEC(0,3)*SUM_QA_WPAIR/DBLE(IGEN)

      END


C=======================================================================
C Find weak open-charm mother index of particle IPART
C=======================================================================

      INTEGER FUNCTION OPENCHARMMOTHERINDEX(IPART)

      IMPLICIT NONE

      INTEGER IPART
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

      INTEGER CUR, MOM, DEPTH, PDG
      LOGICAL ISWEAKOPENCHARM

      CUR = IPART
      OPENCHARMMOTHERINDEX = 0

      DO DEPTH = 1, 80

        IF (CUR.LE.0 .OR. CUR.GT.N) RETURN

        PDG = K(CUR,2)

        IF (ISWEAKOPENCHARM(PDG)) THEN
          OPENCHARMMOTHERINDEX = CUR
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
C Weak open-charm hadrons
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
C Parent species bin: D0, D+, Ds, Lambda_c, other weak charm
C=======================================================================

      INTEGER FUNCTION GETPARENTBIN(PDG)

      IMPLICIT NONE

      INTEGER PDG, APDG

      APDG = ABS(PDG)

      GETPARENTBIN = 5
      IF (APDG.EQ.421)  GETPARENTBIN = 1
      IF (APDG.EQ.411)  GETPARENTBIN = 2
      IF (APDG.EQ.431)  GETPARENTBIN = 3
      IF (APDG.EQ.4122) GETPARENTBIN = 4

      RETURN
      END


C=======================================================================
C Hard-subprocess grouping from PYTHIA6 ISUB.
C This is NOT yet full HF source tagging.
C It groups the hard process that PYTHIA reports through MSTI(1).
C=======================================================================

      INTEGER FUNCTION GETPROCESSBIN(ISUB)

      IMPLICIT NONE

      INTEGER ISUB

C     1 = direct qqbar -> ccbar, massive HF pair creation
C     2 = direct gg    -> ccbar, massive HF pair creation
C     3 = qg -> qg
C     4 = qq / qqbar scattering
C     5 = gg -> qqbar
C     6 = gg -> gg
C     7 = semihard/MPI
C     8 = other

      GETPROCESSBIN = 8

      IF (ISUB.EQ.81) GETPROCESSBIN = 1
      IF (ISUB.EQ.82) GETPROCESSBIN = 2

      IF (ISUB.EQ.28) GETPROCESSBIN = 3

      IF (ISUB.EQ.11) GETPROCESSBIN = 4
      IF (ISUB.EQ.12) GETPROCESSBIN = 4
      IF (ISUB.EQ.13) GETPROCESSBIN = 4

      IF (ISUB.EQ.53) GETPROCESSBIN = 5

      IF (ISUB.EQ.68) GETPROCESSBIN = 6

      IF (ISUB.EQ.96) GETPROCESSBIN = 7

      RETURN
      END

      
C=======================================================================
C Force particle and antiparticle to electron decay channels
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
C Return inclusive electron BR from PYTHIA BRAT table
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
C Kinematics
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


      DOUBLE PRECISION FUNCTION GETMASS(I)

      IMPLICIT NONE

      INTEGER I
      INTEGER N, NPAD, K
      DOUBLE PRECISION P, V
      COMMON /PYJETS/ N, NPAD, K(4000,5), P(4000,5), V(4000,5)

      GETMASS = P(I,5)

      RETURN
      END
