      SUBROUTINE ATGROUP(MY,IV,INERR)
C
C  This subroutine calculates the GROUP fractions in AIR TOXICS
C  using the technology group fractions found in MOBILE6.
C
C  Called by ATCCALX
C
C  Calls None
C
C  Changes: (Last change first)
C
C  15 Apr 02: AIR Task Bug Fix 388, Modification of constants via QUITER.
C  20 Sep 01: AIR Task 15: Comments
C  18 Jun 01: AIR Task 15: New function.
C
C  Input on call:
C
C  argument list: MY,IV,INERR
C    common blocks:
C    /EGSCAL/ TGS
C    /LOOKUP/ IVTAM,IQG
C
C  Output on call:
C
C  argument list: 
C    common blocks:
C    /ATOX3/ GROUP
C
C  Temporary TGS array:
C  T(I)  ---  0: No Air Pump / No Catalyst (Not Used)
C             1: Air Pump Only / No Catalyst (Not Used)
C             2: Air Pump / Catalyst (Ox + 3way)
C             3: No Air Pump / Catalyst (Ox + 3way)
C             4: EGR Only / No 3way Catalyst (Not Used)
C             5: EGR / 3way Catalyst
C             6: No EGR / 3way Catalyst
C             7: No EGR / No 3way Catalyst (Not Used)
C
C  Group Tech Fractions:
C  GROUP(I) ---   1: PFI , 3way   , No Air, EGR
C                 2: PFI , 3way   , No Air, No EGR
C                 3: TBI , 3way   , No Air, EGR
C                 4: PFI , 3way+Ox, Air   , EGR
C                 5: PFI , 3way   , Air   , EGR
C                 6: TBI , 3way   , Air   , EGR
C                 7: TBI , 3way+Ox, Air   , EGR
C                 8: TBI , 3way   , No Air, No EGR
C                 9: CARB, 3way+Ox, Air   , EGR
C                10: All High Emitters
CC
      IMPLICIT NONE
C      
      INCLUDE 'ATOX3.I'   ! GROUP
      INCLUDE 'EGSCAL.I'  ! TGS
      INCLUDE 'INJECT.I'  ! TBI,PFI
      INCLUDE 'IVPCOM.I'  ! IVPTRT
      INCLUDE 'IVTYPE.I'  ! VTGASBUS,VVLDGAS
      INCLUDE 'LOOKUP.I'  ! IVTAN,IQG
C
C  Declare external function.
C
      INTEGER, EXTERNAL :: IEVPTR
      INTEGER, EXTERNAL :: ITAMPT
C
C  Declare parameter list
C
      INTEGER INERR      ! Error
      INTEGER IV         ! Vehicle class
      INTEGER MY         ! Model year
C
C  Declare local variables.
C    
      INTEGER IG1        ! Tampering index
      INTEGER IG5        ! TBI/PFI Fuel delivery group index
      INTEGER ITG        ! Technology group index (see above)
      INTEGER IVTEMP     ! Temporary vehicle class
      REAL    FDSWT(3)   ! Fuel delivery system: 1=Carb, 2=TBI, 3=PFI
      REAL    SUM        ! Summation
      REAL    T(6)       ! Temporary technology group (see above)
C
C  Initialize group array
C
      GROUP = 0.
C
C  Exit if not LDGV/T
C
      IF(VVLDGAS(IV).EQ.0) RETURN
C
C  Determine the fuel delivery fractions
C  (1=Carb, 2=TBI, 3=PFI)
C
      IVTEMP = IVPTRT(IV)
      IF (IV.EQ.VTGASBUS) IVTEMP = 14
      IG5 = IEVPTR(MY,IVTEMP,2,INERR)
      IF(INERR.GT.0) GOTO 99
      FDSWT(2) = TBI(IG5,IVTEMP)            !TBI
      FDSWT(3) = PFI(IG5,IVTEMP)            !PFI
      FDSWT(1) = 1.00 - FDSWT(2) - FDSWT(3) !CARB
C
C  Determine the vehicle class pointer into TGS array.
C
      IVTAM=IV
      IF (IV.EQ.VTGASBUS) IVTAM = 14
C
C  Set up temporary TGS variables.
C  Programming note: TGS must use IQG and IVTAM via LOOKUP.I.
C
      T = 0.
C
      SUM=0.
C
      DO IQG=2,3
        IG1=ITAMPT(MY,1)
        T(IQG)=TGS(IG1,IQG,IVTAM)
        SUM=SUM+T(IQG)
      END DO
C
      DO IQG=5,6
        IG1=ITAMPT(MY,1)
        T(IQG)=TGS(IG1,IQG,IVTAM)
      END DO
C
C  Compute the GROUP fractions, if possible. Split GROUPS 4-7
C  into their components via the Complex Model fractions.
C
      IF(SUM.GT.0.) THEN
        GROUP( 1)=T(3)*T(5)*FDSWT(3)/SUM        !PFI , 3way   , No Air, EGR
        GROUP( 2)=T(3)*T(6)*FDSWT(3)/SUM        !PFI , 3way   , No Air, No EGR
        GROUP( 3)=T(3)*T(5)*FDSWT(2)/SUM        !TBI , 3way   , No Air, EGR
        GROUP( 4)=T(2)*T(5)*FDSWT(3)*0.7262/SUM !PFI , 3way+Ox, Air   , EGR
        GROUP( 5)=T(2)*T(5)*FDSWT(3)*0.2738/SUM !PFI , 3way   , Air   , EGR
        GROUP( 6)=T(2)*T(5)*FDSWT(2)*0.0917/SUM !TBI , 3way   , Air   , EGR
        GROUP( 7)=T(2)*T(5)*FDSWT(2)*0.9083/SUM !TBI , 3way+Ox, Air   , EGR
        GROUP( 8)=T(3)*T(6)*FDSWT(2)/SUM        !TBI , 3way   , No Air, No EGR
        GROUP( 9)=T(2)*T(5)*FDSWT(1)/SUM        !CARB, 3way+Ox, Air   , EGR
        GROUP(10)=1.0                           !All High Emitters
      ELSE
        CALL QUITER(FLOAT(MY),IV,576,INERR)
      ENDIF
C
C  Normalize the Normal Emitter Fractions to 100%
C
      SUM=0.
C
      DO ITG=1,9
        SUM=SUM+GROUP(ITG)
      END DO
C
      DO ITG=1,9
        GROUP(ITG)=GROUP(ITG)/SUM
      END DO
C
  99  RETURN
      END
