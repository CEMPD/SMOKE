      SUBROUTINE PMFRACA(MY,IV,INERR)
C
C  This subroutine calculates the array FRACA in PART5
C  using the technology group fractions found in MOBILE6.
C  FRACA are the fractions of vehicles with the following
C  technology mixes: Ox No Air, Ox Air, 3w No Air, 3w Air
C
C  Called by PMGSO4SO2, PMLEAD
C
C  Calls None
C
C  Changes: (Last change first)
C
C  15 Apr 02: AIR Task Bug Fix: changed call to QUITER so that variable
C             types match.
C  13 Feb 01: AIR Task 11: New function.
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
C    /PART2/ FRACA
C
C  FRACA(I)  ---  1: Ox No Air
C                 2: Ox Air
C                 3: 3w No Air
C                 4: 3w Air
C
C  Local variable dictionary:
C
C   Name   Type              Description
C  ------  ----  -------------------------------------------------------
C  IG1      I    Model year pointer into TGS array, via function ITAMPT
C  SUM23    R    Sum of T(2) and T(3)
C  SUM56    R    Sum of T(5) and T(6)
C  T        R    temporary TGS array, 2=air pump/catalyst, 
C                3=no air pump/catalyst, 5=EGR/catalyst,
C                6=no EGR/catalyst, 1 and 4 not used (yet)
C
      IMPLICIT NONE
C      
      INCLUDE 'EGSCAL.I'
      INCLUDE 'IVTYPE.I'
      INCLUDE 'LOOKUP.I'
      INCLUDE 'PART2.I'
C
C  Declare external function.
C
      INTEGER, EXTERNAL :: ITAMPT
C
C  Declare local variables.
C    
      INTEGER IG1
      INTEGER INERR
      INTEGER IV
      INTEGER MY
      REAL    SUM23
      REAL    SUM56
      REAL    T(6)
C
C  Initialize array FRACA
C
      FRACA=0.0
C
C  If motorcycles, then no catalysts.
C
      IF(IV.EQ.VTMC) RETURN
C
C  Determine the vehicle class pointer into TGS array.
C
      IVTAM=IV
      IF (IV.EQ.VTGASBUS) IVTAM = 14
C
C  Set up temporary TGS variables.
C  Programming note: TGS must use IQG and IVTAM via LOOKUP.I.
C
      SUM23=0.
      DO IQG=2,3
        IG1=ITAMPT(MY,1)
        T(IQG)=TGS(IG1,IQG,IVTAM)
        SUM23=SUM23+T(IQG)
      END DO
C
      SUM56=0.
      DO IQG=5,6
        IG1=ITAMPT(MY,1)
        T(IQG)=TGS(IG1,IQG,IVTAM)
        SUM56=SUM56+T(IQG)
      END DO
C
C  Call QUITER when negative technology fractions occur.
C
      IF(SUM56.GT.SUM23) CALL QUITER(FLOAT(MY),IV,170,INERR)
C
C  Compute the equivalent FRACA, if possible
C
      IF(SUM23.GT.0.) THEN
        FRACA(1)=T(3)*((SUM23-SUM56)/SUM23)  !Ox No Air
        FRACA(2)=T(2)*((SUM23-SUM56)/SUM23)  !Ox Air
        FRACA(3)=T(3)*SUM56/SUM23            !3w No Air
        FRACA(4)=T(2)*SUM56/SUM23            !3w Air
      ENDIF
C
      RETURN
      END
