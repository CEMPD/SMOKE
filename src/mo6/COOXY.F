      REAL FUNCTION COOXY(MY,IV,INH) 
C 
C  COOXY computes the oxy fuel effect of industry average fuel and 
C  reformulated gas adjustments.  Due to the introduction of 
C  oxidation and 3-way catalyst on LDTs, special handling of those 
C  technologies was needed. 
C 
C  Called by FUEL. 
C 
C  Calls ITAMPT, P3WPOX, and NMOGFID 
C 
C  Changes: (Last change first) 
C 
C  04 Dec 00 @EPA-bag Task X6, removed MY, JDX, and IDX from MYCODE.I; removed 
C            MYCODE.I; model year now passed through parameter lists 
C  15 May 00 @EPA-ddj BUG27, New function. 
C 
C  Input on call: 
C 
C    argument list: MY,IV,INH 
C    common blocks: 
C     /FLAGS1/ DTGFLG, OXYFLG 
C     /FLAGS2/ NEWFLG 
C     /OXY1/   SHRMKT,OXYCNT 
C     /RFORM1/ IPHASE,SEAFLG,RFGON 
C     /VVNAMS/ VVDSL, VVHDGV, VVLDGT12, VVLDGT34  (IVTYPE.I) 
C     /VVTYPS/ VTLDGV, VTMC  (IVTYPE.I) 
C 
C  Output on return: 
C 
C    common blocks: 
C              F_CF 
C 
C  Local array subscripts: 
C 
C     COLDGV(10)  - COLDGV(EMY) 
C     COLDGT(10)  - COLDGT(EMY) 
C     COHDGV(10)  - COHDGV(EMY) 
C     COMC(10)    - COMC(EMY) 
C     EMY depends on the emitter type (Normals and High) and 
C          model years  
C 
C 
C  Local variable dictionary: 
C 
C   Name   Type              Description 
C  ------  ----  ------------------------------------------------------- 
C  BENFAC   R    benefit factor for determining oxy adjustment 
C  CAT      R    fraction of vehicles that are catalyst equipped 
C  COLDGV   R    Hardcoded oxygen fuel correction factors for LDGV 
C  COLDGT   R    Hardcoded oxygen fuel correction factors for LDGT 
C  COHDGV   R    Hardcoded oxygen fuel correction factors for HDGV 
C  COMC     R    Hardcoded oxygen fuel correction factors for MC 
C  NOCAT    R    fraction of vehicles that are not catalyst equipped 
C  OXYCT    R    local variable for the oxygen content 
C  P3W      R    percent of vehicles with 3-way catalyst 
C  POX      R    percent of vehicles with ox catalyst 
C 
C 
C  Notes: 
C 
C 
      IMPLICIT NONE 
      INCLUDE 'BASEQ9.I' 
      INCLUDE 'FLAGS1.I' 
      INCLUDE 'FLAGS2.I' 
      INCLUDE 'IVTYPE.I' 
      INCLUDE 'OXY1.I' 
      INCLUDE 'RFORM1.I' 
C 
      INTEGER, INTENT(IN)      :: INH 
      INTEGER, INTENT(IN)      :: IV 
      INTEGER, INTENT(IN)      :: MY 
C 
      REAL NOCAT 
      REAL OXYCT 
      REAL CAT 
      REAL P3W 
      REAL POX 
      REAL BENFAC 
      REAL, PARAMETER :: COMC = -6.6 
C 
      INTEGER IFUEL 
      INTEGER IVFLAG 
C 
C    COLDGV(EMY) 
C       EMY =1 pre-1981 Non catalyst high emitter 
C       EMY =2 pre-1981 Catalyst high emitter 
C       EMY =3 post 1981 high emitter 
C       EMY =4 pre-1981 Non catalyst normal emitter 
C       EMY =5 pre-1981 Catalyst normal emitter 
C       EMY =6 1981-1985 
C       EMY =7 1986-1987 No ADL 
C       EMY =8 1986-1987 ADL 
C       EMY =9 1988-1993 
C       EMY =10 1994+ (not used) 
C 
      REAL, SAVE, DIMENSION(10)   :: COLDGV = 
     &  (/-6.6, -9.4, -5.3, -6.6, -9.4, -4.0, -5.7, -4.8, -3.1, 0.0/) 
C 
C    COLDGT(EMY) 
C       EMY =1 pre-1984 Non catalyst high emitter 
C       EMY =2 pre-1984 Catalyst high emitter 
C       EMY =3 post 1984 high emitter 
C       EMY =4 pre-1984 Non catalyst normal emitter 
C       EMY =5 pre-1984 Catalyst normal emitter 
C       EMY =6 1984-1985 
C       EMY =7 1986-1987 No ADL 
C       EMY =8 1986-1987 ADL 
C       EMY =9 1988-1993 
C       EMY =10 1994+ (not used) 
C 
      REAL, SAVE, DIMENSION(10)   :: COLDGT = 
     &  (/-6.6, -9.4, -5.3, -6.6, -9.4, -4.0, -5.7, -4.8, -3.1, 0.0/) 
C 
C    COHDGV(EMY) 
C       EMY =1 Non catalyst 
C       EMY =2 Oxidation catalyst 
C       EMY =3 three-way catalyst 
C 
      REAL, SAVE, DIMENSION(3)    :: COHDGV = (/-6.6, -9.4, -3.1/) 
C 
      DIMENSION OXYCT(2) 
C 
C  Initialize 
C 
      BENFAC = 0.0 
      COOXY  = 1.0 
C 
C  If vehicle is diesel fueled, no adjustment is necessary. 
C 
      IF(VVDSL(IV) .NE. 1) THEN 
C 
C  Return if no reformulated or oxygenated fuel effects are to be calculated. 
C 
        IF(RFGON .OR. OXYFLG.NE.1) THEN 
C 
C  Local variable for oxygen content. 
C 
          DO IFUEL=1,2 
            OXYCT(IFUEL) = OXYCNT(IFUEL) 
          END DO 
C 
C  Sep-22-1994 @ CSC-bsg Request 426  Modularized the calculation of the 
C  3-Way oxidation values, by calling a function to calculate them. 
C  Compute the number of catalyst and non-catalyst vehicles, 
C 
          CALL P3WPOX(MY,IV,P3W,POX) 
          CAT = P3W + POX 
          NOCAT = 1.0 - CAT 
C 
C Here we determine the vehicle type and we assign the default IVFLAG 
C 
          IF (IV.EQ.1) THEN 
             IVFLAG = 1 
          ELSEIF (VVLDGAS(IV).EQ.1) THEN 
             IVFLAG = 2 
          ELSEIF(VVHDGV(IV).EQ.1) THEN 
             IVFLAG = 3 
          ELSEIF (IV.EQ.VTMC) THEN 
             IVFLAG = 4 
          ELSE 
             IVFLAG = 5 
          END IF 
C 
          SELECT CASE(IVFLAG) 
C 
C  Here we compute the LDGV 
C 
            CASE(1) 
              SELECT CASE(INH) 
C 
C    LDGV oxy effect for HIGH emitters. 
C 
                CASE(2) 
                  BENFAC = COLDGV(2)*POX + COLDGV(1)*NOCAT + 
     &                     P3W*COLDGV(3) 
                CASE DEFAULT 
C 
C      LDGV oxy effect for NORMAL emitters. 
C      NOTE: Normals emitters for post-1995 vehicles do not receive 
C            an OXY benefit. 
C 
                IF(MY.LT.1996) THEN 
                  BENFAC = COLDGV(5)*POX + COLDGV(4)*NOCAT 
                  IF(MY.LT.1981) THEN 
                      BENFAC = BENFAC + P3W*COLDGV(5) 
                  ELSE IF(MY.LE.1985) THEN 
                      BENFAC = BENFAC + P3W*COLDGV(6) 
                  ELSE IF(MY.LE.1987) THEN 
                      BENFAC = BENFAC + P3W*(COLDGV(7)+COLDGV(8))/2 
                  ELSE 
                      BENFAC = BENFAC + P3W*COLDGV(9) 
                  END IF 
C 
                  IF(MY.GE.1994) THEN 
C 
C                    Only Tier0 vehicles affected.  Normal emitting 
C                    Tier1 and LEV vehicles are not affected by different 
C                    oxygen levels. 
C 
                     BENFAC = BENFAC*LEVIMP(1,MY-1993,IV) 
                  END IF 
                END IF 
              END SELECT 
C 
C Here We compute the LDGT 
C 
            CASE(2) 
              SELECT CASE(INH) 
                CASE(2) 
C 
C    LDGT oxy effect for HIGH emitters. 
C 
                  BENFAC = COLDGT(2)*POX + COLDGT(1)*NOCAT + 
     &                     P3W*COLDGT(3) 
C 
                CASE DEFAULT 
C 
C      LDGT oxy effect for NORMAL emitters. 
C 
                  IF(MY.LT.1996) THEN 
                    BENFAC = COLDGT(5)*POX + COLDGT(4)*NOCAT 
                    IF(MY.LT.1981) THEN 
                      BENFAC = BENFAC + P3W*COLDGT(5) 
                    ELSEIF(MY.LE.1985) THEN 
                      BENFAC = BENFAC + P3W*COLDGT(6) 
                    ELSEIF(MY.LE.1987) THEN 
                      BENFAC = BENFAC + P3W*(COLDGT(7)+COLDGT(8))/2 
                    ELSE 
                      BENFAC = BENFAC + P3W*COLDGT(9) 
                    END IF 
C 
                    IF(MY.GE.1994) THEN 
C 
C                  Only Tier0 vehicles affected.  Normal emitting 
C                  Tier1 and LEV vehicles are not affected by different 
C                  oxygen levels. 
C 
                      BENFAC = BENFAC*LEVIMP(1,MY-1993,IV) 
                    END IF 
                  END IF 
              END SELECT 
C 
C Here We compute the HDGV 
C 
              CASE(3) 
                BENFAC = COHDGV(3)*P3W + COHDGV(2)*POX + COHDGV(1)*NOCAT 
C 
C Here We calculate the VTMC 
C 
              CASE (4) 
                BENFAC = COMC 
C 
C Here the case for no oxygen effects 
C 
              CASE DEFAULT 
                CONTINUE 
            END SELECT 
C 
            COOXY = SHRMKT(1) 
     *            + SHRMKT(2)*(1.0+OXYCT(1)*BENFAC) 
     *            + SHRMKT(3)*(1.0+OXYCT(2)*BENFAC) 
          END IF 
      END IF 
C 
      RETURN 
      END 
 
