      REAL FUNCTION ACCF(IV,IP,INH,MY,IR,ISTD,EF,SPD) 
C 
C AC pollution emission based on compressor loading for 
C light duty cars and trucks, each of the three pollutants, 
C normal or high emitters, heat index(temperature & humidity), 
C and A/C unit repair rates(see reports M6.ACE.001 and 
C M6.ACE.002 for details). 
C 
C Changes : 
C 
C  20 Feb 01 @EPA-bag Removed speed parameters, i.e., SPEED9.I 
C  07 Aug 00 @EPA-djb Bug219, Added model year condition to SFTP calculations. 
C  20 Jul 00 @EPA-ddj Bug194 Removed IDX from SFTP_AC parameter list. 
C  16 MAY 00 @EPA-bag    deleted diagnostic print statement 
C  21 Jan 00 @DynTel-ddj 1-020 Added SFTP A/C correction. 
C  14 Jan 00 @DynTel-km  1-013 Cleaned the code 
C  13 Jan 00 @DynTel-ddj 1-013 Debugged non-working version 
C            provided by EPA. 
C  14 Dec 99 @EPA  New subroutine. 
C 
C Called by: SPEED_CF 
C 
C Calls: OFFCYCLE, SFTP_AC, HEAT, QUITER. 
C 
C Input on call: 
C 
C   argument list : IV, IP, INH, MY, IR, ISTD, EF, SPD 
C 
C   common blocks : 
C    /AC/ CLD_CVR, SUNRISE, SUNSET, PEAKSUN(2), 
C   &     VHCLES_AC_EQPD(2,27), FULL_A(10), 
C   &     FULL_B(10), DMND(3,5), NTFCNNG_AC_UNITS(2,20)   
C    /VTNAMS/ VVLDGAS (IVTYPE.I) 
C    /TEMPS/ TEMHRLY(24) 
C 
C Output on return: 
C 
C   argument list : ACCF 
C 
C   common blocks : 
C 
C   function      : ACCF 
C 
C Local variable / array dictionary : 
C 
C Name            Type         Description 
C -----          ------    ----------------------------------------- 
C AC_FCTR          R        A/C full use factor based on compressor loading 
C IP               I        Pollutant ( IP = 1(HC), 2(CO), 3(NOx) ) 
C SPD              R        Speed (mph) 
C VHCLS_FCNNG_AC   R        Fraction of vehicles with functioning AC 
C INERR            I        Input error counter 
C 
C Notes : 
C 
C 
      IMPLICIT NONE 
C 
      INCLUDE 'AC.I' 
      INCLUDE 'FLAGS5.I' 
      INCLUDE 'IVTYPE.I' 
CCC      INCLUDE 'SPEED9.I' 
      INCLUDE 'SFTP.I' 
C 
      REAL, INTENT(IN)       :: EF 
      REAL, INTENT(IN)       :: SPD 
C 
      INTEGER, INTENT(IN)    :: INH 
      INTEGER, INTENT(IN)    :: IP 
      INTEGER, INTENT(IN)    :: IR 
      INTEGER, INTENT(IN)    :: ISTD 
      INTEGER, INTENT(IN)    :: IV 
      INTEGER, INTENT(IN)    :: MY 
C 
      REAL AC_FCTR 
      REAL AC_SFTP 
      REAL FRAC 
      REAL FRACT 
      REAL SFTP_AC 
C 
      INTEGER IMY 
      INTEGER IVI 
C 
C     Initialize function 
C 
      AC_FCTR = 0.0 
C 
C Next  : A/C pollution emission based on compressor loading  
C         for each vehicle and pollutant. Note that for NOx there  
C         is only one emitter type. Also for HC in the low speed 
C         range (or emitter type) the LDV curve is set equal to 
C         the LDT curve for speeds less than or equal to 8 MPH. 
C 
C         Full usage emissions are related to pollutant type and  
C         vehicle class/emitter type. NOx parameters have only a  
C         dependence on vehicle class and only have values for 
C         LDV and LDT. Normal hydrocarbon emitting LDV and LDT  
C         both have individual ac full usage parameters. All  
C         high hydrocarbon vehicle classes use the same set of  
C         parameters. The parameters for carbon monoxide vehicle 
C         class and emitter types follow the same rules as the  
C         hydrocarbons. 
C 
C         IR is the roadway type index 
C            1:  Freeway / Interstate / Expressway 
C            2:  Arterial / Collector 
C            3:  Local / Urban 
C            4:  Ramp 
C 
      IF(EF .GT. 0.0) THEN 
        SELECT CASE (IP) 
C 
C ********** NMHC ********** 
C 
         CASE(1) 
           IF ( INH .EQ. 1 ) THEN 
C 
C  NMHC A/C correction factor for Normal emitters  : 
C 
             IF(IR.EQ.1 .OR. IR.EQ.2 .OR. IR.EQ.4) THEN 
                AC_FCTR = FULL_A(1) * SPD 
             ELSEIF(IR.EQ.3) THEN 
                AC_FCTR = FULL_A(2) * EF 
             ENDIF 
C 
C  NMHC A/C correction factor for High emitters : 
C 
           ENDIF 
C 
C ********** CO  ********** 
C 
         CASE(2) 
C 
           IF(IR.EQ.1 .OR. IR.EQ.2 .OR. IR.EQ.4) THEN  ! Not Local 
             IF(INH.EQ.1) THEN                         ! Normal Emitter 
               IF(IV.EQ.1) THEN                        ! LDGV Only 
                 AC_FCTR = (FULL_A(3) * EF) + (FULL_B(3) * SPD) 
               ELSEIF( VVLDGAS(IV).EQ.1 ) THEN 
                 AC_FCTR = FULL_A(4) * SPD 
               ENDIF                                   ! End Normal 
             ELSEIF(INH.EQ.2 .AND. SPD .LT. 19.0) THEN ! If < 19 mph 
               AC_FCTR = FULL_A(5) * EF 
             ENDIF 
C 
           ELSEIF(IR.EQ.3) THEN                        ! Local 
             IF(INH.EQ.1) THEN                         ! Normal Emitter 
               AC_FCTR = FULL_A(6) * EF 
             ELSE                                      ! High Emitter 
               AC_FCTR = FULL_A(7) * EF 
             ENDIF 
C 
           ENDIF 
C 
C ********** NOx  ********** 
C 
         CASE(3) 
C 
           IF(SPD .LE. 0.0 .OR. EF .LE. 0.0) THEN 
             AC_FCTR = 0.0 
           ELSE 
             IF(IR.EQ.1 .OR. IR.EQ.2 .OR. IR.EQ.3) THEN  ! Not Ramp 
               IF( IV.EQ.1 ) THEN                        ! LDGV Only 
                  AC_FCTR = LOG10(EF + 1.0) * ( FULL_A(8) - 
     *                    (FULL_B(8) * LOG10(SPD)) ) 
               ELSEIF( VVLDGAS(IV).EQ.1 ) THEN           ! LDGTs 
                AC_FCTR = LOG10(EF + 1.0) * ( FULL_A(9) - 
     *                    (FULL_B(9) * LOG10(SPD)) ) 
               ENDIF 
             ELSE IF(IR .EQ. 4) THEN 
                AC_FCTR = FULL_A(10)*LOG10(EF+1.0)         ! Ramps Only 
             ENDIF 
C 
           ENDIF 
C 
        END SELECT 
C 
      END IF               ! EF Greater than 0.0 
C 
C   Apply SFTP_AC effect to only 2001 and later model year 
C   light duty gasoline cars and trucks. 
C 
      IF(MY .GE. 2001) THEN 
        IF(VVLDGAS(IV) .EQ. 1 .AND. SFTP_AC_FLAG .EQ. 1) THEN 
          IF(IP .EQ. 2) THEN 
            AC_SFTP = EF*SFTP_AC(INH,MY,IV,IP,ISTD,SPD,FRAC) - EF 
C 
C     For CO, the adjusted emission values shall not be larger than the 
C     non-SFTP emissions. 
C 
            AC_SFTP = MIN(AC_FCTR, AC_SFTP) 
            IF(MY.GE.2001 .AND. MY.LE.2003) THEN 
               IMY = MY - 2000 
               IF(IMY .GE. 4) THEN      ! If Model Year > 2004 
                  FRACT = 1.0           ! Full Phase-In. 
               ELSE 
                  IF(VVLDGT34(IV) .EQ. 0) THEN 
                     IVI = 1 
                  ELSE               ! LDGT3 or LDGT4 
                     IVI = 2 
                  END IF 
                  FRACT = SFTP_FRAC(IVI,IMY) 
               END IF 
               AC_FCTR = FRACT*AC_SFTP + (1 - FRACT)*AC_FCTR 
            ELSE 
               AC_FCTR = AC_SFTP 
            END IF 
          ELSE       ! IP = HC or NOx 
            AC_FCTR = AC_FCTR*SFTP_AC(INH,MY,IV,IP,ISTD,SPD,FRAC) 
          END IF 
        END IF 
      END IF       ! Model Year Condition 
C 
      ACCF = MAX(AC_FCTR, 0.0) 
      RETURN 
      END 
