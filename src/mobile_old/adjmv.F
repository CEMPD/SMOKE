
        REAL FUNCTION ADJMV( NSIZ, KEY1, KEY2, TABLE1, TABLE2, FACS  )

C***********************************************************************
C  function body starts at line 116
C
C  DESCRIPTION:
C       Assigns a mobile adjustments factor based on the KEYs, which can
C       either be FIPS code (KEY1) and road class code (KEY2), or cell number
C       (KEY1).  The latter situation is assumed when KEY2=0.
C
C  PRECONDITIONS REQUIRED:
C       Keys set to valid values from inventory. KEY1 set to county code for
C       county/road-class-based or cell number for cell-based.  KEY2 set to 
C       zero for cell-based and to road class code for road-class based.
C       Input arrays populated with table read from RDMVADJ that is 
C       consistent with the KEY settings.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       I/O API functions and subroutines.
C
C  REVISION  HISTORY:
C       Written 4/98 by M Houyoux
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Programs Group
C MCNC--North Carolina Supercomputing Center
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C env_progs@mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'MBDIMS3.EXT'     ! mobile-specific information

C...........   FUNCTION ARGUMENTS and their descriptions:

        INTEGER     NSIZ           !  adjustments table size
        INTEGER     KEY1           !  FIPS code or cell number
        INTEGER     KEY2           !  road class code, or zero
        INTEGER     TABLE1( NSIZ ) !  sorted FIP codes or cell codes
        INTEGER     TABLE2( NSIZ ) !  sorted road class codes or nothing
        REAL        FACS  ( NSIZ ) !  adjustement factors

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2 CRLF
        INTEGER     FIND1
        INTEGER     FIND2

        EXTERNAL    CRLF, FIND1, FIND2

C...........   LOCAL PARAMETERS:

        CHARACTER*5  BLANK5
        CHARACTER*10 BLANK10
        PARAMETER  ( BLANK5  = ' '
     &             , BLANK10 = ' ' )

C...........   Grouped tables and indices for FIP/RoadClass-based factors...
C...........   FIPS/state-codes only
        INTEGER   NG1
        INTEGER   FIPG1( MXFIP )
        REAL      FACG1( MXADJ1 )
        SAVE      NG1, FIPG1, FACG1

C...........   Road class codes only
        INTEGER   NG2
        INTEGER   RDTG2( NRCLAS )
        REAL      FACG2( NRCLAs )
        SAVE      NG2, RDTG2, FACG2

C...........   Ultimate default
        INTEGER   NG3
        REAL      FACG3
        SAVE      NG3, FACG3

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER         I                ! Counters and indices

        INTEGER         FIP              ! Temporary FIPS code
        INTEGER         RDT              ! Temporary road class
        INTEGER         STA              ! Temporary state code

        REAL            FAC              ! Temporary adjustment factor

        LOGICAL         CFLAG            ! True: matching is cell-type
        LOGICAL         EFLAG            ! Error flag
        LOGICAL         FIRSTIME         ! True: first time routine called
        DATA            FIRSTIME / .TRUE. /
        SAVE            FIRSTIME

        CHARACTER*256   MESG             ! message buffer

C***********************************************************************
C   begin body of function ADJMV

        EFLAG = .FALSE.

C.........  Initialize the return value in case no match is found
        ADJMV = 1. 

C.........  Determine the type of matching wanted
        CFLAG = .FALSE.
        IF( KEY2 .EQ. 0 ) CFLAG = .TRUE.

C.........  For the first time the routine is called, build tables for
C           FIPS-only and road-class-only matches
        IF( FIRSTIME .AND. .NOT. CFLAG ) THEN

            NG1 = 0
            NG2 = 0
            NG3 = 0
            DO 55 I = 1, NSIZ

                FIP = TABLE1( I )
                RDT = TABLE2( I )
                FAC = FACS  ( I )

                IF( FIP .EQ. 0 ) THEN

                    IF( RDT .EQ. 0 ) THEN
                         NG3   = 1
                         FACG3 = FAC

                    ELSE
                        NG2 = NG2 + 1
                        IF( NG2 .LE. NRCLAS ) THEN
                            RDTG2( NG2 ) = RDT
                            FACG2( NG2 ) = FAC
                        END IF

                    END IF

                ELSEIF( RDT .EQ. 0 ) THEN

                    NG1 = NG1 + 1
                    IF( NG1 .LE. MXFIP ) THEN
                        FIPG1( NG1 ) = FIP
                        FACG1( NG1 ) = FAC
                    END IF

                END IF

55          CONTINUE

            IF( NG1 .GT. MXFIP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 
     &                 'Number of FIP-only entries in mobile ' //
     &                 'adjustments file =', NG1, 'but max (MXFIP) =', 
     &                 MXFIP
                CALL M3MSG2( MESG )
            END IF

            IF( NG2 .GT. NRCLAS ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 
     &                 'Number of road-class-only entries in mobile ' //
     &                 'adjustments file =', NG2, 'but max (NRCLAS) =', 
     &                 NRCLAS
                CALL M3MSG2( MESG )
            END IF

            IF ( EFLAG ) THEN
                MESG = 'Problem processing mobile adjustments file.'
                CALL M3EXIT( 'ADJMV', 0, 0, MESG, 2 )
            ENDIF

            FIRSTIME = .FALSE.

        END IF

C.........  Find the cell-based match
        IF( CFLAG ) THEN

            I = FIND1( KEY1, NSIZ, TABLE1 )

            IF( I .GT. 0 ) ADJMV = FACS( I )

C.........  Otherwise, find the FIPS code and road-class-based match
        ELSE

C.............  FIPS code/road-class complete match
            I = FIND2( KEY1, KEY2, NSIZ, TABLE1, TABLE2 ) 
            IF( I .GT. 0 ) THEN
                ADJMV = FACS( I )
                RETURN
            END IF

C.............  State code/road-class match
            STA = KEY1/1000 * 1000
            I = FIND2( STA, KEY2, NSIZ, TABLE1, TABLE2 ) 
            IF( I .GT. 0 ) THEN
                ADJMV = FACS( I )
                RETURN
            END IF

C.............  FIPS code only match
            I = FIND1( KEY1, NG1, FIPG1 ) 
            IF( I .GT. 0 ) THEN
                ADJMV = FACG1( I )
                RETURN
            END IF
            
C.............  State code only match
            I = FIND1( STA, NG1, FIPG1 ) 
            IF( I .GT. 0 ) THEN
                ADJMV = FACG1( I )
                RETURN
            END IF
            
C.............  Road class code only match
            I = FIND1( KEY2, NG2, RDTG2 ) 
            IF( I .GT. 0 ) THEN
                ADJMV = FACG2( I )
                RETURN
            END IF

C.............  Ultimate default match
            IF( NG3 .GT. 0 ) THEN
                ADJMV = FACG3
                RETURN
            END IF

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

        END

