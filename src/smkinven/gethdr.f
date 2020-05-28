
        SUBROUTINE GETHDR( MXDATA, CFLAG, YFLAG, DFLAG, 
     &                     LINE, ICC, INY, NPOA, EOS )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine interprets the header lines from the inventory files
C
C  PRECONDITIONS REQUIRED:
C      Line must be read in and left-justified. Valid pollutant and activity
C      list must be populated. Country codes list must be populated. The 
C      values of ICC, INY, and NPOA must be initialized.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      created by M. Houyoux (2/2000) 
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C...........   MODULES for public variables
C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: MXIDAT, INVDUNT, INVDCNV, INVDNAM,
     &                      ITCASA, ITNAMA, NINVTBL, ITKEEPA

C.........  This module contains the arrays for state and county summaries
        USE MODSTCY, ONLY: NCOUNTRY, CTRYNAM, CTRYCOD

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: DATPOS, TMPNAM, NPOLID

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)           CRLF
        INTEGER                GETNLIST
        INTEGER                INDEX1
        INTEGER                STR2INT
        REAL                   UNITFAC 
        LOGICAL                USEEXPGEO

        EXTERNAL    CRLF, GETNLIST, INDEX1, STR2INT, UNITFAC, USEEXPGEO

C...........   SUBROUTINE ARGUMENTS
C...........   NOTE that NDROP and EDROP are not used at present
        INTEGER     , INTENT (IN) :: MXDATA  ! max no data variables allowed
        LOGICAL     , INTENT (IN) :: CFLAG   ! true: country header expected
        LOGICAL     , INTENT (IN) :: YFLAG   ! true: year header expected
        LOGICAL     , INTENT (IN) :: DFLAG   ! true: data names header expected
        CHARACTER(*), INTENT(IN OUT):: LINE  ! full record in string
        INTEGER     , INTENT(OUT) :: ICC     ! country code
        INTEGER     , INTENT(OUT) :: INY     ! inventory year
        INTEGER     , INTENT(OUT) :: NPOA    ! no. data names
        INTEGER     , INTENT(OUT) :: EOS     ! error status

C...........   Local allocatable arrays
        CHARACTER(IOULEN3), ALLOCATABLE :: UNITS( : )

C...........   Other local variables
        INTEGER         I, J, K, L, L1, L2, V  ! counters and indices

        INTEGER         COD       !  tmp data variable position in valid list
        INTEGER         IOS       !  i/o status
        INTEGER         NFINAL    !  final count of all "kept" pollutants
        INTEGER      :: NUNIT = 0 !  i/o status

        LOGICAL, SAVE:: ACT_FLAG  = .FALSE. ! true: activities in data file
        LOGICAL, SAVE:: FIRSTIME  = .TRUE.  ! first time subroutine is called
        LOGICAL, SAVE:: UNITFLAG  = .FALSE. ! true: units have already been set

        REAL            DFAC    !  denominator units conversion factor
        REAL            NFAC    !  numerator units conversion factor

        CHARACTER       CBUF    !  single char buffer
        CHARACTER(20)   CNTRY   !  country name
        CHARACTER(300)  MESG    !  message buffer
        CHARACTER(600)  BUFFER  !  tmp line buffer, upper case

        CHARACTER(IOVLEN3) CVAR      ! tmp variable name
        CHARACTER(IOULEN3) DBUF      ! tmp units denominator
        CHARACTER(IOULEN3) NBUF      ! tmp units numerator
        CHARACTER(IOULEN3) UBUF      ! tmp units

        CHARACTER(16) :: PROGNAME = 'GETHDR' ! Program name

C***********************************************************************
C   begin body of subroutine GETHDR

C.........  First time the routine is called, initialize values
        IF( FIRSTIME ) THEN

            IF( CFLAG ) ICC  = -1
            IF( YFLAG ) INY  = 0
            IF( DFLAG ) NPOA = 0
            FIRSTIME = .FALSE.

        END IF

C.........  Initialize i/o status
        EOS = 0

C.............  Scan for header lines
        IF( LINE( 1:1 ) .EQ. CINVHDR ) THEN
            L1 = INDEX( LINE, '=' ) 
            L2 = LEN_TRIM( LINE )

C.............  Convert to upper case
            BUFFER = LINE
            CALL UPCASE( BUFFER )

C.............  Check for inventory type
            IF ( BUFFER(2:5) .EQ. 'TYPE' ) THEN

C.................  Try to find activity in name of data, otherwise, assume
C                   emissions
                I = INDEX( BUFFER, 'ACTIVITY' ) 
                IF( I .GT. 0 ) THEN
                    ACT_FLAG = .TRUE.
                END IF

C.............  Check for country name
            ELSE IF ( BUFFER(2:8) .EQ. 'COUNTRY' ) THEN 

C.................  Skip country header if using expanded geographic codes
                IF ( USEEXPGEO() ) RETURN

                IF( L1 < 1 ) L1 = 8 

                CNTRY = ADJUSTL( LINE( L1+1:L2 ) )
                I   = INDEX1( CNTRY, NCOUNTRY, CTRYNAM )             

                IF( NCOUNTRY .LE. 0 ) THEN
                    EOS = 2
                    MESG = 'INTERNAL ERROR: Valid country list is ' //
     &                     'uninitialized.'
                    CALL M3MSG2( MESG )

                ELSE IF ( I .LE. 0 ) THEN
                    EOS = 1
                    L = LEN_TRIM( CNTRY )
                    MESG = 'ERROR: Country name "' // CNTRY( 1:L ) //
     &                     '" from inventory is not in country names '//
     &                     'file.'
                    CALL M3MSG2( MESG )

                END IF

                ICC   = STR2INT( CTRYCOD( I ) ) / 100000

C.............  Check for inventory year
            ELSE IF ( BUFFER(2:5) .EQ. 'YEAR' ) THEN 

                IF( L1 < 1 ) L1 = 5

                INY = STR2INT( BUFFER( L1+1:L2 ) )

                IF ( INY .LT. 1971 .OR. INY .GT. MXIYEAR ) THEN
                    WRITE( MESG, 94010 ) 'WARNING: Inventory year', INY,
     &                     'is outside of expected range.'
                    CALL M3MSG2( MESG )                    
                END IF

C.............  Check for units field to be used and compare to valid units
C               from master input list
C.............  Compute conversion factors for adjusting emissions in reader
C               routines
            ELSE IF ( BUFFER(2:6) .EQ. 'UNITS' ) THEN ! read in units                

C.................  Get the number of units fields and allocate memory
                BUFFER = BUFFER( 7:L2 )
                L = LEN_TRIM( BUFFER )
                NUNIT = GETNLIST( L, BUFFER )

                IF( UNITFLAG ) THEN
                    MESG = 'ERROR: UNITS header encountered again in '//
     &                     'input. This is not allowed.'
                    CALL M3MSG2( MESG )
                    EOS = 5
                    RETURN

                ELSE IF( .NOT. ALLOCATED( DATPOS ) ) THEN
                    MESG = 'ERROR: UNITS header must be after DATA ' //
     &                     'or POLID field.'
                    CALL M3MSG2( MESG )
                    EOS = 5
                    RETURN

                ELSE IF( .NOT. ACT_FLAG ) THEN
                    MESG = 'ERROR: UNITS header only valid for ' //
     &                     'activity data.'
                    CALL M3MSG2( MESG )
                    EOS = 5
                    RETURN

                ELSE IF( NUNIT .NE. NPOA ) THEN

                    WRITE( MESG,94010 ) 'ERROR: number of UNITS ' //
     &                     'header entries (', NUNIT, ') is not ' //
     &                     CRLF() // BLANK10 // 'consistent with ' //
     &                     'number of DATA or POLID entries (', NPOA,
     &                     ')'
                    CALL M3MSG2( MESG )
                    EOS = 5
                    RETURN

                END IF

C.................  Allocate memory for local units field
                ALLOCATE( UNITS( NUNIT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'UNITS', PROGNAME )
  
C.................  Parse the header line into the unit names
                CALL PARSLINE( LINE( 7:L2 ), NUNIT, UNITS )

C.................  Post-process units and store in final arrays
                K = 0
                DO V = 1, NPOA

C.....................  Get position in master pollutants/activities list
                    J  = DATPOS( V )

C.....................  Convert input units to valid SMOKE output units
                    CALL UNITMATCH( UNITS( V ) )

C.....................  Set conversion factors for numerator and denominator,
C                       if any
                    NFAC = UNITFAC( UNITS( V ), INVDUNT( J ), .TRUE.  )
                    DFAC = UNITFAC( UNITS( V ), INVDUNT( J ), .FALSE. )

                    IF( NFAC .LT. 0. ) NFAC = 1.
                    IF( DFAC .LT. 0. ) DFAC = 1.

C.....................  Error if units are not the same and the conversion
C                       factor is equal to one
                    IF( UNITS( V ) .NE. INVDUNT( J ) .AND.
     &                  NFAC .EQ. 1. .AND. DFAC .EQ. 1.      ) THEN

                        MESG = 'ERROR: Units conversion could not be '//
     &                         'made - use different input units.'
                        CALL M3MSG2( MESG )
                        EOS = 5

                        RETURN

                    END IF

C.....................  Store current values of conversion factor
                    INVDCNV( J ) = NFAC * DFAC

                END DO           ! No. variables

                DEALLOCATE( UNITS )

                UNITFLAG = .TRUE.

            ELSE IF ( BUFFER(2:6) .EQ. 'POLID' .OR.
     &                BUFFER(2:6) .EQ. 'DATA '       ) THEN ! read in data names

C..................... Deallocate names for pollutant, if must
                IF( ALLOCATED( TMPNAM ) ) DEALLOCATE( TMPNAM,DATPOS )

C.................  Allocate memory for current file for reading pol names
C                   and storing positions in master list
                BUFFER = BUFFER( 7:L2 )
                L = LEN_TRIM( BUFFER )
                NPOA = GETNLIST( L, BUFFER )

                ALLOCATE( TMPNAM( NPOA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'TMPNAM', PROGNAME )
                ALLOCATE( DATPOS( NPOA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DATPOS', PROGNAME )
                TMPNAM = ' '
                DATPOS = 0

C.................  Set special value of error status when the number of
C                   data variables attempted is greater than max allowed
                IF( NPOA .GT. MXDATA .AND. DFLAG ) THEN
                    EOS = 4
                END IF

C.................  Parse the header line into the pollutant names
                CALL PARSLINE( BUFFER, NPOA, TMPNAM )

C.................  Store the position in master list of each pollutant
C.................  Write error if pollutant is not found.
                NPOLID = NPOA    ! Store original no of pollutants
                NFINAL = 0
                DO V = 1, NPOA

                    CVAR = TMPNAM( V )
                    L = LEN_TRIM( CVAR )

C.....................  Look for variable in SMOKE variable names
                    COD = INDEX1( CVAR, MXIDAT, INVDNAM )
                    IF( COD .LE. 0 ) THEN

C........................  Look for variable in Inventory Pollutant codes                       
                        COD = INDEX1( CVAR, NINVTBL, ITCASA )

C.......................   If pollutant is not "kept", then take it
C                          out of the count and the list
                        IF( COD .GT. 0 ) THEN
                            IF( .NOT. ITKEEPA( COD ) ) CYCLE
                            COD = INDEX1( ITNAMA(COD), MXIDAT, INVDNAM )
                            DATPOS( V ) = COD

                            NFINAL = NFINAL + 1
                            TMPNAM( NFINAL ) = CVAR

C.........................  If not found in list of names or codes, then error
                        ELSE
                            EOS = 1
                            MESG = 'ERROR: Data variable "' // 
     &                           CVAR( 1:L )// '" not in master '//
     &                           'data variable list!'
                            CALL M3MSG2( MESG )
                        END IF

C.....................  Variable found in SMOKE names
                    ELSE
                        DATPOS( V ) = COD
                        NFINAL = NFINAL + 1

                    END IF

                END DO
C.................  Reset NPOA with final count that drops unkept variables
                NPOA = NFINAL

            END IF

C.............  If this is a header line, return here without checking that
C               all of the valid header fields are set properly
            RETURN

        END IF

        EOS = -1

C.........  If the line is not a header line, make sure that all of the 
C           important header lines have been read in...

C.........  Check for country header
        IF( CFLAG .AND. ICC .LT. 0 ) THEN
            EOS = 1
            ICC = 11         ! to turn off error message
            MESG = 'ERROR: Country name was not set with ' //
     &             '#COUNTRY header before first data line.'
            CALL M3MSG2( MESG )
        END IF

C.........  Check for inventory year header
        IF( YFLAG .AND. INY .LE. 0 ) THEN
            EOS = 1
            INY = 1       ! to turn off error message
            MESG = 'ERROR: Inventory year was not set with ' //
     &             '#YEAR header before first data line.'
            CALL M3MSG2( MESG )
        END IF

C.........  Check for pollutant names header
C.........  This check is coordinated with initialized value
        IF( DFLAG .AND. NPOA .EQ. 0 ) THEN
            EOS = 1
            NPOA = -1       ! to turn off error message
            MESG = 'ERROR: Data variable names were not set with ' //
     &             '#POLID or #DATA header before first data line.'
            CALL M3MSG2( MESG )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GETHDR
