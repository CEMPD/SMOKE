
        SUBROUTINE GETHDR( MXDATA, MXIPPA, CFLAG, YFLAG, DFLAG, 
     &                     INVDNAM, LINE, ICC, INY, NPOA, EOS )

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
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C***************************************************************************

C...........   MODULES for public variables
C.........  This module contains the arrays for state and county summaries
        USE MODSTCY

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

         INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        INTEGER                GETNLIST
        INTEGER                INDEX1
        INTEGER                STR2INT

        EXTERNAL    CRLF, GETNLIST, INDEX1, STR2INT

C...........   SUBROUTINE ARGUMENTS
C...........   NOTE that NDROP and EDROP are not used at present
        INTEGER     , INTENT (IN) :: MXDATA  ! max no data variables allowed
        INTEGER     , INTENT (IN) :: MXIPPA  ! max no recognized data variables
        LOGICAL     , INTENT (IN) :: CFLAG   ! true: country header okay
        LOGICAL     , INTENT (IN) :: YFLAG   ! true: year header okay
        LOGICAL     , INTENT (IN) :: DFLAG   ! true: data names header okay
        CHARACTER(*), INTENT (IN) :: INVDNAM( MXIPPA ) ! valid data var names
        CHARACTER(*), INTENT(IN OUT):: LINE  ! full record in string
        INTEGER     , INTENT(OUT) :: ICC     ! country code
        INTEGER     , INTENT(OUT) :: INY     ! inventory year
        INTEGER     , INTENT(OUT) :: NPOA    ! no. data names
        INTEGER     , INTENT(OUT) :: EOS     ! error status

C...........   Local allocatable arrays
        CHARACTER(LEN=IOULEN3), ALLOCATABLE :: UNITS( : )

C...........   Other local variables
        INTEGER         I, J, K, L, L2, V  ! counters and indices

        INTEGER         COD       !  tmp data variable position in valid list
        INTEGER         IOS       !  i/o status
        INTEGER      :: NUNIT = 0 !  i/o status

        LOGICAL, SAVE:: FIRSTIME  = .TRUE.  !  first time subroutine is called
        LOGICAL, SAVE:: UNITFLAG  = .FALSE. !  true: units have already been set

        CHARACTER*1     CBUF    !  single char buffer
        CHARACTER*20    CNTRY   !  country name
        CHARACTER*300   MESG    !  message buffer

        CHARACTER(LEN=IOVLEN3) CVAR      ! tmp variable name
        CHARACTER(LEN=IOULEN3) UBUF      ! tmp units

        CHARACTER*16 :: PROGNAME = 'GETHDR' ! Program name

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

            L2 = LEN_TRIM( LINE )

C.............  Check for country name
            IF ( LINE(2:8) .EQ. 'COUNTRY' ) THEN 
                CNTRY = ADJUSTL( LINE( 9:L2 ) )
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

                ICC   = CTRYCOD( I )

C.............  Check for inventory year
            ELSE IF ( LINE(2:5) .EQ. 'YEAR' ) THEN 

                INY = STR2INT( LINE( 6:L2 ) )

                IF ( INY .LT. 1971 .OR. INY .GT. MXIYEAR ) THEN
                    WRITE( MESG, 94010 ) 'WARNING: Inventory year', INY,
     &                     'is outside of expected range.'
                    CALL M3MSG2( MESG )                    
                END IF

            ELSE IF ( LINE(2:6) .EQ. 'UNITS' ) THEN ! read in units                

                IF( UNITFLAG ) THEN
                    MESG = 'ERROR: UNITS field encountered again in ' //
     &                     'input. This is not allowed.'
                    CALL M3MSG2( MESG )
                    EOS = 5
                    RETURN

                ELSE IF( .NOT. ALLOCATED( POLPOS ) ) THEN
                    MESG = 'ERROR: UNITS field must be after DATA ' //
     &                     'or POLID field.'
                    CALL M3MSG2( MESG )
                    EOS = 5
                    RETURN

                END IF

C.................  Allocate memory based on the number of allowed data vars
                ALLOCATE( INVDUNT( MXIPPA,NPPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVDUNT', PROGNAME )
                ALLOCATE( INVDCNV( MXIPPA,NPPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVDCNV', PROGNAME )
                INVDUNT = ' '
                INVDCNV = 1.
                  
C.................  Get the number of units fields and allocate memory
                LINE = LINE( 7:L2 )
                L = LEN_TRIM( LINE )
                NUNIT = GETNLIST( L, LINE )

C.................  Error if units list is inconsistent with data fields
        	J = NUNIT/NPPOL
        	IF( NPOA  .GT. 0    .AND.
     &              NUNIT .GT. 0    .AND. 
     &              J     .NE. NPOA       ) THEN
                    WRITE( MESG,94010 ) 'ERROR:', J, 'group of units' //
     &                     'fields found, but', NPOA, 'variable ' //
     &                     'fields found in input file.'
                    CALL M3MSG2( MESG )
                    EOS = 6
                END IF

C.................  Allocate memory for local units field
                ALLOCATE( UNITS( NUNIT ), STAT=IOS )
                CALL CHECKMEM( IOS, 'UNITS', PROGNAME )
  
C.................  Parse the header line into the unit names
                CALL PARSLINE( LINE, NUNIT, UNITS )

C.................  Post-process units and store in final arrays
                K = 0
                DO V = 1, NPOA

                    I = POLPOS( V )
                    DO J = 1, NPPOL

                        K    = K + 1
                        UBUF = UNITS( K )

C.........................  Remove plural               
                	L  = INDEX( UBUF, '/' )
                	CBUF = UBUF( L-1:L-1 )
                	CALL UPCASE( CBUF )
                	IF( CBUF .EQ. 'S' ) THEN
                            L2 = LEN_TRIM( UBUF )
                            UBUF = UBUF(1:L-2) // UBUF(L:L2)
                	END IF

C.........................  Special case for units in 10E6
                	L = INDEX( UBUF, '10E6' )
                	IF( L .GT. 0 ) THEN
                            L2 = LEN_TRIM( UBUF )
                            INVDCNV( I,J ) = 1000000
                            UBUF = ADJUSTL( UBUF(L+4:L2) )
                	END IF

C.........................  Convert mile to mi
                	L  = INDEX( UBUF, 'mile' )
                	IF( L .GT. 0 ) THEN
                            L2 = LEN_TRIM( UBUF )
                            UBUF = UBUF(1:L+1) // UBUF(L+4:L2)
                	END IF

                        IF( I .GT. 0 ) INVDUNT( I,J ) = UBUF

                    END DO   ! No. variables per data variable
                END DO       ! No. variables

                DEALLOCATE( UNITS )

                UNITFLAG = .TRUE.

            ELSE IF ( LINE(2:6) .EQ. 'POLID' .OR.
     &                LINE(2:5) .EQ. 'DATA'       ) THEN ! read in data names

C..................... Deallocate names for pollutant, if must
                IF( ALLOCATED( TMPNAM ) ) DEALLOCATE( TMPNAM,POLPOS )

C.................  Allocate memory for current file for reading pol names
C                   and storing positions in master list
                LINE = LINE( 7:L2 )
                L = LEN_TRIM( LINE )
                NPOA = GETNLIST( L, LINE )

                ALLOCATE( TMPNAM( NPOA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'TMPNAM', PROGNAME )
                ALLOCATE( POLPOS( NPOA ), STAT=IOS )
                CALL CHECKMEM( IOS, 'POLPOS', PROGNAME )
                TMPNAM = ' '
                POLPOS = 0

C.................  Set special value of error status when the number of
C                   data variables attempted is greater than max allowed
                IF( NPOA .GT. MXDATA ) THEN
                    EOS = 4                  
                END IF

C.................  Parse the header line into the pollutant names
                CALL PARSLINE( LINE, NPOA, TMPNAM )

C.................  Store the position in master list of each pollutant
C.................  Write error if pollutant is not found.
                DO V = 1, NPOA

                    CVAR = TMPNAM( V )
                    L = LEN_TRIM( CVAR )
                    COD = INDEX1( CVAR, MXIPPA, INVDNAM )
                    IF( COD .LE. 0 ) THEN
                        EOS = 1
                        MESG = 'ERROR: Data variable "' // CVAR( 1:L )//
     &                         '" not in master data variable list!'
                        CALL M3MSG2( MESG )
                    ELSE
                        POLPOS( V ) = COD
                    END IF

                END DO

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
