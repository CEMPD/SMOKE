
        SUBROUTINE RDTZONE( FDEV, NDIM, NZS, NZF, TZONE0,
     &                      TZONST, TFIPST, TZONEF, TFIPEF )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the time zones file, sorts it, and returns
C      the sorted data
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
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

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)   CRLF

        EXTERNAL CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER      FDEV              !  time zones file unit number
        INTEGER      NDIM              !  dimension for time zone arrays
        INTEGER      NZS               !  no of state-specific
        INTEGER      NZF               !  no of county-specific
        INTEGER      TZONE0            !  fallback zone
        INTEGER      TZONST( NDIM )    !  state-specific time zones
        INTEGER      TFIPST( NDIM )    !  state FIPS codes (2 digit)
        INTEGER      TZONEF( NDIM )    !  FIPS-specific time zones
        INTEGER      TFIPEF( NDIM )    !  state/county FIPS codes (5 digits)

C...........   Unsorted time zone records
        INTEGER      INDXSA( NDIM )
        INTEGER      TZONSA( NDIM )
        INTEGER      TFIPSA( NDIM )
        INTEGER      INDXFA( NDIM )
        INTEGER      TZONFA( NDIM )
        INTEGER      TFIPFA( NDIM )

C...........   Other local variables
        INTEGER         FIP, TZONE          ! tmp fips code and time zone
        INTEGER         I, J                ! counters and indices
        INTEGER         IOS                 ! i/o status
        INTEGER         IREC                ! record counter

        LOGICAL      :: EFLAG = .FALSE.     ! error flag

        CHARACTER(300)  LINE    !  Input line from POINT file
        CHARACTER(300)  MESG    !  message buffer

        CHARACTER(16) :: PROGNAME = 'RDTZONE' ! program name

C***********************************************************************
C   begin body of subroutine RDTZONE

        TZONE0 = 5      !  default:  EST
        NZS    = 0
        NZF    = 0
        IREC   = 0
        DO      

            READ( FDEV, *, END=12, IOSTAT=IOS ) FIP, TZONE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 
     &                 'I/O error', IOS,  'reading time zones file ' //
     &                 'at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF ( FIP .EQ. 0 ) THEN              !  fallback -- all sources

                TZONE0 = TZONE

            ELSE IF ( MOD( FIP,100 ) .EQ. 0 ) THEN     !  state-specific zone

                NZS = NZS + 1
                IF ( NZS .LE. NDIM ) THEN
                    INDXSA( NZS ) = NZS
                    TFIPSA( NZS ) = FIP / 1000
                    TZONSA( NZS ) = TZONE
                END IF

            ELSE                                        !  county-specific zone

                NZF = NZF + 1
                IF ( NZF .LE. NDIM ) THEN
                    INDXFA( NZF ) = NZF
                    TFIPFA( NZF ) = FIP
                    TZONFA( NZF ) = TZONE
                END IF

            END IF      !  if fip zero, or nn000, or not.

        ENDDO

12      CONTINUE        !  exit from loop reading ZDEV

        IF ( NZS .GT. NDIM .OR. NZF .GT. NDIM ) THEN
            WRITE( MESG,94010 ) 
     &        'Number of state-only records  :', NZS, CRLF()// BLANK5//
     &        'Number of state&county records:', NZF, CRLF()// BLANK5//
     &        'Memory allocated:', NDIM
            CALL M3MSG2( MESG )

            MESG = 'INTERNAL ERROR: Insufficient memory allocated ' //
     &             'for time zones tables'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

        CALL SORTI1( NZS, INDXSA, TFIPSA )
        DO I = 1, NZS
            J = INDXSA( I )
            TZONST( I ) = TZONSA( J )
            TFIPST( I ) = TFIPSA( J )
        ENDDO

        CALL SORTI1( NZF, INDXFA, TFIPFA )
        DO I = 1, NZF
            J = INDXFA( I )
            TZONEF( I ) = TZONFA( J )
            TFIPEF( I ) = TFIPFA( J )
        ENDDO

C.........  Rewind file

        REWIND( FDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I5.5 )

        END
