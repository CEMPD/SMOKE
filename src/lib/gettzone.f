
        INTEGER FUNCTION GETTZONE( FIP, NZS, NZF, TZONE0,
     &                             TZONST, TFIPST, TZONEF, TFIPEF )

C***********************************************************************
C  function body starts at line 76
C
C  DESCRIPTION: 
C     Returns the time zone of a FIPs code based on the input tables
C
C  PRECONDITIONS REQUIRED:
C     FIP set to state and county FIPS code of interest
C     Memory allocated for arrays in argument list and arrays populated
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C     Subroutines: I/O API subroutines
C     Functions: I/O API functions
C
C  REVISION  HISTORY:
C     Created 10/98 by M. Houyoux
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

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER      FIND1   !  returns -1 for failure
        EXTERNAL     FIND1

C...........   ARGUMENTS and their descriptions:

        INTEGER      FIP               !  input FIPS code to assign zone
        INTEGER      NZS               !  no of state-specific
        INTEGER      NZF               !  no of county-specific
        INTEGER      TZONE0            !  fallback zone
        INTEGER      TZONST( NZS )     !  state-specific time zones
        INTEGER      TFIPST( NZS )     !  state FIPS codes (2 digit)
        INTEGER      TZONEF( NZF )     !  FIPS-specific time zones
        INTEGER      TFIPEF( NZF )     !  state/county FIPS codes (5 digits)

C...........   LOCAL VARIABLES their descriptions:

        INTEGER       K            ! indice

        CHARACTER*1   BUFFER       ! ASCII LINE from X-ref file
        CHARACTER*300 MESG         ! Message buffer

        CHARACTER*16 :: PROGNAME = 'GETTZONE' ! Program name

C***********************************************************************
C   begin body of subroutine GETTZONE

C.........  Search for FIPS code in county-specific table
        K = FIND1( FIP, NZF, TFIPEF )

        IF ( K .GT. 0 ) THEN
            GETTZONE = TZONEF( K )

C.........  Search for FIPS code in state-specific table
        ELSE
            K = FIND1( FIP/1000, NZS, TFIPST )

            IF ( K .GT. 0 ) THEN
                GETTZONE = TZONST( K )

C.............  Apply default
            ELSE
                GETTZONE = TZONE0

                WRITE( MESG,94040 ) 
     &                'WARNING: Applying default time zone of', TZONE0,
     &                'to state/county FIPS code', FIP
                CALL M3MESG( MESG )

            END IF

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************C...........   Internal buffering formats.............94xxx
 
C...........   Internal buffering formats............ 94xxx

94040   FORMAT( A, 1X, I2.2, 1X, A, 1X, I6.6 )

        END
