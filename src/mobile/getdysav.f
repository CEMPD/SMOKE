
        SUBROUTINE GETDYSAV( NSRC, IFIP, LDAYSAV )

C***********************************************************************
C  subroutine GETDYSAV body starts at line < >
C
C  DESCRIPTION:
C
C  PRECONDITIONS REQUIRED:
C      Creates an array of sources that are affected by daylight savings time.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C 
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER*2  CRLF
        LOGICAL      ENVYN
        INTEGER      FIND1
        INTEGER      GETFLINE
        INTEGER      PROMPTFFILE

        EXTERNAL     CRLF, ENVYN, GETFLINE, FIND1, PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NSRC             ! no. sources
        INTEGER     , INTENT (IN) :: IFIP   ( NSRC )  ! country/state/co codes
        LOGICAL     , INTENT(OUT) :: LDAYSAV( NSRC )  ! true: source gets DLS

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE :: EXMPTCSC( : )

C...........   Logical names and unit numbers
        INTEGER       FDEV   ! unit no. for daylight savings exemption file

C...........   Other local variables
        INTEGER       K, N, S       ! counters and indices

        INTEGER       CNY           ! tmp country code
        INTEGER       CSC           ! tmp country/state/county code
        INTEGER       IOS           ! i/o status
        INTEGER       IREC          ! record counter
        INTEGER       NEXEMPT       ! no. entries in the exemptions file
        INTEGER       NLINES        ! no. lines in input file
        INTEGER       STA           ! tmp country/state code

        LOGICAL    :: EFLAG = .TRUE. ! true: error found 
        LOGICAL    :: FFLAG = .TRUE. ! true: use daylight time exemption file 

        CHARACTER*300 TRAILER       ! ending part of output log message
        CHARACTER*300 MESG          ! message buffer

        CHARACTER*16 :: PROGNAME = 'GETDYSAV' ! program name

C***********************************************************************
C   begin body of subroutine GETDYSAV

C.........  Initialize daylight savings indicator array
        LDAYSAV = .TRUE.  ! array

C.........  Get environment variable to decide whether file prompt is
C           needed
        MESG = 'Daylight time exemption switch'
        FFLAG = ENVYN( 'DAYLIGHT_EXEMPT', MESG, .FALSE., IOS )

C.........  Prompt for and open daylight savings indicator file
        IF( FFLAG ) THEN
            FDEV = PROMPTFFILE( 
     &             'Enter logical name for DAYLIGHT TIME EXEMPT file',
     &             .TRUE., .TRUE., 'DLTEXEMPT', PROGNAME )

C.........  Return if no file was needed
        ELSE
            RETURN

        END IF

C.........  Get number of lines in the exemptions file
        NLINES = GETFLINE( FDEV, 'Daylight exemptions file' )

C.........  Allocate memory for the temporary arrays for reading in file
        ALLOCATE( EXMPTCSC( NLINES ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EXMPTCSC', PROGNAME )

C.........  Loop through file and read it
        DO N = 1, NLINES

            READ( FDEV, *, END=999, IOSTAT = IOS ) CSC

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 'reading daylight ' //
     &              'savings exemptions file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF
                        
            EXMPTCSC( N ) = CSC

        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem reading daylight savings exemptions file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        NEXEMPT = N

        TRAILER = 'exempted from daylight savings time'

C.........  Loop through sources and see if any of the countries, states, or
C           counties do not use daylight time
        DO S = 1, NSRC

            CSC =  IFIP( S )

C.............  Search for county
            K = FIND1( CSC,  NEXEMPT, EXMPTCSC )

            IF( K .GT. 0 ) THEN

                LDAYSAV( S ) = .FALSE.
                WRITE( MESG,94010 ) 'County', CSC, TRAILER
                CALL M3MESG( MESG )
                CYCLE

            END IF

C.............  Search for state
            STA = 1000* ( CSC / 1000 )
            K = FIND1( STA,  NEXEMPT, EXMPTCSC )

            IF( K .GT. 0 ) THEN

                LDAYSAV( S ) = .FALSE.
                WRITE( MESG,94010 ) 'State', STA, TRAILER
                CALL M3MESG( MESG )
                CYCLE

            END IF

C.............  Search for country
            CNY = 100000 * ( CSC / 100000 )
            K = FIND1( CNY, NEXEMPT, EXMPTCSC )

            IF( K .GT. 0 ) THEN

                LDAYSAV( S ) = .FALSE.
                WRITE( MESG,94010 ) 'Country', CNY, TRAILER
                CALL M3MESG( MESG )
                CYCLE

            END IF

        END DO

C.........  Deallocate local memory
        DEALLOCATE( EXMPTCSC )

        RETURN

C.........  Error message for reaching the end of file too soon
999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of daylight' // CRLF() // BLANK5 //
     &         'savings exemptions file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94020   FORMAT( 10 ( A, 1X, I5, 2( 1X, A, 1X, F8.2 ) ) )
 
        END SUBROUTINE GETDYSAV
