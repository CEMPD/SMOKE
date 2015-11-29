
        SUBROUTINE GETDYSAV( NSRC, CIFIP, LDAYSAV )

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
C****************************************************************************

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS 
        CHARACTER(2)  CRLF
        LOGICAL      ENVYN
        INTEGER      FINDC
        INTEGER      GETFLINE
        INTEGER      PROMPTFFILE

        EXTERNAL     CRLF, ENVYN, GETFLINE, FINDC, PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS
        INTEGER           , INTENT (IN) :: NSRC             ! no. sources
        CHARACTER(FIPLEN3), INTENT (IN) :: CIFIP( NSRC )    ! country/state/co codes
        LOGICAL           , INTENT(OUT) :: LDAYSAV( NSRC )  ! true: source gets DLS

C...........   Local allocatable arrays
        CHARACTER(FIPLEN3), ALLOCATABLE :: EXMPTCSC( : )

C...........   Logical names and unit numbers
        INTEGER       FDEV   ! unit no. for daylight savings exemption file

C...........   Other local variables
        INTEGER       K, N, S       ! counters and indices

        INTEGER       IOS           ! i/o status
        INTEGER       IREC          ! record counter
        INTEGER       NEXEMPT       ! no. entries in the exemptions file
        INTEGER       NLINES        ! no. lines in input file

        LOGICAL    :: EFLAG = .TRUE. ! true: error found 
        LOGICAL    :: FFLAG = .TRUE. ! true: use daylight time exemption file 

        CHARACTER(FIPLEN3) CSC       ! tmp country/state/county code
        CHARACTER(FIPLEN3) STA       ! tmp country/state code
        CHARACTER(FIPLEN3) CNY       ! tmp country code
        CHARACTER(300) TRAILER       ! ending part of output log message
        CHARACTER(300) MESG          ! message buffer

        CHARACTER(16) :: PROGNAME = 'GETDYSAV' ! program name

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

        TRAILER = ' exempted from daylight savings time'

C.........  Loop through sources and see if any of the countries, states, or
C           counties do not use daylight time
        DO S = 1, NSRC

            CSC =  CIFIP( S )

C.............  Search for county
            K = FINDC( CSC,  NEXEMPT, EXMPTCSC )

            IF( K .GT. 0 ) THEN

                LDAYSAV( S ) = .FALSE.
                MESG = 'County ' // CSC // TRAILER
                CALL M3MESG( MESG )
                CYCLE

            END IF

C.............  Search for state
            STA = CSC
            STA( STALEN3+1:FIPLEN3 ) = REPEAT( '0', FIPLEN3-STALEN3+1 )
            K = FINDC( STA,  NEXEMPT, EXMPTCSC )

            IF( K .GT. 0 ) THEN

                LDAYSAV( S ) = .FALSE.
                MESG = 'State ' // STA // TRAILER
                CALL M3MESG( MESG )
                CYCLE

            END IF

C.............  Search for country
            CNY = CSC
            CNY( FIPEXPLEN3+2:FIPLEN3 ) = REPEAT( '0', FIPLEN3-FIPEXPLEN3+2 )
            K = FINDC( CNY, NEXEMPT, EXMPTCSC )

            IF( K .GT. 0 ) THEN

                LDAYSAV( S ) = .FALSE.
                MESG = 'Country ' // CNY // TRAILER
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
