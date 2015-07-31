
        SUBROUTINE FMTCSRC( INSTRING, NCIN, OUTBUFF, LENOUT )

C***********************************************************************
C  subroutine body starts at line 78
C
C  DESCRIPTION:
C      This subroutine formats the INSTRING assuming that is a string of 
C      source characteristics that has proper segment lengths.  It adds
C      labels and newlines as necessary and returns the formatted OUTBUFF
C      and its length.
C
C  PRECONDITIONS REQUIRED:
C      INSTRING with correct segement lengths
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Functions: I/O API functions
C
C  REVISION  HISTORY:
C      Created 11/98 by M. Houyoux
C
C**************************************************************************
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

C...........   Modules for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NCHARS, SC_BEGP, SC_ENDP

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS:
        CHARACTER(2)    CRLF

        EXTERNAL        CRLF

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: INSTRING ! Input string
        INTEGER     , INTENT (IN) :: NCIN     ! No. of chars to format (1 to 8)
        CHARACTER(*), INTENT(OUT) :: OUTBUFF  ! Formatted output string
        INTEGER     , INTENT(OUT) :: LENOUT   ! Length of output string

C...........   Local paramaters
        INTEGER, PARAMETER :: MXPTLBL = 8
        INTEGER, PARAMETER :: MXARLBL = 3
        INTEGER, PARAMETER :: MXMBLBL = 6

C.........  Output labels (note: these could be dynamic if in MODINFO)
        CHARACTER(6) :: PTLABEL( MXPTLBL ) =  !  message buffer
     &                ( / 'Region', 'Plant ', 'Char1 ', 'Char2 ',
     &                    'Char3 ', 'Char4 ', 'Char5 ', 'Data  ' / )
 
        CHARACTER(6) :: ARLABEL( MXARLBL ) =  !  message buffer
     &                ( / 'Region', 'SCC   ', 'Data  ' / )
 
        CHARACTER(9) :: MBLABEL( MXMBLBL ) =  !  message buffer
     &                ( / 'Region   ', 'Road Type', 'Link     ', 
     &                    'Vtype    ', 'SCC      ', 'Data     ' / )

        CHARACTER(9), SAVE :: LABEL( MXPTLBL )
 
C...........   Other local variables
        INTEGER         I, L, K, L1, L2, L3, L4  ! counters and indices

        INTEGER, SAVE :: TMPNUM      ! temporary number of chars
        INTEGER          NLOOP       ! number of loop iterations

        LOGICAL, SAVE :: FIRSTIME = .TRUE.

        CHARACTER(300)  BUFFER  !  string buffer

        CHARACTER(16) :: PROGNAME = 'FMTCSRC'  ! program name

C***********************************************************************
C   begin body of subroutine FMTCSRC

        IF( FIRSTIME ) THEN

            FIRSTIME = .FALSE.

            SELECT CASE( CATEGORY )
            CASE( 'AREA' )
                TMPNUM = MXARLBL
                LABEL( 1:MXARLBL ) = ARLABEL  ! array

            CASE( 'MOBILE' )
                TMPNUM = MXMBLBL
                LABEL( 1:MXMBLBL ) = MBLABEL  ! array

            CASE( 'POINT' )
                TMPNUM = MXPTLBL
                LABEL( 1:MXPTLBL ) = PTLABEL  ! array

            CASE DEFAULT
                TMPNUM = NCIN

            END SELECT

        END IF

C.........  Make sure not to exceed NCHARS legal value
        NLOOP = MIN( NCIN, TMPNUM )

C.........  Initialize output buffer
        OUTBUFF = ' '
        IF( NCHARS .GE. 1 ) THEN
            L1 = SC_BEGP( 1 )
            L2 = SC_ENDP( 1 )
            L4 = LEN_TRIM( LABEL( 1 ) )

            WRITE( OUTBUFF, 94900 ) LABEL( 1 )(1:L4), INSTRING( L1:L2 )

        ENDIF

C.........  Loop through the remaining source chars, writing the output string
C           for each populated source characterstic
        K = LEN_TRIM( OUTBUFF )
        DO I = 2, NLOOP

            L = LEN_TRIM( OUTBUFF )

            L1 = SC_BEGP( I )    ! retrieve stored field lengths
            L2 = SC_ENDP( I )
 
            BUFFER = ADJUSTL( INSTRING( L1:L2 ) )   ! Left-justifying
            L3 = LEN_TRIM( BUFFER )
            K = K + 8 + L3

C.............  Continue to contribute to buffer if not blank
            IF( BUFFER .NE. ' ' ) THEN

                L4 = LEN_TRIM( LABEL( I ) )

C.................  Include a return and indent if line gets too long
                IF( K .GT. EMOUTLN3 ) THEN
                    K = 18 + L3
                    WRITE( OUTBUFF, 94900 ) OUTBUFF( 1:L ) //
     &                   CRLF() // BLANK10 // LABEL( I )( 1:L4 ) ,
     &                     BUFFER( 1:L3 )
                ELSE
                    WRITE( OUTBUFF, 94900 ) OUTBUFF( 1:L ) // ' ' //
     &                     LABEL( I )( 1:L4 ) , BUFFER( 1:L3 )

                ENDIF
            ENDIF

        ENDDO

        LENOUT = LEN_TRIM( OUTBUFF )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94900   FORMAT( A, ': ', A )

        END





