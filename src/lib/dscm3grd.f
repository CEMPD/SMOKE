
        LOGICAL FUNCTION DSCM3GRD( GNAME, GDESC, CNAME, CTYPE, PUNIT, 
     &                             P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &                             XORIG, YORIG, XCELL, YCELL, NCOLS, 
     &                             NROWS, NTHIK )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function returns true it successfully reads the properties of the
C      models-3 grid, and false otherwise.  It also returns those properties
C      using the subroutine arguments.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Copied and modified from DSCGRID, version 12/16/97
C
C**************************************************************************
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
C***************************************************************************

        IMPLICIT NONE
        
        INCLUDE 'EMCNST3.EXT'               
        INCLUDE 'PARMS3.EXT'               
        
C...........   ARGUMENTS and their descriptions.  All are output variables:
        
        CHARACTER*(*) GNAME	!  grid name
        CHARACTER*(*) GDESC     !  grid description
        CHARACTER*(*) CNAME	!  coord sys name
        INTEGER       CTYPE	!  coord sys type (I/O API code number)
        CHARACTER*(*) PUNIT     !  projection units (e.g., meters or degrees)
        REAL*8        P_ALP	!  first, second, third map
        REAL*8        P_BET	!  projection descriptive
        REAL*8        P_GAM	!  parameters
        REAL*8        XCENT	!  lon for coord-system X=0
        REAL*8        YCENT	!  lat for coord-system Y=0
        REAL*8        XORIG	!  X-coordinate origin of grid (map units)
        REAL*8        YORIG	!  Y-coordinate origin of grid
        REAL*8        XCELL	!  X-coordinate cell dimension
        REAL*8        YCELL	!  Y-coordinate cell dimension
        INTEGER       NCOLS	!  number of grid columns
        INTEGER       NROWS	!  number of grid rows
        INTEGER       NTHIK	!  BOUNDARY:  perimeter thickness (cells)
                    
C...........   EXTERNAL FUNCTIONS:
        
        CHARACTER*2  CRLF
        INTEGER      GETEFILE
        INTEGER      INDEX1
        INTEGER      TRIMLEN

        EXTERNAL     CRLF, GETEFILE, INDEX1, TRIMLEN
        
C...........   Local parameters
        INTEGER, PARAMETER :: MXGRDTYP = 5

C...........   Grid types and names arrays
        INTEGER      :: GRDTYPES( MXGRDTYP ) = ( / LATGRD3
     &                                           , LAMGRD3
     &                                           , MERGRD3
     &                                           , STEGRD3
     &                                           , UTMGRD3 / )

        CHARACTER*15 :: GRDNAMES( MXGRDTYP ) = ( / 'GEOGRAPHIC     '
     &                                           , 'LAMBERT        '
     &                                           , 'MERCATOR       '
     &                                           , 'STEREOGRAPHIC  '
     &                                           , 'UTM            ' / )

C...........   File units and logical/physical names:
        INTEGER         IDEV    !  unit number of grid information file
        CHARACTER*16    LNAME   !  logical name for grid information file

C...........   Scratch local variables and their descriptions:
            
        INTEGER         I, J, L  !  indexes for names in lists
        INTEGER         IOS      !  I/O status return
        INTEGER         LEGAL    !  actual length of GNAME or GDESC from caller

        CHARACTER*16    CDUM     !  character dummy buffer
        CHARACTER*300   BUFFER   !  multi-purpose buffer
        CHARACTER*300   GNBUF    !  grid name buffer
        CHARACTER*300   GDBUF    !  grid description buffer
        CHARACTER*300   MESG     !  message buffer

        CHARACTER*16 :: PROGNAME = 'DSCM3GRD' ! program name

C***********************************************************************
C   begin body of function DSCM3GRD

C.........  Get value of Models-3 environment variable for grid and check
C           exit status

        LNAME = 'G_GRIDPATH'

C.........  Open grid file
        IDEV = GETEFILE( LNAME, .TRUE., .TRUE., PROGNAME )

        IF( IDEV .LE. 0 ) THEN

            L = TRIMLEN( LNAME )
            MESG = 'Could not open file "' // LNAME( 1:L ) // '"!'
            DSCM3GRD = .FALSE.
            RETURN

        ENDIF

C.........  Read grid information file

        READ( IDEV, 93010, END=111, ERR=222 ) 
     &        GNBUF, GDBUF, CNAME, PUNIT, P_ALP, P_BET, P_GAM, 
     &        XCENT, YCENT, XORIG, YORIG, CDUM, XCELL, YCELL,
     &        NCOLS, NROWS, NTHIK

C.........  Check length of input strings (that matter) with allowable length 
C           as defined in the calling program.

        BUFFER = ADJUSTL( GNBUF )
        L = TRIMLEN( BUFFER )
        LEGAL = LEN( GNAME )
        IF ( L .GT. LEGAL ) THEN

            WRITE( MESG,94010 ) 'Grid name"', BUFFER( 1:L ), 
     &             '" Has maximum allowable length of', LEGAL, '.' //
     &             CRLF() // BLANK5 // 'Truncating to "' // 
     &             BUFFER( 1:LEGAL ) // '"'
            CALL M3WARN( 'DSCM3GRD', 0, 0, MESG )

        END IF

        GNAME = BUFFER( 1:LEGAL )
                
        BUFFER = ADJUSTL( GDBUF )
        L = TRIMLEN( BUFFER )
        LEGAL = LEN( GDESC )
        IF ( L .GT. LEGAL ) THEN

            WRITE( MESG,94010 ) 'Grid description being truncated to',
     &              LEGAL, 'allowable characters.'
            CALL M3WARN( 'DSCM3GRD', 0, 0, MESG )

        END IF
                
        GDESC = BUFFER( 1:LEGAL )

C.........  Convert the coordinate system type to I/O API code

        CALL UPCASE( CNAME )

        I = INDEX1( CNAME, MXGRDTYP, GRDNAMES )

        IF( I .GT. 0 ) THEN
            CTYPE = GRDTYPES( I )

        ELSE
            L = TRIMLEN( CNAME )
            MESG = 'INTERNAL ERROR: Coordinate system type "' //
     &             CNAME( 1:L ) // '" is not regonized by ' // CRLF() //
     &             BLANK5// 'library function "'// PROGNAME(1:16) // '"'
            CALL M3MSG2( MESG )
            DSCM3GRD = .FALSE.
            RETURN

        ENDIF

        DSCM3GRD = .TRUE.

        CLOSE( IDEV )

        RETURN

111     MESG = 'ERROR: End of file "' // LNAME( 1:TRIMLEN( LNAME ) ) //
     &         '" reached too early.'
        CALL M3MSG2( MESG )

222     DSCM3GRD = .FALSE.

        CLOSE( IDEV )

        RETURN
        
C******************  FORMAT  STATEMENTS   ******************************
                
C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( 4(A/), 7(E10.2/), A/ 2(E10.2/), 2(I8/), I8 )
                
C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )


        END

