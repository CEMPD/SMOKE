
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
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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
        
        INCLUDE 'IOCNST3.EXT'
        INCLUDE 'IOSTRG3.EXT'
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
        LOGICAL      CHKINT
        LOGICAL      CHKREAL
        CHARACTER*2  CRLF
        INTEGER      GETEFILE
        INTEGER      GETNLIST
        INTEGER      INDEX1
        INTEGER      STR2INT
        REAL         STR2REAL

        EXTERNAL     CHKINT, CHKREAL, CRLF, GETEFILE, GETNLIST, INDEX1, 
     &               STR2INT, STR2REAL
        
C...........   Local parameters
        INTEGER, PARAMETER :: MXGRDTYP = 5
        INTEGER, PARAMETER :: M3CHAR = 7

C...........   Grid types and names arrays
        INTEGER      :: GRDTYPES( MXGRDTYP ) = ( / LATGRD3
     &                                           , LAMGRD3
     &                                           , MERGRD3
     &                                           , STEGRD3
     &                                           , UTMGRD3 / )

        CHARACTER*15 :: GRDNAMES( MXGRDTYP ) = ( / 'LAT-LON        '
     &                                           , 'LAMBERT        '
     &                                           , 'MERCATOR       '
     &                                           , 'STEREOGRAPHIC  '
     &                                           , 'UTM            ' / )

C...........   Local arrays (note- tried to make these allocatable, but
C              this caused unexplainable crashing on SGI).
        CHARACTER*32 :: SEGMENT( 32 )
        CHARACTER*32 :: UPCSGMT( 32 )

C...........   File units and logical/physical names:
        INTEGER, SAVE :: IDEV    !  unit number of grid information file
        CHARACTER*16     LNAME   !  logical name for grid information file

C...........   Scratch local variables and their descriptions:
            
        INTEGER         I, J, L, L2, N  !  indices and string lengths
        INTEGER         IOS      !  I/O status return
        INTEGER      :: IDUM = 0 !  integer dummy variable
        INTEGER         IREC     !  record counter
        INTEGER         LEGAL    !  valid length of string

        REAL*8       :: RDUM = 0 !  double precision dummy variable
        REAL*8          DMISS3   !  double precision missing value

        LOGICAL      :: CFLAG = .FALSE.   ! true: time to cycle in loop
        LOGICAL      :: EFLAG = .FALSE.   ! true: error detected
        LOGICAL      :: FIRSTIME = .TRUE. ! true: first time routine called

        CHARACTER*16    CDUM     !  character dummy buffer
        CHARACTER*300   BUFFER   !  multi-purpose buffer
        CHARACTER*300   GNBUF    !  grid name buffer
        CHARACTER*300   GDBUF    !  grid description buffer
        CHARACTER*300   LINE     !  line of file
        CHARACTER*300   MESG     !  message buffer
        CHARACTER*300   UPCLINE  !  upper case line of file

        CHARACTER*16 :: PROGNAME = 'DSCM3GRD' ! program name

C***********************************************************************
C   begin body of function DSCM3GRD

C.........  Get value of Models-3 environment variable for grid and check
C           exit status
        IF( FIRSTIME ) THEN

            FIRSTIME = .FALSE.
            LNAME = 'G_GRIDPATH'

C.............  Open grid file
            IDEV = GETEFILE( LNAME, .TRUE., .TRUE., PROGNAME )

            IF( IDEV .LE. 0 ) THEN

                L = LEN_TRIM( LNAME )
                MESG = 'Could not open file "' // LNAME( 1:L ) // '"!'
                DSCM3GRD = .FALSE.
                RETURN

            END IF
        END IF

C.........  Initialize grid setting to missing
        DMISS3  = DBLE( AMISS3 )
        GNAME	= ' '
        GDESC   = ' '
        CNAME	= ' '
        CTYPE	= IMISS3
        PUNIT   = ' '
        P_ALP	= DMISS3
        P_BET	= DMISS3
        P_GAM	= DMISS3
        XCENT	= DMISS3
        YCENT	= DMISS3
        XORIG	= DMISS3
        YORIG	= DMISS3
        XCELL	= DMISS3
        YCELL	= DMISS3
        NCOLS	= IMISS3
        NROWS	= IMISS3
        NTHIK	= IMISS3

C.........  Make sure read it at the start of the file
        REWIND( IDEV )

C.........  Read grid information file
        IREC = 0
        DO 

C.............  Read whole line from the file
            READ( IDEV, 93000, END=111, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS,
     &                 'reading grid information file at line', IREC

                CALL M3MESG( MESG )
                CYCLE
            ENDIF

C.............  Adjust line to left and create upper case line
            LINE = ADJUSTL( LINE )
            UPCLINE = LINE
            CALL UPCASE( UPCLINE )
            L2 = LEN_TRIM( LINE )

C.............  Skip line if it starts with a number
            IF( UPCLINE( 1:1 ) .GE. '0' .AND.
     &          UPCLINE( 1:1 ) .LE. '9'       ) CYCLE

C.............  Parse the lower and upper case lines
            N = GETNLIST( L2, UPCLINE )
            IF( N .GT. 32 ) THEN
                WRITE( MESG,94010 ) 'WARNING: Line', IREC, 
     &                 'in G_GRIDPATH file skipped because more ' //
     &                 'than 32 fields were found.'
                CALL M3MSG2( MESG )
                CYCLE 
            END IF

            SEGMENT = ' '            
            UPCSGMT = ' '            
            CALL PARSLINE( LINE   , L2, SEGMENT )
            CALL PARSLINE( UPCLINE, L2, UPCSGMT )

C.............  Search for keywords.  If the keyword exists, extract the
C               value(s) and cycle to the next line...

C.............  Grid name
            SELECT CASE ( UPCSGMT( 1 ) )
            CASE( 'GDNAME_GD' )
                GNBUF = UPCSGMT( 3 )

C.............  Grid description
            CASE( 'GDDESC_GD' )
                I = INDEX( LINE, '0' )
                GDBUF = ADJUSTL( LINE( I+1:L2 ) )

C.............  Coordinate system type
            CASE( 'GDTYP_GD' )
                CNAME = UPCSGMT( 3 )

C.................  Look for coordinate system type in known types
                J = INDEX1( CNAME, MXGRDTYP, GRDNAMES )
                IF( J .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Coordinate system type "' // 
     &                      CNAME // '" is not known in grid ' //
     &                      'information reader.'
                    CALL M3MSG2( MESG )
                ELSE
                    CTYPE = GRDTYPES( J )
                END IF

C.............  Grid units
            CASE( 'GDUNT_GD' )
                PUNIT = SEGMENT( 3 )

C.............  Alpha value
            CASE( 'P_ALP_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, P_ALP, IDUM )

C.............  Beta value
            CASE( 'P_BET_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, P_BET, IDUM )

C.............  Gamma value
            CASE( 'P_GAM_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, P_GAM, IDUM )

C.............  X-center value
            CASE( 'XCENT_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, XCENT, IDUM )

C.............  Y-center value
            CASE( 'YCENT_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, YCENT, IDUM )

C.............  X-origin value
            CASE( 'XORIG_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, XORIG, IDUM )

C.............  Y-origin value
            CASE( 'YORIG_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, YORIG, IDUM )

C.............  Delta X value
            CASE( 'XCELL_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, XCELL, IDUM )

C.............  Delta Y value
            CASE( 'YCELL_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, YCELL, IDUM )

C.............  Number of columns
            CASE( 'NCOLS' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3INT, RDUM, NCOLS )

C.............  Number of rows
            CASE( 'NROWS' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3INT, RDUM, NROWS )

C.............  Number of boundary cells
            CASE( 'NTHIK' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3INT, RDUM, NTHIK )

            END SELECT

        END DO

C.........  Exit from read loop
111     CONTINUE

C.........  Check length of character variables that matter
        L = LEN_TRIM( GNBUF )
        IF ( L .GT. IOVLEN3 ) THEN
 
            WRITE( MESG,94010 ) 'Grid name "' // GNBUF( 1:L ) //
     &             '" has maximum allowable length of', IOVLEN3, '.' //
     &             CRLF() // BLANK10 // 'Truncating to "' // 
     &             GNBUF( 1:IOVLEN3 ) // '".'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

        END IF
        GNAME = GNBUF( 1:IOVLEN3 )

C.........  Set grid description, and truncate
        GDESC = GDBUF( 1:LEN( GDESC ) )

C.........  Make sure everything important has been defined
        CALL GRDINFO_DEFINED( 'GDNAME_GD', GNAME, RDUM, IDUM, M3CHAR )
        CALL GRDINFO_DEFINED( 'GDDESC_GD', GDESC, RDUM, IDUM, M3CHAR )
        CALL GRDINFO_DEFINED( 'GDTYP_GD', CNAME, RDUM, IDUM, M3CHAR )
        CALL GRDINFO_DEFINED( 'P_ALP_GD', CDUM, P_ALP, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'P_BET_GD', CDUM, P_BET, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'P_GAM_GD', CDUM, P_GAM, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'XCENT_GD', CDUM, XCENT, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'YCENT_GD', CDUM, YCENT, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'XORIG_GD', CDUM, XORIG, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'YORIG_GD', CDUM, YORIG, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'XCELL_GD', CDUM, XCELL, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'YCELL_GD', CDUM, YCELL, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'NCOLS', CDUM, RDUM, NCOLS, M3INT )
        CALL GRDINFO_DEFINED( 'NROWS', CDUM, RDUM, NROWS, M3INT )
        CALL GRDINFO_DEFINED( 'NTHIK', CDUM, RDUM, NTHIK, M3INT )
 
        DSCM3GRD = ( .NOT. EFLAG )

        RETURN
        
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx                

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

C********************** INTERNAL SUBPROGRAMS ****************************

        CONTAINS

            SUBROUTINE GRDINFO_CHECKER( STRINGS, OUTTYPE,
     &                                  ROUT, IOUT )

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN) :: STRINGS( 3 )
            INTEGER     , INTENT (IN) :: OUTTYPE
            REAL*8      , INTENT(OUT) :: ROUT     ! real output value
            INTEGER     , INTENT(OUT) :: IOUT     ! integer output value

C.............  Local variables
            INTEGER   I, L

            CHARACTER*300 MESG

C..........................................................................

            IF( OUTTYPE .EQ. M3INT ) THEN

                IF( .NOT. CHKINT( STRINGS( 3 ) ) ) THEN
                    EFLAG = .TRUE.
                    L = LEN_TRIM( STRINGS( 1 ) )
                    MESG = 'ERROR: ' // STRINGS( 1 )( 1:L ) // 
     &                     ' value is invalid in grid information' // 
     &                     ' file.'
                    CALL M3MSG2( MESG )
                ELSE
                    IOUT = STR2INT( STRINGS( 3 ) )
                END IF

            ELSE IF( OUTTYPE .EQ. M3DBLE ) THEN

                IF( .NOT. CHKREAL( STRINGS( 3 ) ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: ' // STRINGS( 1 )( 1:L ) // 
     &                     ' value is invalid in grid information' // 
     &                     ' file.'
                    CALL M3MSG2( MESG )
                ELSE
                    ROUT = DBLE( STR2REAL( STRINGS( 3 ) ) )
                END IF

            END IF

            RETURN

            END SUBROUTINE GRDINFO_CHECKER

C..........................................................................
C..........................................................................

            SUBROUTINE GRDINFO_DEFINED( KEYWORD, CVAL, RVAL, 
     &                                  IVAL, INTYPE )

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN) :: KEYWORD
            CHARACTER(*), INTENT (IN) :: CVAL
            REAL*8      , INTENT (IN) :: RVAL
            INTEGER     , INTENT (IN) :: IVAL
            INTEGER     , INTENT (IN) :: INTYPE

C.............  Local variables
            CHARACTER*300 MESG

C..........................................................................

            IF( ( INTYPE .EQ. M3CHAR .AND. CVAL .EQ. ' '    ) .OR.
     &          ( INTYPE .EQ. M3DBLE .AND. RVAL .EQ. DMISS3 ) .OR.
     &          ( INTYPE .EQ. M3INT  .AND. IVAL .EQ. IMISS3 ) ) THEN

                EFLAG = .TRUE.
                MESG = 'ERROR: Keyword "' // KEYWORD // 
     &                 '" not found in grid information file.'
                CALL M3MSG2( MESG )

            END IF

            END SUBROUTINE GRDINFO_DEFINED

        END FUNCTION DSCM3GRD

