
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
        LOGICAL      CHKINT
        LOGICAL      CHKREAL
        CHARACTER*2  CRLF
        INTEGER      GETEFILE
        INTEGER      INDEX1
        INTEGER      STR2INT
        REAL         STR2REAL

        EXTERNAL     CHKINT, CHKREAL, CRLF, GETEFILE, INDEX1, 
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

C...........   File units and logical/physical names:
        INTEGER, SAVE :: IDEV    !  unit number of grid information file
        CHARACTER*16     LNAME   !  logical name for grid information file

C...........   Scratch local variables and their descriptions:
            
        INTEGER         I, J, L, L2  !  indices and string lengths
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

            L2 = LEN_TRIM( LINE )

C.............  Search for keywords.  If the keyword exists, extract the
C               value(s) and cycle to the next line...

C.............  Grid name
            I = INDEX( LINE, 'GDNAME_GD' )
            IF( I .GT. 0 ) THEN
                I = INDEX( LINE, '0' )
                GNBUF = ADJUSTL( LINE( I+1:L2 ) )
                CALL UPCASE( GNBUF )
                CYCLE
            END IF

C.............  Grid description
            I = INDEX( LINE, 'GDDESC_GD' )
            IF( I .GT. 0 ) THEN
                I = INDEX( LINE, '0' )
                GDBUF = ADJUSTL( LINE( I+1:L2 ) )
                CALL UPCASE( GDBUF )
                CYCLE
            END IF

C.............  Coordinate system type
            I = INDEX( LINE, 'GDTYP_GD' )
            IF( I .GT. 0 ) THEN
                I = INDEX( LINE, '0' )
                CNAME = ADJUSTL( LINE( I+1:L2 ) )
                CALL UPCASE( CNAME )

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
                CYCLE
            END IF

C.............  Grid units
            I = INDEX( LINE, 'GDUNT_GD' )
            IF( I .GT. 0 ) THEN
                I = INDEX( LINE, '0' )
                PUNIT = ADJUSTL( LINE( I+1:L2 ) )
                CALL UPCASE( PUNIT )
                CYCLE
            END IF

C.............  Alpha value
            CALL GRDINFO_CHECKER( 'P_ALP_GD', M3DBLE, P_ALP, 
     &                            IDUM, CFLAG )
            IF( CFLAG ) CYCLE

C.............  Beta value
            CALL GRDINFO_CHECKER( 'P_BET_GD', M3DBLE, P_BET, 
     &                            IDUM, CFLAG )
            IF( CFLAG ) CYCLE

C.............  Gamma value
            CALL GRDINFO_CHECKER( 'P_GAM_GD', M3DBLE, P_GAM, 
     &                            IDUM, CFLAG )
            IF( CFLAG ) CYCLE

C.............  X-center value
            CALL GRDINFO_CHECKER( 'XCENT_GD', M3DBLE, XCENT, 
     &                            IDUM, CFLAG )
            IF( CFLAG ) CYCLE

C.............  Y-center value
            CALL GRDINFO_CHECKER( 'YCENT_GD', M3DBLE, YCENT, 
     &                            IDUM, CFLAG )
            IF( CFLAG ) CYCLE

C.............  X-origin value
            CALL GRDINFO_CHECKER( 'XORIG_GD', M3DBLE, XORIG, 
     &                            IDUM, CFLAG )
            IF( CFLAG ) CYCLE

C.............  Y-origin value
            CALL GRDINFO_CHECKER( 'YORIG_GD', M3DBLE, YORIG, 
     &                            IDUM, CFLAG )
            IF( CFLAG ) CYCLE

C.............  Delta X value
            CALL GRDINFO_CHECKER( 'XCELL_GD', M3DBLE, XCELL, 
     &                            IDUM, CFLAG )
            IF( CFLAG ) CYCLE

C.............  Delta Y value
            CALL GRDINFO_CHECKER( 'YCELL_GD', M3DBLE, YCELL, 
     &                            IDUM, CFLAG )
            IF( CFLAG ) CYCLE

C.............  Number of columns
            CALL GRDINFO_CHECKER( 'NCOLS', M3INT, RDUM, 
     &                            NCOLS, CFLAG )
            IF( CFLAG ) CYCLE

C.............  Number of rows
            CALL GRDINFO_CHECKER( 'NROWS', M3INT, RDUM, 
     &                            NROWS, CFLAG )
            IF( CFLAG ) CYCLE

C.............  Number of boundary cells
            CALL GRDINFO_CHECKER( 'NTHIK', M3INT, RDUM, 
     &                            NTHIK, CFLAG )
            IF( CFLAG ) CYCLE

        END DO

C.........  Exit from read loop
111     CONTINUE

C.........  Check length of character variables that matter
        BUFFER = ADJUSTL( GNBUF )
        L = LEN_TRIM( BUFFER )
        LEGAL = LEN( GNAME )
        IF ( L .GT. LEGAL ) THEN
 
            WRITE( MESG,94010 ) 'Grid name "', BUFFER( 1:L ), 
     &             '" Has maximum allowable length of', LEGAL, '.' //
     &             CRLF() // BLANK10 // 'Truncating to "' // 
     &             BUFFER( 1:LEGAL ) // '"'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

        END IF
        GNAME = BUFFER( 1:LEGAL )

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

            SUBROUTINE GRDINFO_CHECKER( KEYWORD, OUTTYPE,
     &                                  ROUT, IOUT, CFLAG )

C.............  Subroutine arguments
            CHARACTER(*), INTENT (IN) :: KEYWORD
            INTEGER     , INTENT (IN) :: OUTTYPE
            REAL*8      , INTENT(OUT) :: ROUT     ! real output value
            INTEGER     , INTENT(OUT) :: IOUT     ! integer output value
            LOGICAL     , INTENT(OUT) :: CFLAG    ! true: keyword was found

C.............  Local variables
            INTEGER   I, L

            CHARACTER*300 MESG

C..........................................................................

            CFLAG = .FALSE.
            I = INDEX( LINE, KEYWORD )

            IF( I .GT. 0 ) THEN

                I = INDEX( LINE, '0' )
                IF( OUTTYPE .EQ. M3INT ) THEN

                    IF( .NOT. CHKINT( LINE( I+1:L2 ) ) ) THEN
                        EFLAG = .TRUE.
                        MESG = 'ERROR: ' // KEYWORD // ' value is ' // 
     &                         'invalid in grid information file.'
                        CALL M3MSG2( MESG )
                    ELSE
                        IOUT = STR2INT( LINE( I+1:L2 ) )
                    END IF

                ELSE IF( OUTTYPE .EQ. M3DBLE ) THEN

                    IF( .NOT. CHKREAL( LINE( I+1:L2 ) ) ) THEN
                        EFLAG = .TRUE.
                        MESG = 'ERROR: ' // KEYWORD // ' value is ' // 
     &                         'invalid in grid information file.'
                        CALL M3MSG2( MESG )
                    ELSE
                        ROUT = DBLE( STR2REAL( LINE( I+1:L2 ) ) )
                    END IF

                END IF

                CFLAG = .TRUE.

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

