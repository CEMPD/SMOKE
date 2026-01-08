
        LOGICAL FUNCTION DSCM3LAY( NLAYS, VGTYP, VGTOP, VGLVS )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function returns true it successfully reads the properties of the
C      models-3 layer structure, and false otherwise.  It also returns those 
C      properties using the subroutine arguments.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Copied and modified from DSCM3LAY, version 1.10
C
C**************************************************************************
C
C Project Title: EDSS Tools Library
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************
        USE M3UTILIO

        IMPLICIT NONE
        
        INCLUDE 'IOCNST3.EXT'
C        INCLUDE 'PARMS3.EXT'               
        
C...........   ARGUMENTS and their descriptions.  All are output variables:
        
        INTEGER, INTENT( OUT ) :: NLAYS      ! number of layers
        INTEGER, INTENT( OUT ) :: VGTYP      ! type of vertical coordinates
        REAL   , INTENT( OUT ) :: VGTOP      ! model-top, for sigma coord types
        REAL   , INTENT( OUT ) :: VGLVS( MXLAYS3 + 1 ) ! vertical coord values
                    
C...........   EXTERNAL FUNCTIONS:
        LOGICAL      CHKINT
        LOGICAL      CHKREAL
C       CHARACTER(2) CRLF
C       INTEGER      GETEFILE
        INTEGER      GETNLIST
C       INTEGER      INDEX1
C       INTEGER      STR2INT
C       REAL         STR2REAL

C        EXTERNAL     CHKINT, CHKREAL, CRLF, GETEFILE, GETNLIST, 
C     &               INDEX1, STR2INT, STR2REAL
        EXTERNAL     CHKINT, CHKREAL, GETNLIST
        
C...........   Local parameters
        INTEGER, PARAMETER :: MXLAYTYP = 6

C...........   Layer structure names arrays
        CHARACTER(7):: LAYNAMES( MXLAYTYP ) = ( / 'VGSGPH3'
     &                                          , 'VGSGPN3'
     &                                          , 'VGSIGZ3'
     &                                          , 'VGPRES3'
     &                                          , 'VGZVAL3'
     &                                          , 'VGHVAL3' / )

C...........   Local arrays (note- tried to make these allocatable, but
C              this caused unexplainable crashing on SGI).
        CHARACTER(32) :: SEGMENT( 32 )
        CHARACTER(32) :: UPCSGMT( 32 )

C...........   File units and logical/physical names:
        INTEGER, SAVE :: IDEV    !  unit number of grid information file
        CHARACTER(16)    LNAME   !  logical name for grid information file

C...........   Scratch local variables and their descriptions:
            
        INTEGER         I, J, L, L2, LB, N  !  indices and string lengths
        INTEGER         IOS      !  I/O status return
        INTEGER         IREC     !  record counter
        INTEGER         NLAYREAD !  number of layers for reading sigma levels

        REAL(8)         DMISS3   !  double precision missing value

        LOGICAL      :: EFLAG = .FALSE.   ! true: error detected
        LOGICAL      :: FIRSTIME = .TRUE. ! true: first time routine called

        CHARACTER(300)  BUFFER   !  multi-purpose buffer
        CHARACTER(300)  LINE     !  line of file
        CHARACTER(300)  MESG     !  message buffer
        CHARACTER(300)  UPCLINE  !  upper case line of file

        CHARACTER(16) :: PROGNAME = 'DSCM3LAY' ! program name

C***********************************************************************
C   begin body of function DSCM3LAY

C.........  Get value of Models-3 environment variable for grid and check
C           exit status
        FIRSTIME = .FALSE.
        LNAME = 'G_GRIDPATH'

C.........  Open grid file
        IDEV = GETEFILE( LNAME, .TRUE., .TRUE., PROGNAME )

        IF( IDEV .LE. 0 ) THEN

            L = LEN_TRIM( LNAME )
            MESG = 'Could not open file "' // LNAME( 1:L ) // '"!'
            DSCM3LAY = .FALSE.
            RETURN

        END IF

C.........  Initialize grid setting to missing
        NLAYS   = IMISS3
        VGTYP   = IMISS3
        DMISS3  = DBLE( AMISS3 )

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
            END IF

C.............  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

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
     &                 'in G_GRIDPATH file skipped because ' //
     &                 CRLF() // BLANK10 // 'more than 32' //
     &                 'fields were found.'
                CALL M3MSG2( MESG )
                CYCLE 
            END IF
c note: add allocatable array
            SEGMENT = ' '            
            UPCSGMT = ' '            
            CALL PARSLINE( LINE   , N, SEGMENT )
            CALL PARSLINE( UPCLINE, N, UPCSGMT )

C.............  Search for keywords.  If the keyword exists, extract the
C               value(s) and cycle to the next line...

C.............  Grid name
            SELECT CASE ( UPCSGMT( 1 ) )
            CASE( 'NLAYS' )

C.................  Check to ensure the layer number is an integer
                IF ( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                     'ERROR: NLAYS value is not an integer '//
     &                     'in grid description ' // CRLF() // BLANK10//
     &                     'file at line', IREC
                    CALL M3MSG2( MESG )

C.................  For integer input...
                ELSE
                    NLAYS = STR2INT( SEGMENT( 3 ) )

                END IF

C.............  Grid description
            CASE( 'VG_TYP_GD' )

C.................  Check to ensure the layer number is an integer, and if
C                   not, see if layer type name is being given
                IF ( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN

                    VGTYP = INDEX1( SEGMENT( 3 ), MXLAYTYP, LAYNAMES )

                    IF ( VGTYP .LE. 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                     'ERROR: VG_TYP_GD value is not an integer '//
     &                     'or does not match the '// CRLF()// BLANK10//
     &                     'I/O API keywords in grid description ' // 
     &                     'file at line', IREC
                        CALL M3MSG2( MESG )
                    END IF

C.................  For integer input...
                ELSE
                    VGTYP = STR2INT( SEGMENT( 3 ) )

                END IF

C.............  Coordinate system type
            CASE( 'VGTOP_GD' )

C.................  Check to ensure the layer number is an floating point value
                IF ( .NOT. CHKREAL( SEGMENT( 3 ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                     'ERROR: VGTOP_GD value cannot be read '//
     &                     'from grid description'// CRLF()// BLANK10//
     &                     'file at line', IREC
                    CALL M3MSG2( MESG )

C.................  For real input...
                ELSE
                    VGTOP = STR2REAL( SEGMENT( 3 ) )
                END IF

C.............  Layer levels
            CASE( 'VGLVS_GD' )

C.................  Check to ensure the layer number is an integer
                IF ( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                     'ERROR: VGLVS_GD layer number is not '//
     &                     'and integer in grid '// CRLF()// BLANK10//
     &                     'description file ' //
     &                     'at line', IREC
                    CALL M3MSG2( MESG )

C.................  For integer input...
                ELSE

C.....................  Store the number of levels in the file
                    NLAYREAD = STR2INT( SEGMENT( 3 ) )

C.....................  Check to ensure that the number of layers available
C                       is consistent with the number of layers needed in file
                    IF ( NLAYREAD .GT. MXLAYS3 + 1 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                     'ERROR: Number of layers (', NLAYREAD,
     &                     ') exceeds maximum (', MXLAYS3+1, 
     &                     ') at line', IREC
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF

C.....................  Check to ensure that the number of layers available
C                       is consistent with the number specified elsewhere in the
C                       file
                    IF ( NLAYREAD-1 .NE. NLAYS ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                     'ERROR: Number of layers (', NLAYREAD-1,
     &                     ') differs from NLAYS in file (', NLAYS, 
     &                     ') at line', IREC
                        CALL M3MSG2( MESG )
                        CYCLE

                    END IF

C.....................  Read the layer structure
                    IF( LB .EQ. 0 ) 
     &                  READ( IDEV, *, ERR=999 ) 
     &                      ( VGLVS( I ), I=1, NLAYREAD )

                END IF

            END SELECT

        END DO

C.........  Exit from read loop
111     CONTINUE
 
        DSCM3LAY = ( .NOT. EFLAG )

C.........  Close grid file
        CLOSE( IDEV )

        RETURN

999     DSCM3LAY = .FALSE.
        MESG = 'ERROR: Could not read vertical layer structure ' //
     &         'from grid description file.'
        CALL M3MSG2( MESG )

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
            REAL(8)     , INTENT(OUT) :: ROUT     ! real output value
            INTEGER     , INTENT(OUT) :: IOUT     ! integer output value

C.............  Local variables
            INTEGER   I, L

            CHARACTER(300) MESG

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
            REAL(8)     , INTENT (IN) :: RVAL
            INTEGER     , INTENT (IN) :: IVAL
            INTEGER     , INTENT (IN) :: INTYPE

C.............  Local variables
            CHARACTER(300) MESG

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

        END FUNCTION DSCM3LAY

