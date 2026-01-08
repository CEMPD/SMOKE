
        SUBROUTINE RDSRGHDR( VFLAG, FDEV, SRGFMT )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the header of the spatial surrogates file.
C     
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Added POLGRD3 as a supported coord sys type - E. Giroux CNRC 03/2004
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures

C...........   EXTERNAL FUNCTIONS and their descriptions:

C       CHARACTER(2)   CRLF
        LOGICAL        DSCM3GRD
C       INTEGER        INDEX1
C       INTEGER        STR2INT
C       REAL           STR2REAL
C       REAL*8         STR2DBLE
        LOGICAL        CHKINT
        LOGICAL        CHKREAL

C        EXTERNAL       CRLF, DSCM3GRD, INDEX1, STR2INT, STR2REAL, 
C     &                 CHKINT, CHKREAL, STR2DBLE
        EXTERNAL     DSCM3GRD, CHKINT, CHKREAL

C...........   Subroutine arguments
        LOGICAL      , INTENT  (IN) :: VFLAG      ! true: using variable grid
        INTEGER      , INTENT  (IN) :: FDEV       ! File unit number
        CHARACTER(*) , INTENT  (OUT):: SRGFMT     ! Format of surrogates file

C...........   Local parameters

        INTEGER, PARAMETER :: MXSEG = 16          ! # of potential line segments
        INTEGER, PARAMETER :: MXGRDTYP = 13

C...........   Grid types and names arrays
        INTEGER      :: GRDTYPES( MXGRDTYP ) = ( / LATGRD3
     &                                           , LATGRD3
     &                                           , LATGRD3
     &                                           , LAMGRD3
     &                                           , LAMGRD3
     &                                           , MERGRD3
     &                                           , MERGRD3
     &                                           , STEGRD3
     &                                           , STEGRD3
     &                                           , POLGRD3
     &                                           , POLGRD3
     &                                           , UTMGRD3
     &                                           , UTMGRD3 / )

        CHARACTER(15) :: GRDNAMES( MXGRDTYP ) = ( / 'LAT-LON        '
     &                                           , 'GEOGRAPHIC     '
     &                                           , 'LATGRD3        '
     &                                           , 'LAMBERT        '
     &                                           , 'LAMGRD3        '
     &                                           , 'MERCATOR       '
     &                                           , 'MERGRD3        '
     &                                           , 'STEREOGRAPHIC  '
     &                                           , 'STEGRD3        '
     &                                           , 'POLAR          '
     &                                           , 'POLGRD3        '
     &                                           , 'UTM            '
     &                                           , 'UTMGRD3        ' / )

C...........   Other arrays

        CHARACTER(20) SEGMENT( MXSEG )             ! Segments of parsed lines

C...........   Local variables

        INTEGER       I, J, L    ! counters and indices
        INTEGER       IOS        ! i/o status
        INTEGER       IREC       ! record counter
        INTEGER       NTHIK      ! boundary thickness

        LOGICAL    :: EFLAG = .FALSE.           ! Error flag

        CHARACTER(7)       PROJUNIT   ! Projection units
        CHARACTER(80)      MSGEND     ! tmp end of message
        CHARACTER(320)     LINE       ! Read buffer for a line
        CHARACTER(300)     MESG       ! Message buffer
        CHARACTER(IOVLEN3) COORD3D    ! coordinate system name 
        CHARACTER(IOVLEN3) COORUN3D   ! coordinate system projection units
        CHARACTER(IOVLEN3) PROJTYPE   ! coordinate system projection name

        CHARACTER(16) :: PROGNAME = 'RDSRGHDR'     ! Program name

C***********************************************************************
C   Begin body of subroutine RDSRGHDR

        IREC  = 0
        EFLAG = .FALSE.

        REWIND( FDEV )

C.........  Loop through lines of file until the header line is encountered
        DO

            READ ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

            IREC = IREC + 1
             
            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &            'I/O error', IOS, 'reading gridding information '//
     &            'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL UPCASE( LINE )

            IF ( LINE .EQ. ' ' ) CYCLE  ! skip all blank lines
  
C.............  Determine if current line is the header
            IF( VFLAG ) THEN
                I = INDEX( LINE, '#VARIABLE_GRID' )
            ELSE
                I = INDEX( LINE, '#GRID' )
            END IF
 
            IF ( I .GT. 0 ) THEN
                SRGFMT = 'MODELS3'  ! set format to 'MODELS3'
                IF( IREC .NE. 1 ) THEN   
                    MESG = 'First line of the file did not '//
     &                     'contain the header line.'
                    CALL M3WARN( PROGNAME, 0, 0, MESG )
                END IF

                EXIT                ! exit read loop with header found

            END IF 

        END DO

        SELECT CASE( SRGFMT )

        CASE ( 'MODELS3' )

C.............  Parse the line of data into segments based on the rules
C               for "list-formatted" in fortran, but not requiring 
C               quotes around the text strings

            CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.............  Make sure appropriate segments of header line appear to
C               be either INTEGER or REAL as expected 

            IF ( .NOT. CHKREAL( SEGMENT( 3 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 4 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 5 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 6 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKINT ( SEGMENT( 7 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKINT ( SEGMENT( 8 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKINT ( SEGMENT( 9 ) ) )  EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 12 ) ) ) EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 13 ) ) ) EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 14 ) ) ) EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 15 ) ) ) EFLAG = .TRUE.
            IF ( .NOT. CHKREAL( SEGMENT( 16 ) ) ) EFLAG = .TRUE.

            IF ( EFLAG ) THEN

                WRITE( MESG, 94010 )
     &            'Unexpected data type encountered in header of '//
     &            'the file'
                CALL M3MESG( MESG )

            ELSE

                GDNAM3D =            SEGMENT ( 2 )
                XORIG3D  = STR2DBLE( SEGMENT ( 3 ) )
                YORIG3D  = STR2DBLE( SEGMENT ( 4 ) )
                XCELL3D  = STR2DBLE( SEGMENT ( 5 ) )
                YCELL3D  = STR2DBLE( SEGMENT ( 6 ) )
                NCOLS3D  = STR2INT ( SEGMENT ( 7 ) )
                NROWS3D  = STR2INT ( SEGMENT ( 8 ) )
                NTHIK3D  = STR2INT ( SEGMENT ( 9 ) )
                PROJTYPE =           SEGMENT ( 10 )
                PROJUNIT =           SEGMENT ( 11 )
                P_ALP3D  = STR2DBLE( SEGMENT ( 12 ) )
                P_BET3D  = STR2DBLE( SEGMENT ( 13 ) )
                P_GAM3D  = STR2DBLE( SEGMENT ( 14 ) )
                XCENT3D  = STR2DBLE( SEGMENT ( 15 ) )
                YCENT3D  = STR2DBLE( SEGMENT ( 16 ) )

            END IF

        END SELECT

C.........  Set project code based on projection type
        J =  INDEX1( PROJTYPE, MXGRDTYP, GRDNAMES )
        IF ( J .LE. 0 ) THEN
            EFLAG = .TRUE.
            L = LEN_TRIM( PROJTYPE )
            MESG = 'Projection type "' // PROJTYPE( 1:L ) //
     &             '" is not recognized.'
            CALL M3MSG2( MESG )

C.........  Otherwise, set the grid type code number
C.........  Initialize grid information based on the surrogates file
        ELSE

            GDTYP3D = GRDTYPES( J )
            CALL CHKGRID( GDNAM3D, 'SURROGATES', 0, EFLAG )

        END IF

C.........  Abort if an error is found
        IF ( EFLAG ) THEN

            MESG = 'Problem with processing grid information'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        RETURN

C.........  If get here, then the header was never found!
999     MESG = 'Surrogates file is missing a header line'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDSRGHDR

        



