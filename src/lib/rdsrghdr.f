
        SUBROUTINE RDSRGHDR( FDEV, SRGFMT )

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
C
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C
C See file COPYRIGHT for conditions of use.
C
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures
        INCLUDE 'FLTERR.EXT'    !  error filter statement function

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2    CRLF
        LOGICAL        DSCM3GRD
        INTEGER        INDEX1
        INTEGER        STR2INT
        REAL           STR2REAL
        LOGICAL        CHKINT
	LOGICAL        CHKREAL

        EXTERNAL       CRLF, DSCM3GRD, INDEX1, STR2INT, STR2REAL, 
     &                 CHKINT, CHKREAL

C...........   Subroutine arguments

        INTEGER      , INTENT  (IN) :: FDEV       ! File unit number
        CHARACTER(*) , INTENT  (OUT):: SRGFMT     ! Format of surrogates file

C...........   Local parameters

        INTEGER, PARAMETER :: MXSEG = 16          ! # of potential line segments
        INTEGER, PARAMETER :: MXGRDTYP = 11

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
     &                                           , UTMGRD3
     &                                           , UTMGRD3 / )

        CHARACTER*15 :: GRDNAMES( MXGRDTYP ) = ( / 'LAT-LON        '
     &                                           , 'GEOGRAPHIC     '
     &                                           , 'LATGRD3        '
     &                                           , 'LAMBERT        '
     &                                           , 'LAMGRD3        '
     &                                           , 'MERCATOR       '
     &                                           , 'MERGRD3        '
     &                                           , 'STEREOGRAPHIC  '
     &                                           , 'STEGRD3        '
     &                                           , 'UTM            '
     &                                           , 'UTMGRD3        ' / )

C...........   Other arrays

        CHARACTER*20 SEGMENT( MXSEG )             ! Segments of parsed lines

C...........   Local variables

        INTEGER       I, J, L    ! counters and indices
        INTEGER       IOS        ! i/o status
        INTEGER       IREC       ! record counter
        INTEGER       NTHIK      ! boundary thickness

        LOGICAL    :: EFLAG = .FALSE.           ! Error flag

        CHARACTER*7            PROJUNIT   ! Projection units
        CHARACTER*80           MSGEND     ! tmp end of message
        CHARACTER*200          LINE       ! Read buffer for a line
        CHARACTER*300          MESG       ! Message buffer
        CHARACTER(LEN=IOVLEN3) COORD3D    ! coordinate system name 
        CHARACTER(LEN=IOVLEN3) COORUN3D   ! coordinate system projection units
        CHARACTER(LEN=IOVLEN3) PROJTYPE   ! coordinate system projection name

        CHARACTER*16 :: PROGNAME = 'RDSRGHDR'     ! Program name

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
     &            'I/O error', IOS, 'reading gridding surrogates '//
     &            'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL UPCASE( LINE )

            IF ( LINE .EQ. ' ' ) CYCLE  ! skip all blank lines
  
            I = INDEX( LINE, '#GRID' ) ! determine if current line is the header
 
            IF ( I .GT. 0 ) THEN
                SRGFMT = 'MODELS3'  ! set format to 'MODELS3'
                IF( IREC .NE. 1 ) THEN   
                    MESG = 'Line 1 of the spatial surrogates file '//
     &                     'did not contain the header line.'
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
     &            'spatial surrogates file'
                CALL M3MESG( MESG )

            ELSE

                GDNAM3D =            SEGMENT ( 2 )
                XORIG3D  = STR2REAL( SEGMENT ( 3 ) )
                YORIG3D  = STR2REAL( SEGMENT ( 4 ) )
                XCELL3D  = STR2REAL( SEGMENT ( 5 ) )
                YCELL3D  = STR2REAL( SEGMENT ( 6 ) )
                NCOLS3D  = STR2INT ( SEGMENT ( 7 ) )
                NROWS3D  = STR2INT ( SEGMENT ( 8 ) )
                NTHIK3D  = STR2INT ( SEGMENT ( 9 ) )
                PROJTYPE =           SEGMENT ( 10 )
                PROJUNIT =           SEGMENT ( 11 )
                P_ALP3D  = STR2REAL( SEGMENT ( 12 ) )
                P_BET3D  = STR2REAL( SEGMENT ( 13 ) )
                P_GAM3D  = STR2REAL( SEGMENT ( 14 ) )
                XCENT3D  = STR2REAL( SEGMENT ( 15 ) )
                YCENT3D  = STR2REAL( SEGMENT ( 16 ) )

            END IF

        END SELECT

C.........  Set project code based on projection type
        J =  INDEX1( PROJTYPE, MXGRDTYP, GRDNAMES )
        IF ( J .LE. 0 ) THEN
            EFLAG = .TRUE.
            L = LEN_TRIM( PROJTYPE )
            MESG = 'Projection type "' // PROJTYPE( 1:L ) //
     &             '" is not recognized by the rdsrghdr.f routine.'
            CALL M3MSG2( MESG )

C.........  Otherwise, set the grid type code number
C.........  Initialize grid information based on the surrogates file
        ELSE

            GDTYP3D = GRDTYPES( J )
            CALL CHKGRID( GDNAM3D, 'SURROGATES', 0, EFLAG )

        END IF

C.........  Abort if an error is found
        IF ( EFLAG ) THEN

            MESG = 'Problem with surrogate file'
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

        



