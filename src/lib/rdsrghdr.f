
        SUBROUTINE RDSRGHDR( FDEV, FORMAT, GDNAME, PROJTYPE, PROJUNIT,
     &                       P_ALP, P_BET, P_GAM, XCENT, YCENT, XORIG,
     &                       YORIG, XCELL, YCELL, NCOLS, NROWS, NTHIK )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the header of the spatial surrogates file.
C     
C     
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
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

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER        STR2INT
        REAL           STR2REAL
        LOGICAL        CHKINT
	LOGICAL        CHKREAL

        EXTERNAL       STR2INT, STR2REAL, CHKINT, CHKREAL

C...........   Subroutine arguments

        INTEGER      , INTENT  (IN) :: FDEV       ! File unit number
        INTEGER      , INTENT  (OUT):: NCOLS      ! # cells in X direction
        INTEGER      , INTENT  (OUT):: NROWS      ! # cells in Y direction
        INTEGER      , INTENT  (OUT):: NTHIK      ! Boundary thickness
        CHARACTER(*) , INTENT  (OUT):: FORMAT     ! Format of surrogates file
        CHARACTER(*) , INTENT  (OUT):: GDNAME     ! Grid name
        CHARACTER(*) , INTENT  (OUT):: PROJTYPE   ! Project type
        CHARACTER(*) , INTENT  (OUT):: PROJUNIT   ! Projection units
        REAL         , INTENT  (OUT):: P_ALP      ! First map projection
                                                  ! description parameter
        REAL         , INTENT  (OUT):: P_BET      ! Second map projection
                                                  ! description parameter
        REAL         , INTENT  (OUT):: P_GAM      ! Third map projection
                                                  ! description parameter
        REAL         , INTENT  (OUT):: XCENT      ! Center of coordinate system
        REAL         , INTENT  (OUT):: YCENT      ! Center of coordinate system
        REAL         , INTENT  (OUT):: XORIG      ! X origin
        REAL         , INTENT  (OUT):: YORIG      ! Y origin
        REAL         , INTENT  (OUT):: XCELL      ! Cell size, X direction
        REAL         , INTENT  (OUT):: YCELL      ! Cell size, Y direction

C...........   Local parameters

        INTEGER, PARAMETER :: MXSEG = 16          ! # of potential line segments

C...........   Other arrays

        CHARACTER*20 SEGMENT( MXSEG )             ! Segments of parsed lines

C...........   Local variables

        INTEGER         I                         ! Model format flag
        INTEGER         IOS                       ! i/o status
        INTEGER         IREC                      ! Record counter
        LOGICAL         EFLAG                     ! Error flag
        CHARACTER*80    LINE                      ! Read buffer for a line
        CHARACTER*300   MESG                      ! Text for M3EXIT()
        CHARACTER*16 :: PROGNAME = 'RDSRGHDR'     ! Program name

C***********************************************************************
C   Begin body of subroutine RDSRGHDR

        IREC  = 0
        EFLAG = .FALSE.

        REWIND( FDEV )

        DO

        READ ( FDEV, 93000, END=12, IOSTAT=IOS ) LINE

        IREC = IREC + 1
             
        IF ( IOS .GT. 0 ) THEN
             WRITE( MESG, 94010)
     &            'I/O error', IOS, 'reading speciation profile '//
     &            'file at line', IREC
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL UPCASE( LINE )

        IF ( LINE .EQ. ' ' .AND. IREC .EQ. 1) THEN   ! skip if blank
             MESG = 'Line 1 of the spatial surrogates file did not '//
     &              'contain the header line'
             CALL M3WARN( PROGNAME, 0, 0, MESG )
             CYCLE
        ELSE
        END IF

        I = INDEX( LINE, '#GRID' )   ! determine if current line is the header

        IF ( I .GT. 0 ) THEN
             FORMAT = 'MODELS3'  ! set format to 'MODELS3'
             EXIT                ! exit read loop
        ELSE
        END IF 

        END DO

12      CONTINUE    ! end of read on input file

        SELECT CASE( FORMAT )

        CASE ( 'MODELS3' )

C.........  Parse the line of data into segments based on the rules
C           for "list-formatted" in fortran, but not requiring 
C           quotes around the text strings

        CALL PARSLINE( LINE, MXSEG, SEGMENT )

C.........  Make sure appropriate segments of header line appear to
C           be either INTEGER or REAL as expected 

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
     &       'Unexpected data type encountered in header of spatial '//
     &       'surrogates file'
             CALL M3MESG( MESG )

        ELSE

             GDNAME   =           SEGMENT ( 2 )
             XORIG    = STR2REAL( SEGMENT ( 3 ) )
             YORIG    = STR2REAL( SEGMENT ( 4 ) )
             XCELL    = STR2REAL( SEGMENT ( 5 ) )
             YCELL    = STR2REAL( SEGMENT ( 6 ) )
             NCOLS    = STR2INT ( SEGMENT ( 7 ) )
             NROWS    = STR2INT ( SEGMENT ( 8 ) )
             NTHIK    = STR2INT ( SEGMENT ( 9 ) )
             PROJTYPE =           SEGMENT ( 10 )
             PROJUNIT =           SEGMENT ( 11 )
             P_ALP    = STR2REAL( SEGMENT ( 12 ) )
             P_BET    = STR2REAL( SEGMENT ( 13 ) )
             P_GAM    = STR2REAL( SEGMENT ( 14 ) )
             XCENT    = STR2REAL( SEGMENT ( 15 ) )
             YCENT    = STR2REAL( SEGMENT ( 16 ) )

        END IF

        END SELECT

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END

        



