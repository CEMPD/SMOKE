
        SUBROUTINE RDSRG( FDEV, FORMAT, XCENT, YCENT, XORIG, YORIG,
     &                    XCELL, YCELL, NCOLS, NROWS ) 

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for the spatial surrogate
C      arrays, reads the the spatial surrogates file, and then fills
C      the spatial surrogate arrays.
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

C...........   Modules for public variables

        USE MODSURG

        IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER        STR2INT
        REAL           STR2REAL

        EXTERNAL       STR2INT, STR2REAL

C...........   Subroutine arguments

        INTEGER      , INTENT  (IN) :: FDEV       ! File unit number
        CHARACTER(*) , INTENT  (IN) :: FORMAT     ! Format of surrogates file
        REAL         , INTENT  (IN) :: XCENT      ! Center of coordinate system
        REAL         , INTENT  (IN) :: YCENT      ! Center of coordinate system
        REAL         , INTENT  (IN) :: XORIG      ! X origin
        REAL         , INTENT  (IN) :: YORIG      ! Y origin
        REAL         , INTENT  (IN) :: XCELL      ! Cell size, X direction
        REAL         , INTENT  (IN) :: YCELL      ! Cell size, Y direction
        INTEGER      , INTENT  (IN) :: NCOLS      ! # cells in X direction
        INTEGER      , INTENT  (IN) :: NROWS      ! # cells in Y direction

C...........   Local parameters

        INTEGER, PARAMETER :: MXSEG = 5           ! # of potential line segments

C...........   Other arrays

        CHARACTER*20 SEGMENT( MXSEG )             ! Segments of parsed lines

C...........   Local variables

        INTEGER         COL                   ! Temp grid column number (x-axis)
        INTEGER         COUNTY                ! Temp Country/State/County code
        INTEGER         I                     ! Header flag
        INTEGER         IOS                   ! i/o status
        INTEGER         IREC                  ! Record counter
        INTEGER         J                     ! Counter for surrogate entries
        INTEGER         ROW                   ! Temp grid row number (y-axis)
        INTEGER         SSC                   ! Temp spatial surrogate code
        REAL            RATIO                 ! Temp spatial surrogate Ratio
        CHARACTER*80    LINE                  ! Read buffer for a line
        CHARACTER*300   MESG                  ! Text for M3EXIT()
        CHARACTER*16 :: PROGNAME = 'RDSRG'    !  program name

C***********************************************************************
C   Begin body of subroutine RDSRG

        IREC    = 0
	J       = 0
        NSRGREC = 0

C......... Determine the number surrogate file entries

        REWIND( FDEV )

        SELECT CASE( FORMAT )
	
        CASE( 'MODELS3' )

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

        I = INDEX( LINE, '#GRID' )    ! determine if current line is the header

        IF ( I .GT. 0 ) THEN          ! skip if current line is header
             CYCLE
        ELSEIF ( LINE .EQ. ' ' ) THEN ! skip if current line is blank
             CYCLE
        ELSE
             NSRGREC = NSRGREC + 1
        END IF

        END DO

12      CONTINUE    ! end of read on input file

        END SELECT

C......... Allocate memory for surrogate arrays

        ALLOCATE( IDXSRGA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRGA', PROGNAME )

        ALLOCATE( IDXSRGB( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRGB', PROGNAME )

        ALLOCATE( SCELLA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SCELLA', PROGNAME )

        ALLOCATE( SFIPSA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SFIPSA', PROGNAME )

        ALLOCATE( SSRGIDA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SSRGIDA', PROGNAME )

        ALLOCATE( SFRACA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SFRACA', PROGNAME )

C......... Fill surrogate arrays

        REWIND( FDEV )
	IREC    = 0

        SELECT CASE( FORMAT )
	
        CASE( 'MODELS3' )

        DO

        READ ( FDEV, 93000, END=24, IOSTAT=IOS ) LINE

        IREC = IREC + 1
             
        IF ( IOS .GT. 0 ) THEN
             WRITE( MESG, 94010)
     &            'I/O error', IOS, 'reading speciation profile '//
     &            'file at line', IREC
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL UPCASE( LINE )

        I = INDEX( LINE, '#GRID' )    ! determine if current line is the header

        IF ( I .GT. 0 ) THEN          ! skip if current line is header
             CYCLE
        ELSEIF ( LINE .EQ. ' ' ) THEN ! skip if current line is blank
             CYCLE
        ELSE
             J = J + 1 
        END IF

C.........  Parse the line of data into segments based on the rules
C           for "list-formatted" in fortran, but not requiring 
C           quotes around the text strings

        CALL PARSLINE( LINE, 5, SEGMENT )

        IDXSRGA( J ) = J
        IDXSRGB( J ) = J
        SSC    = STR2INT ( SEGMENT( 1 ) )
        COUNTY = STR2INT ( SEGMENT( 2 ) )
        COL    = STR2INT ( SEGMENT( 3 ) )
        ROW    = STR2INT ( SEGMENT( 4 ) )
        RATIO  = STR2REAL( SEGMENT( 5 ) )

        SCELLA ( J ) = (ROW-1)*NCOLS + COL
        SFIPSA ( J ) = COUNTY
        SSRGIDA( J ) = SSC
        SFRACA ( J ) = RATIO     


        END DO

24      CONTINUE    ! end of read on input file

        END SELECT

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
