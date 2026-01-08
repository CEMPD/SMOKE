
        SUBROUTINE GENPROJ( PDEV, PYEAR, ENAME )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine processes the projection data and writes out
C      the projection matrix.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     
C
C************************************************************************
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
C*************************************************************************

C.........  MODULES for public variables
C.........  This module contains the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: CSOURC, CISIC, CMACT

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: POLSFLAG, NVPROJ, FACTOR, RPTDEV,
     &                      PNAMPROJ, PRJFC

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: MXCHRS, BYEAR, CATDESC, NCHARS, NSRC,
     &                     SC_BEGP, SC_ENDP, CRL, CATEGORY

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  i/o api parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS and their descriptions:
C       INTEGER         GETEFILE
C       INTEGER         PROMPTFFILE

C        EXTERNAL   GETEFILE, PROMPTFFILE

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN OUT) :: PDEV   ! temporary file
        INTEGER     , INTENT (IN) :: PYEAR  ! projection year
        CHARACTER(*), INTENT (IN) :: ENAME  ! emission inventory file name

C...........  Local static arrays
        LOGICAL          LF   ( MXCHRS )      !  true: column should be output
        CHARACTER(20)    CHARS( MXCHRS )      !  source fields for output

C...........   Local variables that depend on module variables
        CHARACTER(10) CHRHDRS( NCHARS )  ! Source characteristics headers

C...........   Logical names and unit numbers
c        INTEGER          ODEV       ! unit number of output tmp file
        CHARACTER(16)    PNAME      ! logical name for projection matrix

C...........   Other local variables

        INTEGER          J, K, L, N, S, V    ! counters and indices
        INTEGER          IDUM          ! dummy integer
        INTEGER          IOS           ! i/o error status
        INTEGER          NC            ! local number src chars
        INTEGER          RDEV          ! Report unit number

        LOGICAL       :: EFLAG    = .FALSE.  ! true: error has occurred
        LOGICAL, SAVE :: FIRSTOUT = .TRUE.   ! true: output header line
        LOGICAL, SAVE :: APPLFLAG = .FALSE.  ! true: something has been applied
        LOGICAL       :: MACTFLAG = .FALSE.  ! true: MACT field is allocated

        CHARACTER(SICLEN3) :: CSIC    ! SIC code
        CHARACTER(200)     :: PATHNM  ! path name for tmp file
        CHARACTER(220)        FILENM  ! file name
        CHARACTER(512)        MESG    ! message buffer
        CHARACTER(512)        HDRBUF  ! report header buffer
        CHARACTER(IOVLEN3) :: PNAM    ! tmp pol/act name

        CHARACTER(16) :: PROGNAME = 'GENPROJ' ! program name

C***********************************************************************
C   begin body of subroutine GENPROJ

C.........  Get path for temporary files
        MESG = 'Path where temporary control files will be written'
        CALL ENVSTR( 'SMK_TMPDIR', MESG, '.', PATHNM, IOS )

C.........  Open reports file
        RPTDEV( 4 ) = PROMPTFFILE( 
     &                'Enter logical name for PROJECTION REPORT', 
     &                .FALSE., .TRUE., CRL // 'PROJREP', PROGNAME )
        RDEV = RPTDEV( 4 )

C.........  Open *output* temporary file
c NOTE: This is commented out because output tmp file is not used in wcntlrep.f
c        FILENM = TRIM( PATHNM ) // '/cntlmat_tmp_proj_rep'
c        ODEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

C.........  Set up and open output projection matrices
        CALL OPENPMAT( ENAME, BYEAR, PYEAR, PNAME )

C..........  Write header for report.
        WRITE( RDEV, 93000 ) '# Processed as '// CATDESC// ' sources'
        WRITE( RDEV, 93000 ) 
     &         '# Projection factors applied with /PROJECTION/ packet'

        WRITE( RDEV, 93390 ) '#      from base year    ', BYEAR
        WRITE( RDEV, 93390 ) '#      to projected year ', PYEAR

        IF( POLSFLAG ) THEN
            WRITE( RDEV, 93000 ) '#      using pollutant-specific ' //
     &                           'assignments'
        ELSE
            WRITE( RDEV, 93000 ) '#      to all pollutants uniformly'
        END IF

C.........  Define source-category specific header
C.........  NOTE that (1) will not be used and none will be for area sources
        CHRHDRS( 1 ) = 'Region'
        SELECT CASE( CATEGORY )
        CASE( 'AREA' )
            CHRHDRS( 2 ) = 'SCC'

        CASE( 'MOBILE' )
            CHRHDRS( 2 ) = 'Road '
            CHRHDRS( 3 ) = 'Link'
            CHRHDRS( 4 ) = 'Veh Type'
            CHRHDRS( 5 ) = 'SCC'

        CASE( 'POINT' )
            CHRHDRS( 2 ) = 'Plant ID'
            IF ( NCHARS .GE. 3 ) CHRHDRS( 3 ) = 'Char 1'
            IF ( NCHARS .GE. 4 ) CHRHDRS( 4 ) = 'Char 2'
            IF ( NCHARS .GE. 5 ) CHRHDRS( 5 ) = 'Char 3'
            IF ( NCHARS .GE. 6 ) CHRHDRS( 6 ) = 'Char 4'
            IF ( NCHARS .GE. 7 ) CHRHDRS( 7 ) = 'Char 5'
            CHRHDRS( NCHARS ) = 'SCC'   

        END SELECT

C.........  Check if MACT field is available from inventory
        MACTFLAG = ( ASSOCIATED( CMACT ) )

C.........  Loop through all sources and store projection information for
C           those that have it.  Otherwise, set projection factor=1.

C.........  Initialize valid columns
        LF = .FALSE.  ! array
        DO J = 1, NCHARS
            LF( J ) = .TRUE.
        END DO

C.........  If no pollutant-specific assignments, write out single pfac var
        IF ( .NOT. POLSFLAG ) NVPROJ = 1

C.........  Loop through pollutants that are getting projections
        DO V = 1, NVPROJ

            IF( POLSFLAG ) THEN
                PNAM = PNAMPROJ( V )
            ELSE
                PNAM = 'pfac'
            END IF

C...........  Loop through sources, retrieve projection packet index,
C             set projection factor, and write out new tmp file info.
            DO S = 1, NSRC

                READ( PDEV, * ) K

C................  If source has projection info...
                IF ( K .GT. 0 ) THEN

C...................  Store projection factor for current pollutant/act
                    FACTOR( S ) = PRJFC( K )

c                    WRITE( ODEV,93300 ) 1, PNAM, FACTOR( S )
                    APPLFLAG = .TRUE.

C....................  Write report
C....................  Format source characteristic information
                    CALL PARSCSRC( CSOURC(S), MXCHRS, SC_BEGP, SC_ENDP, 
     &                             LF, NC, CHARS )
                    NC = MIN( NC, NCHARS )

C.................  Write out projection information for all sources
C                   that are getting projected
                    IF( ASSOCIATED( CISIC ) ) THEN
                        CSIC = CISIC( S )
                    ELSE
                        CSIC = ''
                    END IF

                    IF( POLSFLAG ) THEN

C.........................  Write header and use same header names as Svmkreport
C                           to facilitate EMF importing
                        WRITE( HDRBUF, 94015 ) '#       Variable', 
     &                        (CHRHDRS(N), N=1, NC), 'SIC'
                        WRITE( MESG, 94015 ) PNAM, 
     &                   ( CHARS(J)(1:SC_ENDP(J)-SC_BEGP(J)+1),J=1,NC ),
     &                   CSIC
                    ELSE
                        WRITE( HDRBUF, 94016 ) 
     &                       ( CHRHDRS(N), N=1, NC ), 'SIC'
                        WRITE( MESG, 94016 ) 
     &                   ( CHARS(J)(1:SC_ENDP(J)-SC_BEGP(J)+1),J=1,NC ),
     &                   CSIC
                    END IF

C.....................  Add more to header and output if MACT is present
                    IF( MACTFLAG ) THEN
                        L = LEN_TRIM( HDRBUF )
                        HDRBUF = HDRBUF( 1:L ) // '; MACT'

                        L = LEN_TRIM( MESG )
                        MESG = MESG( 1:L ) // '; ' // CMACT( S )
                    ENDIF

C.....................  Complete header
                    L = LEN_TRIM( HDRBUF )
                    HDRBUF = HDRBUF( 1:L ) // '; Data value'
                    
                    IF( FIRSTOUT ) THEN
                        WRITE( RDEV, '(A)' ) TRIM( HDRBUF )
                        FIRSTOUT = .FALSE.
                    ENDIF

                    WRITE( RDEV, 94020 ) TRIM( MESG ), FACTOR( S )

C...............  If source does not have projection info., set to 1.
                ELSE

                    FACTOR( S ) = 1.0

C..................  Add to tmp file line
c                    WRITE( ODEV, 93300 ) 0, 'D', 1.0

                END IF

            END DO      ! End loop through sources

C.............  Write projection factors
            IF( .NOT. WRITESET( PNAME,PNAM,ALLFILES,0,0,FACTOR )) THEN

                MESG = 'Problem writing "'// TRIM( PNAM )// '" to '//
     &                 'output file "' // TRIM( PNAME )// '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

        END DO          ! End loop through pol/act

        IF( .NOT. APPLFLAG ) THEN

            MESG = 'WARNING: No PROJECTION packet entries match ' //
     &             'inventory.'
            CALL M3MSG2( MESG )

            MESG = 'WARNING: Projection matrix will not be created!'
            CALL M3MSG2( MESG )

C.............  Write not into report file
            WRITE( RDEV, 93000 ) 
     &             '# No projection packet entries matched the inventory.'

            RETURN

        END IF

C.........  Reset tmp file to be file just output
c        PDEV = ODEV

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx
94015   FORMAT( A16, ';', 1X, 10( A, :,';', 1X ) )

94016   FORMAT( 10( A, :, ';', 1X ) )

94020   FORMAT( A, ';', 1X, E13.5 )

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

c93300   FORMAT( I2, 1X, '"', A, '"', 3( 1X, E12.5 ) )

93390   FORMAT( A, I4.4 )

C******************  INTERNAL SUBPROGRAMS  *****************************

        END SUBROUTINE GENPROJ
