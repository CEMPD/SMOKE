
        SUBROUTINE RPNTSCHR( INFILE, FDEV, NPSRC, NVARS, VNAMES, NCHARS)

C***********************************************************************
C  program body starts at line 
CC
C  DESCRIPTION:
C     Reads in SMOKE point inventory characteristics
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

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2            CRLF
        INTEGER                INDEX1

        EXTERNAL               CRLF, INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: INFILE          ! inven logical file name
        INTEGER     , INTENT (IN) :: FDEV            ! inven ASCII file unit no.
        INTEGER     , INTENT (IN) :: NPSRC           ! no. of point sources
        INTEGER     , INTENT (IN) :: NVARS           ! no. of inven vars to read
        CHARACTER(*), INTENT (IN) :: VNAMES( NVARS ) ! variable names 
        INTEGER     , INTENT(OUT) :: NCHARS          ! no. src characterstics

C...........   Index for unread variables
        INTEGER        UNREAD
        INTEGER        UNRIDX( NVARS )

C...........   Error message strings
        CHARACTER*23 ,PARAMETER :: PART1 = 'Error reading variable '
        CHARACTER*24 ,PARAMETER :: PART3 = ' from POINT SOURCES file'

C...........   Other local variables
        INTEGER          J, N, S    ! counters and indices

        INTEGER          ID         ! tmp smoke ID
        INTEGER          IOS        ! i/o status
        INTEGER          NCOL       ! number of columns in FDEV

        LOGICAL       :: BLRIN   = .FALSE.  ! True: Boiler is in input file
        LOGICAL       :: CSRFLAG = .FALSE.  ! True: Source chars requested
        LOGICAL       :: EFLAG   = .FALSE.  ! True: error
        LOGICAL       :: BLRFLAG = .FALSE.  ! True: Boilers requested
        LOGICAL       :: PDSFLAG = .FALSE.  ! True: Plant description requested
        LOGICAL       :: SCCFLAG = .FALSE.  ! True: SCC requested

        CHARACTER*20           HEADER( 20 ) ! header fields
        CHARACTER*300          FILFMT ! ASCII file format after header
        CHARACTER*300          LINE   ! line of ASCII point source file
        CHARACTER*300          MESG   ! message buffer
        CHARACTER(LEN=BLRLEN3) CBLR   ! temporary boiler name
        CHARACTER(LEN=FIPLEN3) CFIP   ! temporary character FIPs code
        CHARACTER(LEN=CHRLEN3) CHARS( 5 ) ! temporary character FIPs code
        CHARACTER(LEN=DSCLEN3) CPDS   ! temporary plant description
        CHARACTER(LEN=SCCLEN3) CS     ! temporary scc
        CHARACTER(LEN=PLTLEN3) FCID   ! temporary facility code
        CHARACTER(LEN=IOVLEN3) INVAR  ! tmp inventory pollutant name

        CHARACTER*16 :: PROGNAME = 'RPNTSCHR'   !  program name

C***********************************************************************
C   begin body of program RPNTSCHR

C.........  Loop through input variables
C.........  Allocate memory and read the ones that are needed from I/O API file

        UNREAD = 0
        DO N = 1, NVARS

            INVAR = VNAMES( N )

            MESG = PART1 // INVAR( 1:LEN_TRIM( INVAR ) ) // PART3

            SELECT CASE( INVAR )

            CASE( 'IFIP' )
                ALLOCATE( IFIP( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IFIP', PROGNAME )

                IF( .NOT. READ3(INFILE,'IFIP',ALLAYS3,0,0,IFIP ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'ISIC' )
                ALLOCATE( ISIC( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'ISIC', PROGNAME )

                IF( .NOT. READ3(INFILE,'ISIC',ALLAYS3,0,0,ISIC ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'IORIS' )
                ALLOCATE( IORIS( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IORIS', PROGNAME )

                IF( .NOT. READ3(INFILE,'IORIS',ALLAYS3,0,0,IORIS)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'TZONES' )
                ALLOCATE( TZONES( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'TZONES', PROGNAME )

                IF(.NOT. READ3(INFILE,'TZONES',ALLAYS3,0,0,TZONES)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'TPFLAG' )
                ALLOCATE( TPFLAG( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'TPFLAG', PROGNAME )

                IF(.NOT. READ3(INFILE,'TPFLAG',ALLAYS3,0,0,TPFLAG)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'INVYR' )
                ALLOCATE( INVYR( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVYR', PROGNAME )

                IF( .NOT. READ3(INFILE,'INVYR',ALLAYS3,0,0,INVYR)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
           
            CASE( 'IDIU' )
                ALLOCATE( IDIU( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IDIU', PROGNAME )

                IF( .NOT. READ3(INFILE,'IDIU',ALLAYS3,0,0,IDIU ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'IWEK' )
                ALLOCATE( IWEK( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IWEK', PROGNAME )

                IF( .NOT. READ3(INFILE,'IWEK',ALLAYS3,0,0,IWEK ) ) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'XLOCA' )
                ALLOCATE( XLOCA( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )

                IF( .NOT. READ3(INFILE,'XLOCA',ALLAYS3,0,0,XLOCA)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'YLOCA' )
                ALLOCATE( YLOCA( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )

                IF( .NOT. READ3(INFILE,'YLOCA',ALLAYS3,0,0,YLOCA)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'STKHT' )
                ALLOCATE( STKHT( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'STKHT', PROGNAME )

                IF( .NOT. READ3(INFILE,'STKHT',ALLAYS3,0,0,STKHT)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'STKDM' )
                ALLOCATE( STKDM( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'STKDM', PROGNAME )

                IF( .NOT. READ3(INFILE,'STKDM',ALLAYS3,0,0,STKDM)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'STKTK' )
                ALLOCATE( STKTK( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'STKTK', PROGNAME )

                IF( .NOT. READ3(INFILE,'STKTK',ALLAYS3,0,0,STKTK)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

            CASE( 'STKVE' )
                ALLOCATE( STKVE( NPSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'STKVE', PROGNAME )

                IF( .NOT. READ3(INFILE,'STKVE',ALLAYS3,0,0,STKVE)) THEN
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF
            
            CASE DEFAULT
                UNREAD = UNREAD + 1
                UNRIDX( UNREAD ) = N

            END SELECT

        ENDDO  ! End loop on variables indicated for reading
            
C.........  Section if the ASCII file is also to be read in
        IF( FDEV .LE. 0 .AND. UNREAD .GT. 0 ) THEN
            EFLAG = .TRUE.

            DO J = 1, UNREAD

                N = UNRIDX( J )
                INVAR = VNAMES( N )

                MESG = 'INTERNAL ERROR: Program "' // 
     &                 PROGNAME( 1:LEN_TRIM( PROGNAME ) ) //
     &                 '" found unread variable "' //
     &                 INVAR( 1:LEN_TRIM( INVAR ) ) //
     &                 '." NOTE: ASCII PNTS file is not opened.'
                CALL M3MSG2( MESG )

            ENDDO

        ELSEIF( FDEV .GT. 0 ) THEN

C.............  Allocate memory for the data that are needed from the ASCII file

            DO J = 1, UNREAD

                N = UNRIDX( J )

                INVAR = VNAMES( N )

                SELECT CASE( INVAR )

                CASE( 'CSCC' )
                    SCCFLAG = .TRUE. 
                    ALLOCATE( CSCC( NPSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CSCC', PROGNAME )

                CASE( 'CBLRID' )
                    BLRFLAG = .TRUE. 
                    ALLOCATE( CBLRID( NPSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CBLRID', PROGNAME )

                CASE( 'CPDESC' )
                    PDSFLAG = .TRUE. 
                    ALLOCATE( CPDESC( NPSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CPDESC', PROGNAME )

                CASE( 'CSOURC' )
                    CSRFLAG = .TRUE. 
                    ALLOCATE( CSOURC( NPSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CSOURC', PROGNAME )

                CASE DEFAULT
                    EFLAG = .TRUE.
                    MESG = 'INTERNAL ERROR: Program "' // 
     &                     PROGNAME( 1:LEN_TRIM( PROGNAME ) ) //
     &                     '" does not know about variable "' // 
     &                     INVAR( 1:LEN_TRIM( INVAR ) ) // '"'
                    CALL M3MSG2( MESG )

                END SELECT

            ENDDO  ! End loop of unread variables

C.............  Read in and store data from ASCII file...

C.............  Read in number of header lines
            READ( FDEV, * ) NCOL, FILFMT

C.............  Read past header
            DO J = 1, NCOL
                READ( FDEV, '(A)' ) HEADER( J )
            ENDDO

C.............  Determine number of source characteristics. Assumes SCC is
C               the first variable after the final source characteristic 
            J = INDEX1( 'SCC', NCOL, HEADER )
            NCHARS = J - 4  ! Four because Source ID, FIPS, Plant ID, & SCC 

C.............  Determine if boiler is present
            J = INDEX1( 'Boiler code', NCOL, HEADER )
            BLRIN = ( J .GT. 0 )

C.............  If boiler not present but has been requested, then internal err
            IF( .NOT. BLRIN .AND. BLRFLAG ) THEN

                MESG = 'INTERNAL ERROR: Boiler request, but is ' //
     &                 'not present in ASCII inventory file'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            ENDIF

C.............  Initialize arrays of source characteristics
            CHARS = ' '  ! array

            DO S = 1, NPSRC

C.................  Initialize temporary characteristics
                CHARS = ' '  ! array
                CBLR  = ' '
                CPDS  = ' '

C.................  Read in line of character data

                IF( BLRIN ) THEN
                    READ( FDEV, FILFMT, END=999 ) ID, CFIP, FCID, 
     &                  ( CHARS( J ), J=1, NCHARS ), CS, CBLR, CPDS

                ELSE

                    READ( FDEV, FILFMT, END=999 ) ID, CFIP, FCID, 
     &                  ( CHARS( J ), J=1, NCHARS ), CS, CPDS

                ENDIF

                IF( SCCFLAG ) CSCC  ( S ) = CS

                IF( BLRFLAG ) CBLRID( S ) = CBLR

                IF( PDSFLAG ) CPDESC( S ) = CPDS

                IF( CSRFLAG ) 
     &              CALL BLDCSRC( CFIP, FCID, CHARS(1), CHARS(2),
     &                            CHARS(3), CHARS(4), CHARS(5),
     &                            POLBLNK3, CSOURC( S ) )

                IF( ID .NE. S ) THEN

                    WRITE( MESG,94010 ) 'SMOKE ASCII inventory ' //
     &                     'file has been corrupted at source ID', ID
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF

            ENDDO  ! End loop on sources

            REWIND( FDEV )

        ENDIF

        IF( EFLAG ) CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C.........  For output, the number of characteristics must include the FIPS
C           code and the facility code
        NCHARS = NCHARS + 2 

        RETURN

999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of temporal' // CRLF() // BLANK5 //
     &         'cross reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )


        END


