
        SUBROUTINE RDINVCHR( CATEGORY, INFILE, FDEV, NSRC, 
     &                       NVARS, VNAMES )

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C     Reads in SMOKE inventory characteristics for any source category
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

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2            CRLF
        INTEGER                INDEX1

        EXTERNAL               CRLF, INDEX1

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY        ! source category
        CHARACTER(*), INTENT (IN) :: INFILE          ! inven logical file name
        INTEGER     , INTENT (IN) :: FDEV            ! inven ASCII file unit no.
        INTEGER     , INTENT (IN) :: NSRC            ! no. of sources
        INTEGER     , INTENT (IN) :: NVARS           ! no. of inven vars to read
        CHARACTER(*), INTENT (IN) :: VNAMES( NVARS ) ! variable names 

C...........   Index for unread variables
        INTEGER        UNREAD
        INTEGER        UNRIDX( NVARS )

C...........   Error message strings
        CHARACTER*23 ,PARAMETER :: PART1 = 'Error reading variable '
        CHARACTER*24 ,PARAMETER :: PART3 = ' from INVENTORY file'

C...........   Other local variables
        INTEGER          J, N, S    ! counters and indices

        INTEGER          ID         ! tmp smoke ID
        INTEGER          IOS        ! i/o status
        INTEGER          NCOL       ! number of columns in FDEV
        INTEGER          NC         ! number of plant characteristics

        LOGICAL       :: BLRIN   = .FALSE.  ! True: boiler is in input file
        LOGICAL       :: CSRFLAG = .FALSE.  ! True: source chars requested
        LOGICAL       :: EFLAG   = .FALSE.  ! True: error
        LOGICAL       :: BLRFLAG = .FALSE.  ! True: boilers requested
        LOGICAL       :: LNKFLAG = .FALSE.  ! True: link ID requested
        LOGICAL       :: ORSIN   = .FALSE.  ! Ture: DOE plant ID in input file
        LOGICAL       :: ORSFLAG = .FALSE.  ! True: DOE plant ID requested
        LOGICAL       :: PDSIN   = .FALSE.  ! True: plant desc in input file
        LOGICAL       :: PDSFLAG = .FALSE.  ! True: plant description requested
        LOGICAL       :: SCCFLAG = .FALSE.  ! True: SCC requested
        LOGICAL       :: VTPFLAG = .FALSE.  ! True: vehicle type requested

        CHARACTER*20           HEADER( 20 ) ! header fields
        CHARACTER*256          FILFMT ! ASCII file format after header
        CHARACTER*256          MESG   ! message buffer
        CHARACTER(LEN=BLRLEN3) CBLR   ! temporary boiler name
        CHARACTER(LEN=FIPLEN3) CFIP   ! temporary character FIPs code
        CHARACTER(LEN=CHRLEN3) CHARS( 5 ) ! temporary plant characteristics
        CHARACTER(LEN=LNKLEN3) CLNK   ! temporary link ID
        CHARACTER(LEN=ORSLEN3) CORS   ! temporary DOE plant ID
        CHARACTER(LEN=DSCLEN3) CPDS   ! temporary plant description
        CHARACTER(LEN=RWTLEN3) CRWT   ! temporary roadway type
        CHARACTER(LEN=SCCLEN3) CS     ! temporary scc
        CHARACTER(LEN=VIDLEN3) CVID   ! temporary vehicle type code
        CHARACTER(LEN=VTPLEN3) CVTP   ! tmp vehicle type
        CHARACTER(LEN=PLTLEN3) FCID   ! temporary facility code
        CHARACTER(LEN=IOVLEN3) INVAR  ! tmp inventory pollutant name

        CHARACTER*16 :: PROGNAME = 'RDINVCHR'   !  program name

C***********************************************************************
C   begin body of program RDINVCHR

        MESG = 'Reading source data from inventory file...'
        CALL M3MSG2( MESG )

C.........  Loop through input variables
C.........  Allocate memory and read the ones that are needed from I/O API file

        UNREAD = 0
        DO N = 1, NVARS

            INVAR = VNAMES( N )

            MESG = PART1 // INVAR( 1:LEN_TRIM( INVAR ) ) // PART3

            SELECT CASE( INVAR )

            CASE( 'IFIP' )
              ALLOCATE( IFIP( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'IFIP', PROGNAME )

              IF(.NOT. READSET(INFILE,'IFIP',ALLAYS3,1,0,0,IFIP)) THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'IRCLAS' )
              ALLOCATE( IRCLAS( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'IRCLAS', PROGNAME )

              IF(.NOT. READSET(INFILE,'IRCLAS',ALLAYS3,1,0,0,IRCLAS)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'ISIC' )
              ALLOCATE( ISIC( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'ISIC', PROGNAME )

              IF( .NOT. READSET(INFILE,'ISIC',ALLAYS3,1,0,0,ISIC )) THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'IVTYPE' )
              ALLOCATE( IVTYPE( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'IVTYPE', PROGNAME )

              IF(.NOT. READSET(INFILE,'IVTYPE',ALLAYS3,1,0,0,IVTYPE)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'CELLID' )
              ALLOCATE( CELLID( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'CELLID', PROGNAME )

              IF(.NOT. READSET(INFILE,'CELLID',ALLAYS3,1,0,0,CELLID)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'TZONES' )
              ALLOCATE( TZONES( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'TZONES', PROGNAME )

              IF(.NOT. READSET(INFILE,'TZONES',ALLAYS3,1,0,0,TZONES)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'TPFLAG' )
              ALLOCATE( TPFLAG( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'TPFLAG', PROGNAME )

              IF(.NOT. READSET(INFILE,'TPFLAG',ALLAYS3,1,0,0,TPFLAG)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'INVYR' )
              ALLOCATE( INVYR( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'INVYR', PROGNAME )

              IF( .NOT. READSET(INFILE,'INVYR',ALLAYS3,1,0,0,INVYR)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF
           
            CASE( 'IDIU' )
              ALLOCATE( IDIU( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'IDIU', PROGNAME )

              IF( .NOT. READSET(INFILE,'IDIU',ALLAYS3,1,0,0,IDIU )) THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'IWEK' )
              ALLOCATE( IWEK( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'IWEK', PROGNAME )

              IF( .NOT. READSET(INFILE,'IWEK',ALLAYS3,1,0,0,IWEK )) THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'XLOCA' )
              ALLOCATE( XLOCA( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )

              IF( .NOT. READSET(INFILE,'XLOCA',ALLAYS3,1,0,0,XLOCA)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'YLOCA' )
              ALLOCATE( YLOCA( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )

              IF( .NOT. READSET(INFILE,'YLOCA',ALLAYS3,1,0,0,YLOCA)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'XLOC1' )
              ALLOCATE( XLOC1( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'XLOC1', PROGNAME )

              IF( .NOT. READSET(INFILE,'XLOC1',ALLAYS3,1,0,0,XLOC1)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'YLOC1' )
              ALLOCATE( YLOC1( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'YLOC1', PROGNAME )

              IF( .NOT. READSET(INFILE,'YLOC1',ALLAYS3,1,0,0,YLOC1)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'XLOC2' )
              ALLOCATE( XLOC2( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'XLOC2', PROGNAME )

              IF( .NOT. READSET(INFILE,'XLOC2',ALLAYS3,1,0,0,XLOC2)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'YLOC2' )
              ALLOCATE( YLOC2( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'YLOC2', PROGNAME )

              IF( .NOT. READSET(INFILE,'YLOC2',ALLAYS3,1,0,0,YLOC2)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'SPEED' )
              ALLOCATE( SPEED( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'SPEED', PROGNAME )

              IF( .NOT. READSET(INFILE,'SPEED',ALLAYS3,1,0,0,SPEED)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF
            
            CASE( 'STKHT' )
              ALLOCATE( STKHT( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'STKHT', PROGNAME )

              IF( .NOT. READSET(INFILE,'STKHT',ALLAYS3,0,0,STKHT)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'STKDM' )
              ALLOCATE( STKDM( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'STKDM', PROGNAME )

              IF( .NOT. READSET(INFILE,'STKDM',ALLAYS3,0,0,STKDM)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'STKTK' )
              ALLOCATE( STKTK( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'STKTK', PROGNAME )

              IF( .NOT. READSET(INFILE,'STKTK',ALLAYS3,1,0,0,STKTK)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'STKVE' )
              ALLOCATE( STKVE( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'STKVE', PROGNAME )

              IF( .NOT. READSET(INFILE,'STKVE',ALLAYS3,1,0,0,STKVE)) 
     &            THEN
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
     &                 '." NOTE: ASCII inventory file is not opened.'
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
                    ALLOCATE( CSCC( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CSCC', PROGNAME )

                CASE( 'CORIS' )
                    ORSFLAG = .TRUE.
                    ALLOCATE( CORIS( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CORIS', PROGNAME )

                CASE( 'CBLRID' )
                    BLRFLAG = .TRUE. 
                    ALLOCATE( CBLRID( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CBLRID', PROGNAME )

                CASE( 'CLINK' )
                    LNKFLAG = .TRUE. 
                    ALLOCATE( CLINK( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CLINK', PROGNAME )

                CASE( 'CPDESC' )
                    PDSFLAG = .TRUE. 
                    ALLOCATE( CPDESC( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CPDESC', PROGNAME )

                CASE( 'CSOURC' )
                    CSRFLAG = .TRUE. 
                    ALLOCATE( CSOURC( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CSOURC', PROGNAME )

                CASE( 'CVTYPE' )
                    VTPFLAG = .TRUE. 
                    ALLOCATE( CVTYPE( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CVTYPE', PROGNAME )

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

CC.............  Depending on source category, read in and store ASCII source
C               characteristics

            SELECT CASE ( CATEGORY )
            CASE ( 'AREA' )

                DO S = 1, NSRC

C.....................  Read in line of character data
                    READ( FDEV, FILFMT, END=999 ) ID, CFIP, CS

                    IF( SCCFLAG ) CSCC( S ) = CS

                    IF( CSRFLAG ) 
     &                  CALL BLDCSRC( CFIP, CS, CHRBLNK3, CHRBLNK3,
     &                                CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                                POLBLNK3, CSOURC( S ) )

                    CALL CHECK_CORRUPTED

                END DO  ! End loop on sources

            CASE ( 'MOBILE' )

                DO S = 1, NSRC

C.....................  Initialize temporary characteristics
                    CVID  = ' '
                    CLNK  = ' '

C.....................  Read in line of character data
                    READ( FDEV, FILFMT, END=999 ) ID, CFIP, CRWT, 
     &                                            CLNK, CVID, CS, CVTP

                    IF( SCCFLAG ) CSCC  ( S ) = CS

                    IF( VTPFLAG ) CVTYPE( S ) = CVTP

                    IF( LNKFLAG ) CLINK ( S ) = CLNK

                    IF( CSRFLAG ) 
     &                  CALL BLDCSRC( CFIP, CRWT, CLNK, CVID,
     &                                CS, CHRBLNK3, CHRBLNK3,
     &                                POLBLNK3, CSOURC( S ) )

                    CALL CHECK_CORRUPTED

                END DO  ! End loop on sources

            CASE ( 'POINT' )

C.................  Determine number of plant characteristics. Assumes SCC is
C                   the first variable after the final source characteristic 
                J = INDEX1( 'SCC', NCOL, HEADER )
                NC = J - 4  ! Four because Source ID, FIPS, Plant ID, & SCC 

C.................  Determine if DOE plant ID is present
                J = INDEX1( 'DOE plant ID', NCOL, HEADER )
                ORSIN = ( J .GT. 0 )

C.................  Determine if boiler is present
                J = INDEX1( 'Boiler code', NCOL, HEADER )
                BLRIN = ( J .GT. 0 )

C.................  Determine if plant description is present
                J = INDEX1( 'Facility description', NCOL, HEADER )
                PDSIN = ( J .GT. 0 )

C.................  If DOE plant ID not present but has been requested, then 
C                   internal err
                IF( .NOT. ORSIN .AND. ORSFLAG ) THEN

                    MESG = 'WARNING: ORIS ID requested, but ' //
     &                     'is not present in ASCII inventory file'
                    CALL M3MSG2( MESG )

                END IF

C.................  If boiler not present but has been requested, then 
C                   internal err
                IF( .NOT. BLRIN .AND. BLRFLAG ) THEN

                    MESG = 'WARNING: Boiler requested, but ' //
     &                     'is not present in ASCII inventory file'
                    CALL M3MSG2( MESG )

                END IF

                IF( EFLAG ) THEN

                    MESG = 'ERROR: Problem reading inventory file(s)'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF
                
                CORS = ' '
                CBLR  = ' '
                CPDS  = ' '

                DO S = 1, NSRC

C.....................  Initialize temporary characteristics
                    CHARS = ' '  ! array

C.....................  Read in line of character data

                    IF( ORSIN .AND. BLRIN .AND. PDSIN ) THEN
                        READ( FDEV, FILFMT, END=999 ) ID, CFIP, FCID, 
     &                      ( CHARS( J ), J=1,NC ), CS, CORS, CBLR, CPDS

                    ELSE IF( ORSIN .AND. PDSIN ) THEN
                        READ( FDEV, FILFMT, END=999 ) ID, CFIP, FCID, 
     &                      ( CHARS( J ), J=1, NC ), CS, CORS, CPDS

                    ELSE IF( ORSIN .AND. BLRIN ) THEN
                        READ( FDEV, FILFMT, END=999 ) ID, CFIP, FCID, 
     &                      ( CHARS( J ), J=1,NC ), CS, CORS, CBLR
     
                    ELSE IF( PDSIN ) THEN
                        READ( FDEV, FILFMT, END=999 ) ID, CFIP, FCID,
     &                      ( CHARS( J ), J=1, NC ), CS, CPDS      

                    ELSE 
                        READ( FDEV, FILFMT, END=999 ) ID, CFIP, FCID, 
     &                      ( CHARS( J ), J=1, NC ), CS

                    ENDIF

                    IF( SCCFLAG ) CSCC  ( S ) = CS

                    IF( ORSFLAG ) CORIS ( S ) = ADJUSTR( CORS )

                    IF( BLRFLAG ) CBLRID( S ) = ADJUSTR( CBLR )

                    IF( PDSFLAG ) CPDESC( S ) = CPDS

                    IF( CSRFLAG ) 
     &                  CALL BLDCSRC( CFIP, FCID, CHARS(1), CHARS(2),
     &                                CHARS(3), CHARS(4), CHARS(5),
     &                                POLBLNK3, CSOURC( S ) )

                    CALL CHECK_CORRUPTED

                END DO  ! End loop on sources

            END SELECT

            REWIND( FDEV )

        ENDIF

        IF( EFLAG ) CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        RETURN

999     MESG = 'End of file reached unexpectedly. ' //
     &         'Check format of ASCII' // CRLF() // BLANK5 //
     &         'inventory file.'

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram checks to make sure the ASCII inventory
C               file has not been corrupted, and if it has, stops the program
            SUBROUTINE CHECK_CORRUPTED

C......................................................................

            IF( ID .NE. S ) THEN

                WRITE( MESG,94010 ) 'SMOKE ASCII inventory ' //
     &                      'file has been corrupted at source ID', ID
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

C--------------------  FORMAT  STATEMENTS   ----------------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I10, :, 1X ) )

            END SUBROUTINE CHECK_CORRUPTED

C----------------------------------------------------------------------

        END SUBROUTINE RDINVCHR


