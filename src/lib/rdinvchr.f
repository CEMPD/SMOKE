
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

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: CIFIP, IRCLAS, IVTYPE, CELLID, TZONES,
     &                      TPFLAG, INVYR, IDIU, IWEK, XLOCA, YLOCA, 
     &                      XLOC1, YLOC1, XLOC2, YLOC2, SPEED, STKHT,
     &                      STKDM, STKTK, STKVE, CSCC, CORIS, CBLRID,
     &                      CLINK, CPDESC, CSOURC, CVTYPE, CMACT,
     &                      CNAICS, CSRCTYP, CERPTYP, CNEIUID, CEXTORL,
     &                      CINTGR, CISIC, FUGHGT, FUGWID, FUGLEN, FUGANG

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
C       CHARACTER(2)           CRLF
C       INTEGER                INDEX1

C        EXTERNAL               CRLF, INDEX1

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
        CHARACTER(23), PARAMETER :: PART1 = 'Error reading variable '
        CHARACTER(24), PARAMETER :: PART3 = ' from INVENTORY file'

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
        LOGICAL       :: ERPIN   = .FALSE.  ! True: emission release type in input file
        LOGICAL       :: ERPFLAG = .FALSE.  ! True: emission release point type requested
        LOGICAL       :: FIPFLAG = .FALSE.  ! True: geographic code requested
        LOGICAL       :: LNKFLAG = .FALSE.  ! True: link ID requested
        LOGICAL       :: EXTIN   = .FALSE.  ! True: additional extended orl vars in input file
        LOGICAL       :: EXTFLAG = .FALSE.  ! True: additional extended orl vars requested
        LOGICAL       :: ITGIN   = .FALSE.  ! True: INTG(integrate) code in input file
        LOGICAL       :: ITGFLAG = .FALSE.  ! True: INTG(integrate) code requested
        LOGICAL       :: MCTIN   = .FALSE.  ! True: MACT code in input file
        LOGICAL       :: MACFLAG = .FALSE.  ! True: MACT code requested
        LOGICAL       :: NAIIN   = .FALSE.  ! True: NAICS code in input file
        LOGICAL       :: NAIFLAG = .FALSE.  ! True: NAICS code requested
        LOGICAL       :: NEIIN   = .FALSE.  ! True: NEI Unique ID in input file
        LOGICAL       :: NEIFLAG = .FALSE.  ! True: NEI Unique ID requested
        LOGICAL       :: ORSIN   = .FALSE.  ! Ture: DOE plant ID in input file
        LOGICAL       :: ORSFLAG = .FALSE.  ! True: DOE plant ID requested
        LOGICAL       :: PDSIN   = .FALSE.  ! True: plant desc in input file
        LOGICAL       :: PDSFLAG = .FALSE.  ! True: plant description requested
        LOGICAL       :: SCCFLAG = .FALSE.  ! True: SCC requested
        LOGICAL       :: SICFLAG = .FALSE.  ! True: SIC requested
        LOGICAL       :: SICIN   = .FALSE.  ! True: SIC code in input file
        LOGICAL       :: STPIN   = .FALSE.  ! True: source type code in input file
        LOGICAL       :: STPFLAG = .FALSE.  ! True: source type code requested
        LOGICAL       :: VTPFLAG = .FALSE.  ! True: vehicle type requested

        CHARACTER(20)      HEADER( 20 ) ! header fields
        CHARACTER(256)     FILFMT ! ASCII file format after header
        CHARACTER(256)     MESG   ! message buffer
        CHARACTER(CHRLEN3) CHARS( 5 ) ! temporary plant characteristics
        CHARACTER(BLRLEN3) :: CBLR = ' '   ! temporary boiler name
        CHARACTER(FIPLEN3) :: CFIP = ' '   ! temporary character FIPs code
        CHARACTER(LNKLEN3) :: CLNK = ' '   ! temporary link ID
        CHARACTER(NEILEN3) :: CNEI = ' '   ! NEI unique ID
        CHARACTER(EXTLEN3) :: CEXT = ' '   ! Extended ORL vars
        CHARACTER(ORSLEN3) :: CORS = ' '   ! temporary DOE plant ID
        CHARACTER(DSCLEN3) :: CPDS = ' '   ! temporary plant description
        CHARACTER(RWTLEN3) :: CRWT = ' '   ! temporary roadway type
        CHARACTER(SCCLEN3) :: CS   = ' '   ! temporary scc
        CHARACTER(SICLEN3) :: CSIC = ' '   ! temporary SIC
        CHARACTER(VIDLEN3) :: CVID = ' '   ! temporary vehicle type code
        CHARACTER(VTPLEN3) :: CVTP = ' '   ! tmp vehicle type
        CHARACTER(MACLEN3) :: CMT  = ' '   ! tmp MACT code
        CHARACTER(NAILEN3) :: CNAI = ' '   ! tmp NAICS code
        CHARACTER(STPLEN3) :: CSTP = ' '   ! tmp source type code
        CHARACTER(ERPLEN3) :: CERP = ' '   ! tmp emission release point code
        CHARACTER(INTLEN3) :: CINT = ' '   ! tmp integrate code
        CHARACTER(PLTLEN3) :: FCID = ' '   ! temporary facility code
        CHARACTER(IOVLEN3) INVAR  ! tmp inventory pollutant name

C.........  File format handling
        INTEGER, PARAMETER :: MXITEMS = 40
        INTEGER FMTITEM
        CHARACTER(8)  FMTSEGS( MXITEMS )
        CHARACTER(10) FMTSEG

        CHARACTER(16) :: PROGNAME = 'RDINVCHR'   !  program name

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

            CASE( 'IRCLAS' )
              ALLOCATE( IRCLAS( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'IRCLAS', PROGNAME )

              IF(.NOT. READSET(INFILE,'IRCLAS',ALLAYS3,1,0,0,IRCLAS)) 
     &            THEN
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

              IF( .NOT. READSET(INFILE,'STKHT',ALLAYS3,1,0,0,STKHT)) 
     &            THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'STKDM' )
              ALLOCATE( STKDM( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'STKDM', PROGNAME )

              IF( .NOT. READSET(INFILE,'STKDM',ALLAYS3,1,0,0,STKDM)) 
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

            CASE( 'FUGHGT' )
              ALLOCATE( FUGHGT( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'FUGHGT', PROGNAME )

              IF( .NOT. READSET(INFILE,'FUG_HEIGHT',ALLAYS3,1,0,0,FUGHGT)) THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'FUGWID' )
              ALLOCATE( FUGWID( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'FUGWID', PROGNAME )

              IF( .NOT. READSET(INFILE,'FUG_WIDTH',ALLAYS3,1,0,0,FUGWID)) THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'FUGLEN' )
              ALLOCATE( FUGLEN( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'FUGLEN', PROGNAME )

              IF( .NOT. READSET(INFILE,'FUG_LENGTH',ALLAYS3,1,0,0,FUGLEN)) THEN
                  CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
              ENDIF

            CASE( 'FUGANG' )
              ALLOCATE( FUGANG( NSRC ), STAT=IOS )
              CALL CHECKMEM( IOS, 'FUGANG', PROGNAME )

              IF( .NOT. READSET(INFILE,'FUG_ANGLE',ALLAYS3,1,0,0,FUGANG)) THEN
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
                
                CASE( 'CIFIP' )
                    FIPFLAG = .TRUE.
                    ALLOCATE( CIFIP( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CIFIP', PROGNAME )

                CASE( 'CSCC' )
                    SCCFLAG = .TRUE. 
                    ALLOCATE( CSCC( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CSCC', PROGNAME )
                
                CASE( 'CISIC' )
                    SICFLAG = .TRUE.
                    ALLOCATE( CISIC( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CISIC', PROGNAME )
                    
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
                    
                CASE( 'CMACT' )
                    MACFLAG = .TRUE.
                    ALLOCATE( CMACT( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CMACT', PROGNAME )
                    
                CASE( 'CNAICS' )
                    NAIFLAG = .TRUE.
                    ALLOCATE( CNAICS( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CNAICS', PROGNAME )
                    
                CASE( 'CSRCTYP' )
                    STPFLAG = .TRUE.
                    ALLOCATE( CSRCTYP( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CSRCTYP', PROGNAME )
                
                CASE( 'CERPTYP' )
                    ERPFLAG = .TRUE.
                    ALLOCATE( CERPTYP( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CERPTYP', PROGNAME )

                CASE( 'CNEIUID' )
                    NEIFLAG = .TRUE.
                    ALLOCATE( CNEIUID( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CNEIUID', PROGNAME )

                CASE( 'CINTGR' )
                    ITGFLAG = .TRUE.
                    ALLOCATE( CINTGR( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CINTGR', PROGNAME )

                CASE( 'CEXTORL' )
                    EXTFLAG = .TRUE.
                    ALLOCATE( CEXTORL( NSRC ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CEXTORL', PROGNAME )

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

C.............  Read in number of header lines and file format
            READ( FDEV, * ) NCOL, FILFMT
            CALL PARSLINE( FILFMT( 2:LEN_TRIM( FILFMT )-1 ), MXITEMS, FMTSEGS )

C.............  Read past header
            DO J = 1, NCOL
                READ( FDEV, '(A)' ) HEADER( J )
            ENDDO

CC............  Read in and store common ASCII source characteristics 
C               over all source categories

C..............  Determine if source type code is present
            J = INDEX1( 'Source type code', NCOL, HEADER )
            STPIN = ( J > 0 )

C.............  Determine if plant description is present
            J = INDEX1( 'Integrate flag', NCOL, HEADER )
            ITGIN = ( J .GT. 0 )

C.............  Determine if plant description is present
            J = INDEX1( 'Additional extended', NCOL, HEADER )
            EXTIN = ( J .GT. 0 )

C.............  If source type code not present but has been requested,
C               deallocate memory for array
            IF( .NOT. STPIN .AND. STPFLAG ) THEN

                MESG = 'WARNING: Source type code requested, but '//
     &                 'is not present in ASCII inventory file'
c                CALL M3MSG2( MESG )

                DEALLOCATE( CSRCTYP )
                NULLIFY( CSRCTYP )
            END IF

C.............  If integrate columns not present but has been requested,
C               deallocate memory for array
            IF( .NOT. ITGIN .AND. ITGFLAG ) THEN

                MESG = 'WARNING: Integrate flag requested, ' //
     &                 'but is not present in ASCII inventory file'
c                CALL M3MSG2( MESG )

                DEALLOCATE( CINTGR )
                NULLIFY( CINTGR )
            END IF

C.............  If extended columns not present but has been requested,
C               deallocate memory for array
            IF( .NOT. EXTIN .AND. EXTFLAG ) THEN

                MESG = 'WARNING: Additional extended requested, ' //
     &                 'but is not present in ASCII inventory file'
c                CALL M3MSG2( MESG )

                DEALLOCATE( CEXTORL )
                NULLIFY( CEXTORL )
            END IF

CC............  Depending on source category, read in and store ASCII source
C               characteristics

            SELECT CASE ( CATEGORY )
            CASE ( 'AREA' )

C.................  Determine if MACT code is present
                J = INDEX1( 'MACT code', NCOL, HEADER )
                MCTIN = ( J > 0 )
                
C.................  Determine if NAICS code is present
                J = INDEX1( 'NAICS code', NCOL, HEADER )
                NAIIN = ( J > 0 )
                
C.................  Determine if DOE plant ID is present
                J = INDEX1( 'DOE plant ID', NCOL, HEADER )
                ORSIN = ( J .GT. 0 )

C.................  Determine if SIC code is present
                J = INDEX1( 'SIC', NCOL, HEADER )
                SICIN = ( J > 0 )

C.................  If MACT not present but has been requested, then 
C                   internal err
                IF( .NOT. MCTIN .AND. MACFLAG ) THEN

                    MESG = 'WARNING: MACT requested, but ' //
     &                     'is not present in ASCII inventory file'
c                    CALL M3MSG2( MESG )

                    DEALLOCATE( CMACT )
                    NULLIFY( CMACT )

                END IF

C.................  If boiler not present but has been requested, then 
C                   internal err
                IF( .NOT. NAIIN .AND. NAIFLAG ) THEN

                    MESG = 'WARNING: NAICS requested, but ' //
     &                     'is not present in ASCII inventory file'
c                    CALL M3MSG2( MESG )

                    DEALLOCATE( CNAICS )
                    NULLIFY( CNAICS )

                END IF

C.................  If DOE plant ID not present but has been requested, then
C                   internal err
                IF( .NOT. ORSIN .AND. ORSFLAG ) THEN

                    MESG = 'WARNING: ORIS ID requested, but ' //
     &                     'is not present in ASCII inventory file'
C                    CALL M3MSG2( MESG )

                    DEALLOCATE( CORIS )
                END IF

C.................  If SIC not present but has been requested, then
C                   internal err
                IF( .NOT. SICIN .AND. SICFLAG ) THEN

                    MESG = 'WARNING: SIC requested, but ' //
     &                     'is not present in ASCII inventory file'
C                    CALL M3MSG2( MESG )

                    DEALLOCATE( CISIC )
                    NULLIFY( CISIC )

                END IF

                DO S = 1, NSRC

C.....................  Read source information from record of inventory file
                    FMTITEM = 1
                    CALL BUILD_FMTSEG( FMTSEGS( FMTITEM ) )
                    READ( FDEV, FMTSEG, ADVANCE='NO', END=999 ) ID

                    CALL READ_NEXT_VAL( CFIP )
                    CALL READ_NEXT_VAL( CS )
                    IF( STPIN ) CALL READ_NEXT_VAL( CSTP )
                    IF( MCTIN ) CALL READ_NEXT_VAL( CMT )
                    IF( NAIIN ) CALL READ_NEXT_VAL( CNAI )
                    IF( ITGIN ) CALL READ_NEXT_VAL( CINT )
                    IF( SICIN ) CALL READ_NEXT_VAL( CSIC )
                    IF( EXTIN ) CALL READ_NEXT_VAL( CEXT )

C.....................  Advance to next line
                    READ( FDEV, *, END=999 )

                    IF( FIPFLAG ) CIFIP( S ) = CFIP

                    IF( SCCFLAG ) CSCC( S ) = CS

                    IF( MACFLAG .AND. MCTIN ) CMACT( S ) = CMT
                    
                    IF( NAIFLAG .AND. NAIIN ) CNAICS( S ) = CNAI
                    
                    IF( STPFLAG .AND. STPIN ) CSRCTYP( S ) = CSTP

                    IF( ITGFLAG .AND. ITGIN ) CINTGR( S ) = CINT

                    IF( EXTFLAG .AND. EXTIN ) CEXTORL( S ) = 
     &                                             ADJUSTL( CEXT )

                    IF( SICFLAG .AND. SICIN ) CISIC( S ) = CSIC

                    IF( CSRFLAG ) 
     &                  CALL BLDCSRC( CFIP, CS, CHRBLNK3, CHRBLNK3,
     &                                CHRBLNK3, CHRBLNK3, CHRBLNK3,
     &                                POLBLNK3, CSOURC( S ) )

                    CALL CHECK_CORRUPTED

                END DO  ! End loop on sources

            CASE ( 'MOBILE' )

                DO S = 1, NSRC

C.....................  Read source information from record of inventory file
                    FMTITEM = 1
                    CALL BUILD_FMTSEG( FMTSEGS( FMTITEM ) )
                    READ( FDEV, FMTSEG, ADVANCE='NO', END=999 ) ID
                    
                    CALL READ_NEXT_VAL( CFIP )
                    CALL READ_NEXT_VAL( CRWT )
                    CALL READ_NEXT_VAL( CLNK )
                    CALL READ_NEXT_VAL( CVID )
                    CALL READ_NEXT_VAL( CS )
                    CALL READ_NEXT_VAL( CVTP )
                    IF( STPIN ) CALL READ_NEXT_VAL( CSTP )
                    IF( ITGIN ) CALL READ_NEXT_VAL( CINT )
                    IF( EXTIN ) CALL READ_NEXT_VAL( CEXT )

C.....................  Advance to next line
                    READ( FDEV, *, END=999 )
                    
                    IF( FIPFLAG ) CIFIP ( S ) = CFIP

                    IF( SCCFLAG ) CSCC  ( S ) = CS

                    IF( VTPFLAG ) CVTYPE( S ) = CVTP

                    IF( LNKFLAG ) CLINK ( S ) = CLNK

                    IF( STPFLAG .AND. STPIN ) CSRCTYP( S ) = CSTP

                    IF( ITGFLAG .AND. ITGIN ) CINTGR( S ) = CINT

                    IF( EXTFLAG .AND. EXTIN ) CEXTORL( S ) = 
     &                                             ADJUSTL( CEXT )

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

C.................  Determine if MACT code is present
                J = INDEX1( 'MACT code', NCOL, HEADER )
                MCTIN = ( J > 0 )
                
C.................  Determine if NAICS code is present
                J = INDEX1( 'NAICS code', NCOL, HEADER )
                NAIIN = ( J > 0 )
                
C.................  Determine if emission release code is present
                J = INDEX1( 'Emission release pt', NCOL, HEADER )
                ERPIN = ( J > 0 )
                
C.................  Determine if plant description is present
                J = INDEX1( 'Facility description', NCOL, HEADER )
                PDSIN = ( J .GT. 0 )

C.................  Determine if plant description is present
                J = INDEX1( 'NEI unique ID', NCOL, HEADER )
                NEIIN = ( J .GT. 0 )

C.................  Determine if SIC code is present
                J = INDEX1( 'SIC', NCOL, HEADER )
                SICIN = ( J > 0 )

C.................  If MACT not present but has been requested, then 
C                   internal err
                IF( .NOT. MCTIN .AND. MACFLAG ) THEN

                    MESG = 'WARNING: MACT requested, but ' //
     &                     'is not present in ASCII inventory file'
C                    CALL M3MSG2( MESG )

                    DEALLOCATE( CMACT )
                    NULLIFY( CMACT )

                END IF

C.................  If boiler not present but has been requested, then 
C                   internal err
                IF( .NOT. NAIIN .AND. NAIFLAG ) THEN

                    MESG = 'WARNING: NAICS requested, but ' //
     &                     'is not present in ASCII inventory file'
C                    CALL M3MSG2( MESG )

                    DEALLOCATE( CNAICS )
                    NULLIFY( CNAICS )

                END IF

C.................  If DOE plant ID not present but has been requested, then 
C                   internal err
                IF( .NOT. ORSIN .AND. ORSFLAG ) THEN

                    MESG = 'WARNING: ORIS ID requested, but ' //
     &                     'is not present in ASCII inventory file'
C                    CALL M3MSG2( MESG )

                END IF

C.................  If boiler not present but has been requested, then 
C                   internal err
                IF( .NOT. BLRIN .AND. BLRFLAG ) THEN

                    MESG = 'WARNING: Boiler requested, but ' //
     &                     'is not present in ASCII inventory file'
C                    CALL M3MSG2( MESG )

                END IF

C.................  If NEI unique ID not present but has been requested, then 
C                   internal err
                IF( .NOT. NEIIN .AND. NEIFLAG ) THEN

                    MESG = 'WARNING: NEI unique ID requested, but ' //
     &                     'is not present in ASCII inventory file'
C                    CALL M3MSG2( MESG )

                    DEALLOCATE( CNEIUID )

                END IF

C.................  If emission release code not present but has been requested,
C                   deallocate memory for array
                IF( .NOT. ERPIN .AND. ERPFLAG ) THEN

                    MESG = 'WARNING: Emission release point requested, '
     &                  // 'but is not present in ASCII inventory file'
C                    CALL M3MSG2( MESG )

                    DEALLOCATE( CERPTYP )

                END IF

C.................  If SIC not present but has been requested, then 
C                   internal err
                IF( .NOT. SICIN .AND. SICFLAG ) THEN

                    MESG = 'WARNING: SIC requested, but ' //
     &                     'is not present in ASCII inventory file'
C                    CALL M3MSG2( MESG )

                    DEALLOCATE( CISIC )
                    NULLIFY( CISIC )

                END IF
                
                IF( EFLAG ) THEN

                    MESG = 'ERROR: Problem reading inventory file(s)'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF
                
                DO S = 1, NSRC

C.....................  Initialize temporary characteristics
                    CHARS = ' '  ! array

C.....................  Read source information from record of inventory file
                    FMTITEM = 1
                    CALL BUILD_FMTSEG( FMTSEGS( FMTITEM ) )
                    READ( FDEV, FMTSEG, ADVANCE='NO', END=999 ) ID
                    
                    CALL READ_NEXT_VAL( CFIP )
                    CALL READ_NEXT_VAL( FCID )
                    DO J = 1, NC
                        CALL READ_NEXT_VAL( CHARS( J ) )
                    END DO
                    CALL READ_NEXT_VAL( CS )
                    IF( ORSIN ) CALL READ_NEXT_VAL( CORS )
                    IF( BLRIN ) CALL READ_NEXT_VAL( CBLR )
                    IF( MCTIN ) CALL READ_NEXT_VAL( CMT )
                    IF( NAIIN ) CALL READ_NEXT_VAL( CNAI )
                    IF( STPIN ) CALL READ_NEXT_VAL( CSTP )
                    IF( ERPIN ) CALL READ_NEXT_VAL( CERP )
                    IF( PDSIN ) CALL READ_NEXT_VAL( CPDS )
                    IF( NEIIN ) CALL READ_NEXT_VAL( CNEI )
                    IF( ITGIN ) CALL READ_NEXT_VAL( CINT )
                    IF( SICIN ) CALL READ_NEXT_VAL( CSIC )
                    IF( EXTIN ) CALL READ_NEXT_VAL( CEXT )

C.....................  Advance to next line
                    READ( FDEV, *, END=999 )

                    IF( FIPFLAG ) CIFIP ( S ) = CFIP

                    IF( SCCFLAG ) CSCC  ( S ) = CS

                    IF( ORSFLAG .AND. ORSIN ) CORIS ( S ) = ADJUSTR( CORS )

                    IF( BLRFLAG .AND. BLRIN ) CBLRID( S ) = ADJUSTR( CBLR )

                    IF( ITGFLAG .AND. ITGIN ) CINTGR( S ) = CINT

                    IF( NEIFLAG .AND. NEIIN ) CNEIUID( S ) = 
     &                                             ADJUSTR( CNEI )

                    IF( EXTFLAG .AND. EXTIN ) CEXTORL( S ) = 
     &                                             ADJUSTL( CEXT )

                    IF( MACFLAG .AND. MCTIN ) CMACT( S ) = CMT

                    IF( NAIFLAG .AND. NAIIN ) CNAICS( S ) = CNAI

                    IF( STPFLAG .AND. STPIN ) CSRCTYP( S ) = CSTP

                    IF( ERPFLAG .AND. ERPIN ) CERPTYP( S ) = CERP

                    IF( PDSFLAG .AND. PDSIN ) CPDESC( S ) = CPDS

                    IF( SICFLAG .AND. SICIN ) CISIC( S ) = CSIC

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

C.............  This internal subprogram build a format specification string
C               from the given format fragment
            SUBROUTINE BUILD_FMTSEG( FMTFRAGMENT )

C.................  Subprogram arguments
                CHARACTER(*) FMTFRAGMENT

C......................................................................

                WRITE( FMTSEG, '(A1, A, A1)' ) '(', TRIM( FMTFRAGMENT ), ')'
            
            END SUBROUTINE BUILD_FMTSEG

C----------------------------------------------------------------------

C.............  This internal subprogram reads the next string value from
C               the current record in the inventory file
            SUBROUTINE READ_NEXT_VAL( VARIABLE )

C.................  Subprogram arguments
                CHARACTER(*) VARIABLE

C......................................................................

C.................  Skip blank space
                FMTITEM = FMTITEM + 1
                CALL BUILD_FMTSEG( FMTSEGS( FMTITEM ) )
                READ( FDEV, FMTSEG, ADVANCE='NO', END=998 )

C.................  Extract string format and read variable            
                FMTITEM = FMTITEM + 1
                CALL BUILD_FMTSEG( FMTSEGS( FMTITEM ) )
                READ( FDEV, FMTSEG, ADVANCE='NO', END=998 ) VARIABLE
                
                RETURN

998             MESG = 'End of file reached unexpectedly. ' //
     &                 'Check format of ASCII' // CRLF() // BLANK5 //
     &                 'inventory file.'

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            
            END SUBROUTINE READ_NEXT_VAL

        END SUBROUTINE RDINVCHR


