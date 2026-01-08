
        SUBROUTINE GETSINFO( ENAME )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION: 
C     This subroutine get the source-category-specific source information from
C     the inventory header, which should be read prior to calling this
C     routine.
C
C  PRECONDITIONS REQUIRED:
C     SMK_AVEDAY_YN e.v. to have been checked to set the value of INVPIDX
C
C  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
C
C  REVISION  HISTORY:
C       Created 3/99 by M Houyoux
C       09/2025 by HT UNC-IE:  Use M3UTILIO
C
C***********************************************************************
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
C****************************************************************************
        USE M3UTILIO

C.........  MODULES for public variables

C.........  This module is required by the FileSetAPI
        USE MODFILESET

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, NCHARS, JSCC, JSTACK, CATEGORY,
     &                     LSCCEND, RSCCBEG, NPPOL, PLTIDX, MXCHRS,
     &                     SCCLEV1, SCCLEV2, SCCLEV3, SCCLEV4,
     &                     NIPOL, NIACT, NPACT, NIPPA, BYEAR,
     &                     ATTRUNIT, EANAM, EAREAD, EAUNIT, EADESC,
     &                     NMAP, ACTVTY, SC_BEGP, SC_ENDP, EINAM,
     &                     INVPIDX, MAPFIL, MAPNAM

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C       INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C...........   EXTERNAL FUNCTIONS:
c       CHARACTER(2)    CRLF
c       LOGICAL         ENVYN
c       INTEGER         GETIFDSC
c       LOGICAL         SETENVVAR

c       EXTERNAL        CRLF, ENVYN, GETIFDSC, SETENVVAR
        INTEGER, EXTERNAL :: GETIFDSC

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME   ! inventory file logical file name

C...........   LOCAL VARIABLES their descriptions:
        INTEGER       I, J, K     ! counters and indices 
        INTEGER       IOS         ! memory allocation status
        INTEGER       NVAR        ! number of non-pollutant variables 

        LOGICAL    :: EFLAG = .FALSE.  ! true: error found
        LOGICAL       LAVEDAY     ! true: use average day emissions

        CHARACTER(IOFLEN3)  TMPNAME ! tmp logical file name for map data files
        CHARACTER(256) MESG         ! Message buffer

        CHARACTER(16) :: PROGNAME = 'GETSINFO'    ! Program name

C***********************************************************************
C   begin body of subroutine GETSINFO

C.........  Get header description of inventory file, error if problem
        IF( .NOT. DESCSET( ENAME, ALLFILES ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set the number of sources 
        NSRC = NROWS3D

C.........  Set the number of source characteristics and allocate memory for
C           the field positions
        NCHARS   = GETIFDSC( FDESC3D, '/NUMBER CHARS/', .TRUE. )
        JSCC     = GETIFDSC( FDESC3D, '/SCC POSITION/', .TRUE. )
        JSTACK   = GETIFDSC( FDESC3D, '/STACK POSITION/', .FALSE. )

        IF( .NOT. ALLOCATED( SC_BEGP ) ) THEN
            ALLOCATE( SC_BEGP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SC_BEGP', PROGNAME )
            ALLOCATE( SC_ENDP( NCHARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SC_ENDP', PROGNAME )
        END IF

C.........  Set source-category-specific information
        SELECT CASE( CATEGORY )

        CASE( 'AREA' )

            LSCCEND  = SCCEXPLEN3 + 7
            RSCCBEG  = SCCEXPLEN3 + 8
            NPPOL    = NARPPOL3
            PLTIDX   = MXARCHR3
            MXCHRS   = MXARCHR3

            DO I = 1, NCHARS
                SC_BEGP( I ) = ARBEGL3( I )
                SC_ENDP( I ) = ARENDL3( I )
            END DO

            SCCLEV1 = SCCEXPLEN3 + 2
            SCCLEV2 = SCCEXPLEN3 + 4
            SCCLEV3 = SCCEXPLEN3 + 7
            SCCLEV4 = SCCEXPLEN3 + 10

        CASE( 'MOBILE' )

            LSCCEND  = SCCLEN3 - VIDLEN3   ! For reset SCCs CRWT // VCID
            RSCCBEG  = LSCCEND + 1
            NPPOL    = NMBPPOL3
            PLTIDX   = 2
            MXCHRS   = MXMBCHR3

            DO I = 1, NCHARS
                SC_BEGP( I ) = MBBEGL3( I )
                SC_ENDP( I ) = MBENDL3( I )
            END DO

            SCCLEV1 = SCCEXPLEN3 + 2
            SCCLEV2 = SCCEXPLEN3 + 4
            SCCLEV3 = SCCEXPLEN3 + 7
            SCCLEV4 = SCCEXPLEN3 + 10
            
        CASE( 'POINT' )

            LSCCEND  = SCCEXPLEN3 + 5
            RSCCBEG  = SCCEXPLEN3 + 6
            PLTIDX   = 2
            NPPOL    = NPTPPOL3
            MXCHRS   = MXPTCHR3

            DO I = 1, NCHARS
                SC_BEGP( I ) = PTBEGL3( I )
                SC_ENDP( I ) = PTENDL3( I )
            END DO
            
            SCCLEV1 = SCCEXPLEN3 + 3   ! Assumes right-justified 10-digit w/ first 2 zero
            SCCLEV2 = SCCEXPLEN3 + 5
            SCCLEV3 = SCCEXPLEN3 + 8
            SCCLEV4 = SCCEXPLEN3 + 10
            
        CASE DEFAULT
            MESG = 'Category ' // CATEGORY // ' not known in program.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

C.............  Get sizes of inventory from the header

        NVAR     = GETIFDSC( FDESC3D, '/NON POLLUTANT/', .TRUE. )
        NIPOL    = GETIFDSC( FDESC3D, '/POLLUTANTS/', .FALSE. )
        NPPOL    = GETIFDSC( FDESC3D, '/PER POLLUTANT/', .FALSE. )
        NIACT    = GETIFDSC( FDESC3D, '/ACTIVITIES/', .FALSE. )
        NPACT    = GETIFDSC( FDESC3D, '/PER ACTIVITY/', .FALSE. )

        NIPOL = MAX( 0, NIPOL )
        NPPOL = MAX( 0, NPPOL )
        NIACT = MAX( 0, NIACT )
        NPACT = MAX( 0, NPACT )
        NIPPA = NIPOL + NIACT

C............. Retrieve base year information
        BYEAR = GETIFDSC( FDESC3D, '/BASE YEAR/', .FALSE. )

C............. Allocate memory and store the source attributes units
        IF( .NOT. ALLOCATED( ATTRUNIT ) ) THEN

            ALLOCATE( ATTRUNIT( NVAR ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ATTRUNIT', PROGNAME )

        END IF

        DO I = 1, NVAR
            ATTRUNIT( I ) = VUNITSET( I )
        END DO

        IF( .NOT. ALLOCATED( EANAM ) ) THEN

C............. Allocate memory for the array that stores pollutant
C............. names and activity names, units, and descriptions. Populate 
C              this array in the loops below that fill EINAM and NIACT
            ALLOCATE( EANAM( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
            ALLOCATE( EAREAD( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EAREAD', PROGNAME )
            ALLOCATE( EAUNIT( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EAUNIT', PROGNAME )
            ALLOCATE( EADESC( NIPPA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EADESC', PROGNAME )
            EANAM  = ' '  ! array
            EAREAD = ' '  ! array
            EAUNIT = ' '  ! array
            EADESC = ' '  ! array            

C.............  Allocate memory for and store pollutant names 
            ALLOCATE( EINAM( NIPOL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EINAM', PROGNAME )

        END IF

C.........  Ensure that the number of mapped variables is
C           consistent with the inventory file header
        IF ( NMAP .GT. 0 .AND. NMAP .LT. NIPPA ) THEN
            WRITE( MESG,94010 ) 'WARNING: Number of map-formatted ' //
     &             'variables is less than number in header of '//
     &             'inventory file.' // CRLF()// BLANK10 // 
     &             'I/O API header has', NIPPA, 'variables indicated'//
     &              ', but map-formated inventory file has', NMAP
            CALL M3MSG2( MESG )

            MESG = 'This occurs when the inventory contained ' //
     &             'pollutants with all zero values.' // CRLF() // 
     &             BLANK5//'You should confirm this in your '//
     &             'Smkinven log file.'
            CALL M3MSG2( MESG )
            NIPOL = NIPOL - ( NIPPA - NMAP )
            NIPPA = NIPOL + NIACT

        ELSE IF ( NMAP .GT. 0 .AND. NMAP .GT. NIPPA ) THEN
            WRITE( MESG,94010 ) 'Number of map-formatted ' //
     &             'variables is greater than number in header of '//
     &             'inventory file.' // CRLF()// BLANK10 // 
     &             'I/O API header has', NIPPA, 'variables indicated'//
     &              ', but map-formated inventory file has', NMAP
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )


        END IF

        K = 0
        J = NVAR + 1
        
        DO I = 1, NIPOL

            K = K + 1

C............  For map-formatted inventory
            IF( NMAP .GT. 0 ) THEN
                TMPNAME = 'TMP_POL_FILE'
                IF( .NOT. SETENVVAR( TMPNAME, MAPFIL( K ) ) ) THEN
                    MESG = 'ERROR: Could not set logical file name '
     &                          // CRLF() // BLANK10 // 'for file ' //
     &                          TRIM( MAPFIL( K ) )
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

                ELSE IF( .NOT. OPENSET( TMPNAME,FSREAD3,PROGNAME )) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not open "'// TRIM(MAPNAM(I))//
     &                      '" file listed in ' // CRLF() // BLANK10 // 
     &                      'map-formatted inventory for file name: ' //
     &                      CRLF() // BLANK10 // TRIM( MAPFIL( K ) )
                    CALL M3MSG2( MESG )
                    CYCLE

                ELSE IF( .NOT. DESCSET( TMPNAME,ALLFILES ) ) THEN
                    MESG = 'Could not get description of file "' //
     &                     TRIM( TMPNAME ) // '".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C................  Ensure that map variable name and file variable 
C                  name are the same
                IF( MAPNAM( K ) .NE. VNAMESET( 2 ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Variable name "'// TRIM(MAPNAM(I))//
     &                    '" in map-formatted inventory is mapped to '//
     &                     CRLF()// BLANK10// 'a file with variable "'//
     &                     TRIM( VNAMESET( 2 ) )//'". Map-formatted '//
     &                     'inventory has been corrupted.'
                    CALL M3MSG2( MESG )
                END IF

                EINAM ( I ) = VNAMESET( 2 )
                EANAM ( K ) = VNAMESET( 2 )
                EAREAD( K ) = VNAMESET( 2 + INVPIDX )
                EAUNIT( K ) = VUNITSET( 2 + INVPIDX )
                EADESC( K ) = VDESCSET( 2 + INVPIDX )

                IF( .NOT. CLOSESET( TMPNAME ) ) THEN
                    MESG = 'Could not close file:'//CRLF()//BLANK10//
     &                     TRIM( MAPFIL( K ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C............  For old style inventory
            ELSE

                EINAM ( I ) = VNAMESET( J )
                EANAM ( K ) = VNAMESET( J )
                EAREAD( K ) = VNAMESET( J + INVPIDX )
                EAUNIT( K ) = VUNITSET( J + INVPIDX )
                EADESC( K ) = VDESCSET( J + INVPIDX )

                J = J + NPPOL   ! skip over other pollutant-spec variables

            END IF

        END DO

C.........  Allocate memory for and store activity names
        IF( .NOT. ALLOCATED( ACTVTY ) ) THEN
            ALLOCATE( ACTVTY( NIACT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ACTVTY', PROGNAME )
        ENDIF

        DO I = 1, NIACT

            K = K + 1

C............  For map-formatted inventory
            IF( NMAP .GT. 0 ) THEN
                TMPNAME = 'TMP_POL_FILE'
                IF( .NOT. SETENVVAR( TMPNAME, MAPFIL( K ) ) ) THEN
                    MESG = 'ERROR: Could not set logical file name '
     &                          // CRLF() // BLANK10 // 'for file ' //
     &                          TRIM( MAPFIL( K ) )
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

                ELSE IF( .NOT. OPENSET( TMPNAME,FSREAD3,PROGNAME )) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Could not open "'// TRIM(MAPNAM(I))//
     &                      '" file listed in ' // CRLF() // BLANK10 // 
     &                      'map-formatted inventory for file name: ' //
     &                      CRLF() // BLANK10 // TRIM( MAPFIL( K ) )
                    CALL M3MSG2( MESG )
                    CYCLE

                ELSE IF( .NOT. DESCSET( TMPNAME,ALLFILES ) ) THEN
                    MESG = 'Could not get description of file "' //
     &                     TRIM( TMPNAME ) // '".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

C................  Ensure that map variable name and file variable 
C                  name are the same
                IF( MAPNAM( K ) .NE. VNAMESET( 2 ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Variable name "'// TRIM(MAPNAM(I))//
     &                    '" in map-formatted inventory is mapped to '//
     &                     CRLF()// BLANK10// 'a file with variable "'//
     &                     TRIM( VNAMESET( 2 ) )//'". Map-formatted '//
     &                     'inventory has been corrupted.'
                    CALL M3MSG2( MESG )
                END IF

                ACTVTY( I ) = VNAMESET( 2 )
                EANAM ( K ) = VNAMESET( 2 )
                EAREAD( K ) = VNAMESET( 2 + INVPIDX )
                EAUNIT( K ) = VUNITSET( 2 + INVPIDX )
                EADESC( K ) = VDESCSET( 2 + INVPIDX )

                IF( .NOT. CLOSESET( TMPNAME ) ) THEN
                    MESG = 'Could not close file:'//CRLF()//BLANK10//
     &                     TRIM( MAPFIL( K ) )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

C............  For old style inventory
            ELSE
                ACTVTY( I ) = VNAMESET( J )
                EANAM ( K ) = VNAMESET( J )
                EAREAD( K ) = VNAMESET( J )
                EAUNIT( K ) = VUNITSET( J )
                EADESC( K ) = VDESCSET( J )
                J = J + NPACT   ! skip over other activity-spec variables

            END IF

        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem with pollutant files'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Get header description of inventory file, so that
C           the I/O API header settings passed by the include files
C           will be set with this header.  Give error if problem.
        IF( .NOT. DESCSET( ENAME, ALLFILES ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx
 
94010   FORMAT( 10( A, :, I6, :, 2X ) )

        END SUBROUTINE GETSINFO
