
        SUBROUTINE OPENSEF( SRCCT, DESC, SDATE, EDATE, EMISDIR, FNAME )

C***********************************************************************
C  subroutine body starts at line 89
C
C  DESCRIPTION:
C       Opens output emission factor files. Creates file name based on 
C       SMK_EMISPATH directory, current averaging type, and current start
C       date. Stores end date in file header.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C     10/01: Created by C. Seppanen
C
C***********************************************************************
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
C***********************************************************************

C...........   MODULES for public variables
C...........   This module contains the information about the source category
        USE MODINFO

C.........  This module contains emission factor tables and related
        USE MODEMFAC
        
        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures
        INCLUDE 'M6CNST3.EXT'   !  Mobile6 constants

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(LEN=IODLEN3) GETCFDSC
        INTEGER                INDEX1
        CHARACTER*16           VERCHAR
        LOGICAL                SETENVVAR

        EXTERNAL     GETCFDSC, INDEX1, VERCHAR, SETENVVAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT    (IN) :: SRCCT    ! total number of sources
        CHARACTER(*), INTENT    (IN) :: DESC     ! description of file type
        INTEGER     , INTENT    (IN) :: SDATE    ! julian start date
        INTEGER     , INTENT    (IN) :: EDATE    ! julian end date
        CHARACTER(*), INTENT    (IN) :: EMISDIR  ! directory of output files
        CHARACTER(*), INTENT(IN OUT) :: FNAME    ! name output emission factors file

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Other local variables
        INTEGER     I, J, K, L, LJ

        CHARACTER*300   MESG      ! message buffer
        CHARACTER*256   FULLNAME  ! full file name
        CHARACTER*16    CURRVNAME ! current variable name
        CHARACTER*16    POLNAME   ! current pollutant name

        CHARACTER*16 :: PROGNAME = 'OPENSEF' ! program name

C***********************************************************************
C   begin body of subroutine OPENSEF

C.........  Initialize I/O API output file headers
        CALL HDRMISS3

        FDESC3D( 1 ) = CATEGORY( 1:LEN_TRIM( CATEGORY ) ) //
     &                 ' emission factors file'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/NOTE/ Time 000000 in file represents ' //
     &                 '6 AM in local time zone'
        WRITE( FDESC3D( 5 ), 94010 ) '/END DATE/ ', EDATE

C.........  Set header values that cannot be default

        SDATE3D = SDATE
        STIME3D = 0
        TSTEP3D = 10000
        NVARS3D = MXETYPE + 1
        NROWS3D = SRCCT
        NLAYS3D = 1
 
        I = 1
        VNAME3D( I ) = 'SOURCES'
        UNITS3D( I ) = 'n/a'
        VDESC3D( I ) = 'Source number'
        VTYPE3D( I ) = M3INT
        
C.........  Loop through emission process/pollutant combos 
C           (EMTNAM is created from MEPROC file and contains only the ones we want, 
C           EFS* contains all possible MOBILE6 outputs with units and descriptions)        
        DO I = 1, SIZE( EMTNAM,1 )

            CURRVNAME = TRIM( EMTNAM( I,1 ) )

C.............  Stop when there are no variable names left            
            IF( CURRVNAME == ' ' ) EXIT
            
C.............  Make sure there aren't too many variables
            IF( I > MXVARS3 ) THEN
                MESG = 'Too many emission factors for single ' //
     &                 'I/O API file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            
C.............  Check to see if current pollutant is the output hydrocarbon
C               If so, change name back to input name to match M6 names
            K = 0
            IF( OUTPUTHC /= ' ' ) THEN
                K = INDEX( CURRVNAME, TRIM( OUTPUTHC ) )
                IF( K > 0 ) THEN
                    CURRVNAME = CURRVNAME( 1:K-1 ) // INPUTHC
                END IF
            END IF

C.............  Look for name in master list of emission factor names            
            L = INDEX1( CURRVNAME, NEFS, EFSNAM )

C.............  If this is a non-HAP HC name or user-defined HAP, 
C               use name from MEPROC file; otherwise use master list name            
            IF( K > 0 .OR. L <= 0 ) THEN
                VNAME3D( I+1 ) = EMTNAM( I,1 )
            ELSE
                VNAME3D( I+1 ) = EFSNAM( L )
            END IF

C.............  Use master list description and units for intrinsic pollutants            
            IF( L > 0 ) THEN
                UNITS3D( I+1 ) = EFSUNT( L )
                VDESC3D( I+1 ) = EFSDSC( L )
            ELSE
            	
C.................  Otherwise need to build info for user-defined HAPS
                UNITS3D( I+1 ) = M6UNIT

C.................  Extract pollutant name from current variable name
                L = INDEX( CURRVNAME, ETJOIN )
                LJ = LEN_TRIM( ETJOIN )
                POLNAME = CURRVNAME( L+LJ:LEN_TRIM( CURRVNAME ) )

C.................  Determine emission process and create description
                DO J = 1, MXM6EPR
                    L = INDEX( CURRVNAME, TRIM( M6PROCS( J ) ) )
                    IF( L > 0 ) THEN
                        VDESC3D( I+1 ) = 'EFs for ' // 
     &                      TRIM( M6PRCDSC( J ) ) // ' ' // 
     &                      POLNAME
                        EXIT
                    END IF
                END DO
            END IF
            
            VTYPE3D( I+1 ) = M3REAL
            
C.............  Append to description of hydrocarbon pollutant if necessary
            IF( K > 0 ) THEN
                VDESC3D( I+1 ) = TRIM( VDESC3D( I+1 ) ) // 
     &                           ' (minus HAPS)'
            END IF            
        END DO

C.........  Create full file name
        WRITE( FULLNAME,94010 ) EMISDIR( 1:LEN_TRIM( EMISDIR ) ) //
     &                          '/emisfacs.' // 
     &                          DESC( 1:LEN_TRIM( DESC ) ) // '.', 
     &                          SDATE, '.ncf'

C.........  Set logical file name
        IF( .NOT. SETENVVAR( FNAME, FULLNAME ) ) THEN
            MESG = 'Could not set logical file name for file ' //
     &             TRIM( FULLNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF       

C.........  Open new file
        IF( .NOT. OPEN3( FNAME, FSUNKN3, PROGNAME ) ) THEN
            MESG = 'Could not create new output file ' // 
     &             TRIM( FULLNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, I7, A )

94020   FORMAT( A, I7, A1, I7, A )

        END SUBROUTINE OPENSEF


