
        SUBROUTINE OPENTMP( II, ENAME, SDATE, STIME, TSTEP, NSTEPS,
     &                      TZONE, NPELV, TNAME, PDEV, PFLAG )

C***********************************************************************
C  subroutine body starts at line 103
C
C  DESCRIPTION:
C      This subroutine opens the file or files for output from the tmppoint
C      program
C 
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 1/99 by M. Houyoux
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
C***************************************************************************

C........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: ACTVTY, BYEAR, CATDESC, CRL, EANAM, EADESC, 
     &                     EAREAD, EAUNIT, EINAM, INVPIDX, NIACT, 
     &                     NIPPA, NIPOL

C.........  This module contains the temporal profile tables
C        USE MODTMPRL, ONLY: STDATE

C.........  This module is required by the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(2)       CRLF
        INTEGER            INDEX1
        INTEGER            IOAPI_GRD_SIZE
        CHARACTER(IODLEN3) GETCFDSC
        INTEGER            GETIFDSC
        CHARACTER(IOULEN3) MULTUNIT
        INTEGER            PROMPTFFILE
        CHARACTER(16)      VERCHAR
        LOGICAL            SETENVVAR

        EXTERNAL        CRLF, INDEX1, IOAPI_GRD_SIZE, GETCFDSC, 
     &                  GETIFDSC, MULTUNIT, VERCHAR, SETENVVAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: II     ! episode time preriod index 
        CHARACTER(*), INTENT (IN) :: ENAME  ! emissions inven logical name
        INTEGER     , INTENT (IN) :: SDATE  ! episode start date 
        INTEGER     , INTENT (IN) :: STIME  ! episode start time
        INTEGER     , INTENT (IN) :: TSTEP  ! episode time step
        INTEGER     , INTENT (IN) :: NSTEPS ! number of time steps
        INTEGER     , INTENT (IN) :: TZONE  ! zone used for hours in output files
        INTEGER     , INTENT (IN) :: NPELV  ! number of elevated sources
        CHARACTER(*), INTENT(OUT) :: TNAME  ! lay-1 (or all) hourly logical name 
        INTEGER     , INTENT(OUT) :: PDEV   ! unit number of temporal supmtl file
        LOGICAL     , INTENT (IN) :: PFLAG  ! true: episode time periods needed

C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.8_Jun2020$'  ! CVS revision tag

C...........   Other local variables

        INTEGER         I, J, K, V     ! counters and indices

        INTEGER         IOS         ! i/o status
        INTEGER         FILESIZE    ! approximate size of emission factors file
        INTEGER         NINVVAR     ! number of inventory variables
        INTEGER         NVARFILE    ! number of variables per file
        INTEGER         PYEAR       ! projected year from inventory file (or -1)

        CHARACTER(5)    CTZONE      ! string of time zone
        CHARACTER(300)  MESG        ! message buffer 

        CHARACTER(NAMLEN3)  NAMBUF       ! [A|M|P]TMP file name buffer
        CHARACTER(512)  NAMBUFT          ! [A|M|P]TMPNAME file name buffer
        CHARACTER(512)  NAMBUFS          ! [A|M|P]TSUPNAME file name buffer
        CHARACTER(IODLEN3)  IFDESC2, IFDESC3 ! fields 2 & 3 from PNTS FDESC

        CHARACTER(16) :: PROGNAME = 'OPENTMP' ! program name

C***********************************************************************
C   begin body of subroutine OPENTMP

C.........  Write time zone to character string
        WRITE( CTZONE,94000 ) TZONE
 
C.........  Set up file header(s) for opening I/O API output(s). Base this on
C           inventory header...

C.........  Get header information from inventory file
        IF( .NOT. DESCSET( ENAME,-1 ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )
        NINVVAR = NVARSET
        BYEAR   = GETIFDSC( FDESC3D, '/BASE YEAR/', .TRUE. )
        PYEAR   = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )

        NVARSET = NIPPA
        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = TSTEP

        FDESC3D = ' '   ! array

        FDESC3D( 1 ) = CATDESC // ' source hourly emissions data'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/TZONE/ '   // CTZONE
        WRITE( FDESC3D( 5 ),94010 ) '/BASE YEAR/ ', BYEAR 
        IF( PYEAR .GT. 0 ) 
     &      WRITE( FDESC3D( 6 ),94010 ) '/PROJECTED YEAR/ ', PYEAR
        WRITE( FDESC3D( 7 ),94010 ) '/AVERAGE DAY/', INVPIDX

        FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Allocate memory for output arrays
        DEALLOCATE( VNAMESET, VUNITSET, VTYPESET, VDESCSET )
        DEALLOCATE( VARS_PER_FILE )
        ALLOCATE( VNAMESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNAMESET', PROGNAME )
        ALLOCATE( VUNITSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VUNITSET', PROGNAME )
        ALLOCATE( VTYPESET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VTYPESET', PROGNAME )
        ALLOCATE( VDESCSET( NVARSET ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VDESCSET', PROGNAME )

C.........  Check file size and adjust number of files to avoid 2 GB limit
        NFILESET = 1
        DO
            NVARFILE = ( NVARSET + NFILESET - 1 ) / NFILESET
            FILESIZE = IOAPI_GRD_SIZE( NCOLS3D, NROWS3D, NLAYS3D, 
     &                                 NVARFILE, NSTEPS )
            
            IF( FILESIZE > 1500 ) THEN
                NFILESET = NFILESET + 1
            ELSE
                EXIT
            END IF
        END DO

        IF( NFILESET > 1 ) THEN
            ALLOCATE( VARS_PER_FILE( NFILESET ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VARS_PER_FILE', PROGNAME )
            
            DO I = 1, NFILESET - 1
                VARS_PER_FILE( I ) = NVARFILE
            END DO

            VARS_PER_FILE( NFILESET ) = 
     &            NVARSET - ( NVARFILE*( NFILESET - 1 ) )
        END IF

C.........  Set variable names and characteristics from the pollutants
        K = 0
        DO V = 1, NIPOL

C.............  Double check that pollutant is in the inventory file
            I = INDEX1( EINAM( V ), NIPPA, EANAM )
            IF( I .LE. 0 ) THEN
                MESG='INTERNAL ERROR: inventory file variables changed!'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF

            K = K + 1
            IF( K .GT. NVARSET ) THEN
                MESG = 'INTERNAL ERROR: Memory overflow building '//
     &                 'I/O API output variables'
                CALL M3MSG2( MESG )
                CYCLE
            END IF

            VNAMESET( K ) = EINAM ( V )
            VUNITSET( K ) = EAUNIT( I )
            VDESCSET( K ) = EADESC( I )
            VTYPESET( K ) = M3REAL

        END DO  ! End loop on pollutants for output

C.........  Prompt for and open I/O API output file(s)...
        
        CALL GETENV( CRL // 'TMPNAME', NAMBUFT )
        WRITE( NAMBUFT, '( A, I7, A )' ) TRIM( NAMBUFT ), SDATE, '.ncf'

C..........  Set logical file name
        IF( .NOT. SETENVVAR( CRL // 'TMP', NAMBUFT ) ) THEN
            MESG = 'Could not set logical file name for ' //
     &             'file ' // TRIM( NAMBUFT )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MESG = 'Enter name for output HOURLY EMISSIONS file'
        NAMBUF = PROMPTSET( MESG, FSUNKN3, CRL // 'TMP', PROGNAME ) 
        TNAME = NAMBUF
     
C.........  Open supplemental speciation file
        CALL GETENV( CRL // 'TSUPNAME', NAMBUFS )
        WRITE( NAMBUFS, '( A, I7, A )' ) TRIM( NAMBUFS ), SDATE, '.txt'
        
C..........  Set logical file name
        IF( .NOT. SETENVVAR( CRL // 'TSUP', NAMBUFS ) ) THEN
            MESG = 'Could not set logical file name for ' //
     &             'file ' // TRIM( NAMBUFS )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MESG = 'Enter logical name for the TEMPORAL SUPPLEMENTAL '//
     &         'file'
        PDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 
     &                      CRL // 'TSUP', PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I5 )
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( A, :, I8, :, A )

        END SUBROUTINE OPENTMP

