
        SUBROUTINE OPENTMP_U( ENAME, UINAME,SDATE, STIME, TSTEP, TZONE,
     &                        RLZN, MAXLEN, UNAME, UDEV )

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
C      Created 2/02 by G. Cano
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$ $
C
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C Pathname: $ $
C Last updated: $ $ 
C
C***************************************************************************

C...........   MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

C.........  This module contains the information about the source category
        USE MODINFO

C...........   This module contains the uncertainty arrays and variables
        USE MODUNCERT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2            CRLF
        INTEGER                INDEX1
        CHARACTER(LEN=IODLEN3) GETCFDSC
        INTEGER                GETIFDSC
        CHARACTER(LEN=IOULEN3) MULTUNIT
        INTEGER                PROMPTFFILE
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL        CRLF, INDEX1, GETCFDSC, GETIFDSC, MULTUNIT,
     &                  PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME  ! emissions inven logical name
        CHARACTER(*), INTENT (IN) :: UINAME ! emissions uncertainty inven logical name
        INTEGER     , INTENT (IN) :: SDATE  ! episode start date 
        INTEGER     , INTENT (IN) :: STIME  ! episode start time
        INTEGER     , INTENT (IN) :: TSTEP  ! episode time step
        INTEGER     , INTENT (IN) :: TZONE  ! zone used for hours in output files
        INTEGER     , INTENT (IN) :: RLZN   ! realization number
        INTEGER     , INTENT (IN) :: MAXLEN ! the length of the number of realizations 
        CHARACTER(*), INTENT(OUT) :: UNAME  ! lay-1 (or all) hourly logical name 
        INTEGER     , INTENT(OUT) :: UDEV   ! unit number of temporal supmtl file

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$'  ! CVS revision tag
        CHARACTER*300              TMPBUF              ! temporay name label

C...........   Other local variables

        INTEGER         I, J, K, V               ! counters and indices

        INTEGER                      NINVVAR     ! number of inventory variables
        INTEGER                      PYEAR       ! projected year from inventory file (or -1)
        INTEGER, SAVE             :: NLEN = 0    ! temporary string length

        CHARACTER*5     CTZONE                   ! string of time zone
        CHARACTER*300   MESG                     ! message buffer 

        CHARACTER(LEN=NAMLEN3)  NAMBUF           ! file name buffer
        CHARACTER(LEN=IODLEN3)  IFDESC2, IFDESC3 ! fields 2 & 3 from PNTS FDESC

        CHARACTER*16 :: PROGNAME = 'OPENTMP_U' ! program name

C***********************************************************************
C   begin body of subroutine OPENTMP

C.........  Write time zone to character string
        WRITE( CTZONE,94000 ) TZONE
 
C.........  Set up file header(s) for opening I/O API output(s). Base this on
C           inventory header...

c BUG may need to add uncertainty file info gets
C.........  Get header information from inventory file
c        IF( .NOT. DESC3( UINAME ) ) THEN
c            MESG = 'Could not get description of file "' 
c     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
c            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
c        END IF

C.........  Get header information from inventory file
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )
        NINVVAR = NVARS3D
        BYEAR   = GETIFDSC( FDESC3D, '/BASE YEAR/', .TRUE. )
        PYEAR   = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )

        NROWS3D = UNSRC
        NVARS3D = UNIPOL
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
	WRITE( FDESC3D( 7 ),94010 ) '/OZONE SEASON/', INVPIDX

        FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Set variable names and characteristics from the emission types
        K = 0


c        DO J = 1, NIACT
C.............  Double check that pollutant is in the inventory file
C.............  Use EAREAD, because for mobile sources, EANAM has been
C               expanded to contain the emission types
c            I = INDEX1( ACTVTY( J ), NIACT + NIPOL, EAREAD )
c            IF( I .LE. 0 ) THEN
c                MESG='INTERNAL ERROR: inventory file variables changed!'
c                CALL M3MSG2( MESG )
c                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
c            END IF

c            DO V = 1, NETYPE( I-NIPOL )

c        	K = K + 1
c       	VNAME3D( K ) = EMTNAM( V,J )
c        	UNITS3D( K ) = EAUNIT( I )   
c        	VDESC3D( K ) = EMTDSC( V,J )
c        	VTYPE3D( K ) = M3REAL

c            END DO  ! End loop on emission types for output
c        END DO      ! End loop on activities output

C.........  Set variable names and characteristics from the pollutants
        DO V = 1, NIPOL

C.............  Double check that uncertainty pollutant is in the inventory file
            I = INDEX1( EINAM( V ), NIPPA, EANAM )
            IF( .NOT. ( UEINAM( I ) ) ) CYCLE

            K = K + 1
            VNAME3D( K ) = EINAM ( V )
            UNITS3D( K ) = EAUNIT( I )
            VDESC3D( K ) = EADESC( I )
            VTYPE3D( K ) = M3REAL

        END DO  ! End loop on pollutants for output

C.........  Set up name for output file
        WRITE( TMPBUF,94010 ) '', RLZN
        TMPBUF = ADJUSTL( TMPBUF )
        IF (NLEN .EQ. 0) NLEN = MAXLEN + 1
        CALL PADNZERO( NLEN, TMPBUF )
        WRITE( UNAME,94010 ) CRL // 'TMP_U' // TMPBUF(1:NLEN)

C.........  Prompt for and open I/O API output file(s)...
        MESG = 'Enter name for UNCERTAINTY output HOURLY EMISSIONS file'
        NAMBUF = PROMPTMFILE( MESG, FSUNKN3, UNAME, PROGNAME )
        UNAME = NAMBUF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE OPENTMP_U
