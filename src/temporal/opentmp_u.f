
        SUBROUTINE OPENTMP_U( ENAME, UINAME,SDATE, STIME, TSTEP, TZONE,
     &                        RLZN, UNAME, UDEV )

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
        LOGICAL                SETENVVAR

        EXTERNAL        CRLF, INDEX1, GETCFDSC, GETIFDSC, MULTUNIT,
     &                  PROMPTMFILE, VERCHAR, SETENVVAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME  ! emissions inven logical name
        CHARACTER(*), INTENT (IN) :: UINAME ! emissions uncertainty inven logical name
        INTEGER     , INTENT (IN) :: SDATE  ! episode start date 
        INTEGER     , INTENT (IN) :: STIME  ! episode start time
        INTEGER     , INTENT (IN) :: TSTEP  ! episode time step
        INTEGER     , INTENT (IN) :: TZONE  ! zone used for hours in output files
        INTEGER     , INTENT (IN) :: RLZN   ! realization number
        CHARACTER(*), INTENT(IN OUT) :: UNAME  ! lay-1 (or all) hourly logical name 
        INTEGER     , INTENT(OUT) :: UDEV   ! unit number of temporal supmtl file

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$'  ! CVS revision tag

C...........   Other local variables

        INTEGER         I, J, K, V               ! counters and indices

        INTEGER                      NINVVAR     ! number of inventory variables
        INTEGER                      PYEAR       ! projected year from inventory file (or -1)

        CHARACTER*5     CTZONE                   ! string of time zone
        CHARACTER*16    FNAM                     ! file name buffer 
        CHARACTER*300   DIRBUF                   ! directory buffer 
        CHARACTER*300   FNAMX                    ! file name extension 
        CHARACTER*300   LNAME                    ! full path and name for uncertainty files
        CHARACTER*300   MESG                     ! message buffer 

        CHARACTER(LEN=NAMLEN3)  NAMBUF           ! file name buffer
        CHARACTER(LEN=IODLEN3)  IFDESC2, IFDESC3 ! fields 2 & 3 from PNTS FDESC

        CHARACTER*16 :: PROGNAME = 'OPENTMP_U' ! program name

C***********************************************************************
C   begin body of subroutine OPENTMP

        IF( RLZN .EQ. 1) THEN

            UNAME = UCAT( 1:1 ) // 'TMPU'

        ELSE IF( RLZN .GT. 1) THEN
            IF( .NOT.( CLOSE3( UNAME ) ) ) THEN
                WRITE( MESG,94010 ) 
     &                 'WARNING: Unable to close file' // UNAME
                CALL M3MSG2( MESG )
            END IF
        END IF

C.........  Set up file header(s) for opening I/O API output(s). Base this on
C           inventory header...
        CALL GRFNAME
        IF( .NOT.( SETENVVAR( UNAME, LNAME ) ) )
     &       CALL M3EXIT( PROGNAME, 0, 0, 'Unable to assign setenv', 2 )

C.........  Write time zone to character string
        WRITE( CTZONE,94000 ) TZONE

C.........  Get header information from inventory file
        IF( .NOT. DESC3( UINAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // UINAME( 1:LEN_TRIM( UINAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

c        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
c        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )
c        NINVVAR = NVARS3D
c        BYEAR   = GETIFDSC( FDESC3D, '/BASE YEAR/', .TRUE. )
c        PYEAR   = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )

        NLAYS3D = 1
        NCOLS3D = 1
        NROWS3D = UNSRC
        NVARS3D = UNIPOL
        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = TSTEP
        FTYPE3D = GRDDED3
c        GDNAM3D = ' '
c        FDESC3D = ' '   ! array
 
        FDESC3D( 1 ) = CATDESC // ' source hourly emissions data'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4 ) = '/TZONE/ '   // CTZONE

C.........  Set variable names and characteristics from the emission types
        K = 0

C.........  Set variable names and characteristics from the pollutants
        DO V = 1, NIPOL

C.............  Double check that uncertainty pollutant is in the inventory file
            I = INDEX1( EINAM( V ), NIPPA, EANAM )
            IF( I .LE. 0 ) THEN
                MESG='INTERNAL ERROR: inventory file variables changed!'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF
            IF( .NOT. ( UEINAM( I ) ) ) CYCLE

            K = K + 1
            VNAME3D( K ) = EINAM ( V )
            UNITS3D( K ) = EAUNIT( I )
            VDESC3D( K ) = EADESC( I )
            VTYPE3D( K ) = M3REAL

        END DO  ! End loop on pollutants for output

C.........  Prompt for and open I/O API output file(s)...
        MESG = 'Enter name for UNCERTAINTY output HOURLY EMISSIONS file'
        NAMBUF = PROMPTMFILE( MESG, FSUNKN3, 
     &                        UCAT( 1:1 ) // 'TMPU', PROGNAME )
        UNAME = NAMBUF

C.........  Open supplemental speciation file
c        MESG = 'Enter logical name for the TEMPORAL SUPPLEMENTAL '//
c     &         'file'
c        PDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 
c     &                      CRL // 'TSUP', PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )
 
94010   FORMAT( 10( A, :, I8, :, 1X ) )

        CONTAINS

            SUBROUTINE GRFNAME

            INTEGER         L1, L2, L3          ! temporary numbers

            CHARACTER*300   TMPBUF              ! temporay name label
C**********************************************************************

C.........  Set up name for output file
        CALL NAMEVAL( UCAT( 1:1 ) // 'TMPU_ODIR', DIRBUF )
        L1 = LEN_TRIM( DIRBUF )
        CALL NAMEVAL( UCAT( 1:1 ) // 'TMPU_FNAM', FNAM )
        L2 = LEN_TRIM( FNAM )

        WRITE( TMPBUF,94010 ) '', RLZN
        TMPBUF = ADJUSTL( TMPBUF )
        CALL PADNZERO( RMXLEN, TMPBUF )
        L3 = LEN_TRIM( TMPBUF )

        LNAME = DIRBUF( 1:L1 ) // '/' // 
     &          FNAM( 1: L2 )  // '_' //
     &          TMPBUF( 1: L3 )  // '.ncf' 


C...........   Internal buffering formats............ 94xxx
94010   FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE GRFNAME

        END SUBROUTINE OPENTMP_U
