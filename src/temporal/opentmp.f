
        SUBROUTINE OPENPTMP( ENAME, UFLAG, SDATE, STIME, TSTEP, TZONE, 
     &                       NPELV, NIPOL, EINAM, TNAME, UNAME )

C***********************************************************************
C  subroutine body starts at line
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

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        INTEGER                INDEX1
        CHARACTER(LEN=IODLEN3) GETCFDSC
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE
        CHARACTER*16           VERCHAR

        EXTERNAL        INDEX1, GETCFDSC, PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(LEN=NAMLEN3) ENAME ! emissions inven logical name (in)
        LOGICAL       UFLAG          ! true if elevated output (in)
        INTEGER       SDATE          ! episode start date (in)
        INTEGER       STIME          ! episode start time (in)
        INTEGER       TSTEP          ! episode time step (in)
        INTEGER       TZONE          ! zone used for hours in output files (in)
        INTEGER       NPELV          ! number of elevated sources (in)
        INTEGER       NIPOL          ! number of inventory pollutants (in)
        CHARACTER(LEN=IOVLEN3) EINAM( NIPOL ) ! Name of actual pollutants (in)
        CHARACTER(LEN=NAMLEN3) TNAME ! lay-1 (or all) hourly logical name (out)
        CHARACTER(LEN=NAMLEN3) UNAME ! elevated hourly logical name (out)

C...........   LOCAL PARAMETERS
        CHARACTER*50  SCCSW          ! SCCS string with version number at end

        PARAMETER   ( SCCSW   = '@(#)$Id$'
     &              )

C...........   Other local variables

        INTEGER         I, J, V     ! counters and indices

        INTEGER         NINVVAR     ! number of inventory variables

        CHARACTER*5     CTZONE      ! string of time zone
        CHARACTER*300   MESG        ! message buffer 

        CHARACTER(LEN=IODLEN3)  IFDESC2, IFDESC3 ! fields 2 & 3 from PNTS FDESC

        CHARACTER*16 :: PROGNAME = 'OPENPTMP' ! program name

C***********************************************************************
C   begin body of subroutine OPENPTMP

C.........  Write time zone to character string
        WRITE( CTZONE,94000 ) TZONE
 
C.........  Set up file header(s) for opening I/O API output(s). Base this on
C           inventory header...

C.........  Get header information from inventory file

        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not get description of file "' 
     &             // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/' )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/' )
        NINVVAR = NVARS3D

        NVARS3D = NIPOL
        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = TSTEP

        DO I = 1, MXDESC3
            FDESC3D( I ) = ' '
        ENDDO

        FDESC3D( 1 ) = 'Point source hourly emissions data'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )
        FDESC3D( 4 ) = '/TZONE/ '   // CTZONE

        FDESC3D( 11 ) = '/PNTS FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/PNTS VERSION/ ' // IFDESC3

C.........  Define source characteristic variables that are not strings

        DO V = 1, NIPOL

            I = INDEX1( EINAM( V ), NINVVAR, VNAME3D )
            IF( I .LE. 0 ) THEN
                MESG='INTERNAL ERROR: inventory file variables changed!'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF

            J = INDEX( UNITS3D( I ), '/' )

            VNAME3D( V ) = EINAM( V )
            UNITS3D( V ) = UNITS3D( I )( 1:J ) // 'hour'
            VDESC3D( V ) = VDESC3D( I )
            VTYPE3D( V ) = M3REAL

        ENDDO  ! End loop on pollutants for output

C.........  Prompt for and open I/O API output file(s)...

C.........  If elevated sources are being written to a separate file, reset the
C           applicable header entries and prompt for both layer-1 and upper-
C           level point sources files
        IF( UFLAG ) THEN 
        
            FDESC3D( 1 ) = 'Point source layer-1 hourly emissions data'

            MESG = 'Enter name for output LAYER-1 HOURLY EMISSIONS file'
            TNAME= PROMPTMFILE( MESG, FSUNKN3, 'PTMP', PROGNAME )

            NVARS3D = NIPOL + 1
            NROWS3D = NPELV
            VNAME3D( NVARS3D ) = 'INDXE'
            UNITS3D( NVARS3D ) = 'n/a'
            VDESC3D( NVARS3D ) = 'SMOKE source ID number'
            VTYPE3D( NVARS3D ) = M3INT

            FDESC3D( 1 ) = 'Point source elevated hourly emissions data'

            MESG ='Enter name for output ELEVATED HOURLY EMISSIONS file'
            UNAME= PROMPTMFILE( MESG, FSUNKN3, 'ETMP', PROGNAME )

C.........  Without elevated sources, the header is already set, and just need
C           to prompt for one file.
        ELSE

            MESG = 'Enter name for output HOURLY EMISSIONS file'
            TNAME = PROMPTMFILE( MESG, FSUNKN3, 'PTMP', PROGNAME )

            UNAME = ' '

        ENDIF      !  if generating upper-level point sources file or not

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I2.2 )
 
        END SUBROUTINE OPENPTMP

