
        SUBROUTINE OPENGRWOUT( ENAME, FYEAR, ONAME, FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine 
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 2/2000 by M. Houyoux
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
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
C Pathname: $Source$
C Last updated: $Date$ 
C
C***************************************************************************

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        CHARACTER*2             CRLF
        INTEGER                 GETIFDSC
        INTEGER                 PROMPTFFILE
        CHARACTER(LEN=NAMLEN3)  PROMPTMFILE
        CHARACTER*16            VERCHAR

        EXTERNAL CRLF, PROMPTFFILE, PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME ! emis input inven logical name
        INTEGER     , INTENT (IN) :: FYEAR ! future year or 0 for no projection
        CHARACTER(*), INTENT(OUT) :: ONAME ! emis output inven logical name
        INTEGER     , INTENT(OUT) :: FDEV  ! IDA output unit number

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: SCCSW  = '@(#)$Id$'  ! SCCS string with vers no.

C...........   Other local variables

        INTEGER         IOS       ! i/o status
        INTEGER         L         ! counters and indices

        LOGICAL      :: IFLAG = .TRUE.   ! true: output IDA file
        LOGICAL      :: SFLAG = .TRUE.   ! true: output SMOKE file
 
        CHARACTER*16    INAME     ! tmp output IDA file name
        CHARACTER*30    EVNAME    ! tmp environment variable name
        CHARACTER*300   MESG      ! message buffer 
        CHARACTER(LEN=IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from inven FDESC

        CHARACTER*16 :: PROGNAME = 'OPENGRWOUT' ! program name

C***********************************************************************
C   begin body of subroutine OPENGRWOUT

C.........  Evaluate environment variables that control the output files
C           used
        EVNAME = 'SMK_GRWSMKOUT_YN'
        MESG   = 'Output SMOKE inventory file'
        SFLAG = ENVYN( EVNAME, MESG, .TRUE., IOS )

        EVNAME = 'SMK_GRWIDAOUT_YN'
        MESG   = 'Output IDA inventory file'
        IFLAG = ENVYN( EVNAME, MESG, .TRUE., IOS )

C.........  Initialize outputs
        ONAME = 'NONE'
        FDEV  = 0

C.........  Re-read header for input inventory file
        IF( .NOT. DESC3( ENAME ) ) THEN
            L = LEN_TRIM( ENAME )
            MESG = 'Could not read description for "' //
     &             ENAME( 1:L ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Store non-category-specific header information
        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )

C.........  If projection year is non-zero, store in FDESC3D
        IF( FYEAR .NE. 0 ) THEN
             WRITE( FDESC3D( 12 ),94010 ) '/PROJECTED YEAR/ ', FYEAR
        END IF

        FDESC3D( 13 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 14 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Prompt file output SMOKE file (or NONE)
        MESG= 'Enter logical name for the I/O API INVENTORY output file'
        ONAME = ENAME // '_OUT'
        ONAME = PROMPTMFILE( MESG, FSUNKN3, ONAME, PROGNAME )

C.........  Prompt file output IDA file (or NONE)
        MESG  = 'Enter logical name for the IDA INVENTORY output file'
        INAME = CRL // 'IDA_OUT'

        DDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., INAME, PROGNAME )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENGRWOUT

