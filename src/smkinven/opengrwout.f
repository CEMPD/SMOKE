
        SUBROUTINE OPENGRWOUT( ENAME, FYEAR, ONAME, DDEV, VDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 2/2000 by M. Houyoux
C
C***************************************************************************
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
        LOGICAL                 ENVYN
        CHARACTER(LEN=IODLEN3)  GETCFDSC
        INTEGER                 INDEX1
        INTEGER                 PROMPTFFILE
        CHARACTER(LEN=NAMLEN3)  PROMPTMFILE
        CHARACTER*16            VERCHAR

        EXTERNAL CRLF, ENVYN, GETCFDSC, INDEX1, PROMPTFFILE, 
     &           PROMPTMFILE, VERCHAR

C...........   SUBROUTINE ARGUMENTS
        CHARACTER*16, INTENT (IN) :: ENAME ! emis input inven logical name
        INTEGER     , INTENT (IN) :: FYEAR ! future year or 0 for no projection
        CHARACTER*16, INTENT(OUT) :: ONAME ! emis output inven logical name
        INTEGER     , INTENT(OUT) :: DDEV  ! IDA output emissions file number
        INTEGER     , INTENT(OUT) :: VDEV  ! IDA output acitivty file number

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Other local variables

        INTEGER         CATIDX    ! index for source category
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
        IFLAG = ENVYN( EVNAME, MESG, .FALSE., IOS )

C.........  Initialize outputs
        ONAME = 'NONE'
        DDEV  = 0
        VDEV  = 0

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
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

C.........  If projection year is non-zero, store in FDESC3D
        IF( FYEAR .NE. 0 ) THEN
             WRITE( FDESC3D( 14 ),94010 ) '/PROJECTED YEAR/ ', FYEAR
        END IF

        FDESC3D( 15 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 16 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Prompt file output SMOKE file
        IF( SFLAG ) THEN
            MESG= 'Enter logical name for the I/O API INVENTORY ' //
     &            'output file'
            L = LEN_TRIM( ENAME )
            ONAME = ENAME( 1:L ) // '_O'
            ONAME = PROMPTMFILE( MESG, FSUNKN3, ONAME, PROGNAME )
        END IF

C.........  Get index for source category to use for output file names
        CATIDX = INDEX1( CATEGORY, NCAT, CATLIST )

C.........  Prompt file emissions IDA file
        IF( IFLAG .AND. NIPOL .GT. 0 ) THEN
            MESG  = 'Enter logical name for the IDA EMISSIONS ' //
     &              'output file'
            INAME = ANAMLIST( CATIDX ) // '_O'

            DDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., INAME, PROGNAME )
        END IF

C.........  Prompt file activity IDA file
        IF( IFLAG .AND. NIACT .GT. 0 ) THEN
            MESG  = 'Enter logical name for the IDA ACTIVITY ' //
     &              'output file'
            INAME = ANAMLIST( CATIDX ) // '_AO'

            VDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., INAME, PROGNAME )
        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE OPENGRWOUT

