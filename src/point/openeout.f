
        SUBROUTINE OPENEOUT( NGROUP, SDATE, STIME, ENAME, PDEV, MNAME )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C       Opens the output files for the Elevpoint program
C       
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Create 8/99 by M Houyoux
C
C************************************************************************
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
C Last updated: %G 
C  
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        LOGICAL         DSCM3GRD
        CHARACTER*50    GETCFDSC
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        CHARACTER*16    VERCHAR

        EXTERNAL        CRLF, DSCM3GRD, GETCFDSC, PROMPTFFILE, 
     &                  PROMPTMFILE, VERCHAR

C..........    Subroutine arguments and their descriptions
        INTEGER     , INTENT (IN) :: NGROUP  ! number of ping groups
        INTEGER     , INTENT (IN) :: SDATE   ! start date of episode
        INTEGER     , INTENT (IN) :: STIME   ! start time of episode
        CHARACTER(*), INTENT (IN) :: ENAME   ! i/o api inventory file
        INTEGER     , INTENT (OUT):: PDEV    ! ASCII file for major/ping src IDs
        CHARACTER(*), INTENT (OUT):: MNAME   ! logical name of ping srcs groups

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: SCCSW  = '@(#)$Id$'  ! SCCS string with vers no.

C...........   Other local variables
        INTEGER         J      ! indices and counters
         
        LOGICAL      :: PINGFLAG = .FALSE.  !  true: there are ping sources
        LOGICAL      :: EFLAG    = .FALSE.  !  true: error found

        CHARACTER*80    GDESC               !  grid description
        CHARACTER*300   MESG

        CHARACTER(LEN=NAMLEN3) NAMBUF           ! file name buffer
        CHARACTER(LEN=IOVLEN3) COORD3D
        CHARACTER(LEN=IOVLEN3) COORUN3D
        CHARACTER(LEN=IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from inven FDESC

        CHARACTER*16 :: PROGNAME = 'OPENEOUT'   !  subroutine name

C***********************************************************************
C   begin body of subroutine OPENEOUT

C.........  Open ASCII output file
        PDEV = PROMPTFFILE(  
     &        'Enter name for ELEVATED POINT SOURCE output file',
     &        .FALSE., .TRUE., CRL // 'ELV', PROGNAME )

        PINGFLAG = ( NGROUP .GT. 0 )

C.........  If needed, set up and open plume-in-grid i/o api output file
        IF( PINGFLAG ) THEN

C.........  Get header information from inventory file

            IF ( .NOT. DESC3( ENAME ) ) THEN
        	MESG = 'Could not get description of file "' 
     &                 // ENAME( 1:LEN_TRIM( ENAME ) ) // '".'
        	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
            IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

C.............  Set up for opening I/O API output file header
            CALL HDRMISS3  ! Initialize for emissions 

C.............  Get grid description for setting most grid parameters for
C               the output file
            IF( .NOT. DSCM3GRD( 
     &                GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,
     &                P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

        	MESG = 'Could not get Models-3 grid description.'
        	CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )            
            END IF            

C............. Finalize i/o api header fields
C............. NOTE - this is a time-independent file, but the Plume Dynamics
C              Model that reads this file needs to have the start date and time
C              of the Met file in here.
            FDESC3D( 1 ) = CATDESC // ' source stack groups file'
            FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
            FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( SCCSW )
            WRITE( FDESC3D(5), 94010 ) '/NCOLS3D/ ', NCOLS3D
            WRITE( FDESC3D(6), 94010 ) '/NROWS3D/ ', NROWS3D
            FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
            FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

            NVARS3D = 13
            NROWS3D = NGROUP
            NCOLS3D = 1
            NLAYS3D = 1
            SDATE3D = SDATE
            STIME3D = STIME
            TSTEP3D = 10000

C.............  Set the file variables
            J = 1
            VNAME3D( J ) = 'ISTACK'
            VTYPE3D( J ) = M3INT
            UNITS3D( J ) = 'none'
            VDESC3D( J ) = 'Stack group number'
            J = J + 1

            VNAME3D( J ) = 'LATITUDE'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'Latitude'
            J = J + 1

            VNAME3D( J ) = 'LONGITUDE'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'Longitude'
            J = J + 1

            VNAME3D( J ) = 'STKDM'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'm'
            VDESC3D( J ) = 'Inside stack diameter'
            J = J + 1

            VNAME3D( J ) = 'STKHT'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'm'
            VDESC3D( J ) = 'Stack height above ground surface'
            J = J + 1

            VNAME3D( J ) = 'STKTK'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees K'
            VDESC3D( J ) = 'Stack exit temperature'
            J = J + 1

            VNAME3D( J ) = 'STKVE'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'm/s'
            VDESC3D( J ) = 'Stack exit velocity'
            J = J + 1

            VNAME3D( J ) = 'STKFLW'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'm**3/s'
            VDESC3D( J ) = 'Stack exit flow rate'
            J = J + 1

            VNAME3D( J ) = 'STKCNT'
            VTYPE3D( J ) = M3INT
            UNITS3D( J ) = 'none'
            VDESC3D( J ) = 'Number of stacks in group'
            J = J + 1

            VNAME3D( J ) = 'ROW'
            VTYPE3D( J ) = M3INT
            UNITS3D( J ) = 'none'
            VDESC3D( J ) = 'Grid row number'
            J = J + 1

            VNAME3D( J ) = 'COL'
            VTYPE3D( J ) = M3INT
            UNITS3D( J ) = 'none'
            VDESC3D( J ) = 'Grid column number'
            J = J + 1

            VNAME3D( J ) = 'XLOCA'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = COORUN3D
            VDESC3D( J ) = 'Projection x coordinate'
            J = J + 1

            VNAME3D( J ) = 'YLOCA'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = COORUN3D
            VDESC3D( J ) = 'Projection y coordinate'
 
            MESG = 'Enter logical name for ELEVATED STACK GROUPS ' //
     &              'file'

C.............  Open file. Use NAMBUF for HP.
            MNAME  = 'STACK_GROUPS'
            NAMBUF = PROMPTMFILE( MESG, FSUNKN3, MNAME, PROGNAME )
            MNAME  = NAMBUF

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( 4 I10, F10.2 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 2X  ) )

94015   FORMAT( A, 1X, I6, 1X, A, 1X, I5.5, 1X, A, 1X, I8.8, 1X,
     &          A, 1X, I6, 1X, A, 1X, I6,   1X, A )

94020   FORMAT( 5X, 'H[m]:', 1X, F6.2, 1X, 'D[m]:'  , 1X, F4.2, 1X,
     &              'T[K]:', 1X, F7.1, 1X, 'V[m/s]:', 1X, F10.1 )

        END SUBROUTINE OPENEOUT

