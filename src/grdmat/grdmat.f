
        PROGRAM GRDPMAT

C***********************************************************************
C  program body starts at line 
C
C  NOTE to mrh: should there be explicit prompting of the grid definition
C   file?
C
C  DESCRIPTION:
C     Creates the point source gridding matrix
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

C...........   MODULES for public variables
C...........   This module is the source inventory arrays
        USE MODSOURC

        IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        CHARACTER*2            CRLF
        LOGICAL                DSCM3GRD
        CHARACTER(LEN=IODLEN3) GETCFDSC
        CHARACTER*16           PROMPTMFILE
        INTEGER                TRIMLEN
        CHARACTER*16           VERCHAR
   
        EXTERNAL  CRLF, DSCM3GRD, GETCFDSC, PROMPTMFILE, 
     &            TRIMLEN, VERCHAR

C...........   LOCAL PARAMETERS
        CHARACTER*50  SCCSW          ! SCCS string with version number at end

        PARAMETER   ( SCCSW   = '@(#)$Id$'
     &              )

C...........   LOCAL VARIABLES and their descriptions:

C...........   Gridding Matrix

        INTEGER, ALLOCATABLE :: GMAT( : ) ! Contiguous gridding matrix

C...........   Indicator for which public inventory arrays need to be read
        INTEGER               , PARAMETER :: NINVARR = 3
        CHARACTER(LEN=IOVLEN3), PARAMETER :: IVARNAMS( NINVARR ) = 
     &                                 ( / 'IFIP           '
     &                                   , 'XLOCA          '
     &                                   , 'YLOCA          ' / )

C...........   File units and logical/physical names
        INTEGER         LDEV    !  log-device
        CHARACTER*16    ENAME   !  logical name for point inventory input file
        CHARACTER*16    GNAME   !  logical name for grid matrix output file

C...........   Other local variables
        
        INTEGER         L1, K     !  indices and counters.

        INTEGER         CMAX  ! max number srcs per cell
        INTEGER         CMIN  ! min number srcs per cell
        INTEGER         IOS   ! i/o status
        INTEGER         NK    ! Number of coeficients
        INTEGER         NPSRC ! No of point sources
        INTEGER         NGRID ! No of cells
                                 
        REAL            CAVG  ! average number sources per cell

        CHARACTER*16            COORD    !  coordinate system name
        CHARACTER*16            COORUNIT !  coordinate system projection units
        CHARACTER*16            GRDNM    !  grid name
        CHARACTER*80            GDESC    !  grid description
        CHARACTER*300           MESG     !  message buffer

        CHARACTER(LEN=IODLEN3)  IFDESC2, IFDESC3 !  fields 2 & 3 from PNTS FDESC

        CHARACTER*16 :: PROGNAME = 'GRDPMAT'   !  program name

C***********************************************************************
C   begin body of program GRDPMAT
        
        LDEV = INIT3()
        
C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, SCCSW, PROGNAME )

C.........  Get file name; open inventory file
        ENAME = PROMPTMFILE( 
     &        'Enter logical name for the POINT I/O API INVENTORY file',
     &        FSREAD3, 'PNTS', PROGNAME )

C.........  Try to get header information from inventory file
        IF( .NOT. DESC3( ENAME ) ) THEN
            MESG = 'Could not read description for "' //
     &             ENAME( 1:TRIMLEN( ENAME ) ) // '"'
 
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
 
C.........  Store header information that will be needed later
        ELSE
            NPSRC   = NROWS3D
            IFDESC2 = GETCFDSC( FDESC3D, '/FROM/' )
            IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/' )

        ENDIF

C.........  Get grid name from the environment and read grid parameters
        IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP3D, COORUNIT,
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'ERROR: Could not get Models-3 grid description'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C.........  Store grid parameters for later processing
        ELSE

            NGRID = NCOLS3D * NROWS3D

        ENDIF

C.........  Write message stating grid name and description
        L1 = TRIMLEN( GRDNM )
        MESG = 'Grid "' // GRDNM( 1:L1 ) // '" set; defined as' // 
     &         CRLF() // BLANK5 // GDESC
        CALL M3MSG2( MESG )

C.........  Allocate memory for gridding matrix
        ALLOCATE( GMAT( NGRID + NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GMAT', PROGNAME )

C.........  Read point source characteristics
        CALL M3MSG2( 'Reading in POINT SOURCES file...' )

C.........  Allocate memory for and read in required inventory characteristics
        CALL RPNTSCHR( ENAME, 0, NPSRC, NINVARR, IVARNAMS )

C.........  Call subroutine to convert grid coordinates from lat-lon to
C           coordinate system of the destination grid
        CALL CONVRTXY( NPSRC, XLOCA, YLOCA, GDTYP3D, P_ALP3D, P_BET3D,
     &                 P_GAM3D, XCENT3D, YCENT3D )

C.........  Get file name; open output gridding matrix file
C.........      with grid characteristics from DSCM3GRD() above        
        FTYPE3D = SMATRX3
        SDATE3D = 0
        STIME3D = 0
        TSTEP3D = 0
        NVARS3D = 0
        NCOLS3D = NPSRC
        NROWS3D = NGRID
        NLAYS3D = 1
        NTHIK3D = NPSRC
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        GDNAM3D = GRDNM         

        DO K = 1, MXDESC3
            FDESC3D( K ) = ' '
        ENDDO

        FDESC3D( 1  ) = 'Point source gridding matrix'
        FDESC3D( 2  ) = '/FROM/ ' // PROGNAME
        FDESC3D( 3  ) = '/VERSION/ ' // VERCHAR( SCCSW )
        FDESC3D( 4  ) = '/GDESC/ ' // GDESC
        FDESC3D( 11 ) = '/PNTS FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/PNTS VERSION/ ' // IFDESC3

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRIDDING MATRIX output file',
     &          FSUNKN3, 'PGMAT', PROGNAME )

C.........  Generate point source gridding matrix
   
        CALL GENPGMAT( GNAME, NPSRC, NGRID, XLOCA, YLOCA, 
     &                 GMAT( 1 ), GMAT( NGRID+1 ), NK, CMAX, CMIN )

C.........  Report statistics:

        CAVG = FLOAT( NK ) / FLOAT( NGRID )
        WRITE( LDEV,92010 ) 
     &      'Total number of coefficients   ', NK   ,
     &      'Max  number of sources per cell', CMAX,
     &      'Min  number of sources per cell', CMIN
        WRITE( LDEV,92020 ) 
     &      'Mean number of sources per cell', CAVG,
     &      ' '

C...............   End of program
      
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, :, I12 )

92020   FORMAT( 5X, A, :, F17.4 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )


        END


