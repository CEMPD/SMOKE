
        SUBROUTINE OPENGMAT( NMATX, RLZN, INVPROG, INVVERS, GNAME, 
     &                       UNAME, FDEV ) 

C***********************************************************************
C  subroutine body starts at line 95
C
C  DESCRIPTION:
C      This subroutine sets up the header and variables for the gridding matrix
C      and ungridding matrix for mobile sources.  It prompts for the file
C      names and opens the files.
C
C  PRECONDITIONS REQUIRED:
C      Correct number of sources and matrix dimensions are set.  Grid parameters
C      are defined and passed the FDESC3D include file.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 5/99 by M. Houyoux
C      Modified 5/02 by Gabe Cano - deterministic/stochastic mode
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

C.........  This module contains the uncertainty information
        USE MODUNCERT

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptionsNRAWIN
        CHARACTER*2            CRLF
        LOGICAL                DSCM3GRD
        INTEGER                PROMPTFFILE 
        CHARACTER(LEN=NAMLEN3) PROMPTMFILE
        CHARACTER*16           VERCHAR
        LOGICAL                SETENVVAR

        EXTERNAL CRLF, DSCM3GRD, PROMPTFFILE, PROMPTMFILE, VERCHAR,
     &           SETENVVAR

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN)    :: NMATX   ! no. of source-cell intersections
        INTEGER     , INTENT (IN)    :: RLZN    ! no. of uncertainty realizations
        CHARACTER(*), INTENT (IN)    :: INVPROG ! inventory program
        CHARACTER(*), INTENT (IN)    :: INVVERS ! inventory program version
        CHARACTER(*), INTENT(IN OUT) :: GNAME   ! gridding matrix logical name
        CHARACTER(*), INTENT(IN OUT) :: UNAME   ! ungridding matrix logical name
        INTEGER     , INTENT(IN OUT) :: FDEV    ! report file

C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name$' ! CVS release tag

C...........   Other local variables

        INTEGER       I, J, L1, L2, V     ! counter and indices
        INTEGER       NCOLS     ! number of columns
        INTEGER       NIOVARS   ! number of I/O API file non-emis variables
        INTEGER       NPOLMAX   ! max no of pollutants, based on I/O API
        INTEGER       NNPVAR    ! no. non-pollutant inventory variables
        INTEGER       NROWS     ! number of rows

        CHARACTER*16  COORD3D   ! coordinate system name
        CHARACTER*16  COORUN3D  ! coordinate system projection units
        CHARACTER*16  SNAME     ! logical name for the supplemental file
        CHARACTER*80  GDESC     ! grid description
        CHARACTER*300 LMNAME    ! full path and file name for ?GMATU
        CHARACTER*300 LSNAME    ! full path and file name for ?GSUPU
        CHARACTER*300 LUNAME    ! full path and file name for ?UMATU
        CHARACTER*300 MESG      ! message buffer 

        CHARACTER(LEN=NAMLEN3) NAMBUF   ! file name buffer

        CHARACTER*16 :: PROGNAME = 'OPENGMAT' ! program name

C***********************************************************************
C   begin body of subroutine OPENGMAT

        IF( RLZN .GT. 0) THEN

            IF( .NOT.( CLOSE3( GNAME ) ) ) THEN
                WRITE( MESG,94010 ) 
     &                 'WARNING: Unable to close file' // GNAME
                CALL M3MSG2( MESG )
            END IF

            IF( FDEV .GT. 0 ) THEN
                CLOSE( FDEV )
                WRITE( MESG,94010 ) 'NOTE: Supplemental file device ',
     &                 FDEV, 'has been closed'
                CALL M3MSG2( MESG )
            END IF

            IF( .NOT.( CLOSE3( UNAME ) ) .AND.
     &                                  CATEGORY .EQ. 'MOBILE' ) THEN
                WRITE( MESG,94010 ) 
     &                 'WARNING: Unable to close file' // UNAME
                CALL M3MSG2( MESG )
            END IF

            IF( RLZN .EQ. 1 ) THEN

                GNAME = CRL // 'GMATU'
                SNAME = CRL // 'GSUPU'
                UNAME = CRL // 'UMATU'

            END IF
        END IF

C.........  Initialize header by getting the Models-3 grid information file
        IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,
     &                      P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &                      XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                      NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Set up I/O API output file header for gridding matrix
        NCOLS   = NCOLS3D
        NROWS   = NROWS3D
        FTYPE3D = SMATRX3
        NCOLS3D = NMATX
        NROWS3D = NCOLS * NROWS
        NTHIK3D = NSRC

        IF( CATEGORY .EQ. 'POINT' ) THEN
            NVARS3D = 0
        ELSE
            NVARS3D = 1
            VNAME3D( 1 ) = CRL // 'GRDMAT'
            UNITS3D( 1 ) = 'n/a'
            VDESC3D( 1 ) = CATDESC // ' source gridding coefficients'
            VTYPE3D( 1 ) = M3REAL
        END IF

        FDESC3D( 1  ) = CATDESC // 'source gridding matrix'
        FDESC3D( 2  ) = '/FROM/ ' // PROGNAME
        FDESC3D( 3  ) = '/VERSION/ ' // VERCHAR( CVSW )
        FDESC3D( 4  ) = '/GDESC/ ' // GDESC
        WRITE( FDESC3D(5), 94010 ) '/NCOLS3D/ ', NCOLS
        WRITE( FDESC3D(6), 94010 ) '/NROWS3D/ ', NROWS

        FDESC3D( 11 ) = '/INVEN FROM/ ' // INVPROG
        FDESC3D( 12 ) = '/INVEN VERSION/ ' // INVVERS

C.........  Get name of gridding matrix
        IF ( RLZN .GT. 0 ) THEN

            CALL GRFNAME
            IF( .NOT.( SETENVVAR( GNAME, LMNAME ) ) )
     &           CALL M3EXIT( PROGNAME, 0, 0, 
     &                        'Unable to assign setenv for' // 
     &                         LMNAME, 2 )

            NAMBUF = PROMPTMFILE( 
     &               'Enter logical name for GRIDDING MATRIX ' //
     &               'unceratainty output file',
     &               FSUNKN3, GNAME, PROGNAME )
            GNAME = NAMBUF

        ELSE 

            NAMBUF = PROMPTMFILE( 
     &               'Enter logical name for GRIDDING MATRIX' // 
     &               'output file', FSUNKN3, CRL // 'GMAT', PROGNAME )
            GNAME = NAMBUF

        END IF

C.........  Set up I/O API output file header for ungridding matrix
C.........  Leave everything the same as for the gridding matrix, but a couple
C           of header items

        UNAME  = 'NONE'
        IF( CATEGORY .EQ. 'MOBILE' ) THEN

            NROWS3D = NSRC
            NTHIK3D = NCOLS * NROWS
            FDESC3D( 1 ) = CATDESC // ' source ungridding matrix'

            VNAME3D( 1 ) = CRL // 'UGRDMAT'
            UNITS3D( 1 ) = 'n/a'
            VDESC3D( 1 ) = CATDESC // ' source ungridding coefficients'
            VTYPE3D( 1 ) = M3REAL

            IF( RLZN .GT.0 ) THEN

                IF( .NOT.( SETENVVAR( UNAME, LUNAME ) ) )
     &               CALL M3EXIT( PROGNAME, 0, 0, 
     &                            'Unable to assign setenv for '
     &                            // LUNAME, 2 )

                NAMBUF = PROMPTMFILE( 
     &               'Enter logical name for UNGRIDDING MATRIX ' //
     &               'output file', FSUNKN3, UNAME, PROGNAME )
                UNAME = NAMBUF

            ELSE 

                NAMBUF = PROMPTMFILE( 
     &               'Enter logical name for UNGRIDDING MATRIX ' //
     &               'output file', FSUNKN3, CRL // 'UMAT', PROGNAME )
                UNAME = NAMBUF

            END IF

        END IF

C.........  Open report file
        IF ( ( CATEGORY .EQ. 'AREA' .OR. CATEGORY .EQ. 'MOBILE' )
     &       .AND. RLZN .GT. 0 ) THEN

            IF( .NOT.( SETENVVAR( SNAME, LSNAME ) ) )
     &           CALL M3EXIT( PROGNAME, 0, 0, 
     &                        'Unable to assign setenv for ' //
     &                        LSNAME, 2 )
            MESG = 'Enter logical name for the GRIDDING SUPPLEMENTAL '//
     &             'uncertainty file'
            FDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 
     &                          SNAME, PROGNAME )

        ELSE IF ( ( CATEGORY .EQ. 'AREA' .OR. 
     &           CATEGORY .EQ. 'MOBILE' ) ) THEN

            MESG = 'Enter logical name for the GRIDDING SUPPLEMENTAL '//
     &             'file'
            FDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 
     &                          CRL // 'GSUP', PROGNAME )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        CONTAINS

            SUBROUTINE GRFNAME

            INTEGER        L1, L2, L3 

            CHARACTER*16  FMNAM     ! GMATU file name buffer 
            CHARACTER*16  FSNAM     ! GSUPU file name buffer 
            CHARACTER*16  FUNAM     ! UMATU file name buffer 
            CHARACTER*16  TMPBUF    ! temporary string buffer   
            CHARACTER*300 DIRBUF    ! directory buffer 

C**********************************************************************

C.............  Set up name for uncertainty output files
            CALL NAMEVAL( CRL // 'GMATU_ODIR', DIRBUF )
            L1 = LEN_TRIM( DIRBUF )

            WRITE( TMPBUF,94010 ) '', RLZN
            TMPBUF = ADJUSTL( TMPBUF )
            CALL PADNZERO( RMXLEN, TMPBUF )
            L3 = LEN_TRIM( TMPBUF )

            CALL NAMEVAL( CRL // 'GMATU_FNAM', FMNAM )
            L2 = LEN_TRIM( FMNAM )
            LMNAME = DIRBUF( 1:L1 ) // '/' // 
     &               FMNAM( 1: L2 )  // '_' //
     &               TMPBUF( 1: L3 )  // '.ncf' 

            CALL NAMEVAL( CRL // 'GSUPU_FNAM', FSNAM )
            L2 = LEN_TRIM( FSNAM )
            LSNAME = DIRBUF( 1:L1 ) // '/' // 
     &               FSNAM( 1: L2 )  // '_' //
     &               TMPBUF( 1: L3 )  // '.txt' 

            IF( CATEGORY .EQ. 'MOBILE' ) THEN

                CALL NAMEVAL( CRL // 'UMATU_FNAM', FUNAM )
                L2 = LEN_TRIM( FUNAM )
                LUNAME = DIRBUF( 1:L1 ) // '/' // 
     &                   FUNAM( 1: L2 )  // '_' //
     &                   TMPBUF( 1: L3 )  // '.ncf' 

            END IF

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE GRFNAME
 
        END SUBROUTINE OPENGMAT

