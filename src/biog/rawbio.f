
        PROGRAM RAWBIO

C***********************************************************************
C  program body starts at line  116
C
C  DESCRIPTION:
C       Computes normalized gridded biogenic emissions in terms of 
C       county level biomass, land use, emissions factors, and 
C       surrogate factors.  The FIPS codes to use from the county
C       level biomass file are obtained from the surrogate factors
C       file.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C        RDSRGHDR, RDSRG, RDBEFAC, CYBIO, GRDBIO 
C
C  REVISION  HISTORY:
C       Prototype 11/99 by JMV from version 4.2 of RAWBIO SMOKE prototype 
C       02/00 JMV: added option of using gridded landuse and reorganized
C                  by creating CYBIO and GRDBIO routines. 
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

C...........   This module contains the gridding surrogates tables
        USE MODSURG

C...........   This module contains biogenic variables
        USE MODBIOG
 
C.........  This module contains the global variables for the 3-d grid
        USE MODGRID

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
        INCLUDE 'EMCNST3.EXT'     !
        INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants

C...........   PARAMETERS and their descriptions:

        REAL        MICR2G      !  conversion factor:  ug~~>g
        REAL        HA2MSQ      !  hectares to square meters

C...........   LOCAL PARAMETERS

        CHARACTER*50  CVSW          ! CVS release tag

        PARAMETER ( MICR2G    = 1.0E-6 ,
     &              HA2MSQ    = 1.0E4  , 
     &              CVSW      = '$Name$' )

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2     CRLF
        LOGICAL         DSCM3GRD
        LOGICAL         ENVYN
        INTEGER         GETFLINE
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        CHARACTER*16    VERCHAR

        EXTERNAL        CRLF, DSCM3GRD, ENVYN, GETFLINE, PROMPTFFILE, 
     &                  PROMPTMFILE, VERCHAR

C...........   File names and unit numbers
        INTEGER         FDEV    !  unit number for emissions factor file:
        INTEGER         GDEV    !  unit number for gridded land use file
        INTEGER         LDEV    !  unit number for log file
        INTEGER         SDEV    !  unit number for surrogate factors
        INTEGER         UDEV    !  unit number for county land use file

        CHARACTER*16    ENAME   !  logical name for emissions output

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER         B, C, R, I, J, K, L, M, N ! loop counters and subscripts

        INTEGER         DATNROWS            ! no. rows in input data file
        INTEGER         DATNCOLS            ! no. cols in input data file
        INTEGER         IOS                 ! i/o status result

        LOGICAL ::      EFLAG = .FALSE.     ! true: error found
        LOGICAL ::      GLUSE_YN = .TRUE.   ! gridded or county luse to be used
        LOGICAL ::      LGRID = .FALSE.     ! true: use gridded land use

        CHARACTER*16     COORUNIT           ! coordinate system projection units
        CHARACTER*16  :: DATGRDNM = ' '     ! input data grid name
        CHARACTER*16     SRGFMT             ! surrogates format
        CHARACTER*80     GDESC              ! grid description
        CHARACTER*256    MESG               ! message buffer

        CHARACTER*16 :: PROGNAME = 'RAWBIO'   !  program name

C***********************************************************************
C   begin body of program RAWBIO

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, CVSW, PROGNAME )
   
C.......   Get file name; open emission factors file

        FDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION FACTORS file',
     &           .TRUE., .TRUE., 'BFAC', PROGNAME )

C.......   Determine if gridded landuse is to be used

        MESG = 'Use gridded land use or not'
        LGRID = ENVYN ( 'GLUSE_YN', MESG, .FALSE., IOS )

C.......   If county landuse is to be used open surrogates and 
C          county landuse file

        IF ( .NOT. LGRID ) THEN

C.......   Get file name; open county landuse file

           UDEV = PROMPTFFILE( 
     &           'Enter logical name for COUNTY LANDUSE file',
     &           .TRUE., .TRUE., 'BCUSE', PROGNAME )

C.......   Get file name; open surrogates fractions file

           SDEV = PROMPTFFILE( 
     &           'Enter logical name for GRIDDING SURROGATES file',
     &           .TRUE., .TRUE., 'BGPRO', PROGNAME )

           CALL M3MSG2( 'Reading gridding surrogates header...' )

C.............  Read the surrogates header, obtain the format of the file, 
C               and save the name of the input grid
           CALL RDSRGHDR( SDEV, SRGFMT )
           DATGRDNM = GRDNM
           DATNCOLS = NCOLS
           DATNROWS = NROWS

C.........  Use gridded landuse file...
        ELSE

C.......   Get file name; open county landuse file

           GDEV = PROMPTFFILE(
     &           'Enter logical name for GRIDDED LANDUSE file',
     &           .TRUE., .TRUE., 'BGUSE', PROGNAME )

C.............  Read the header of landuse file, obtain the format of the file, 
C               and save the name of the input grid
           CALL RDSRGHDR(  GDEV, SRGFMT )
           DATGRDNM = GRDNM
           DATNCOLS = NCOLS
           DATNROWS = NROWS

C.............   Rewind gridded landuse file
           REWIND ( GDEV )

        END IF

C.........  Get grid name from the environment and read grid parameters
        IF ( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD, GDTYP3D, COORUNIT,
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, 
     &                     YCENT3D, XORIG3D, YORIG3D, XCELL3D,
     &                     YCELL3D, NCOLS3D, NROWS3D, NTHIK3D)) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C.........  Check the output grid grid settings as compared to the input
C           surrogates or land use data
        CALL CHKGRID( GDNAM3D, 'GRIDDESC', 1, EFLAG )

C.........  Abort if error
        IF ( EFLAG ) THEN
            MESG = 'Inconsistency between grid selected and ' //
     &             'grid in input data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Write message stating grid name and description
        L = LEN_TRIM( GRDNM )
        MESG = 'NOTE: Output grid "' // GRDNM( 1:L ) // 
     &         '" set; described as' // CRLF() // BLANK10 // GDESC
        CALL M3MSG2( MESG )

C.........  If output grid is different from surrogates, write message
        IF ( OFFLAG ) THEN
            L = LEN_TRIM( DATGRDNM )
            MESG = 'NOTE: gridding surrogates extracted for output '//
     &             'grid from grid "' // DATGRDNM( 1:L ) // '"'
            CALL M3MSG2( MESG )
        END IF

C.........  If surrogates are needed, allocate memory for and read the 
C           gridding surrogates
        IF( .NOT. LGRID ) THEN

            CALL RDSRG( SDEV, SRGFMT, DATNROWS, DATNCOLS )

        END IF

C............. Set up header variables for output file BGRD...

C............. Some grid description info obtained from call to RDSRG above

        FTYPE3D = GRDDED3
        NROWS3D = NROWS
        NCOLS3D = NCOLS
        GDNAM3D = GRDNM
        SDATE3D = 0       !  n/a
        STIME3D = 0       !  n/a
        TSTEP3D = 0       !  time independent
        NVARS3D = BTYPES * ( BSPCS - 1 ) + LUSES + 1
        NLAYS3D = 1
        NTHIK3D = 1
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3

        FDESC3D = ' '
        FDESC3D( 1 ) = 'Biogenic Source normalized emissions values.'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        IF ( .NOT. LGRID ) THEN
          FDESC3D( 4 ) = '/LANDUSE/ COUNTY '
        ELSE
          FDESC3D( 4 ) = '/LANDUSE/ GRIDDED ' 
        END IF

        I = 0
        DO  M = 1, BSPCS - 1
          DO B = 1, BTYPES
            I = I + 1
            VNAME3D( I ) = BIOLTYPE( B ) // BIOSPC( M )
            VDESC3D( I ) = 'Normalized emissions--forest land use'
            UNITS3D( I ) = 'grams/hour' 
            VTYPE3D( I ) = M3REAL

          END DO
        END DO

        I = I + 1
        VNAME3D( I ) = 'AVLAI'
        VDESC3D( I ) = 'Average leaf area index'
        UNITS3D( I ) = 'index'
        VTYPE3D( I ) = M3REAL

        DO  L = 1, LUSES

            I = I + 1
            VNAME3D( I ) = BIOLUSE( L )( 1:LEN_TRIM( BIOLUSE(L)))//'NO'
            VDESC3D( I ) = 'Normalized emissions--nonforest land use'
            UNITS3D( I ) = 'grams/hour' 
            VTYPE3D( I ) = M3REAL

        END DO

        ENAME = PROMPTMFILE(  
     &          'Enter logical name for NORMALIZED BIO output file',
     &          FSUNKN3, 'BGRD', PROGNAME )

C.........  Get number of records of BFAC file
        NVEG = GETFLINE( FDEV, 'Emissions factor file' )

C.........  Allocate memory for emission factor variables   
        ALLOCATE( VEGID ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VEGID', PROGNAME )
        ALLOCATE( EMFAC ( NVEG, NSEF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMFAC', PROGNAME )
        ALLOCATE( LAI ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LAI', PROGNAME )

C.........  Read emissions factor file
        CALL M3MSG2( 'Reading EMISSIONS FACTOR file...' )

        CALL RDBEFAC( FDEV, NVEG, VEGID, EMFAC, LAI ) 

C.........  Fold ug~~>g, hectare~~>m^2 factors into emfac:

        DO  J = 1, NSEF
          DO  I = 1, NVEG
            EMFAC( I, J ) = MICR2G * HA2MSQ * EMFAC( I, J )
          END DO
        END DO

C.........  Allocate memory and initialize variables for normalized
C           emissions categories

        ALLOCATE( PINE ( NGRID, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PINE', PROGNAME )
        ALLOCATE( DECD ( NGRID, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DECD', PROGNAME )
        ALLOCATE( CONF ( NGRID, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CONF', PROGNAME )
        ALLOCATE( AGRC ( NGRID, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AGRC', PROGNAME )
        ALLOCATE( LEAF ( NGRID, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LEAF', PROGNAME )
        ALLOCATE( OTHR ( NGRID, BSPCS-1  ), STAT=IOS )
        CALL CHECKMEM( IOS, 'OTHR', PROGNAME )

        ALLOCATE( AVLAI ( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVLAI', PROGNAME )
        ALLOCATE( GRASNO ( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRASNO', PROGNAME )
        ALLOCATE( FORENO ( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FORENO', PROGNAME )
        ALLOCATE( WETLNO ( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WETLNO', PROGNAME )
        ALLOCATE( AGRINO ( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AGRINO', PROGNAME )

        PINE = 0.0      ! array
        DECD = 0.0      ! array
        CONF = 0.0      ! array
        AGRC = 0.0      ! array
        LEAF = 0.0      ! array
        OTHR = 0.0      ! array

        AVLAI  = 0.0    ! array
        GRASNO = 0.0    ! array
        FORENO = 0.0    ! array
        WETLNO = 0.0    ! array
        AGRINO = 0.0    ! array

C............  Calculate normalized biogenic emissions using
C............  either county or gridded landuse

        IF ( .NOT. LGRID ) THEN

          CALL CYBIO( UDEV, NGRID )

        ELSE

          CALL GRDBIO( GDEV, DATNROWS, DATNCOLS )

        END IF
 
C...............   Write output file:

        I = 0        
        DO  M = 1, BSPCS - 1

            I = I + 1
            IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                         PINE( 1,M ) ) ) THEN
                MESG = 'Could not write "' //
     &                  VNAME3D( I )( 1: LEN_TRIM( VNAME3D( I ) ) ) //
     &                  '" to ' // ENAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            I = I + 1
            IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                         DECD( 1,M ) ) ) THEN
                MESG = 'Could not write "' //
     &                  VNAME3D( I )( 1: LEN_TRIM( VNAME3D( I ) ) ) //
     &                  '" to ' // ENAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = I + 1
            IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                         CONF( 1,M ) ) ) THEN
                MESG = 'Could not write "' //
     &                  VNAME3D( I )( 1: LEN_TRIM( VNAME3D( I ) ) ) //
     &                  '" to ' // ENAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = I + 1
            IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                         AGRC( 1,M ) ) ) THEN
                MESG = 'Could not write "' //
     &                  VNAME3D( I )( 1: LEN_TRIM( VNAME3D( I ) ) ) //
     &                  '" to ' // ENAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = I + 1
            IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                         LEAF( 1,M ) ) ) THEN
                MESG = 'Could not write "' //
     &                  VNAME3D( I )( 1: LEN_TRIM( VNAME3D( I ) ) ) //
     &                  '" to ' // ENAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            I = I + 1
            IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                         OTHR( 1,M ) ) ) THEN
                MESG = 'Could not write "' //
     &                  VNAME3D( I )( 1: LEN_TRIM( VNAME3D( I ) ) ) //
     &                  '" to ' // ENAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END DO        !  end loop on VOC species M

        IF ( .NOT. WRITE3( ENAME, 'AVLAI', 0, 0,
     &                     AVLAI ) ) THEN
            MESG = 'Could not write "AVLAI" to ' // ENAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'GRASNO', 0, 0,
     &                     GRASNO ) ) THEN
            MESG = 'Could not write "GRASNO"to ' // ENAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'FORENO', 0, 0,
     &                     FORENO ) ) THEN
            MESG = 'Could not write "FORENO"to ' // ENAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'WETLNO', 0, 0,
     &                     WETLNO ) ) THEN
            MESG = 'Could not write "WETLNO"to ' // ENAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'AGRINO', 0, 0,
     &                     AGRINO ) ) THEN
            MESG = 'Could not write "AGRINO"to ' // ENAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF


C.........   End of program:

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )

C...........   Formatted file I/O formats............ 93xxx
                                   
93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx
94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

        END PROGRAM  RAWBIO 

