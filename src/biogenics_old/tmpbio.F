
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM TMPBIO

C***********************************************************************
C  program body starts at line  111
C
C  DESCRIPTION:
C       Computes time stepped gridded biogenic emissions in terms of 
C       normalized gridded emissions from RAWBIO and EDSS meteorology.
C
C  PRECONDITIONS REQUIRED:
C       Standard EDSS meteorology GRID_CRO_2D, MET_CRO_2D, MET_CRO_3D
C       Normalized gridded emissions BGRD from RAWBIO
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC adapted from UAM BEIS2 subroutine HRLYEM
C
C***********************************************************************

        IMPLICIT NONE

C...........   INCLUDES:

      INCLUDE 'PARMS3.EXT'      ! I/O API constants
      INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
      INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
      INCLUDE 'GRDIMS3.EXT'     ! grid-related constants
      INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants


C...........   PARAMETERS and their descriptions:

        REAL        MICR2G      !  conversion factor:  ug~~>g
        REAL        HA2MSQ      !  hectares to square meters
        PARAMETER ( MICR2G    = 1.0E-6 ,
     &              HA2MSQ    = 1.0E4 )

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         INDEX1
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN

        EXTERNAL        GETNUM, GETYN, INDEX1, PROMPTMFILE, TRIMLEN

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER         LDEV    !  unit number for log device
        INTEGER         SDATE   !  run starting date (YYYYDDD)
        INTEGER         STIME   !  run starting time (HHMMSS)
        INTEGER         NHRS    !  run duration (hours)
        INTEGER         JDATE   !  current simulation date (YYYYDDD)
        INTEGER         JTIME   !  current simulation time (HHMMSS)
        INTEGER         HR      !  current simulation hour
        INTEGER         I, J, K, L !  loop counters and subscripts
        INTEGER         B, M    !  counters for biogenic, model species
        CHARACTER*16    ENAME   !  logical name for emissions output
        CHARACTER*16    NNAME   !  logical name for normalized-emissions input
        CHARACTER*16    GNAME   !  logical name for GRID_CRO_2D
        CHARACTER*16    MNAME   !  logical name for MET_CRO_3D

        CHARACTER*256   MESG    !  message buffer for M3EXIT()

C.......   Input met and grid variables:
                
        REAL LAT  ( NCOLS, NROWS )      !  grid lat (deg) -90 <= LAT <= 90
        REAL LON  ( NCOLS, NROWS )      !  grid lon (deg) -180 <= LON <= 180 
        REAL TASFC( NCOLS, NROWS )      !  level-1 air  temperature (K)
        REAL PRSFC( NCOLS, NROWS )      !  pressure (Pa)

C.......   Gridded normalized emissions description input from file BGRD

        REAL NVOC( NCOLS, NROWS, BSPCS-1, BTYPES ) !  normalized VOC emissions
        REAL PINE( NCOLS, NROWS, BSPCS-1 )         !  " for pine
        REAL DECD( NCOLS, NROWS, BSPCS-1 )         !  " for deciduous forest
        REAL CONF( NCOLS, NROWS, BSPCS-1 )         !  " for coniferous forest
        REAL AGRC( NCOLS, NROWS, BSPCS-1 )         !  " for agricultural land
        REAL LEAF( NCOLS, NROWS, BSPCS-1 )         !  " for leaf area
        REAL OTHR( NCOLS, NROWS, BSPCS-1 )         !  " for other land

        EQUIVALENCE( PINE, NVOC( 1,1,1,1 ) )
        EQUIVALENCE( DECD, NVOC( 1,1,1,2 ) )
        EQUIVALENCE( CONF, NVOC( 1,1,1,3 ) )
        EQUIVALENCE( AGRC, NVOC( 1,1,1,4 ) )
        EQUIVALENCE( LEAF, NVOC( 1,1,1,5 ) )
        EQUIVALENCE( OTHR, NVOC( 1,1,1,6 ) )

        REAL AVLAI( NCOLS, NROWS )              !  average LAI
        REAL NORNO( NCOLS, NROWS, LUSES )       !  normalized NO emissions

C.......   BEIS2 internal, output species

        REAL EMIS( NCOLS, NROWS, MSPCS )

C.......   Name tables for file NNAME

        CHARACTER*16    NORMV( BTYPES * ( BSPCS - 1 ) )  !  names for VOC-emission vbles
        CHARACTER*16    NORMN( LUSES )                  !  names for  NO-emission vbles


C***********************************************************************
C   begin body of program TMPBIO
        open(unit = 98, file = 'sgrid.txt')

        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program TMPBIO to take EDSS meteorology and normalized',
     &  'gridded emissions from RAWBIO, and produce time stepped',
     &  'gridded speciated biogenic emissions.',
     &  ' ',
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").',
     &  ' ',
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'TMPBIO', 0, 0, 'Ending program TMPBIO', 0 )
        END IF

C.......   Get cross-point grid parameters file and its description

        WRITE( LDEV,92000 ) ' ', 'Reading GRID and LAND USE data', ' '

        NNAME = PROMPTMFILE( 
     &          'Enter name for NORMALIZED EMISSIONS input file',
     &          FSREAD3, 'BNORM', 'TMPBIO' )

        GNAME = PROMPTMFILE( 
     &          'Enter name for GRID_CRO_2D input file',
     &          FSREAD3, 'GRID_CRO_2D', 'TMPBIO' )

        MNAME = PROMPTMFILE( 
     &          'Enter name for MET_CRO_3D input file',
     &          FSREAD3, 'MET_CRO_3D', 'TMPBIO' )

        IF ( .NOT. DESC3( MNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             MNAME( 1:TRIMLEN( MNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
        END IF

        JDATE = SDATE3D
        JTIME = STIME3D


C...........   Get starting date, time, run duration

        CALL NEXTIME( JDATE, JTIME, 10000*MXREC3D )
        SDATE = GETNUM( SDATE3D, JDATE, SDATE3D, 
     &                  'Enter starting date (YYYYDDD)' )
        STIME = GETNUM( 0, 235959, STIME3D, 
     &                  'Enter starting time  (HHMMSS)' )
        NHRS  = GETNUM( 1, MXREC3D, MXREC3D, 
     &                  'Enter run duration    (hours)' )


C.......   Build description for, and create/open output file
C.......   (all but variables-table in description is borrowed from M2NAME)

        NVARS3D = MSPCS

        DO  11  M = 1, MSPCS
            VNAME3D( M ) = EMSPC( M )
            UNITS3D( M ) = 'moles/hour'
            VDESC3D( M ) = 'biogenic emissions of the indicated species'
            VTYPE3D( M ) = M3REAL
11      CONTINUE

        FDESC3D( 1 ) = 'Gridded biogenic emissions from SMOKE-BEIS2'
        DO  12  M = 2, MXDESC3
            FDESC3D( M ) = ' '
12      CONTINUE
        
        ENAME = PROMPTMFILE( 
     &          'Enter name for BGTS output file',
     &          FSUNKN3, 'BGTS', 'TMPBIO' )


C.......   Build name table for variables in normalized emissions file"

        I = 0
        DO  23  K = 1, BTYPES
        DO  22  B = 1, BSPCS - 1

            I = I + 1
            NORMV( I ) = BIOLTYPE( K ) // BIOSPC( B )

22      CONTINUE
23      CONTINUE

        DO  24  L = 1, LUSES

            NORMN( L ) = BIOLUSE( L )( 1:TRIMLEN( BIOLUSE( L )))//'NO'

24      CONTINUE

C...........   Read time-independent data:

        IF ( .NOT. READ3( GNAME, 'LAT', 1, 0, 0, LAT ) ) THEN
             MESG = 'Could not read LAT from file "' //
     &              GNAME( 1:TRIMLEN( GNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. READ3( GNAME, 'LON', 1, 0, 0, LON ) ) THEN
             MESG = 'Could not read LON from file "' //
     &              GNAME( 1:TRIMLEN( GNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
        END IF

C.......   Loops reading the various categories of normalized emissions:

        I = 0
        DO  34  K = 1, BTYPES
        DO  33  J = 1, BSPCS - 1

            I = I + 1
            IF ( .NOT. READ3( NNAME, NORMV( I ), 1, 0, 0, 
     &                      NVOC( 1,1,J,K ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  NORMV( I )( 1 : TRIMLEN( NORMV( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

33      CONTINUE        !  end loop on VOC's
34      CONTINUE        !  end loop on VOC-emitting land types

        IF ( .NOT. READ3( NNAME, 'AVLAI', 1, 0, 0, AVLAI ) ) THEN
             MESG = 'Could not read AVLAI from file "' //
     &              NNAME( 1:TRIMLEN( NNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
        END IF

        I = 0
        DO  35  J = 1, BSPCS

            I = I + 1
            IF ( .NOT. READ3( NNAME, NORMN( I ), 1, 0, 0, 
     &                      NORNO( 1,1,J ) ) ) THEN
                 MESG = 'Could not read "' // 
     &                  NORMN( I )( 1 : TRIMLEN( NORMN( I ) ) ) //
     &                  '" from file "' //
     &                  NNAME( 1:TRIMLEN( NNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF

35      CONTINUE        !  end loop reading normalized NO's


C.......   Loop on time steps:

        JDATE = SDATE
        JTIME = STIME

        DO  899  HR = 1, NHRS

            IF ( .NOT. READ3( MNAME, 'TA', 1, 
     &                      JDATE, JTIME, TASFC ) ) THEN
                MESG = 'Could not read TA from file "' //
     &                 MNAME( 1:TRIMLEN( MNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', JDATE, JTIME, MESG, 2 )
            END IF

            IF ( .NOT. READ3( MNAME, 'PRES', 1, 
     &                      JDATE, JTIME, PRSFC ) ) THEN
                MESG = 'Could not read PRES from file "' //
     &                 MNAME( 1:TRIMLEN( MNAME ) ) // '"'
                CALL M3EXIT( 'TMPBIO', JDATE, JTIME, MESG, 2 )
            END IF

            CALL HRBIO( JDATE, JTIME, LAT, LON, PRSFC, TASFC,
     &                  PINE, DECD, CONF, LEAF, OTHR, AVLAI, NORNO,
     &                  EMIS )

            IF ( .NOT. WRITE3( ENAME, 'ALL', 
     &                         JDATE, JTIME, EMIS ) ) THEN
                CALL M3EXIT( 'TMPBIO', JDATE, JTIME, 
     &                       'Error writing BIO OUTPUT file' , 2 )
            END IF                              !  if write3 failed

            CALL NEXTIME( JDATE, JTIME, 10000 )

899     CONTINUE                !  end loop on hours HR


C.......   End of program:

      CALL M3EXIT( 'TMPBIO', 0, 0, 
     &             'Successful completion of program TMPBIO', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx

91000   FORMAT ( //5X , '*** ERROR ABORT in program TMPBIO ***',
     &            /5X , A ,
     &           // )        !  generic error message format


C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )


C...........   Formatted file I/O formats............ 93xxx
                                   
93010   FORMAT( 1X, A4, 3 F9.0, F6.0, I2 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

