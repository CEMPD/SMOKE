
        PROGRAM TMPBIO

C***********************************************************************
C  program body starts at line  153
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
C       m3io
C       HRBIO
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC adapted from UAM BEIS2 subroutine HRLYEM
C
C       Revised    8/96 by SL and CJC: algorithm change to match OTAG BEIS2
C             algorithm
C
C       Revised    1/97 by SL: integrated tmpbios.F's functionality into this
C             program. added flexibility in calculating PAR by using
C             1. observational cloud, 2. RGRND from MM5, 3, KUO cloud, 
C             4. KF cloud
C
C***********************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File Version: @(#)$Id$
C  
C COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
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

        INTEGER         GETDATE
        INTEGER         GETMENU
        INTEGER         GETNUM
        LOGICAL         GETYN
        CHARACTER*16    PROMPTMFILE
        INTEGER         PROMPTFFILE
        INTEGER         TRIMLEN
        EXTERNAL        GETDATE, GETMENU,GETNUM, GETYN, PROMPTMFILE, 
     &                  PROMPTFFILE, TRIMLEN

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER         LDEV    !  unit number for log device
        INTEGER         SDATE   !  run starting date (YYYYDDD)
        INTEGER         STIME   !  run starting time (HHMMSS)
        INTEGER         NHRS    !  run duration (hours)
        INTEGER         JDATE   !  current simulation date (YYYYDDD)
        INTEGER         JTIME   !  current simulation time (HHMMSS)
        INTEGER         PARTYPE !  method number to calculate PAR
        INTEGER         HR      !  current simulation hour
        INTEGER         I, J, K, L, C, R  !  loop counters and subscripts
        INTEGER         B, M    !  counters for biogenic, model species
        CHARACTER*16    ENAME   !  logical name for emissions output
        CHARACTER*16    NNAME   !  logical name for normalized-emissions input
        CHARACTER*16    GNAME   !  logical name for GRID_CRO_2D
        CHARACTER*16    M3NAME  !  logical name for MET_CRO_3D
        CHARACTER*16    M2NAME  !  logical name for MET_CRO_2D
        CHARACTER*256   MESG    !  message buffer for M3EXIT()

C.......   Input met and grid variables:
                
        REAL LAT    ( NCOLS, NROWS )    !  grid lat (deg) -90 <= LAT <= 90
        REAL LON    ( NCOLS, NROWS )    !  grid lon (deg) -180 <= LON <= 180 
        REAL TASFC  ( NCOLS, NROWS )    !  level-1 air  temperature (K)
        REAL PRSFC  ( NCOLS, NROWS )    !  pressure (Pa)
        REAL TSOLAR ( NCOLS, NROWS)     !  Photosynthetic Active Radiation (PAR)
        REAL WM2LY, LY2UE, SOL2PAR      !  coefficients to convert RGRND to PAR
        PARAMETER   ( WM2LY  = 0.001433,
     &                LY2UE  = 2916.0
     &              )

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

        CHARACTER*16    NORMV( BTYPES * ( BSPCS - 1 ) )   ! names for VOC vbles
        CHARACTER*16    NORMN( LUSES )          !  names for  NO-emission vbles
        CHARACTER*72    PARMENU( 4 )            ! Methods to calc. PAR
        DATA     PARMENU
     &           / 'Use observational clouds',
     &             'Use MM5 generated RGND',
     &             'Use KUO cloud attenuation',
     &             'Use KF  cloud attenuation'   /
            

C***********************************************************************
C   begin body of program TMPBIO

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

C.......   Get the method to calculate PAR

         PARTYPE = GETMENU( 4, 2, 'How will PAR be calculated?',
     &                      PARMENU )

C.......   Get input files

C        WRITE( LDEV,92000 ) ' ', 'Reading GRID and LAND USE data', ' '

        NNAME = PROMPTMFILE( 
     &          'Enter name for NORMALIZED EMISSIONS input file',
     &          FSREAD3, 'BGRD', 'TMPBIO' )

        M2NAME = PROMPTMFILE( 
     &          'Enter name for 2D METEOROLOGY input file',
     &          FSREAD3, 'SFCCRO', 'TMPBIO' )

        IF ( .NOT. DESC3( M2NAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             M2NAME( 1:TRIMLEN( M2NAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
        END IF

        IF ( PARTYPE .NE. 2 ) THEN 

           M3NAME = PROMPTMFILE( 
     &             'Enter name for 3D METEOROLOGY input file',
     &             FSREAD3, 'LAYCRO', 'TMPBIO' )
           GNAME = PROMPTMFILE( 
     &             'Enter name for 2D GRID PARAMETERS input file',
     &             FSREAD3, 'SFCGRD', 'TMPBIO' )

        END IF    ! Use cloud to calculate PAR


        JDATE = SDATE3D
        JTIME = STIME3D


C...........   Get starting date, time, run duration

        CALL NEXTIME( JDATE, JTIME, 10000*MXREC3D )

        SDATE = GETDATE( SDATE3D, 
     &                  'Enter starting date (YYYYDDD|YYYYMMDD)' )
        STIME = GETNUM( 0, 235959, STIME3D, 
     &                  'Enter starting time  (HHMMSS)' )
        NHRS  = GETNUM( 1, MXREC3D-1, MXREC3D-1, 
     &                  'Enter run duration    (hours)' )


C.......   Build description for, and create/open output file
C.......   (all but variables-table in description is borrowed from M2NAME)

        NVARS3D = MSPCS
        NLAYS3D = 1

        DO  11  M = 1, MSPCS
            VNAME3D( M ) = EMSPC( M )
            UNITS3D( M ) = 'moles/hour'
            VDESC3D( M ) = 'biogenic emissions of the indicated species'
            VTYPE3D( M ) = M3REAL
11      CONTINUE

        FDESC3D( 1 ) = 'Gridded biogenic emissions from SMOKE-BEIS2'

        DO  12  M = 4, MXDESC3
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

        IF ( PARTYPE .NE. 2 ) THEN     ! Use clouds to calc. PAR
           IF ( .NOT. READ3( GNAME, 'LAT', 1, 0, 0, LAT ) ) THEN
              MESG = 'Could not read LAT from file "' //
     &                GNAME( 1:TRIMLEN( GNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
           END IF

           IF ( .NOT. READ3( GNAME, 'LON', 1, 0, 0, LON ) ) THEN
              MESG = 'Could not read LON from file "' //
     &                GNAME( 1:TRIMLEN( GNAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
           END IF
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

           IF ( .NOT. READ3( M2NAME, 'TGD', 1, 
     &          JDATE, JTIME, TASFC ) ) THEN
              MESG = 'Could not read sfc. temp. from file "' //
     &                M2NAME( 1:TRIMLEN( M2NAME ) ) // '"'
              CALL M3EXIT( 'TMPBIO', JDATE, JTIME, MESG, 2 )
           END IF

           IF ( PARTYPE .EQ. 2 ) THEN      ! Use RGRND
 
              IF ( .NOT. READ3( M2NAME, 'RGRND', ALLAYS3, JDATE,
     &             JTIME, TSOLAR(1,1) ) ) THEN
                 MESG = 'Could not read from file MET_CRO_2D'
                 CALL M3EXIT( 'TMPBIO', JDATE, JTIME, MESG, 2 )
              END IF

              SOL2PAR = 0.5 * WM2LY * LY2UE  ! RGRND to PAR conversion factor

              DO 50 I = 1, NCOLS
                 DO 40 J = 1, NROWS
                    TSOLAR( I, J ) = TSOLAR( I, J ) * SOL2PAR ! Calc. PAR
40               END DO
50            END DO
 
              CALL HRBIOS( JDATE, JTIME,
     &                   PINE, DECD, CONF, AGRC, LEAF, OTHR, AVLAI,
     &                   NORNO, TASFC, TSOLAR, EMIS )

           ELSE        ! Use clouds to calculate PAR

              IF ( .NOT. READ3( M3NAME, 'PRES', 1, 
     &                        JDATE, JTIME, PRSFC ) ) THEN
                  MESG = 'Could not read PRES from file "' //
     &                   M3NAME( 1:TRIMLEN( M3NAME ) ) // '"'
                  CALL M3EXIT( 'TMPBIO', JDATE, JTIME, MESG, 2 )
              END IF

              IF ( (PARTYPE.EQ.3) .OR. (PARTYPE.EQ.4) ) THEN
                 DO 70 R = 1, NROWS
                    DO 60 C = 1, NCOLS
                       PRSFC( C, R ) = PRSFC( C, R ) / 100.0  ! Pa to mb
60                  END DO
70               END DO
              END IF

              CALL HRBIO( JDATE, JTIME, LAT, LON, PRSFC, TASFC,
     &                    PINE, DECD, CONF, LEAF, OTHR, AVLAI, NORNO,
     &                    PARTYPE, EMIS )
           END IF 

           IF ( .NOT. WRITE3( ENAME, 'ALL', 
     &                        JDATE, JTIME, EMIS ) ) THEN
               CALL M3EXIT( 'TMPBIO', JDATE, JTIME, 
     &                      'Error writing BIO OUTPUT file' , 2 )
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

