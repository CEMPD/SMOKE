C ****************************************************************************
C *
C * Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C *                System
C * File Version: @(#)$Id$
C *
C * COPYRIGHT (C) 1996, MCNC--North Carolina Supercomputing Center
C * All Rights Reserved
C *
C * See file COPYRIGHT for conditions of use.
C *
C * Environmental Programs Group
C * MCNC--North Carolina Supercomputing Center
C * P.O. Box 12889
C * Research Triangle Park, NC  27709-2889
C *
C * env_progs@mcnc.org
C *
C * Pathname: $Source$
C * Last updated: $Date$ 
C *
C ****************************************************************************/

        PROGRAM TMPBIOS

C***********************************************************************
C  program body starts at line  134
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
C       HRBIOS
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC adapted from UAM BEIS2 subroutine HRLYEM
C
C       Version 8/96 by SL and CJC to use Dyntel-style ASCII solar radiation
C       and surface temperature files instead of MCIP meteorology.
C       Name change HRBIO ~~> HRBIOS to accompany arglist change.
C   
C       Version 8/96 by SL to use PAR converted from RGRND in MCIP, and 
C       surface temperature (at 1.5 m) calculated by METCDF using Daewon 
C       Byun's Similarity theory
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
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN

        EXTERNAL        GETNUM, GETYN, INDEX1, PROMPTFFILE,PROMPTMFILE, 
     &                  TRIMLEN

C...........    Hardwired mod for read in gridded solar radiation .. sl

c       INTEGER*4       METHR, METREC, MCOL, MROW, INDEV, JMETHR
c       INTEGER*4       INDEV
        REAL*4          TSFC(NCOLS, NROWS), TSOLAR(NCOLS, NROWS)


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
        CHARACTER*16    INDEV   !  logical name for PAR-TEMP input
        CHARACTER*256   MESG    !  message buffer for M3EXIT()

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
C   begin body of program TMPBIOS

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
     &  ' ',
     &  'COPYRIGHT(C)1996, MCNC--North Carolina Supercomputing Center',
     &  'See file COPYRIGHT for conditions of use.',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'TMPBIO', 0, 0, 'Ending program TMPBIO', 0 )
        END IF

C.......   Get input files

c       INDEV = PROMPTFFILE(
c    &          'Enter logical name for SOLAR RADIATION file',
c    &          .TRUE., .TRUE., 'SOLAR', 'TMPBIO' )

        INDEV = PROMPTMFILE( 
     &          'Enter name for PAR and SFC. TEMP input file',
     &          FSREAD3, 'SOLAR', 'TMPBIO' )

        WRITE( LDEV,92000 ) ' ', 'Reading Normalized Emission data', ' '

        NNAME = PROMPTMFILE( 
     &          'Enter name for NORMALIZED EMISSIONS input file',
     &          FSREAD3, 'BNORM', 'TMPBIO' )

        IF ( .NOT. DESC3( INDEV ) ) THEN
            MESG = 'Could not get description of file "' //
     &             NNAME( 1:TRIMLEN( NNAME ) ) // '"'
            CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
        END IF

        JDATE = SDATE3D
        JTIME = STIME3D
        TSTEP3D = 10000

C...........   Get starting date, time, run duration

        CALL NEXTIME( JDATE, JTIME, 10000*MXREC3D )
        SDATE = GETNUM( SDATE3D, JDATE, SDATE3D, 
     &                  'Enter starting date (YYYYDDD)' )
        STIME = GETNUM( 0, 235959, STIME3D, 
     &                  'Enter starting time  (HHMMSS)' )
c       NHRS  = GETNUM( 1, 999999, MXREC3D, 
c    &                  'Enter run duration    (hours)' )
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

        JDATE  = SDATE
        JTIME  = STIME
c       METREC = 1

        DO  899  HR = 1, NHRS

c          DO 200 K = 1, NGRID

c             READ (INDEV, *) METHR, MCOL, MROW, 
c    &               TSFC(MCOL, MROW), TSOLAR(MCOL, MROW)

c             JMETHR = MOD ( (METHR+5), 24 )  ! EST Convert to GMT

c             IF ( JMETHR .NE. (JTIME/10000) ) THEN 
c                WRITE (6, 150) METHR, JTIME/10000
c150              FORMAT (1X, 'METHR =', I4, 'JTIME =', I4)
c                CALL M3EXIT( 'TMPBIO', 0, 0,
c    &                         'METHR and JTIME are not the same', 0 )
c             END IF
 
c             METREC = METREC + 1
 
c200        END DO              !  end loop on this hour's met input.

C.......... Convert from degree C to degree K

c           DO 400 I = 1, NCOLS
c              DO 350 J = 1, NROWS
c                 TSFC(I, J) = TSFC(I, J) + 273.15 
c350            END DO
c400         END DO

            IF ( .NOT. READ3( INDEV, 'PAR', ALLAYS3, JDATE,
     &         JTIME, TSOLAR(1,1) ) ) THEN
               MESG = 'Could not read from file SOLAR'
               CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
 
            IF ( .NOT. READ3( INDEV, 'TA', ALLAYS3, JDATE,
     &         JTIME, TSFC(1,1) ) ) THEN
               MESG = 'Could not read from file SOLAR'
               CALL M3EXIT( 'TMPBIO', 0, 0, MESG, 2 )
            END IF
          

            CALL HRBIOS( JDATE, JTIME,
     &                   PINE, DECD, CONF, AGRC, LEAF, OTHR, AVLAI, 
     &                   NORNO, TSFC, TSOLAR, EMIS )

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

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )

        END

