
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM GRDBIO

C***********************************************************************
C  program body starts at line  131
C
C  DESCRIPTION:
C       Computes normalized gridded biogenic emissions in terms of gridded
C       or county level biomass, land use, and emissions factors, and 
C       surrogate factors if biomass is at county level.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC adapted from UAM BEIS2 subroutine NORMALEM
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

        LOGICAL         DSCGRID
        INTEGER         FIND1
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         INDEX1
        INTEGER         LBLANK
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN

        EXTERNAL        DSCGRID, FIND1, GETNUM, GETYN, INDEX1, LBLANK, 
     &                  PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL,
     &                  TRIMLEN


C...........   LOCAL VARIABLES and their descriptions:

        INTEGER         B, C, R, I, J, K, L, M, N  !  loop counters and subscripts
        INTEGER         ISTAT   !  I/O status result
        INTEGER         LINE    !  input line number
        LOGICAL         EFLAG   !  error flag
        INTEGER         CDEV    !  unit number for FIP CODES file
        INTEGER         FDEV    !  unit number for emissions factor file:
        INTEGER         SDEV    !  unit number for surrogate factors
        INTEGER         UDEV    !  unit number for county land use file
        INTEGER         LDEV    !  unit number for log file:
        CHARACTER*16    ENAME   !  logical name for emissions output
        CHARACTER*16    GNAME   !  scratch copy of grid name
        CHARACTER*16    CNAME   !  coordinate system name from DSCGRID()

        CHARACTER*80    INBUF   !  input buffer
        CHARACTER*256   MESG    !  message buffer for M3EXIT()

        REAL            AREA    !  land-use area
        REAL            ATYPE   !  area for this land use type
        INTEGER         NTYPE   !  number of land use types
        CHARACTER*4     TYPE    !  land use type

C.......   Source-level variables

        REAL    EMIS, AFOR, ENOFOR
        REAL    AAG, ENOAG 
        REAL    EOTH ( BSPCS )
        REAL    EPINE( BSPCS )
        REAL    EDECD( BSPCS )
        REAL    ECONF( BSPCS )
        REAL    EAG  ( BSPCS )
        REAL    ELAI ( BSPCS ), ENOGRS, ENOWTF 
        REAL    AOTH, AGRS, AWTF, ALAI, SUMLAI, AVGLAI

C.......   Gridded normalized emissions description output to file BGRD

        REAL PINE( NCOLS, NROWS, BSPCS-1 )         !  pine forest
        REAL DECD( NCOLS, NROWS, BSPCS-1 )         !  deciduous forest
        REAL CONF( NCOLS, NROWS, BSPCS-1 )         !  other coniferous forest
        REAL AGRC( NCOLS, NROWS, BSPCS-1 )         !  grasslands
        REAL LEAF( NCOLS, NROWS, BSPCS-1 )         !  leaf area
        REAL OTHR( NCOLS, NROWS, BSPCS-1 )         !  other biogenic area
        REAL NORM( NCOLS, NROWS, BSPCS-1, BTYPES ) !  aggregate VOC emis

        EQUIVALENCE( PINE, NORM( 1,1,1,1 ) )
        EQUIVALENCE( DECD, NORM( 1,1,1,2 ) )
        EQUIVALENCE( CONF, NORM( 1,1,1,3 ) )
        EQUIVALENCE( AGRC, NORM( 1,1,1,4 ) )
        EQUIVALENCE( LEAF, NORM( 1,1,1,5 ) )
        EQUIVALENCE( OTHR, NORM( 1,1,1,6 ) )

        REAL AVLAI ( NCOLS, NROWS )

        REAL GRASNO( NCOLS, NROWS )
        REAL FORENO( NCOLS, NROWS )
        REAL WETLNO( NCOLS, NROWS )
        REAL AGRINO( NCOLS, NROWS )

C...........   emission factor, vegetation types tables:

        INTEGER      NVEG
        INTEGER      LAI  ( MXID )
        REAL         EMFAC( BSPCS, MXID )
        CHARACTER*4  VEGID( MXID )


C***********************************************************************
C   begin body of program GRDBIO

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program GRDBIO to take the county level biomass, the gridded',
     &  'land use, the emissions factors, and the surrogate factors,',
     &  'and produce gridded normalized biogenic emissions.',
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
            CALL M3EXIT( 'GRDBIO', 0, 0, 'Ending program GRDBIO', 0 )
        END IF

C.......   Get file name; open input raw area source file

        CDEV = PROMPTFFILE( 
     &           'Enter logical name for FIP CODES file',
     &           .TRUE., .TRUE., 'BFIP', 'GRDBIO' )

        FDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION FACTORS file',
     &           .TRUE., .TRUE., 'BFAC', 'GRDBIO' )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for SURROGATE FACTORS file',
     &           .TRUE., .TRUE., 'BGPRO', 'GRDBIO' )

        UDEV = PROMPTFFILE( 
     &           'Enter logical name for GRIDDED LANDUSE file',
     &           .TRUE., .TRUE., 'BCUSE', 'GRDBIO' )


C.......   Set up file description and open emissions file:

        IF ( .NOT. DSCGRID( GRDNM, CNAME, GDTYP3D, 
     &             P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &             XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &             NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            GNAME = GRDNM
            MESG  = 'Could not get description of grid "' //
     &              GNAME( 1:TRIMLEN( GNAME ) ) // 
     &              '" from GRIDDESC file.'
            CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )

        END IF		!  if dscgrid() failed

        GDNAM3D = GRDNM
        FTYPE3D = GRDDED3
        SDATE3D = 0             !  n/a
        STIME3D = 0             !  n/a
        TSTEP3D = 0             !  time independent
        NVARS3D = BTYPES * ( BSPCS - 1 ) + LUSES + 1
        NLAYS3D = 1
        NTHIK3D = 1
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3
        FDESC3D( 1 ) = 'Biogenic Source normalized emissions values.'
        FDESC3D( 2 ) = 'Inventory  type COUNTY_BASED'

        DO  11  I = 3, MXDESC3
            FDESC3D( I ) = ' '
11      CONTINUE

        I = 0
        DO  13  B = 1, BTYPES
        DO  12  M = 1, BSPCS - 1

            I = I + 1
            VNAME3D( I ) = BIOLTYPE( B ) // BIOSPC( M )
            VDESC3D( I ) = 'Normalized emissions--forest land use'
            UNITS3D( I ) = CMISS3
            VTYPE3D( I ) = M3REAL

12      CONTINUE
13      CONTINUE

        I = I + 1
        VNAME3D( I ) = 'AVLAI'
        VDESC3D( I ) = 'Average leaf area index'
        UNITS3D( I ) = CMISS3
        VTYPE3D( I ) = M3REAL

        DO  14  L = 1, LUSES

            I = I + 1
            VNAME3D( I ) = BIOLUSE( L )( 1:TRIMLEN( BIOLUSE( L )))//'NO'
            VDESC3D( I ) = 'Normalized emissions--nonforest land use'
            UNITS3D( I ) = CMISS3
            VTYPE3D( I ) = M3REAL

14      CONTINUE

        ENAME = PROMPTMFILE(  
     &          'Enter logical name for NORMALIZED BIO output file',
     &          FSUNKN3, 'BNORM', 'GRDBIO' )


C.......   Loop:  read emission factor file

        WRITE( LDEV,92000 ) ' ', 'Reading EMISSIONS FACTOR file', ' '
        I     = 0
        EFLAG = .FALSE.

22      CONTINUE		!  head of loop reading EMISSIONS FACTOR file

            I = I + 1
            READ( FDEV, 93010, END=23, IOSTAT=ISTAT )
     &          VEGID( I ),
     &        ( EMFAC( I, J ) , J = 1, NSEF ),
     &          LAI( I )
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'Error', ISTAT, 
     &              'reading EMISSION FACTOR file at line', I
                CALL M3MESG( MESG )
            END IF

            IF ( I .LT. MXID ) GO TO  22 !  end loop reading emissions factor file

        CALL M3EXIT( 'GRDBIO', 0, 0, 
     &               'EMISSIONS FACTOR table overflow', 2 )

23      CONTINUE

        NVEG = I
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRDBIO', 0, 0,
     &                   'Error reading EMISSIONS FACTOR file', 2 )
        END IF

C...........   Fold ug~~>g, hectare~~>m^2 factors into emfac:

        DO  25  J = 1, NSEF
        DO  24  I = 1, NVEG
            EMFAC( I, J ) = MICR2G * HA2MSQ * EMFAC( I, J )
24      CONTINUE
25      CONTINUE


C.......   Initialize normalized-emissions arrays:

        DO 33 I = 1, BSPCS
        DO 32 R = 1, NROWS
        DO 31 C = 1, NCOLS
            PINE( C, R, I ) = 0.0
            DECD( C, R, I ) = 0.0
            CONF( C, R, I ) = 0.0
            AGRC( C, R, I ) = 0.0
            LEAF( C, R, I ) = 0.0
            OTHR( C, R, I ) = 0.0
31      CONTINUE
32      CONTINUE
33      CONTINUE

        DO 35 R = 1, NROWS
        DO 34 C = 1, NCOLS
            AVLAI ( C, R ) = 0.0
            GRASNO( C, R ) = 0.0
            FORENO( C, R ) = 0.0
            WETLNO( C, R ) = 0.0
            AGRINO( C, R ) = 0.0
34      CONTINUE
35      CONTINUE

C.......   Loops:  read and process gridded LANDUSE file:

        WRITE( LDEV,92000 ) ' ', 'Reading LANDUSE file', ' '
        LINE  = 0
        EFLAG = .FALSE.

44      CONTINUE                !  head of loop:  read file header

            READ( UDEV, 93000, IOSTAT=ISTAT ) INBUF
            LINE =  LINE + 1
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading header of LANDUSE file at line', LINE
                CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
            ELSE IF ( INBUF( 1:10 ) .NE. 'ENDHEADER ' ) THEN
                GO TO  44	!  to head of loop reading file header
            END IF


55      CONTINUE                !  head of loop:  read data recordsh

            READ( UDEV, 93000, END = 122, IOSTAT=ISTAT ) C, R, AREA
            LINE =  LINE + 1
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading LANDUSE file at line', LINE
                CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
            END IF

            READ( UDEV, *, IOSTAT=ISTAT ) ATYPE, NTYPE
            LINE =  LINE + 1
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading FOREST AREA from ' //
     &          'LANDUSE at line', LINE
                CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
            END IF

C...........   Initialize accumulator variables for this cell:

            DO  56  M = 1, BSPCS - 1

                EPINE( M ) = 0.0
                EDECD( M ) = 0.0
                ECONF( M ) = 0.0
                EOTH ( M ) = 0.0
                ELAI ( M ) = 0.0 
                EAG  ( M ) = 0.0

56          CONTINUE

            AFOR   = 0.0
            AAG    = 0.0
            AGRS   = 0.0
            AWTF   = 0.0
            AOTH   = 0.0
            ALAI   = 0.0
            ENOFOR = 0.0
            ENOAG  = 0.0
            ENOGRS = 0.0
            ENOWTF = 0.0
            SUMLAI = 0.0

C...........   Land use type:  (rural) forest.  Process subtypes:

            DO  69  N = 1, NTYPE

                READ( UDEV, *, IOSTAT=ISTAT ) TYPE, AREA
                LINE =  LINE + 1
                IF ( ISTAT .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'LANDUSE at line', LINE
                    CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
                END IF

                K     = INDEX1( TYPE, NVEG, VEGID )

                IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LANDUSE file in VEGID at line', LINE
                    CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
                END IF

                L = INDEX1( TYPE, SPTREE, SPFORID )

                IF ( L .GT. 0 ) THEN
                    AFOR   = AFOR   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    ALAI   = ALAI   + AREA
                    DO  61  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( M , K )
61                  CONTINUE
                    ENOFOR = ENOFOR + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Wdcp' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AAG    = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  62  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( M , K )
62                  CONTINUE
                    ENOAG  = ENOAG + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AGRS   = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  63  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( M , K )
63                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.2*AREA
                    AGRS   = AAG    + 0.2*AREA
                    AOTH   = AAG    + 0.6*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  64  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( M , K )
64                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( NO, K )

                ELSE    !  Add area to the forest area, add NO to the forest NO

                    AFOR = AFOR + AREA

                    IF     ( LAI( K ) .EQ. 3 ) THEN
                        DO  65  M = 1, BSPCS - 1
                            EPINE( M ) = EPINE( M ) + AREA*EMFAC( M,K )
65                      CONTINUE

                    ELSE IF( LAI( K ) .EQ. 5 ) THEN

                        DO  66  M = 1, BSPCS - 1
                            EDECD( M ) = EDECD( M ) + AREA*EMFAC( M,K )
66                      CONTINUE

                    ELSE IF( LAI( K ) .EQ. 7 ) THEN

                        DO  67  M = 1, BSPCS - 1
                            ECONF( M ) = ECONF( M ) + AREA*EMFAC( M,K )
67                      CONTINUE

                    ELSE

                        SUMLAI = SUMLAI + AREA * LAI( K )
                        ALAI   = ALAI   + AREA
                        DO  68  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA*EMFAC( M,K )
68                      CONTINUE
                        ENOFOR = ENOFOR + AREA * EMFAC( NO, K )

                    END IF      !  if lai is 3,5,7, or otherwise

                END IF  !  if some spforid, or 'Wdcp' or 'Scwd' or 'Urba or not

69          CONTINUE


C...........   Land use type:  urban forest.  Process subtypes:

            READ( UDEV, *,IOSTAT=ISTAT )  ATYPE, NTYPE
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading URBAN FOREST from ' //
     &          'LANDUSE at line', LINE
                CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
            END IF

            DO  79  N = 1, NTYPE

                READ( UDEV, *, IOSTAT=ISTAT ) TYPE, AREA
                LINE =  LINE + 1
                IF ( ISTAT .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'LANDUSE at line', LINE
                    CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
                END IF

                K     = INDEX1( TYPE, NVEG, VEGID )

                IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LANDUSE file in VEGID at line', LINE
                    CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
                END IF

                L = INDEX1( TYPE, SPTREE, SPFORID )

                IF ( L .GT. 0 ) THEN

                    AFOR   = AFOR   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    ALAI   = ALAI   + AREA
                    DO  71  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( M , K )
71                  CONTINUE
                    ENOFOR = ENOFOR + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Wdcp' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AAG    = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  72  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( M , K )
72                  CONTINUE
                    ENOAG  = ENOAG + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AGRS   = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  73  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( M , K )
73                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.2*AREA
                    AGRS   = AAG    + 0.2*AREA
                    AOTH   = AAG    + 0.7*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  74  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( M , K )
74                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( NO, K )

                ELSE    !  Add area to the forest area, add NO to the forest NO

                    AFOR = AFOR + AREA

                    IF     ( LAI( K ) .EQ. 3 ) THEN
                        DO  75  M = 1, BSPCS - 1
                            EPINE( M ) = EPINE( M ) + AREA*EMFAC( M,K )
75                      CONTINUE

                    ELSE IF( LAI( K ) .EQ. 5 ) THEN

                        DO  76  M = 1, BSPCS - 1
                            EDECD( M ) = EDECD( M ) + AREA*EMFAC( M,K )
76                      CONTINUE

                    ELSE IF( LAI( K ) .EQ. 7 ) THEN

                        DO  77  M = 1, BSPCS - 1
                            ECONF( M ) = ECONF( M ) + AREA*EMFAC( M,K )
77                      CONTINUE

                    ELSE

                        SUMLAI = SUMLAI + AREA * LAI( K )
                        ALAI   = ALAI   + AREA
                        DO  78  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA*EMFAC( M,K )
78                      CONTINUE
                        ENOFOR = ENOFOR + AREA * EMFAC( NO, K )

                    END IF      !!  if lai is 3,5,7, or otherwise

                END IF  !  if some spforid, or 'Wdcp' or 'Scwd' or 'Urba or not

79          CONTINUE


C...........   Land use type:  agriculture.  Process subtypes:

            READ( UDEV, *,IOSTAT=ISTAT )  ATYPE, NTYPE
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading AGRICULTURE from ' //
     &          'LANDUSE at line', LINE
                CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
            END IF

            DO  82  N = 1, NTYPE

                READ( UDEV, *, IOSTAT=ISTAT ) TYPE, AREA
                LINE =  LINE + 1
                IF ( ISTAT .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'LANDUSE at line', LINE
                    CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
                END IF

                K     = INDEX1( TYPE, NVEG, VEGID )

                IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LANDUSE file in VEGID at line', LINE
                    CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
                END IF

                AAG = AAG + AREA
                DO  81  M = 1, BSPCS - 1
                    EAG( M ) = EAG( M ) + AREA * EMFAC( M , K )
81              CONTINUE
                ENOAG  = ENOAG + AREA * EMFAC( NO, K )

82          CONTINUE


C...........   Land use type:  remaining/other.  Process subtypes:

            READ( UDEV, *,IOSTAT=ISTAT )  ATYPE, NTYPE
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading OTHER from ' //
     &          'LANDUSE at line', LINE
                CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
            END IF

            DO  101  N = 1, NTYPE

                READ( UDEV, *, IOSTAT=ISTAT ) TYPE, AREA
                LINE =  LINE + 1
                IF ( ISTAT .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'LANDUSE file at line', LINE
                    CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
                END IF

                K     = INDEX1( TYPE, NVEG, VEGID )

                IF ( K .EQ. 0 ) THEN

                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LANDUSE file in VEGID at line', LINE
                    CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )

                END IF		!  if k = 0

                L = INDEX1( TYPE, RMTREE, OTHERID )

                IF ( L .GT. 0 ) THEN

                    ALAI   = ALAI + AREA
                    AFOR   = AFOR + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  91  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( M , K )
91                  CONTINUE
                    ENOFOR = ENOFOR + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Pacp' ) THEN

                    AAG = AAG + AREA
                    DO  92  M = 1, BSPCS - 1
                        EAG( M ) = EAG( M ) + AREA * EMFAC( M , K )
92                  CONTINUE
                    ENOAG = ENOAG + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Gras' .OR.

     &                    TYPE .EQ. 'Scru' .OR.
     &                    TYPE .EQ. 'Ugra' .OR.
     &                    TYPE .EQ. 'Othe' ) THEN
                    AGRS = AGRS + AREA
                    DO  93  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( M , K )
93                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Wetf' ) THEN

                    AWTF = AWTF + AREA
                    ALAI = ALAI + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  94  M = 1, BSPCS - 1
                        EMIS = AREA * EMFAC( M , K )
                        EOTH( M ) = EOTH( M ) + EMIS
                        ELAI( M ) = ELAI( M ) + EMIS
94                  CONTINUE
                    ENOWTF = ENOWTF + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Wate' .OR.
     &                    TYPE .EQ. 'Barr' .OR.
     &                    TYPE .EQ. 'Uoth' ) THEN

                        AOTH = AOTH + AREA

                ELSE IF ( TYPE .EQ. 'Wdcp' ) THEN

                    ALAI   = ALAI   + AREA
                    AFOR   = AFOR   + 0.5 * AREA
                    AAG    = AAG    + 0.5 * AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  95  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( M , K )
95                  CONTINUE
                    ENOAG = ENOAG + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Scwd' ) THEN

                    ALAI   = ALAI   + AREA
                    AFOR   = AFOR   + 0.5 * AREA
                    AGRS   = AGRS   + 0.5 * AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  96  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( M , K )
96                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2 * AREA
                    AGRS   = AGRS   + 0.2 * AREA
                    AOTH   = AOTH   + 0.6 * AREA
                    DO  97  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( M , K )
97                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( NO, K )

                ELSE IF ( TYPE .EQ. 'Desh' ) THEN

                    AGRS   = AGRS   + 0.5 * AREA
                    AOTH   = AOTH   + 0.5 * AREA
                    DO  98  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( M , K )
98                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( NO, K )

                ELSE

                    AOTH = AOTH + AREA

                    IF ( LAI( K ) .GT. 0 ) THEN

                        ALAI   = ALAI   + AREA
                        SUMLAI = SUMLAI + AREA * LAI( K )
                        DO  99  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA * EMFAC( M,K )
99                      CONTINUE
                        ENOFOR = ENOFOR + AREA * EMFAC( NO, K )

                    ELSE

                        DO  100  M = 1, BSPCS - 1
                            EOTH( M ) = EOTH( M ) + AREA * EMFAC( M,K )
100                      CONTINUE
                        ENOGRS = ENOGRS + AREA * EMFAC( NO, K )

                    END IF

                END IF

101        CONTINUE


C...........   Accumulate gridded results

            DO  111  M = 1, BSPCS - 1
                PINE( C,R,M ) = PINE( C,R,M ) + EPINE( M )
                DECD( C,R,M ) = DECD( C,R,M ) + EDECD( M )
                CONF( C,R,M ) = CONF( C,R,M ) + ECONF( M )
                AGRC( C,R,M ) = AGRC( C,R,M ) + EAG  ( M )
                LEAF( C,R,M ) = LEAF( C,R,M ) + ELAI ( M )
                OTHR( C,R,M ) = OTHR( C,R,M ) + EOTH ( M )
111         CONTINUE
            GRASNO( C,R ) = GRASNO( C,R ) + ENOGRS
            FORENO( C,R ) = FORENO( C,R ) + ENOFOR
            WETLNO( C,R ) = WETLNO( C,R ) + ENOWTF
            AGRINO( C,R ) = AGRINO( C,R ) + ENOAG
            AVLAI ( C,R ) = AVLAI ( C,R ) + SUMLAI / ALAI

           GO TO  55            !  to head of loop:  read LANDUSE data records


122     CONTINUE                !  exit from loop:  read LANDUSE data records


C...............   Write output file:

        I = 0
        DO  312  B = 1, BTYPES
        DO  311  M = 1, BSPCS - 1

            I = I + 1
            IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), ALLAYS3, 0, 0,
     &                         NORM( 1,1,M,B ) ) ) THEN
                MESG = 'Could not write "' //
     &                  VNAME3D( I )( 1: TRIMLEN( VNAME3D( I ) ) ) //
     &                  '" to ' // ENAME
                CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
            END IF

311     CONTINUE        !  end loop on VOC species M
312     CONTINUE        !  end loop on VOC-emitting types B

        IF ( .NOT. WRITE3( ENAME, 'AVLAI', ALLAYS3, 0, 0,
     &                     AVLAI ) ) THEN
            MESG = 'Could not write "AVLAI" to ' // ENAME
            CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'GRASNO', ALLAYS3, 0, 0,
     &                     GRASNO ) ) THEN
            MESG = 'Could not write "GRASNO"to ' // ENAME
            CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'FORENO', ALLAYS3, 0, 0,
     &                     FORENO ) ) THEN
            MESG = 'Could not write "FORENO"to ' // ENAME
            CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'WETLNO', ALLAYS3, 0, 0,
     &                     WETLNO ) ) THEN
            MESG = 'Could not write "WETLNO"to ' // ENAME
            CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'AGRINO', ALLAYS3, 0, 0,
     &                     AGRINO ) ) THEN
            MESG = 'Could not write "AGRINO"to ' // ENAME
            CALL M3EXIT( 'GRDBIO', 0, 0, MESG, 2 )
        END IF


C.......   End of program:

      CALL M3EXIT( 'GRDBIO', 0, 0, 
     &             'Successful completion of program GRDBIO', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx
C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )


C...........   Formatted file I/O formats............ 93xxx
                                   
93000   FORMAT( A )

93010   FORMAT( 1X, A4, 3 F9.0, F6.0, I2 )

93020   FORMAT( 4 F10.0, I10, 2F10.0 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )


C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END

