
        PROGRAM RAWBIO

C***********************************************************************
C  program body starts at line  206
C
C  DESCRIPTION:
C       Computes normalized gridded biogenic emissions in terms of gridded
C       or county level biomass, land use, and emissions factors, and 
C       surrogate factors if biomass is at county level.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       m3io
C
C  REVISION  HISTORY:
C       Prototype 12/95 by CJC adapted from UAM BEIS2 subroutine NORMALEM
C
C       Revised 8/96 by SL and CJC: algorithm change to match UAM BEIS2
C       algorithm
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File Version: @(#)$Id$
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
C************************************************************************

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
        INTEGER         FINDFIP    ! sl For finding unsorted FIP CODE file

        EXTERNAL        DSCGRID, FIND1, GETNUM, GETYN, INDEX1, LBLANK, 
     &                  PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL,
     &                  TRIMLEN,FINDFIP


C...........   LOCAL VARIABLES and their descriptions:

        INTEGER         B, C, R, I, J, K, L, M, N ! loop counters and subscripts
        INTEGER         FIP     !  fip code
        INTEGER         CTY     !  subscript into FIPS table
        INTEGER         ISTAT   !  I/O status result
        INTEGER         LINE    !  I/O status result
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

        CHARACTER*20    CO_NAME !  sl add for reading FIP CODES file right
        REAL            DUMMY1, DUMMY2, DUMMY3

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
        REAL    BUTME, BUTMN, EUTME, EUTMN
        REAL    XCELSZ, YCELSZ, XCOOR, YCOOR, XX, YY
        REAL    FRAC( 100 ), DWIDTH, DHEIGHT
        INTEGER UTMZONE, ZONE


C.......   County-level variables:

        REAL  CTYPINE(  BSPCS , NBFIP )
        REAL  CTYDECD(  BSPCS , NBFIP )
        REAL  CTYCONF(  BSPCS , NBFIP )
        REAL  CTYAGRI(  BSPCS , NBFIP )
        REAL  CTYLEAF(  BSPCS , NBFIP )
        REAL  CTYOTHR(  BSPCS , NBFIP )
        REAL  CTYNOFOR( NBFIP )
        REAL  CTYNOGRS( NBFIP )
        REAL  CTYNOAG ( NBFIP )
        REAL  CTYNOWTF( NBFIP )
        REAL  CTYAVLAI( NBFIP )
        REAL  CTYAFOR ( NBFIP )
        REAL  CTYAAG  ( NBFIP )
        REAL  CTYAGRS ( NBFIP )
        REAL  CTYAWTF ( NBFIP )
        REAL  CTYAOTH ( NBFIP )

C.......   Surrogate indices

        INTEGER         PSRG    !  for pine forest
        INTEGER         DSRG    !  for deciduous forest
        INTEGER         CSRG    !  for coniferous nonpine forest
        INTEGER         ASRG    !  for agriculture
        INTEGER         OSRG    !  for other land uses
        INTEGER         FSRG    !  for forest
        INTEGER         GSRG    !  for grasslands
        INTEGER         WSRG    !  for grasslands
        INTEGER         MXSRG   !  maximum of all of the above

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

C.......   County FIP CODE table

        INTEGER  NFIP                   !  number of actual counties
        INTEGER  FIPS( NBFIP )          !  county FIPS codes


C...........   emission factor, vegetation types tables:

        INTEGER      NVEG
        INTEGER      LAI  ( MXID )
        REAL         EMFAC( MXID,BSPCS )
        CHARACTER*4  VEGID( MXID )

C...............   STATEMENT FUNCTIONS:

        REAL    PP
        REAL    QQ
        LOGICAL FLTERR          !  pp, qq _NOT_ approx equal
        FLTERR( PP, QQ ) = 
     &  ( ( PP - QQ )**2 .GT. 1.0E-10 * ( PP*PP + QQ*QQ + 1.0E-5 ) )
        
C***********************************************************************
C   begin body of program RAWBIO

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program RAWBIO to take the county level biomass, the gridded',
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
            CALL M3EXIT( 'RAWBIO', 0, 0, 'Ending program RAWBIO', 0 )
        END IF

C.......   Get file name; open input raw area source file

        CDEV = PROMPTFFILE( 
     &           'Enter logical name for FIP CODES file',
     &           .TRUE., .TRUE., 'BFIP', 'RAWBIO' )

        FDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION FACTORS file',
     &           .TRUE., .TRUE., 'BFAC', 'RAWBIO' )

        SDEV = PROMPTFFILE( 
     &           'Enter logical name for SURROGATE FACTORS file',
     &           .TRUE., .TRUE., 'BGPRO', 'RAWBIO' )

        UDEV = PROMPTFFILE( 
     &           'Enter logical name for COUNTY LANDUSE file',
     &           .TRUE., .TRUE., 'BCUSE', 'RAWBIO' )


C.......   Set up file description and open emissions file:

        IF ( .NOT. DSCGRID( GRDNM, CNAME, GDTYP3D, 
     &             P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &             XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &             NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            GNAME = GRDNM
            MESG  = 'Could not get description of grid "' //
     &              GNAME( 1:TRIMLEN( GNAME ) ) // 
     &              '" from GRIDDESC file.'
            CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )

        END IF		!  if dscgrid() failed

        GDNAM3D = GRDNM
        FTYPE3D = GRDDED3
        SDATE3D = 1988167       !  n/a
        STIME3D = 050000        !  n/a
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
        DO  23  B = 1, BTYPES
        DO  22  M = 1, BSPCS - 1

            I = I + 1
            VNAME3D( I ) = BIOLTYPE( B ) // BIOSPC( M )
            VDESC3D( I ) = 'Normalized emissions--forest land use'
            UNITS3D( I ) = CMISS3
            VTYPE3D( I ) = M3REAL

22      CONTINUE
23      CONTINUE

        I = I + 1
        VNAME3D( I ) = 'AVLAI'
        VDESC3D( I ) = 'Average leaf area index'
        UNITS3D( I ) = CMISS3
        VTYPE3D( I ) = M3REAL

        DO  24  L = 1, LUSES

            I = I + 1
            VNAME3D( I ) = BIOLUSE( L )( 1:TRIMLEN( BIOLUSE( L )))//'NO'
            VDESC3D( I ) = 'Normalized emissions--nonforest land use'
            UNITS3D( I ) = CMISS3
            VTYPE3D( I ) = M3REAL

24      CONTINUE

        ENAME = PROMPTMFILE(  
     &          'Enter logical name for NORMALIZED BIO output file',
     &          FSUNKN3, 'BGRD', 'RAWBIO' )


C.......   Get surrogate categories for each emissions type:

        PSRG = GETNUM( 1, 15, 1, 
     &         'Enter surrogate index for pine forest' )

        DSRG = GETNUM( 1, 15, 1, 
     &         'Enter surrogate index for deciduous forest' )

        MXSRG = MAX( PSRG, DSRG )

        CSRG = GETNUM( 1, 15, 1, 
     &         'Enter surrogate index for coniferous forest' )

        MXSRG = MAX( MXSRG, CSRG )

        FSRG = GETNUM( 1, 15, 1, 
     &         'Enter surrogate index for all forest' )

        MXSRG = MAX( MXSRG, FSRG )

        ASRG = GETNUM( 1, 15, 1, 
     &         'Enter surrogate index for agriculture' )

        MXSRG = MAX( MXSRG, ASRG )

        GSRG = GETNUM( 1, 15, 1, 
     &         'Enter surrogate index for grasslands' )

        MXSRG = MAX( MXSRG, GSRG )

        WSRG = GETNUM( 1, 15, 1, 
     &         'Enter surrogate index for wetlands' )

        MXSRG = MAX( MXSRG, WSRG )

        OSRG = GETNUM( 1, 15, 1, 
     &         'Enter surrogate index for other land uses' )

        MXSRG = MAX( MXSRG, OSRG )


C.......   Loop:  read emission factor file

        WRITE( LDEV,92000 ) ' ', 'Reading EMISSIONS FACTOR file', ' '
        I     = 0
        EFLAG = .FALSE.
33      CONTINUE

            I = I + 1
            READ( FDEV, 93010, END=34, IOSTAT=ISTAT )
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

            IF ( I .LT. MXID ) GO TO 33 ! end loop reading emissions factor file

        CALL M3EXIT( 'RAWBIO', 0, 0, 
     &               'EMISSIONS FACTOR table overflow', 2 )

34      CONTINUE

        NVEG = I
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'RAWBIO', 0, 0,
     &                   'Error reading EMISSIONS FACTOR file', 2 )
        END IF

C...........   Fold ug~~>g, hectare~~>m^2 factors into emfac:

        DO  36  J = 1, NSEF
        DO  35  I = 1, NVEG
            EMFAC( I, J ) = MICR2G * HA2MSQ * EMFAC( I, J )
35      CONTINUE
36      CONTINUE


C.......   Loop:  read FIPS list from FIPS CODE file

        WRITE( LDEV,92000 ) ' ', 'Reading FIPS CODE file', ' '
        I     = 0
        EFLAG = .FALSE.

C.......   Skip the first line in FIPS CODE file
        READ( CDEV,'(A20)', END = 45, IOSTAT= ISTAT ) CO_NAME

44      CONTINUE                !  head of loop:  read FIPS CODE file

            READ( CDEV,92010, END = 45, IOSTAT= ISTAT ) CO_NAME, FIP,
     &                                  DUMMY1,DUMMY2,DUMMY3    ! sl
            IF ( ISTAT .EQ. 0 ) THEN    ! sl modify .NE. to .EQ.
                I = I + 1
                IF ( I .LE. NBFIP ) THEN
                    FIPS( I ) = FIP
                ELSE
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'FIPS buffer overflow at line', I
                    CALL M3MESG( MESG )
                END IF
                GO TO 44
            END IF

45      CONTINUE                !  exit from loop:  read FIPS CODE file
        NFIP = I                ! sl NFIP was not assigned anywhere    

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'RAWBIO', 0, 0, 
     &                   'Error reading FIPS CODE file', 2 )
        END IF


C.......   Loop:  read county land use file:

        WRITE( LDEV,92000 ) ' ', 'Reading COUNTY LAND USE file', ' '
        LINE  = 0
        EFLAG = .FALSE.

55      CONTINUE                !  head of loop:  read COUNTY LAND USE file

            READ( UDEV, 93000, END = 122, IOSTAT=ISTAT ) INBUF
            LINE =  LINE + 1
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading FIP from COUNTY LAND USE at line', LINE
                CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
            END IF

            FIP = STR2INT( INBUF )
c sl            CTY = FIND1( FIP, NFIP, FIPS )
            CTY = FINDFIP( FIP, NFIP, FIPS )

C...........   Initialize accumulator variables for this county:

            IF ( CTY .GT. 0 ) THEN

                DO  56  M = 1, BSPCS - 1

                    EPINE( M ) = 0.0
                    EDECD( M ) = 0.0
                    ECONF( M ) = 0.0
                    EOTH ( M ) = 0.0
                    ELAI ( M ) = 0.0 
                    EAG  ( M ) = 0.0

56              CONTINUE

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

            END IF

C...........   Land use type:  (rural) forest.  Process subtypes:

            READ( UDEV, 93000, IOSTAT=ISTAT ) INBUF
            LINE =  LINE + 1
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading FOREST AREA from ' //
     &          'COUNTY LAND USE at line', LINE
                CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
            END IF

            I     = INDEX( INBUF, ',' )
            ATYPE = STR2REAL( INBUF(   1 : I-1 ) )
            NTYPE = STR2INT ( INBUF( I+1 :  80 ) )

            DO  69  N = 1, NTYPE

                READ( UDEV, 93000, IOSTAT=ISTAT ) INBUF
                LINE =  LINE + 1
                IF ( ISTAT .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'COUNTY LAND USE at line', LINE
                    CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
                END IF

                J     = LBLANK( INBUF )
                I     = INDEX( INBUF( J+1:80 ) , ',' )
                TYPE  = INBUF( J+1 : J+I-1 )
                AREA  = STR2REAL( INBUF( J+I+1 : 80 ) )
                K     = INDEX1( TYPE, NVEG, VEGID )

                IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LU FILE in VEGID at line', LINE
                    CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
                END IF

                L = INDEX1( TYPE, SPTREE, SPFORID )

                IF ( L .GT. 0 ) THEN
                    AFOR   = AFOR   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    ALAI   = ALAI   + AREA
                    DO  61  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
61                  CONTINUE
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Wdcp' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AAG    = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  62  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
62                  CONTINUE
                    ENOAG  = ENOAG + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AGRS   = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  63  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
63                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2*AREA
                    AGRS   = AAG    + 0.2*AREA
                    AOTH   = AAG    + 0.6*AREA
                    DO  64  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
64                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                ELSE    !  Add area to the forest area, add NO to the forest NO

                    AFOR = AFOR + AREA
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                    IF     ( LAI( K ) .EQ. 3 ) THEN
                        DO  65  M = 1, BSPCS - 1
                            EPINE( M ) = EPINE( M ) + AREA*EMFAC( K,M )
65                      CONTINUE

                    ELSE IF( LAI( K ) .EQ. 5 ) THEN

                        DO  66  M = 1, BSPCS - 1
                            EDECD( M ) = EDECD( M ) + AREA*EMFAC( K,M )
66                      CONTINUE

                    ELSE IF( LAI( K ) .EQ. 7 ) THEN

                        DO  67  M = 1, BSPCS - 1
                            ECONF( M ) = ECONF( M ) + AREA*EMFAC( K,M )
67                      CONTINUE

                    ELSE

                        SUMLAI = SUMLAI + AREA * LAI( K )
                        ALAI   = ALAI   + AREA
                        DO  68  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA*EMFAC( K,M )
68                      CONTINUE

                    END IF      !  if lai is 3,5,7, or otherwise

                END IF  !  if some spforid, or 'Wdcp' or 'Scwd' or 'Urba or not

69          CONTINUE


C...........   Land use type:  urban forest.  Process subtypes:

            READ( UDEV, 93000, IOSTAT=ISTAT ) INBUF
            LINE =  LINE + 1
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading URBAN FOREST from ' //
     &          'COUNTY LAND USE at line', LINE
                CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
            END IF

            I     = INDEX( INBUF, ',' )
            ATYPE = STR2REAL( INBUF(   1 : I-1 ) )
            NTYPE = STR2INT ( INBUF( I+1 :  80 ) )

            DO  79  N = 1, NTYPE

                READ( UDEV, 93000, IOSTAT=ISTAT ) INBUF
                LINE =  LINE + 1
                IF ( ISTAT .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'COUNTY LAND USE at line', LINE
                    CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
                END IF

                J     = LBLANK( INBUF )
                I     = INDEX( INBUF( J+1:80 ) , ',' )
                TYPE  = INBUF( J+1 : J+I-1 )
                AREA  = STR2REAL( INBUF( J+I+1 : 80 ) )

                K     = INDEX1( TYPE, NVEG, VEGID )

                IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LU FILE in VEGID at line', LINE
                    CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
                END IF

                L = INDEX1( TYPE, SPTREE, SPFORID )

                IF ( L .GT. 0 ) THEN

                    AFOR   = AFOR   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    ALAI   = ALAI   + AREA
                    DO  71  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
71                  CONTINUE
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Wdcp' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AAG    = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  72  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
72                  CONTINUE
                    ENOAG  = ENOAG + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AGRS   = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  73  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
73                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2*AREA
                    AGRS   = AAG    + 0.2*AREA
                    AOTH   = AAG    + 0.6*AREA
                    DO  74  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
74                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                ELSE    !  Add area to the forest area, add NO to the forest NO

                    AFOR = AFOR + AREA
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                    IF     ( LAI( K ) .EQ. 3 ) THEN
                        DO  75  M = 1, BSPCS - 1
                            EPINE( M ) = EPINE( M ) + AREA*EMFAC( K,M )
75                      CONTINUE

                    ELSE IF( LAI( K ) .EQ. 5 ) THEN

                        DO  76  M = 1, BSPCS - 1
                            EDECD( M ) = EDECD( M ) + AREA*EMFAC( K,M )
76                      CONTINUE

                    ELSE IF( LAI( K ) .EQ. 7 ) THEN

                        DO  77  M = 1, BSPCS - 1
                            ECONF( M ) = ECONF( M ) + AREA*EMFAC( K,M )
77                      CONTINUE

                    ELSE

                        SUMLAI = SUMLAI + AREA * LAI( K )
                        ALAI   = ALAI   + AREA
                        DO  78  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA*EMFAC( K,M )
78                      CONTINUE

                    END IF      !!  if lai is 3,5,7, or otherwise

                END IF  !  if some spforid, or 'Wdcp' or 'Scwd' or 'Urba or not

79          CONTINUE


C...........   Land use type:  agriculture.  Process subtypes:

            READ( UDEV, 93000, IOSTAT=ISTAT ) INBUF
            LINE =  LINE + 1
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading AGRICULTURE from ' //
     &          'COUNTY LAND USE at line', LINE
                CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
            END IF

            I     = INDEX( INBUF, ',' )
            ATYPE = STR2REAL( INBUF(   1 : I-1 ) )
            NTYPE = STR2INT ( INBUF( I+1 :  80 ) )

            DO  82  N = 1, NTYPE

                READ( UDEV, 93000, IOSTAT=ISTAT ) INBUF
                LINE =  LINE + 1
                IF ( ISTAT .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'COUNTY LAND USE at line', LINE
                    CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
                END IF

                J     = LBLANK( INBUF )
                I     = INDEX( INBUF( J+1:80 ) , ',' )
                TYPE  = INBUF( J+1 : J+I-1 )
                AREA  = STR2REAL( INBUF( J+I+1 : 80 ) )

                K     = INDEX1( TYPE, NVEG, VEGID )

                IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LU FILE in VEGID at line', LINE
                    CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
                END IF

                AAG = AAG + AREA
                DO  81  M = 1, BSPCS - 1
                    EAG( M ) = EAG( M ) + AREA * EMFAC( K , M )
81              CONTINUE
                ENOAG  = ENOAG + AREA * EMFAC( K, NO )

82          CONTINUE


C...........   Land use type:  remaining/other.  Process subtypes:

            READ( UDEV, 93000, IOSTAT=ISTAT ) INBUF
            LINE =  LINE + 1
            IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading OTHER AREA from ' //
     &          'COUNTY LAND USE at line', LINE
                CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
            END IF

            I     = INDEX( INBUF, ',' )
            ATYPE = STR2REAL( INBUF(   1 : I-1 ) )
            NTYPE = STR2INT ( INBUF( I+1 :  80 ) )

            DO  101  N = 1, NTYPE

                READ( UDEV, 93000, IOSTAT=ISTAT ) INBUF
                LINE =  LINE + 1
                IF ( ISTAT .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'COUNTY LAND USE at line', LINE
                    CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
                END IF

                J     = LBLANK( INBUF )
                I     = INDEX( INBUF( J+1:80 ) , ',' )
                TYPE  = INBUF( J+1 : J+I-1 )
                AREA  = STR2REAL( INBUF( J+I+1 : 80 ) )

                K     = INDEX1( TYPE, NVEG, VEGID )

                IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LU FILE in VEGID at line', LINE
                    CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
                END IF

                L = INDEX1( TYPE, RMTREE, OTHERID )

                IF ( L .GT. 0 ) THEN

                    ALAI   = ALAI + AREA
                    AFOR   = AFOR + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  91  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
91                  CONTINUE
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Pacp' ) THEN

                    AAG = AAG + AREA
                    DO  92  M = 1, BSPCS - 1
                        EAG( M ) = EAG( M ) + AREA * EMFAC( K , M )
92                  CONTINUE
                    ENOAG = ENOAG + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Gras' .OR.

     &                    TYPE .EQ. 'Scru' .OR.
     &                    TYPE .EQ. 'Ugra' .OR.
     &                    TYPE .EQ. 'Othe' ) THEN
                    AGRS = AGRS + AREA
                    DO  93  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
93                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Wetf' ) THEN

                    AWTF = AWTF + AREA
                    ALAI = ALAI + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  94  M = 1, BSPCS - 1
                        EMIS = AREA * EMFAC( K , M )
                        ELAI( M ) = ELAI( M ) + EMIS
94                  CONTINUE
                    ENOWTF = ENOWTF + AREA * EMFAC( K, NO )

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
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
95                  CONTINUE
                    ENOAG = ENOAG + AREA * EMFAC( K, NO )
                ELSE IF ( TYPE .EQ. 'Scwd' ) THEN
                    ALAI   = ALAI   + AREA
                    AFOR   = AFOR   + 0.5 * AREA
                    AGRS   = AGRS   + 0.5 * AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  96  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
96                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2 * AREA
                    AGRS   = AGRS   + 0.2 * AREA
                    AOTH   = AOTH   + 0.6 * AREA
                    DO  97  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
97                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                ELSE IF ( TYPE .EQ. 'Desh' ) THEN

                    AGRS   = AGRS   + 0.5 * AREA
                    AOTH   = AOTH   + 0.5 * AREA
                    DO  98  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
98                  CONTINUE
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                ELSE

                    AOTH = AOTH + AREA

                    IF ( LAI( K ) .GT. 0 ) THEN

                        ALAI   = ALAI   + AREA
                        SUMLAI = SUMLAI + AREA * LAI( K )
                        DO  99  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA * EMFAC( K,M )
99                      CONTINUE
                        ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                    ELSE

                        DO  100  M = 1, BSPCS - 1
                            EOTH( M ) = EOTH( M ) + AREA * EMFAC( K,M )
100                      CONTINUE
                        ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                    END IF

                END IF

101        CONTINUE


C...........   Save only counties found in the county-list:

           IF ( CTY .GT. 0 ) THEN

               IF ( ALAI .GT. 0.0 ) THEN
                   AVGLAI = SUMLAI / ALAI 
               ELSE
                   AVGLAI = 0.0
               END IF

               DO 111 M = 1, BSPCS - 1

                   CTYPINE( M,CTY ) = EPINE( M )
                   CTYDECD( M,CTY ) = EDECD( M )
                   CTYCONF( M,CTY ) = ECONF( M )
                   CTYAGRI( M,CTY ) = EAG  ( M )
                   CTYLEAF( M,CTY ) = ELAI ( M )
                   CTYOTHR( M,CTY ) = EOTH ( M )

111            CONTINUE           

               CTYNOFOR( CTY ) = ENOFOR
               CTYNOGRS( CTY ) = ENOGRS
               CTYNOAG ( CTY ) = ENOAG
               CTYNOWTF( CTY ) = ENOWTF
               CTYAVLAI( CTY ) = AVGLAI
               CTYAFOR ( CTY ) = AFOR
               CTYAAG  ( CTY ) = AAG
               CTYAGRS ( CTY ) = AGRS
               CTYAWTF ( CTY ) = AWTF
               CTYAOTH ( CTY ) = AOTH

           END IF

           GO TO  55            !  to head of loop:  read COUNTY LAND USE file


122     CONTINUE                !  exit from loop:  read COUNTY LAND USE file



C.......   Read and process SURROGATE FACTORS file:

        WRITE( LDEV,92000 ) ' ', 'Reading SURROGATE file', ' '

        DO 203 I = 1, BSPCS-1
        DO 202 R = 1, NROWS
        DO 201 C = 1, NCOLS
            PINE( C, R, I ) = 0.0
            DECD( C, R, I ) = 0.0
            CONF( C, R, I ) = 0.0
            AGRC( C, R, I ) = 0.0
            LEAF( C, R, I ) = 0.0
            OTHR( C, R, I ) = 0.0
201     CONTINUE
202     CONTINUE
203     CONTINUE
 
        DO 205 R = 1, NROWS
        DO 204 C = 1, NCOLS
            AVLAI ( C, R ) = 0.0
            GRASNO( C, R ) = 0.0
            FORENO( C, R ) = 0.0
            WETLNO( C, R ) = 0.0
            AGRINO( C, R ) = 0.0
204     CONTINUE
205     CONTINUE

C.......   Read file header for spatial allocation factor file:

        READ( SDEV, 93020, IOSTAT = ISTAT )
     &          BUTME, BUTMN, EUTME, EUTMN, UTMZONE, XCELSZ, YCELSZ

        XX = ( XORIG3D - BUTME ) / XCELSZ
        YY = ( YORIG3D - BUTMN ) / YCELSZ
        IF ( ISTAT .NE. 0 ) THEN
            WRITE( MESG,94010 ) 
     &      'Error', ISTAT, 'reading header of SURROGATE file'
            CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
        ELSE IF ( FLTERR( XCELSZ, SNGL( XCELL3D ) ) ) THEN
            CALL M3EXIT( 'RAWBIO', 0, 0,
     &                   'Bad XCELSZ in SURROGATE file', 2 )
        ELSE IF ( FLTERR( YCELSZ, SNGL( YCELL3D ) ) ) THEN
            CALL M3EXIT( 'RAWBIO', 0, 0,
     &                   'Bad YCELSZ in SURROGATE file', 2 )
        ELSE IF ( FLTERR( XX, FLOAT( NINT( XX ) ) ) ) THEN
            CALL M3EXIT( 'RAWBIO', 0, 0,
     &                   'Bad BUTME alignment in SURROGATE file', 2 )
        ELSE IF ( FLTERR( YY, FLOAT( NINT( YY ) ) ) ) THEN
            CALL M3EXIT( 'RAWBIO', 0, 0,
     &                   'Bad BUTMN alignment in SURROGATE file', 2 )
        ELSE IF ( XX .LT. -1.0E-4 * XCELSZ ) THEN
            CALL M3EXIT( 'RAWBIO', 0, 0,
     &                   'Bad BUTME placement in SURROGATE file', 2 )
        ELSE IF ( YY .LT. -1.0E-4 * YCELSZ ) THEN
            CALL M3EXIT( 'RAWBIO', 0, 0,
     &                   'Bad BUTMN placement in SURROGATE file', 2 )
c sl        ELSE IF ( EUTME .LT. 
c    &            XORIG3D + ( NCOLS3D + 1.0E-4 ) * XCELSZ ) THEN
c           CALL M3EXIT( 'RAWBIO', 0, 0,
c    &                   'Bad EUTME placement in SURROGATE file', 2 )
c       ELSE IF ( EUTMN .LT.
c    &            YORIG3D + ( NROWS3D + 1.0E-4 ) * YCELSZ ) THEN
c           CALL M3EXIT( 'RAWBIO', 0, 0,
c sl     &                   'Bad EUTMN placement in SURROGATE file', 2 )
        END IF

        DWIDTH  = 1.0 / XCELSZ
        DHEIGHT = 1.0 / YCELSZ
        LINE    = 1
        EFLAG   = .FALSE.

211     CONTINUE                !  head of loop:  read SURROGATE FACTORS file

            READ( SDEV, 93030, IOSTAT = ISTAT, END = 299 )
     &          FIP,
     &          XCOOR,
     &          YCOOR,
     &          ZONE,
     &        ( FRAC(K), K = 1, MXSRG )
            LINE = LINE + 1

            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 
     &          'Error ', ISTAT, 
     &          'reading SURROGATE FACTORS file at line', LINE
                CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
            END IF
            
c sl            K = FIND1( FIP, NFIP, FIPS )
            K = FINDFIP( FIP, NFIP, FIPS )
            IF ( K .LE. 0 ) THEN
                EFLAG = .FALSE.
                WRITE( MESG, 94010 ) 
     &          'FIPS code ', FIP, 
     &           'at line', LINE, 'not found in FIP list'
                CALL M3MESG( MESG )
                GO TO  211
            END IF

            C = NINT( ( XCOOR - SNGL(XORIG3D) ) * DWIDTH  ) + 1
            R = NINT( ( YCOOR - SNGL(YORIG3D) ) * DHEIGHT ) + 1

            IF ( C .LT. 1  .OR.  C .GT. NCOLS  .OR.
     &           R .LT. 1  .OR.  R .GT. NROWS ) THEN
                WRITE( MESG, 94010 )
     &          'Cell (C,R)=', C, ',', R, 'outside grid at line', LINE
                CALL M3MESG( MESG )
            END IF

            DO  222  M = 1, BSPCS - 1
                PINE( C,R,M ) = PINE( C,R,M ) + 
     &                          CTYPINE( M,K ) * FRAC( PSRG )
                DECD( C,R,M ) = DECD( C,R,M ) + 
     &                          CTYDECD( M,K ) * FRAC( DSRG )
                CONF( C,R,M ) = CONF( C,R,M ) + 
     &                          CTYCONF( M,K ) * FRAC( CSRG )
                AGRC( C,R,M ) = AGRC( C,R,M ) + 
     &                          CTYAGRI( M,K ) * FRAC( ASRG )
                LEAF( C,R,M ) = LEAF( C,R,M ) + 
     &                          CTYLEAF( M,K ) * FRAC( FSRG )
                OTHR( C,R,M ) = OTHR( C,R,M ) + 
     &                          CTYOTHR( M,K ) * FRAC( OSRG )
222         CONTINUE
            GRASNO( C,R ) = GRASNO( C,R ) + CTYNOGRS( K ) * FRAC( GSRG )
            FORENO( C,R ) = FORENO( C,R ) + CTYNOFOR( K ) * FRAC( FSRG )
            WETLNO( C,R ) = WETLNO( C,R ) + CTYNOWTF( K ) * FRAC( WSRG )
            AGRINO( C,R ) = AGRINO( C,R ) + CTYNOAG ( K ) * FRAC( ASRG )

            AVLAI ( C,R ) = AVLAI ( C,R ) + CTYAVLAI( K ) * FRAC( FSRG )
                
            GO TO 211           ! end of loop reading SURROGATE FACTORS file


299     CONTINUE                ! exit from loop reading SURROGATE FACTORS file

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'RAWBIO', 0, 0,
     &                   'Error reading SURROGATE FACTORS file', 2 )
        END IF


C...............   Write output file:

        I = 0        
        DO  312  B = 1, BTYPES
        DO  311  M = 1, BSPCS - 1

            I = I + 1
            IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                         NORM( 1,1,M,B ) ) ) THEN
                MESG = 'Could not write "' //
     &                  VNAME3D( I )( 1: TRIMLEN( VNAME3D( I ) ) ) //
     &                  '" to ' // ENAME
                CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
            END IF

311     CONTINUE        !  end loop on VOC species M
312     CONTINUE        !  end loop on VOC-emitting types B

        IF ( .NOT. WRITE3( ENAME, 'AVLAI', 0, 0,
     &                     AVLAI ) ) THEN
            MESG = 'Could not write "AVLAI" to ' // ENAME
            CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'GRASNO', 0, 0,
     &                     GRASNO ) ) THEN
            MESG = 'Could not write "GRASNO"to ' // ENAME
            CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'FORENO', 0, 0,
     &                     FORENO ) ) THEN
            MESG = 'Could not write "FORENO"to ' // ENAME
            CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'WETLNO', 0, 0,
     &                     WETLNO ) ) THEN
            MESG = 'Could not write "WETLNO"to ' // ENAME
            CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
        END IF

        IF ( .NOT. WRITE3( ENAME, 'AGRINO', 0, 0,
     &                     AGRINO ) ) THEN
            MESG = 'Could not write "AGRINO"to ' // ENAME
            CALL M3EXIT( 'RAWBIO', 0, 0, MESG, 2 )
        END IF


C.......   End of program:

      CALL M3EXIT( 'RAWBIO', 0, 0, 
     &             'Successful completion of program RAWBIO', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx
C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )

92010   FORMAT ( A20, I5, 3F11.3 )

C...........   Formatted file I/O formats............ 93xxx
                                   
93000   FORMAT( A )

93010   FORMAT( 1X, A4, 3 F9.0, F6.0, I2 )

93020   FORMAT( 4 F10.0, I10, 2F10.0 )

93030   FORMAT( I5, 2F10.0, I10, 100F10.4 )    ! sl For surrogate data


C...........   Internal buffering formats............ 94xxx
94010   FORMAT( 10 ( A, :, I5, :, 2X ) )


C...........   Miscellaneous formats................. 95xxx

95000   FORMAT ( /5X , A , $ )          !  generic prompt format.


        END


        INTEGER FUNCTION FINDFIP( K, N, LIST)

C***********************************************************************
C  RETURNS:
C       subscript at which the targeted key appears, or
C       -1 in case of failure.
C***********************************************************************

        INTEGER        K
        INTEGER        N
        INTEGER        LIST( N )

        DO I = 1, N
           IF ( K .EQ. LIST( I ) ) THEN
              FINDFIP = I
              RETURN
           END IF
        END DO 
       
        FINDFIP = -1
        RETURN
        END

