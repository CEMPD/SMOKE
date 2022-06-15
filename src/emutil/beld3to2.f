
        PROGRAM BELD3TO2

C***********************************************************************
C
C  DESCRIPTION: Converts IOAPI netCDF landuse files created from SMOKE
C               tool into BELD2 gridded ASCII format.  The SMOKE tool
C               takes BELD3 1km resolution data and creates landuse
C               files for user-defined domain for use with BEIS3.  BELD3
C               databases have over 230 landuse types.  This database
C               cannot be used for SMOKE-BEIS2 modeling since this model
C               uses BELD2 format data that has 126 landuse types.  This
C               program maps the 230 BELD3 types to the 126 BELD2 landuse 
C               types and outputs an ASCII format file that SMOKE-BEIS2 
C               can use.  There are 4 BELD2 categories that are output :
C                     1: Urban Forest
C                     2: Rural Forest
C                     3: Agricultural 
C                     4: Other
C
C               SMOKE-BEIS2 treats Urban and Rural Forest the same.  So,
C               the program outputs all forest land use types for a grid
C               cell into the first category (Urban Forest).  The BELD3
C               to BELD2 xref file includes the BELD2 category for each
C               BELD3 landuse types. 
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C              02/01 : prototype created by J. Vukovich
C                      tested on Lambert projection
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
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
        INCLUDE 'EMCNST3.EXT'     !

C...........   PARAMETERS and their descriptions:

C...........   LOCAL PARAMETERS

        CHARACTER(50), PARAMETER :: 
     &  CVSW = '$Name SMOKEv4.9_Jun2022$' ! CVS release tag
        REAL,          PARAMETER :: TOHECT = 0.0001  ! factor to convert to hectares
        CHARACTER,     PARAMETER :: AQUT = "'"  
        INTEGER,       PARAMETER :: IZERO = 0

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL         ENVYN
        INTEGER         GETFLINE
        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE
        INTEGER         TRIMLEN
        CHARACTER(16)   VERCHAR

        EXTERNAL        ENVYN, GETFLINE, PROMPTFFILE, PROMPTMFILE,  
     &                  TRIMLEN, VERCHAR

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER         B, C, R, I, J, K, L, M, N ! loop counters and subscripts
        INTEGER         IOS     !  I/O status result
        INTEGER         JCAT    !  temporary BELD2 category
        INTEGER         FDEV    !  unit number for BELD3 to BELD2 xref file
        INTEGER         LDEV    !  unit number for log file:
        INTEGER         XDEV    !  BELD2 output file

        CHARACTER(256)  MESG    !  message buffer for M3EXIT()
        CHARACTER(16), ALLOCATABLE  :: LUNAMA( : )   ! Landuse file A var. names
        CHARACTER(16), ALLOCATABLE  :: LUNAMB( : )   ! Landuse file B var. names
        CHARACTER(16), ALLOCATABLE  :: VEG3( : )     ! BELD3 landuse types
        CHARACTER(16)   LUNAM   !  temporary landuse name
        CHARACTER(16)   GRDNM   !  grid name
        CHARACTER(16)   GNAMEA  !  unit number for gridded land use file A
        CHARACTER(16)   GNAMEB  !  unit number for gridded land use file B
        CHARACTER(16)   GNAMET  !  unit number for gridded land use totals file
        CHARACTER(16)   COORD   !  coordinate system used
        CHARACTER(8)    CUNIT   !  units

        CHARACTER(4), ALLOCATABLE :: VEG2( : )      ! BELD2 xref record
        CHARACTER(4), ALLOCATABLE :: B2TYPES ( : )  ! unsorted BELD2 types
        CHARACTER(4), ALLOCATABLE :: B2SORT( : )    ! sorted BELD2 types
        CHARACTER(4)                 TMPVEG         ! temporary veg type

        INTEGER         NCOLS   ! no. of grid columns
        INTEGER         NGRID   ! no. of grid cells
        INTEGER         NROWS   ! no. of grid rows
        INTEGER         NVARSA  ! no. of variables in A file
        INTEGER         NVARSB  ! no. of variables in B file
        INTEGER         NVARST  ! total of variables in A and B
        INTEGER         IFOUND  ! temporary variable to help with indexing
        INTEGER         IXREFS  ! no. of records in xref file
        INTEGER         IB2     ! no. of BELD2 landuse types to output
        INTEGER         IUSDEC, IUSBRD, IUSCON, IUSMIX  
        INTEGER         JFLAG   ! flag for counting BELD2 types

        INTEGER, ALLOCATABLE :: ICAT  ( : )      ! BELD2 category flag
        INTEGER, ALLOCATABLE :: IFLAG ( : )      ! flag for aggregration
        INTEGER, ALLOCATABLE :: IB2CAT( : )      ! 
        INTEGER, ALLOCATABLE :: LUINDX( : )      ! landuse name index
        INTEGER, ALLOCATABLE :: NLUS( :, :, : )  ! no. BELD2 landuse per category
        INTEGER, ALLOCATABLE :: INDXB2( : )      ! index sort of BELD2 types

        REAL, ALLOCATABLE    :: LUSE ( :, :, : ) ! BELD3 landuse
        REAL, ALLOCATABLE    :: LUB2 ( :, :, : ) ! BELD2 landuse
        REAL, ALLOCATABLE    :: TCAT ( :, :, : ) ! Total landuse per BELD2 category
        REAL, ALLOCATABLE    :: TOTAL( :, : )    ! Total landuse per grid cell
        REAL, ALLOCATABLE    :: FIA ( :, : )     ! Forest inv % in each cell
        REAL  XRES, YRES, BTMP, TTMP

        LOGICAL         EFLAG                    ! error flag

        CHARACTER(16) :: PROGNAME = 'BELD3TO2'    !  program name

C***********************************************************************
C   begin body of program NORMBIO

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.

        CALL INITEM( LDEV, CVSW, PROGNAME )

C.......   Open BELD3 to BELD2 xref file and BELD2 landuse names file

        FDEV = PROMPTFFILE( 
     &           'Enter logical name for Land Use xref file',
     &           .TRUE., .TRUE., 'B3XRF', PROGNAME )


C.......   Loop:  read BELD3 to BELD2 xref file

        CALL M3MSG2( 'Reading BELD3to2 XREF file' )

        IXREFS = GETFLINE( FDEV, 'BELD3 to BELD2 xref file ' )

C.......  Allocate memory 

        ALLOCATE( VEG2 ( IXREFS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VEG2', PROGNAME )
        ALLOCATE( VEG3 ( IXREFS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VEG3', PROGNAME )
        ALLOCATE( ICAT ( IXREFS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ICAT', PROGNAME )
        ALLOCATE( B2TYPES ( IXREFS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'B2TYPES', PROGNAME )
        ALLOCATE( IB2CAT ( IXREFS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IB2CAT', PROGNAME )
  
C........   Read in BELD3 to BELD2 xref and keep
C........   track of number of BELD2 types to output.

        IB2 = 0  
        JFLAG = 0

        DO I = 1, IXREFS

            READ( FDEV, 93010 ) VEG2( I ), JCAT, VEG3( I )

C.........  Assign BELD2 category for each output type

            IF ( JCAT .GT. 4 ) THEN
                ICAT( I ) = 1
            ELSE
                ICAT( I ) = JCAT
            ENDIF

C..........  Check to see if this is new BELD2 type

            DO J = 1, IB2
                IF ( VEG2( I ) .EQ. B2TYPES( J ) ) JFLAG = 1
            ENDDO

C..........  If it is a new type assign values to appropriate arrays

            IF ( JFLAG .EQ. 0 ) THEN
               IB2 = IB2 + 1
               B2TYPES( IB2 ) = VEG2( I ) 
               IB2CAT( IB2 ) = ICAT( I )
            ELSE
               JFLAG = 0
            ENDIF

        ENDDO

C........... Allocate memory and initialize

        ALLOCATE( INDXB2 ( IB2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXB2', PROGNAME )
        ALLOCATE( B2SORT ( IB2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXB2', PROGNAME )
        INDXB2 = 0 ! array

C........... Setup for alphabetically sorting BELD2 types

        DO I = 1, IB2
            B2SORT ( I ) = B2TYPES( I ) 
            INDXB2 ( I ) = I
        ENDDO

C...........  Sort BELD2 types alphabetically

        CALL SORTIC ( IB2, INDXB2 , B2SORT ) 

C.......   Open gridded landuse files 

        GNAMET = PROMPTMFILE(
     &           'Enter logical name for GRIDDED LANDUSE totals file',
     &           FSREAD3, 'BELD3_TOT', PROGNAME )

        IF ( .NOT. DESC3( GNAMET ) ) THEN

            MESG = 'Could not get description of file "' //
     &             GNAMET( 1:TRIMLEN( GNAMET ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

C.......   Initialize grid definition

        CALL CHKGRID( GNAMET, 'GRID' , 0, EFLAG )

        GNAMEA = PROMPTMFILE(
     &           'Enter logical name for GRIDDED LANDUSE A file',
     &           FSREAD3, 'BELD3_A', PROGNAME )

        IF ( .NOT. DESC3( GNAMEA ) ) THEN

            MESG = 'Could not get description of file "' //
     &             GNAMEA( 1:TRIMLEN( GNAMEA ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

C.......   Check grid definition

        CALL CHKGRID( GNAMEA, 'GRID' , 0, EFLAG )

C........  If grid definition does not match first landuse file then stop

        IF ( EFLAG ) THEN
            MESG = 'Problems opening input files. See ERROR(S) above.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF


        NVARSA = NVARS3D

C.......   Resolution variables assumed to be in meters


        NCOLS = NCOLS3D
        NROWS = NROWS3D
        XRES  = XCELL3D 
        YRES  = YCELL3D  

        ALLOCATE( LUNAMA( NVARSA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUNAMA', PROGNAME )

        DO I = 1, NVARSA

            LUNAMA ( I ) = VNAME3D ( I ) 

        ENDDO
 
        GNAMEB = PROMPTMFILE(
     &           'Enter logical name for GRIDDED LANDUSE B file',
     &           FSREAD3, 'BELD3_B', PROGNAME )

        IF ( .NOT. DESC3( GNAMEB ) ) THEN

            MESG = 'Could not get description of file "' //
     &             GNAMEB( 1:TRIMLEN( GNAMEB ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

C.......   Check grid definition

        CALL CHKGRID( GNAMEB, 'GRID' , 0, EFLAG )

C........  If grid definition does not match first landuse file then stop

        IF ( EFLAG ) THEN
            MESG = 'Problems opening input files. See ERROR(S) above.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF


        NVARSB = NVARS3D

        ALLOCATE( LUNAMB( NVARSB ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUNAMB', PROGNAME )

        DO I = 1, NVARSB
            
            LUNAMB ( I ) = VNAME3D ( I ) 

        ENDDO
     
C............. some grid description info obtained from call to RDSRG above

        XDEV = PROMPTFFILE(
     &           'Enter logical name for output BELD2 LANDUSE file',
     &           .FALSE., .TRUE., 'BGUSE', PROGNAME )


        NVARST = NVARSA + NVARSB


C.............. prepare info to output for header of BELD2 output file
C.............. this header will contain grid information
 
        CUNIT = 'meters  '
        IF ( GDTYP3D .EQ. 1 ) THEN 
           COORD = 'LAT-LON '
           CUNIT = 'degrees '
        ELSE IF ( GDTYP3D .EQ. 2 ) THEN
           COORD = 'LAMBERT '
        ELSE
           MESG = 'Code not tested for this grid type!'
           CALL M3WARN( PROGNAME, 0, 0, MESG )          
        ENDIF

        WRITE( XDEV, 93100 ) GDNAM3D,XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &                       NCOLS3D, NROWS3D, IZERO, COORD, CUNIT,
     &                       P_ALP3D, P_BET3D,
     &                       P_GAM3D, XCENT3D, YCENT3D 

C.............  Allocate memory for BELD3 land use and index

        ALLOCATE( LUSE ( NCOLS, NROWS, NVARST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUSE', PROGNAME )
        ALLOCATE( LUINDX ( NVARST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUINDX', PROGNAME )

        LUINDX = -9   ! array

C.......  Check to see if all BELD3 types in xref are available in
C.......  xref file

        DO I = 1, NVARST

           IFOUND = 0

           IF ( I .LE. NVARSA ) THEN
              LUNAM = LUNAMA( I ) 
           ELSE
              K = I - NVARSA
              LUNAM = LUNAMB( K ) 
           ENDIF 
 
           DO J = 1, IXREFS

             IF ( VEG3 ( J ) .EQ. LUNAM ) THEN
                IFOUND = 1  
                LUINDX( I ) = J
             ENDIF

C.........  Assign landuse type value for 4 USGS classes to be used later

             IF ( VEG3( J ) .EQ. 'USGS_decidforest' ) IUSDEC = J
             IF ( VEG3( J ) .EQ. 'USGS_evbrdleaf  ' ) IUSBRD = J
             IF ( VEG3( J ) .EQ. 'USGS_coniferfor ' ) IUSCON = J
             IF ( VEG3( J ) .EQ. 'USGS_mxforest   ' ) IUSMIX = J

           ENDDO

           IF ( IFOUND .EQ. 0 ) THEN

               MESG =   LUNAM // 
     &                ' does NOT have xref record in B3TO2 xref file'
               CALL M3WARN( PROGNAME, 0, 0, MESG )

           ENDIF 
         
        ENDDO


C.......  Read the gridded landuse from the landuse files
       
        DO M = 1, NVARST 
            N = LUINDX ( M )
            IF ( N .GT. 0 ) THEN
               IF ( M .LE. NVARSA ) THEN
        
                  MESG = ' Reading ' // LUNAMA( M )
                  CALL M3MESG( MESG )
                  IF ( .NOT. READ3( GNAMEA, LUNAMA( M ) , 1, 0, 0 , 
     &                LUSE( 1, 1, N ) ) ) THEN
 
                    MESG = 'Could not find ' // LUNAMA( M ) // 
     &                  ' in file ' // GNAMEA( 1:TRIMLEN( GNAMEA ) )

                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                  ENDIF          !  if read3() failed
 
               ELSE

                  K = M - NVARSA
                  MESG = ' Reading ' // LUNAMB( K )
                  CALL M3MESG( MESG )

                  IF ( .NOT. READ3( GNAMEB, LUNAMB( K ) , 1, 0, 0 , 
     &                LUSE( 1, 1, N ) ) ) THEN

                     MESG = 'Could not find ' // LUNAMB( K ) // 
     &                   ' in file '// GNAMEB( 1:TRIMLEN( GNAMEB ) )

                     CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                  ENDIF          !  if read3() failed

               ENDIF
            ENDIF 

        ENDDO

C......... Read FIA percentage from landuse totals file

        ALLOCATE( FIA( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FIA', PROGNAME )

        MESG = ' Reading FIA from totals file ' 
        CALL M3MESG( MESG )
        IF ( .NOT. READ3( GNAMET, 'FOREST', 1, 0, 0, FIA ) ) THEN
             MESG = 'Could not read FOREST from file "' //
     &                GNAMET( 1:TRIMLEN( GNAMET ) ) // '"'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF


C.......   Allocate memory for output purposes

        ALLOCATE( LUB2 ( NCOLS, NROWS, IB2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUB2', PROGNAME )
        ALLOCATE( TCAT ( NCOLS, NROWS, 4 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TCAT', PROGNAME )
        ALLOCATE( NLUS ( NCOLS, NROWS, 4 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NLUS', PROGNAME )
        ALLOCATE( IFLAG ( IB2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IFLAG', PROGNAME )
        ALLOCATE( TOTAL ( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TOTAL', PROGNAME )

C.......  Initialize arrays

        LUB2 =  0.000 ! array
        TCAT =  0.000 ! array
        NLUS =  0     ! array
        TOTAL = 0.000 ! array

C.......   Calculate category totals, landuse type totals for
C.......   BELD2 output file 

        DO I = 1, NCOLS
           DO J = 1, NROWS
              TTMP = 0
              IFLAG = 0   !array

C.......  If Forest inventory used greater than 0
             
              IF ( FIA( I, J ) .GT. 0.000 ) THEN

C.......   Convert from percent to hectares
                          
                 DO N = 1, IXREFS

                    BTMP = LUSE( I, J, N ) * 0.01 * XRES * YRES * TOHECT

C.......  Assume these 4 USGS classes are Scrub

                    IF ( N .EQ. IUSDEC .OR. N .EQ. IUSBRD .OR.
     &                   N .EQ. IUSCON .OR. N .EQ. IUSMIX )  THEN
   
                         TMPVEG = 'Scru'
                         K = 4
                    ELSE
                         TMPVEG = VEG2( N )
                         K = ICAT( N )
                    ENDIF

                    IF ( BTMP .GT. 0.00049 ) THEN
 
                       TTMP = TTMP + BTMP

                       DO M = 1, IB2
                          L = INDXB2 ( M )
                          IF ( TMPVEG .EQ. B2SORT( L ) ) THEN

                             LUB2( I, J, L ) = LUB2( I, J, L ) + BTMP
                             TCAT( I, J, K ) = TCAT( I, J, K ) + BTMP

                             IF ( IFLAG( L ) .EQ. 0 ) THEN
                                NLUS( I, J, K ) = NLUS( I, J, K ) + 1
                                IFLAG( L ) = 1
                             ENDIF

                             TOTAL( I, J ) = TOTAL( I, J ) + BTMP
                          ENDIF
                       ENDDO
                    ENDIF
                 ENDDO

C.........  If FIA less than or equal 0 (Canada) then use USGS xrefs as is

              ELSE
                          
                 DO N = 1, IXREFS
                    BTMP = LUSE( I, J, N ) * 0.01 * XRES * YRES * TOHECT

                    TTMP = TTMP + BTMP
                    
                    IF ( BTMP .GT. 0.00049 ) THEN
                       K = ICAT( N )
                       TMPVEG = VEG2( N )

                       DO M = 1, IB2
                          L = INDXB2 ( M ) 
                          IF ( TMPVEG .EQ. B2SORT( L ) ) THEN
                              
                             LUB2( I, J, L ) = LUB2( I, J, L ) + BTMP
                             TCAT( I, J, K ) = TCAT( I, J, K ) + BTMP

                             IF ( IFLAG( L ) .EQ. 0 ) THEN
                                NLUS( I, J, K ) = NLUS( I, J, K ) + 1
                                IFLAG( L ) = 1
                             ENDIF

                             TOTAL( I, J ) = TOTAL( I, J ) + BTMP
                          ENDIF
                       ENDDO
                    ENDIF
                 ENDDO
              ENDIF
              
           ENDDO
        ENDDO


C...............   Write BELD2 output file

        DO I = 1, NCOLS
           DO J = 1, NROWS
 
              WRITE( XDEV, 93030 ) I, J , TOTAL( I, J ) 
              DO K = 1, 4
                 WRITE( XDEV, 93040 ) TCAT( I ,J, K ), NLUS( I, J, K )
                 DO M = 1, IB2
                    L = INDXB2( M )
                    IF ( IB2CAT( L ) .EQ. K ) THEN
                       IF ( LUB2( I, J, L ) .GT. 0.00049 ) THEN
                     
                          WRITE( XDEV, 93050 ) AQUT, B2SORT( L ),
     &                                         AQUT, LUB2( I, J, L ) 

                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO

           ENDDO
        ENDDO

C.........   End of program:

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )

C...........   Formatted file I/O formats............ 93xxx
                                   
93000   FORMAT( A )
93010   FORMAT( A4, 1x, I3, 1x , A16 )
93020   FORMAT( 1x, A4 )
93030   FORMAT( 2(I5,1x), F13.3 )
93040   FORMAT( F13.3, 1x, I3 )
93050   FORMAT( 11x, A, A4, A, F13.3 )
93100   FORMAT('#GRID ',A, 2F9.0,2F7.0,3I4, 2(1x, A), 5F5.0 )

C...........   Internal buffering formats............ 94xxx
94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

        END PROGRAM  BELD3TO2 

