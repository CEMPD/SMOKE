
        SUBROUTINE GRDBIO( GDEV, LUNROWS, LUNCOLS ) 

C***********************************************************************
C  subroutine body starts at line 111
C
C  DESCRIPTION:
C       Computes normalized gridded biogenic emissions in terms of
C       biomass, gridded land use, and emissions factors.d
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       02/00 JMV : prototype 
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
C*************************************************************************

C...........   MODULES for public variables
C...........   This module contains the biogenic variables
        USE MODBIOG, ONLY: NVEG, VEGID, PINE, DECD, CONF, AGRC, LEAF,
     &                     OTHR, GRASNO, FORENO, WETLNO, AGRINO, AVLAI,
     &                     LAI, EMFAC   

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: NCOLS, NROWS, XOFF, YOFF

        IMPLICIT NONE

        INCLUDE 'EMCNST3.EXT'     ! emissions constants
        INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         GETFLINE
        INTEGER         INDEX1
        REAL            STR2REAL

        EXTERNAL        GETFLINE, INDEX1, STR2REAL

C...........   SUBROUTINE ARGUMENTS

        INTEGER, INTENT (IN)  ::  GDEV    ! unit no. for county landuse file
        INTEGER, INTENT (IN)  ::  LUNROWS ! no. rows in land use header
        INTEGER, INTENT (IN)  ::  LUNCOLS ! no. columns in land use header

C...........   Allocatable, source-level variables

        REAL, ALLOCATABLE :: EOTH ( : )   ! other landuse
        REAL, ALLOCATABLE :: EPINE( : )   ! pine
        REAL, ALLOCATABLE :: EDECD( : )   ! deciduous
        REAL, ALLOCATABLE :: ECONF( : )   ! coniferous
        REAL, ALLOCATABLE :: EAG  ( : )   ! agriculture
        REAL, ALLOCATABLE :: ELAI ( : )   ! leaf area index

C...........   Other local variables

        INTEGER         B, I, J, K, L, M, N ! loop counters and subscripts

        INTEGER         CELL      !  tmp cell no. in grid
        INTEGER         COL       !  tmp column no.
        INTEGER         FLINES    !  number of lines in BGUSE
        INTEGER         IOS       !  i/o status
        INTEGER         LC, LR    !  length of COLRANGE & ROWRANGE
        INTEGER         NLINE     !  line number counter for BGUSE
        INTEGER         NTYPE     !  number of land use types
        INTEGER         ROW       !  tmp row no.

        REAL            AOTH, AGRS, AWTF, ALAI, SUMLAI, AVGLAI
        REAL            EMIS, AFOR, ENOFOR
        REAL            AAG, ENOAG, ENOGRS, ENOWTF
        REAL            AREA      !  land-use area
        REAL            ATYPE     !  area for this land use type

        CHARACTER(4)    LUTYPE    !  land use type
        CHARACTER(20)   COLRANGE  !  buffer w/ column range
        CHARACTER(20)   ROWRANGE  !  buffer w/ row range
        CHARACTER(80)   INBUF     !  input buffer
        CHARACTER(300)  MESG      !  message buffer 

        LOGICAL      :: EFLAG    = .FALSE.    ! true: error found
        LOGICAL      :: INHEADER = .TRUE.     ! true: in BEIS header section
        LOGICAL      :: OFLAG    = .FALSE.    ! true: at least 1 warning found
        LOGICAL      :: WFLAG    = .FALSE.    ! true: warning on current record

        CHARACTER(16) :: PROGNAME = 'GRDBIO'   !  program name

C***********************************************************************
C   begin body of subroutine  GRDBIO

C.........  Allocate memory for reading and storing county land use
        ALLOCATE( EPINE ( BSPCS-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EPINE', PROGNAME )
        ALLOCATE( EDECD ( BSPCS-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EDECD', PROGNAME )
        ALLOCATE( ECONF ( BSPCS-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ECONF', PROGNAME )
        ALLOCATE( EOTH ( BSPCS-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EOTH', PROGNAME )
        ALLOCATE( ELAI ( BSPCS-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ELAI', PROGNAME )
        ALLOCATE( EAG ( BSPCS-1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EAG', PROGNAME )

C.........  Create message fields for errors
        WRITE( COLRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, LUNCOLS
        WRITE( ROWRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, LUNROWS

        LC = LEN_TRIM( COLRANGE )
        LR = LEN_TRIM( ROWRANGE )

C.........  Read and process gridded land use file:

        MESG = 'Reading GRIDDED LAND USE file'
        CALL M3MSG2( MESG )

C.........  Get length of BGUSE file

        FLINES = GETFLINE( GDEV, ' Gridded land use file' )  

C.........  Read until end of the BGUSE file
        NLINE = 0
        DO WHILE ( NLINE .LT. FLINES )  

            READ( GDEV, 93000, IOSTAT=IOS ) INBUF
            NLINE =  NLINE + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS,
     &                 'reading gridded land use file at line', NLINE
                CALL M3MESG( MESG )
                CYCLE
            END IF

            INBUF = ADJUSTL( INBUF )

C.............  Skip all lines in file until ENDHEADER found
            IF( INHEADER ) THEN
                J = INDEX( INBUF, 'ENDHEADER' )
                IF( J .GT. 0 ) INHEADER = .FALSE.
            END IF

C.............  If SMOKE header line found cycle
            IF ( INBUF(1:1) .EQ. CINVHDR ) THEN
                INHEADER = .FALSE.
                CYCLE 
            END IF

C.............  Skip blank lines
            IF( INBUF .EQ. ' ' ) CYCLE

C.............  Initialize warning flag
            WFLAG = .FALSE.

C.............  Read column, row and total area in cell

            READ ( INBUF, * , IOSTAT=IOS ) COL, ROW, AREA

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS,
     &                 'reading gridded land use file at line', NLINE
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.................  Check the value of the column number
            IF( COL .LE. 0 .OR. COL .GT. LUNCOLS ) THEN
                WFLAG = .TRUE.
                OFLAG = .TRUE.
                WRITE( MESG,94010 ) 'WARNING: Column value ', COL,
     &                     'is outside range ' // COLRANGE( 1:LC ) // 
     &                     ' at line', NLINE
                CALL M3MESG( MESG )
            END IF

C.................  Check the value of the row number
            IF( ROW .LE. 0 .OR. ROW .GT. LUNROWS ) THEN
                WFLAG = .TRUE.
                OFLAG = .TRUE.
                WRITE( MESG,94010 ) 'WARNING: Row value ', ROW,
     &                 'is outside range ' // ROWRANGE( 1:LR ) // 
     &                 ' at line', NLINE
                CALL M3MESG( MESG )                    
            END IF

C.............  Adjust column and row for subgrid
            COL = COL - XOFF
            ROW = ROW - YOFF

C.............  Skip entry after subgrid adjustment
            IF( COL .LE. 0 .OR. COL .GT. NCOLS .OR.
     &          ROW .LE. 0 .OR. ROW .GT. NROWS       ) WFLAG = .TRUE.

C.............   If landuse file has grid cells outside the valid ranges
C                given in the header, skip them
            IF( WFLAG ) THEN

               DO B = 1, 4    ! Skip all four blocks in landuse file
                  READ( GDEV, *, IOSTAT=IOS ) ATYPE, NTYPE
                  NLINE =  NLINE + 1

                  IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'I/O error', IOS, 'while '//
     &                'skipping line', NLINE, 'of gridded land use'
                    CALL M3MESG( MESG )
                    CYCLE
                  END IF

                  DO N = 1, NTYPE
                     READ( GDEV, 93000 , IOSTAT=IOS ) INBUF
                     NLINE = NLINE + 1
                     IF ( IOS .NE. 0 ) THEN
                       EFLAG = .TRUE.
                       WRITE( MESG,94010 ) 'I/O error', IOS, 'while '//
     &                   'skipping line', NLINE, 'of gridded land use'
                       CALL M3MESG( MESG )
                       CYCLE
                    END IF

                  END DO
               END DO

C.............  If no warning, continue with read...
            ELSE

C................  Land use type: (rural) forest.  Process subtypes:
               READ( GDEV, *, IOSTAT=IOS ) ATYPE, NTYPE
               NLINE =  NLINE + 1
               IF ( IOS .NE. 0 ) THEN
                 EFLAG = .TRUE.
                 WRITE( MESG,94010 )
     &             'I/O error', IOS, 'reading FOREST AREA from ' //
     &             'gridded land use at line', NLINE
                 CALL M3MESG( MESG )
                 CYCLE
               END IF

               EPINE = 0.0    ! array
               EDECD = 0.0    ! array
               ECONF = 0.0    ! array
               EOTH  = 0.0    ! array
               ELAI  = 0.0    ! array
               EAG   = 0.0    ! array

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

               DO  N = 1, NTYPE

                 READ( GDEV, 93000, IOSTAT=IOS ) INBUF
                 NLINE =  NLINE + 1
                 IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                'I/O error', IOS, 'reading AREA TYPE from ' //
     &                'gridded land use at line', NLINE
                    CALL M3MESG( MESG )
                    CYCLE
                 END IF

                 INBUF = ADJUSTL( INBUF )  
                 LUTYPE  = INBUF( 2 : 5 )
                 AREA  = STR2REAL( INBUF( 7 : 80 ) )
                 K     = INDEX1( LUTYPE, NVEG, VEGID )

                 IF ( K .EQ. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                'Could not find "' // LUTYPE //
     &                '" from LU FILE in VEGID at line', NLINE
                    CALL M3MESG( MESG )
                    CYCLE
                 END IF

                 L = INDEX1( LUTYPE, SPTREE, SPFORID )

                 IF ( L .GT. 0 ) THEN
                    AFOR   = AFOR   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    ALAI   = ALAI   + AREA
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Wdcp' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AAG    = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOAG  = ENOAG + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AGRS   = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2*AREA
                    AGRS   = AAG    + 0.2*AREA
                    AOTH   = AAG    + 0.6*AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE    !  Add area to the forest area, add NO to the forest NO


                    AFOR = AFOR + AREA
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                    IF     ( LAI( K ) .EQ. 3 ) THEN
                        DO  M = 1, BSPCS - 1
                            EPINE( M ) = EPINE( M ) + AREA*EMFAC( K,M )
                        END DO

                    ELSE IF( LAI( K ) .EQ. 5 ) THEN

                        DO  M = 1, BSPCS - 1
                            EDECD( M ) = EDECD( M ) + AREA*EMFAC( K,M )
                        END DO

                    ELSE IF( LAI( K ) .EQ. 7 ) THEN

                        DO  M = 1, BSPCS - 1
                            ECONF( M ) = ECONF( M ) + AREA*EMFAC( K,M )
                        END DO

                    ELSE

                        SUMLAI = SUMLAI + AREA * LAI( K )
                        ALAI   = ALAI   + AREA
                        DO  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA*EMFAC( K,M )
                        END DO

                    END IF      !  if lai is 3,5,7, or otherwise

                  END IF  !  if some spforid, or 'Wdcp' or 'Scwd' or 'Urba or not
               END DO

C...........   Land use type:  urban forest.  Process subtypes:

               READ( GDEV, 93000, IOSTAT=IOS ) INBUF
               NLINE =  NLINE + 1
               IF ( IOS .NE. 0 ) THEN
                 EFLAG = .TRUE.
                 WRITE( MESG,94010 )
     &                 'I/O error', IOS, 'reading URBAN FOREST from ' //
     &                 'gridded land use at line', NLINE
                 CALL M3MESG( MESG )
                 CYCLE
               END IF

               INBUF = ADJUSTL ( INBUF ) 
               READ( INBUF, *, IOSTAT=IOS ) ATYPE, NTYPE

               IF ( IOS .NE. 0 ) THEN
                 EFLAG = .TRUE.
                 WRITE( MESG,94010 )
     &                 'I/O error', IOS, 'reading URBAN FOREST from ' //
     &                 'gridded land use at line', NLINE
                 CALL M3MESG( MESG )
                 CYCLE
              END IF


               DO  N = 1, NTYPE

                 READ( GDEV, 93000, IOSTAT=IOS ) INBUF
                 NLINE =  NLINE + 1
                 IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                'I/O error', IOS, 'reading AREA TYPE from ' //
     &                'gridded land use at line', NLINE
                    CALL M3MESG( MESG )
                    CYCLE
                 END IF

                 INBUF = ADJUSTL( INBUF )
                 LUTYPE  = INBUF( 2 : 5 )
                 AREA  = STR2REAL( INBUF( 7 : 80 ) )

                 K     = INDEX1( LUTYPE, NVEG, VEGID )

                 IF ( K .EQ. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // LUTYPE // 
     &              '" from LU FILE in VEGID at line', NLINE
                    CALL M3MESG( MESG )
                    CYCLE
                 END IF

                 L = INDEX1( LUTYPE, SPTREE, SPFORID )

                 IF ( L .GT. 0 ) THEN

                    AFOR   = AFOR   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    ALAI   = ALAI   + AREA
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Wdcp' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AAG    = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOAG  = ENOAG + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AGRS   = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2*AREA
                    AGRS   = AAG    + 0.2*AREA
                    AOTH   = AAG    + 0.6*AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE    !  Add area to the forest area, add NO to the forest NO

                    AFOR = AFOR + AREA
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                    IF     ( LAI( K ) .EQ. 3 ) THEN
                        DO  M = 1, BSPCS - 1
                            EPINE( M ) = EPINE( M ) + AREA*EMFAC( K,M )
                        END DO

                    ELSE IF( LAI( K ) .EQ. 5 ) THEN

                        DO  M = 1, BSPCS - 1
                            EDECD( M ) = EDECD( M ) + AREA*EMFAC( K,M )
                        END DO

                    ELSE IF( LAI( K ) .EQ. 7 ) THEN

                        DO  M = 1, BSPCS - 1
                            ECONF( M ) = ECONF( M ) + AREA*EMFAC( K,M )
                        END DO

                    ELSE

                        SUMLAI = SUMLAI + AREA * LAI( K )
                        ALAI   = ALAI   + AREA
                        DO  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA*EMFAC( K,M )
                        END DO

                    END IF      !!  if lai is 3,5,7, or otherwise

                 END IF  !  if some spforid, or 'Wdcp' or 'Scwd' or 'Urba or not

               END DO


C...........   Land use type:  agriculture.  Process subtypes:

               READ( GDEV, 93000, IOSTAT=IOS ) INBUF
               NLINE =  NLINE + 1
               IF ( IOS .NE. 0 ) THEN
                 EFLAG = .TRUE.
                 WRITE( MESG,94010 )
     &                'I/O error', IOS, 'reading AGRICULTURE from ' //
     &                'gridded land use at line', NLINE
                 CALL M3MESG( MESG )
                 CYCLE
               END IF

               INBUF = ADJUSTL ( INBUF ) 
               READ( INBUF, *, IOSTAT=IOS ) ATYPE, NTYPE

               IF ( IOS .NE. 0 ) THEN
                 EFLAG = .TRUE.
                 WRITE( MESG,94010 )
     &                'I/O error', IOS, 'reading AGRICULTURE from ' //
     &                'gridded land use at line', NLINE
                 CALL M3MESG( MESG )
                 CYCLE
               END IF

               DO  N = 1, NTYPE

                 READ( GDEV, 93000, IOSTAT=IOS ) INBUF
                 NLINE =  NLINE + 1
                 IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                'I/O error', IOS, 'reading AREA TYPE from ' //
     &                'gridded land use at line', NLINE
                    CALL M3MESG( MESG )
                    CYCLE
                 END IF

                 INBUF = ADJUSTL( INBUF )
                 LUTYPE  = INBUF( 2 : 5 )
                 AREA  = STR2REAL( INBUF( 7 : 80 ) )

                 K     = INDEX1( LUTYPE, NVEG, VEGID )

                 IF ( K .EQ. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // LUTYPE // 
     &              '" from LU FILE in VEGID at line', NLINE
                    CALL M3MESG( MESG )
                 END IF

                 AAG = AAG + AREA
                 DO  M = 1, BSPCS - 1
                    EAG( M ) = EAG( M ) + AREA * EMFAC( K , M )
                 END DO
                 ENOAG  = ENOAG + AREA * EMFAC( K, NO )

               END DO


C...........   Land use type:  remaining/other.  Process subtypes:

               READ( GDEV, 93000, IOSTAT=IOS ) INBUF
               NLINE =  NLINE + 1
               IF ( IOS .NE. 0 ) THEN
                 EFLAG = .TRUE.
                 WRITE( MESG,94010 ) 
     &                  'I/O error', IOS, 'reading OTHER AREA from ' //
     &                 'gridded land use at line', NLINE
                 CALL M3MESG( MESG )
                 CYCLE
               END IF

               INBUF = ADJUSTL ( INBUF ) 
               READ( INBUF, *, IOSTAT=IOS ) ATYPE, NTYPE

               IF ( IOS .NE. 0 ) THEN
                 EFLAG = .TRUE.
                 WRITE( MESG,94010 ) 
     &                  'I/O error', IOS, 'reading OTHER AREA from ' //
     &                  'gridded land use at line', NLINE
                 CALL M3MESG( MESG )
                 CYCLE
               END IF

               DO  N = 1, NTYPE

                 READ( GDEV, 93000, IOSTAT=IOS ) INBUF
                 NLINE =  NLINE + 1
                 IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                'I/O error', IOS, 'reading AREA TYPE from ' //
     &                'gridded land use at line', NLINE
                    CALL M3MESG( MESG )
                    CYCLE
                 END IF

                 INBUF = ADJUSTL( INBUF )
                 LUTYPE  = INBUF( 2 : 5 )
                 AREA  = STR2REAL( INBUF( 7 : 80 ) )

                 K     = INDEX1( LUTYPE, NVEG, VEGID )

                 IF ( K .EQ. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &                  'Could not find "' // LUTYPE // 
     &                  '" from LU FILE in VEGID at line', NLINE 
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                 END IF

                 L = INDEX1( LUTYPE, RMTREE, OTHERID )

                 IF ( L .GT. 0 ) THEN

                    ALAI   = ALAI + AREA
                    AFOR   = AFOR + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Pacp' ) THEN

                    AAG = AAG + AREA
                    DO  M = 1, BSPCS - 1
                        EAG( M ) = EAG( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOAG = ENOAG + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Gras' .OR.

     &                    LUTYPE .EQ. 'Scru' .OR.
     &                    LUTYPE .EQ. 'Ugra' .OR.
     &                    LUTYPE .EQ. 'Othe' ) THEN
                    AGRS = AGRS + AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Wetf' ) THEN

                    AWTF = AWTF + AREA
                    ALAI = ALAI + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        EMIS = AREA * EMFAC( K , M )
                        ELAI( M ) = ELAI( M ) + EMIS
                    END DO
                    ENOWTF = ENOWTF + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Wate' .OR.
     &                    LUTYPE .EQ. 'Barr' .OR.
     &                    LUTYPE .EQ. 'Uoth' ) THEN

                        AOTH = AOTH + AREA

                 ELSE IF ( LUTYPE .EQ. 'Wdcp' ) THEN

                    ALAI   = ALAI   + AREA
                    AFOR   = AFOR   + 0.5 * AREA
                    AAG    = AAG    + 0.5 * AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOAG = ENOAG + AREA * EMFAC( K, NO )
                 ELSE IF ( LUTYPE .EQ. 'Scwd' ) THEN
                    ALAI   = ALAI   + AREA
                    AFOR   = AFOR   + 0.5 * AREA
                    AGRS   = AGRS   + 0.5 * AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2 * AREA
                    AGRS   = AGRS   + 0.2 * AREA
                    AOTH   = AOTH   + 0.6 * AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( LUTYPE .EQ. 'Desh' ) THEN

                    AGRS   = AGRS   + 0.5 * AREA
                    AOTH   = AOTH   + 0.5 * AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    END DO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE

                    AOTH = AOTH + AREA

                    IF ( LAI( K ) .GT. 0 ) THEN

                        ALAI   = ALAI   + AREA
                        SUMLAI = SUMLAI + AREA * LAI( K )
                        DO  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA * EMFAC( K,M )
                        END DO
                        ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                    ELSE

                        DO  M = 1, BSPCS - 1
                            EOTH( M ) = EOTH( M ) + AREA * EMFAC( K,M )
                        END DO
                        ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                    END IF

                 END IF

               END DO

               IF ( ALAI .GT. 0.0 ) THEN
                   AVGLAI = SUMLAI / ALAI 
               ELSE
                   AVGLAI = 0.0
               END IF

C..............  Assign normalized emissions to appropriate
C..............  grid cell

               CELL = ( ( ROW - 1 ) * NCOLS )  + COL
               DO M = 1, BSPCS - 1
                  PINE( CELL, M ) = EPINE( M )
                  DECD( CELL, M ) = EDECD( M ) 
                  CONF( CELL, M ) = ECONF( M )
                  AGRC( CELL, M ) = EAG ( M )
                  LEAF( CELL, M ) = ELAI( M ) 
                  OTHR( CELL, M ) = EOTH( M )
               END DO

               GRASNO( CELL ) = ENOGRS
               FORENO( CELL ) = ENOFOR
               WETLNO( CELL ) = ENOWTF
               AGRINO( CELL ) = ENOAG
               AVLAI ( CELL ) = AVGLAI

          END IF

        END DO    ! End read loop for gridded land use

C.........  Error found: abort run
        IF( EFLAG ) THEN

            MESG = 'Problem reading gridding land use file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx
C...........   Informational (LOG) message formats... 92xxx

92020   FORMAT ( 5X , A, I4, A, I4, A )

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I8, :, 1X ) )

        END SUBROUTINE GRDBIO
 
