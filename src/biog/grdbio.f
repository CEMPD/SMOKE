       SUBROUTINE  GRDBIO ( GDEV, NCOLS, NROWS) 

C***********************************************************************
C  subroutine body starts at line 103
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
C*************************************************************************

C...........   MODULES for public variables
C...........   This module contains the biogenic variables

        USE MODBIOG

        IMPLICIT NONE

        INCLUDE 'EMCNST3.EXT'     ! emissions constants
        INCLUDE 'BIODIMS3.EXT'    ! biogenic-related constants

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         FIND1
        INTEGER         GETFLINE
        INTEGER         INDEX1
        INTEGER         STR2INT
        REAL            STR2REAL

        EXTERNAL        FIND1, GETFLINE, INDEX1, 
     &                  STR2INT, STR2REAL

        INTEGER, INTENT (IN)  ::  GDEV    !  unit number for county landuse file
        INTEGER, INTENT (IN)  ::  NCOLS   !  no. cols
        INTEGER, INTENT (IN)  ::  NROWS   !  no. rows

C.......   Source-level variables

        REAL    EMIS, AFOR, ENOFOR
        REAL    AAG, ENOAG, ENOGRS, ENOWTF

        REAL, ALLOCATABLE :: EOTH ( : )   ! other landuse
        REAL, ALLOCATABLE :: EPINE( : )   ! pine
        REAL, ALLOCATABLE :: EDECD( : )   ! deciduous
        REAL, ALLOCATABLE :: ECONF( : )   ! coniferous
        REAL, ALLOCATABLE :: EAG  ( : )   ! agriculture
        REAL, ALLOCATABLE :: ELAI ( : )   ! leaf area index

        REAL    AOTH, AGRS, AWTF, ALAI, SUMLAI, AVGLAI

        INTEGER         NLINE   !  Line number counter for BGUSE
        INTEGER         FLINES  !  Number of lines in BGUSE

        INTEGER         IOS     !  I/O status result

        INTEGER         B, C, R, I, J, K, L, M, N ! loop counters and subscripts
        INTEGER         CELL    ! cell no. in grid

        REAL            AREA    !  land-use area
        REAL            ATYPE   !  area for this land use type
        INTEGER         NTYPE   !  number of land use types
        CHARACTER*4     TYPE    !  land use type
        CHARACTER*80    INBUF   !  input buffer
        CHARACTER*256   MESG    !  message buffer for M3EXIT()

        LOGICAL         EFLAG   !  error flag

        CHARACTER*16 :: PROGNAME = 'GRDBIO'   !  program name

C***********************************************************************
C   begin body of subroutine  GRDBIO

C.......   Allocate memory for reading and storing county land use

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

C.......   Read and process gridded land use file:

        WRITE( MESG, 93000 ) 'Reading GRIDDED LAND USE file'
        CALL M3MSG2( MESG )

C.......  Get length of BGUSE file

        FLINES = GETFLINE( GDEV, ' Gridded land use file' )  

        EFLAG = .FALSE.

        NLINE = 0

C....... Read until end of the BGUSE file

        DO WHILE ( NLINE .LT. FLINES )  

            READ( GDEV, 93000, IOSTAT=IOS ) INBUF
            NLINE =  NLINE + 1
            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &          'Error reading GRIDDED LANDUSE file at line', NLINE
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            INBUF = ADJUSTL( INBUF )

C........  If header line found cycle

            IF ( INBUF(1:1) .EQ. CINVHDR ) CYCLE 

C........  Read column, row and total area in cell

            READ ( INBUF, * , IOSTAT=IOS ) C, R, AREA

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &          'Error reading GRIDDED LANDUSE file at line', NLINE
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            
C...........   If landuse file provides more grid cells than the test domain
C              skip them.

            IF ( (C .GT. NCOLS) .OR. (R .GT. NROWS) ) THEN
               WRITE( MESG,92020 )
     &               'Grid ( ', C, ', ', R, ' ) ' //
     &               'of landuse file is outside the domain'
               CALL M3WARN( PROGNAME, 0, 0, MESG )

               DO B = 1, 4    ! Skip all four blocks in landuse file
                  READ( GDEV, *, IOSTAT=IOS ) ATYPE, NTYPE
                  NLINE =  NLINE + 1
                  IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &             'Error reading GRIDDED ' //
     &             'LANDUSE at line', NLINE
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                  END IF

                  DO N = 1, NTYPE
                     READ( GDEV, 93000 , IOSTAT=IOS ) INBUF
                     NLINE = NLINE + 1
                     IF ( IOS .NE. 0 ) THEN
                       EFLAG = .TRUE.
                       WRITE( MESG,94010 )
     &                'Error reading GRIDDED ' //
     &                'LANDUSE at line', NLINE
                       CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                  END DO
               END DO

            ELSE

C...........   Land use type:  (rural) forest.  Process subtypes:

               READ( GDEV, *, IOSTAT=IOS ) ATYPE, NTYPE
               NLINE =  NLINE + 1
               IF ( IOS .NE. 0 ) THEN
                 EFLAG = .TRUE.
                 WRITE( MESG,94010 )
     &           'Error reading FOREST AREA from ' //
     &           'LANDUSE at line', NLINE
                 CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
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
     &              'Error reading AREA TYPE from ' //
     &              'GRIDDED LAND USE at line', NLINE
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF

                 INBUF = ADJUSTL( INBUF )  
                 TYPE  = INBUF( 2 : 5 )
                 AREA  = STR2REAL( INBUF( 7 : 80 ) )
                 K     = INDEX1( TYPE, NVEG, VEGID )

                 IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'Could not find "' // TYPE //
     &              '" from LU FILE in VEGID at line', NLINE
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF

                 L = INDEX1( TYPE, SPTREE, SPFORID )

                 IF ( L .GT. 0 ) THEN
                    AFOR   = AFOR   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    ALAI   = ALAI   + AREA
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Wdcp' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AAG    = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOAG  = ENOAG + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AGRS   = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2*AREA
                    AGRS   = AAG    + 0.2*AREA
                    AOTH   = AAG    + 0.6*AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE    !  Add area to the forest area, add NO to the forest NO


                    AFOR = AFOR + AREA
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                    IF     ( LAI( K ) .EQ. 3 ) THEN
                        DO  M = 1, BSPCS - 1
                            EPINE( M ) = EPINE( M ) + AREA*EMFAC( K,M )
                        ENDDO

                    ELSE IF( LAI( K ) .EQ. 5 ) THEN

                        DO  M = 1, BSPCS - 1
                            EDECD( M ) = EDECD( M ) + AREA*EMFAC( K,M )
                        ENDDO

                    ELSE IF( LAI( K ) .EQ. 7 ) THEN

                        DO  M = 1, BSPCS - 1
                            ECONF( M ) = ECONF( M ) + AREA*EMFAC( K,M )
                        ENDDO

                    ELSE

                        SUMLAI = SUMLAI + AREA * LAI( K )
                        ALAI   = ALAI   + AREA
                        DO  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA*EMFAC( K,M )
                        ENDDO

                    END IF      !  if lai is 3,5,7, or otherwise

                  END IF  !  if some spforid, or 'Wdcp' or 'Scwd' or 'Urba or not
               ENDDO

C...........   Land use type:  urban forest.  Process subtypes:

               READ( GDEV, 93000, IOSTAT=IOS ) INBUF
               NLINE =  NLINE + 1
               IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading URBAN FOREST from ' //
     &          'GRIDDED LAND USE at line', NLINE
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
               END IF

               INBUF = ADJUSTL ( INBUF ) 
               READ( INBUF, *, IOSTAT=IOS ) ATYPE, NTYPE

               IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &          'Error reading URBAN FOREST from ' //
     &          'GRIDDED LAND USE at line', NLINE
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
               END IF


               DO  N = 1, NTYPE

                 READ( GDEV, 93000, IOSTAT=IOS ) INBUF
                 NLINE =  NLINE + 1
                 IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'GRIDDED LAND USE at line', NLINE
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF

                 INBUF = ADJUSTL( INBUF )
                 TYPE  = INBUF( 2 : 5 )
                 AREA  = STR2REAL( INBUF( 7 : 80 ) )

                 K     = INDEX1( TYPE, NVEG, VEGID )

                 IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LU FILE in VEGID at line', NLINE
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF

                 L = INDEX1( TYPE, SPTREE, SPFORID )

                 IF ( L .GT. 0 ) THEN

                    AFOR   = AFOR   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    ALAI   = ALAI   + AREA
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Wdcp' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AAG    = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOAG  = ENOAG + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Scwd' ) THEN

                    AFOR   = AFOR   + 0.5*AREA
                    AGRS   = AAG    + 0.5*AREA
                    ALAI   = ALAI   + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2*AREA
                    AGRS   = AAG    + 0.2*AREA
                    AOTH   = AAG    + 0.6*AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE    !  Add area to the forest area, add NO to the forest NO

                    AFOR = AFOR + AREA
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                    IF     ( LAI( K ) .EQ. 3 ) THEN
                        DO  M = 1, BSPCS - 1
                            EPINE( M ) = EPINE( M ) + AREA*EMFAC( K,M )
                        ENDDO

                    ELSE IF( LAI( K ) .EQ. 5 ) THEN

                        DO  M = 1, BSPCS - 1
                            EDECD( M ) = EDECD( M ) + AREA*EMFAC( K,M )
                        ENDDO

                    ELSE IF( LAI( K ) .EQ. 7 ) THEN

                        DO  M = 1, BSPCS - 1
                            ECONF( M ) = ECONF( M ) + AREA*EMFAC( K,M )
                        ENDDO

                    ELSE

                        SUMLAI = SUMLAI + AREA * LAI( K )
                        ALAI   = ALAI   + AREA
                        DO  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA*EMFAC( K,M )
                        ENDDO

                    END IF      !!  if lai is 3,5,7, or otherwise

                 END IF  !  if some spforid, or 'Wdcp' or 'Scwd' or 'Urba or not

               ENDDO


C...........   Land use type:  agriculture.  Process subtypes:

               READ( GDEV, 93000, IOSTAT=IOS ) INBUF
               NLINE =  NLINE + 1
               IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading AGRICULTURE from ' //
     &          'COUNTY LAND USE at line', NLINE
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
               END IF

               INBUF = ADJUSTL ( INBUF ) 
               READ( INBUF, *, IOSTAT=IOS ) ATYPE, NTYPE

               IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &          'Error reading URBAN FOREST from ' //
     &          'GRIDDED LAND USE at line', NLINE
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
               END IF

               DO  N = 1, NTYPE

                 READ( GDEV, 93000, IOSTAT=IOS ) INBUF
                 NLINE =  NLINE + 1
                 IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'COUNTY LAND USE at line', NLINE
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF

                 INBUF = ADJUSTL( INBUF )
                 TYPE  = INBUF( 2 : 5 )
                 AREA  = STR2REAL( INBUF( 7 : 80 ) )

                 K     = INDEX1( TYPE, NVEG, VEGID )

                 IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LU FILE in VEGID at line', NLINE
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF

                 AAG = AAG + AREA
                 DO  M = 1, BSPCS - 1
                    EAG( M ) = EAG( M ) + AREA * EMFAC( K , M )
                 ENDDO
                 ENOAG  = ENOAG + AREA * EMFAC( K, NO )

               ENDDO


C...........   Land use type:  remaining/other.  Process subtypes:

               READ( GDEV, 93000, IOSTAT=IOS ) INBUF
               NLINE =  NLINE + 1
               IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &          'Error reading OTHER AREA from ' //
     &          'COUNTY LAND USE at line', NLINE
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
               END IF

               INBUF = ADJUSTL ( INBUF ) 
               READ( INBUF, *, IOSTAT=IOS ) ATYPE, NTYPE

               IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &          'Error reading URBAN FOREST from ' //
     &          'GRIDDED LAND USE at line', NLINE
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
               END IF

               DO  N = 1, NTYPE

                 READ( GDEV, 93000, IOSTAT=IOS ) INBUF
                 NLINE =  NLINE + 1
                 IF ( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 
     &              'Error reading AREA TYPE from ' //
     &              'COUNTY LAND USE at line', NLINE
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF

                 INBUF = ADJUSTL( INBUF )
                 TYPE  = INBUF( 2 : 5 )
                 AREA  = STR2REAL( INBUF( 7 : 80 ) )

                 K     = INDEX1( TYPE, NVEG, VEGID )

                 IF ( K .EQ. 0 ) THEN
                    WRITE( MESG,94010 ) 
     &              'Could not find "' // TYPE // 
     &              '" from LU FILE in VEGID at line', NLINE 
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                 END IF

                 L = INDEX1( TYPE, RMTREE, OTHERID )

                 IF ( L .GT. 0 ) THEN

                    ALAI   = ALAI + AREA
                    AFOR   = AFOR + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Pacp' ) THEN

                    AAG = AAG + AREA
                    DO  M = 1, BSPCS - 1
                        EAG( M ) = EAG( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOAG = ENOAG + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Gras' .OR.

     &                    TYPE .EQ. 'Scru' .OR.
     &                    TYPE .EQ. 'Ugra' .OR.
     &                    TYPE .EQ. 'Othe' ) THEN
                    AGRS = AGRS + AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Wetf' ) THEN

                    AWTF = AWTF + AREA
                    ALAI = ALAI + AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        EMIS = AREA * EMFAC( K , M )
                        ELAI( M ) = ELAI( M ) + EMIS
                    ENDDO
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
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOAG = ENOAG + AREA * EMFAC( K, NO )
                 ELSE IF ( TYPE .EQ. 'Scwd' ) THEN
                    ALAI   = ALAI   + AREA
                    AFOR   = AFOR   + 0.5 * AREA
                    AGRS   = AGRS   + 0.5 * AREA
                    SUMLAI = SUMLAI + AREA * LAI( K )
                    DO  M = 1, BSPCS - 1
                        ELAI( M ) = ELAI( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Urba' ) THEN

                    AFOR   = AFOR   + 0.2 * AREA
                    AGRS   = AGRS   + 0.2 * AREA
                    AOTH   = AOTH   + 0.6 * AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE IF ( TYPE .EQ. 'Desh' ) THEN

                    AGRS   = AGRS   + 0.5 * AREA
                    AOTH   = AOTH   + 0.5 * AREA
                    DO  M = 1, BSPCS - 1
                        EOTH( M ) = EOTH( M ) + AREA * EMFAC( K , M )
                    ENDDO
                    ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                 ELSE

                    AOTH = AOTH + AREA

                    IF ( LAI( K ) .GT. 0 ) THEN

                        ALAI   = ALAI   + AREA
                        SUMLAI = SUMLAI + AREA * LAI( K )
                        DO  M = 1, BSPCS - 1
                            ELAI( M ) = ELAI( M ) + AREA * EMFAC( K,M )
                        ENDDO
                        ENOFOR = ENOFOR + AREA * EMFAC( K, NO )

                    ELSE

                        DO  M = 1, BSPCS - 1
                            EOTH( M ) = EOTH( M ) + AREA * EMFAC( K,M )
                        ENDDO
                        ENOGRS = ENOGRS + AREA * EMFAC( K, NO )

                    END IF

                 END IF

               ENDDO

               IF ( ALAI .GT. 0.0 ) THEN
                   AVGLAI = SUMLAI / ALAI 
               ELSE
                   AVGLAI = 0.0
               END IF

C..............  Assign normalized emissions to appropriate
C..............  grid cell

               CELL = ( ( R - 1 ) * NCOLS )  + C
               DO M = 1, BSPCS - 1
                  PINE( CELL, M ) = EPINE( M )
                  DECD( CELL, M ) = EDECD( M ) 
                  CONF( CELL, M ) = ECONF( M )
                  AGRC( CELL, M ) = EAG ( M )
                  LEAF( CELL, M ) = ELAI( M ) 
                  OTHR( CELL, M ) = EOTH( M )
               ENDDO

               GRASNO( CELL ) = ENOGRS
               FORENO( CELL ) = ENOFOR
               WETLNO( CELL ) = ENOWTF
               AGRINO( CELL ) = ENOAG
               AVLAI ( CELL ) = AVGLAI

          ENDIF

        ENDDO

C******************  FORMAT  STATEMENTS   ******************************

C...........   Error and warning message formats..... 91xxx
C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT ( 5X , A )

92020   FORMAT ( 5X , A, I4, A, I4, A )

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )


  
        RETURN

        END
 
