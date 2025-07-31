
        PROGRAM NORMBEIS4

C***********************************************************************
C
C  DESCRIPTION:  Produces normalized biogenic emissions for use with
C                SMOKE-BEIS4.   Inputs include emissions factor file,
C                gridded biomass file (new), and gridded land use 
C                file (Biogenic Emissions Landcover Database).
C
C  SUBROUTINES AND FUNCTIONS CALLED: RDB3FAC
C
C  REVISION  HISTORY: 3/2022 Initial version   J Vukovich and J Bash  
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2022, Environmental Modeling for Policy Development
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

C.........  MODULES for public variables
C.........  This module contains biogenic variables
C.......... MODBEIS3 can be used for BEIS4; BEIS4 uses less of the
C.......... module since LFBIO and SLW no longer needed.
        USE MODBEIS3, ONLY: NVEG, VEGID, AVGEMIS, AVGLAI, NOEMIS, 
     &                      EMFAC, LAI, WFAC
        USE MODGRDLIB
 
        IMPLICIT NONE

C.........  INCLUDES
        INCLUDE 'PARMS3.EXT'      ! I/O API constants
        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations
C.......... Same include file for BEIS3 can be use for BEIS4
        INCLUDE 'B3V14DIMS3.EXT'     ! BEIS3-related and BEIS4 declarations
  
C.........  EXTERNAL FUNCTIONS and their descriptions
        LOGICAL         ENVYN
        INTEGER         GETFLINE
        INTEGER         PROMPTFFILE
        CHARACTER(16)   PROMPTMFILE

        EXTERNAL        ENVYN, GETFLINE, PROMPTFFILE, PROMPTMFILE

C.........  LOCAL VARIABLES and their descriptions
        INTEGER         B, C, R, I, J, K, L, M, N ! loop counters and subscripts
        INTEGER         IOS     !  I/O status result

        INTEGER         FDEV    !  unit number for emissions factor file
        INTEGER         LDEV    !  unit number for log file

        CHARACTER(16)   ENAME   !  logical name for normalized emissions output
        CHARACTER(16)   BNAME   !  logical name for gridded biomass file 
        CHARACTER(16)   GRDNM   !  grid name
        CHARACTER(16), ALLOCATABLE  :: VGLIST( : )  ! land use type names	
        CHARACTER(16), ALLOCATABLE  :: BIOM_ID( : )  ! biomass veg names 
        CHARACTER(16), ALLOCATABLE  :: BUNITS( : )  ! biomass veg names        
        LOGICAL, ALLOCATABLE  :: VAG_YN ( : )   ! ag veg type or not

        CHARACTER(16)   LUNAME  !  logical name for gridded land use totals file
        CHARACTER(16)   BEISVER !  version of BEIS4 to use

        CHARACTER(256)  MESG    !  message buffer for M3EXIT()
        CHARACTER(5)    BTMP    ! temporary tag for naming output variables
C.........  CVS release tag
        CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv5.2_Jul2025$'

        INTEGER         NCOLS   ! no. of grid columns
        INTEGER         NROWS   ! no. of grid rows

        INTEGER         NBELD   ! total number variables in BELD file 
        INTEGER         NBIOM   ! total number variables in BIOMASS file
        INTEGER         IFOUND  ! used for checking land use vs. emis facs
        INTEGER         IGRASS  ! Grass

        INTEGER         MODIS14  ! modis mosaic 1/3 (grass + mxforest + drycrop)
        INTEGER     ::    LAI_SAVE_INDEX(3)
        INTEGER :: NLINES

        REAL, ALLOCATABLE    :: LUSE ( :, :, :  )  ! BELD land use data
        REAL, ALLOCATABLE    :: BIOM ( :, :, :  )  ! Gridded biomass data

        REAL, ALLOCATABLE    ::   SUMEM( : )       ! Summer emissions 
        REAL, ALLOCATABLE    ::   SUMEMW( : )      ! Winter emissions
        REAL, ALLOCATABLE    ::   NOEM( : )        ! NO emissions 
        REAL, ALLOCATABLE    ::   SUMLAI( : )      ! Summer LAIs
        REAL, ALLOCATABLE    ::   SUMLAIW( : )     ! Winter LAIs
        REAL  VEGAREA                              ! Veg. area tmp
        REAL  VEGBIOM                              ! Veg biomass tmp 
        REAL  EFTMP                                ! Emission factors
        REAL  SVEGA 
        DOUBLE PRECISION PRCNT2KM2     ! Prcnt to km**2
        DOUBLE PRECISION BIOMASS2GM    ! factor to assist getting grams 
        LOGICAL CHK_SUM_LUSE
 
        LOGICAL       :: EFLAG = .FALSE.           ! Error flag
        CHARACTER(16) :: PROGNAME = 'NORMBEIS4'  ! Program name

C***********************************************************************
C   begin body of program NORMBEIS4

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and
C           prompt to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get the BEIS4 model version to use
        MESG = 'Version of BEIS4 to use'
        CALL ENVSTR( 'BEIS_VERSION', MESG, '4.0', BEISVER, IOS )

C.........  Get file name; open emission factors file
        FDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION FACTORS file',
     &           .TRUE., .TRUE., 'BEISFAC', PROGNAME )

C.........  Open gridded landuse file 
        LUNAME = PROMPTMFILE(
     &           'Enter logical name for GRIDDED LANDUSE totals file',
     &           FSREAD3, 'BELD6', PROGNAME )

        IF ( .NOT. DESC3( LUNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( LUNAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        NBELD = NVARS3D

C.........  Initialize grid definition
        CALL CHKGRID( LUNAME, 'GRID' , 0, EFLAG )

C.........  Grid cell resolution assumed to be in meters
C           Input land use will be in percentages
C           Compute conversion factor to go from percentages
C           to km**2

        NCOLS = NCOLS3D
        NROWS = NROWS3D
        GRDNM = GDNAM3D

C.......  factors to arrive at desired units
C.......  for BELD landuse convert from percent to km**2
C.......  for BIOMASS  g/m**2 multiplied by emis fac in ug/(g-hr)

        PRCNT2KM2  = XCELL3D * YCELL3D * 1E-08
        BIOMASS2GM = XCELL3D * YCELL3D * 1E-06


C.........  Store landuse variable names from first file

        ALLOCATE( VGLIST( NBELD ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VGLIST', PROGNAME )

        NVEG = NBELD

        ALLOCATE( BUNITS( NBELD ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BUNITS', PROGNAME )
        
        DO I = 1, NBELD
           BUNITS ( I ) = TRIM( UNITS3D ( I ) )
        END DO

C.........  Allocate memory for emission factor variables   

        ALLOCATE( VEGID ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VEGID', PROGNAME )

        DO I = 1, NBELD
           VGLIST ( I ) = TRIM( VNAME3D ( I ) )
        END DO
        DO I = 1, NBELD
           VEGID ( I ) = VGLIST ( I )
        END DO

C.........  Open gridded biomass file
        BNAME = PROMPTMFILE(
     &           'Enter logical name for GRIDDED LANDUSE totals file',
     &           FSREAD3, 'BIOMASS', PROGNAME )

        IF ( .NOT. DESC3( BNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             TRIM( BNAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        NBIOM = NVARS3D

C.............  Check grid definition
        CALL CHKGRID( BNAME, 'GRID' , 0 , EFLAG )
        EFLAG = .FALSE.

        IF ( EFLAG ) THEN
           MESG = 'Grid in file "' // TRIM( BNAME ) //
     &            '" does not match previously set grid.'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for biomass variable names from BIOMASS file
        ALLOCATE( BIOM_ID( NBIOM ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BIOM_ID', PROGNAME )

        DO I = 1, NBIOM
           BIOM_ID ( I ) = TRIM( VNAME3D ( I ) )
        END DO

C.............. Check variables names BELD vs BIOMASS files

        IF ( NBELD .NE. NBIOM ) THEN
           WRITE( MESG,94010 ) 'Number of veg types in BELD ', 
     &       NBELD,
     &       ' does not match number of veg types in BIOMASS ', 
     &       NBIOM
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C............  Ensure variable names the same in BELD and BIOMASS file

        EFLAG = .FALSE.
        DO I = 1, NBELD
           IFOUND = 0 
           DO N = 1, NBIOM
             IF ( VEGID( I ) .EQ. BIOM_ID( N ) ) IFOUND = 1
           ENDDO
           IF ( IFOUND .EQ. 0 ) THEN
             MESG = VEGID ( I ) // 'not found in BIOMASS file' 
             CALL M3MESG( MESG )
             EFLAG = .TRUE.
           ENDIF
        ENDDO

        IF ( EFLAG ) THEN
           MESG = "'ERROR: Variable names not identical in BIOMASS and ',
     &            'BELD files.'"
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
 
C.. Add check for units of Landuse and Biomass files to avoid using the same
        EFLAG = .FALSE.
        DO I = 1, NBELD
                IFOUND=0
                IF ( UNITS3D( I ) .NE. "g m-2           " ) IFOUND =1
                        IF ( IFOUND .EQ. 1 ) THEN
                MESG = VEGID ( I ) // 'Biomass file has incorrect units'
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                ENDIF
        ENDDO   
        
        IF ( EFLAG ) THEN
           MESG = "'ERROR: Incorrect Biomass File',
     &            'BELD files.'"
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.. Add check for units of Landuse and Biomass files to avoid using the same
        EFLAG = .FALSE.
        DO I = 1, NBIOM
                IFOUND=0
                IF ( BUNITS( I ) .NE. "percent         " ) IFOUND =1
                        IF ( IFOUND .EQ. 1 ) THEN
                MESG = VEGID ( I ) // 'Landuse file has incorrect units'
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                ENDIF
        ENDDO

        IF ( EFLAG ) THEN
           MESG = "'ERROR: Incorrect Landuse File',
     &            'BELD files.'"
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF


C.........  Set up header variables for output file BEIS_NORM_EMIS
        NROWS3D = NROWS
        NCOLS3D = NCOLS
        GDNAM3D = GRDNM
        FTYPE3D = GRDDED3

        SDATE3D = 0       !  n/a
        STIME3D = 0       !  n/a
        TSTEP3D = 0       !  time independent
        NVARS3D = ( (NSEF-1) + NLAI ) * NSEASONS + NNO    ! special treatment of NO
        NLAYS3D = 1
        NTHIK3D = 1
        VGTYP3D = IMISS3
        VGTOP3D = AMISS3

        FDESC3D = ' '   ! array

        FDESC3D( 1 ) = 'BEIS4 normalized emissions values.'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // BEISVER 
        FDESC3D( 4 ) = '/LANDUSE/ BELD6 '                 
C........... BELD6 only works with BEISv4 or later

        I = 0

C.........  Set up variable names and output units
        DO M = 1, NSEASONS 
            DO B = 1, NSEF
            
                BTMP = BIOTYPES( B ) 

C.................  Handle types except NO
                IF( TRIM( BTMP ) /= 'NO' ) THEN
                    I = I + 1
                    
                    VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // SEASON( M )
                    VDESC3D( I ) = 'normalized emissions'
                    UNITS3D( I ) = 'gramsC/hour'
                    VTYPE3D( I ) = M3REAL
                END IF
            END DO
            
            DO N = 1, NLAI
                
                BTMP = LAITYPES( N )
                I = I + 1
                
                VNAME3D( I ) = 'LAI_' // TRIM( BTMP ) // SEASON( M )
                VDESC3D( I ) = 'normalized emissions'
                UNITS3D( I ) = 'index'
                VTYPE3D( I ) = M3REAL
            END DO
        END DO

C.............  Handle NO types (not dependent on season)            
        DO B = 1, NSEF
        
            BTMP = BIOTYPES( B )
            
            IF( TRIM( BTMP ) == 'NO' ) THEN
            
                I = I + 1
                VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // 'AG_GROW'
                VDESC3D( I ) = 'normalized emissions for NO AG_GROW'
                UNITS3D( I ) = 'gramsN/hour'
                VTYPE3D( I ) = M3REAL
                
                I = I + 1
                VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // 'AG_NONGROW'
                VDESC3D( I ) = 'normalized emissions for NO AG_NONGROW'
                UNITS3D( I ) = 'gramsN/hour'
                VTYPE3D( I ) = M3REAL

                I = I + 1
                VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // 'NONAG'
                VDESC3D( I ) = 'normalized emissions for NO NONAG'
                UNITS3D( I ) = 'gramsN/hour'
                VTYPE3D( I ) = M3REAL

            END IF
        END DO

C.........  Open output file
        ENAME = PROMPTMFILE(  
     &        'Enter logical name for NORMALIZED emissions output file',
     &        FSUNKN3, 'BEIS_NORM_EMIS', PROGNAME )

C.........  Get length of BFAC file
        NLINES = GETFLINE( FDEV, 'Emissions factor file' )

C.........  Allocate memory for emission factor variables   
        ALLOCATE( EMFAC ( NVEG, NSEF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMFAC', PROGNAME )

        ALLOCATE( LAI ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LAI', PROGNAME )

        ALLOCATE( WFAC ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WFAC', PROGNAME )

        ALLOCATE( VAG_YN ( NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VAG_YN', PROGNAME )

C...... iniital values for checking
        EMFAC(1:NVEG,1:NSEF) = -99.0
        LAI(1:NVEG) = -99
        WFAC(1:NVEG) = -99.0
        VAG_YN(1:NVEG) = .TRUE.
 
C.........  Read emissions factor file
        MESG = 'Reading emissions factor file...'
        CALL M3MSG2( MESG )
        
        WRITE( MESG,94010 ) 'Number of landuse types in factor file: ',
     &                      NVEG
        CALL M3MSG2( MESG )

C.........  This routine reads in emission factors for BEIS4 normalized
C.........  emissions calculations

        CALL RDNORMBEIS4_EFS( NLINES, NSEF, FDEV, NBELD, 
     &               VGLIST, BIOTYPES, LAI, WFAC, EMFAC, VAG_YN ) 
     
        DO I = 1, NVEG
           IF (LAI(I) .eq. -99) THEN
              MESG = 'ERROR: MISSING LAI FOR VEG TYPE: '//VGLIST(I)
              CALL M3MSG2( MESG )
              EFLAG = .TRUE.
           ENDIF

           IF (WFAC(I) .eq. -99.0) THEN
              MESG = 'ERROR: MISSING WINTER FAC FOR VEG TYPE: 
     &        '//VGLIST(I)
              CALL M3MSG2( MESG )
              EFLAG = .TRUE.
           ENDIF
           IF ( VAG_YN( I ) ) THEN
              MESG = VGLIST(I) // 'being treated as agricultural type'
              CALL M3MSG2( MESG )
           ENDIF
           IF ( VEGID( I ) .EQ. 'MODIS_14        ' ) MODIS14 = I
           IF ( VEGID( I ) .EQ. 'Other_Grass     ' )  IGRASS = I
        ENDDO

        DO J = 1, NSEF
           DO I = 1, NVEG
              IF (EMFAC(I,J) .eq. -99.0 ) THEN
                 MESG = 
     &            'ERROR: MISSING EMISSION FACTOR FOR VEG TYPE:'
     &            //VGLIST(I)//
     &            'AND SPECIES: '//BIOTYPES(J)
                 CALL M3MSG2( MESG )
                 EFLAG = .TRUE.    
              ENDIF
           ENDDO
        ENDDO

        IF( EFLAG ) THEN
            MESG = 'ERRORs are occurred above. Check the messages!'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        ALLOCATE( LUSE ( NCOLS, NROWS, NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUSE', PROGNAME )
        ALLOCATE( BIOM ( NCOLS, NROWS, NVEG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BIOM', PROGNAME )

C.........  Read the gridded landuse and biomass from the input files
        DO M = 1, NVEG
     
           MESG = 'Reading in landuse and biomass ' // VEGID( M )
           CALL M3MESG( MESG )

           IF( .NOT. READ3( LUNAME, VEGID( M ), 1, 0, 0, 
     &                               LUSE( 1, 1, M ) ) ) THEN
                MESG = 'Could not find "' // VEGID( M ) // 
     &                 '" in file "' // TRIM( LUNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
           END IF

           IF( .NOT. READ3( BNAME, VEGID( M ), 1, 0, 0,
     &                               BIOM( 1, 1, M ) ) ) THEN
                MESG = 'Could not find "' // VEGID( M ) //
     &                 '" in file "' // TRIM( BNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
           END IF

        END DO

C.........  Allocate memory for output normalized fluxes
        ALLOCATE( AVGEMIS( NCOLS, NROWS, NSEF, NSEASONS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGEMIS', PROGNAME )

        ALLOCATE( NOEMIS( NCOLS, NROWS, NNO ), STAT=IOS )
        CALL CHECKMEM( IOS, ' NOEMIS', PROGNAME )

        ALLOCATE( AVGLAI( NCOLS, NROWS, NLAI, NSEASONS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AVGLAI', PROGNAME )

        ALLOCATE( SUMEM( NSEF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMEM', PROGNAME )
       
        ALLOCATE( NOEM( NNO ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NOEM', PROGNAME )

        ALLOCATE( SUMEMW( NSEF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMEMW', PROGNAME )

        ALLOCATE( SUMLAI( NLAI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMLAI', PROGNAME )

        ALLOCATE( SUMLAIW( NLAI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUMLAIW', PROGNAME )

        AVGEMIS = 0.000  !  array
        AVGLAI  = 0.000  !  array
        NOEMIS  = 0.000  !  array
        LAI_SAVE_INDEX(1:3) = 0

        MESG = 'Check sum of landuse?'
        CHK_SUM_LUSE = ENVYN ( 'CHECK_LUSE', MESG, .TRUE., IOS )

        IF ( CHK_SUM_LUSE ) THEN

         MESG = 'Checking landuse sums...'
         CALL M3MESG( MESG )

         EFLAG = .FALSE. 
         DO J = 1, NROWS
          DO I = 1, NCOLS
            SVEGA = SUM( LUSE( I, J, : ) )
            IFOUND=0
            IF ( SVEGA  < 99.99  ) IFOUND =1
            IF ( SVEGA  > 100.01 ) IFOUND =1
            IF ( IFOUND .EQ. 1 ) THEN
              WRITE( MESG,94020 ) 'Incorrect BELD sums at ncol',
     &           I, ' and nrow ', J, ' sum= ', SVEGA
              CALL M3MESG( MESG )
              EFLAG = .TRUE.
            ENDIF
 
          ENDDO
         ENDDO

         IF ( EFLAG ) THEN
           MESG = 'BELD landuse did not sum to 100%'
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
         ENDIF

        ENDIF

C.........  Calculate normalized fluxes 
        DO J = 1, NROWS
           DO I = 1, NCOLS

C.................  Initialize variables
              SUMEM = 0.000   ! array
              SUMEMW = 0.000  ! array
              NOEM  = 0.000   ! array
              SUMLAI = 0.000  ! array
              SUMLAIW = 0.000 ! array

              DO M = 1, NVEG

C.....................  Assuming that land use is in percentages
C                       Converting to units of area (km**2)
                    VEGAREA = LUSE( I, J, M ) * PRCNT2KM2
                    VEGBIOM = BIOM( I, J, M ) * BIOMASS2GM  

                    DO N = 1, NSEF
                        BTMP = BIOTYPES( N ) 

C.........................  Special handling for NO emissions
                        IF ( TRIM( BTMP ) == 'NO' ) THEN

                            IF ( VEGAREA > 0.0000 ) THEN

                               IF ( VAG_YN( M ) ) THEN
     
C............... Compute NO emissions for agriculture regions 
C                during growing season

C............... Mosaics of small-scale cultivation 40-60% with
C................natural tree, shrub, or herbaceous vegetation.
C................(from latest MODIS doc for MODIS_14 at 
C   https://lpdaac.usgs.gov/documents/101/MCD12_User_Guide_V6.pdf 
C................ using 50% below

                                 IF( M .EQ. MODIS14 ) THEN
                                   NOEM( 1 ) = NOEM( 1 ) 
     &                                + VEGAREA * EMFAC(M,N)*0.500
                                   NOEM( 2 ) = NOEM( 2 )
     &                                + VEGAREA * EMFAC(IGRASS,N)*
     &                                    0.500
                                   NOEM( 3 ) = NOEM( 3 ) 
     &                                + VEGAREA * EMFAC(M,N)*0.500
                                 ELSE
                                   NOEM( 1 ) = NOEM( 1 )
     &                                + VEGAREA * EMFAC(M,N)
C...................  Compute NO emissions for agriculture regions 
C                                       outside growing season
                                   NOEM( 2 ) = NOEM( 2 )
     &                                + VEGAREA * EMFAC(IGRASS,N) 
                                 ENDIF

                                ELSE
C.....................................  Compute NO emissions for Non-Agriculture regions 
                                 NOEM( 3 ) = NOEM( 3 ) 
     &                               + VEGAREA * EMFAC(M,N)

                                END IF
                            END IF

                        ELSE
C........................Normalized BVOC calculations
                            EFTMP = EMFAC( M, N )

                            IF ( VEGBIOM > 0.000 ) THEN

C.................................  Compute summer emissions
                                SUMEM( N ) = SUMEM( N ) 
     &                               + VEGBIOM * EFTMP

C.................................  Compute winter emissions
                                SUMEMW( N ) = SUMEMW( N )
     &                               + VEGBIOM * EFTMP * WFAC( M )
                            END IF

C.............................  Compute LAI on ISOP and MBO and METH

                            IF ( TRIM( BTMP ) == 'ISOP') THEN
                               LAI_SAVE_INDEX(1) = N
                                SUMLAI( 1 ) = SUMLAI( 1 ) + VEGBIOM
     &                                * LAI( M ) * EFTMP
                                SUMLAIW( 1 ) = SUMLAIW( 1 ) + VEGBIOM
     &                                * LAI( M ) * EFTMP * WFAC( M )
                            END IF

                            IF ( TRIM( BTMP ) == 'MBO' ) THEN
                                LAI_SAVE_INDEX(2) = N

                                SUMLAI( 2 ) = SUMLAI( 2 ) + VEGBIOM
     &                                * LAI( M ) * EFTMP
                                SUMLAIW( 2 ) = SUMLAIW( 2 ) + VEGBIOM
     &                                * LAI( M ) * EFTMP * WFAC( M )
                            END IF

                            IF ( TRIM( BTMP ) == 'METH' ) THEN
                                LAI_SAVE_INDEX(3) = N

                                SUMLAI( 3 ) = SUMLAI( 3 ) + VEGBIOM
     &                                * LAI( M ) * EFTMP
                                SUMLAIW( 3 ) = SUMLAIW( 3 ) + VEGBIOM
     &                                * LAI( M ) * EFTMP * WFAC( M )
                            END IF
                            
                        END IF  ! check if NO emissions
                    END DO  ! end of emis fac loop
                END DO  ! end of veg land use loop2


                DO K = 1, NLAI

                    IF ( SUMLAI( K ) <= 1E-06 ) THEN
                        AVGLAI( I, J, K, 1 ) = 0.
                    ELSE IF ( SUMEM(  LAI_SAVE_INDEX(K) ) == 0. ) THEN
                        AVGLAI(  I, J, K, 1 ) = 0.
                    ELSE
                        AVGLAI( I, J, K, 1 ) =
     &                  SUMLAI( K )/SUMEM( LAI_SAVE_INDEX(K) )
                    ENDIF
                    
                    IF ( SUMLAIW( K ) <= 1E-06 ) THEN
                        AVGLAI( I, J, K, 2 ) = 0.
                    ELSE IF ( SUMEMW(  LAI_SAVE_INDEX(K) ) == 0. ) THEN
                        AVGLAI(  I, J, K, 2 ) = 0.
                    ELSE
                        AVGLAI( I, J, K, 2) =
     &                  SUMLAIW( K ) /SUMEMW( LAI_SAVE_INDEX(K) )
                    ENDIF

                END DO

                DO N = 1, NSEF

                    BTMP = BIOTYPES( N ) 

C.....................  Check for NO emissions
                    IF ( TRIM( BTMP ) == 'NO' ) THEN    
                        NOEMIS( I, J, 1 ) = NOEM( 1 ) 
                        NOEMIS( I, J, 2 ) = NOEM( 2 )
                        NOEMIS( I, J, 3 ) = NOEM( 3 )
                    ELSE
                        AVGEMIS( I, J, N, 1 ) = SUMEM( N ) 
                        AVGEMIS( I, J, N, 2 ) = SUMEMW( N )
                    END IF

                END DO  ! end loop over emission factors
                END DO  ! end loop over rows
        END DO  ! end loop over columns



C.........  Write output file
        I = 0
        DO M = 1, NSEASONS
            DO B = 1, NSEF
                BTMP = BIOTYPES( B ) 

C.................  Handle types other than NO
                IF ( TRIM( BTMP ) /= 'NO' ) THEN  

                    I = I + 1 
                    IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                                 AVGEMIS(1, 1, B, M ) ) ) THEN
                        MESG = 'Could not write "' //
     &                          TRIM( VNAME3D( I ) ) //
     &                          '" to "' // TRIM( ENAME ) // '"'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF
                END IF

            END DO

            DO N = 1, NLAI

                I = I + 1
                IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                             AVGLAI(1, 1, N, M ) ) ) THEN
                    MESG = 'Could not write "' //
     &                      TRIM( VNAME3D( I ) ) //
     &                      '" to "' // TRIM( ENAME ) // '"'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END DO

        END DO

        I = I + 1 
        IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                     NOEMIS(1, 1, 1 ) ) ) THEN
            MESG = 'Could not write "' //
     &              TRIM( VNAME3D( I ) ) //
     &              '" to "' // TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        I = I + 1 
        IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                     NOEMIS(1, 1, 2 ) ) ) THEN
            MESG = 'Could not write "' //
     &              TRIM( VNAME3D( I ) ) //
     &              '" to "' // TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        I = I + 1 
        IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0,
     &                     NOEMIS(1, 1, 3 ) ) ) THEN
            MESG = 'Could not write "' //
     &              TRIM( VNAME3D( I ) ) //
     &              '" to "' // TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  End of program
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

C.........  End of program
        STOP

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )
94020   FORMAT( 2 ( A, :, I5, :, 2X ), A, F10.5 )

C-----------------------------------------------------------------------------
        
        END PROGRAM NORMBEIS4
