
        PROGRAM LAYALLOC

C***********************************************************************
C  program body starts at line 
C
C  DESCRIPTION:
C
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       CRLF, DSCGRID, INDEX1, 
C       PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL,
C       TRIMLEN, ENVINT
C
C  REVISION  HISTORY:
C
C***************************************************************************
C
C COPYRIGHT (C) 2008, Institute for the Environment, UNC-CH
C All Rights Reserved
C
C
C*************************************************************************


        IMPLICIT NONE

C...........   INCLUDES:
        

        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        LOGICAL         DSCGRID
        INTEGER         ENVINT
        REAL            ENVREAL
        INTEGER         FIND1
        INTEGER         FINDC
        INTEGER         GETFLINE
        INTEGER         INDEX1
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN

        EXTERNAL        CRLF, DSCGRID, ENVINT, FIND1, FINDC, 
     &                  GETFLINE, INDEX1, 
     &                  PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL,
     &                  TRIMLEN, WKDAY
C
C...........   PARAMETERS:


C...........   LOCAL PARAMETERS
        CHARACTER*50, PARAMETER :: CVSW = '$Name SMOKEv5.2.1_Sep2025$' ! CVS release tag
        REAL,         PARAMETER :: CONVPA = 1.0E-2  ! conversion factor for Pa to mb

C...........   LOCAL VARIABLES and their descriptions:


        CHARACTER*80            GDESC    !  grid description 

C...........   Logical file names and unit numbers
   
        INTEGER         LDEV    !  log-device
        INTEGER         FDEV       ! input file for user-defined vertical distribution

        CHARACTER*16    SNAME      ! grid-point src met file
        CHARACTER*16    XNAME      ! cross-point layered met file
        CHARACTER*16    ONAME      ! Output factor file
        CHARACTER*16    CNAME      ! Input emission file

C...........   Allocatable arrays
        REAL     , ALLOCATABLE :: UBOT( : )      ! user-defined layer bottom
        REAL     , ALLOCATABLE :: UTOP( : )      ! user-defined layer top
        REAL     , ALLOCATABLE :: UFRAC( : )     ! user-defined layer fraction
        REAL     , ALLOCATABLE :: LFRAC( : )     ! model layer fractions 
        REAL     , ALLOCATABLE :: ZZF( :,:,: )   ! layer's full height (m)
        REAL     , ALLOCATABLE :: VGLVLS( : )    ! gridded mask values to be output
        REAL     , ALLOCATABLE :: TMPBUF( :,:,: )  ! Tmp buffer
        REAL     , ALLOCATABLE :: TMP3D ( :,:,: )
        
C...........   Other local variables
 
        INTEGER         I, J, K, M, N, C, F, S, R, NG, L, NL ! counters, subscripts
        INTEGER         IOS                     !  I/O status
        INTEGER         KREC,IREC                    !  input line (record) number
        INTEGER         JHR, JMAX 
        INTEGER         IHR
        INTEGER         NVARS
        INTEGER         NFILES
        INTEGER         NSTEPS
        INTEGER         NCOLS
        INTEGER         NROWS
        INTEGER         NLAYS
        INTEGER         ULAYS           ! no of user-defined layers
        INTEGER         VGTYP
        INTEGER         LTOP
        INTEGER         LBOT, DDP

        INTEGER      :: JDATE = 0 ! Julian start date (YYYYDDD)
        INTEGER      :: MDATE = 0 ! Julian start date (YYYYDDD) for met file
        INTEGER      :: JTIME = 0 ! start time (HHMMSS)
        INTEGER      :: MTIME = 0 ! start time (HHMMSS) for met file
        INTEGER         TSTEP     ! output time step

        REAL            VGTOP
        REAL            ZBOT, ZTOP, ZFRAC, TFRAC
        REAL            PFRAC, PDIFF
	REAL            LTOT

        LOGICAL         EFLAG      ! Input error flat
        LOGICAL    ::   FIRSTIME = .TRUE.  ! first time flag

        CHARACTER*5     FFORMAT !  temporary indicator for input formats
        CHARACTER*16    CVAR, TFILE, TLAY, VNM
        CHARACTER*100   LINE
        CHARACTER*256   MESG
        CHARACTER*30    SEGMENT( 5 )

        CHARACTER*16 :: PROGNAME = 'LAYALLOC'   !  program name

***********************************************************************
C   begin body of program

        LDEV = INIT3()

C.........  Write out copywrite, version, web address, header info, and prompt
C           to continue running the program.
        CALL INITEM( LDEV, CVSW, PROGNAME )

C.........  Get input files
        MESG = 'Enter name for CROSS-POINT LAYERED MET file'
        XNAME = PROMPTMFILE( MESG, FSREAD3, 'MET_CRO_3D', PROGNAME )
     
C.........  Read description of 3d file for defining layer structure
        IF ( .NOT. DESC3( XNAME ) ) THEN
            MESG = 'Could not get description of file "' //
     &             XNAME( 1:TRIMLEN( XNAME ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF
        
C.........  Store local layer information
        ALLOCATE( VGLVLS( 0:MXLAYS3 ), STAT= IOS) 
        CALL CHECKMEM( IOS, 'VGLVLS', PROGNAME )

        MDATE  = SDATE3D
        MTIME  = STIME3D
        NSTEPS = MXREC3D
        NLAYS  = NLAYS3D
        VGTYP  = VGTYP3D
        VGTOP  = VGTOP3D
        NCOLS  = NCOLS3D
        NROWS  = NROWS3D
        VGLVLS = VGLVS3D   ! array

C.........  Allocate arrays
        ALLOCATE( ZZF( NCOLS,NROWS,NLAYS ), STAT= IOS) 
        CALL CHECKMEM( IOS, 'ZZF', PROGNAME )
        ALLOCATE( LFRAC( NLAYS ), STAT= IOS) 
        CALL CHECKMEM( IOS, 'LFRAC', PROGNAME )

        ZZF    = 0.0
        LFRAC  = 0.0

C.........  Get input layer fractions file
        MESG = 'Enter logical name for gridded input file'
        FDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'LAYER_FRACTION', 
     &                      PROGNAME )
        
        ULAYS = GETFLINE( FDEV, 'User-defined Vertical Layer fraction file' ) 

C........  Allocate memory for emissions fields
        ALLOCATE( UTOP( ULAYS ), STAT= IOS) 
        CALL CHECKMEM( IOS, 'UTOP', PROGNAME )
        ALLOCATE( UBOT( ULAYS ), STAT= IOS) 
        CALL CHECKMEM( IOS, 'UBOT', PROGNAME )
        ALLOCATE( UFRAC( ULAYS ), STAT= IOS) 
        CALL CHECKMEM( IOS, 'UFRAC', PROGNAME )
        
        UTOP  = 0.0
        UBOT  = 0.0
        UFRAC = 0.0
C.........  Read the layers file and store in array
        NL = 0
        IREC = 0

        DO L = 1, ULAYS

            READ ( FDEV, 93000, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94020)
     &             'I/O error', IOS, 'reading LAYERS file at line',IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Left adjust line
            LINE = TRIM( LINE )

C.............  Skip blank and comment lines
            IF( LINE(1:1) == '#' ) CYCLE

C.............  Get line
            CALL PARSLINE( LINE, 5, SEGMENT )
            NL = NL + 1
            UBOT ( NL ) = STR2REAL( SEGMENT( 2 ) )
            UTOP ( NL ) = STR2REAL( SEGMENT( 3 ) )
            UFRAC( NL ) = STR2REAL( SEGMENT( 4 ) )

        END DO
        
        ULAYS = NL

C.........  Read 2-D emissions file
        CNAME = PROMPTMFILE( 
     &          'Enter name for netCDF 2-d emissions file',
     &          FSREAD3, 'INFILE', PROGNAME )
          
        IF ( .NOT. DESC3( CNAME ) ) THEN
              MESG = 'Could not get description of file "' //
     &        CNAME( 1:TRIMLEN( CNAME ) ) // '"'
              CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

C.........  Allocate arrays
        ALLOCATE( TMP3D( NCOLS3D,NROWS3D,NLAYS ), STAT= IOS) 
        CALL CHECKMEM( IOS, 'TMP3D', PROGNAME )
        ALLOCATE( TMPBUF( NCOLS3D,NROWS3D,1 ), STAT= IOS) 
        CALL CHECKMEM( IOS, 'TMPBUF', PROGNAME )

        TMP3D  = 0.0
        TMPBUF = 0.0

C.........  Multiply layer fracitons by 2-d emissions
        JDATE = SDATE3D
        JTIME = STIME3D
        DO IHR = 1, NSTEPS
        
C.............  Read layer's top height (meter)
            IF( .NOT. READ3( XNAME, 'ZF', -1,
     &                       JDATE, JTIME, ZZF ) ) THEN
                MESG = 'ERROR : Could not read ZF from file '
     &                  // XNAME( 1:TRIMLEN( XNAME ) )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

C.............  Loop through variables
            DO K = 1,NVARS3D

                TMPBUF = 0.0    ! array
                VNM = VNAME3D( K )

                IF( .NOT. READ3( CNAME, VNM, 1, JDATE, JTIME,
     &              TMPBUF ) ) THEN
                     MESG = 'ERROR : Could not read ' //VNAME3D( K ) //
     &                      ' from file ' // CNAME( 1:TRIMLEN( CNAME ) )
                     CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                DO R = 1, NROWS3D

                    DO C = 1, NCOLS3D
                        
C.........................  Compute layer fractions given top an bottom
C                           height of the user-defined layers
                        LFRAC = 0.0
                        TFRAC = 0.0
        
                        DO I = 1, ULAYS

C.............................  User-define layers
                            ZTOP  = UTOP( I )
                            ZBOT  = UBOT( I )
                            ZFRAC = UFRAC( I )
            
C.............................  Sum of user-defined layer fractions
                            TFRAC = TFRAC + ZFRAC
            
                            DO L = 1, NLAYS - 1
                                IF ( ZBOT <= ZZF( C,R,L ) ) THEN
                                    LBOT = L
                                    GO TO  111   ! end loop and skip reset of LBOT
                                END IF
                            END DO
                            LBOT = NLAYS           !  fallback

111                         CONTINUE                !  loop exit:  bottom found at LBOT
 
                            IF ( ZTOP <= ZZF( C,R,LBOT ) ) THEN  !  plume in this layer
 
                                PFRAC = 1.0 * ZFRAC
                                LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
                                LTOP = LBOT

                            ELSE IF( LBOT == NLAYS ) THEN    ! plume above top layer
 
                                PFRAC = 1.0 * ZFRAC
                                LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
                                LTOP = NLAYS
                
                            ELSE                               ! plume crosses layers
 
                                DO L = LBOT + 1, NLAYS
                                    IF ( ZTOP <= ZZF( C,R,L ) ) THEN
                                        LTOP = L
                                        GO TO 222  ! end loop and skip reset of LTOP
                                    END IF
                                END DO
                                LTOP = NLAYS

222                             CONTINUE
 
C.................................  Calculate between layer 
                                PDIFF = ZTOP - ZBOT
                    
C.................................  Calculate a fraction for the bottom layer
                                PFRAC = ( ( ZZF( C,R,LBOT ) - ZBOT )
     &                                  / PDIFF ) * ZFRAC
                                LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
                     
C.................................  Calculate a fraction for the top layer
                                PFRAC = ( ( ZTOP - ZZF( C,R,LTOP-1 ) )
     &                                  / PDIFF ) * ZFRAC
                                LFRAC( LTOP ) = LFRAC( LTOP ) + PFRAC

                                DDP = LTOP - LBOT

                                IF( DDP >= 2 ) THEN

                                    DO L = LBOT+1, LTOP-1 !  layers in plume

                                       PFRAC = 
     &                                    ( ( ZZF(C,R,L) -ZZF(C,R,L-1) )
     &                                     / PDIFF ) * ZFRAC
                                       LFRAC( L ) = LFRAC( L ) + PFRAC

                                    END DO

                                ENDIF

                            END IF          !  if ztop in same layer as zbot, or not

                        ENDDO     ! loop over user-define layers
        
                        IF( TFRAC > 1.001 ) THEN
                            MESG = 'ERROR: The total user-defined layer fractions ' //
     &                             ' can not be greater than 1.0'
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        ENDIF
                       
c......................... Before applying layer fractions make sure that they add to 1.0
                       LTOT = 0.0

                       DO L = 1, NLAYS
                          LTOT = LTOT + LFRAC( L )
                       ENDDO

                       IF ( LTOT < 0.999 ) THEN
                           MESG = 'ERROR: Calculated layer fractions '//
     &                         'are less than 1.0. Emissions will be '//
     &                         'dropped.'
                           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )                  
                       ELSE IF ( LTOT > 1.001 ) THEN
                           MESG = 'ERROR: Calculated layer fractions '//
     &                            'are greater than 1.0. Emissions '//
     &                            'will be increased.'
                           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )                  
                       ENDIF

C.........................  Apply layered fraction values to emissions
                        DO L = 1, NLAYS
                          TMP3D( C,R,L ) = TMPBUF( C,R,1 ) * LFRAC( L )
                        ENDDO

                    ENDDO

                ENDDO

C.............  Get output file name using environment variable
                IF( FIRSTIME ) THEN

C.................  Define top layer for output file
                    NLAYS3D = LTOP
                    VGTYP3D = VGTYP
                    VGTOP3D = VGTOP
                    VGLVS3D = VGLVLS
                    MESG = 'Enter logical name for output file'
                    ONAME = PROMPTMFILE( MESG, FSUNKN3, 'OUTFILE', 
     &                                   PROGNAME )
                    FIRSTIME = .FALSE.

                ENDIF

                IF ( .NOT. WRITE3( ONAME, VNM, JDATE, JTIME, TMP3D ) )
     &                                                              THEN
                    WRITE( MESG, 93000 ) 'Could not write to "'
     &                     // ONAME( 1:TRIMLEN( ONAME ) ) // '".'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ENDIF

C.................  Re-initialize tmp array
                TMP3D = 0.0

            ENDDO

            CALL NEXTIME( JDATE, JTIME, 10000 )

        ENDDO

        IF ( .NOT. CLOSE3( CNAME ) ) THEN
            CALL M3ERR( PROGNAME, 0, 0,
     &                  CNAME( 1 : TRIMLEN( CNAME ) ) // '".', .TRUE. )
        END IF      !  if close3() failed 

        IF ( .NOT. CLOSE3( ONAME ) ) THEN
            CALL M3ERR( PROGNAME, 0, 0,
     &                  ONAME( 1 : TRIMLEN( ONAME ) ) // '".', .TRUE. )
        END IF      !  if close3() failed 

C.........   End of program:
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )


C******************  FORMAT  STATEMENTS   ******************************
   
C...........   Formatted file I/O formats............ 93xxx

93100   FORMAT( I6, 1x, A16, F10.3 )

C...........   Internal buffering formats............ 94xxx

94000   FORMAT( I7 ) 
94020   FORMAT( 10( A, I3.0 ) )
93000   FORMAT( A )
94030   FORMAT( 2( A, F10.3 ), 2( A, I3.0 ), A )

C-----------------------------------------------------------------------

        END PROGRAM LAYALLOC
