
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM TMPAREA

C***********************************************************************
C  program body starts at line 177
C
C  DESCRIPTION:
C       Produces source-level inventory-speciated temporally allocated
C       emissions, using the emissions-source vector produced by GRDAMAT
C       and the temporal profiles and cross-reference tables for area
C       sources (extracted as the areas-only portions of the AIRS tables
C       with the "stack-ID field removed).
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level area source emissions in canonical order
C       produced by GRDAMAT
C       Sorted, cut-down, only-active-area-sources XREF file.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       FIND1, FIND2, FIND3, GETEFILE, GETNUM, GETYN, DAYMON, 
C       TIME2SEC, TRIMLEN, WKDAY
C
C  REVISION  HISTORY:
C       Prototype  3/95 by CJC.
C       Modified   8/95 by CJC -- proper handling of time zones, 
C       daylight saving time, annualized vs. day-specific inventories
C       Modified  10/95 by CJC -- Uses TPFLAG, to determine source-dependent
C       applicability of monthly, weekly profiles (applicable iff 
C	mod( tflag( source ), [M|W]PRFAC ) is zero
C***********************************************************************

      IMPLICIT NONE

C...........   INCLUDES:
        
        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        INTEGER         FIND1, FIND2, FIND3
        INTEGER         GETEFILE
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         INDEX1
        LOGICAL         ISDSTIME
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDTPROF
        INTEGER         STR2INT
        INTEGER         TIME2SEC
        INTEGER         TRIMLEN
        INTEGER         WKDAY           !  day of week (1...7)
        
        EXTERNAL        FIND1, FIND2, FIND3, GETEFILE, GETNUM, GETYN, 
     &                  INDEX1, ISDSTIME, PROMPTFFILE, PROMPTMFILE, 
     &                  RDTPROF, STR2INT, TIME2SEC, TRIMLEN, WKDAY

                                            
C...........   PARAMETERS and their descriptions

        CHARACTER*16    BLANK16
        INTEGER		MXREF
        PARAMETER     ( BLANK16 = ' ' ,
     &                  MXREF   = NTRFA * NIPOL )
        
                
C...........   LOCAL VARIABLES and their descriptions:
C...........   Area Sources input and output arrays
        
        INTEGER         IFIP ( NASRC )         ! FIP codes for sources
        INTEGER         ASC7 ( NASRC )         ! 7 leading digits of ASC
        INTEGER         ASC3 ( NASRC )         ! 3 trailing digits of ASC
        INTEGER         ZONES( NASRC )         ! time zones
        INTEGER         TFLAG( NASRC )         ! governs profile type usage
        INTEGER         MDEX ( NASRC, NIPOL )  ! monthly       profile subscript
        INTEGER         WDEX ( NASRC, NIPOL )  ! weeky         profile subscript
        INTEGER         DDEX ( NASRC, NIPOL )  ! wkday diurnal profile subscript
        INTEGER         EDEX ( NASRC, NIPOL )  ! wkend-diurnal profile subscript
        REAL            EMIS ( NASRC, NIPOL )  ! input inventory emissions.
        REAL            TOUT ( NASRC, NIPOL )  ! timestepped output emssions
        
        INTEGER         MREF ( NIPOL )  ! monthly       profile subscript
        INTEGER         WREF ( NIPOL )  ! weeky         profile subscript
        INTEGER         DREF ( NIPOL )  ! wkday diurnal profile subscript

C.......   Actually-occurring ASC table

        INTEGER     NASC
        INTEGER     ASCA7( NASCC )
        INTEGER     ASCA3( NASCC )


C...........   Time-stepping Matrix.  
C.......   TMAT( S,V,H ) is profile coeffs for source S, pollutant V, hour H
        
        REAL            TMAT( NASRC, NIPOL, 24 )    ! coeff profiles
      
        
C...........   Temporal Profile XREF tables:  give monthly, weekly, hourly 
C...........   profile-indexes for lookup into temporal profile tables.
C.......   Before sorting:

        INTEGER         NPRFA
        INTEGER         INDXA( MXREF )
        INTEGER         FIPFA( MXREF )
        INTEGER         ID7FA( MXREF )
        INTEGER         ID3FA( MXREF )
        INTEGER         MPRFA( MXREF )
        INTEGER         WPRFA( MXREF )
        INTEGER         DPRFA( MXREF )
        INTEGER         CODEA( MXREF )

C.......   After sorting:
C.......   Default no-FIP, no-ACC profile
        
        INTEGER         MPR00( NIPOL )  !  month-in-year
        INTEGER         WPR00( NIPOL )  !  day-in-week
        INTEGER         DPR00( NIPOL )  !  hour-in-day (diurnal)
        
C.......   No-FIP, by-ASC profiles
        
        INTEGER         NPRNF
        INTEGER         ID7NF( NASCC )
        INTEGER         ID3NF( NASCC )
        INTEGER         MPRNF( NASCC, NIPOL )
        INTEGER         WPRNF( NASCC, NIPOL )
        INTEGER         DPRNF( NASCC, NIPOL )
        
C.......   No-county, (by-STATE FIP code-and-ASC) profiles
        
        INTEGER         NPRNC
        INTEGER         FIPNC( NASCC )
        INTEGER         ID7NC( NASCC )
        INTEGER         ID3NC( NASCC )
        INTEGER         MPRNC( NASCC, NIPOL )
        INTEGER         WPRNC( NASCC, NIPOL )
        INTEGER         DPRNC( NASCC, NIPOL )
        
C.......   FIP-ASC-specific profiles

        INTEGER         NPRFC
        INTEGER         FIPFC( NTRFA )
        INTEGER         ID7FC( NTRFA )
        INTEGER         ID3FC( NTRFA )
        INTEGER         MPRFC( NTRFA, NIPOL )
        INTEGER         WPRFC( NTRFA, NIPOL )
        INTEGER         DPRFC( NTRFA, NIPOL )

C...........   Temporal Profile arrays

        INTEGER         MONREF( NMCOD )     ! Monthly         XREF table
        INTEGER         WEKREF( NWCOD )     ! Weekly          XREF table
        INTEGER         DIUREF( NDCOD )     ! Weekday-diurnal XREF table
        INTEGER         ENDREF( NDCOD )     ! Weekend-diurnal XREF table

        REAL            MONFAC( 12, NMCOD ) ! Monthly         profile coeffs
        REAL            WEKFAC(  7, NWCOD ) ! Weekly          profile coeffs
        REAL            DIUFAC( 24, NDCOD ) ! Weekday-diurnal profile coeffs
        REAL            ENDFAC( 24, NDCOD ) ! Weekend-diurnal profile coeffs

C.......   Numbers of month, week, day, weekday, weekend profiles read

        INTEGER         NMON, NWEK, NDIU, NEND

C...........   Other local variables

        LOGICAL         EFLAG   !  error-flag
        LOGICAL         MFLAG   !  use monthly profiles
        LOGICAL         WFLAG   !  use weekly profiles
        INTEGER         TZONE
        INTEGER         JDATE
        INTEGER         JTIME
        INTEGER         TSTEP
        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         ADEV    ! unit number for actual-ASC file
        INTEGER         XDEV	! ...for cross-reference input file
        INTEGER         YDEV	! ...cross-reference output file
        INTEGER         PDEV	! ...temporal profile input file
        INTEGER         IOS
        INTEGER         IREC
        INTEGER         B, H, I, II, J, K, C0, C1, F, S, T, V
        INTEGER         FIP, PLANT, ID7, ID3, CODE, IMON, IWEK, IDIU
        INTEGER         LID7, LID3
        INTEGER         MONTH, DAY, LDATE

        INTEGER         IBADC       !  Bad pollutant code counter
        INTEGER         BPCODE(20)  !  Bad pollutant code list 
 
        CHARACTER*7     ID7C
        CHARACTER*3     ID3C
        CHARACTER*16    ANAME   !  logical name for area-source      input file
        CHARACTER*16    TNAME   !  logical name for timestepped     output file
        CHARACTER*256   MESG    !  buffer for M3EXIT()


C***********************************************************************
C   begin body of program TMPAREA
        
        LDEV = INIT3()
        
        WRITE( *,92000 ) 
     &  ' ',
     &  'Program TMPAREA to take a sorted source level area source ',
     &  'emissions file, the ASCII temporal profile file, and the ',
     &  'ASCII temporal cross-reference file, and produce temporally',
     &  'allocated (time-stepped) emissions for the user-requested ',
     &  'time period.', 
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>"). The XREF file must',
     &  'have been sorted as indicated, prior to program execution.',
     &  'You will be prompted for which species to allocate -- these',
     &  'should only be the chemical species, and not the FIP-related',
     &  'or ASC-related variables.)', 
     &  ' '
        WRITE( *,92000 ) 
     &  'You will also be prompted for the TIME PERIOD to be covered',
     &  'by this run of the program, as specified by STARTING DATE',
     &  'in format YYYYDDD (= 1000*year + Julian day n0 (1...365,6) )',
     &  'STARTING TIME in format HHMMSS (= 10000*hr + 100*min + sec) )',
     &  'and DURATION in format HHMMSS.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 'Ending program', 2 )
        END IF
        
        DO  2  S = 1, NMCOD
        DO  1  I = 1, 12
            MONFAC( I,S ) = 1.0
1       CONTINUE
2       CONTINUE

        DO  4  S = 1, NWCOD
        DO  3  I = 1, 7
            WEKFAC( I,S ) = 1.0
3       CONTINUE
4       CONTINUE

        DO  6  S = 1, NDCOD
        DO  5  I = 1, 24
            DIUFAC( I,S ) = 1.0
            ENDFAC( I,S ) = 1.0
5       CONTINUE
6       CONTINUE

        DO  9  I = 1, NIPOL
            MPR00( I ) = IMISS3
            WPR00( I ) = IMISS3
            DPR00( I ) = IMISS3
            DO  7  S = 1, NASCC
                MPRNF( S,I ) = IMISS3
                WPRNF( S,I ) = IMISS3
                DPRNF( S,I ) = IMISS3
                MPRNC( S,I ) = IMISS3
                WPRNC( S,I ) = IMISS3
                DPRNC( S,I ) = IMISS3
7           CONTINUE
            DO  8  S = 1, NTRFA
                MPRFC( S,I ) = IMISS3
                WPRFC( S,I ) = IMISS3
                DPRFC( S,I ) = IMISS3
8           CONTINUE
9       CONTINUE


C.......   Get file name; open input area sources, temporal cross-reference,
C.......   and temporal profiles files

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for AREA SOURCE input file',
     &          FSREAD3, 'AREA', 'TMPAREA' )

        ADEV = PROMPTFFILE( 
     &           'Enter logical name for ACTUAL ASC file',
     &           .TRUE., .TRUE., 'ASCC', 'TMPAREA' )

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL XREF file',
     &           .TRUE., .TRUE., 'ATREF', 'TMPAREA' )

        PDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL PROFILES file',
     &           .TRUE., .TRUE., 'ATPRO', 'TMPAREA' )

        YDEV = PROMPTFFILE( 
     &  'Enter logical name for XREF OUTPUT file, or "NONE"',
     &           .FALSE., .TRUE., 'NONE', 'TMPAREA' )


C.......   Get list of output variables (shrink list in file description
C.......   data structure to just those variables desired):
        
        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                  'Could not get description of file "' 
     &                  // ANAME( 1:TRIMLEN( ANAME ) ) // '"',
     &                  2 )
        END IF
        
        
C.......   Determine whether input represents peak-day emissions or
C.......   annual total emissions:

        C0 =   1 + INDEX( FDESC3D( 6 ),           '"' )
        IF ( C0 .EQ. 0 ) THEN
            MESG = 'Bad field [' // FDESC3D( 6 )( 1:24 ) //
     &             '] in file description -- no leading quote'
            CALL M3EXIT( 'TMPAREA', 0, 0, MESG, 2 )
        END IF
        C1 =  INDEX( FDESC3D( 6 )( C0 : 24 ) , '"' )
        IF ( C1 .EQ. 0 ) THEN
            MESG = 'Bad field [' // FDESC3D( 6 )( 1:24 ) //
     &             '] in file description -- no trailing quote'
            CALL M3EXIT( 'TMPAREA', 0, 0, MESG, 2 )
        END IF

C.......   Get time period for output file:

        IF ( ISDSTIME( SDATE3D ) ) THEN
            TZONE = 4
        ELSE 
            TZONE = 5
        END IF
        TZONE  = GETNUM( -12, 12, TZONE, 
     &                   'Enter time zone   (0 for GMT)' )
        JDATE  = GETNUM( 0, 9999999, 1988200, 
     &                   'Enter starting date (YYYYDDD)' )
        JTIME  = GETNUM( 0, 999999, 0, 
     &                   'Enter starting time  (HHMMSS)' )
        TSTEP  = 10000
        NSTEPS = GETNUM( 1, 999999, 24,
     &                   'Enter output duration (hours)' )

C...........   Compute offset of pollutant codes within input VNAME3D:
                
        DO  11  V = 1, NIPOL

            I = INDEX1( EINAM( V ), NVARS3D, VNAME3D )
            IF( I .LE. 0 ) THEN
                CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                       'Could not find variable "' // 
     &                       EINAM( V )( 1:TRIMLEN( EINAM( V ) ) )
     &                       // '" in "' // 
     &                       ANAME( 1:TRIMLEN( ANAME ) ) //  '"',
     &                       2 )
            END IF

            VNAME3D( V ) = EINAM( V )
            UNITS3D( V ) = 'tons/hour'
            VDESC3D( V ) = VDESC3D( I )
            VTYPE3D( V ) = M3REAL

11      CONTINUE

        NVARS3D = NIPOL
        SDATE3D = JDATE
        STIME3D = JTIME
        TSTEP3D = 10000
        WRITE( FDESC3D( 9 ),94010 ) 'Time zone ', TZONE
        
        
C.......   Get file name; open output time-stepped area sources file

        TNAME = PROMPTMFILE( 
     &  'Enter logical name for output TIME-STEPPED EMIS file',
     &          FSUNKN3, 'ATMP', 'TMPAREA' )


C.......   Read the ACTUAL ASC FILE
        
        WRITE( *,92000 ) 
     &       ' ', 'Reading ACTUAL ASC file...', ' '

        IREC  =  0
        I     =  0
        NASC  =  0
        LID3  = -1
        LID7  = -1
        EFLAG = .FALSE.
        
13      CONTINUE        !  head of the SDEV-read loop
        
            READ( ADEV, 93015, END=14, IOSTAT=IOS ) ID7, ID3

            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( *,92010 ) 
     &              'Unit # ', ADEV, 
     &              'IOSTAT ', IOS, 
     &              'Line # ', IREC, ' '
                CALL M3WARN( 'TMPAREA', 0, 0, 
     &              'I/O error reading ACTUAL ASC file.' )

            ELSE IF ( I .GE. NASCC ) THEN 

                EFLAG = .TRUE.
                WRITE( *,92010 ) 
     &              'Record number  ', IREC,
     &              'Dimensioned max', NASCC, ' '
                CALL M3WARN( 'TMPAREA', 0, 0, 
     &              'Max ASC TABLE size NASCC exceeded' )

            ELSE IF ( ( LID7 .GT. ID7 ) .OR.
     &                ( LID7 .EQ. ID7 .AND. LID3 .GE. ID3 ) ) THEN 

                EFLAG = .TRUE.
                WRITE( *,92010 ) 
     &              'Record number  ', IREC
                CALL M3WARN( 'TMPAREA', 0, 0, 
     &              'ASC table out of order' )

            ELSE

                I          = I + 1
                LID7       = ID7
                LID3       = ID3
                ASCA7( I ) = ID7
                ASCA3( I ) = ID3

            END IF

            GO TO  13  !  to head of loop
      
14      CONTINUE        !  end of the SDEV-read loop

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                   'Error reading ACTUAL ASC file.', 2 )
        END IF

        NASC = I


C.......   Read temporal-profile cross-reference file.

        WRITE( *,92000 ) ' ', 'Reading temporal-profile XREF file.'

        EFLAG = .FALSE.
        IREC  = 0
        K     = 0
        IBADC = 0
        
21      CONTINUE
            
            IREC = IREC + 1
            READ( XDEV, 93020, END=22, IOSTAT=IOS ) 
     &              FIP, PLANT, ID7C, ID3C, CODE, IMON, IWEK, IDIU
            
            IF ( IOS .NE. 0 ) THEN

                WRITE( *,92010 ) 'Line number ', IREC, 
     &                           'I/O status  ', IOS
                CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                       'Error reading CROSS-REFERENCE file', 2 )
            ENDIF

C.............  Define pollutant for this record
            IF ( CODE .GT. 0 ) THEN
                I = FIND1( CODE, NIPOL, PCODE )
            ELSE
                I = 0
            END IF

C.............  Check for all reasons to move on to the next record
            IF( PLANT .GT. 0 ) THEN   ! Point source x-ref entry

                GO TO 21

            ELSEIF( ID7C(1:1) .GT. '9' ) THEN  ! Motor vehicle x-ref entry

                GO TO 21

            ELSEIF( I .LT. 0 ) THEN   ! Pollutant code not found

                II = -1
                DO B = 1,IBADC
                    IF( CODE .EQ. BPCODE( B ) ) II = B
                ENDDO

                IF( II .LE. 0 ) THEN

                    IBADC = IBADC + 1
                    BPCODE( IBADC ) = CODE

                    WRITE( MESG,94030 ) 
     &                     'Pollutant code', CODE, 'not used'
                    CALL M3MESG( MESG )

                ENDIF

                GO TO 21

C.............  For AREA x-ref entries, check if needed
            ELSE
           
                ID7 = STR2INT( ID7C ) 
                ID3 = STR2INT( ID3C ) 

                IF( ID7 .GT. 0 .AND. ID3 .GE. 0 ) THEN
                    J = FIND2( ID7, ID3, NASC, ASCA7, ASCA3 )
                ELSE
                    J = 1
                ENDIF

                IF( J .LE. 0 ) THEN   ! ASC code not found in ASCC LIST file

                    GO TO 21

                ELSE IF ( IMON .LE. 0 ) THEN  ! Bad monthly profile code
 
                    WRITE( MESG,94020 ) 
     &                  'Bad IMON=', IMON, 
     &                  'at record', IREC,
     &                  'for FIP', FIP, 
     &                  'and ASC', ID7, ID3
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  21
 
                ELSE IF ( IWEK .LE. 0 ) THEN  ! Bad weekly profile code
 
                    WRITE( MESG,94020 ) 
     &                  'Bad IWEK=', IWEK, 
     &                  'at record', IREC,
     &                  'for FIP', FIP, 
     &                  'and ASC', ID7, ID3
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  21
 
                ELSE IF ( IDIU .LE. 0 ) THEN  ! Bad diurnal profile code
 
                    WRITE( MESG,94020 ) 
     &                  'Bad IDIU=', IDIU, 
     &                  'at record', IREC,
     &                  'for FIP', FIP, 
     &                  'and ASC', ID7, ID3
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  21

C.....................  For correct x-ref records, store if fits
                ELSEIF( K .GE. MXREF ) THEN  !  X-ref storage overflow

                    WRITE( MESG,94030 ) 
     &                     'XREF table overflow at record', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.

                ELSE

                    K          = K + 1
                    INDXA( K ) = K
                    FIPFA( K ) = FIP
                    ID7FA( K ) = ID7
                    ID3FA( K ) = ID3
                    MPRFA( K ) = IMON
                    WPRFA( K ) = IWEK
                    DPRFA( K ) = IDIU
                    CODEA( K ) = I

                ENDIF !  No overflow, ASC found in list, no errors

            ENDIF     !  Area records and pollutant code found

            GO TO 21  !  To head of read-XREF loop
                            
22      CONTINUE
                    
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 'Bad XREF file.', 2 )
        END IF
        
        NPRFA = K


C.......   Sort and process temporal-profile cross-reference file.

        WRITE( *,92000 ) ' ', 'Processing temporal-profile XREF file.'

        CALL SORTI3( NPRFA, INDXA, FIPFA, ID7FA, ID3FA )
        
        EFLAG = .FALSE.
        IREC  = 0
        I     = 0
        NPRNF = 0
        NPRNC = 0
        NPRFC = 0
        
        DO  33  K = 1, NPRFA
        
            J    = INDXA( K )
            FIP  = FIPFA( J )
            ID7  = ID7FA( J )
            ID3  = ID3FA( J )
            IMON = MPRFA( J )
            IWEK = WPRFA( J )
            IDIU = DPRFA( J )
            I    = CODEA( J )

            IF ( FIP .EQ. 0 ) THEN !  FIP-independent fallback profiles
                
                IF ( ID7 .EQ. 0 ) THEN  !  ultimate fallback profiles
                    
                    IF ( I .EQ. 0 ) THEN
                        DO 23  I = 1, NIPOL
                            MPR00( I ) = IMON
                            WPR00( I ) = IWEK
                            DPR00( I ) = IDIU
23                      CONTINUE
                    ELSE
                        MPR00( I ) = IMON
                        WPR00( I ) = IWEK
                        DPR00( I ) = IDIU
                        DO 24  V = 1, NIPOL
                            IF ( MPR00( V ) .EQ. IMISS3 ) 
     &                           MPR00( V ) = IMON
                            IF ( WPR00( V ) .EQ. IMISS3 ) 
     &                           WPR00( V ) = IWEK
                            IF ( DPR00( V ) .EQ. IMISS3 ) 
     &                           DPR00( V ) = IDIU
24                      CONTINUE
                    END IF
                    
                ELSE                    !  ASC-dependent-only profiles
                
                    IF ( NPRNF .EQ. 0   .OR.
     &                   ID7 .NE. ID7NF( NPRNF )  .OR.
     &                   ID3 .NE. ID3NF( NPRNF ) ) THEN
                        
                        NPRNF = NPRNF + 1
                        IF ( NPRNF .GT. NASCC ) THEN
                            EFLAG = .TRUE.
                            WRITE( *,92010 )
     &                              'Line number ', IREC,
     &                              'FIP         ', FIP,
     &                              'ID7         ', ID7,
     &                              'ID3         ', ID3
                            CALL M3MESG( 
     &                           'Subscript NPRNF out of range' )
                            GO TO 33
                        END IF

                        ID7NF( NPRNF ) = ID7
                        ID3NF( NPRNF ) = ID3

                    END IF
                    
                    IF ( I .EQ. 0 ) THEN
                        DO  25  I = 1, NIPOL
                            MPRNF( NPRNF,I ) = IMON
                            WPRNF( NPRNF,I ) = IWEK
                            DPRNF( NPRNF,I ) = IDIU
25                      CONTINUE
                    ELSE
                        MPRNF( NPRNF,I ) = IMON
                        WPRNF( NPRNF,I ) = IWEK
                        DPRNF( NPRNF,I ) = IDIU
                        DO  26  V = 1, NIPOL
                            IF ( MPRNF( NPRNF,V ) .EQ. IMISS3 ) 
     &                           MPRNF( NPRNF,V ) = IMON
                            IF ( WPRNF( NPRNF,V ) .EQ. IMISS3 ) 
     &                           WPRNF( NPRNF,V ) = IWEK
                            IF ( DPRNF( NPRNF,V ) .EQ. IMISS3 )
     &                           DPRNF( NPRNF,V ) = IDIU
26                      CONTINUE
                    END IF

                END IF
            
            ELSE IF ( MOD( FIP,1000 ) .EQ. 0 ) THEN     !  state-ASC dependent
            
                IF ( NPRNC .EQ. 0            .OR.
     &               FIP/1000 .NE. FIPNC( NPRNC ) .OR.
     &               ID7 .NE. ID7NC( NPRNC ) .OR.
     &               ID3 .NE. ID3NC( NPRNC ) ) THEN

                    NPRNC = NPRNC + 1
                    IF ( NPRNC .GT. NASCC ) THEN
                        EFLAG = .TRUE.
                        WRITE( *,92010 )
     &                          'Line number ', IREC,
     &                          'FIP         ', FIP,
     &                          'ID7         ', ID7,
     &                          'ID3         ', ID3
                        CALL M3MESG(
     &                       'Subscript NPRNC out of range' )
                        GO TO 33
                    END IF

                    FIPNC( NPRNC ) = FIP / 1000
                    ID7NC( NPRNC ) = ID7
                    ID3NC( NPRNC ) = ID3

                END IF
                    
                IF ( I .EQ. 0 ) THEN
                    DO 27  I = 1, NIPOL
                        MPRNC( NPRNC,I ) = IMON
                        WPRNC( NPRNC,I ) = IWEK
                        DPRNC( NPRNC,I ) = IDIU
27                  CONTINUE                    
                ELSE
                    MPRNC( NPRNC,I ) = IMON
                    WPRNC( NPRNC,I ) = IWEK
                    DPRNC( NPRNC,I ) = IDIU
                    DO 28  V = 1, NIPOL
                        IF ( MPRNC( NPRNC,V ) .EQ. IMISS3 ) 
     &                       MPRNC( NPRNC,V ) = IMON
                        IF ( WPRNC( NPRNC,V ) .EQ. IMISS3 ) 
     &                       WPRNC( NPRNC,V ) = IWEK
                        IF ( DPRNC( NPRNC,V ) .EQ. IMISS3 )
     &                       DPRNC( NPRNC,V ) = IDIU
28                  CONTINUE                    
                END IF

            ELSE                !  FIP-and-ASC dependent profiles
            
                IF ( NPRFC .EQ. 0             .OR.
     &               FIP .NE. FIPFC( NPRFC )  .OR.
     &               ID7 .NE. ID7FC( NPRFC )  .OR.
     &               ID3 .NE. ID3FC( NPRFC ) ) THEN
                                           
                    NPRFC = NPRFC + 1
                    IF ( NPRFC .GT. NASCC ) THEN
                        EFLAG = .TRUE.
                        WRITE( *,92010 )
     &                          'Line number ', IREC,
     &                          'FIP         ', FIP,
     &                          'ID7         ', ID7,
     &                          'ID3         ', ID3
                        CALL M3MESG(
     &                      'Subscript NPRFC out of range' )
                        GO TO 33
                    END IF

                    FIPFC( NPRFC ) = FIP
                    ID7FC( NPRFC ) = ID7
                    ID3FC( NPRFC ) = ID3

                END IF
                    
                IF ( I .EQ. 0 ) THEN
                    DO 29  I = 1, NIPOL
                        MPRFC( NPRFC,I ) = IMON
                        WPRFC( NPRFC,I ) = IWEK
                        DPRFC( NPRFC,I ) = IDIU
29                  CONTINUE
                ELSE
                    MPRFC( NPRFC,I ) = IMON
                    WPRFC( NPRFC,I ) = IWEK
                    DPRFC( NPRFC,I ) = IDIU
                    DO 30  V = 1, NIPOL
                        IF ( MPRFC( NPRFC,V ) .EQ. IMISS3 )
     &                       MPRFC( NPRFC,V ) = IMON
                        IF ( WPRFC( NPRFC,V ) .EQ. IMISS3 ) 
     &                       WPRFC( NPRFC,V ) = IWEK
                        IF ( DPRFC( NPRFC,V ) .EQ. IMISS3 )
     &                       DPRFC( NPRFC,V ) = IDIU
30                  CONTINUE
                END IF

            END IF      !  if i=0; else if fip=0; else...

33      CONTINUE
                    
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 'Bad XREF file.', 2 )
        END IF


C...........   If fixed-up output file desired, write it:

        IF ( YDEV .GT. 0 ) THEN

            WRITE( YDEV,93030 )
     &      ( 0, 0, 0, 
     &        PCODE( V ), MPR00( V ), WPR00( V ), DPR00( V ), 
     &        V = 1, NIPOL ),
     &      ( ( 0, ID7NF( S ), ID3NF( S ), 
     &          PCODE( V ), MPRNF( S,V ), WPRNF( S,V ),DPRNF( S,V ),
     &          V = 1, NIPOL ), S = 1, NPRNF ),
     &      ( ( 1000 * FIPNC( S ), ID7NC( S ), ID3NC( S ), 
     &          PCODE( V ), MPRNC( S,V ), WPRNC( S,V ),DPRNC( S,V ),
     &          V = 1, NIPOL ), S = 1, NPRNC ),
     &      ( ( FIPFC( S ), ID7FC( S ), ID3FC( S ), 
     &          PCODE( V ), MPRFC( S,V ), WPRFC( S,V ),DPRFC( S,V ),
     &          V = 1, NIPOL ), S = 1, NPRFC )
     
        END IF	!  if ydev>0:  optional temporal xref file being produced


C...........   Read temporal-profiles file:  3 parts (monthly, weekly, diurnal)

        WRITE( *,92000 ) ' ', 'Reading temporal-profile file.'

        NMON = RDTPROF( PDEV, 'MONTHLY', 12, NMCOD, MONREF, MONFAC )
        NWEK = RDTPROF( PDEV, 'WEEKLY',   7, NWCOD, WEKREF, WEKFAC )
        NDIU = RDTPROF( PDEV, 'WEEKDAY', 24, NDCOD, DIUREF, DIUFAC )
        NEND = RDTPROF( PDEV, 'WEEKEND', 24, NDCOD, ENDREF, ENDFAC )

     
C...........   Read in emissions indexing data IFIP, ASC7, ASC3.
        
        WRITE( *,92000 ) ' ', 'Reading source list from AREA file.'

        IF ( .NOT. READ3( ANAME, 'FIP', ALLAYS3, 0, 0, IFIP ) ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                   'Error reading FIP from file' // ANAME, 2)
        END IF

        IF ( .NOT. READ3( ANAME, 'ASC7', ALLAYS3, 0, 0, ASC7 ) ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                   'Error reading ASC7 from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'ASC3', ALLAYS3, 0, 0, ASC3 ) ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                   'Error reading ASC3 from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'ZONES', ALLAYS3, 0, 0, ZONES ) ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                   'Error reading ZONES from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'TPFLAG', ALLAYS3,0,0, TFLAG ) ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                   'Error reading TFLAG from file' // ANAME, 2 )
        END IF


C...........   Construct  -- mapping source number to effective
C...........   composite profile number.
              
        WRITE( *,92000 ) ' ', 'Building source indexes.'
        EFLAG = .FALSE.

        DO  79  S = 1, NASRC
        
            MFLAG = ( MOD( TFLAG( S ), MTPRFAC ) .EQ. 0 )
            WFLAG = ( MOD( TFLAG( S ), WTPRFAC ) .EQ. 0 )
            
            IF ( .NOT. MFLAG ) THEN
                DO  71  V = 1, NIPOL
                    MDEX( S,V ) = 0
71              CONTINUE
            END IF

            IF ( .NOT. WFLAG ) THEN
                DO  72  V = 1, NIPOL
                    WDEX( S,V ) = 0
72              CONTINUE
            END IF

            FIP = IFIP( S )
            ID7 = ASC7( S )
            ID3 = ASC3( S )

C...........   Search for FIP-ASC match; else for STATE-ASC match; 
C...........   else for ASC-only match; else use defaults.
            
            F = FIND3( FIP, ID7, ID3, NPRFC, FIPFC, ID7FC, ID3FC )  ! FIP-ASC

            IF ( F .GT. 0 ) THEN

                DO  73  V = 1, NIPOL
                    MREF( V ) = MPRFC( F,V )
                    WREF( V ) = WPRFC( F,V )
                    DREF( V ) = DPRFC( F,V )
73              CONTINUE

                GO TO  77       !  to end of sources-loop

            END IF      !  if find3() FIP-and-ASC succeeded

            F = FIND3( FIP/1000, ID7, ID3, NPRNC, FIPNC, ID7NC, ID3NC )  ! STATE-ASC
  
            IF ( F .GT. 0 ) THEN

                DO  74  V = 1, NIPOL
                    MREF( V ) = MPRNC( F,V )
                    WREF( V ) = WPRNC( F,V )
                    DREF( V ) = DPRNC( F,V )
74              CONTINUE

                GO TO  77       !  to end of sources-loop

            END IF      !  if find3() STATE-and-ASC succeeded

            F = FIND2( ID7, ID3, NPRNF, ID7NF, ID3NF )    ! ASC search
  
            IF ( F .GT. 0 ) THEN

                DO  75  V = 1, NIPOL
                    MREF( V ) = MPRNF(F,V)
                    WREF( V ) = WPRNF(F,V)
                    DREF( V ) = DPRNF( F,V )
75              CONTINUE

            ELSE                        !  fallback to ASC-independent default

                DO  76  V = 1, NIPOL
                    MREF( V ) = MPR00( V )
                    WREF( V ) = WPR00( V )
                    DREF( V ) = DPR00( V )
76              CONTINUE

            END IF      !  if find2() ASC (only) succeeded, or not

77          CONTINUE    !  check index results:

            DO  78  V = 1, NIPOL

                IF ( MFLAG ) THEN
                    MDEX( S,V ) = FIND1( MREF( V ), NMON, MONREF )
                    IF( MDEX( S,V ) .LT. 0 ) THEN
                        WRITE( MESG,94010 )
     &                  'Cross-reference index not found for MREF=',
     &                   MREF( V )
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF

                IF ( WFLAG ) THEN
                    WDEX( S,V ) = FIND1( WREF( V ), NWEK, WEKREF )
                    IF( WDEX( S,V ) .LT. 0 ) THEN
                        WRITE( MESG,94010 )
     &                  'Cross-reference index not found for WREF=',
     &                   WREF( V )
                        CALL M3MESG( MESG )
                        EFLAG = .TRUE.
                    END IF
                END IF

                DDEX( S,V ) = FIND1( DREF( V ), NDIU, DIUREF )
                IF( DDEX( S,V ) .LT. 0 .AND. NDIU .GT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'Cross-reference index not found for DREF=',
     &               DREF( V )
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF

                EDEX( S,V ) = FIND1( DREF( V ), NEND, ENDREF )
                IF( EDEX( S,V ) .LT. 0 .AND. NEND .GT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'Cross-reference index not found for EREF=',
     &               DREF( V )
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF

78          CONTINUE	!  end index lookup loop on species V

79      CONTINUE        !  end loop on sources S, constructing TNDX
        
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 'Bad XREF::PROF match.', 2 )
        END IF


C...........   Construct TMAT -- array of effective composite 
C...........   profile coefficients
        
C...........   Read in emissions values.
        
        WRITE( *,92000 ) ' ', 'Reading emissions from AREA file.'

        DO  99  V = 1, NIPOL
        
            IF ( .NOT. READ3( ANAME, EINAM( V ), ALLAYS3, 0, 0, 
     &                        EMIS( 1,V ) ) ) THEN
                CALL M3EXIT( 'TMPAREA', 0, 0, 
     &                       'Error reading "' // 
     &                       EINAM( V )( 1:TRIMLEN( EINAM( V ) ) )//
     &                       '" from file "' // 
     &                       ANAME( 1:TRIMLEN( ANAME ) ) // '".', 2 )
            END IF

99      CONTINUE


C...........   For each time step, write out the area source emissions values:
        
        WRITE( *,92000 ) ' ', 'Writing out emissions values...'
        
        LDATE = -1   !  initialization forces TMAT calculation on first iter.
        
        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or transitions between
C...............   weekday and weekend, calculate a new TMAT

         
            IF ( LDATE .NE. JDATE ) THEN

                LDATE = JDATE

                CALL DAYMON( JDATE, MONTH, DAY )  !  get month
                DAY = WKDAY( JDATE )              !  get day-of-week

                IF ( NEND .EQ. 0 ) THEN  !  no weekend packet

                    WRITE( *,92011 ) 'Processing day', DAY
                    CALL MKTMAT( NASRC,  JDATE, MONTH,  DAY, 
     &                           TZONE,  ZONES,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           DIUFAC, DDEX, 0,0, TMAT )

                ELSE IF ( DAY .GE. 6 ) THEN     !  day is Saturday or Sunday
                                                !  in SAI's brain-damaged
                                                !  day-of-week conventions.

                    WRITE( *,92011 ) 'Processing weekend day', DAY
                    CALL MKTMAT( NASRC,  JDATE, MONTH, DAY, 
     &                           TZONE,  ZONES,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           ENDFAC, EDEX, 0,0, TMAT )

                ELSE

                    WRITE( *,92011 ) 'Processing midweek day', DAY
                    CALL MKTMAT( NASRC,  JDATE, MONTH, DAY, 
     &                           TZONE,  ZONES,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           DIUFAC, DDEX, 0,0, TMAT )

                END IF

            END IF      ! if jdate not ldate

            H = 1 + MOD( JTIME / 10000 , 24 )

            DO  198  V = 1, NIPOL
            DO  197  S = 1, NASRC

                TOUT( S,V ) = EMIS( S,V ) * TMAT( S,V,H )

197         CONTINUE
198         CONTINUE

            IF ( .NOT.WRITE3( TNAME, 'ALL', JDATE, JTIME, TOUT ) ) THEN

                CALL M3EXIT( 'TMPAREA', JDATE, JTIME,
     &                       'Could not write timestep to "' // 
     &                       TNAME( 1 : TRIMLEN( TNAME ) ) // '".', 2 )

            END IF

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE

999   CONTINUE          !  exit program
      
      CALL M3EXIT( 'TMPAREA', 0, 0, 'Normal completion of TMPAREA', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, I10 )

92011   FORMAT( /5X, A, I10 )

92020   FORMAT( 5X, A, I9 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )

93015   FORMAT ( I7, I3 )

cmh 93020   FORMAT( I5, I8, I3, I6, 3 I4 )
93020   FORMAT( I5, 1X, I5, 9X, A7, A3, I6, 3 I4 )

93030   FORMAT( I5.5, 15X, I7.7, I3.3, 1X, I5, 3 I4 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, I5, :, 1X, A, :, '  FIP', I7.5, '  ASC', I9.7, I3.3 )

94020   FORMAT( A, I5, 1X, A, I7, 1X, A, I7, 1X, A, I9.7, I3.3 )

94030   FORMAT( 10( A, :, I7, :, 1X ) )

        END

