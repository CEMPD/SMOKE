            
        PROGRAM TMPAREA

C***********************************************************************
C  program body starts at line 254
C
C  DESCRIPTION:
C       Produces source-level inventory-speciated temporally allocated
C       emissions, using the emissions-source vector produced by EMSAREA
C       and the EPS2-formatted temporal profiles and cross-reference 
C       tables for area sources.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level area source emissions in canonical order
C       produced by EMSAREA
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       RDTPROF(), RDASCC(), MKTMAT()
C       Models-3 I/O.
C
C  REVISION  HISTORY:
C       Prototype  3/95 by CJC.
C
C       Modified   8/95 by CJC -- proper handling of time zones, 
C       daylight saving time, annualized vs. day-specific inventories
C
C       Modified  10/95 by CJC -- Uses TPFLAG, to determine source-dependent
C       applicability of monthly, weekly profiles (applicable iff 
C	mod( tflag( source ), [M|W]PRFAC ) is zero
C
C	Modified 8/96 by CJC:  style changes.  Now uses AIRS-format
c  	ATREF, rather than filtered, cut-down ATREF.
C
C	Modified 12/96 by CJC for OTAG:  now assumes that temporal profiles
C  	are independent of pollutant.  New lookup algorithm goes from XREF
C       to profiles.  Uses RDASCC().
C
C       Modified  2/97 by MRH for additional checks related to having
C       taken out pollutant-specific X-refs in the previous revision.
C       Also modified output formatting in many places.
C
C***********************************************************************
C
C PROJECT TITLE: 
C	Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling System
C
C File Version @(#)$Id$
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

        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         CHKXLIN
        CHARACTER*2     CRLF
        INTEGER         FIND1, FIND2, FIND3
        INTEGER         GETDATE
        INTEGER         GETEFILE
        INTEGER         GETNUM
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        LOGICAL         ISDSTIME
        INTEGER         JULIAN
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDTPROF
        INTEGER         STR2INT
        INTEGER         TIME2SEC
        INTEGER         TRIMLEN
        INTEGER         WKDAY           !  day of week (1...7)

        EXTERNAL  CHKXLIN, CRLF, FIND1, FIND2, FIND3, GETDATE, GETEFILE, 
     &            GETNUM, GETYN, HHMMSS, INDEX1, ISDSTIME, JULIAN,
     &            MMDDYY, PROMPTFFILE, PROMPTMFILE, RDTPROF, STR2INT,
     &            TIME2SEC, TRIMLEN, WKDAY


C...........   PARAMETERS and their descriptions

        CHARACTER*5     BLANK5 
        CHARACTER*16    BLANK16
        PARAMETER     ( BLANK5  = ' ' ,
     &                  BLANK16 = ' '  )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Area Sources input and output arrays

        INTEGER         IFIP ( NASRC )         ! FIP codes for sources
        INTEGER         ASC7 ( NASRC )         ! 7 leading digits of ASC
        INTEGER         ASC3 ( NASRC )         ! 3 trailing digits of ASC
        INTEGER         ZONES( NASRC )         ! time zones
        INTEGER         TFLAG( NASRC )         ! governs profile type usage
        INTEGER         MDEX ( NASRC )  ! monthly       profile subscript
        INTEGER         WDEX ( NASRC )  ! weeky         profile subscript
        INTEGER         DDEX ( NASRC )  ! wkday diurnal profile subscript
        INTEGER         EDEX ( NASRC )  ! wkend-diurnal profile subscript
        REAL            EMIS ( NASRC, NIPOL )  ! input inventory emissions.
        REAL            TOUT ( NASRC, NIPOL )  ! timestepped output emssions

        INTEGER         MREF 	  	! monthly       profile subscript
        INTEGER         WREF 	  	! weeky         profile subscript
        INTEGER         DREF 	  	! wkday diurnal profile subscript

C.......   Actually-occurring ASC table

        INTEGER     NASC
        INTEGER     ASCA7( NASCC )
        INTEGER     ASCA3( NASCC )


C...........   Time-stepping Matrix.  
C.......   TMAT( S,V,H ) is profile coeffs for source S, pollutant V, hour H

        REAL            TMAT( NASRC, 24 )    ! coeff profiles

C...........   Temporal Profile XREF tables:  give monthly, weekly, hourly 
C...........   profile-indexes for lookup into temporal profile tables.
C.......   Before sorting:

        INTEGER         NPRFA
        INTEGER         INDXA( NTRFA )
        INTEGER         FIPFA( NTRFA )
        INTEGER         ID7FA( NTRFA )
        INTEGER         ID3FA( NTRFA )
        INTEGER         MPRFA( NTRFA )
        INTEGER         WPRFA( NTRFA )
        INTEGER         DPRFA( NTRFA )

C.......   After sorting:
C.......   Default no-FIP, no-ACC profile

        INTEGER         MPR00   !  month-in-year
        INTEGER         WPR00   !  day-in-week
        INTEGER         DPR00   !  hour-in-day (diurnal)
        INTEGER         EPR00   !  hour-in-day (diurnal)

C.......   No-FIP, by-ASC profiles

        INTEGER         NPRNF
        INTEGER         ID7NF( NASCC )
        INTEGER         ID3NF( NASCC )
        INTEGER         MPRNF( NASCC )
        INTEGER         WPRNF( NASCC )
        INTEGER         DPRNF( NASCC )
        INTEGER         EPRNF( NASCC )

C.......   No-county, (by-STATE FIP code-and-ASC) profiles

        INTEGER         NPRNC
        INTEGER         FIPNC( NASCC )		!  state code (FIP/1000)
        INTEGER         ID7NC( NASCC )
        INTEGER         ID3NC( NASCC )
        INTEGER         MPRNC( NASCC )
        INTEGER         WPRNC( NASCC )
        INTEGER         DPRNC( NASCC )
        INTEGER         EPRNC( NASCC )

C.......   FIP-ASC-specific profiles

        INTEGER         NPRFC
        INTEGER         FIPFC( NTRFA )
        INTEGER         ID7FC( NTRFA )
        INTEGER         ID3FC( NTRFA )
        INTEGER         MPRFC( NTRFA )
        INTEGER         WPRFC( NTRFA )
        INTEGER         DPRFC( NTRFA )
        INTEGER         EPRFC( NTRFA )

C...........   Temporal Profile and XREF-code arrays

        INTEGER         MONREF( NMCOD )     ! Month:Year   XREF-code table
        INTEGER         WEKREF( NWCOD )     ! Day:Week     XREF-code table
        INTEGER         DIUREF( NDCOD )     ! Hour:Weekday XREF-code table
        INTEGER         ENDREF( NDCOD )     ! Hour:Weekend XREF-code table

        INTEGER         NMON		    ! number of month:year profiles
        INTEGER         NWEK		    ! number of day:week   profiles
        INTEGER         NDIU		    ! number of hr:weekday profiles
        INTEGER         NEND  		    ! number of hr:weekend profiles
        REAL            MONFAC( 12, NMCOD ) ! Month:Year   profile coeffs
        REAL            WEKFAC(  7, NWCOD ) ! Day:Week     profile coeffs
        REAL            DIUFAC( 24, NDCOD ) ! Hour:Weekday profile coeffs
        REAL            ENDFAC( 24, NDCOD ) ! Hour:Weekend profile coeffs

C...........   Other local variables

        LOGICAL         EFLAG   !  error-flag
        LOGICAL         DAYLIT  !  iff ZONES changed to daylight time
        LOGICAL         DSTNOW  !  iff currently daylight time now
        INTEGER         TZONE
        INTEGER         DCNT    !  count of uses for default profiles
        INTEGER         JDATE
        INTEGER         JTIME
        INTEGER         TSTEP
        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         ADEV    ! unit number for actual-ASC file
        INTEGER         XDEV	! ...for cross-reference input file
        INTEGER         PDEV	! ...temporal profile input file
        INTEGER         IOS
        INTEGER         IREC
        INTEGER         H, I, J, K, C0, C1, F, S, T, V
        INTEGER         FIP, ID7, ID3, CODE
        INTEGER         IMON, IWEK, IDIU, IEND
        INTEGER         YEAR, MONTH, DAY, LDATE
        INTEGER         LFIP, LID7, LID3, LCOD   ! From previous line of X-REF
        INTEGER         LMON, LWEK, LDIU         ! From previous line of X-REF
        CHARACTER*14    DTBUF   !  buffer for MMDDYY
        CHARACTER*16    ANAME   !  logical name for area-source  input file
        CHARACTER*16    TNAME   !  logical name for timestepped output file
        CHARACTER*128   LINE    !  buffer for input of X-ref file
        CHARACTER*256   MESG    !  message buffer

        LOGICAL         FIRSTCOD  !  iff haven't yet encounter CODE > 0
        DATA            FIRSTCOD / .TRUE. /

        CHARACTER*10    DAYS( 7 )
        DATA            DAYS
     &          /
     &          'Monday', 'Tuesday', 'Wednesday', 'Thursday', 
     &          'Friday', 'Saturday', 'Sunday'
     &          /

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
     &  'or in format YYYYMMDD, STARTING TIME in format HHMMSS',
     &  '(= 10000*hr + 100*min + sec) ) and DURATION in format HHMMSS.',
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

        MPR00 = IMISS3
        WPR00 = IMISS3
        DPR00 = IMISS3
        DO  7  S = 1, NASCC
            MPRNF( S ) = IMISS3
            WPRNF( S ) = IMISS3
            DPRNF( S ) = IMISS3
            MPRNC( S ) = IMISS3
            WPRNC( S ) = IMISS3
            DPRNC( S ) = IMISS3
7       CONTINUE
        DO  8  S = 1, NTRFA
            MPRFC( S ) = IMISS3
            WPRFC( S ) = IMISS3
            DPRFC( S ) = IMISS3
8       CONTINUE


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


C.......   Get list of output variables (shrink list in file description
C.......   data structure to just those variables desired):

        IF ( .NOT. DESC3( ANAME ) ) THEN
            MESG = 'Could not get description of file "'  // 
     &              ANAME( 1:TRIMLEN( ANAME ) ) // '"'
            CALL M3EXIT( 'TMPAREA', 0, 0, MESG, 2 )
        ELSE IF ( NROWS3D .NE. NASRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  AREA SOURCES file:', NROWS3D,
     &      'program:', NASRC
            CALL M3EXIT( 'TMPAREA', 0, 0, MESG, 2 )
        END IF


C.......   Determine whether input represents peak-day emissions or
C.......   annual total emissions:

        C0 =   1 + INDEX( FDESC3D( 6 ),        '"' )
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
     &           'Enter time zone for output emissions (0 for GMT)' )
        JDATE  = GETDATE( 1988200, 
     &           'Enter simulation starting date (YYYYDDD)|(YYYYMMDD)' )

        JTIME  = GETNUM( 0, 999999, 0, 
     &           'Enter simulation starting time  (HHMMSS)' )
        TSTEP  = 10000
        NSTEPS = GETNUM( 1, 999999, 24,
     &                   'Enter output duration (hours)' )

        DTBUF = MMDDYY( JDATE )
        WRITE( MESG,94050 ) 
     &  'Output Time Zone :', TZONE,           CRLF() // BLANK5 //
     &  '       Start Date:', DTBUF( 1:TRIMLEN( DTBUF ) ) // 
     &                                         CRLF() // BLANK5 //
     &  '       Start Time:', JTIME,'HHMMSS'// CRLF() // BLANK5 //
     &  '       Time Step :', 1    ,'hour'  // CRLF() // BLANK5 //
     &  '       Duration  :', NSTEPS, 'hours'

        CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )

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

        CALL M3MSG2( 'Reading ACTUAL ASC file...' )
        CALL RDASCC ( ADEV, NASCC, NASC, ASCA7, ASCA3 )

        WRITE( MESG,94010 ) 
     &      'Actual-ASC table max:', NASCC, 'used:', NASC
        CALL M3MSG2( MESG )


C...........   Read temporal-profiles file:  3 parts (monthly, weekly, diurnal)

        CALL M3MSG2( 'Reading temporal-profile file...' )

        NMON = RDTPROF( PDEV, 'MONTHLY', 12, NMCOD, MONREF, MONFAC )
        NWEK = RDTPROF( PDEV, 'WEEKLY',   7, NWCOD, WEKREF, WEKFAC )
        NDIU = RDTPROF( PDEV, 'WEEKDAY', 24, NDCOD, DIUREF, DIUFAC )
        NEND = RDTPROF( PDEV, 'WEEKEND', 24, NDCOD, ENDREF, ENDFAC )


C.......   Read temporal-profile cross-reference file.

        CALL M3MSG2( 'Reading temporal-profile XREF file...' )

        EFLAG = .FALSE.
        IREC  =  0
        K     =  0
        LFIP  = IMISS3
        LID7  = IMISS3
        LID3  = IMISS3
        LCOD  = IMISS3
        LMON  = IMISS3
        LWEK  = IMISS3
        LDIU  = IMISS3

21      CONTINUE

            IREC = IREC + 1

C.............  Read in X-REF using LINE variable, so that blank
C.............  values will be converted to negative values by STR2INT
C.............  Also enables the use of CHKXLIN function, and is 
C.............  consistent with the other SMOKE programs which require
C.............  this. MRH.
            READ( XDEV, 93000, END=22, IOSTAT=IOS ) LINE

C............. Check for bad IOS value and elimiate non-area records
            IF( CHKXLIN( LINE, IOS, IREC, 'AREA' ) .NE. 0 ) THEN

                GO TO 21 

            ELSE

                FIP  = STR2INT( LINE(  1: 5 ) )
                ID7  = STR2INT( LINE( 21:27 ) )
                ID3  = STR2INT( LINE( 28:30 ) ) 
                CODE = STR2INT( LINE( 32:36 ) )
                IMON = STR2INT( LINE( 38:40 ) )
                IWEK = STR2INT( LINE( 42:44 ) )
                IDIU = STR2INT( LINE( 45:48 ) )

                IF( ID7 .GT. 0  .AND.  ID3 .GE. 0 ) THEN

                    J = FIND2( ID7, ID3, NASC, ASCA7, ASCA3 )
                    IF( J .LE. 0 ) THEN   ! ASC code not found in ASCC LIST file
                        GO TO 21
                    END IF

                ELSE IF( ID7 .LT. 0  .OR.  ID3 .LT. 0 ) THEN
                    
                    WRITE( MESG,94010 ) 
     &                  'Bad ASC at record', IREC,
     &                  'for FIP', FIP, 
     &                  'and ASC ', ID7, ID3
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                ENDIF

C.................  Check each possible fault with input record in
C.................  an individual IF statement, because multiple problems
C.................  can exist for each line.  Slower yes, but more helpful
C.................  to the users. MRH.
 
                IF ( IMON .LE. 0 ) THEN  ! Bad monthly profile code
 
                    WRITE( MESG,94020 ) 
     &                  'Bad IMON=', IMON, 
     &                  'at record', IREC,
     &                  'for FIP', FIP, 
     &                  'and ASC', ID7, ID3
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF
 
                IF ( IWEK .LE. 0 ) THEN  ! Bad weekly profile code
 
                    WRITE( MESG,94020 ) 
     &                  'Bad IWEK=', IWEK, 
     &                  'at record', IREC,
     &                  'for FIP', FIP, 
     &                  'and ASC', ID7, ID3
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF
 
                IF ( IDIU .LE. 0 ) THEN  ! Bad diurnal profile code
 
                    WRITE( MESG,94020 ) 
     &                  'Bad IDIU=', IDIU, 
     &                  'at record', IREC,
     &                  'for FIP', FIP, 
     &                  'and ASC', ID7, ID3
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF

C.................  If pollutant CODE value is positive, then report that
C.................  it is being ignored. MRH.
                IF ( FIRSTCOD .AND. CODE .GT. 0 ) THEN
                    FIRSTCOD = .FALSE.
                    CALL M3MSG2( 
     &                   'Polluant code(s) in XREF file are ignored' )
                END IF

C.....................  For incorrect x-ref records, go on to next line
                IF( EFLAG ) GO TO 21

C.................  Since not using CODE, there could be multiple x-ref
C.................  entries with same FIP, ID7, and ID3.  If this happens
C.................  and there are difference IMON, IWEK, or IDIU: problem. MRH.
                IF( FIP .EQ. LFIP  .AND.  
     &              ID3 .EQ. LID3  .AND. 
     &              ID7 .EQ. LID7        ) THEN

                    IF( IMON .NE. LMON ) THEN
                        WRITE( MESG,94020 ) 
     &                      'Pollutant-specific monthly prof', IMON, 
     &                      'at line', IREC,
     &                      'not used. Profile', LMON,
     &                      'used instead.'
                        CALL M3MESG( MESG )
                    ENDIF

                    IF( IWEK .NE. LWEK ) THEN
                        WRITE( MESG,94020 ) 
     &                      'Pollutant-spe  weekly profile', IWEK, 
     &                      'at record', IREC,
     &                      'not used. Profile', LWEK,
     &                      'used instead.'
                        CALL M3MESG( MESG )
                    ENDIF

                    IF( IDIU .NE. LDIU ) THEN
                        WRITE( MESG,94020 ) 
     &                      'Pollutant-specific diurnal profile', IDIU, 
     &                      'at record', IREC,
     &                      'not used. Profile', LDIU,
     &                      'used instead.'
                        CALL M3MESG( MESG )
                    ENDIF

                    GO TO 21   ! Don't want to restorw for same FIP, ID7, ID3

                ENDIF

C.....................  For correct x-ref records, store if fits

                K = K + 1
                IF( K .LE. NTRFA ) THEN  !  X-ref storage bounds check

                    INDXA( K ) = K
                    FIPFA( K ) = FIP
                    ID7FA( K ) = ID7
                    ID3FA( K ) = ID3
                    MPRFA( K ) = IMON
                    WPRFA( K ) = IWEK
                    DPRFA( K ) = IDIU

                    LFIP       = FIP
                    LID7       = ID7
                    LID3       = ID3
                    LCOD       = CODE
                    LMON       = IMON
                    LWEK       = IWEK
                    LDIU       = IDIU

                END IF 		!  bounds check OK

            END IF     		!  Area records and pollutant code found

            GO TO 21  		!  To head of read-XREF loop

22      CONTINUE

        NPRFA = K

        WRITE( MESG,94010 ) 
     &      'Number of XREF entries--dimensioned', NTRFA,
     &      'actual', NPRFA
        CALL M3MSG2( MESG )
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 'Bad XREF file', 2 )
        ELSE IF ( K .GT. NTRFA ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 'XREF table overflow', 2 )
        END IF


C.......   Sort and process temporal-profile cross-reference file.

        CALL M3MSG2( 'Processing TEMPORAL XREF file...' )
        CALL SORTI3( NPRFA, INDXA, FIPFA, ID7FA, ID3FA )

        IREC  =  0
        NPRNF =  0
        NPRNC =  0
        NPRFC =  0

        DO  33  K = 1, NPRFA

            J    = INDXA( K )
            FIP  = FIPFA( J )
            ID7  = ID7FA( J )
            ID3  = ID3FA( J )
            MREF = MPRFA( J )
            WREF = WPRFA( J )
            DREF = DPRFA( J )

C...........   Look up profiles for these XREFs:

            IF ( NMON .GT. 0 ) THEN
                IMON = FIND1( MREF, NMON, MONREF )
                IF( IMON .LT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'Cross-reference index not found for MREF=', MREF
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF
            ELSE
                IMON = 0
            END IF

            IF ( NWEK .GT. 0  ) THEN
                IWEK = FIND1( WREF, NWEK, WEKREF )
                IF( IWEK .LT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'Cross-reference index not found for WREF=', WREF
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF
            ELSE
                IWEK = 0
            END IF

            IF( NDIU .GT. 0 ) THEN
                IDIU = FIND1( DREF, NDIU, DIUREF )
                IF( IDIU .LT. 0 .AND. NDIU .GT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'Cross-reference index not found for DREF=', DREF
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF
            ELSE
                IDIU = 0
            END IF

            IF( NEND .GT. 0 ) THEN
                IEND = FIND1( DREF, NEND, ENDREF )
                IF( IEND .LT. 0 ) THEN
                    WRITE( MESG,94010 )
     &              'Cross-reference index not found for EREF=', DREF
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF
            ELSE
                IEND = IDIU
            END IF

            IF ( FIP .EQ. 0 ) THEN !  FIP-independent fallback profiles

                IF ( ID7 .EQ. 0 ) THEN  !  ultimate fallback profiles

                    MPR00 = IMON
                    WPR00 = IWEK
                    DPR00 = IDIU
                    EPR00 = IEND

                ELSE                    !  ASC-dependent-only profiles

                    NPRNF = NPRNF + 1
                    IF ( NPRNF .LE. NASCC ) THEN
                        ID7NF( NPRNF ) = ID7
                        ID3NF( NPRNF ) = ID3
                        MPRNF( NPRNF ) = IMON
                        WPRNF( NPRNF ) = IWEK
                        DPRNF( NPRNF ) = IDIU
                        EPRNF( NPRNF ) = IEND
                    END IF

                END IF

            ELSE IF ( MOD( FIP,1000 ) .EQ. 0 ) THEN     !  state-ASC dependent

                NPRNC = NPRNC + 1
                IF ( NPRNC .LE. NASCC ) THEN
                    FIPNC( NPRNC ) = FIP / 1000
                    ID7NC( NPRNC ) = ID7
                    ID3NC( NPRNC ) = ID3
                    MPRNC( NPRNC ) = IMON
                    WPRNC( NPRNC ) = IWEK
                    DPRNC( NPRNC ) = IDIU
                    EPRNC( NPRNC ) = IEND
                END IF

            ELSE                !  FIP-and-ASC dependent profiles

                NPRFC = NPRFC + 1
                IF ( NPRFC .LE. NASCC ) THEN
                    FIPFC( NPRFC ) = FIP
                    ID7FC( NPRFC ) = ID7
                    ID3FC( NPRFC ) = ID3
                    MPRFC( NPRFC ) = IMON
                    WPRFC( NPRFC ) = IWEK
                    DPRFC( NPRFC ) = IDIU
                    EPRFC( NPRFC ) = IEND
                END IF

            END IF      !  if i=0; else if fip=0; else...

33      CONTINUE

        WRITE( MESG,94010 ) 
     &      'ASC-only XREF table max:', NASCC, 'used:', NPRNF
        CALL M3MSG2( MESG )
        WRITE( MESG,94010 ) 
     &      'By-state XREF table max:', NASCC, 'used:', NPRNC
        CALL M3MSG2( MESG )
        WRITE( MESG,94010 ) 
     &      'by-FIP   XREF table max:', NTRFA, 'used:', NPRFC
        CALL M3MSG2( MESG )
        IF ( NPRNF .GT. NASCC  .OR.
     &       NPRNC .GT. NASCC  .OR.
     &       NPRFC .GT. NTRFA ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 'XREF table overflow', 2 )
        END IF

        IF ( MPR00 .LT. 0 ) THEN
            CALL M3WARN( 'TMPAREA', 0, 0, 
     &                   'XREF:  Missing fallback MONTHLY profile' )
        ELSE
            WRITE( MESG,94010 ) 
     &      'Fallback (default)  MONTHLY  profile', MONREF( MPR00 )
            CALL M3MSG2( MESG )
        END IF

        IF ( WPR00 .LT. 0 ) THEN
            CALL M3WARN( 'TMPAREA', 0, 0, 
     &                   'XREF:  Missing fallback WEEKLY profile' )
        ELSE
            WRITE( MESG,94010 ) 
     &      'Fallback (default)   WEEKLY  profile', WEKREF( WPR00 )
            CALL M3MSG2( MESG )
        END IF

        IF ( DPR00 .LT. 0 ) THEN
            CALL M3WARN( 'TMPAREA', 0, 0, 
     &                   'XREF:  Missing fallback DIURNAL profile' )
        ELSE
            WRITE( MESG,94010 ) 
     &      'Fallback (default)  DIURNAL  profile', DIUREF( DPR00 )
            CALL M3MSG2( MESG )
        END IF


C...........   Read in emissions indexing data IFIP, ASC7, ASC3.

        CALL M3MSG2( 'Reading source list from AREA file...' )

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
C...........   Search for FIP-ASC match; else for STATE-ASC match; 
C...........   else for ASC-only match; else use defaults.

        CALL M3MSG2( 'Building source indexes...' )
        EFLAG = .FALSE.
        DCNT  = 0

        DO  88  S = 1, NASRC

            FIP = IFIP( S )
            ID7 = ASC7( S )
            ID3 = ASC3( S )

            F = FIND3( FIP, ID7, ID3, NPRFC, FIPFC, ID7FC, ID3FC )  ! FIP-ASC

            IF ( F .GT. 0 ) THEN

                MDEX( S ) = MPRFC( F )
                WDEX( S ) = WPRFC( F )
                DDEX( S ) = DPRFC( F )
                EDEX( S ) = EPRFC( F )

            ELSE	!  FIP-ASC failed

                T = FIP / 1000
                F = FIND3( T, ID7, ID3, NPRNC, FIPNC, ID7NC, ID3NC ) ! state-ASC

                IF ( F .GT. 0 ) THEN

                    MDEX( S ) = MPRNC( F )
                    WDEX( S ) = WPRNC( F )
                    DDEX( S ) = DPRNC( F )
                    EDEX( S ) = EPRNC( F )

                ELSE		!  state-ASC failed

                    F = FIND2( ID7, ID3, NPRNF, ID7NF, ID3NF ) ! ASC only

                    IF ( F .GT. 0 ) THEN

                        MDEX( S ) = MPRNF( F )
                        WDEX( S ) = WPRNF( F )
                        DDEX( S ) = DPRNF( F )
                        EDEX( S ) = EPRNF( F )

                    ELSE          !  fallback to ASC-independent default

                        MDEX( S ) = MPR00
                        WDEX( S ) = WPR00
                        DDEX( S ) = DPR00
                        EDEX( S ) = EPR00

                        DCNT      = DCNT + 1
                        
                        WRITE( MESG,94040 )
     &                  'Fallback temporal profile used for', 
     &                  FIP, ID7, ID3
                        CALL M3MESG( MESG )

                    END IF      !  if find2() ASC (only) succeeded, or not

                END IF      !  if find3() STATE-and-ASC succeeded

            END IF      !  if find3() FIP-and-ASC succeeded

88      CONTINUE        !  end loop on sources S, constructing TNDX

        IF ( DCNT .GT. 0 ) THEN
            WRITE( MESG, 94010 ) 
     &      'Number of times default profiles used', DCNT
            CALL M3MSG2( MESG )

            IF ( MPR00 .LT. 0  .OR.	!  defaults are missing!
     &           WPR00 .LT. 0  .OR.
     &           DPR00 .LT. 0  ) THEN
                CALL M3MESG( 'Default profiles MISSING !!' )
                EFLAG = .TRUE.
            END IF			!  if defaults missing

        END IF				!  if default-count positive

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPAREA', 0, 0, 'Bad XREF::PROF match.', 2 )
        END IF


C...........   Read in emissions values.

        CALL M3MSG2( 'Reading emissions from AREA file...' )

        DO  99  V = 1, NIPOL

            IF ( .NOT. READ3( ANAME, EINAM( V ), ALLAYS3, 0, 0, 
     &                        EMIS( 1,V ) ) ) THEN
                MESG = 'Error reading "' // 
     &                 EINAM( V )( 1:TRIMLEN( EINAM( V ) ) ) //
     &                 '" from file "' //
     &                 ANAME( 1:TRIMLEN( ANAME ) ) // '".'
                CALL M3EXIT( 'TMPAREA', 0, 0, MESG, 2 )
            END IF

99      CONTINUE


C...........   For each time step, write out the area source emissions values:

        CALL M3MSG2( 'Writing out hourly emissions values...' )

        LDATE  = -1   !  initialization forces TMAT calculation on first iter.
        DAYLIT = .FALSE.

        DO  199  T = 1, NSTEPS

C...........   Construct TMAT -- array of effective composite 
C...........   profile coefficients as appropriate:

            IF ( LDATE .NE. JDATE ) THEN

C...............   If this is a new month, or transitions between
C...............   weekday and weekend, calculate a new TMAT
C...............   Adjust sources' time zones to account for daylight time
C..................   Subtract 1 if date is daylight and ZONES is not
C..................   Add 1 if date is not daylight and ZONES is daylight

                DSTNOW = ISDSTIME( JDATE )
                IF ( DSTNOW .AND. .NOT. DAYLIT ) THEN
 
                    DAYLIT = .TRUE.
                    LDATE  = -1

                    DO 101 S = 1, NASRC
                        ZONES( S ) = ZONES( S ) - 1
  101               CONTINUE
 
                ELSE IF( DAYLIT .AND. .NOT. DSTNOW ) THEN
 
                    DAYLIT = .FALSE.
                    LDATE  = -1
 
                    DO 103 S = 1, NASRC
                        ZONES( S ) = ZONES( S ) + 1
  103               CONTINUE
 
                END IF

                CALL DAYMON( JDATE, MONTH, DAY )  !  get month
                DAY = WKDAY( JDATE )              !  get day-of-week

                IF ( NEND .EQ. 0 .AND. LDATE .LT. 0 ) THEN !  no weekend packet

                    WRITE( *,92030 ) 
     &                  'Processing ', DAYS( DAY ), MMDDYY( JDATE )
                    CALL MKTMAT( NASRC,  JDATE, MONTH,  DAY, 
     &                           TZONE,  ZONES, TFLAG,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           DIUFAC, DDEX, 0,0, TMAT )

                ELSE IF ( DAY .GE. 6 ) THEN     !  day is Saturday or Sunday
                                                !  in SAI's brain-damaged
                                                !  day-of-week conventions.

                    WRITE( *,92030 ) 
     &                  'Processing weekend day ', 
     &                  DAYS( DAY ), MMDDYY( JDATE )
                    CALL MKTMAT( NASRC,  JDATE, MONTH, DAY, 
     &                           TZONE,  ZONES, TFLAG,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           ENDFAC, EDEX, 0,0, TMAT )

                ELSE

                    WRITE( *,92030 ) 
     &                  'Processing midweek day ', 
     &                  DAYS( DAY ), MMDDYY( JDATE )
                    CALL MKTMAT( NASRC,  JDATE, MONTH, DAY, 
     &                           TZONE,  ZONES, TFLAG,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           DIUFAC, DDEX, 0,0, TMAT )

                END IF

                LDATE = JDATE		!  reset LDATE

            END IF      ! if jdate not ldate

            WRITE( *, 94030 ) HHMMSS( JTIME )

            H = 1 + MOD( JTIME / 10000 , 24 )

            DO  198  S = 1, NASRC
            DO  197  V = 1, NIPOL

                TOUT( S,V ) = EMIS( S,V ) * TMAT( S,H )

197         CONTINUE
198         CONTINUE

            IF ( .NOT. WRITE3( TNAME, ALLVAR3, 
     &                         JDATE, JTIME, TOUT ) ) THEN

                MESG = 'Could not write time step to "' //
     &                 TNAME( 1 : TRIMLEN( TNAME ) ) // '".'
                CALL M3EXIT( 'TMPAREA', JDATE, JTIME, MESG, 2 )

            END IF		! if write3() failed

            CALL NEXTIME( JDATE, JTIME, TSTEP )

199     CONTINUE	!  end loop on time steps T

999   CONTINUE          !  exit program

      CALL M3EXIT( 'TMPAREA', 0, 0, 'Normal completion of TMPAREA', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, I10 )

92011   FORMAT( /5X, A, I10, / )

92020   FORMAT( 5X, A, I9 )

92030   FORMAT( /5X, 3 A )

92040   FORMAT( /5X, A, I9,
     &          /5X, A, I9.7, ':', I6.6, ' (', A, ')',
     &          /5X, A, I9  , 2X, A,
     &          /5X, A, I9  , 2X, A, / )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )

93015   FORMAT ( I7, I3 )

93020   FORMAT( I5, 1X, I5, 9X, A7, A3, I6, 3 I4 )

93030   FORMAT( I5.5, 15X, I7.7, I3.3, 1X, I5, 3 I4 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I7, :, 1X ) )

94020   FORMAT( A, I5, 1X, A, I7, 1X, A, I7, 1X, A, I9.7, I3.3 )

94030   FORMAT( 8X, 'at time ', A8 )

94040   FORMAT( A, ' FIP', I7.5, '  ASC', I9.7, I3.3 )

94050   FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,
     &          A, 1X, I2.2, 1X, A, 1X, I2.2, 1X, A   )

        END

