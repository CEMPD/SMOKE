
C Version "@(#)$Id$ $Source$ $Date$ 

        PROGRAM TMPMOBIL

C***********************************************************************
C  program body starts at line 188
C
C  DESCRIPTION:
C       Produces source-level temporally allocated VMT using the 
C       mobile-source inventory vector produced by RAWMOBIL
C       and the temporal profiles and cross-reference tables for mobile
C       sources (extracted as the mobile-only portions of the AIRS tables
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level mobile source emissions in canonical order
C       produced by GRDAMAT
C       Sorted, cut-down, only-active-mobile-sources XREF file.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; 
C       FIND1, FIND2, FIND3 GETEFILE, GETNUM, GETYN, DAYMON, 
C       TIME2SEC, TRIMLEN, WKDAY
C
C  REVISION  HISTORY:
C       Prototype  1/96 by CJC.
C***********************************************************************

      IMPLICIT NONE


C...........   INCLUDES:

        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'TMDIMS3.EXT'   !  emis temporal modeling parms
        INCLUDE 'MBDIMS3.EXT'   !  mobile-source dimensioning parameters
        INCLUDE 'GRDIMS3.EXT'   !  grid parms
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'CONST3.EXT'    !  physical and mathematical constants


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         FIND1, FIND2, FIND3
        INTEGER         GETEFILE
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         INDEX1
        LOGICAL         ISDSTIME
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         RDMPLIST
        INTEGER         RDTPROF
        INTEGER         STR2INT
        INTEGER         SEC2TIME
        INTEGER         TIME2SEC
        INTEGER         TRIMLEN
        INTEGER         WKDAY           !  day of week (1...7)

        EXTERNAL    FIND1, FIND2, FIND3, GETEFILE, GETNUM, GETYN, 
     &              INDEX1, ISDSTIME, PROMPTFFILE, PROMPTMFILE,
     &              RDMPLIST, RDTPROF, STR2INT, SEC2TIME, TIME2SEC, 
     &              TRIMLEN, WKDAY


C...........   PARAMETERS and their descriptions

        CHARACTER*16    BLANK16
        REAL*8          FTOC
        PARAMETER     ( BLANK16 = ' ' ,
     &                  FTOC    = 5.0D0 / 9.0D0 )


C...........   LOCAL VARIABLES and their descriptions:
C...........   Mobile Sources input and output arrays

        INTEGER         IFIP  ( NMSRC )        ! FIP codes for sources
        INTEGER         IRCLAS( NMSRC )        ! EPA road classification code
        INTEGER         ILINK ( NMSRC )        ! link number or 0 for non-link
        REAL            VMT   ( NMSRC )        ! vehicle miles traveled

        INTEGER         TZONES( NMSRC )        ! time zones
        INTEGER         TFLAG( NMSRC )         ! governs profile type usage
        INTEGER         MDEX ( NMSRC, NVTYPE, NETYPE ) ! monthly       profile subscript
        INTEGER         WDEX ( NMSRC, NVTYPE, NETYPE ) ! weeky         profile subscript
        INTEGER         DDEX ( NMSRC, NVTYPE, NETYPE ) ! wkday diurnal profile subscript
        INTEGER         EDEX ( NMSRC, NVTYPE, NETYPE ) ! wkend-diurnal profile subscript
        REAL            TOUT ( NMSRC, NVTYPE, NETYPE ) ! timestepped output emssions

C...........   Vehicle Mix transform matrix

        REAL            MXMAT( NMSRC, NVTYPE )

C...........   MPLIST emissions factor cross reference tables:

        INTEGER   NMPLIST
        INTEGER   MPFIP( NLIST )
        INTEGER   MPRDT( NLIST )
        INTEGER   MPLNK( NLIST )
        INTEGER   MPREF( 0:23, NLIST )    !  hour&source-specific xrefs
        INTEGER   REFCNT
        INTEGER   REFLIST( NLIST )

        INTEGER   EFDEX( 0:23, NMSRC )        ! subscript into emissions factor tables

C...........   Emissions factor tables:

        REAL    TEMP0		!  starting em factor table temperature
        REAL    TEMP1		!  ending em factor table temperature
        REAL    DTEMP		!  em factor table temperature increment (K)
        REAL    DDTEMP		!  inverse temperature increment (K)
        INTEGER NTEMP		!  number of active elements in temp grid

        REAL	EXH_CO ( MAXTMP, NVTYPE, NLIST )	!  emissions factor
        REAL	EXH_NOX( MAXTMP, NVTYPE, NLIST )	!  emissions factor
        REAL	EXH_VOC( MAXTMP, NVTYPE, NLIST )	!  emissions factor
        REAL	EVP_VOC( MAXTMP, NVTYPE, NLIST )	!  emissions factor
        REAL	RNL_VOC( MAXTMP, NVTYPE, NLIST )	!  emissions factor
        REAL	RST_VOC( MAXTMP, NVTYPE, NLIST )	!  emissions factor

C.......   Temporal profiles:

        INTEGER         NMON, NWEK, NDIU, NEND	!  numbers of profiles

        INTEGER         MONREF( NMCOD )     ! Monthly         XREF table
        INTEGER         WEKREF( NWCOD )     ! Weekly          XREF table
        INTEGER         DIUREF( NDCOD )     ! Weekday-diurnal XREF table
        INTEGER         ENDREF( NDCOD )     ! Weekend-diurnal XREF table

        REAL            MONFAC( 12, NMCOD ) ! Monthly         profile coeffs
        REAL            WEKFAC(  7, NWCOD ) ! Weekly          profile coeffs
        REAL            DIUFAC( 24, NDCOD ) ! Weekday-diurnal profile coeffs
        REAL            ENDFAC( 24, NDCOD ) ! Weekend-diurnal profile coeffs

C...........   Time-stepping Matrix.  
C.......   TMAT( S,V,H ) is profile coeffs for source S, pollutant V, hour H

        REAL            TMAT( NMSRC, NVTYPE, NETYPE, 24 )    ! coeff profiles


C...........   Ungridding matrix:

        INTEGER     NU( NMSRC )
        INTEGER     IU( MXGCOF )
        REAL        CU( MXGCOF )

        COMMON  / UGRIDMAT / NU, IU, CU


C...........   Other local variables

        LOGICAL         EFLAG   !  error-flag
        INTEGER         TZONE
        INTEGER         JYEAR
        INTEGER         FDATE
        INTEGER         JDATE	!  for output file
        INTEGER         JTIME	!  for output file
        INTEGER         MDATE	!  for met input file
        INTEGER         MTIME	!  for met input file
        INTEGER         TSTEP
        INTEGER         NSTEPS
        INTEGER         LDEV
        INTEGER         FDEV    ! unit number for emissions factors xref file
        INTEGER         VDEV    ! unit number for VMT mix file
        INTEGER         XDEV    ! ...for cross-reference input file
        INTEGER         PDEV    ! ...temporal profile input file
        INTEGER         H, HCORR, II, J, KU, M, C0, C1, F, S, T, V
        INTEGER         FIP, RDT, LNK
        INTEGER         MONTH, DAY, LDATE
        REAL            SCR, PP, QQ

        REAL            TA( NGRID )	!  one-layer slice of temperature
        REAL            TASRC           !  temperature, interpolated to srcs

        CHARACTER*16    ANAME   !  logical name for mobile-source     input file
        CHARACTER*16    FNAME   !  logical name for emissions factors input file
        CHARACTER*16    MNAME   !  logical name for met_cro_3d        input file
        CHARACTER*16    UNAME   !  logical name for ungridding-matrix input file
        CHARACTER*16    TNAME   !  logical name for timestepped      output file
        CHARACTER*256   MESG    !  buffer for M3EXIT()


C***********************************************************************
C   begin body of program TMPMOBIL

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program TMPMOBIL to take a sorted source level mobile source',
     &  'emissions file, the ASCII temporal profile file, and the ',
     &  'ASCII temporal cross-reference file, and produce temporally',
     &  'allocated (time-stepped) VMT for each vehicle class for the',
     &  'user-requested time period.',
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
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 'Ending program', 2 )
        END IF

C.......   Get file name; open input mobile sources, temporal cross-reference,
C.......   and temporal profiles files

        ANAME = PROMPTMFILE( 
     &          'Enter logical name for mobile SOURCE input file',
     &          FSREAD3, 'MOBL', 'TMPMOBIL' )

        FNAME = PROMPTMFILE( 
     &          'Enter logical name for EMISSIONS FACTOR input file',
     &          FSREAD3, 'MEFS', 'TMPMOBIL' )

        MNAME = PROMPTMFILE( 
     &          'Enter logical name for MET_CRO_3D input file',
     &          FSREAD3, 'MET_CRO_3D', 'TMPMOBIL' )

        UNAME = PROMPTMFILE( 
     &          'Enter logical name for UNGRIDDING MATRIX input file',
     &          FSREAD3, 'MUGMAT', 'TMPMOBIL' )

        VDEV = PROMPTFFILE( 
     &           'Enter logical name for VMT MIX file',
     &           .TRUE., .TRUE., 'MVMTM', 'TMPMOBIL' )

        FDEV = PROMPTFFILE( 
     &           'Enter logical name for EMISSION FACTORS XREF file',
     &           .TRUE., .TRUE., 'MPLIST', 'TMPMOBIL' )

        XDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL XREF file',
     &           .TRUE., .TRUE., 'MTREF', 'TMPMOBIL' )

        PDEV = PROMPTFFILE( 
     &           'Enter logical name for TEMPORAL PROFILES file',
     &           .TRUE., .TRUE., 'MTPRO', 'TMPMOBIL' )


C.......   Get list of output variables (shrink list in file description
C.......   data structure to just those variables desired):

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
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
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
        END IF
        C1 =  INDEX( FDESC3D( 6 )( C0 : 24 ) , '"' )
        IF ( C1 .EQ. 0 ) THEN
            MESG = 'Bad field [' // FDESC3D( 6 )( 1:24 ) //
     &             '] in file description -- no trailing quote'
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
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
        NSTEPS = GETNUM( 1, 999999, 24,
     &                   'Enter output duration (hours)' )

        TSTEP  = 10000
        JYEAR  = JDATE / 1000

        MDATE = JDATE
        MTIME = JTIME
        CALL NEXTIME( MDATE, MTIME, -TZONE )	!  GMT for met file


C.......   Get file name; open output time-stepped area sources file

        SDATE3D = JDATE
        STIME3D = JTIME
        TSTEP3D = TSTEP
        NROWS3D = NMSRC
        NCOLS3D = 1
        NTHIK3D = IMISS3
        NLAYS3D = NVTYPE
        NVARS3D = 1

        VTYPE3D( 1 ) = M3REAL
        VNAME3D( 1 ) = 'VMT'
        UNITS3D( 1 ) = 'miles/hour'
        VDESC3D( 1 ) = 'Hourly miles traveled by vehicle types'

        FDESC3D( 1 ) = 
     &  'Temporally allocated source level mobile source VMT values.'
        DO  20 V = 1, NVTYPE
            WRITE( FDESC3D( V+1 ), 94010 )
     &      'Layer', V, 'represents vehicle type ' // VTYPE3( V )
20      CONTINUE
        DO  21 V = NVTYPE+2, MXDESC3
            FDESC3D( V ) = ' '
21      CONTINUE

        TNAME = PROMPTMFILE( 
     &  'Enter logical name for output TIME-STEPPED EMIS file',
     &          FSUNKN3, 'MTMP', 'TMPMOBIL' )


C...........   Read temporal-profiles file:  3 parts (monthly, weekly, diurnal)

        WRITE( *,92000 ) ' ', 'Reading temporal-profile file.'

        NMON = RDTPROF( PDEV, 'MONTHLY', 12, NMCOD, MONREF, MONFAC )
        NWEK = RDTPROF( PDEV, 'WEEKLY',   7, NWCOD, WEKREF, WEKFAC )
        NDIU = RDTPROF( PDEV, 'WEEKDAY', 24, NDCOD, DIUREF, DIUFAC )
        NEND = RDTPROF( PDEV, 'WEEKEND', 24, NDCOD, ENDREF, ENDFAC )


C...........   Read temporal XREF file; Construct mapping of source number 
C...........   to effective composite temporal profile number.


C...........   Read in emissions indexing data IFIP, IRCLAS, ILINK.

        WRITE( *,92000 ) ' ', 'Reading source list from MOBILE file.'

        IF ( .NOT. READ3( ANAME, 'FIP', ALLAYS3, 0,0, IFIP ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading FIP from file' // ANAME, 2)
        END IF

        IF ( .NOT. READ3( ANAME, 'IRCLAS', ALLAYS3, 0,0, IRCLAS ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading IRCLAS from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'ILINK', ALLAYS3, 0,0, ILINK ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading ILINK from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'TZONES', ALLAYS3, 0,0, TZONES ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading TZONES from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'TPFLAG', ALLAYS3,0,0, TFLAG ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading TFLAG from file' // ANAME, 2 )
        END IF

        IF ( .NOT. READ3( ANAME, 'VMT', ALLAYS3, 0, 0, VMT ) ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Error reading "VMT" from file "' // 
     &                   ANAME( 1:TRIMLEN( ANAME ) ) // '".', 2 )
        END IF


C...........   Read cross references and map to sources  (RDVMIX also
C...........   multiplies vehicle mix fractions by the VMT totals):

        CALL RDMTREF( XDEV, IFIP, IRCLAS, ILINK,
     &                      NMON, NWEK, NDIU, NEND,
     &                      MONREF, WEKREF, DIUREF, ENDREF,
     &                      MDEX, WDEX, DDEX )

        CALL RDVMIX( VDEV, JYEAR, IFIP, IRCLAS, ILINK, VMT, MXMAT )


C...........   Read and process emissions factors cross reference file MPLIST

        NMPLIST = RDMPLIST( FDEV, MPFIP, MPRDT, MPLNK,
     &                      MPREF, REFCNT, REFLIST )
        IF ( NMPLIST .LE. 0 ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0,0, 'Error reading MPLIST', 2 )
        END IF

C...........   Map emissions factor indexes onto sources:

        EFLAG = .FALSE.
        DO  355  S = 1, NMSRC

            FIP = IFIP  ( S )
            RDT = IRCLAS( S )
            LNK = ILINK ( S )

C............   Link-specific matches:

            F = FIND3( FIP, RDT, LNK, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  301  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
301             CONTINUE
                GO TO  355
            END IF

C............   FIP-RDT matches:

            F = FIND3( FIP, RDT, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  311  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
311             CONTINUE
                GO TO  355
            END IF

C............   state--RDT matches:

            FIP = 1000 * ( FIP / 1000 )
            F = FIND3( FIP, RDT, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  322  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
322             CONTINUE
                GO TO  355
            END IF

C............   RDT-only matches:

            F = FIND3( 0, RDT, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  333  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
333             CONTINUE
                GO TO  355
            END IF

C............   fallback

            FIP = 1000 * ( FIP / 1000 )
            F = FIND3( 0, 0, 0, NMPLIST, MPFIP, MPRDT, MPLNK )
            IF ( F .GT. 0 ) THEN
                DO  344  T = 0, 23
                    EFDEX( T,S ) = FIND1( MPREF( T,F ), 
     &                                    REFCNT, REFLIST )
344             CONTINUE
                GO TO  355
            END IF

C...........   If you get to here:  error

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &      'No emission factor xref for FIP', FIP,
     &      'RCLASS', RDT,
     &      'Link', LNK
            CALL M3MESG( MESG )

355     CONTINUE

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'TMPMOBIL', 0, 0, 
     &                   'Bad emission factor matchup', 2 )
        END IF                                               


C...........   Read relevant emissions factors from MEFS

        IF ( .NOT. DESC3( FNAME ) ) THEN
            MESG = 'Could not get description for ' // FNAME
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
        END IF


C.......   Get starting point, increment, number of items for
C.......   emissions factor temperature tables:

        FDATE  = SDATE3D
        NTEMP  = STR2INT( FDESC3D( 1 ) )
        TEMP0  = CTOK  +  SNGL( FTOC * ( XORIG3D - 32.0D0 ) )
        DTEMP  = SNGL( FTOC * XCELL3D )
        TEMP1  = TEMP0 +  DTEMP * FLOAT( NTEMP - 1 )
        DDTEMP = 1.0 / DTEMP


C.......   Read emissions factor tables:

        DO  366  T = 1, REFCNT

            IF ( .NOT. READ3( FNAME, 'EXH_CO', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        EXH_CO( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable EXH_CO from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'EXH_NOX', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        EXH_NOX( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable EXH_NOX from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'EXH_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        EXH_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable EXH_VOC from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'EVP_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        EVP_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable EVP_VOC from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'RNL_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        RNL_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable RNL_VOC from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. READ3( FNAME, 'RST_VOC', ALLAYS3, 
     &                        FDATE, SEC2TIME( REFLIST( T ) ),
     &                        RST_VOC( 1,1,T ) ) ) THEN
                WRITE( MESG,94010 ) 
     &          'Could not read record', REFLIST( T ),
     &          'of variable RST_VOC from ' // FNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF

366     CONTINUE



C.......   Read ungridding matrix:

        IF ( .NOT. READ3( UNAME, ALLVAR3, ALLAYS3, 0, 0, NU ) ) THEN
            MESG = 'Could not read ungridding matrix from ' // UNAME
            CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
        END IF


C...........   Construct TMAT -- array of effective composite 
C...........   profile coefficients

C...........   For each time step, write out the area source emissions values:

        WRITE( *,92000 ) ' ', 'Writing out emissions values...'

        LDATE = -1   !  initialization forces TMAT calculation on first iter.

        DO  899  T = 1, NSTEPS

            IF ( .NOT. READ3( MNAME, 'TA', 1, MDATE, MTIME, TA ) ) THEN
                MESG = 'Could not read TA from ' // MNAME
                CALL M3EXIT( 'TMPMOBIL', 0, 0, MESG, 2 )
            END IF
C...............   If this is a new month, or transitions between
C...............   weekday and weekend, calculate a new TMAT


            IF ( LDATE .NE. JDATE ) THEN

                LDATE = JDATE

                CALL DAYMON( JDATE, MONTH, DAY )  !  get month
                DAY = WKDAY( JDATE )              !  get day-of-week

                IF ( NEND .EQ. 0 ) THEN  !  no weekend packet

                    WRITE( *,92011 ) 'Processing day', DAY
                    CALL MKVMAT( NMSRC,  JDATE, MONTH,  DAY, 
     &                           TZONE,  TZONES,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           DIUFAC, DDEX, MXMAT, VMT, TMAT )

                ELSE IF ( DAY .GE. 6 ) THEN     !  day is Saturday or Sunday
                                                !  in SAI's brain-damaged
                                                !  day-of-week conventions.

                    WRITE( *,92011 ) 'Processing weekend day', DAY
                    CALL MKVMAT( NMSRC,  JDATE, MONTH, DAY, 
     &                           TZONE,  TZONES,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           ENDFAC, EDEX, MXMAT, VMT, TMAT )

                ELSE

                    WRITE( *,92011 ) 'Processing midweek day', DAY
                    CALL MKVMAT( NMSRC,  JDATE, MONTH, DAY, 
     &                           TZONE,  TZONES,
     &                           MONFAC, MDEX,
     &                           WEKFAC, WDEX,
     &                           DIUFAC, DDEX, MXMAT, VMT, TMAT )

                END IF

            END IF      ! if jdate not ldate

            H  = 1 + MOD( JTIME / 10000 , 24 )
            KU = 0

            DO  888  S = 1, NMSRC

C...............   Ungrid temperature:

                TASRC = 0.0
                DO  811  J = 1, NU( S )
                    KU = KU + 1
                    TASRC = TASRC + TA( IU( KU ) ) * CU( KU )
811             CONTINUE

                IF ( TASRC .LT. TEMP0 ) THEN	!  trap TASRC against TEMP0

                    WRITE( MESG,94030 )
     &                  'Temperature TEMP0', 
     &                  'at source', S,
     &                  'below MEFS table min', TEMP0
                    CALL M3MESG( MESG )
                    PP = 1.0
                    QQ = 0.0
                    II = 1

                ELSE IF ( TASRC .GT. TEMP1 ) THEN		!  and against TEMP1

                    WRITE( MESG,94030 )
     &                  'Temperature TEMP1', 
     &                  'at source', S,
     &                  'above MEFS table max', TEMP1
                    CALL M3MESG( MESG )
                    PP = 0.0
                    QQ = 1.0
                    II = NTEMP - 1

                ELSE			!  get table entry, interp coeffs

                    SCR = DDTEMP * ( TASRC - TEMP0 )
                    II  = 1 + INT( SCR )
                    PP  = AMOD( SCR, 1.0 )
                    QQ  = 1.0 - PP

                END IF		!  tasrc trapped; coeffs & index found.

                HCORR = MOD( H + TZONE + 25 - TZONES( S ), 24 )
                M     = EFDEX( HCORR, S )

                DO  818  V = 1, NVTYPE

                    TOUT( S,V,1 ) = TMAT( S,V,1,H ) *
     &                              ( PP * EXH_CO( II  ,V,M ) +
     &                                QQ * EXH_CO( II+1,V,M ) )

                    TOUT( S,V,2 ) = TMAT( S,V,2,H ) *
     &                              ( PP * EXH_NOX( II  ,V,M ) +
     &                                QQ * EXH_NOX( II+1,V,M ) )

                    TOUT( S,V,3 ) = TMAT( S,V,3,H ) *
     &                              ( PP * EXH_VOC( II  ,V,M ) +
     &                                QQ * EXH_VOC( II+1,V,M ) )

                    TOUT( S,V,4 ) = TMAT( S,V,4,H ) *
     &                              ( PP * EVP_VOC( II  ,V,M ) +
     &                                QQ * EVP_VOC( II+1,V,M ) )

                    TOUT( S,V,5 ) = TMAT( S,V,5,H ) *
     &                              ( PP * RNL_VOC( II  ,V,M ) +
     &                                QQ * RNL_VOC( II+1,V,M ) )

                    TOUT( S,V,6 ) = TMAT( S,V,6,H ) *
     &                              ( PP * RST_VOC( II  ,V,M ) +
     &                                QQ * RST_VOC( II+1,V,M ) )

818             CONTINUE	!  end loop calculating emissions for S

888         CONTINUE		!  end loop on sources S

            IF ( .NOT.WRITE3( TNAME, 'ALL', JDATE, JTIME, TOUT ) ) THEN

                CALL M3EXIT( 'TMPMOBIL', JDATE, JTIME,
     &                       'Could not write timestep to "' // 
     &                       TNAME( 1 : TRIMLEN( TNAME ) ) // '".', 2 )

            END IF

            CALL NEXTIME( JDATE, JTIME, TSTEP )
            CALL NEXTIME( MDATE, MTIME, TSTEP )

899     CONTINUE

999   CONTINUE          !  exit program

      CALL M3EXIT( 'TMPMOBIL',0,0, 'Normal completion of TMPMOBIL', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( 5X, A, I10 )

92011   FORMAT( /5X, A, I10 )

92020   FORMAT( 5X, A, I9 )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( A16 )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I5, :, 2X ) )

94020   FORMAT( A, I5, 1X, A, I7, 1X, 10( A, I7, :, 1X ) )

94030   FORMAT( 10 ( A, :, 1PG12.5, :, 2X , A, :, I5, :, 2X ) )

        END

