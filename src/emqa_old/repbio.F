
        PROGRAM REPBIO

C***********************************************************************
C  program body starts at line 174
C
C  DESCRIPTION:
C       Produces state and county hourly totals for hourly biogenic  
C       emission species.
C
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for surrogate coeffs.
C       EPS2 Formatted Surrogate file (will not work with EMS95 surgs)
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       CRLF, DSCGRID, INDEX1,
C       PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL,
C       TRIMLEN, ENVINT
C
C  REVISION  HISTORY:
C       Prototype  4/97 by JMV.
C
C***************************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
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
C****************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'BIODIMS3.EXT'  !  biogenic-source dimensioning parameters 
        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal dimensioning parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        
C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        LOGICAL         DSCGRID
        INTEGER         ENVINT
        LOGICAL         ENVYN
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN
        INTEGER         WKDAY           !  day of week (1...7)

        EXTERNAL        CRLF, DSCGRID, ENVINT, ENVYN, INDEX1, MMDDYY,
     &                  PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL,
     &                  TRIMLEN, WKDAY
C
C...........   PARAMETERS:

        INTEGER         NRMAX                   !  max # of surrogate records
        CHARACTER*5     BLANK5

        PARAMETER     ( NRMAX = 275 000 ,
     &                  BLANK5= ' '     )       

        REAL GM2TON
        PARAMETER ( GM2TON = 1.102311E-06 )

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER     OLDFIP, MREC
        INTEGER     KFIP ( NBFIP )
        INTEGER     JSRG 
        INTEGER     SSCT( NASRG )
  
C...........   Surrogate-cell::FIPS tables (unsorted; sorted and processed)

        INTEGER     INDXA( NRMAX )        !  subscripts for sorti1()
        INTEGER     COLA( NRMAX )         !  column cell for source
        INTEGER     ROWA( NRMAX )         !  row cell for source
        INTEGER     FIPSA( NRMAX )        !  FIP code,           before sorting
        REAL        FRACA( NRMAX )        !  surrogate fraction, before sorting

        INTEGER     STATES( NASID )       !  state ids 
  
C...........   Logical file names and unit numbers
   
        INTEGER         RDEV    !  output hourly report file
        INTEGER         SDEV    !  for surrogate coeff file
        INTEGER         TDEV    !  output daily report file  
        INTEGER         LDEV    !  log-device
        CHARACTER*16    ANAME   !  logical name for emission source gridded file
        CHARACTER*25    CTIME                   ! date and hour

C.......   Name tables for file ANAME

        REAL            EMISG ( NGRID, MSPCS )  !  gridded emissions
        REAL            CTOT( NCOLS, NROWS )    !  normalized grid fraction  

        REAL            CYEM( NBFIP, MSPCS )    !  summed county emissions
        REAL            STEM( NASID, MSPCS )    !  hourly state emissions 
        REAL            STDAY( NASID, MSPCS )   !  daily totals for states

C...........   Other local variables

        REAL            DDX, DDY                !  inverse of cell widths
        REAL            SRGT( NASRG )           !  surrogate fractions
        REAL            X1, X2, Y1, Y2, XC, YC  !  surrogate file header
        REAL            XOFF, YOFF              !  surrogate offsets from grid
        REAL            XX, YY                  !  cell number or l-l position
        REAL            SNORM, TTMP             !  conversion variables
        REAL            RWT( MSPCS )            !  molcular wt of BEIS2 species 
        REAL            FAC     !  temporary conversion factor

        INTEGER         BTIME              ! beginning time of day in run
        INTEGER         LDATE, LTIME       
        INTEGER         JDATE, JTIME, TSTEP     !  time variables
        INTEGER         NSTEPS
        INTEGER         S, I, J, K, M, N, T, V ! counters, subscripts
        INTEGER         COL,  ROW                 
        INTEGER         FIP, SID, CID
        INTEGER         HCOL, HROW              !  surrogate hdr cols and rows
        INTEGER         IOS                     !  I/O status
        INTEGER         IREC                    !  input line (record) number
        INTEGER         MFIP, LSTA, STA 
        INTEGER         NSRG                 ! Count of srgts in EPS2 srgts
        INTEGER         B, IJ
        INTEGER         NSTA
        INTEGER         KS
        INTEGER         RTIME                ! time for reporting state totals

        LOGICAL         EFLAG      ! Input error flat
        LOGICAL         ST_REPORT  ! Output state-based emissions totals
        LOGICAL         BSPCFLAG( NMPOL )  ! Flag for BEIS2 species

        CHARACTER*5     FFORMAT !  temporary indicator for input formats
        CHARACTER*16    SCRBUF
        CHARACTER*256   LINE
        CHARACTER*256   MESG

C...........   STATEMENT FUNCTIONS:
C.......   floating point "unequal" -- true iff
C.......   | P - Q | > 1e-5 * sqrt( p*p + q*q + 1e-5 )

        REAL            P, Q
        LOGICAL         FLTERR
        FLTERR( P, Q ) =
     &  ( (P - Q)**2  .GT.  1.0E-10*( P*P + Q*Q + 1.0E-5 ) )

***********************************************************************
C   begin body of program REPBIO

        LDEV = INIT3()

        CALL INITEM( LDEV )
C
C  Read in hourly biogenic emission file
C
        ANAME = PROMPTMFILE(
     &  'Enter logical name for Biogenic Gridded Emissions file',
     &    FSREAD3, 'BGTS', 'REPBIO' )
C
C  Read in gridding surrogate file and open output report file
C
        SDEV = PROMPTFFILE(
     &  'Enter logical name for GRIDDING SURROGATE COEFF file',
     &           .TRUE., .TRUE., 'AGPRO', 'REPBIO' )

        RDEV  = PROMPTFFILE(
     &             'Enter logical name for HOURLY REPORT or "NONE"',
     &             .FALSE., .TRUE., 'BHRREP' , 'REPBIO' )

        TDEV  = PROMPTFFILE(
     &             'Enter logical name for DAILY REPORT or "NONE"',
     &             .FALSE., .TRUE., 'BDAYREP' , 'REPBIO' )

C
C  Write title to report file
C
        MESG = 'Biogenic hourly emissions ' //
     &           '(tons/hr) within grid ' //
     &             GRDNM

        IF( RDEV .GT. 0 ) WRITE( RDEV,93000 ) MESG

C.......   Get environment variable settings

        ST_REPORT = ENVYN( 'ST_REPORT', 'Per-state report writing flag',
     &                      .TRUE., IOS )

        RTIME = ENVINT( 'ST_REPTIME', 'Per-State summation reporting '
     *                  // 'time', 230000, IOS )

        JSRG = ENVINT( 'JSRG', 
     &             'Record number for area surrogate',
     &             3 , IOS)

        IF ( IOS .GT. 0 ) THEN
           MESG = ' Bad integer for record number'
           CALL M3EXIT( 'REPBIO', 0, 0, MESG, 2 ) 
        END IF


C.......   Read the grid parameters from GRIDDESC for use in checking
C.......   EPS2.0-formatted surrogate files

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF, GDTYP3D,
     &              P_ALP3D, P_BET3D,P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D,
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'REPBIO', 0, 0, MESG, 2 )

        END IF                  ! if dscgrid() failed


C...........   Read and process the surrogate coefficient file:
C...............   Must be sorted by FIP (state ID, county ID), surrogate ID,
C...............   col, and row.
C...............   All actually-occurring area source FIP codes must be present.

        CALL M3MSG2( 'Reading GRIDDING SURROGATE COEFF file...' )

        FFORMAT = 'EPS2'

C........  Can read EMS-95 _or_ EPS2.0 formatted gridding surrogate
C........  files.  The distinguishing factor is that EMS95 files do not
C........  have a blank space in column 6

C...........   Initialize parameters for reading surrogates data
        IREC    =  0
        N       =  0
        EFLAG   = .FALSE.
C........  Determine input format of gridding surrogates file
C........  For EPS2 format, treat 1st line as header line, and QA file
C........  For Lat-lon grids, header will be irrelevant, but can't check either

        READ( SDEV, 93000, END=77, IOSTAT=IOS ) LINE

        IF( TRIMLEN( LINE ) .LT. 60 ) THEN  ! EMS-95 format

            FFORMAT = 'EMS95'
            NSRG    = 1
            REWIND( SDEV )

        ELSE                                ! EPS2   format

            IREC = IREC + 1

            IF( GDTYP3D .NE. LATGRD3 ) THEN

                X1  = STR2REAL( LINE(  1:10 ) )
                Y1  = STR2REAL( LINE( 11:20 ) )
                X2  = STR2REAL( LINE( 21:30 ) )
                Y2  = STR2REAL( LINE( 31:40 ) )
                NSRG= STR2INT ( LINE( 41:50 ) )
                XC  = STR2REAL( LINE( 51:60 ) )
                YC  = STR2REAL( LINE( 61:70 ) )

                IF( XC .EQ. 0 ) THEN
                    CALL M3EXIT( 'REPBIO', 0, 0,
     &                   'XCELL is zero in SURROGATE file header', 2 )

                ELSEIF( YC .EQ. 0 ) THEN
                    CALL M3EXIT( 'REPBIO', 0, 0,
     &                   'YCELL is zero in SURROGATE file header', 2 )

                ELSE
                    DDX = 1.0 / XC
                    DDY = 1.0 / YC
                ENDIF

                XOFF= DDX * ( X1 - XORIG3D )
                YOFF= DDY * ( Y1 - YORIG3D )
                HCOL= NINT( DDX * ( X2 - X1 ) )
                HROW= NINT( DDY * ( Y2 - Y1 ) )
             IF ( FLTERR( XC, SNGL( XCELL3D ) ) ) THEN
                    CALL M3EXIT( 'REPBIO', 0, 0,
     &                           'Bad XCELL in SURROGATE file', 2 )
                ELSEIF ( FLTERR( YC, SNGL( YCELL3D ) ) ) THEN
                    CALL M3EXIT( 'REPBIO', 0, 0,
     &                           'Bad YCELL in SURROGATE file', 2 )
                ELSEIF ( FLTERR( XOFF, FLOAT( NINT( XOFF ) ) ) ) THEN
                    CALL M3EXIT( 'REPBIO', 0, 0,
     &                     'Bad X alignment in SURROGATE file', 2 )
                ELSEIF ( FLTERR( YOFF, FLOAT( NINT( YOFF ) ) ) ) THEN
                    CALL M3EXIT( 'REPBIO', 0, 0,
     &                     'Bad Y alignment in SURROGATE file', 2 )
                ELSEIF ( FLTERR( REAL( HCOL ), REAL( NCOLS ) ) ) THEN
                    CALL M3EXIT( 'REPBIO', 0, 0,
     &                     'Bad ending X-coord in SURROGATE file', 2 )
                ELSEIF ( FLTERR( REAL( HROW ), REAL( NROWS ) ) ) THEN
                    CALL M3EXIT( 'REPBIO', 0, 0,
     &                     'Bad ending Y-coord in SURROGATE file', 2 )
                ELSEIF ( NSRG .LE. 0 ) THEN
                    WRITE( MESG, 94010 )
     &               'NOTE: Use ZONE portion of EPS2.0 surrogates '//
     &               'header' // CRLF() // BLANK5 //
     &               '      to specify number of surrogates ' //
     &               'in the file.'
                    CALL M3MSG2( MESG )

                    WRITE( MESG, 94010 )
     &               'Number of surrogates in header was ', NSRG
                    CALL M3EXIT( 'REPBIO', 0, 0, MESG, 2 )

                ELSEIF ( NSRG .GT. NASRC ) THEN
                    WRITE( MESG, 94010 )
     &               'Number of surrogates in file =', NSRG,
     &               'but maximum (NASRG) =', NASRG
                    CALL M3EXIT( 'REPBIO', 0, 0, MESG, 2 )

                END IF

            ELSE  ! Grid type is LAT-LON

                MESG= 'Using LAT-LON grid, with EPS2.0 surrogates file.'
                CALL M3WARN( 'REPBIO', 0, 0, MESG )

            END IF

        ENDIF

C...........   Start loop for actually reading surrogates data

66      CONTINUE                !  head of the SDEV-read loop

C...............  Read EMS-95 formatted file
            IF( FFORMAT .EQ. 'EMS95' ) THEN

                READ( SDEV, *, END=77, IOSTAT=IOS )
     &              SSCT( NSRG ),   !  surrogate ID
     &              SID, CID,       !  state, county FIP codes
     &              COL, ROW,       !  grid col and row
     &              SRGT( NSRG )    !  surrogate fraction

                IREC = IREC + 1

                FIP  = 1000 * SID + CID

C...............  Read EPS2 file.  This read is made somewhat more complex
C...............  because we're conforming to EMS95 file structure.  However,
C...............  since the EPS2.0 format does not specify exact column
C...............  numbers for each surrogate class (just must start on
C...............  column 36 and continue to column 186), and the number of
C...............  surrogates is not pre-defined, the EPS2.0 format is ill-
C...............  defined anyway.  Here, we add the constraint that the
C...............  zone field is used to define the number of surrogate
C...............  entries. Also, all file columns are as described in the
C...............  EPS2 documentation, AND the surrogates entries are 10
C...............  columns wide.  These constraints allow us to read in the
C...............  file more assuredly.  MRH.
            ELSEIF( FFORMAT .EQ. 'EPS2' ) THEN

                READ( SDEV, 93010, END=77, IOSTAT=IOS )
     &                FIP, XX, YY, ( SRGT( J ), J = 1, NSRG )

                IREC = IREC + 1

C.................  First try col/row format
                COL = INT( XX )
                ROW = INT( YY )

C.................  If no good, then try coordinates format
                IF( COL .LT. 1  .OR. COL .GT. NCOLS .OR.
     &              ROW .LT. 1  .OR. ROW .GT. NROWS      ) THEN

                    COL = 1 + NINT( DDX * ( XX - XORIG3D ) )
                    ROW = 1 + NINT( DDY * ( YY - YORIG3D ) )

                ENDIF


            END IF         ! end input format type
            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94010 )
     &              'I/O error ', IOS,
     &              'reading SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  66   ! to head of SDEV read loop

C.................  If columns and rows are still a problem, give message,
C.................  but don't error because may want to use bigger
C.................  surrogates file to do smaller domain
            ELSE IF ( COL .LT. 1  .OR. COL .GT. NCOLS ) THEN

                WRITE( MESG,94010 )
     &              'Column ', COL,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                GO TO  66   ! to head of SDEV read loop

            ELSE IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN

                WRITE( MESG,94010 )
     &              'Row ', ROW,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                GO TO  66   ! to head of SDEV read loop

            END IF              !  if IOS bad, or col, or row out of range

              
C...........   Record surrogate fraction for this cell.

            N = N + 1

            IF( N .LE. NRMAX ) THEN

               INDXA( N ) = N
               ROWA( N ) = ROW
               COLA( N ) = COL
               FIPSA( N ) = FIP
               FRACA( N ) = SRGT( JSRG )

            ELSE

              WRITE( MESG, 94010 )
     &        'Number of records in AGPRO exceeds ', NRMAX
              CALL M3EXIT( 'REPBIO', 0, 0, MESG, 2 )
 
            ENDIF

64          CONTINUE

            GO TO  66  !  to head of SDEV read loop

77      CONTINUE       !  end of the SDEV-read loop

        CALL M3MSG2( 'NOTE: File read in as ' // FFORMAT // ' format.' )
C
C Loop thru to add up area fraction totals for each cell.
C Note a cell may be in more than one county.  Compute inverse
C to use for normalization of fraction.
C
        DO 80 I = 1, NCOLS
         DO 82 J = 1, NROWS

           CTOT( I, J ) = 0.0
           TTMP=0.0

           DO 84 MREC = 1, N

             IF (I .EQ. COLA( MREC ) .AND. J .EQ. ROWA( MREC )) THEN
               TTMP = TTMP + FRACA ( MREC )
             ENDIF

84         CONTINUE 

           IF (TTMP .GT. 0.0)THEN
            CTOT( I, J ) = 1.0 / TTMP
           ELSE
            CTOT( I, J ) = 0.0
           ENDIF

82       CONTINUE

80      CONTINUE         

C
C  Read description of biogenic emissions file
C
        IF ( .NOT. DESC3( ANAME ) ) THEN
           MESG = 'Could not get description of file ' // ANAME
           CALL M3EXIT( 'REPBIO', 0, 0, MESG, 2 )
        END IF
C
C  True if a biogenic species
C
        DO  11  V = 1, NMPOL

           BSPCFLAG( V ) = ( INDEX1( EMNAM( V ),
     &                       NVARS3D, VNAME3D ) .GT. 0 )

11      CONTINUE

C  Initialize daily totals

        DO I = 1, NASID

           DO K = 1, MSPCS
                 STDAY( I, K ) = 0.0
           ENDDO

        ENDDO

C
C   Setup local time variables
C
        JDATE   = SDATE3D
        JTIME   = STIME3D
        TSTEP   = TSTEP3D

        IF ( TSTEP .EQ. 0 ) THEN
           NSTEPS  = 1
        ELSE
           NSTEPS  = MXREC3D
        END IF


        WRITE( *,92000 )
     &      ' ',
     &      'Calculating and writing emissions totals...',
     &      ' '
C
C  Begin time steps and initialized county emissions.
C
        BTIME = JTIME
        LDATE = 0

        DO  199  T = 1, NSTEPS
          DO I = 1, NBFIP

             DO J = 1, MSPCS
                CYEM( I, J ) = 0.0
             ENDDO
  
          ENDDO

          KS = 0


             DO  197  K = 1, NMPOL


C...................   Check to see if a biogenic species and store it 
C...................   and its molecular weight 
     
                    IF ( BSPCFLAG( K ) ) THEN

                      KS = KS + 1

                      IF ( KS .GT. MSPCS ) THEN

                       MESG = 'Number of biogenic species from file"' //
     &                 ANAME( 1:TRIMLEN( ANAME ) ) // '" exceeds MSPCS' 
                       CALL M3EXIT( 'REPBIO', 0, 0, MESG, 2 )

                      ENDIF

                      RWT( KS ) =  EMWT( K )

                      IF( .NOT. READ3( ANAME, EMSPC( KS ), ALLAYS3,
     &                    JDATE, JTIME, EMISG( 1, KS ) )) THEN

                       MESG = 'Could not read emissions from file"' //
     &                 ANAME( 1:TRIMLEN( ANAME ) ) // '"'
                       CALL M3EXIT( 'REPBIO', 0, 0, MESG, 2 )

                      ENDIF

                    ENDIF


197             CONTINUE

C 
C  Begin loop to compute county emissions
C 
                OLDFIP = -9
                M = 0

                DO 94 MREC = 1, N

                  I = COLA( MREC )
                  J = ROWA( MREC )
                  IJ = I +  NCOLS * ( J - 1 ) 
C
C  Conversion factor of grams to tons included in normalization
C  factor SNORM
C
                  SNORM = FRACA( MREC ) * CTOT ( I , J ) 

                  IF (FIPSA( MREC ) .GT. OLDFIP ) THEN

                      M = M + 1 
                      OLDFIP = FIPSA( MREC )
                      KFIP( M ) = FIPSA( MREC )

                  ENDIF
C
C  Conversion plus molecular weight conversion (moles)
C 
                  DO  103  B = 1, MSPCS 
                    CYEM( M, B ) = CYEM( M, B ) + ( EMISG( IJ, B)
     &                    * SNORM ) 
103               CONTINUE

94               CONTINUE

C
C  Initialized state emission totals and other variables
C
                 DO I = 1, NASID

                  DO K = 1, MSPCS
                   STEM( I, K ) = 0.0
                  ENDDO

                 ENDDO

                 MFIP = M
                 K = 0
                 LSTA = -9

C 
C  Begin loop to compute state totals
C
                 DO 120 J= 1, MFIP
                    STA = KFIP( J ) / 1000
          
                    IF ( STA .GT. LSTA ) THEN
                       K = K + 1 
                       STATES( K ) = STA
                       LSTA = STA
                    ENDIF
 
                    DO  123  B = 1, MSPCS 

                     STEM( K, B ) = STEM( K, B ) + CYEM( J, B ) 

123                 CONTINUE


120               CONTINUE

C   Total number of states in the domain
C
                  NSTA = K 

C.............  Write report if report flag is switched on

            LDATE = JDATE
            LTIME = JTIME

            CALL NEXTIME( JDATE, JTIME, TSTEP )

            IF( ST_REPORT ) THEN  

               DO 193 S = 1, MSPCS

                   FAC = RWT( S ) * GM2TON 

                   DO 191 I = 1, NSTA

                       STEM( I,S ) = STEM( I,S ) * FAC
                       STDAY( I,S ) = STDAY( I,S ) + STEM( I,S ) 

191                CONTINUE

193            CONTINUE

               WRITE(CTIME,92010)LDATE, LTIME


               MESG = 'For date / hour ' //
     &                CTIME

               IF( RDEV .GT. 0 )
     &         CALL STATELST( RDEV  , MESG , NASID , NSTA ,
     &              MSPCS, 'F8.1', EMSPC , STATES, STEM )


C.............  Write totals report if previous time equals the report time

               IF ( LTIME .EQ. RTIME .OR. T .EQ. NSTEPS ) THEN

                  IF( TDEV .GT. 0 )
     &            WRITE( TDEV, 94010 ) DAYS( WKDAY( LDATE ) ) //
     &                               MMDDYY( LDATE ) // ' from',
     &                               BTIME, 'to', LTIME

                  MESG = 'Biogenic emissions (tons/specific day) ' //
     &                 'within grid ' // GRDNM

                  IF( TDEV .GT. 0 )
     &            CALL STATELST( TDEV  , MESG , NASID , NSTA ,
     &              MSPCS, 'F9.1', EMSPC , STATES, STDAY )

                  DO 298 I = 1, NASID

                    DO 297 S = 1, MSPCS
                        STDAY( I,S ) = 0.
297                 CONTINUE

298               CONTINUE

C.................  Set begining time for next day
                  BTIME = JTIME

               ENDIF 


            ENDIF


199     CONTINUE

C...........   Formatted file I/O formats............ 93xxx

92000   FORMAT( 5X, A )
92010   FORMAT( I10,' / ',I10 ) 
93000   FORMAT( A )
93010   FORMAT( I5, 2F10.0, 10X, 100( F10.8 ) )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )

        STOP
        END
