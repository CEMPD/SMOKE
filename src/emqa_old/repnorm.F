        PROGRAM REPNORM

C***********************************************************************
C  program body starts at line 45
C
C  DESCRIPTION:
C       Produces state and county totals for biogenic normalized 
C       emission categories.
C
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for surrogate coeffs.
C       EPS2 Formatted Surrogate file (will not work with EMS95 surgs)
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       CRLF, DSCGRID,
C       PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL,
C       TRIMLEN, ENVINT
C
C  REVISION  HISTORY:
C       Prototype  4/97 by JMV.
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
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
C****************************************************************************

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'BIODIMS3.EXT'  !  biogenic-source dimensioning parameters 
        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER*2     CRLF
        LOGICAL         DSCGRID
        LOGICAL         ENVYN
        INTEGER         ENVINT
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         STR2INT
        REAL            STR2REAL
        INTEGER         TRIMLEN

        EXTERNAL        CRLF, DSCGRID,
     &                  PROMPTFFILE, PROMPTMFILE, STR2INT, STR2REAL,
     &                  TRIMLEN, ENVINT

C
C...........   PARAMETERS:

        INTEGER         NRMAX                   !  max # of surrogate records
        CHARACTER*5     BLANK5
        CHARACTER*16    CLAI

        PARAMETER     ( NRMAX = 275 000 ,
     &                  BLANK5= ' '     ,       
     &                  CLAI  = 'LAI            ' )

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

        INTEGER         RDEV    !  output report file
        INTEGER         SDEV    !  for surrogate coeff file
        INTEGER         LDEV    !  log-device
        CHARACTER*16    ANAME   !  logical name for emission source gridded file

C.......   Name tables for file ANAME

        CHARACTER*16    NORMV( BSPCS-1, BTYPES ) ! names for VOC vbles
        CHARACTER*16    NORMN( LUSES )          !  names for  NO-emission vbles
 
C.......   Gridded normalized emissions description input from file BGRD

        REAL NVOC( NCOLS, NROWS, BSPCS-1, BTYPES ) !  normalized VOC emissions
        REAL AVLAI( NCOLS, NROWS )              !  average LAI
        REAL NORNO( NCOLS, NROWS, LUSES )       !  normalized NO emissions
        REAL CYVOC( NBFIP , BSPCS-1, BTYPES )   !  county VOC landuse totals
        REAL CYLAI( NBFIP )                     !  county LAI totals
        REAL CYNO( NBFIP, LUSES )               !  county NO totals
        REAL CTOT( NCOLS, NROWS )               !  normalized area fractions
        REAL STLAI( NASID )                     !  state LAI totals
        REAL STNO( NASID, LUSES )               !  state NO totals
        REAL STVOC( NASID, BSPCS-1, BTYPES )    !  state VOC totals
 
C...........   Other local variables

        REAL            DDX, DDY                !  inverse of cell widths
        REAL            SRGT( NASRG )           !  surrogates
        REAL            X1, X2, Y1, Y2, XC, YC  !  surrogate file header
        REAL            XOFF, YOFF              !  surrogate offsets from grid
        REAL            XX, YY                  !  cell number or l-l position
        REAL            SNORM, SLAI, TTMP       !  tmp conversion factors
        INTEGER         C, S, F, I, J, K, L, M, N, T  !  counters, subscripts
        INTEGER         COL,  ROW
        INTEGER         FIP, SID, CID
        INTEGER         HCOL, HROW              !  surrogate hdr cols and rows
        INTEGER         IOS                     !  I/O status
        INTEGER         IREC                    !  input line (record) number
        INTEGER         MFIP, LSTA, STA 
        INTEGER         LO, HI
        INTEGER         NSRG                 ! Count of srgts in EPS2 srgts
        INTEGER         B
        INTEGER         NSTA 
        LOGICAL         EFLAG   !  input error flat
        LOGICAL         ST_REPORT  ! Output state-based emissions totals

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
C   begin body of program REPNORM

        LDEV = INIT3()
C
C   Read in normalized emissions and surrogate factor(s)
C
        ANAME = PROMPTMFILE(
     &  'Enter logical name for Biogenic Gridded Emissions file',
     &    FSREAD3, 'BGRD', 'REPNORM' )


        SDEV = PROMPTFFILE(
     &  'Enter logical name for GRIDDING SURROGATE COEFF file',
     &           .TRUE., .TRUE., 'AGPRO', 'REPNORM' )
C
C   Open output file
C
        RDEV  = PROMPTFFILE(
     &             'Enter logical name for TOTALS REPORT or "NONE"',
     &             .FALSE., .TRUE., 'NORMREP' , 'REPNORM' )

C.......   Get environment variable settings

        ST_REPORT = ENVYN( 'ST_REPORT', 'Per-state report writing flag',
     &                      .TRUE., IOS )

        JSRG = ENVINT( 'JSRG', 
     &             'Record number for area surrogate',
     &             3 , IOS)

        IF ( IOS .GT. 0 ) THEN
           MESG = ' Bad integer for record number'
           CALL M3EXIT( 'REPNORM', 0, 0, MESG, 2 ) 
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
            CALL M3EXIT( 'REPNORM', 0, 0, MESG, 2 )

        END IF                  ! if dscgrid() failed

C.......   Build name table for variables in normalized emissions file"

        DO  23  K = 1, BTYPES
        DO  22  B = 1, BSPCS - 1

            NORMV( B, K ) = BIOLTYPE( K ) // BIOSPC( B )

22      CONTINUE
23      CONTINUE

        DO  24  L = 1, LUSES

            NORMN( L ) = BIOLUSE( L )( 1:TRIMLEN( BIOLUSE( L )))//'NO'

24      CONTINUE


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
                    CALL M3EXIT( 'REPNORM', 0, 0,
     &                   'XCELL is zero in SURROGATE file header', 2 )

                ELSEIF( YC .EQ. 0 ) THEN
                    CALL M3EXIT( 'REPNORM', 0, 0,
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
                    CALL M3EXIT( 'REPNORM', 0, 0,
     &                           'Bad XCELL in SURROGATE file', 2 )
                ELSEIF ( FLTERR( YC, SNGL( YCELL3D ) ) ) THEN
                    CALL M3EXIT( 'REPNORM', 0, 0,
     &                           'Bad YCELL in SURROGATE file', 2 )
                ELSEIF ( FLTERR( XOFF, FLOAT( NINT( XOFF ) ) ) ) THEN
                    CALL M3EXIT( 'REPNORM', 0, 0,
     &                     'Bad X alignment in SURROGATE file', 2 )
                ELSEIF ( FLTERR( YOFF, FLOAT( NINT( YOFF ) ) ) ) THEN
                    CALL M3EXIT( 'REPNORM', 0, 0,
     &                     'Bad Y alignment in SURROGATE file', 2 )
                ELSEIF ( FLTERR( REAL( HCOL ), REAL( NCOLS ) ) ) THEN
                    CALL M3EXIT( 'REPNORM', 0, 0,
     &                     'Bad ending X-coord in SURROGATE file', 2 )
                ELSEIF ( FLTERR( REAL( HROW ), REAL( NROWS ) ) ) THEN
                    CALL M3EXIT( 'REPNORM', 0, 0,
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
                    CALL M3EXIT( 'REPNORM', 0, 0, MESG, 2 )

                ELSEIF ( NSRG .GT. NASRC ) THEN
                    WRITE( MESG, 94010 )
     &               'Number of surrogates in file =', NSRG,
     &               'but maximum (NASRG) =', NASRG
                    CALL M3EXIT( 'REPNORM', 0, 0, MESG, 2 )

                END IF

            ELSE  ! Grid type is LAT-LON

                MESG= 'Using LAT-LON grid, with EPS2.0 surrogates file.'
                CALL M3WARN( 'REPNORM', 0, 0, MESG )

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
              CALL M3EXIT( 'REPNORM', 0, 0, MESG, 2 )
 
            ENDIF

64          CONTINUE

            GO TO  66  !  to head of SDEV read loop

77      CONTINUE       !  end of the SDEV-read loop

        CALL M3MSG2( 'NOTE: File read in as ' // FFORMAT // ' format.' )
C
C Loop thru to add up area fraction totals for each cell.
C Note a cell may be in more than one county
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


C.......   Loops reading the various categories of normalized emissions:

        I = 0

        DO  34  K = 1, BTYPES

         DO  33  J = 1, BSPCS - 1

            IF ( .NOT. READ3( ANAME, NORMV( J, K ), 1, 0, 0,
     &                      NVOC( 1,1,J,K ) ) ) THEN
                 MESG = 'Could not read "' //
     &                NORMV( J, K )( 1 : TRIMLEN( NORMV( J, K ) ) ) //
     &                '" from file "' //
     &                ANAME( 1:TRIMLEN( ANAME ) ) // '"'
                CALL M3EXIT( 'REPNORM', 0, 0, MESG, 2 )
            END IF

33       CONTINUE        !  end loop on VOC's

34      CONTINUE        !  end loop on VOC-emitting land types

        IF ( .NOT. READ3( ANAME, 'AVLAI', 1, 0, 0, AVLAI ) ) THEN
             MESG = 'Could not read AVLAI from file "' //
     &              ANAME( 1:TRIMLEN( ANAME ) ) // '"'
            CALL M3EXIT( 'REPNORM', 0, 0, MESG, 2 )
        END IF

        I = 0

        DO  35  J = 1, LUSES

            I = I + 1
            IF ( .NOT. READ3( ANAME, NORMN( I ), 1, 0, 0,
     &                      NORNO( 1,1,J ) ) ) THEN
                 MESG = 'Could not read "' //
     &                  NORMN( I )( 1 : TRIMLEN( NORMN( I ) ) ) //
     &                  '" from file "' //
     &                  ANAME( 1:TRIMLEN( ANAME ) ) // '"'
                CALL M3EXIT( 'REPNORM', 0, 0, MESG, 2 )
            END IF


35      CONTINUE        !  end loop reading normalized NO's

        DO I = 1, NBFIP

          CYLAI( I ) = 0.0
 
          DO J = 1, BSPCS-1
           DO K = 1, BTYPES
             CYVOC( I, J, K ) = 0.0
           ENDDO
          ENDDO
     
          DO J = 1, LUSES
             CYNO( I, J ) = 0.0
          ENDDO

        ENDDO
              

        OLDFIP = -9
        M = 0


        DO 94 MREC = 1, N

          I = COLA( MREC )
          J = ROWA( MREC )
          
          SNORM = FRACA( MREC ) * CTOT ( I , J ) * GM2TON
          SLAI  = FRACA( MREC ) * CTOT ( I , J ) 
          IF (FIPSA( MREC ) .GT. OLDFIP ) THEN
            M = M + 1 
            OLDFIP = FIPSA( MREC )
            KFIP( M ) = FIPSA( MREC )
          ENDIF

          CYLAI( M ) = CYLAI( M ) + ( AVLAI ( I, J ) * SLAI )  

          DO  104  K = 1, BTYPES

           DO  103  B = 1, BSPCS - 1
             CYVOC( M, B, K ) = CYVOC( M, B, K ) + 
     &                        ( NVOC( I, J, B, K ) * SNORM ) 
103        CONTINUE

104       CONTINUE

          DO 112 K = 1, LUSES

             CYNO( M, K ) = CYNO( M, K ) + 
     &                      ( NORNO( I, J, K) * SNORM )

112       CONTINUE

94      CONTINUE

        DO I = 1, NASID

          STLAI( I ) = 0.0

          DO J = 1, BSPCS-1

           DO K = 1, BTYPES
             STVOC( I, J, K ) = 0.0
           ENDDO

          ENDDO

          DO J = 1, LUSES
             STNO( I, J ) = 0.0
          ENDDO

        ENDDO

        MFIP = M
        K = 0
        LSTA = -9

        DO 120 J= 1, MFIP

          STA = KFIP( J ) / 1000
          
          IF ( STA .GT. LSTA ) THEN
       
             K = K + 1 
             STATES( K ) = STA
             LSTA = STA
      
          ENDIF
 
          STLAI( K ) = STLAI( K ) +CYLAI ( J )           

          DO  124  N = 1, BTYPES

           DO  123  B = 1, BSPCS - 1

             STVOC( K, B, N ) = STVOC( K, B, N ) + CYVOC( J, B, N ) 

123        CONTINUE

124       CONTINUE

          DO 132 N = 1, LUSES

             STNO( K, N ) = STNO( K, N ) + CYNO( J, N ) 

132       CONTINUE

120     CONTINUE

        NSTA = K 
C.............  Write report if report flag is switched on

            IF( ST_REPORT ) THEN 

              MESG = 'Biogenic normalized emissions ' //
     &                '(tons/ hour) within grid ' //
     &                GRDNM
              DO M = 1, BTYPES

               CALL STATELST( RDEV  , MESG , NASID , NSTA , BSPCS-1,
     &              'F8.2', NORMV( 1, M ) , STATES, STVOC( 1, 1, M)) 

              ENDDO


              CALL STATELST( RDEV  , MESG , NASID , NSTA , LUSES,
     &              'F8.2', NORMN , STATES, STNO )

              MESG = 'Biogenic leaf area index ' //
     &                'within grid ' //
     &                GRDNM

              CALL STATELST( RDEV  , MESG , NASID , NSTA , 1 ,
     &              'F8.2', CLAI , STATES, STLAI )
 
            ENDIF

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )
93010   FORMAT( I5, 2F10.0, 10X, 100( F10.8 ) )


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )

        STOP
        END
