
        PROGRAM SRGGRID

C***********************************************************************
C  program body starts at line 138
C
C  DESCRIPTION:
C       Construct by-FIP and sum-total gridded files containing data from
C       EPS-style surrogates file.
C
C  PRECONDITIONS REQUIRED:
C       Sorted, cut-down input data for surrogate coeffs.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND1, FIND2, FIND3, GETYN, PROMPTFFILE, PROMPTMFILE, TRIMLEN
C
C  REVISION  HISTORY:
C       Prototype  9/96 by CJC.
C       NOTE: County totals processing not QA'd
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE)
C                Modeling System
C File Version  @(#)$Id$
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
C Pathname:     $Source$
C Last updated: $Date$ 
C
C***************************************************************************/

      IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'ARDIMS3.EXT'   !  area-source dimensioning parameters
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         FIND1, FIND2, FIND3
        INTEGER         GETMENU, GETNUM
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        LOGICAL         DSCGRID
        INTEGER         TRIMLEN
        
        EXTERNAL  DSCGRID, FIND1, FIND2, FIND3, GETMENU, GETNUM, GETYN,
     &            PROMPTFFILE, PROMPTMFILE, TRIMLEN


C...........   PARAMETERS:

        INTEGER         NRMAX                   !  max # of surrogate records
        PARAMETER     ( NRMAX = 275 000 )


C...........   LOCAL VARIABLES and their descriptions:
C...........   NOTE that ASC (area-source-category) ID's are 10-digit
C...........   FIPS/areas tables

        INTEGER     NSTALIST              !  # of FIPS
        INTEGER     STALIST( NAFIP )      !  STATE-codes list
        INTEGER     STAINDX( NAFIP )      !  1st loc of state in FIPLIST

        INTEGER     NFIPLIST                     !  # of FIPS
        INTEGER     FIPLIST(        NAFIP )      !  FIP-codes list
        REAL        FIPINDX(        NAFIP )      !  1st loc of FIP in FIPS()
        REAL        FIPSIZE( NASRG, NAFIP )      !  1 / sum_cells( srg( fip ) )

C...........   Surrogate-cell::FIPS tables (unsorted; sorted and processed)

        INTEGER     NSREC                 !  number of (unsorted) entries
        INTEGER     INDXA( NRMAX )        !  subscripts for sorti1()
        INTEGER     CELLA( NRMAX )        !  cell number,        before sorting
        INTEGER     FIPSA( NRMAX )        !  FIP code,           before sorting
        INTEGER     SIDSA( NRMAX )        !  surrogate ID,       before sorting
        REAL        FRACA( NRMAX )        !  surrogate fraction, before sorting

        INTEGER     CELL( NRMAX )        !  cell number,        before sorting
        INTEGER     FIPS( NRMAX )        !  FIP code,           before sorting
        INTEGER     SIDS( NRMAX )        !  surrogate ID,       before sorting
        REAL        FRAC( NRMAX )        !  surrogate fraction, before sorting

C...........   Output grids:

        REAL            GRID( NGRID, NASRG )

C...........   Other local variables

        REAL            SRG
        INTEGER         FIP, STA, SSC, SID, CID
        INTEGER         LFIP, LSTA
        INTEGER         COL,  ROW
        INTEGER         LO, HI
        INTEGER         JDATE, JTIME
        INTEGER         CHOICE, LOLIST, HILIST
        INTEGER         C, S, F, J, K, N  !  counters, subscripts
        INTEGER         SDEV    !  for surrogate coeff file
        INTEGER         TDEV    !  for surrogate total file
        INTEGER         LDEV    !  log-device
        INTEGER         IOS     !  I/O status
        INTEGER         IREC    !  input line (record) number
        LOGICAL         EFLAG   !  input error flat
        LOGICAL         TFLAG   !  iff totals file provided
        CHARACTER*16    GNAME   !  logical name for grid matrix     output file
        CHARACTER*16    SCRBUF
        CHARACTER*256   MESG

        CHARACTER*72    MENU( 3 )
     &  /
     &  'Compute gridded surrogates for ALL COUNTIES',	!  choice 1
     &  '... for all counties in a selected STATE',	!  choice 2
     &  '... for a selected COUNTY in a selected STATE'	!  choice 3
     &  /

C***********************************************************************
C   begin body of program SRGGRID

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &  ' ',
     &  'Program SRGGRID to take the EMS-95 SURROGATE COEFFICIENT,',
     &  'file and the SURROGATE TOTAL BY FIP file, and produce gridded',
     &  'SURROGATE TOTAL and SURROGATE BY-FIPS file.', 
     &  ' ',
     &  'NOTE: the EMS95 surrogate cross reference file must have',
     &  '      been sorted by FIP (state/county) and SURROGATE codes',
     &  '      prior to program execution.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  Input files must',
     &  'have been sorted as indicated, prior to program execution.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry. Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'SRGGRID', 0, 0, 'Ending program SRGGRID', 2 )
        END IF


C.......   Get file name; open actual-scc file

        SDEV = PROMPTFFILE( 
     &  'Enter logical name for SURROGATE COEFF file',
     &           .TRUE., .TRUE., 'AGPRO', 'SRGGRID' )

        TDEV = PROMPTFFILE( 
     &  'Enter logical name for SURROGATE TOTAL file or "NONE"',
     &           .TRUE., .TRUE., 'AGTOT', 'SRGGRID' )

        TFLAG = ( TDEV .NE. -2 )

C.......   Get the grid name and parameters
C.......   to put into file header:

        IF ( .NOT. DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,
     &              XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            SCRBUF = GRDNM
            MESG   = '"' // SCRBUF( 1:TRIMLEN( SCRBUF ) ) //
     &               '" not found in GRIDDESC file'
            CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )

        END IF                  ! if dscgrid() failed

        FTYPE3D = GRDDED3
        SDATE3D = 0
        STIME3D = 10000
        TSTEP3D = 10000
        NVARS3D = NASRG
        NLAYS3D = 1
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3
        FDESC3D( 1 ) = 'Area-Source gridding surrogates, by FIP'
        FDESC3D( 2 ) = 
     &  'First record is totals, SUM_fip( surg( fip )*tot( fip ) )'
        FDESC3D( 3 ) = 
     &  'Other records are FIP-specific gridded surrogate fractions'
        FDESC3D( 4 ) = 
     &  'with the FIPS in consecutive sorted order.'

        DO  11   K = 5, MXDESC3
            FDESC3D( K ) = ' '
11      CONTINUE

        VNAME3D( 1 ) = 'AGRI'
        UNITS3D( 1 ) = 'n/a'
        VDESC3D( 1 ) = 'AGRI gridding surrogate, by FIP'
        VTYPE3D( 1 ) = M3REAL

        VNAME3D( 2 ) = 'AIRPORT'
        UNITS3D( 2 ) = 'n/a'
        VDESC3D( 2 ) = 'AIRPORT gridding surrogate, by FIP'
        VTYPE3D( 2 ) = M3REAL

        VNAME3D( 3 ) = 'AREA'
        UNITS3D( 3 ) = 'n/a'
        VDESC3D( 3 ) = 'AREA gridding surrogate, by FIP'
        VTYPE3D( 3 ) = M3REAL

        VNAME3D( 4 ) = 'HOUSING'
        UNITS3D( 4 ) = 'n/a'
        VDESC3D( 4 ) = 'HOUSING gridding surrogate, by FIP'
        VTYPE3D( 4 ) = M3REAL

        VNAME3D( 5 ) = 'INVHOUSING'
        UNITS3D( 5 ) = 'n/a'
        VDESC3D( 5 ) = 'INVERSE HOUSING gridding surrogate, by FIP'
        VTYPE3D( 5 ) = M3REAL

        VNAME3D( 6 ) = 'INVPOP'
        UNITS3D( 6 ) = 'n/a'
        VDESC3D( 6 ) = 'INVERSE POPULATION gridding surrogate, by FIP'
        VTYPE3D( 6 ) = M3REAL

        VNAME3D( 7 ) = 'HWY'
        UNITS3D( 7 ) = 'n/a'
        VDESC3D( 7 ) = 'MAJOR HIGHWAYS gridding surrogate, by FIP'
        VTYPE3D( 7 ) = M3REAL

        VNAME3D( 8 ) = 'POP'
        UNITS3D( 8 ) = 'n/a'
        VDESC3D( 8 ) = 'POPULATION gridding surrogate, by FIP'
        VTYPE3D( 8 ) = M3REAL

        VNAME3D( 9 ) = 'PORTS'
        UNITS3D( 9 ) = 'n/a'
        VDESC3D( 9 ) = 'PORTS gridding surrogate, by FIP'
        VTYPE3D( 9 ) = M3REAL

        VNAME3D( 10 ) = 'RR'
        UNITS3D( 10 ) = 'n/a'
        VDESC3D( 10 ) = 'RAILROADS gridding surrogate, by FIP'
        VTYPE3D( 10 ) = M3REAL

        VNAME3D( 11 ) = 'WATER'
        UNITS3D( 11 ) = 'n/a'
        VDESC3D( 11 ) = 'WATER gridding surrogate, by FIP'
        VTYPE3D( 11 ) = M3REAL

        GNAME = PROMPTMFILE( 
     &      'Enter logical name for GRIDDING SURROGATES output file',
     &      FSUNKN3, 'ASURG', 'SRGGRID' )


C...........   Read and process the surrogate coefficient file:
C...............   Must be sorted by FIP (state ID, county ID), surrogate ID,
C...............   col, and row.  
C...............   All actually-occurring area source FIP codes must be present.

        WRITE( *,92000 ) 
     &      ' ', 'Reading GRIDDING SURROGATE COEFF file...', ' '

        IREC    =  1     
        N       =  0
        EFLAG   = .FALSE.

22      CONTINUE                !  head of the SDEV-read loop

            IREC = IREC + 1
            READ( SDEV, *, END=23, IOSTAT=IOS ) 
     &          SSC,            !  surrogate ID
     &          SID, CID,       !  state, county FIP codes
     &          COL, ROW,       !  grid col and row
     &          SRG             !  surrogate fraction

            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94010 ) 
     &              'I/O error ', IOS, 
     &              'reading SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  22

            ELSE IF ( COL .LT. 1  .OR. COL .GT. NCOLS ) THEN

                WRITE( MESG,94010 )
     &              'Column ', COL,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  22

            ELSE IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN

                WRITE( MESG,94010 )
     &              'Row ', ROW,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  22

            ELSE IF ( SSC .LT. 1  .OR. SSC .GT. NASRG ) THEN

                WRITE( MESG,94010 )
     &              'Surrogate code ', SSC,
     &              'out of range in SURG COEFF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                GO TO  22

            END IF              !  if IOS bad, or col,row, or ssc out of range


C...........   Record surrogate fraction for this cell.

            N = N + 1

            IF ( N .LE. NRMAX ) THEN

                INDXA( N ) = N
                CELLA( N ) = COL  +  NCOLS * ( ROW - 1 )
                SIDSA( N ) = SSC
                FIPSA( N ) = 1000 * SID + CID
                FRACA( N ) = SRG

            END IF                      ! if n in bounds

            GO TO  22                   !  to head of loop

23      CONTINUE                        !  end of the SDEV-read loop


        NSREC = N
        WRITE( MESG,94010 ) 
     &      'Max SURG COEF records - actual:', N,
     &      'dimensioned:', NRMAX
        CALL M3MESG( MESG )

        IF ( N .GT. NRMAX ) THEN
            CALL M3EXIT( 'SRGGRID', 0, 0, 
     &                   'Table overflow reading SURG COEFF file', 2 )
        END IF

        IF ( EFLAG ) THEN
            MESG = 'Error reading SURG COEF file.'
            CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )
        END IF


C...........   Sort and process the surrogates:

        WRITE( *,92000 ) 
     &      ' ', 'Processing GRIDDING SURROGATE COEFF file...', ' '

        CALL SORTI3( NSREC, INDXA, FIPSA, CELLA, SIDSA )

        F    = 0
        S    = 0
        LFIP = 0
        LSTA = 0

        DO  33  J = 1, NSREC            !  shuffle SURROGATE COEFFs

            K   = INDXA( J )
            C   = CELLA( K )
            FIP = FIPSA( K )
            SID = SIDSA( K )
            SRG = FRACA( K )

            CELL( J ) = C
            FIPS( J ) = FIP
            SIDS( J ) = SID
            FRAC( J ) = SRG

            IF ( FIP .NE. LFIP ) THEN

                F = F + 1

                IF ( F .LT. NAFIP ) THEN
                    FIPLIST( F ) = FIP
                    FIPINDX( F ) = J
                    STA  = FIP / 1000
                    IF ( STA .NE. LSTA ) THEN
                        S = S + 1
                        STALIST( S ) = STA
                        STAINDX( S ) = F
                        LSTA = STA
                    END IF
                END IF          !  if fip found, or not

                LFIP = FIP
                
            END IF              !  if new FIP encountered

33      CONTINUE                !  end loop shuffling SURROGATE COEFFs

        FIPINDX( F+1 ) = NSREC + 1
        STAINDX( S+1 ) = F     + 1
        NFIPLIST = F
        NSTALIST = S

        WRITE( MESG,94010 ) 
     &      'SURG COEF table - actual:', NFIPLIST,
     &      'dimensioned:', NAFIP
        CALL M3MESG( MESG )

        IF ( NFIPLIST .GT. NAFIP ) THEN
            CALL M3EXIT( 'SRGGRID', 0, 0, 
     &      'Table overflow processing SURG COEFF file', 2 )
        END IF


C...........  If surrogates file present, process SURROGATE-TOTAL-BY-FIP file:

        IF( TFLAG ) THEN

            DO  44  F = 1, NFIPLIST
            DO  43  S = 1, NASRG        
                FIPSIZE( S,F ) = 0.0
43          CONTINUE                !  end loop on SID's
44          CONTINUE                !  end loop on FIPS

            IREC  = 0
            EFLAG = .FALSE.

55          CONTINUE                !  head of loop reading TDEV

                READ( TDEV,*, END=66, IOSTAT=IOS ) FIP, SSC, SRG

                IF ( IOS .NE. 0 ) THEN

                    WRITE( MESG,94010 ) 
     &                  'I/O error ', IOS, 
     &                  'reading SURG TOTAL file at line', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  55

                ELSE IF ( SSC .LT. 1  .OR. SSC .GT. NASRG ) THEN

                    WRITE( MESG,94010 )
     &                  'Surrogate code ', SSC,
     &                  'out of range in SURG TOTAL file at line', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                    GO TO  55

                END IF              !  if IO error, or if SID out of range

                F = FIND1( FIP, NFIPLIST, FIPLIST )

                IF ( S .GE. 0 ) THEN
                    FIPSIZE( SSC,F ) = SRG
                ELSE
                    WRITE( MESG,94010 )
     &                  'FIP', FIP, 
     &                  'not found in SURG TOTAL file at line', IREC
                    CALL M3MESG( MESG )
                    EFLAG = .TRUE.
                END IF

                GO TO  55           !  to head of loop

66          CONTINUE                !  end of loop reading TDEV

            IF ( EFLAG ) THEN
                MESG = 'Error reading SURG TOTAL file.'
                CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )
            END IF

C...............   Process gridded surrogate-totals:

            DO  122  S = 1, NASRG
            DO  111  C = 1, NGRID
                GRID( C,S ) = 0.0
111         CONTINUE
122         CONTINUE

            LFIP = -1
            DO  133  K = 1, NSREC
                FIP = FIPS( K )
                IF ( FIP .NE. LFIP ) THEN
                    F   = FIND1( FIP, NFIPLIST, FIPLIST )
                    LFIP = FIP
                END IF
                C = CELL( K )
                S = SIDS( K )
                GRID( C,S ) = GRID( C,S )  +  FRAC( K ) * FIPSIZE( S,F )
133         CONTINUE

CXCC        IF ( .NOT. WRITE3( GNAME, 0, 0, 'ALL' ) ) THEN
CXCC            MESG = 'Could not write SURROGATE TOTALS to "' //
CXCC &                 GNAME( 1 : TRIMLEN( GNAME ) ) // '"'
CXCC            CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )
CXCC        END IF

        ENDIF   ! End if surrogates file present or not

C...........   Select sublist to compute:

        CHOICE = GETMENU( 3, 2, 'Enter type of county-selection',
     &                    MENU )

        IF ( CHOICE .EQ. 1 ) THEN

            LOLIST = 1
            HILIST = NFIPLIST

        ELSE IF ( CHOICE .EQ. 2 ) THEN			!  select a state

            WRITE( *,92010 )
     &          'The available state codes are:',
     &          ( S, STALIST( S ), S = 1, NSTALIST )
            S = GETNUM( 1, NSTALIST, 1, 
     &                  'Enter number for state selected' )
            LOLIST = STAINDX( S )
            HILIST = STAINDX( S+1 ) - 1

        ELSE IF ( CHOICE .EQ. 3 ) THEN			!  select a fip

            WRITE( *,92010 )
     &          'The available state codes are:',
     &          ( S, STALIST( S ), S = 1, NSTALIST )
            S = GETNUM( 1, NSTALIST, 1, 
     &                  'Enter number for state selected' )
            LO = STAINDX( S )
            HI = STAINDX( S+1 ) - 1
            WRITE( *,92010 )
     &          'The available county codes for this state are:',
     &          ( F, FIPLIST( F ), F = LO, HI )
            F = STAINDX( S ) + GETNUM( 1, NSTALIST, 1, 
     &                                 'Enter number for state' )
            LOLIST = F
            HILIST = F

        END IF			!  if all counties, else if state, else if fip


C...........   Process gridded county-level surrogates

        WRITE( *,92000 ) 
     &  ' ', 
     &  'Constructinging GRIDDED SURROGATE FRACTIONS output file...',
     &  ' '

        JDATE = 0
        JTIME = 10000

        DO  299  F = LOLIST, HILIST                !  loop on counties

            DO  222  S = 1, NASRG
            DO  211  C = 1, NGRID
                GRID( C,S ) = 0.0
211         CONTINUE
222         CONTINUE

            LO = FIPINDX( F   )
            HI = FIPINDX( F+1 ) - 1

            DO  288  K = LO, HI         !  loop on records K for this fip

                C = CELL( K )
                S = SIDS( K )
                GRID( C,S ) = GRID( C,S )  +  FRAC( K )

288         CONTINUE                    !  end loop on records K for this fip

            IF ( WRITE3( GNAME, 'ALL', JDATE, JTIME, GRID ) ) THEN

                WRITE( MESG,94020 ) 
     &              'Writing', FIPLIST( F ),  
     &              'to "' // GNAME( 1 : TRIMLEN( GNAME ) ) // 
     &              '" as record', F - LOLIST + 1
                CALL M3MESG( MESG )

            ELSE                        !  else write3() failed

                WRITE( MESG,94020 ) 
     &              'Could not write', FIP,  
     &              'to "' // GNAME( 1 : TRIMLEN( GNAME ) ) // '"'
                CALL M3EXIT( 'SRGGRID', 0, 0, MESG, 2 )

            END IF                      !  if write3() succeeded, or not
            
            CALL NEXTIME( JDATE, JTIME, 10000 )

299     CONTINUE                        !  end loop on counties F


C...............   End of program

        CALL M3EXIT( 'SRGGRID', 0, 0, 
     &               'SUCCESSFUL COMPLETION of program SRGGRID', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

92010   FORMAT( /5X, A , /, 10 ( /, 5( 4X, I3, ':', I5, : ) ) )


C...........   Formatted file I/O formats............ 93xxx


C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )

94020   FORMAT( A, :, I7, :, 2X, A, :, I10 )


        END

