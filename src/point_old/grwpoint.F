
        PROGRAM GRWPOINT

C***********************************************************************
C  program body      starts at line  224
C  subroutine NORMID starts at line  946
C
C  DESCRIPTION:
C       Construct monthly point source PNTS file(s) from monthly
C       projection factors by stid cyid pltid and scc.  If a source
C       does not have a projection factor match from this file it
C       uses a default factor from a netCDF input file (e.g. divide
C       1995 PNTS/1990 PNTS to get a default factor).
C
C  PRECONDITIONS REQUIRED:
C       Monthly Growth/control factors ASCII file 
C       Monthly Growth/control factor netCDF file 
C       SMOKE PNTS file
C       SMOKE SRCLIST file
C     
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND4, PROMPTFFILE, PROMPTMFILE, GETMENU, GETNUM, GETYN, TRIMLEN
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

        INCLUDE 'PTDIMS3.EXT'   !  point-source dimensioning parameters
        INCLUDE 'CHDIMS3.EXT'   !  emis chem parms (inventory + model)
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         FINDC
        INTEGER         GETMENU
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         STR2INT
        REAL            STR2REAL
        REAL            YR2DAY

        EXTERNAL        FINDC, GETMENU, GETNUM, GETYN,
     &                  PROMPTFFILE, PROMPTMFILE, TRIMLEN, STR2INT,
     &                  STR2REAL, YR2DAY

C...........   PARAMETERS and their descriptions:

        INTEGER         MXLST               ! max plant|scc replace-list size
        INTEGER         MXMON               ! no. of. months of fuel use input
        INTEGER         MXPLT               ! max XREF fip/scc/plant count
        INTEGER         MXSPPL              ! max no. sources per fip/plant

        PARAMETER     ( MXLST  = 5000 ,
     &                  MXMON  =   12 , 
     &                  MXPLT  =  700 ,
     &                  MXSPPL =   60 )

C...........   LOCAL VARIABLES and their descriptions:

C...........   Growth factors from point Sources Table from GNAME

        REAL          NOXGR ( NPSRC )  !  NOx growth factors 
        REAL          COGR  ( NPSRC )  !  CO growth factors 
        REAL          ROGGR ( NPSRC )  !  ROG growth factors 
                        
C...........   Inventory variables

        INTEGER        IFIP  ( NPSRC )  !  source FIPS (county) ID
        INTEGER        ISCC  ( NPSRC )  !  source SCC
        INTEGER        ISIC  ( NPSRC )  !  source SIC
        INTEGER        IPLANT( NPSRC )  !  Plant ID
        INTEGER        ISTACK( NPSRC )  !  Stack ID
        INTEGER        TZONES( NPSRC )  !  time zones
        INTEGER        TPFLAG( NPSRC )  !  temporal profile types
        INTEGER        INVYR ( NPSRC )  !  inventory year for this record
        REAL           XLOCA ( NPSRC )  !  UTM X-location (m)
        REAL           YLOCA ( NPSRC )  !  UTM Y-location (m)
        REAL           STKHT ( NPSRC )  !  stack height   (m)
        REAL           STKDM ( NPSRC )  !  stack diameter (m)
        REAL           STKTK ( NPSRC )  !  exhaust temperature (deg K)
        REAL           STKVE ( NPSRC )  !  exhaust velocity    (m/s)
        REAL           CTLEFF( NPSRC, NIPOL )  !  control efficiency fraction
        REAL           RULPEN( NPSRC, NIPOL )  !  rule penetration   fraction
        REAL           RULEFF( NPSRC, NIPOL )  !  rule effectiveness fraction
        REAL           EMISV ( NPSRC, NIPOL )  !  emissions values (tons/yr)

        COMMON /EMISREC / IFIP, ISIC, ISCC, IPLANT, ISTACK, TZONES,
     &                    TPFLAG, INVYR, XLOCA, YLOCA,
     &                    STKHT, STKDM, STKTK, STKVE,
     &                    CTLEFF, RULEFF, RULPEN, EMISV

        REAL           EMIS( NPSRC, NIPOL )  !  tmp emis values (tons/yr)

        CHARACTER*15   PLIDS( NPSRC )        !  EMS-95 plant ID for all srcs
        CHARACTER*40   PDESC( NPSRC )        !  plant description

C...........   Input monthly adjustments and associated variables

        INTEGER        NPLTSRC          !  number of source matches to XREF
        INTEGER        NSCCPLT          !  size of fip-scc-plant table
        INTEGER        INDXF ( MXPLT )  !  sorting index for XRFFPA
        INTEGER        INDXP ( MXPLT )  !  sorting index for XRFFSPA
        INTEGER        XRFSCC( MXPLT )  !  identifier for scc matchup

        REAL           RMFAC ( MXMON, MXPLT ) ! monthly factors for fip/scc/plt

        CHARACTER*15   FPLID  ( MXPLT ) !  EMS-95 plant code
        CHARACTER*28   XRFFSPA( MXPLT ) !  identifier for fip-scc-plant matchup
        CHARACTER*28   XRFFSP ( MXPLT ) !  sorted XRFFSPA

C...........   For report by facility

        INTEGER        NPLT             !  size of fip-plant table 
        INTEGER        CNTSP  ( MXPLT ) !  count of sources for each FIP/plant
        INTEGER        SNUMPPL( MXPLT, MXSPPL ) ! sources for each FIP/plant

        REAL           SUMNOX ( MXPLT, MXMON )  ! plant sum of monthly NOx
        REAL           SUMINV1( MXPLT ) !  plant sum of inventory year one NOx
        REAL           SUMINV2( MXPLT ) !  plant sum of inventory year two NOx

        LOGICAL        LPPL   ( MXPLT ) !  true if any srcs matched XRFFSP

        CHARACTER*15   PLTPPL ( MXPLT ) !  full plant ID names
        CHARACTER*20   XRFFPA ( MXPLT ) !  identifier for fip-plant matchup
        CHARACTER*20   XRFFP  ( MXPLT ) !  sorted XRFFPA

C...........   Sources from PSRC in monthly adjustments

        INTEGER        PLTSRC( MXLST ) !  fip/scc/plant record subscript
        INTEGER        SRCPLT( MXLST ) !  matching src  record subscript

C...........   LOGICAL NAMES and UNIT NUMBERS:

        INTEGER         LDEV           !  log-device
        INTEGER         SRDEV          !  sources list file device
        INTEGER         XRDEV          !  x-reference file device
        INTEGER         FDEV           !  facility report-file device
        INTEGER         RDEV           !  sources report-file device

        CHARACTER*16    ENAME          !  point-source input file logical name 
        CHARACTER*16    GNAME          !  growth factors file (netCDF)
        CHARACTER*16    MNAME( MXMON ) !  ctl matrix  output file logical names
        CHARACTER*16    NAME           !  scratch filename 
        CHARACTER*16    VNAME

C...........   OTHER LOCAL VARIABLES :

        INTEGER        I, J, K, L   !  scratch indices
        INTEGER        M, S, V         !  scratch loop variables
        INTEGER        XCNT, ECNT      !  scratch counters

        INTEGER        AMON, ZMON      !  first, last month to be averaged
        INTEGER        IOS             !  I/O status
        INTEGER        IREC            !  record # from file
        INTEGER        JFIP, JSCC      !  scratch input vars from pt src list 
        INTEGER        KMAX            !  
        INTEGER        MONCNT          !  count of months processed
        INTEGER        NEWYR           !  New inventory year
        INTEGER        NSRC            !  # recs, pt src files 
        INTEGER        NXRFOB          !  size of fip-oris-boiler table
        INTEGER        SMKPL, SMKST    !  temporary SMOKE plant and stack IDs

        REAL           ROLDEM          !  temporary old emissions value
        REAL           PERC            !  percent difference
        REAL           AMAT( MXMON, NPSRC ) !  adjustment coefficients
        REAL           YRINVNEW        !  yr 2 day conversion for new year
        REAL           YRINVOLD        !  yr 2 day conversion for old year

        LOGICAL        EFLAG          !  flag:  input error
        LOGICAL        FFLAG          !  flag:  make facility report-file
        LOGICAL        RFLAG          !  flag:  make sources report-file

        CHARACTER*5    CFIP            !  temporary character FIPS code 
        CHARACTER*8    SCC             !  temporary SCC code
        CHARACTER*12   DEVID           !  temporary EMS95  device ID
        CHARACTER*12   STKID           !  temporary EMS95 stack ID 
        CHARACTER*15   PLID,PLID2      !  temporary OTAG plant id on xrf file
        CHARACTER*20   KEY20, SCR20    !  temporary fips//plant
        CHARACTER*28   KEY28, SCR28    !  temporary fips//scc//plant
        CHARACTER*40   DUMNAM          !  temporary plant name
        CHARACTER*90   LINE            !  file line buffer
        CHARACTER*256  MESG            !  error message buffer

C.........  DATA STATEMENTS

        CHARACTER*3    MON( MXMON )    ! month for output file logical names
        DATA           MON / 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
     &                       'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /


C***********************************************************************
C   begin body of program GRWPOINT

        LDEV = INIT3()

        CALL INITEM( LDEV )

        WRITE( *,92000 ) 
     &' ',
     &  'Program GRWPOINT takes a SMOKE POINT INVENTORY file, an ASCII',
     &  'POINT SOURCES LIST file, a monthly adjustment factors file,',
     &  'and a growth factors file, to produce MONTHLY INVENTORY files',
     &  'for each month requested and optional ASCII REPORTS.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to have set the logical names for the input and',
     &  'output files (using "setenv <logicalname> <pathname>") prior',
     &  'to program launch, and will need to enter these logical names',
     &  'at the appropriate prompts.', 
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses to prompts',
     &  'are indicated in brackets [LIKE THIS].',
     &  ' '

        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0, 'Ending program GRWPOINT', 2)
        END IF

C....... Get file name; open input pt src, src list, Xref and fuel use files 

        ENAME = PROMPTMFILE(
     &           'Enter logical name for POINT SOURCES input file',
     &           FSREAD3, 'PNTS', 'GRWPOINT' )

C.......   Read in description of point-source emissions file

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &          'Could not get description of file ' // ENAME, 3 )
        END IF

C.......  Get other input files

        SRDEV = PROMPTFFILE( 
     &           'Enter logical name of the FACILITY/STACK input file',
     &           .TRUE., .TRUE., 'PSRC', 'GRWPOINT' )

        XRDEV = PROMPTFFILE( 
     &          'Enter logical name of MONTHLY ADJUSTMENTS input file',
     &          .TRUE., .TRUE., 'PASCFAC', 'GRWPOINT' )

C.........  Get date for projection and months for processing

        NEWYR = GETNUM( 1980, 2020, 1995 ,
     &         'New inventory year (projection year)' )

        AMON = GETMENU( 12, 5, 
     &         'Enter first month for factors in PASCFAC', MON )

        ZMON = GETNUM( AMON, 12, 12 ,
     &         'Enter last month for factors in PASCFAC' )
C
C.........  Open output inventory files for each month.
C.........  Reuse much of input PNTS description
C
        DO  44  I = AMON, ZMON

            FDESC3D( 1 ) = FDESC3D( 1 )( 1:TRIMLEN( FDESC3D( 1 ) ) ) // 
     &                    ' For ' // MON( I ) // '.'
            NAME = 'PNTS_' // MON( I )
            MESG = 'Enter logical name of ' // MON( I ) // 
     &             'POINT SOURCES output file' 
            MNAME( I ) = PROMPTMFILE( MESG, FSUNKN3, NAME, 'GRWPOINT' )

44      CONTINUE
C
        RDEV  = PROMPTFFILE( 
     &          'Enter logical name of SOURCE REPORT file, or "NONE"',
     &          .FALSE., .TRUE., 'PAREPT', 'GRWPOINT' )

        RFLAG = ( RDEV .GT. 0 )

        FDEV = PROMPTFFILE( 
     &         'Enter logical name of FACILITY REPORT file, or "NONE"',
     &         .FALSE., .TRUE., 'PAFCRPT', 'GRWPOINT' )

        FFLAG = ( FDEV .GT. 0 )

        DO  51   J = 1, MXPLT
        DO  50   M = 1, MXMON
            RMFAC( M,J ) = 0.0
50      CONTINUE
51      CONTINUE

C...........   Read the ASCII factors file
C.......   Format is as follows
C.......       State FIP                I2      cols   1: 2
C.......       County FIP               I3      cols   3: 5
C.......       OTAG Plant ID           A15      cols   6:20
C.......       SCC                      I8      cols  21:28
C.......       monthly factors   NMON*F8.4      cols  29:36 and so on
C............. 

        CALL M3MESG( 'Processing monthly adjustments...' )

        ECNT = 0
        IREC = 0
        J    = 0

77      CONTINUE

            READ ( XRDEV, 93000, END=88, IOSTAT=IOS ) LINE

            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 )
     &              'Error', IOS, 
     &              'reading monthly adjustments file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                ECNT  = ECNT + 1
                IF ( ECNT .GT. 10 ) THEN
                    CALL M3EXIT( 'GRWPOINT', 0, 0, 
     &                           'Too many errors', 2 )
                END IF
                GO TO  77  

            ELSEIF ( LINE(  1: 2 ) .EQ. ' ' .OR.
     &           LINE(  3: 5 ) .EQ. ' ' .OR.
     &           LINE(  6:20 ) .EQ. ' ' .OR.
     &           LINE( 21:28 ) .EQ. ' ' ) THEN
               EFLAG = .TRUE.
               WRITE( MESG, 94010 ) 'Missing ID field in monthly ' //
     &                              'adjustments file at line', IREC
               CALL M3MESG( MESG )
               GO TO 77  

            END IF                      !  if field missing
            
            JFIP = 1000 * STR2INT( LINE( 1:2 ) ) +   ! Convert to contain 0s
     &                    STR2INT( LINE( 3:5 ) )
            WRITE( CFIP,94020 ) JFIP
            CALL NORMID( PLID, LINE( 6:20 ) )
            SCC         = LINE( 21:28 )

            J = J + 1 

            IF( J .LE. MXPLT )THEN

                XRFSCC ( J ) = STR2INT( SCC )
                INDXP  ( J ) = J
                XRFFSPA( J ) = CFIP // SCC // PLID

C.................  Assign factors to specific month
                K = 29
                L = 36
                DO 81 M = AMON, ZMON
                    RMFAC( M,J ) = AMAX1( 0.0, STR2REAL( LINE( K:L ) ) )
                    K = K + 8
                    L = L + 8
81              CONTINUE

            ENDIF
              
            GO TO 77                    !  to head of loop reading xref file
                   
88      CONTINUE                        !  exit from loop reading xref file

        NXRFOB = J

        IF ( NXRFOB .GT. MXPLT ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &        'Monthly adjustments table overflow; actual:', NXRFOB,
     &        'dimensioned MXPLT=', MXPLT
            CALL M3MSG2( MESG )

        ELSE
            WRITE( MESG,94010 ) 
     &          'Number of valid monthly adjustment records:         ',
     &          NXRFOB
            CALL M3MESG( MESG )

        END IF

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading monthly adjustments file', 2 )
        END IF

C...........   Create sorted, unique fip-scc-plant table:

        CALL SORTIC( NXRFOB, INDXP, XRFFSPA )             !  sort on plants
        KEY28 =  CMISS3		!  initialize previous fip/scc/plant
        J     =  0		!  fip/scc/plant count

        DO  101  I = 1, NXRFOB                   ! loop through xref records

            K     = INDXP  ( I )
            SCR28 = XRFFSPA( K )

            IF ( SCR28 .NE. KEY28 ) THEN

                J = J + 1
                XRFFSP( J ) = SCR28 ! J must be < MXPLT b/c NXRFOB is
                KEY28       = SCR28

                INDXF ( J ) = J     ! Store list of FIP // PLID
                XRFFPA( J ) = SCR28( 1:5 ) // SCR28( 14:28 )

            END IF                  !  if new key encountered
            
101     CONTINUE                !  end loop initializing fip-plants table

        NSCCPLT = J

        WRITE( MESG,94010 ) 
     &  'Number of factors fips/scc/plt:                ', NSCCPLT
        CALL M3MESG( MESG )

C...........   Create sorted, unique fip-plant table:

        CALL SORTIC( NSCCPLT, INDXF, XRFFPA )             !  sort on plants
        KEY20 =  CMISS3		!  initialize previous fip/scc/plant
        J     =  0		!  fip/scc/plant count

        DO  103  I = 1, NSCCPLT                  ! loop through xref records

            K     = INDXF ( I )
            SCR20 = XRFFPA( K )

            IF ( SCR20 .NE. KEY20 ) THEN

                J = J + 1
                XRFFP( J ) = SCR20  ! J must be < MXPLT b/c NXRFOB is
                KEY20      = SCR20

            END IF                  !  if new key encountered
            
103     CONTINUE                !  end loop initializing fip-plants table

        NPLT = J

        WRITE( MESG,94010 ) 
     &  'Number of factors fips/plt:                ', NPLT
        CALL M3MESG( MESG )

C.........  Read in netCDF file with growth factors (not MONTHLY)

        GNAME = PROMPTMFILE(
     &           'Enter logical name for GROWTH FACTORS file',
     &           FSREAD3, 'PNTGRW', 'GRWPOINT' )

            IF ( .NOT. READ3( GNAME, 'NOX', ALLAYS3, 0, 0,
     &                        NOXGR ) ) THEN
                CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                       'Could not read NOX from ' // GNAME, 2 )
            END IF

            IF ( .NOT. READ3( GNAME, 'CO', ALLAYS3, 0, 0,
     &                        COGR ) ) THEN
                CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                       'Could not read CO from ' // GNAME, 2 )
            END IF

            IF ( .NOT. READ3( GNAME, 'ROG', ALLAYS3, 0, 0,
     &                        ROGGR ) ) THEN
                CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                       'Could not read ROG from ' // GNAME, 2 )
            END IF


C...........   Read SMOKE::EMS95 point source list file.
C.......   Strategy:  use NOx growth factors if available and
C.......   and override with monthly growth factors from PASCFAC.

        CALL M3MESG( 'Processing SOURCE LIST file...' )


        DO 106 M = AMON, ZMON
            DO 105 K = 1, NPSRC
                AMAT( M, K ) = NOXGR( K )         ! no-match control matrix
105         CONTINUE
106     CONTINUE

C.........  Initialize count of sources per EMS-95 facility code
        DO 109 J = 1, NPLT
            CNTSP( J ) = 0
            LPPL ( J ) = .FALSE.
109     CONTINUE

        EFLAG = .FALSE.
        IREC  = 0
        KMAX  = 0
        S     = 0 ! source counter
        XCNT  = 0
        
111     CONTINUE       !  head of loop reading PSRC file

            READ ( SRDEV, 93070, END=115, IOSTAT=IOS ) 
     &             JFIP, JSCC, SMKPL, SMKST, PLID2, STKID, DEVID, DUMNAM

            IREC = IREC + 1          

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 )
     &              'Error', IOS,
     &              'reading SOURCE LIST file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                ECNT  = ECNT + 1
                IF ( ECNT .GT. 10 ) THEN
                    CALL M3EXIT( 'GRWPOINT', 0, 0, 'Too many errors', 2)
                END IF
                GO TO  111
            END IF 

            S = S + 1

            IF ( S .LE. NPSRC ) THEN

                PDESC( S ) = DUMNAM
                PLIDS( S ) = PLID2

                WRITE( CFIP, 94020 ) JFIP    
                WRITE( SCC , 94015 ) JSCC
                CALL NORMID( PLID, PLID2 )  

C.................  Find source's FIP/PLID in list from monthly adjustments
C.................  and initialize flag 
                KEY20 = CFIP // PLID
                I = FINDC( KEY20, NPLT, XRFFP )    !find fip/plant match

C.................  Determine if SMOKE source is in list of monthly adjustments
                KEY28 = CFIP // SCC // PLID
                J     = FINDC( KEY28, NSCCPLT, XRFFSP )!find fip/scc/plant match

                IF( J .GT. 0 ) THEN

                    XCNT = XCNT + 1
                    IF ( XCNT .LE. MXLST ) THEN
                        PLTSRC( XCNT ) = J
                        SRCPLT( XCNT ) = S
                    END IF

                    FPLID( J ) = PLID2

C.....................  Set flag that says at least one source for this plant
C.....................  matches the monthly adjustments table
                    IF( I .GT. 0 ) LPPL ( I ) = .TRUE.

                END IF                      !  if fip/scc/plant match or not

            END IF      !  if source count S in bounds

            GO TO 111   !  to head of PSRC read loop

115      CONTINUE	!  exit from loop reading pt-source list file
        
        NSRC    = S
        NPLTSRC = XCNT

        WRITE( MESG,94010 ) 
     &  'Number of PSRC POINT SOURCES          ', NSRC
        CALL M3MESG( MESG )

        WRITE( MESG,94010 ) 
     &  'Number of PSRC::XREF PLANT matches    ', XCNT
        CALL M3MESG( MESG )

        IF ( NSRC .GT. NPSRC ) THEN
            EFLAG = .TRUE.
            CALL M3MSG2( 'Overflow reading PSRC file' )
        ELSE IF ( NSRC .LT. NPSRC ) THEN
            EFLAG = .TRUE.
            CALL M3MSG2( 'Underflow reading PSRC file' )
        END IF

        IF ( XCNT .GT. MXLST ) THEN
            EFLAG = .TRUE.
            CALL M3MSG2( 'Overflow building PLANT::SOURCE table' )
        END IF
          
        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0, 
     &                   'Error reading PSRC file', 2 )
        ENDIF

C.......   Write over the emissions where matches were found:

        CALL M3MESG( 'Storing EMISSIONS with matches found...')

        DO 122 M = AMON, ZMON
            DO 121 K = 1, NPLTSRC

                J = PLTSRC( K )
                AMAT( M, SRCPLT( K ) ) = RMFAC( M,J )

121         CONTINUE
122     CONTINUE                !  end loop on months M in the desired season

C
C...........   Read in emissions values.

        DO  133 V = 1, NIPOL

            IF ( .NOT. READ3( ENAME, EINAM( V ), ALLAYS3, 0, 0,
     &                        EMIS( 1,V ) ) ) THEN
                CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                       'Error reading "' //
     &                       EINAM( V )( 1:TRIMLEN( EINAM( V ) ) )//
     &                       '" from file "' //
     &                       ENAME( 1:TRIMLEN( ENAME ) ) // '".', 2 )
            END IF

            VNAME = 'CTLEFF_' // EINAM( V )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3,
     &                        0, 0, CTLEFF( 1,V ) ) ) THEN
                CALL M3EXIT( 'GRWPOINT', 0, 0,
     &               'Error reading CTLEFF from file' // ENAME, 2 )
            END IF

            VNAME = 'RULEFF_' // EINAM( V )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3,
     &                        0, 0, RULEFF( 1,V ) ) ) THEN
                CALL M3EXIT( 'GRWPOINT', 0, 0,
     &               'Error reading RULEFF from file' // ENAME, 2 )
            END IF

            VNAME = 'RULPEN_' // EINAM( V )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3,
     &                        0, 0, RULPEN( 1,V ) ) ) THEN
                CALL M3EXIT( 'GRWPOINT', 0, 0,
     &               'Error reading RULPEN from file' // ENAME, 2 )
            END IF

133     CONTINUE

C.........  Read inventory file characteristics

        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3,0,0, ISCC ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading ISCC from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading IPLANT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3,0,0, ISTACK ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading ISTACK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'TZONES', ALLAYS3,0,0, TZONES ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading TZONES from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'TPFLAG', ALLAYS3,0,0, TPFLAG ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading TPFLAG from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'INVYR', ALLAYS3,0,0, INVYR ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                  'Error reading INVYR from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'XLOCA', ALLAYS3,0,0, XLOCA ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading XLOCA from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'YLOCA', ALLAYS3,0,0, YLOCA ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading YLOCA from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKHT', ALLAYS3,0,0, STKHT ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading STKHT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKDM', ALLAYS3,0,0, STKDM ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading STKDM from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKTK', ALLAYS3,0,0, STKTK ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading STKTK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKVE', ALLAYS3,0,0, STKVE ) ) THEN
            CALL M3EXIT( 'GRWPOINT', 0, 0,
     &                   'Error reading STKVE from file' // ENAME, 2 )
        END IF


C.........  Store the source IDs associated with this FIP // PLID
C.........  if it is found and if there was at least one match to
C.........  the monthly adjustments table for this plant.
C.........  Must do this after IFIP has been read.
        DO 137 S = 1, NPSRC

            WRITE( CFIP, 94020 ) IFIP( S ) 
            CALL NORMID( PLID, PLIDS( S ) )  

            KEY20 = CFIP // PLID
            J     = FINDC( KEY20, NPLT, XRFFP )

            IF( J .GT. 0 ) THEN

                IF( LPPL( J ) ) THEN

                    K = CNTSP( J ) + 1

                    IF( K .LE. MXSPPL ) THEN
                        SNUMPPL( J,K ) = S
                    ELSE
                        KMAX = MAX( KMAX,K )
                    ENDIF

                    CNTSP ( J ) = K
                    PLTPPL( J ) = PLIDS( S )

                END IF 
            END IF      !  if fip/plant match or not

137     CONTINUE
          
        IF ( KMAX .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 
     &        'Sources per FIP/plant table overflow; actual:', KMAX,
     &        'dimensioned MXSPPL=', MXSPPL
            CALL M3MSG2( MESG )

            CALL M3EXIT( 'GRWPOINT', 0, 0, 
     &                   'Error reading PSRC file', 2 )
        ENDIF

C............................................................................
C.........  Calculate inverse of old projection years
        YRINVOLD = YR2DAY( INVYR( 1 ) )  ! Assume sources same year
        YRINVNEW = YR2DAY( NEWYR )

        DO 140 S = 1, NPSRC
            INVYR( S ) = NEWYR
140     CONTINUE

C.........  Use default factors for CO and ROG and default/specific factors
C.........  for NOx

        DO 166 M = AMON, ZMON

C.............  MRH NOTE: This may not be general for inventory (eg., CO must=1)
            DO 144 S = 1, NPSRC
                DO 143 V = 1, NIPOL

                    IF(V .EQ. 1) EMISV( S,V )= EMIS( S,V )* COGR ( S )
                    IF(V .EQ. 2) EMISV( S,V )= EMIS( S,V )* AMAT ( M,S )
                    IF(V .EQ. 3) EMISV( S,V )= EMIS( S,V )* ROGGR( S ) 

143             CONTINUE
144         CONTINUE

C.......   If rflag, write out report by source:
C.......   JV NOTE: This may need some more work later....
C
            IF ( RFLAG ) THEN

                WRITE( RDEV,93000 ) ' '
                WRITE( RDEV,93000 ) 'EMISSIONS BY SOURCE and MONTH'
                WRITE( RDEV,93000 ) 'units are tons/year'
                WRITE( RDEV,93000 ) ' '

                WRITE( RDEV,95100 ) MON( M )

                DO 155 L = 1, NPLTSRC
                    S = SRCPLT( L )
                    J = PLTSRC( L )

                    IF( MOD( L,30 ) .EQ. 1 ) WRITE( RDEV,95105 )

                    ROLDEM = EMIS( S,2 ) * NOXGR( S )
                    IF ( ROLDEM .GT. 0.0 )THEN
                        PERC= 100* ( ( EMISV( S,2 ) - ROLDEM )/ ROLDEM )
                    ELSE
                        PERC= 0.0
                    ENDIF

                    WRITE( RDEV,95106 ) 
     &               IFIP( S ), ISCC( S ), FPLID( J ), ISTACK( S ),
     &               EMIS( S,2 ), ROLDEM, EMISV( S,2 ), PERC, PDESC( S )

155             CONTINUE

                WRITE( RDEV,93000 ) ' '

            END IF   ! if rflag:  write report-file

C.............  If fflag, Sum emissions for all plants for this month and store
            IF ( FFLAG ) THEN

                DO 159 J = 1, NPLT

                    DO 157 K = 1, CNTSP( J )

                        IF( K .EQ. 1 ) SUMNOX( J,M ) = 0.  ! Initialize

                        S = SNUMPPL( J,K )

                        SUMNOX( J,M ) = SUMNOX( J,M ) + 
     &                                  EMISV( S,2 ) * YRINVNEW

157                 CONTINUE
159             CONTINUE

            ENDIF

            IF ( .NOT. WRITE3( MNAME( M ), ALLVAR3, 0, 0,  IFIP ) ) THEN
                CALL M3EXIT( 'GRWPOINT', 0, 0,
     &               'Error writing output file "' //
     &               MNAME( M )( 1:TRIMLEN( MNAME( M )) ) // '"',  2 )
            END IF

166     CONTINUE  

C.........  If FFLAG, write report by facility for all months

        IF( FFLAG ) THEN

C.............  Calculate the number of months for format statement 93100 
            MONCNT = ZMON - AMON + 1

C.............  Sum up emissions by plant for non-adjusted inventories
            DO 188 J = 1, NPLT
                DO 187 K = 1, CNTSP( J )

                    IF( K .EQ. 1 ) THEN  ! Initialize
                        SUMINV1( J ) = 0.
                        SUMINV2( J ) = 0.
                    ENDIF

                    S = SNUMPPL( J,K )

                    SUMINV1( J ) = SUMINV1( J ) + EMIS( S,2 ) * YRINVOLD
                    SUMINV2( J ) = SUMINV2( J ) + 
     &                             EMIS( S,2 ) * NOXGR( S ) * YRINVNEW

187             CONTINUE
188         CONTINUE

C.............  Write the report (no header)
            DO 201 J = 1, NPLT

                IF( CNTSP( J ) .GT. 0 ) THEN  ! Report on this FIP/plant

                    S    = SNUMPPL( J,1 )   ! Use first src per plant for name
                    CFIP = XRFFP  ( J )( 1:5 )
                    PLID = PLTPPL ( J )
                    L    = TRIMLEN( PDESC( S ) )

                    WRITE( FDEV, 93100 ) CFIP, PLID, 
     &                     ( SUMNOX( J,M ), M = AMON, ZMON ),
     &                     SUMINV2( J ), SUMINV1( J ), PDESC( S )( 1:L )
                ENDIF

201         CONTINUE

        ENDIF

        MESG = 'Successful completion of Program GRWPOINT'
        CALL M3EXIT( 'GRWPOINT', 0, 0, MESG, 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93070   FORMAT( I5.5, 1X, I8.8, 1X, 2(I6, 1X), A15, 1X, A12 , 1X, A12,
     &   1x,A40 )

93100   FORMAT( A5, 1X, A15, 1X, <MONCNT>( F8.2, 1X ), 
     &          2( F8.2, 1X ), A ) 

C...........   Internal bufferring  formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) ) 

94015   FORMAT( I8.8 ) 

94020   FORMAT( I5.5 )

94100   FORMAT( A5, A8, A15 )

C...........   Other formats........................... 95xxx

95100   FORMAT( 'for month :', A3, / )

95105   FORMAT('  FIP',4x,'SCC',6x,'PLTID',4x,'90EMIS',4x,
     &         '95t2',5x,'95new',4x,'%DIFF',3x,'NAME')

95106   FORMAT( I5, I9, 1X, A15, I8, 3F10.2, F7.1, 1X, A40 )

        END

C-------------------------------------------------------------------------

        SUBROUTINE NORMID( IDSTR, SRCSTR )

C***********************************************************************
C  FUNCTION:
C	Copies SRCSTR to IDSTR, exclusive of leading blanks and leading zeros 
C
C  PRECONDITIONS REQUIRED:
C	adequate string length for IDSTR (else SRCSTR gets truncated)
C
C  REVISION  HISTORY:
C	Prototype 2/97 by CJC
C
C***********************************************************************

      IMPLICIT NONE

C...........   ARGUMENTS and their descriptions:

        CHARACTER*(*)	IDSTR	!  result string (after copying,stripping
        CHARACTER*(*)	SRCSTR	!  source string

C...........   SCRATCH LOCAL VARIABLES and their descriptions:

        INTEGER		SRCLEN
        INTEGER		C
        CHARACTER*1	CH

C***********************************************************************
C   begin body of subroutine  NORMID

        SRCLEN = LEN( SRCSTR )

        DO  11  C = 1, SRCLEN	

            CH = SRCSTR( C:C )

            IF ( CH .NE. ' '  .AND. CH .NE. '0' ) THEN
                IDSTR = SRCSTR( C : SRCLEN )	!  copy rest of string
                GO TO  12
            END IF

11      CONTINUE

        IDSTR = ' '		!  fallback if SRCSTR is all blanks or zeros

12      CONTINUE		!  exit from loop

        RETURN
        END

