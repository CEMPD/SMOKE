
        PROGRAM GRDMOBIL

C***********************************************************************
C  program body starts at line 144
C
C  DESCRIPTION:
C       Use gridding-transform from program GRDMMAT to grid time-stepped,
C       source level mobile source emission files.
C
C  PRECONDITIONS REQUIRED:  
C       M3IO source-level area source emissions in canonical order used
C       by GRDMMAT
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O; PROMPTMFILE
C
C  REVISION  HISTORY:
C       Prototype  2/96 by CJC.
C
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
        
        INCLUDE 'GRDIMS3.EXT'   !  grid parameters
        INCLUDE 'MBDIMS3.EXT'   !  mobil-source dimensioning parameters
        INCLUDE 'TMDIMS3.EXT'   !  temporal dimensioning parameters
        INCLUDE 'CHDGEN3.EXT'   !  Generalized chemistry parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.


C...........   EXTERNAL FUNCTIONS and their descriptions:
        
        LOGICAL         DSCGRID
        INTEGER         ENVINT
        LOGICAL         ENVYN
        LOGICAL         GETYN
        CHARACTER*10    HHMMSS
        INTEGER         INDEX1
        CHARACTER*14    MMDDYY
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         WKDAY           !  day of week (1...7)

        EXTERNAL        DSCGRID, ENVINT, ENVYN, GETYN, HHMMSS, INDEX1,
     &                  MMDDYY, PROMPTFFILE, PROMPTMFILE, TRIMLEN, WKDAY


C...........   PARAMETER

        CHARACTER*16    NONE
        PARAMETER     ( NONE = 'NONE' )

C...........   LOCAL VARIABLES and their descriptions:
C...........   Mobil Sources input and output arrays (one variable at a time)
        
        REAL            EMISV( NMSRC, NVTYPE )  !  input emissions.
        REAL            EMISG( NGRID, NVTYPE )  !  gridded output emissions
        
C...........   Gridding Matrix, Control Matrix, Speciation Matrix
      
        INTEGER         NX( NGRID )
        INTEGER         IX( NMATX )
        REAL            CX( NMATX )
        
        COMMON  / GRIDMAT / NX, IX, CX

C...........   Inventory arrays
        REAL            EMISR ( NMSRC )   ! Emissions summed over layers

        INTEGER         IFIP  ( NMSRC )   ! FIPs  code for each source
        INTEGER         INVYR ( NMSRC )   ! Inventory year code for each source
        INTEGER         ISPTR ( NMSRC )   ! Pointer for each state
        INTEGER         STATES( MXSTA )   ! Unique list of state IDs
 
C...........   Report arrays
        REAL            EMSUM ( MXSTA, MXEMIS+1 ) ! Summed emissions by state
        INTEGER         VOLIDX( MXEMIS )
           
C...........   Other local variables

        INTEGER         BEGBLK  !  beg of array of blanks in output file desc
        INTEGER         BTIME
        INTEGER         IOS
        INTEGER         I, J, K, L, S, T, V
        INTEGER         JDATE, JTIME, TSTEP
        INTEGER         LDATE, LTIME
        INTEGER         LDEV
        INTEGER         LSTA               ! previous state code
        INTEGER         NSTA               ! actual number of states
        INTEGER         NSTEPS
        INTEGER         NSUM               ! No. of emis types to report
        INTEGER         NVOL               ! No. of emis types w/ volatile pol
        INTEGER         RDEV               ! output report unit number
        INTEGER         RTIME              ! time for reporting state totals
        INTEGER         STA

        LOGICAL         OUTPUTFLAG
        LOGICAL         ST_REPORT  ! Output state-based emissions totals

        CHARACTER*3     INVPOL  !  buffer for inventory pollutant
        CHARACTER*4     OUTDESC !  contains descriptor of output file type
        CHARACTER*4     OUTNAM  !  contains default  logical output name
        CHARACTER*5     OUTFMT  !  contains output report format
        CHARACTER*8     REPNAM  !  contains default logical report output name
        CHARACTER*16    ANAME   !  logical name for area-source input file
        CHARACTER*16    ENAME   !  logical name for emission output file
        CHARACTER*16    GNAME   !  logical name for grid matrix input file
        CHARACTER*16    INAME   !  logical name for inventory input file

        CHARACTER*16    SCRBUF  !  scratch name area
        CHARACTER*256   MESG    !  scratch message area


C***********************************************************************
C   begin body of program GRDMOBIL
        
        LDEV = INIT3()
        
        WRITE( *,92000 )
     &  ' ',
     &  'Program GRDMOBIL to take a sorted source level mobile source ',
     &  'emissions file and a gridding-transform matrix computed by',
     &  'program GRDMMAT, and produce gridded vehicle-typed, mode',
     &  'specific emissions.  The output time step structure of the ',
     &  'gridded file will be the same as the time step structure of',
     &  'the input area source emissions file.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You will need to enter the logical names for the input and',
     &  'output files (and to have set them prior to program launch,',
     &  'using "setenv <logicalname> <pathname>").  Use "NONE" as the',
     &  'name for the control-matrix file if you wish to omit control',
     &  'from the operations performed, or as the name of the output',
     &  'file if you want only to time program performance without the',
     &  'overhead of additional I/O.',
     &  ' '
        WRITE( *,92000 ) 
     &  'You may use END_OF-FILE (control-D) to quit the program',
     &  'during logical-name entry.  Default responses are indicated',
     &  'in brackets [LIKE THIS].',
     &  ' '
        
        IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
            CALL M3EXIT( 'GRDMOBIL', 0, 0, 'Ending Program', 2 )
        END IF
 
C.......   Get environment variable settings
        ST_REPORT = ENVYN( 'ST_REPORT', 'Per-state report writing flag',
     &                      .TRUE., IOS )

        RTIME = ENVINT( 'ST_REPTIME', 
     &                  'Per-State summation reporting time (HHMMSS)',
     &                   230000, IOS )

C.......   Get file names; open input gridding matrix, mobile source, and
C.......   output files.

        GNAME = PROMPTMFILE( 
     &          'Enter logical name for GRIDDING MATRIX file',
     &          FSREAD3, 'MGMAT', 'GRDMOBIL' )

C......... NOTE: Close to 80 character prompt + default display limit
        ANAME = PROMPTMFILE( 
     &          'Enter logical name for ANNUAL VMT / TIME-STEPPED ' // ! 49
     &          'MOBILE EMIS file',                                    ! 16
     &          FSREAD3, 'MTMP', 'GRDMOBIL' )

C.......   Build description of output file, and optionally open it:

        IF ( .NOT. DESC3( ANAME ) ) THEN
            CALL M3EXIT( 'GRDMOBIL', 0, 0, 
     &              'Could not get description of file ' // ANAME, 2 )

        ELSEIF( NROWS3D .NE. NMSRC ) THEN
            WRITE( MESG, 94010 )
     &      'Dimension mismatch.  MOBILE SOURCES file:', NROWS3D,
     &      'program:', NMSRC
            CALL M3EXIT( 'GRDMOBIL', 0, 0, MESG, 2 )
 
        END IF

        JDATE   = SDATE3D
        JTIME   = STIME3D
        TSTEP   = TSTEP3D

C.........  Annual inventory VMT file as input
        IF ( TSTEP .EQ. 0 ) THEN
            REPNAM = 'MGRDREP'
            NSUM   = 1
 
C.........  Hourly inventory file as input
        ELSE
            REPNAM = 'MGRTREP'
            NSUM   = NVARS3D
        END IF

C.......   Read FIPS code for reporting purpose
        IF( ST_REPORT ) THEN
 
C............  If we've read in temporal file, must get inventory file name
            IF ( TSTEP .NE. 0 ) THEN
 
                INAME = PROMPTMFILE(
     &                  'Enter logical name for MOBIL INVENTORY file',
     &                  FSREAD3, 'MOBL', 'GRDMOBIL' )
            ELSE
                INAME = ANAME
 
            ENDIF
 
            IF( .NOT. READ3( INAME, 'IFIP', ALLAYS3, 0, 0, IFIP ) ) THEN
 
                CALL M3EXIT( 'GRDMOBIL', 0, 0,
     &                       'Error reading FIP from file' // INAME, 2 )
 
            END IF
 
            IF( .NOT.
     &          READ3( INAME, 'INVYR', ALLAYS3, 0, 0, INVYR ) ) THEN
 
                CALL M3EXIT( 'GRDMOBIL', 0, 0,
     &               'Error reading INVYR from file' // INAME, 2 )
 
            END IF
 
            RDEV  = PROMPTFFILE(
     &             'Enter logical name for TOTALS REPORT or "NONE"',
     &             .FALSE., .TRUE., REPNAM, 'GRDMOBIL' )
 
C............  Create state codes arrays
            LSTA = -9
            K    = 0
            DO 22 I = 1, NMSRC
 
                STA = IFIP( I ) / 1000
 
                IF( STA .GT. LSTA ) THEN
 
                    K = K + 1
                    STATES( K ) = STA
                    LSTA        = STA
 
                ENDIF
 
                ISPTR ( I ) = K
 
22          CONTINUE
 
            NSTA = K
 
            IF( NSTA .GT. MXSTA ) THEN
                WRITE( MESG,94010 )
     &                'Number of states is ', NSTA,
     &                'but maximum (MXSTA) is', MXSTA
                CALL M3EXIT( 'GRDMOBIL', JDATE, JTIME, MESG, 2 )
            ENDIF
 
C............  Initialize emission summary array
            DO 28 I = 1, MXSTA
                DO 27 S = 1, MXEMIS
                    EMSUM( I,S ) = 0.
27              CONTINUE
28          CONTINUE
 
        ENDIF

        IF ( .NOT.DSCGRID( GRDNM, SCRBUF, GDTYP3D, 
     &                     P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, 
     &                     XORIG3D, YORIG3D, XCELL3D, YCELL3D, 
     &                     NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            CALL M3EXIT( 'GRDMOBIL', 0, 0,
     &              'Could not get description of grid ' // GRDNM, 2 )

        END IF

C......... Store variable properties based on DESC3 call
        FTYPE3D = GRDDED3   !  shares most of file-description with input file.
        GDNAM3D = GRDNM
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3

C.........  Set single output variable if input is inventory VMT file
        IF ( TSTEP .EQ. 0 ) THEN
            NSTEPS  = 1
            FDESC3D( 1 ) = 'Mobile source gridded VMT'
            OUTNAM = 'MGRD'
            OUTDESC= 'VMT'
            BEGBLK= 2

            L = INDEX1( 'VMT', NVARS3D, VNAME3D )

            IF( L .LE. 0 ) THEN
 
                MESG = 'ERROR: Variable "VMT"' //
     &                 ' not found in input file!'
                CALL M3EXIT( 'GRDMOBIL', 0, 0, MESG, 2 )

            ELSE

                NVARS3D      = 1
                VNAME3D( 1 ) = VNAME3D( L )
                UNITS3D( 1 ) = UNITS3D( L )
                VTYPE3D( 1 ) = VTYPE3D( L )
                VDESC3D( 1 ) = VDESC3D( L ) 

            ENDIF
           
C.........  Otherwise, use same output variables as input variables 
        ELSE  
            NSTEPS  = MXREC3D
            FDESC3D( 1 ) = 'Mobile source gridded emission types'

C.............  Copy vehicle type descriptions from MTMP file
            DO 31 I = 1, NVTYPE
                FDESC3D( 1 + I ) = FDESC3D( 3 + I )
31          CONTINUE

            BEGBLK = NVTYPE + 2

            OUTNAM = 'MGRT'
            OUTDESC= 'EMIS'

        END IF

        DO 33 I = BEGBLK, MXDESC3
            FDESC3D( I ) = ' '
33      CONTINUE 
        
        ENAME = PROMPTMFILE( 
     &          'Enter logical name for OUTPUT ' // OUTDESC // 
     &          ' file or "NONE"', FSUNKN3, OUTNAM, 'GRDMOBIL' )

        OUTPUTFLAG = ( ENAME .NE. NONE )

C.......   Create pointer for which input variables are volatile
        I    = 0
        NVOL = 0 
        IF( NVARS3D .GT. 1 ) THEN
            DO 44 V = 1, NVARS3D

C.................  Check if pollutant part of name is in volatile pol list
                L = INDEX  ( VNAME3D( V ), '_' ) 
                J = TRIMLEN( VNAME3D( V ) ) 
                K = INDEX1 ( VNAME3D( V )( L+1: J ), NUMVOL, VOLTYP )

C.................  If so, store pointer
                IF( K .GT. 0 ) THEN
                    I = I + 1
                    VOLIDX( I ) = V
                    INVPOL = VNAME3D( V )( L+1: J )
                ENDIF 

44          CONTINUE
            NVOL = I
            NSUM = NSUM + 1

C..........   Create title for the statelst report, and use VNAME3D array for it
            IF( NVOL .GT. 0 ) VNAME3D( NSUM ) = 'TOT_' // INVPOL

        ENDIF

C.......   Read gridding, speciation, and (optional) control matrices:

        CALL M3MSG2( 'Reading GRIDDING matrix...' )

        IF ( .NOT. READ3( GNAME, 'ALL', 1, 0, 0, NX ) ) THEN

            CALL M3EXIT( 'GRDMOBIL', 0, 0,
     &           'Could not read gridding matrix from file "' //
     &           GNAME( 1 : TRIMLEN( GNAME ) ) // '".', 2 )
            
        END IF      !  if read3() failed for gridding matrix

C.......   Transform and write out area source emissions values:

        CALL M3MSG2( 'Calculating gridded categorized VMT...' )

        BTIME = JTIME
        LDATE = 0
        DO  199  T = 1, NSTEPS

C...............   If this is a new month, or new day, write message
            IF ( JDATE .GT. 0 .AND. LDATE. NE. JDATE ) THEN
 
                MESG = 'Processing ' //
     &                 DAYS( WKDAY( JDATE ) ) // MMDDYY( JDATE )
                CALL M3MSG2( MESG( 1:TRIMLEN( MESG ) ) )
 
            END IF

C.............  Write to screen because WRITE3 only writes to LDEV
            IF( JDATE .GT. 0 ) WRITE( *, 93020 ) HHMMSS( JTIME )

            DO  188  V = 1, NVARS3D

C...................   Read in entire emissions field

                IF( .NOT. READ3( ANAME, VNAME3D( V ), ALLAYS3, 
     &                           JDATE, JTIME, EMISV )         ) THEN

                    MESG = 'Could not read "' //
     &                     VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) ) //
     &                     '" from ' // ANAME
                    CALL M3EXIT( 'GRDMOBIL', JDATE, JTIME, MESG, 2 )

                END IF          !  if read3() succeeds or not

                DO  155  L = 1, NLAYS3D

                    CALL SMATVEC( NMSRC, NGRID, NMATX, NX, IX, CX,
     &                            EMISV( 1,L ), EMISG( 1,L ) )

155             CONTINUE

                IF( OUTPUTFLAG ) THEN

                    IF ( .NOT. WRITE3( ENAME, VNAME3D( V ),
     &                                 JDATE, JTIME, EMISG  ) ) THEN

                        MESG = 'Could not write "' //
     &                         VNAME3D( V )( 1:TRIMLEN( VNAME3D( V ) ) )
     &                         // '" to ' // ENAME
                        CALL M3EXIT( 'GRDMOBIL', JDATE, JTIME, MESG, 2 )

                    END IF          !  if write3() failed

                ENDIF

C.........................  When reporting is switched on, calculate emissions
C.........................  summaries by state and report to output QA file
C.........................  EMSUM must be initialized
                IF( ST_REPORT ) THEN

C.....................  Sum emissions across vehicle types
                    DO 177  L = 1, NLAYS3D
                        DO 175 S = 1, NMSRC
                            EMISR( S ) = EMISR( S ) + EMISV( S,L )
175                     CONTINUE
177                 CONTINUE
 
                    CALL AGRDSUM( NMSRC, NGRID, NMATX, NX, IX,
     &                            CX, ISPTR, EMISR, EMSUM( 1,V ) )
 
                    DO 181 S = 1, NMSRC
                        EMISR( S ) = 0.
181                 CONTINUE

                ENDIF

188         CONTINUE            ! end loop on input variables

            LDATE = JDATE
            LTIME = JTIME

            CALL NEXTIME( JDATE, JTIME, TSTEP )

C.............  Write report if report flag is switched on and previous
C.............  time equals the report time
            IF( ST_REPORT .AND.
     &          ( LTIME .EQ. RTIME .OR. T .EQ. NSTEPS ) ) THEN
 
C.................  Convert summary emissions from moles/day to tons/day
 
                IF( TSTEP .EQ. 0 ) THEN
                    MESG = 'Mobile source VMT (miles/ave ' //
     &                     'inventory day) within grid ' // GRDNM
                    OUTFMT = 'F12.0'

                ELSE 
                    IF( RDEV .GT. 0)
     &              WRITE( RDEV, 94010 ) DAYS( WKDAY( LDATE ) ) //
     &                                   MMDDYY( LDATE ) // ' from',
     &                                   BTIME, 'to', LTIME

                    MESG = 'Mobile source emissions (tons/' //
     &                     'specific day) within grid ' // GRDNM
                    OUTFMT = 'F8.1'

C.................  Sum volatile component
                    IF( NVOL .GT. 1 ) THEN
                        DO 194 I = 1, MXSTA
                            DO 193 J = 1, NVOL
                                EMSUM( I, NSUM ) = EMSUM( I, NSUM ) + 
     &                                             EMSUM( I, VOLIDX(J) )
193                         CONTINUE
194                     CONTINUE
                    ENDIF

                ENDIF
                IF( RDEV .GT. 0)  
     &          CALL STATELST( RDEV  , MESG , MXSTA , NSTA , NSUM,
     &                         OUTFMT, VNAME3D, STATES, EMSUM        )
 
C.................  Reinitialize emissions summary
                DO 198 I = 1, MXSTA
                    DO 197 S = 1, NSUM
                        EMSUM( I,S ) = 0.
197                 CONTINUE 
198             CONTINUE

C.................  Set begining time for next day
                BTIME = JTIME
 
            ENDIF

199     CONTINUE          !  end loop on time steps

999     CONTINUE          !  exit program
      
        CALL M3EXIT( 'GRDMOBIL', 0, 0, 
     &               'Normal completion of PROGRAM GRDMOBIL', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )

C...........   Formatted file I/O formats............ 93xxx
 
93020   FORMAT( 8X, 'at time ', A8 )
 

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT ( 10 ( A, :, I10, :, 2X ) )

        END

