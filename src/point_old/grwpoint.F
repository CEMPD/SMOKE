
        PROGRAM GRWPMAT

C***********************************************************************
C  program body      starts at line  189
C  subroutine NORMID starts at line 1179
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
C COPYRIGHT (C) 1997, MCNC--North Carolina Supercomputing Center
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
        INTEGER         FIND1
        INTEGER         GETMENU
        INTEGER         GETNUM
        LOGICAL         GETYN
        INTEGER         PROMPTFFILE
        CHARACTER*16    PROMPTMFILE
        INTEGER         TRIMLEN
        INTEGER         STR2INT
        REAL            STR2REAL

        EXTERNAL        FINDC, FIND1, GETMENU, GETNUM, GETYN,
     &                  PROMPTFFILE, PROMPTMFILE, TRIMLEN, STR2INT,
     &                  STR2REAL

C...........   PARAMETERS and their descriptions:

        INTEGER         MXMON               ! no. of. months of fuel use input
        INTEGER         MXPLT               ! max XREF fip/scc/plant count
        INTEGER         MXLST               ! max plant|scc replace-list size

        PARAMETER     ( MXMON =   12 , 
     &                  MXPLT =  700 ,
     &                  MXLST = 5000 )

C...........   Growth factors from point Sources Table from GNAME

        REAL          NOXGR ( NPSRC )  !  NOx growth factors 
        REAL          COGR  ( NPSRC )  !  CO growth factors 
        REAL          ROGGR ( NPSRC )  !  ROG growth factors 
                        

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER        AMON, ZMON      !  first, last month to be averaged
        INTEGER        NEWYR           !  New inventory year
        INTEGER        NSRC            !  # recs, pt src files 
        REAL           AMAT  ( MXMON, NPSRC ) !  adjustment coefficients
        INTEGER        JFIP, JSCC      !  scratch input vars from pt src list 
        INTEGER        SMKPL, SMKST    !    "       "     "       "      "  
        INTEGER        NPLTSRC
        INTEGER        PLTSRC( MXLST ) !  fip/scc/plant record subscript
        INTEGER        SRCPLT( MXLST ) !  matching src  record subscript
        INTEGER        SRCSCC( MXLST ) !  matching src  record subscript
        REAL           RMFAC ( MXMON, MXPLT ) !  factors by month for plt/scc
        INTEGER        NPLT            !  size of fip-plant table
        CHARACTER*28   XRFPL ( MXPLT ) !  identifier for fip-scc-plant matchup
        CHARACTER*28   XRFPLS( MXPLT ) !  sorted "
        INTEGER        INDXP ( MXPLT ) !  index array for XRFPL(S
        INTEGER        XRFSCC( MXPLT ) !  identifier for scc matchup
        INTEGER        KFIP( NPSRC ), KSCC( NPSRC )
        CHARACTER*15   FPLID( NPSRC )
        INTEGER        NXRFOB          !  size of fip-oris-boiler table
        CHARACTER*90   LINE            !  xref file line buffer
        CHARACTER*5    FIPS            !  FIPS code for state and county
        CHARACTER*8    SCC             !  SCC code
        CHARACTER*12   FSTKID(NPSRC) , STKID        !  EMS95 stack ID 
        CHARACTER*15   PLID,PLID2      !  scratch OTAG plant id on xrf file
        CHARACTER*12   DEVID           !   "       "   device "     "   "
        CHARACTER*5    ORISID          !  scratch doe ORIS id string  
        CHARACTER*6    BLRID           !  scratch doe boiler id string
        CHARACTER*28   KEY28, SCR28    !  scratch fips//scc//plant
        CHARACTER*40   DUMNAM, PDESC( NPSRC)   
        CHARACTER*4    MON( MXMON )    ! month for output file logical names
        DATA           MON /
     &          'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',  
     &          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /

        INTEGER         C, F, I, J, K, L  !  scratch indices
        INTEGER         M, S, V, Y  !  scratch loop variables
        INTEGER         LDEV           !  log-device
        INTEGER         SRDEV          !  sources list file device
        INTEGER         XRDEV          !  x-reference file device
        INTEGER         FUDEV          !  fuel use file device
        INTEGER         RDEV           !  report-file device
        INTEGER         IOS            !  I/O status
        INTEGER         IREC           !  record # from file
        INTEGER         XCNT, FCNT, ECNT
        INTEGER         MCNT, NCNT
        LOGICAL         EFLAG          !  flag:  input error
        LOGICAL         FFLAG          !  flag:  bad input fields
        LOGICAL         RFLAG          !  flag:  make report-file

        CHARACTER*16    ENAME          !  point-source input file logical name 
        CHARACTER*16    MNAME( MXMON ) !  ctl matrix  output file logical names
        CHARACTER*16    GNAME          !  growth factors file (netCDF)
        CHARACTER*16    NAME           !  scratch filename 
        CHARACTER*16    VNAME
        CHARACTER*256   MESG           !  error message buffer

        REAL            DDD, DIV, LGF1, LGF2, X, X2, LRAT, SUM, SIG
        REAL            ROLDEM, PERC

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
        REAL           RULPEN( NPSRC, NIPOL )  !  rule penetration   fraction
        REAL           RULEFF( NPSRC, NIPOL )  !  rule effectiveness fraction
        REAL           CTLEFF( NPSRC, NIPOL )  !  control efficiency fraction
        REAL           EMISV ( NPSRC, NIPOL )  !  emissions values (tons/yr)
        REAL           EMIS ( NPSRC, NIPOL )  !  tmp emis values (tons/yr)

        COMMON /EMISREC / IFIP, ISIC, ISCC, IPLANT, ISTACK, TZONES,
     &                    TPFLAG, INVYR, XLOCA, YLOCA,
     &                    STKHT, STKDM, STKTK, STKVE,
     &                    CTLEFF, RULEFF, RULPEN, EMISV


C***********************************************************************
C   begin body of program GRWPMAT

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &' ',
     &'Program GRWPMAT takes a SMOKE PNTS inventory file, an ASCII ',
     &'POINT SOURCES LIST file, an ASCII factors file, ',
     &'and netCDF factors file, to produce one GROWTH MATRIX file',
     &'of adjustment factors for each month requested, and an ',
     &'optional ASCII REPORT file.',
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
            WRITE( *,92000 ) ' ', 'Ending program GRWPMAT', ' '
            CALL EXIT( 0 )
        END IF


C....... Get file name; open input pt src, src list, Xref and fuel use files 

        ENAME = PROMPTMFILE(
     &           'Enter logical name for point source LIST file',
     &           FSREAD3, 'PNTS', 'GRWPMAT' )

C.......   Read in description of point-source emissions file

        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'GRWPMAT', 0, 0,
     &          'Could not get description of file ' // ENAME, 3 )
        END IF
C.......

        SRDEV = PROMPTFFILE( 
     &           'Enter logical name for SRCLIST file',
     &           .TRUE., .TRUE., 'PSRC', 'GRWPMAT' )

        XRDEV = PROMPTFFILE( 
     &           'Enter logical name for ASCII factors file',
     &           .TRUE., .TRUE., 'PASCFAC', 'GRWPMAT' )

C.......   Open control matrix file:  reuse much of POINT SOURCE description

C        FTYPE3D = GRDDED3


        RDEV  = PROMPTFFILE( '...for REPORT file, or "NONE"',
     &                       .FALSE., .TRUE., 'PAREPT', 'GRWPMAT' )
        RFLAG = ( RDEV .GT. 0 )

        NEWYR = GETNUM( 1980, 2020, 1995 ,
     &         'New inventory year (projection year) ' )

        AMON = GETMENU( 12, 5, 'Enter first months factors in PASCFAC',
     &         MON )

        ZMON = GETNUM( AMON, 12, 12 ,
     &         'Enter last months factors in PASCFAC ' )
C
C...... growth matrices for months specified above
C
        DO  44  I = AMON, ZMON

         NAME = 'PNTS_'//MON( I )
         MESG = 'Enter logical name for PNTS monthly file for ' 
     &           // MON(I)
         MNAME( I ) = PROMPTMFILE( MESG, FSUNKN3, NAME, 'GRWPMAT' )

44      CONTINUE
C

        DO  143  S = 1, MXPLT
        DO  141  M = 1, MXMON
            RMFAC( M,S ) = 0.0
141     CONTINUE
143     CONTINUE

C...........   Read the ASCII factors file
C.......   Format is as follows
C.......       State FIP                I2      cols   1: 2
C.......       County FIP               I3      cols   3: 5
C.......       OTAG Plant ID           A15      cols   6:20
C.......       SCC                      I8      cols  21:28
C.......       monthly factors   NMON*F8.4      cols  29:36 and so on
C............. 

        CALL M3MESG( 'Processing (Monthly) ASCII factors file' )

        ECNT  = 0
        J     = 0

101     CONTINUE

            READ ( XRDEV, 93000, END=111, IOSTAT=IOS ) LINE
            J = J + 1 
            IF( J .GT. MXPLT )THEN
                WRITE( MESG, 94000 )
     &             'Exceeding MXPLT, currently MXPLT =' , MXPLT
                CALL M3MESG( MESG )
                STOP
            ENDIF


            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94000 )
     &              'Error', IOS, 
     &              'reading Ascii factors file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                ECNT  = ECNT + 1
                IF ( ECNT .GT. 10 ) THEN
                    CALL M3EXIT( 'GRWPMAT', 0, 0, 
     &                           'Too many errors', 2 )
                END IF
                GO TO  101
            END IF 

            IF ( LINE(  1: 2 ) .EQ. ' ' .OR.
     &           LINE(  3: 5 ) .EQ. ' ' .OR.
     &           LINE(  6:20 ) .EQ. ' ' .OR.
     &           LINE( 21:28 ) .EQ. ' ' ) THEN

               WRITE( MESG, 94000 ) 
     &           'Missing ID field in ascii factors file at line' , IREC
               CALL M3MESG( MESG )
               GO TO 101

            END IF                      !  if field missing
            
            JFIP = 1000 * STR2INT( LINE(  1: 2 ) ) + 
     &                    STR2INT( LINE(  3: 5 ) )
            WRITE( FIPS,94020 ) JFIP
            CALL NORMID( PLID, LINE( 6:20 ) )
            SCC         = LINE( 21:28 )
            XRFSCC( J ) = STR2INT( SCC )

            WRITE(XRFPL( J ),94100)FIPS,SCC,PLID

94100       FORMAT(A5,A8,A15)

            INDXP( J )  = J
            K = 29
            L = 36
C
C   assign factor to specific month
C
            DO M = AMON , ZMON
              X = AMAX1( 0.0, STR2REAL( LINE( K:L ) ) )
              RMFAC(M,J) = X
              K = K + 8
              L = L + 8
            ENDDO
              
            GO TO 101                   !  to head of loop reading xref file
                   
111     CONTINUE                        !  exit from loop reading xref file

        NXRFOB = J

        WRITE( MESG,94000 ) 
     &  'Number of valid ASCII factor records:         ', NXRFOB
        CALL M3MESG( MESG )

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                   'Error reading ASCII factors file', 2 )
        ELSE IF ( J .GT. MXPLT ) THEN
            CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                   'Overflow reading ASCII factors file', 2 )
        END IF


C...........   Create fip-scc-plant table:

        CALL SORTIC( NXRFOB, INDXP, XRFPL )             !  sort on plants
        KEY28 =  CMISS3		!  initialize "last fip/scc/plant"
        S     =  0		!  fip/scc/plant count

        DO  122  I = 1, NXRFOB                   ! loop through xref records

            K     = INDXP( I )
            SCR28 = XRFPL( K )
            IF ( SCR28 .NE. KEY28 ) THEN

                S     = S + 1
                KEY28 = SCR28
                IF ( S .LE. MXPLT ) THEN

                    XRFPLS( S ) = SCR28

                END IF                  !  if S in bounds
            END IF                      !  if new key encountered
            
122     CONTINUE                !  end loop initializing fip-plants table

        NPLT = S
        WRITE( MESG,94000 ) 
     &  'Number of factors fips/scc/plt:                ', S
        CALL M3MESG( MESG )

        IF ( NPLT .GT. MXPLT ) THEN
            CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                   'Overflow building fips/scc/plt table', 2 )
        END IF
C
C...... read in netCDF file with growth factors (not MONTHLY)
C
        GNAME = PROMPTMFILE(
     &           'Enter logical name for netCDF factors file',
     &           FSREAD3, 'PNTGRW', 'GRWPMAT' )

C
            IF ( .NOT. READ3( GNAME, 'NOX', ALLAYS3, 0, 0,
     &                        NOXGR ) ) THEN
                CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                       'Could not read NOX from ' // GNAME, 2 )
            END IF

            IF ( .NOT. READ3( GNAME, 'CO', ALLAYS3, 0, 0,
     &                        COGR ) ) THEN
                CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                       'Could not read CO from ' // GNAME, 2 )
            END IF

            IF ( .NOT. READ3( GNAME, 'ROG', ALLAYS3, 0, 0,
     &                        ROGGR ) ) THEN
                CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                       'Could not read ROG from ' // GNAME, 2 )
            END IF


C...........   Read SMOKE::EMS95 point source list file.
C.......   Strategy:  use NOx growth factors if available and
C.......   and override with monthly growth factors from PASCFAC.

        CALL M3MESG( 'Processing SOURCE LIST file' )


        DO M = AMON, ZMON
           DO  505  K = 1, NPSRC
            AMAT( M, K ) = NOXGR( K )         ! no-match control matrix
505        CONTINUE
        ENDDO

        EFLAG = .FALSE.
        IREC  = 0
        FCNT  = 0
        XCNT  = 0

511      CONTINUE		!  head of loop reading pt-source list file

            READ ( SRDEV, 93070, END=522, IOSTAT=IOS ) 
     &          JFIP, JSCC, SMKPL, SMKST, PLID2, STKID, DEVID,DUMNAM

            IREC = IREC + 1          
            PDESC( IREC ) = DUMNAM


            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94000 )
     &              'Error', IOS,
     &              'reading SOURCE LIST file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                ECNT  = ECNT + 1
                IF ( ECNT .GT. 10 ) THEN
                    CALL M3EXIT( 'GRWPMAT', 0, 0, 
     &                           'Too many errors', 2 )
                END IF
                GO TO  511
            END IF 

            IF ( IREC .LE. NPSRC ) THEN

                WRITE( FIPS, 94020 ) JFIP
                WRITE( SCC , 94010 ) JSCC
                CALL NORMID( PLID, PLID2 )  

                KEY28 = FIPS // SCC // PLID
                S     = FINDC( KEY28, NPLT, XRFPLS ) !  find fip/scc/plant match

                IF ( S .GT. 0 ) THEN

                    XCNT = XCNT + 1
                    IF ( XCNT .LE. MXLST ) THEN
                        PLTSRC( XCNT ) = S
                        SRCPLT( XCNT ) = IREC
                    END IF
                    KFIP(S)=JFIP
                    KSCC(S)=JSCC
                    FPLID(S)=PLID2
                    FSTKID(S)=STKID

                ELSE
 
9302               format('SRDEV NOT matched! ',A28,I8)

                END IF                      !  if fip/scc/plant match or not

            END IF                                      !  if irec in bounds

            GO TO 511

522      CONTINUE		!  exit from loop reading pt-source list file
        
        NSRC    = IREC
        NPLTSRC = XCNT

        WRITE( MESG,94000 ) 
     &  'Number of PSRC POINT SOURCES          ', NSRC
        CALL M3MESG( MESG )

        WRITE( MESG,94000 ) 
     &  'Number of PSRC::XREF PLANT matches    ', XCNT
        CALL M3MESG( MESG )


        IF ( EFLAG ) THEN
            CALL M3EXIT( 'GRWPMAT', 0, 0, 
     &                   'Error reading LIST file', 2 )
        ELSE IF ( IREC .GT. NPSRC ) THEN
            CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                   'Overflow reading LIST file', 2 )
        ELSE IF ( IREC .LT. NPSRC ) THEN
            CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                   'Underflow reading LIST file', 2 )
        ELSE IF ( XCNT .GT. MXLST ) THEN
            CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                   'Overflow building PLANT::SOURCE table', 2 )
        END IF                                  !  if errors or over/under-flow
          

C............................................................................
C.......   Construct and write over the control matrix where matches
C.......   were found:

        CALL M3MESG( 'Overwriting GROWTH factors with matches found..' )

        DO  633  M = AMON, ZMON

          DO  611  S = 1, NPLTSRC

              K = PLTSRC( S )
              AMAT( M, SRCPLT( S ) ) = RMFAC( M, K)

611       CONTINUE

633     CONTINUE                !  end loop on months M in the desired season

C
C...........   Read in emissions values.

        DO  88  V = 1, NIPOL

            IF ( .NOT. READ3( ENAME, EINAM( V ), ALLAYS3, 0, 0,
     &                        EMIS( 1,V ) ) ) THEN
                CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                       'Error reading "' //
     &                       EINAM( V )( 1:TRIMLEN( EINAM( V ) ) )//
     &                       '" from file "' //
     &                       ENAME( 1:TRIMLEN( ENAME ) ) // '".', 2 )

            END IF

            VNAME = 'CTLEFF_' // EINAM( V )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3,
     &                        0, 0, CTLEFF( 1,V ) ) ) THEN
                CALL M3EXIT( 'CTLPMAT', 0, 0,
     &               'Error reading CTLEFF from file' // ENAME, 2 )
            END IF

            VNAME = 'RULEFF_' // EINAM( V )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3,
     &                        0, 0, RULEFF( 1,V ) ) ) THEN
                CALL M3EXIT( 'CTLPMAT', 0, 0,
     &               'Error reading RULEFF from file' // ENAME, 2 )
            END IF

            VNAME = 'RULPEN_' // EINAM( V )
            IF ( .NOT. READ3( ENAME, VNAME, ALLAYS3,
     &                        0, 0, RULPEN( 1,V ) ) ) THEN
                CALL M3EXIT( 'CTLPMAT', 0, 0,
     &               'Error reading RULPEN from file' // ENAME, 2 )
            END IF

88      CONTINUE
        IF ( .NOT. READ3( ENAME, 'IFIP', ALLAYS3,0,0, IFIP ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0,
     &                   'Error reading FIP from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'ISCC', ALLAYS3,0,0, ISCC ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0,
     &                   'Error reading ISCC from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'IPLANT', ALLAYS3,0,0, IPLANT ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0,
     &                   'Error reading IPLANT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'ISTACK', ALLAYS3,0,0, ISTACK ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0,
     &                   'Error reading ISTACK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'TZONES', ALLAYS3,0,0, TZONES ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0,
     &                   'Error reading TZONES from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'TPFLAG', ALLAYS3,0,0, TPFLAG ) ) THEN
            CALL M3EXIT( 'TMPPOINT', 0, 0,
     &                   'Error reading TPFLAG from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'INVYR', ALLAYS3,0,0, INVYR ) ) THEN
            CALL M3EXIT( 'CTLPMAT', 0, 0,
     &                  'Error reading INVYR from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'XLOCA', ALLAYS3,0,0, XLOCA ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                   'Error reading XLOCA from file' // ENAME, 2)
        END IF

        IF ( .NOT. READ3( ENAME, 'YLOCA', ALLAYS3,0,0, YLOCA ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                   'Error reading YLOCA from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKHT', ALLAYS3,0,0, STKHT ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                   'Error reading STKHT from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKDM', ALLAYS3,0,0, STKDM ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                   'Error reading STKDM from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKTK', ALLAYS3,0,0, STKTK ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                   'Error reading STKTK from file' // ENAME, 2 )
        END IF

        IF ( .NOT. READ3( ENAME, 'STKVE', ALLAYS3,0,0, STKVE ) ) THEN
            CALL M3EXIT( 'LAYPOINT', 0, 0,
     &                   'Error reading STKVE from file' // ENAME, 2 )
        END IF

        DO K = 1, NPSRC
          INVYR( K ) = NEWYR
        ENDDO
C
C      use default factors for CO and ROG and default/specific factors
C      for NOx

        DO 95 M = AMON, ZMON

          DO 98 K = 1, NPSRC
            DO 97 V = 1, NIPOL
             IF ( V .EQ. 1) EMISV( K, V ) = EMIS( K, V ) * COGR( K )
             IF ( V .EQ. 2) EMISV( K, V ) = EMIS( K, V ) * AMAT( M, K )
             IF ( V .EQ. 3) EMISV( K, V ) = EMIS( K, V ) * ROGGR( K ) 
97          CONTINUE
98        CONTINUE

C............................................................................
C.......   If rflag, write out report:

C   This may need some more work later....
C

        IF ( RFLAG ) THEN


            WRITE( RDEV, * )
            WRITE( RDEV, * ) 'EMISSIONS BY SOURCE and MONTH'
            WRITE( RDEV, * ) 'units are tons/year'
            WRITE( RDEV, * )

             WRITE( RDEV, 95100) MON( M )
95100        FORMAT('for month :', A3)
             WRITE(RDEV, * ) 

             DO  S = 1, NPLTSRC
               K = SRCPLT(S)
               J = PLTSRC( S )
              IF(MOD(S,30).EQ.1)THEN
              WRITE( RDEV, 95105)
95105        FORMAT('  FIP',4x,'SCC',6x,'PLTID',4x,'90EMIS',4x,
     &        '95t2',5x,'95new',4x,'%DIFF',3x,'NAME')
              ENDIF

               ROLDEM = EMIS(K,2) * NOXGR( K )
               IF (ROLDEM .GT. 0.0)THEN
                PERC = 100*((EMISV(K,2) - ROLDEM)/ROLDEM)
               ELSE
                PERC = 0.0
               ENDIF

               WRITE(RDEV, 95106)IFIP(K),ISCC(K),FPLID(J),ISTACK(K),
     &           EMIS(K,2) ,ROLDEM, EMISV(K,2),PERC, PDESC(K)

95106          FORMAT(I5,I9,1x,A15,I8,3F10.2,F7.1,1x,A40)
             ENDDO

            WRITE( RDEV, * )



        END IF                          ! if rflag:  write report-file

          IF ( .NOT. WRITE3( MNAME( M ), ALLVAR3, 0, 0,  IFIP ) ) THEN
            CALL M3EXIT( 'GRWPMAT', 0, 0,
     &                'Error writing output file "' //
     &                 MNAME( M )( 1:TRIMLEN( MNAME( M )) ) // '"',  2 )
          END IF

95      CONTINUE  

        CALL M3EXIT( 'GRWPMAT', 0, 0,
     &             '====  END PROGRAM GRWPMAT  ====', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93070   FORMAT( I5.5, 1X, I8.8, 1X, 2(I6, 1X), A15, 1X, A12 , 1X, A12,
     &   1x,A40 )

93100   FORMAT( /1X, A )

93101   FORMAT( X, A )

93110   FORMAT( 16X,                   12 ( 4X, A4, :, 4X ) )

93120   FORMAT( 1X, I4, ':', 11X,      12 F12.6 )

93121   FORMAT( 1X, I8, 2( ':',  I4 ), 12 F12.0 )

93122   FORMAT( 18X,                   12 F12.6 )

93130   FORMAT( 1X, A16,               12 F12.6 )

93131   FORMAT( 1X, A16,               12 F12.0 )

93140   FORMAT( 1X, 'FIPS ', A5, 5X, 'SCC ', A8, 5X, 'EMS PLT ', A15 )


C...........   Internal bufferring  formats............ 94xxx

94000   FORMAT( 10 ( A, :, I10, :, 2X ) ) 

94010   FORMAT( I8.8 ) 

94020   FORMAT( I5.5 )

94030   FORMAT( I5.5, I8.8, A )

94040   FORMAT( A, :, I10, :, 2X, 5 ( A, :, X, A, :, 2X ) ) 

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

