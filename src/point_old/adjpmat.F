
        PROGRAM ADJPMAT

C***********************************************************************
C  program body      starts at line  189
C  subroutine NORMID starts at line 1179
C
C  DESCRIPTION:
C       Construct monthly point source fuel-adjustment matrices from data 
C       in fuel use data files.  Version aggregates fuel by fip/scc/plant.
C
C  PRECONDITIONS REQUIRED:
C       fuel use file
C       SMOKE PNTS file
C       SMOKE SRCLIST file
C       OIRM Boiler Cross-reference file
C       Fuel use file
C       LEN( FUSTR ) consistent with parameter MXYRS
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C       Models-3 I/O
C       FIND4, PROMPTFFILE, PROMPTMFILE, GETMENU, GETNUM, GETYN, TRIMLEN
C
C***********************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE)
C                Modeling System
C
C File Version @(#)$Id$
C Pathname:    $Source$
C Last updated: $Date$ 
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
C*************************************************************************
C
C  REVISION  HISTORY:
C       Prototype  2/97 by CJC.
C
C*************************************************************************

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

        EXTERNAL        FINDC, FIND1, GETMENU, GETNUM, GETYN,
     &                  PROMPTFFILE, PROMPTMFILE, TRIMLEN, STR2INT


C...........   PARAMETERS and their descriptions:

        INTEGER         MXMON               ! no. of. months of fuel use input
        INTEGER         MXYRS               ! yrs in adj. factor matrix
        INTEGER         MXXRF               ! max XREF ORIS::BLR count
        INTEGER         MXPLT               ! max XREF fip/scc/plant count
        INTEGER         MXSCC               ! max XREF scc count
        INTEGER         MXLST               ! max plant|scc replace-list size

        PARAMETER     ( MXMON =   12, 
     &                  MXYRS =    5 , ! May-Sept, 1990-94  
     &                  MXXRF = 5000 ,
     &                  MXPLT = 2000 ,
     &                  MXSCC =   60 ,
     &                  MXLST = 5000 )

                        

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER        AMON, ZMON      !  first, last month to be averaged
        INTEGER        KMON, LMON      !  first, last months to be normalized
        REAL           DMON, EMON      !  1 / #( months  in avg, norm )
        INTEGER        NYRS            !  number of years to be averaged
        INTEGER        NSRC            !  # recs, pt src files 
        REAL           AMAT  ( NPSRC ) !  adjustment coefficients
        INTEGER        IFIP, ISCC      !  scratch input vars from pt src list 
        INTEGER        SMKPL, SMKST    !    "       "     "       "      "  
        INTEGER        NPLTSRC
        INTEGER        PLTSRC( MXLST ) !  fip/scc/plant record subscript
        INTEGER        SRCPLT( MXLST ) !  matching src  record subscript
        INTEGER        NSCCSRC
        INTEGER        SCCSRC( MXLST ) !  scc           record subscript
        INTEGER        SRCSCC( MXLST ) !  matching src  record subscript
        INTEGER        NSCC            !  size of scc table
        INTEGER        INDXC ( MXXRF ) !  index array for XRFSCC(S)
        INTEGER        XRFSCC( MXXRF ) !  identifier for scc matchup
        INTEGER        SCCLST( MXSCC ) !  identifier for scc matchup
        REAL           CFUSE ( MXMON, MXSCC ) !  sum     by month,scc
        REAL           CFFUSE( MXMON, MXSCC ) !  factors by month,scc
        INTEGER        CNFUSE(        MXSCC ) !  number of fu-records per scc
        INTEGER        CMATCH(        MXSCC ) !  number of src-scc matches
        INTEGER        NPLT            !  size of fip-plant table
        CHARACTER*28   XRFPL ( MXXRF ) !  identifier for fip-scc-plant matchup
        CHARACTER*28   XRFPLS( MXPLT ) !  sorted "
        INTEGER        INDXP ( MXXRF ) !  index array for XRFPL(S)
        INTEGER        NXRFOB          !  size of fip-oris-boiler table
        CHARACTER*24   XRFOB ( MXXRF ) !  Xref-Plant-ORIS-blr matchup
        CHARACTER*24   XRFOBS( MXXRF ) !  sorted "
        INTEGER        INDXO ( MXXRF ) !  index array for XRFOB(S)
        INTEGER        INDXPS( MXXRF ) !  " into plants-table from oris-blr table
        CHARACTER*90   LINE            !  xref file line buffer
        CHARACTER*1404 FUSTR           !  fuel use file string buffer
        REAL           FUSE  ( MXMON, MXYRS ) !  scratch initial copy
        REAL           LFUSE ( MXMON, MXYRS, MXPLT ) !  fuel use
        REAL           FFUSE ( MXMON, MXYRS, MXPLT ) !  log( norm. fueluse )
        REAL           PFUSE ( MXMON,        MXPLT ) !  sum over years for LFUSE
        REAL           SPFUSE( MXMON,        MXPLT ) !  sum-sq over years
        INTEGER        NPFUSE(               MXPLT ) !  # of entries in PFUSE
        REAL           YFUSE ( MXMON, MXYRS ) !  per-year all-plant sum
        REAL           SYFUSE( MXMON, MXYRS ) !  per-year all-plant sum-squares
        INTEGER        NYFUSE(        MXYRS ) !  number of entries in YFUSE
        REAL           MFUSE ( MXMON ) !  all-year all-plant sum
        REAL           SMFUSE( MXMON ) !  all-year all-plant sum-sq
        INTEGER        NMFUSE          !  number of entries in MFUSE
        REAL           TFUSE ( MXMON ) !  all-year all-plant sum
        REAL           TTFUSE( MXMON ) !  full-season sum
        REAL           UFUSE ( MXMON ) !  Unmatched fuel use
        REAL           WFUSE ( MXMON ) !  Matched fuel use
        CHARACTER*5    FIPS            !  FIPS code for state and county
        CHARACTER*8    SCC             !  SCC code
        CHARACTER*12   STKID           !  EMS95 stack ID 
        CHARACTER*15   PLID            !  scratch OTAG plant id on xrf file
        CHARACTER*12   DEVID           !   "       "   device "     "   "
        CHARACTER*5    ORISID          !  scratch doe ORIS id string  
        CHARACTER*6    BLRID           !  scratch doe boiler id string
        CHARACTER*24   KEY24           !  scratch fips//scc//oris
        CHARACTER*28   KEY28, SCR28    !  scratch fips//scc//plant
     
        INTEGER        YR ( MXYRS )    ! inventory yr for adj factors
        CHARACTER*4    MON( MXMON )    ! month for output file logical names
        DATA           MON /
     &          'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',  
     &          'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /

        INTEGER         C, F, I, J, K, L  !  scratch indices
        INTEGER         M, Y, S, V, Y  !  scratch loop variables
        INTEGER         LDEV           !  log-device
        INTEGER         SRDEV          !  sources list file device
        INTEGER         XRDEV          !  x-reference file device
        INTEGER         FUDEV          !  fuel use file device
        INTEGER         RDEV           !  report-file device
        INTEGER         IOS            !  I/O status
        INTEGER         IREC           !  record # from file
        INTEGER         XCNT, FCNT, ECNT, SCNT
        INTEGER         MCNT, NCNT
        LOGICAL         EFLAG          !  flag:  input error
        LOGICAL         FFLAG          !  flag:  bad input fields
        LOGICAL         RFLAG          !  flag:  make report-file

        CHARACTER*16    ENAME          !  point-source input file logical name 
        CHARACTER*16    MNAME( MXMON ) !  ctl matrix  output file logical names
        CHARACTER*16    NAME           !  scratch filename 
        CHARACTER*256   MESG           !  error message buffer

        REAL            DDD, DIV, LGF, X, LRAT, SUM, SIG


C***********************************************************************
C   begin body of program ADJPMAT

        LDEV = INIT3()

        WRITE( *,92000 ) 
     &' ',
     &'Program ADJPMAT takes a SMOKE PNTS inventory file, an ASCII ',
     &'POINT SOURCES LIST file, an ASCII BOILER CROSS_REFERENCE file,',
     &'and an ASCII FUEL USE file, to produce one CONTROL MATRIX file',
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
            WRITE( *,92000 ) ' ', 'Ending program ADJPMAT', ' '
            CALL EXIT( 0 )
        END IF


C....... Get file name; open input pt src, src list, Xref and fuel use files 

        ENAME = PROMPTMFILE(
     &           'Enter logical name for point source LIST file',
     &           FSREAD3, 'PNTS', 'ADJPMAT' )

        SRDEV = PROMPTFFILE( 
     &           'Enter logical name for SRCLIST file',
     &           .TRUE., .TRUE., 'PSRC', 'ADJPMAT' )

        XRDEV = PROMPTFFILE( 
     &           'Enter logical name for XREF file',
     &           .TRUE., .TRUE., 'FREF', 'ADJPMAT' )

        FUDEV = PROMPTFFILE( 
     &           'Enter logical name for FUELUSE file',
     &           .TRUE., .TRUE., 'FUSE', 'ADJPMAT' )

        NYRS    = GETNUM( 1, MXYRS, 5, 
     &            'Enter number of years in this file' )

        YR( 1 ) = GETNUM( 0, 9999, 1990, 
     &            'Enter starting year for this file' )

        DO  11  I = 2, NYRS
            YR( I ) = YR( I-1 ) + 1
11      CONTINUE

        AMON = GETMENU( 12, 5, 
     &         'Enter first month for season', MON )

        ZMON = GETNUM( AMON, 12, MIN( AMON+4, 12 ) ,
     &         'Enter last month for season' )

        KMON = GETNUM( 1, 12, ( AMON + ZMON ) / 2,
     &         'Enter first normalization month' )

        LMON = GETNUM( KMON, 12, KMON,
     &         'Enter last normalization month' )

        DMON = 1.0 / FLOAT( ZMON - AMON + 1 )
        EMON = 1.0 / FLOAT( LMON - KMON + 1 )

C.......   Read in description of point-source emissions file
 
        IF ( .NOT. DESC3( ENAME ) ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0, 
     &          'Could not get description of file ' // ENAME, 3 )
        END IF
 
C.......   Open control matrix file:  reuse much of POINT SOURCE description
 
        NVARS3D = 3
        FTYPE3D = GRDDED3

        DO  22  V = 1, NIPOL
            VNAME3D( V ) = 'CONTROL_' // EINAM( V )
            UNITS3D( V ) = 'n/a'
            VDESC3D( V ) ='Adjustment factor for this source, variable ' 
     &                     // EINAM( V )( 1:TRIMLEN( EINAM( V ) ) )
            VTYPE3D( V ) = M3REAL
22     CONTINUE

        FDESC3D( 1 ) = 'Matrix of adjustment factors, by source'
        DO  33  I = 2, MXDESC3
            FDESC3D( I ) = ' '
33       CONTINUE

C.......   Get output file name; open control file, and report file if any

        DO  44  I = AMON, ZMON

            NAME = 'PAMAT_'//MON( I )
            MESG = 'Enter logical name for CONTROL file for ' // MON(I)
            MNAME( I ) = PROMPTMFILE( MESG, FSUNKN3, NAME, 'ADJPMAT' )

44      CONTINUE

        RDEV  = PROMPTFFILE( '...for REPORT file, or "NONE"',
     &                       .FALSE., .TRUE., 'REPORT', 'ADJPMAT' )
        RFLAG = ( RDEV .GT. 0 )


C...........   Read and process the ORIS::EMS95 XREF file
C.......   Format is as follows
C.......       State FIP        I2      cols   1: 2
C.......       County FIP       I3      cols   3: 5
C.......       OTAG Plant ID    A15     cols   6:20
C.......       Point ID         A15     cols  21:35
C.......       SCC              I8      cols  36:43
C.......       plant name       A40     cols  44:83
C.......       ORIS ID          A5      cols  84:88

        CALL M3MESG( 'Processing ORIS XREF file' )

        ECNT  = 0
        IREC  = 0
        J     = 0

101     CONTINUE

            READ ( XRDEV, 93000, END=111, IOSTAT=IOS ) LINE
            IREC = IREC + 1 

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94000 )
     &              'Error', IOS, 
     &              'reading ORIS XREF file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                ECNT  = ECNT + 1
                IF ( ECNT .GT. 10 ) THEN
                    CALL M3EXIT( 'ADJPMAT', 0, 0, 
     &                           'Too many errors', 2 )
                END IF
                GO TO  101
            END IF 

            IF ( LINE(  1: 2 ) .EQ. ' ' .OR. 
     &           LINE(  3: 5 ) .EQ. ' ' .OR.
     &           LINE(  6:20 ) .EQ. ' ' .OR.
     &           LINE( 36:43 ) .EQ. ' ' .OR.
     &           LINE( 84:88 ) .EQ. ' ' ) THEN

               WRITE( MESG, 94000 ) 
     &             'Missing field in ORIS XREF file at line' , IREC
               CALL M3MESG( MESG )
               GO TO 101

            END IF                      !  if field missing
            
            CALL NORMID( BLRID, LINE( 21:35 ) )
            IF ( BLRID .EQ. ' ' ) THEN
               GO TO 101
            END IF

            CALL NORMID( ORISID, LINE( 84:88 ) )
            IF ( ORISID .EQ. ' ' ) THEN
               GO TO 101
            END IF

            J = J + 1

            IF ( J .LE. MXXRF ) THEN
                INDXP( J )  = J
                INDXO( J )  = J
                INDXC( J )  = J
                IFIP = 1000 * STR2INT( LINE(  1: 2 ) ) + 
     &                        STR2INT( LINE(  3: 5 ) )
                WRITE( FIPS,94020 ) IFIP
                CALL NORMID( PLID, LINE( 6:20 ) )
                SCC         = LINE( 36:43 )
                XRFSCC( J ) = STR2INT( SCC )
                XRFPL( J )  = FIPS // SCC // PLID
                XRFOB( J )  = FIPS // SCC // 
     &                        ORISID( 1:TRIMLEN( ORISID ) ) //
     &                        BLRID                 
            END IF                      !  if J in bounds

            GO TO 101                   !  to head of loop reading xref file
                   
111     CONTINUE                        !  exit from loop reading xref file

        NXRFOB = J

        WRITE( MESG,94000 ) 
     &  'Number of valid XREF records:         ', NXRFOB
        CALL M3MESG( MESG )

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0,
     &                   'Error reading cross reference file', 2 )
        ELSE IF ( J .GT. MXXRF ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0,
     &                   'Overflow reading cross reference file', 2 )
        END IF


C...........   Create SCC-only table, fip-scc-plant tables:

        CALL SORTI1( NXRFOB, INDXC, XRFSCC )
            
        CALL SORTIC( NXRFOB, INDXP, XRFPL )             !  sort on plants

        KEY28 =  CMISS3		!  initialize "last fip/scc/plant"
        L     = -1		!  initialize "last scc"
        S     =  0		!  fip/scc/plant count
        C     =  0		!  scc count

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
            
            K = INDXC( I )
            M = XRFSCC( K )
            IF ( M .NE. L ) THEN
                C           = C + 1
                L           = M
                IF ( C .LE. MXSCC ) THEN
                    SCCLST( C ) = M
                    CNFUSE( C ) = 0
                    CMATCH( C ) = 0
                END IF			!  if C in bounds
            END IF			!  if new SCC encountered

122     CONTINUE                !  end loop initializing fip-plants table

        NSCC = C
        WRITE( MESG,94000 ) 
     &  'Number of XREF   SCCs:                ', C
        CALL M3MESG( MESG )

        NPLT = S
        WRITE( MESG,94000 ) 
     &  'Number of XREF plants:                ', S
        CALL M3MESG( MESG )

        IF ( NSCC .GT. MXSCC ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0,
     &                   'Overflow building XREF SCCs table', 2 )
        ELSE IF ( NPLT .GT. MXPLT ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0,
     &                   'Overflow building XREF PLANTS table', 2 )
        END IF


C...........   fip-oris-boiler tables (with index INDXPS into fip-plant table)

        CALL SORTIC( NXRFOB, INDXO, XRFOB )             !  sort on oris-blr

        DO  133  I = 1, NXRFOB                   ! loop through xref records

            K           = INDXO( I )
            XRFOBS( I ) = XRFOB( K )
            INDXPS( I ) = FINDC( XRFPL( K ), NPLT, XRFPLS )

133     CONTINUE                        !  end loop sorting by fip-oris-boiler


C...........   Read, process, and sort fuel use data file:
C.......   Fuel use records are formatted ( A83, MXMON*NYRS( I22 ) )
C.......   where the A83 field contains FIPS in cols 1:5, 
C.......   ORISID in cols 66:70, boiler id in cols 71:75,
C.......   and SCC in cols 76:83

        CALL M3MESG( 'Processing FUEL USE file' )

        DO  143  S = 1, MXPLT
        DO  142  Y = 1, MXYRS
        DO  141  M = 1, MXMON
            LFUSE( M,Y,S ) = 0.0
141     CONTINUE
142     CONTINUE
143     CONTINUE

        DO  145  C = 1, NSCC
        DO  144  M = 1, MXMON
            CFUSE( M,C ) = 0.0
144     CONTINUE
145     CONTINUE

        DO  146  M = 1, MXMON
            TTFUSE( M ) = 0.0
            TFUSE ( M ) = 0.0
            UFUSE ( M ) = 0.0
            WFUSE ( M ) = 0.0
146     CONTINUE

        EFLAG = .FALSE.
        IREC  = 0
        FCNT  = 0
        MCNT  = 0			!  number of XREF match-failures
        NCNT  = 0			!  number of  SCC match-failures

155     CONTINUE                        !  head of loop reading fuel use file

            READ ( FUDEV, 93000, END=177, IOSTAT=IOS )  FUSTR
            IREC = IREC + 1 

            IF ( IOS .NE. 0 ) THEN

                WRITE( MESG,94000 ) 
     &          'Error', IOS, 'reading FUEL USE at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                ECNT  = ECNT + 1
                IF ( ECNT .GT. 10 ) THEN
                    CALL M3EXIT( 'ADJPMAT', 0, 0, 
     &                           'Too many errors', 2 )
                END IF
               GO TO 155                !  to head of loop reading fuel use file

            ELSE IF ( FUSTR( 1:83 ) .EQ. ' ') THEN

               WRITE( MESG, 94000 ) 
     &              'Missing field in fuel use file at line', IREC
               CALL M3MESG( MESG )
               GO TO 155                !  to head of loop reading fuel use file

            END IF

            FIPS = FUSTR(  1: 5 )
            SCC  = FUSTR( 76:83 )
            CALL NORMID( ORISID , FUSTR( 66:70 ) )
            CALL NORMID( BLRID  , FUSTR( 71:75 ) )

            IF      ( FIPS .EQ. ' ' ) THEN
               WRITE( MESG, 94000 ) 'Missing FIPS at line' , IREC
               CALL M3MESG( MESG )
               GO TO 155                !  to head of loop reading fuel use file
            ELSE IF ( ORISID .EQ. ' ' ) THEN
               WRITE( MESG, 94000 ) 'Missing ORISID at line' , IREC
               CALL M3MESG( MESG )
               GO TO 155                !  to head of loop reading fuel use file
            ELSE IF ( SCC  .EQ. ' ' ) THEN
               WRITE( MESG, 94000 ) 'Missing SCC at line' , IREC
               CALL M3MESG( MESG )
               GO TO 155                !  to head of loop reading fuel use file
            END IF

            DO  162  Y = 1, NYRS
            DO  161  M = MIN( AMON, KMON ), MAX( ZMON, LMON )
                        
                K = 84 + 22 * ( M - 1 + 12 * ( Y - 1 ) )
                L = K  + 21

                X = FLOAT( MAX( 0, STR2INT( FUSTR( K:L ) ) ) )

                FUSE ( M,Y ) = X
                TFUSE( M   ) = TFUSE( M ) + X

161         CONTINUE                ! end  loop on months M in this season
162         CONTINUE                !  end loop on years Y

            KEY24  = FIPS // SCC // 
     &               ORISID( 1:TRIMLEN( ORISID ) ) // BLRID
            J      = FINDC( KEY24, NXRFOB, XRFOBS )

            IF ( J .GT. 0 ) THEN

                FCNT = FCNT + 1
                S    = INDXPS( J )       !  plant index for this oris:boiler

                DO  164  Y = 1, NYRS
                DO  163  M = MIN( AMON, KMON ), MAX( ZMON, LMON )
                    X              = FUSE( M,Y ) 
                    LFUSE( M,Y,S ) = X
                    WFUSE( M     ) = WFUSE( M ) + X
163             CONTINUE                ! end  loop on months M in this season
164             CONTINUE                !  end loop on years Y

            ELSE

                DO  166  Y = 1, NYRS
                DO  165  M = MIN( AMON, KMON ), MAX( ZMON, LMON )
                    UFUSE( M     ) = UFUSE( M ) + FUSE( M,Y )
165             CONTINUE                ! end  loop on months M in this season
166             CONTINUE                !  end loop on years Y

                MCNT = MCNT + 1
                WRITE( MESG, 94040 ) 
     &              'No match at line', IREC,
     &              'fips', FIPS, 'oris', ORISID, 'boiler', BLRID
                CALL M3MESG( MESG )

            END IF                      !  if ORIS:BOILER found
            
            ISCC = STR2INT( SCC )
            C    = FIND1( ISCC, NSCC, SCCLST )
            IF ( C .GT. 0 ) THEN

                CNFUSE( C ) = CNFUSE( C ) + 1
                
                DO  168  Y = 1, NYRS
                DO  167  M = MIN( AMON, KMON ), MAX( ZMON, LMON )
                    CFUSE( M,C ) = CFUSE( M,C ) + FUSE( M,Y )
167             CONTINUE                ! end  loop on months M in this season
168             CONTINUE                !  end loop on years Y

            ELSE
                NCNT = NCNT + 1
            END IF			!  if SCC found

            GO TO 155                   !  to head of loop reading fuel-use file

177     CONTINUE                        !  end of loop reading fuel-use

        WRITE( MESG,94000 ) 
     &  'Number of FUSE::XREF PLANT    matches:', FCNT
        CALL M3MESG( MESG )            

        WRITE( MESG,94000 ) 
     &  'Number of FUSE::XREF PLANT mismatches:', MCNT
        CALL M3MESG( MESG )

        WRITE( MESG,94000 ) 
     &  'Number of FUSE::XREF  SCC  mismatches:', NCNT
        CALL M3MESG( MESG )

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0,
     &              'Error reading fuel use datafile', 2 )
        END IF                                  !  if eflag or overflow


C...........   Collapse the SCC-list (remove unused entries) and 
C...........   normalize each month's by-SCC factors

        K = 0
        DO  188  C = 1, NSCC

            SUM = 0.0
            DO  181  M = KMON, LMON
                SUM = SUM + CFUSE( M,C )
181         CONTINUE

            IF ( SUM .GT. 0.0 ) THEN

                K   = K + 1
                SUM = 1.0 / ( EMON * SUM )
                DO  182  M = AMON, ZMON
                    CFFUSE( M,K ) = SUM * CFUSE( M,C )
                    CFUSE ( M,K ) =       CFUSE( M,C )
182             CONTINUE
                SCCLST( K ) = SCCLST( C )
                CNFUSE( K ) = CNFUSE( C )

            ELSE				!  no fuel use for this SCC:

                WRITE( MESG , 94000 ) 'Unused SCC', SCCLST( M )
                CALL M3MESG( MESG )

            END IF			!  if sum>0 for this SCC, or not
            
188     CONTINUE			!  end loop on XREF SCC's C
        
        NSCC = K
        WRITE( MESG,94000 ) 
     &  'Number of FUSE::SCC matched SCCs:     ', K
        CALL M3MESG( MESG )


C...........   Accumulate sums for statistics:

        NMFUSE = 0
        DO  191  J = 1, MXMON
            MFUSE ( J ) = 0.0
            SMFUSE( J ) = 0.0
191     CONTINUE

        DO  193  J = 1, MXPLT
            NPFUSE( J ) = 0
            DO  192  I = 1, MXMON
                PFUSE ( I,J ) = 0.0
                SPFUSE( I,J ) = 0.0
192         CONTINUE
193     CONTINUE

        DO  195  J = 1, MXYRS
            NYFUSE( J ) = 0
            DO  194  I = 1, MXMON
                YFUSE ( I,J ) = 0.0
                SYFUSE( I,J ) = 0.0
194         CONTINUE
195     CONTINUE

        DO  244  F = 1, NPLT
        DO  243  Y = 1, NYRS

            FFLAG = .FALSE.

            DO  202  M = MIN( AMON, KMON ), MAX( ZMON, LMON )
                X = LFUSE( M,Y,F )
                IF ( X .LE. 0.0 ) THEN
                    FFLAG = .TRUE.
                END IF
202         CONTINUE

            IF ( .NOT. FFLAG ) THEN   !  no bad/zero entries for this src

                LGF = 0.0

                DO  222   M = KMON, LMON
                    LGF = LGF + LFUSE( M,Y,F )
222             CONTINUE        ! end loop on normalization season

                LGF = LOG( EMON * LGF ) !  log( norm-season mean fuel use )

                DO  233   M = AMON, ZMON

                    X               = LFUSE( M,Y,F )
                    LRAT            = LOG( X ) - LGF

                    TTFUSE( M     ) = TTFUSE( M     ) + X
                    MFUSE ( M     ) = MFUSE ( M     ) + LRAT
                    SMFUSE( M     ) = SMFUSE( M     ) + LRAT*LRAT

                    YFUSE ( M,Y   ) = YFUSE ( M,Y   ) + LRAT
                    SYFUSE( M,Y   ) = SYFUSE( M,Y   ) + LRAT*LRAT

                    PFUSE ( M,  F ) = PFUSE ( M,  F ) + LRAT
                    SPFUSE( M,  F ) = SPFUSE( M,  F ) + LRAT*LRAT

                    FFUSE ( M,Y,F ) = EXP( LRAT )

233             CONTINUE        !  end loop on months M in projected season

                NPFUSE( F ) = NPFUSE( F ) + 1
                NYFUSE( Y ) = NYFUSE( Y ) + 1
                NMFUSE      = NMFUSE      + 1

            END IF                      	!  if bad/zero entries or not

243     CONTINUE                        	!  end loop on  years Y
244     CONTINUE                        	!  end loop on plants F


C.......   Sorting and Statistical analysis:

        DO  313  J = 1,    NPLT   !  all-years, per-source means and deviations
        
            K = NPFUSE( J )

            IF ( K .GT. 0 ) THEN

                DIV = 1.0 / FLOAT( K )
            
                DO  311  M = AMON, ZMON
                    SUM = DIV * PFUSE ( M,J )
                    SIG = DIV * SPFUSE( M,J ) - SUM * SUM
                    PFUSE ( M,J ) = EXP( SUM )
                    SPFUSE( M,J ) = EXP( 0.5 * SIG )
311             CONTINUE                        !  end loop on  months  M

            ELSE                                !  reset pfuse for AMAT

                DO  312  M = AMON, ZMON
                    PFUSE ( M,J ) = 1.0
                    SPFUSE( M,J ) = 0.0
312             CONTINUE                        !  end loop on  months  M

            END IF                      !  if num of entries > 0 for source S

313     CONTINUE                        !  end loop on sources  S

        DO  322  Y = 1,    NYRS         !  per-year, all-source means and devs

            K = NYFUSE( Y )

            IF ( K .GT. 0 ) THEN

                DIV = 1.0 / FLOAT( K )

                DO  321  M = AMON, ZMON
                        SUM = DIV * YFUSE ( M,Y )
                        SIG = DIV * SYFUSE( M,Y ) - SUM * SUM
                        YFUSE ( M,Y ) = EXP( SUM )
                        SYFUSE( M,Y ) = EXP( 0.5 * SIG )
321             CONTINUE                        !  end loop on months  M

            END IF                      !  if num of entries > 0 for year Y

322     CONTINUE                        !  end loop on  years  Y

        IF ( NMFUSE .GT. 0 ) THEN

            DIV = 1.0 / FLOAT( NMFUSE )
        
            DO  333  M = AMON, ZMON     !  all source, all-year means and devs
                    SUM = DIV * MFUSE ( M )
                    SIG = DIV * SMFUSE( M ) - SUM * SUM
                    MFUSE ( M ) = EXP( SUM )
                    SMFUSE( M ) = EXP( 0.5 * SIG )
333         CONTINUE                        !  end loop on months  M

        END IF                              !  if number of entries > 0


C...........   Read SMOKE::EMS95 point source list file.
C.......   Create xref-ID, fueluse-ID strings and look up sources.
C.......   Strategy:  use xref file lookup to get ORIS ID and boiler ID.
C.......   Then use these to construct a key into the fuel-use ID table.

        CALL M3MESG( 'Processing SOURCE LIST file' )

        DO  505  K = 1, NPSRC
            AMAT( K ) = 1.0                    ! no-match control matrix
505     CONTINUE

        EFLAG = .FALSE.
        IREC  = 0
        FCNT  = 0
        XCNT  = 0
        SCNT  = 0

511      CONTINUE		!  head of loop reading pt-source list file

            READ ( SRDEV, 93070, END=522, IOSTAT=IOS ) 
     &          IFIP, ISCC, SMKPL, SMKST, PLID, STKID, DEVID
            IREC = IREC + 1          

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94000 )
     &              'Error', IOS,
     &              'reading SOURCE LIST file at line', IREC
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                ECNT  = ECNT + 1
                IF ( ECNT .GT. 10 ) THEN
                    CALL M3EXIT( 'ADJPMAT', 0, 0, 
     &                           'Too many errors', 2 )
                END IF
                GO TO  511
            END IF 

            IF ( IREC .LE. NPSRC ) THEN

                WRITE( FIPS, 94020 ) IFIP
                WRITE( SCC , 94010 ) ISCC

                KEY28 = FIPS // SCC // PLID
                S     = FINDC( KEY28, NPLT, XRFPLS ) !  find fip/scc/plant match

                IF ( S .GT. 0 ) THEN

                    XCNT = XCNT + 1
                    IF ( XCNT .LE. MXLST ) THEN
                        PLTSRC( XCNT ) = S
                        SRCPLT( XCNT ) = IREC
                    END IF

                ELSE 	!  no plant match; try to find scc match

                    C = FIND1( ISCC, NSCC, SCCLST )

                    IF ( C .GT. 0 ) THEN

                        SCNT        = SCNT + 1
                        CMATCH( C ) = CMATCH( C ) + 1
                        IF ( SCNT .LT. MXLST ) THEN
                            SCCSRC( SCNT ) = C
                            SRCSCC( SCNT ) = IREC
                        END IF

                    END IF                                  !  if scc match

                END IF                      !  if fip/scc/plant match or not

            END IF                                      !  if irec in bounds

            GO TO 511

522      CONTINUE		!  exit from loop reading pt-source list file
        
        NSRC    = IREC
        NPLTSRC = XCNT
        NSCCSRC = SCNT

        WRITE( MESG,94000 ) 
     &  'Number of PSRC POINT SOURCES          ', NSRC
        CALL M3MESG( MESG )

        WRITE( MESG,94000 ) 
     &  'Number of PSRC::XREF PLANT matches    ', XCNT
        CALL M3MESG( MESG )

        WRITE( MESG,94000 ) 
     &  'Number of PSRC::XREF  SCC  matches    ', SCNT
        CALL M3MESG( MESG )

        IF ( EFLAG ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0, 
     &                   'Error reading LIST file', 2 )
        ELSE IF ( IREC .GT. NPSRC ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0,
     &                   'Overflow reading LIST file', 2 )
        ELSE IF ( IREC .LT. NPSRC ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0,
     &                   'Underflow reading LIST file', 2 )
        ELSE IF ( XCNT .GT. MXLST ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0,
     &                   'Overflow building PLANT::SOURCE table', 2 )
        ELSE IF ( SCNT .GT. MXLST ) THEN
            CALL M3EXIT( 'ADJPMAT', 0, 0,
     &                   'Overflow building SCC::SOURCE table', 2 )
        END IF                                  !  if errors or over/under-flow
          

C............................................................................
C.......   Construct and write out the control matrix:

        CALL M3MESG( 'Writing out CONTROL MATRIX...' )

        DO  633  M = AMON, ZMON

          DO  611  S = 1, NPLTSRC
              AMAT( SRCPLT( S ) ) = PFUSE( M, PLTSRC( S ) )
611       CONTINUE

          DO  612  S = 1, NPLTSRC
              AMAT( SRCSCC( S ) ) = CFFUSE( M, SCCSRC( S ) )
612       CONTINUE

          NAME = MNAME( M )
          DO  622  V = 1, NIPOL

              IF ( .NOT.WRITE3( NAME, VNAME3D( V ), 0, 0, AMAT ) ) THEN
                  CALL M3EXIT( 'ADJPMAT', 0, 0, 
     &                       'Error writing CONTROL MATRIX file.', 2 )
              END IF

622       CONTINUE              !  end loop writing AMAT by pollutant V

633     CONTINUE                !  end loop on months M in the desired season


C............................................................................
C.......   If rflag, write out report:

        IF ( RFLAG ) THEN

            WRITE( RDEV, * )
            WRITE( RDEV, 93110 )  ( MON( M ), M = AMON, ZMON )

            WRITE( RDEV, 93100 ) 'RAW TOTALS (all plants, all years)'
            WRITE( RDEV, * )

            WRITE( RDEV, 93131 ) 
     &          'XREF-Match  tot',
     &          ( WFUSE( M ) , M = AMON, ZMON )

            WRITE( RDEV, 93131 ) 
     &          'Full-season tot',
     &          ( TTFUSE( M ) , M = AMON, ZMON )

            WRITE( RDEV, 93131 ) 
     &          'Unmatched   tot',
     &          ( UFUSE( M ) , M = AMON, ZMON )

            WRITE( RDEV, 93131 ) 
     &          'Grand Total    ',
     &          ( TFUSE( M ) , M = AMON, ZMON )

            WRITE( RDEV, 93130 ) 
     &          'Frac full-season',
     &          ( TTFUSE( M ) / TFUSE( M ), M = AMON, ZMON )

            WRITE( RDEV, 93130 ) 
     &          'Frac Matched    ',
     &          ( 1.0 - UFUSE( M ) / TFUSE( M ), M = AMON, ZMON )

 
            DIV  = 0.0
            DDD  = 0.0
            SUM  = 0.0
            DO  707  M = KMON, LMON

                X = TFUSE( M )
                IF ( X .GT. 0.0 ) THEN
                    DDD  = DDD  + X
                END IF

                X = WFUSE( M )
                IF ( X .GT. 0.0 ) THEN
                    DIV  = DIV  + X
                END IF

                X = TTFUSE( M )
                IF ( X .GT. 0.0 ) THEN
                    SUM  = SUM  + X
                END IF

707         CONTINUE
            DIV = 1.0 / ( EMON * DIV )          ! 1/(avg monthly fuel-use)
            DDD = 1.0 / ( EMON * DDD )          ! 1/(avg monthly fuel-use)
            SUM = 1.0 / ( EMON * SUM )          ! 1/(avg monthly fuel-use)

            WRITE( RDEV, * )
            WRITE( RDEV, 93100 ) 'RATIOS (all plants, all years)'
            WRITE( RDEV, * )

            WRITE( RDEV, 93130 ) 
     &          'Match   total ',
     &          ( SUM * TTFUSE( M ) , M = AMON, ZMON )

            WRITE( RDEV, 93130 ) 
     &          'Full-Seas tot ',
     &          ( DIV * WFUSE( M ) , M = AMON, ZMON )

            WRITE( RDEV, 93130 ) 
     &          'File    total ',
     &          ( DDD * TFUSE( M ),  M = AMON, ZMON )

            WRITE( RDEV, * )
            WRITE( RDEV, 93100 ) 
     &      'By-year source-averaged GEOMETRIC MEAN ADJ. FACTORS'
            WRITE( RDEV, * )

            DO  711  Y = 1, NYRS

                WRITE( RDEV, 93120 ) 
     &                YR( Y ),
     &              ( YFUSE ( M, Y ), M = AMON, ZMON )

711         CONTINUE			!  end loop reporting factors by year

            WRITE( RDEV, 93130 )  
     &          'All years      ', 
     &          ( MFUSE( M ) , M = AMON, ZMON )

            WRITE( RDEV, 93100 ) 
     &      'By-year source-averaged GEOMETRIC DEVIATION'

            DO  722  Y = 1, NYRS

                WRITE( RDEV, 93120 ) 
     &                YR( Y ),
     &              ( SYFUSE( M, Y ), M = AMON, ZMON )

722         CONTINUE			!  end loop reporting devs by year

            WRITE( RDEV, 93130 ) 
     &          'All years      ', 
     &          ( SMFUSE( M ) , M = AMON, ZMON )

            WRITE( RDEV, * ) 
            WRITE( RDEV, 93100 ) 'BY-SCC FUEL USE totals'
            WRITE( RDEV, 93101 ) 'and composite ADJUSTMENT FACTORS'
            WRITE( RDEV, * )

            DO  733  C = 1, NSCC

                WRITE( RDEV, 93121 )
     &              SCCLST( C ),
     &              CNFUSE( C ),
     &              CMATCH( C ),
     &            ( CFUSE( M,C ) , M = AMON, ZMON )

                WRITE( RDEV, 93122 )
     &            ( CFFUSE( M,C ), M = AMON, ZMON )


                WRITE( RDEV, * )

733         CONTINUE			!  end loop reporting totals by SCC

            WRITE( RDEV, * ) 
            WRITE( RDEV, 93100 ) 'BY-SCC composite ADJUSTMENT FACTORS'
            WRITE( RDEV, * )


            DO  744  C = 1, NSCC

744         CONTINUE			!  end loop reporting factors by SCC

            WRITE( RDEV, * ) 
            WRITE( RDEV, 93100 ) 'Adjustment factors BY PLANT AND YEAR'
            WRITE( RDEV, * ) 
            WRITE( RDEV, 93110 ) ( MON( M ), M = AMON, ZMON )
            WRITE( RDEV, * ) 

            DO  799  S = 1, NPLT	!  report factors and devs by plant

                WRITE( RDEV, 93140 ) 
     &                  XRFPLS( S )(  1: 5 ),	!  fips code
     &                  XRFPLS( S )(  6:13 ),	!  scc  code
     &                  XRFPLS( S )( 14:28 )	!  EMS-95 plant id
                            
                FCNT = 0
                DO  788  Y = 1, NYRS

                    IF ( FFUSE( AMON,Y,S ) .GT. 0.0 ) THEN
                        FCNT = FCNT + 1
                        WRITE( RDEV, 93120 ) 
     &                    YR( Y ),  
     &                  ( FFUSE( M,Y,S ), M = AMON, ZMON )
                    END IF                      !  if valid year for this plant
     
788             CONTINUE                        !  end loop on years Y

                IF ( FCNT .GT. 1 ) THEN

                    WRITE( RDEV, 93130 ) 
     &                  'Geometric mean  ',
     &                  ( PFUSE ( M,S ),  M = AMON, ZMON )

                    WRITE( RDEV, 93130 ) 
     &                  'Geometric dev   ',
     &                  ( SPFUSE( M,S ),  M = AMON, ZMON )

                END IF                  !  if any year valid for this plant

                WRITE( RDEV, * ) 

799         CONTINUE            !  end loop reporting factors and devs by plant

            WRITE( RDEV, * )

        END IF                          ! if rflag:  write report-file


        CALL M3EXIT( 'ADJPMAT', 0, 0,
     &             '====  END PROGRAM ADJPMAT  ====', 0 )


C******************  FORMAT  STATEMENTS   ******************************

C...........   Informational (LOG) message formats... 92xxx

92000   FORMAT( 5X, A )


C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93070   FORMAT( I5.5, 1X, I8.8, 1X, 2(I6, 1X), A15, 1X, A12 , 1X, A12 )

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

