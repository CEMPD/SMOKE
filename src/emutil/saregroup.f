
      PROGRAM REGROUP

C***********************************************************************
C  program REGROUP body starts at line 118
C
C  DESCRIPTION:
C       The purpose of this program is to merge multiple source-groups
C       from program SMKMERGE, with optional renaming, and then output
C       the set of new groups.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Initial version 10/29/2019 by Carlie J. Coats, Jr., UNC IE.
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2019, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C****************************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER, EXTERNAL :: GETFLINE
        LOGICAL, EXTERNAL :: BLKORCMT

C.........  LOCAL PARAMETERS and their descriptions:

        INTEGER,       PARAMETER :: SRTLEN3  = FIPLEN3+SCCLEN3+PLTLEN3+CHRLEN3+24

        CHARACTER(1) , PARAMETER :: BLANK    = ' '
        CHARACTER(16), PARAMETER :: PROGNAME = 'SAREGROUP' ! program name
        CHARACTER(64), PARAMETER ::
     &  BAR = '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
        CHARACTER(50), PARAMETER ::
     &  CVSW = '$Name SMOKEv4.7_Nov2019$' ! CVS release tag

C...........   LOCAL VARIABLES and their descriptions:

        INTEGER             LDEV, GDEV, ISTAT
        INTEGER             GRPDEV              !! unit number for REGROUP_GROUPS file
        INTEGER             INRPRT              !! unit number for IN_REPORT file
        INTEGER             OUTRPT              !! unit number for OUT_REPORT file
        INTEGER             F, K, L, M, N, V, LO, HI
        INTEGER             IGRP, ICOL, IROW, GDEX
        INTEGER             SDATE, STIME, TSTEP, NRECS, JDATE, JTIME
        LOGICAL             EFLAG

        INTEGER             MXGID        !! number of lines  in REGROUP_GROUPS
        INTEGER             NGRPS        !! number of grpIDs in REGROUP_GROUPS
        INTEGER             NRPTLINE     !! number of lines in IN_REPORT
        INTEGER             NSRCALL      !! number of srcs  in IN_REPORT
        INTEGER             NVECSSRC     !! number of srcs  in  IN_SRCS,  IN_EMIS
        INTEGER             NVECOUT      !! number of srcs  in OUT_SRCS, OUT_EMIS

        INTEGER           , ALLOCATABLE :: GID_IN( : )      !!  group-ID rename table
        INTEGER           , ALLOCATABLE :: GIDOUT( : )      !!  from REGROUP_GROUPS

        INTEGER           , ALLOCATABLE :: COLSRC( : )      !!  COL    from IN_SRCS
        INTEGER           , ALLOCATABLE :: ROWSRC( : )      !!  ROW    from IN_SRCS
        INTEGER           , ALLOCATABLE :: GRPSIN( : )      !!  IGROUP from IN_SRCS
        INTEGER           , ALLOCATABLE :: GRPOUT( : )      !!  IGROUP for OUT_SRCS
        INTEGER           , ALLOCATABLE :: INDXIN( : )      !!  sort-index
        INTEGER           , ALLOCATABLE :: VECCNT( : )      !!  counts  for aggregation-matrix
        INTEGER           , ALLOCATABLE :: VECDEX( : )      !!  indices for aggregation-matrix
        
        CHARACTER(NAMLEN3)  SECTOR
        CHARACTER(1)        CATEGORY
        CHARACTER(256)      MESG
        CHARACTER(256)      LINE
        CHARACTER(1024)     CBUF
        CHARACTER(32)       SEGMENT( 64 )

        INTEGER           , ALLOCATABLE :: IGROUP( : )      !! (NRPTLINE):  input-groups group-id
        INTEGER           , ALLOCATABLE :: SGROUP( : )      !! (NRPTLINE): output-groups group-id
        CHARACTER(FIPLEN3), ALLOCATABLE :: CFIPIN( : )      !! (NRPTLINE):  input FIP list
        CHARACTER(SCCLEN3), ALLOCATABLE :: CSCCIN( : )      !! (NRPTLINE):  input SCC list
        CHARACTER(PLTLEN3), ALLOCATABLE :: CPLTIN( : )      !! (NRPTLINE):  input plant list
        CHARACTER(CHRLEN3), ALLOCATABLE :: CPNTIN( : )      !! (NRPTLINE):  input stack list

        INTEGER             NVARSSRC                !! number of data variables for IN_SRCS
        CHARACTER(NAMLEN3)  VARSSRC( MXVARS3 )
        INTEGER             TYPESRC( MXVARS3 )

        INTEGER             NVAREMIS                !! number of data variables for IN_EMIS
        CHARACTER(NAMLEN3)  VAREMIS( MXVARS3 )
        INTEGER             TYPEMIS( MXVARS3 )

C***********************************************************************
C   begin body of program REGROUP

        EFLAG = .FALSE.
        LDEV  = INIT3()
        WRITE( LDEV, '(5X, A )' ) BAR,
     &  'SMOKE ---------------',
     &  'Program ' // PROGNAME,
     &  'Copyright (c) 2019 Environmental Modeling for Policy Development',
     &  'All rights reserved',
     &  '',
     &  'Online documentation available at:',
     &  '    https://www.cmascenter.org/smoke/',
     &  '',
     &  'Program REGROUP to merge/rename multiple source-groups from ,',
     &  'program SMKMERGE, and then output the set of new groups.',
     &  '',
     &  'PRECONDITIONS REQUIRED:',
     &  '   setenv  MRG_SOURCE      <A|P|B|M> (one letter only)',
     &  '   setenv  MRG_SECTOR      <sector',
     &  '',
     &  '   setenv  REGROUP_GROUPS  <ASCII re-group definitions, with',
     &  '   with sets of lines formatted',
     &  '',
     &  '       <input-file group-ID>  <output group-ID>',
     &  '',
     &  '   setenv  IN_REPORT       <ASCII SRCGRP_REPORT for output groups>',
     &  '   setenv  IN_SRCS         <I/O API output source-attribute file>',
     &  '   setenv  IN_EMIS         <I/O API output timestepped data file>',
     &  '',
     &  '   setenv  OUT_REPORT      <ASCII SRCGRP_REPORT for output groups>',
     &  '   setenv  OUT_SRCS        <I/O API output source-attribute file>',
     &  '   setenv  OUT_EMIS        <I/O API output timestepped data file>',
     &  ''

        CALL ENVSTR( 'MRG_SOURCE', 'Source category:  A, B, M, or P>', 'A', CATEGORY, ISTAT  )
        IF ( ISTAT .GT. 0 ) THEN
            CALL M3MESG( 'Bad environment variable "MRG_SOURCE"' )
            EFLAG = .TRUE.
        ELSE IF ( INDEX( 'ABMPabmp', CATEGORY ) .EQ. 0 ) THEN
            CALL M3MESG( 'Bad MRG_SOURCE:  not A,B,M,or P' )
            EFLAG = .TRUE.
        END IF
        CALL UPCASE( CATEGORY )


        CALL ENVSTR( 'MRG_SECTOR', 'Sector', BLANK, SECTOR, ISTAT  )
        IF ( ISTAT .GT. 0 ) THEN
            CALL M3MESG( 'Bad environment variable "MRG_SECTOR"' )
            EFLAG = .TRUE.
        ELSE IF ( SECTOR .EQ. BLANK ) THEN
            CALL M3MESG( 'Bad/missing MRG_SECTOR' )
            EFLAG = .TRUE.
        END IF
        CALL UPCASE( SECTOR )

        GRPDEV = GETEFILE( 'REGROUP_GROUPS', .TRUE., .TRUE., PROGNAME )
        IF ( GRPDEV .LT. 0 ) THEN
            CALL M3MESG( 'Could not open "REGROUP_GROUPS"' )
            EFLAG = .TRUE.
        ELSE
            MXGID = GETFLINE( GRPDEV, 'ASCII re-group definitions' )
        END IF

        IF ( EFLAG ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Environment error(s) ', 2 )
        END IF


C.........  Process REGROUP_GROUPS:  First read and parse header-line
C.........  then read group-renaming data:

        ALLOCATE( GID_IN( MXGID ),
     &            GIDOUT( MXGID ), STAT = ISTAT )
        CALL CHECKMEM( ISTAT, 'GID_IN,GIDOUT', PROGNAME )

        READ( GRPDEV, '(A)', IOSTAT=ISTAT, END=199 ) CBUF       !!  process header-line
        IF ( ISTAT .NE. 0 ) THEN
            WRITE( MESG, '( A, I9, 2X, A  )' )
     &              'ERROR=', ISTAT, 'reading "REGROUP_GROUPS" at line 1'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF

        CALL PARSLINE( CBUF, 64, SEGMENT )
        GDEX = IMISS3
        DO N = 1, 64
            CALL UPCASE( SEGMENT( N ) )
            IF ( SEGMENT( N ) .EQ. SECTOR ) THEN
                GDEX = N
                EXIT
            END IF
        END DO

        IF ( GDEX .EQ. IMISS3 ) THEN
            MESG = 'SECTOR "' // TRIM( SECTOR) // '" not found in REGROUP_GROUPS header'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )            
        END IF

        N = 0
        DO L = 2, MXGID     !!  process data lines

            READ( GRPDEV, '(A)', IOSTAT=ISTAT, END=199 ) LINE
            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( 2( A, I9, :, 2X ) )' )
     &              'ERROR=', ISTAT, 'reading "REGROUP_GROUPS" at line', L
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                CYCLE
            ELSE IF (  BLKORCMT( LINE ) ) THEN
                CYCLE
            END IF
            
            CALL PARSLINE( LINE, 64, SEGMENT )

            N = N + 1
            GIDOUT( N ) = STR2INT( SEGMENT( GDEX ) )
            GID_IN( N ) = STR2INT( SEGMENT( 1 ) )
            
            IF ( GIDOUT( N ) .EQ. IMISS3 ) THEN
                WRITE( MESG, '( A, I9, :, 2X )' )
     &              'Bad GIDOUT in "REGROUP_GROUPS" at line', L
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF
            
            IF ( GID_IN( N ) .EQ. IMISS3 ) THEN
                WRITE( MESG, '( A, I9 )' )
     &              'Bad GID_IN in "REGROUP_GROUPS" at line', L
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF

        END DO

199     CONTINUE

        NGRPS = N

        IF ( EFLAG ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error(s) processing "REGROUP_GROUPS"', 2 )
        END IF


C.........  Process input REPORT file:

        INRPRT = GETEFILE( 'IN_REPORT', .TRUE., .TRUE., PROGNAME )
        IF ( INRPRT .LT. 0 ) THEN
            CALL M3MESG( 'Could not open "IN_REPORT"' )
            EFLAG = .TRUE.
        ELSE
            NRPTLINE = GETFLINE( INRPRT, 'IN_REPORT' )
        END IF

        ALLOCATE( CFIPIN( NRPTLINE ),
     &            CSCCIN( NRPTLINE ),
     &            CPLTIN( NRPTLINE ),
     &            CPNTIN( NRPTLINE ),
     &            IGROUP( NRPTLINE ),
     &            SGROUP( NRPTLINE ), STAT = ISTAT )
        CALL CHECKMEM( ISTAT, 'CFIPIN...SGROUP', PROGNAME )
        SGROUP = 0

        DO L = 1, NRPTLINE

            READ( INRPRT, '(A)', IOSTAT=ISTAT,END=399 ) LINE
            IF ( ISTAT .NE. 0 ) THEN
                WRITE( MESG, '( A, I9)' )
     &              'ERROR=', ISTAT, 'reading "IN_REPORT" at line', L
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
                CYCLE
            END IF

            IF ( BLKORCMT( LINE ) )  CYCLE

            SELECT CASE ( CATEGORY )
                CASE( 'A' )
                    READ( LINE, '(I8,1X,A12,1X,A20,1X,I8)')
     &                  M, CFIPIN(M), CSCCIN(M), IGROUP(M)
                    CPLTIN( M ) = BLANK
                    CPNTIN( M ) = BLANK
                CASE( 'B' )
                    READ( LINE, '(A12,1X,I8)')
     &                  CFIPIN(M), IGROUP(M)
                        M = L
                    CSCCIN( M ) = BLANK
                    CPLTIN( M ) = BLANK
                    CPNTIN( M ) = BLANK
                CASE( 'M' )
                    READ( LINE, '(I8,1X,A12,1X,A20,1X,I8)')
     &                  M, CFIPIN(M), CSCCIN(M), IGROUP(M)
                    CPLTIN( M ) = BLANK
                    CPNTIN( M ) = BLANK
                CASE( 'P' )
                    READ( LINE, '(I8,1X,A12,1X,A20,1X,A20,1X,A20,1X,I8)')
     &                  M, CFIPIN(M), CSCCIN(M), CPLTIN(M), CPNTIN(M), IGROUP(M)
            END SELECT

            SGROUP( M ) = IGROUP( M )
            DO K = 1, NGRPS
                IF ( IGROUP( M ) .EQ. GID_IN( K ) ) THEN
                    SGROUP( M ) = GIDOUT( K )
                    EXIT
                END IF
            END DO

        END DO      !!  end loop on lines for this report-file

399     CONTINUE
        NSRCALL = M

        IF ( EFLAG ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error reading "IN_REPORT"', 2 )
        END IF
        WRITE( MESG, '( A, I10 )' ) 'File IN_REPORT:  Number of sources', NSRCALL
        CALL M3MESG( MESG )


CC.........  Check consistency for regrouping:
C
C        DO N = 1, NGRPS
C            IGRP = GIDOUT( N )
C            IF ( INDEXINT1( IGRP, NGRPS,   GID_IN ) .GT. 0 )  CYCLE
C            IF ( INDEXINT1( IGRP, NSRCALL, IGROUP ) .GT. 0 )  THEN
C                EFLAG = .TRUE.
C                WRITE( MESG, '(A,1X,I8)' ) 'Invalid output group ID ', IGRP
C                CALL M3MESG( MESG )
C            END IF
C        END DO
C
C        IF ( EFLAG ) THEN
C            CALL M3EXIT( PROGNAME, 0,0, 'Invalid output group ID(s)', 2 )
C        END IF


C.........  Write the OUT_REPORT

        CALL WRRPT_GRPS()


C.........  Open and check the IN_SRCS file

        IF ( .NOT.OPEN3( 'IN_SRCS', FSREAD3, PROGNAME ) ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error opening "IN_SRCS"', 2 )
        ELSE IF ( .NOT.DESC3( 'IN_SRCS' ) ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error in  DESC3( IN_SRCS )', 2 )
        ELSE IF ( NLAYS3D .NE. 1 .OR. NCOLS3D .NE. 1 ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Unexpected COL, LAY dims for "IN_SRCS"', 2 )
        END IF

        SDATE    = SDATE3D
        STIME    = STIME3D
        NVARSSRC = NVARS3D
        NVECSSRC = NROWS3D
        VARSSRC( 1:NVARSSRC ) = VNAME3D( 1:NVARSSRC )
        TYPESRC( 1:NVARSSRC ) = VTYPE3D( 1:NVARSSRC )


C.........  allocate / calculate output-group arrays

        ALLOCATE( COLSRC( NVECSSRC ),
     &            ROWSRC( NVECSSRC ),
     &            GRPSIN( NVECSSRC ),
     &            GRPOUT( NVECSSRC ),
     &            INDXIN( NVECSSRC ),
     &            VECCNT( NVECSSRC ),
     &            VECDEX( NVECSSRC ), STAT = ISTAT )
        CALL CHECKMEM( ISTAT, 'COLSRC...VECDEX', PROGNAME )

        IF ( .NOT.READ3( 'IN_SRCS', 'COL',    1,SDATE,STIME, COLSRC ) ) THEN
            EFLAG = .TRUE.
        END IF

        IF ( .NOT.READ3( 'IN_SRCS', 'ROW',    1,SDATE,STIME, ROWSRC ) ) THEN
            EFLAG = .TRUE.
        END IF

        IF ( .NOT.READ3( 'IN_SRCS', 'IGROUP', 1,SDATE,STIME, GRPSIN ) ) THEN
            EFLAG = .TRUE.
        END IF

        IF ( EFLAG ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error reading "IN_SRCS"', 2 )
        END IF


C.........  Compute the input-vector :: output-group mapping

        DO M = 1, NVECSSRC

            INDXIN( M ) = M

            GRPOUT( M ) = GRPSIN( M )
            DO K = 1, NGRPS
                IF ( GRPSIN( M ) .EQ. GID_IN( K ) ) THEN
                    GRPOUT( M ) = GIDOUT( K )
                    EXIT
                END IF
            END DO      !!  end loop on input-IDs for this group

        END DO

        CALL SORTI( NVECSSRC, INDXIN, GRPOUT, ROWSRC, COLSRC )

        VECCNT = 0      !!  array
        VECDEX = 0      !!  array
        L      = 0
        M      = 0
        IGRP   = -9999
        ICOL   = -9999
        IROW   = -9999

        DO N = 1, NVECSSRC

            K = INDXIN( N )
            M = M + 1
            IF ( IGRP .EQ. GRPOUT(K) .AND. 
     &           IROW .EQ. ROWSRC(K) .AND. 
     &           ICOL .EQ. COLSRC(K) ) THEN     !!  repeated (grpID,row,col)
                VECCNT(L) = VECCNT(L) + 1
                VECDEX(M) = K
            ELSE                                !!  new (grpID,row,col)
                L    = L + 1
                IGRP = GRPOUT(K)
                IROW = ROWSRC(K)
                ICOL = COLSRC(K)
                VECCNT(L) = 1
                VECDEX(M) = K
            END IF

        END DO

C.........  No group is not allowed for point sources.
        IF( CATEGORY == 'P' ) THEN
            L = NVECSSRC
            VECCNT = 1
            VECDEX = INDXIN
        END IF

        NVECOUT = L


C.........  Open the OUT_SRCS, borrowing much of description from IN_SRCS

        NROWS3D = NVECOUT

        IF ( .NOT.OPEN3( 'OUT_SRCS', FSUNKN3, PROGNAME ) ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error(s) opening OUT_SRCS file', 2 )
        END IF
        CALL M3MESG( BLANK )


C.........  Write the OUT_SRCS

        DO V = 1, NVARSSRC

            IF ( TYPESRC(V) .EQ. M3INT ) THEN
                IF ( .NOT.WRINT_AGGGRPS( 'IN_SRCS', 'OUT_SRCS', VARSSRC(V), SDATE, STIME ) ) THEN
                    EFLAG = .TRUE.
                END IF
            ELSE IF ( TYPESRC(V) .EQ. M3REAL ) THEN
                IF ( .NOT.WRREAL_AGGGRPS( 'IN_SRCS', 'OUT_SRCS', VARSSRC(V), SDATE, STIME ) ) THEN
                    EFLAG = .TRUE.
                END IF
            ELSE IF ( TYPESRC(V) .EQ. M3DBLE ) THEN
                IF ( .NOT.WRDBLE_AGGGRPS( 'IN_SRCS', 'OUT_SRCS', VARSSRC(V), SDATE, STIME ) ) THEN
                    EFLAG = .TRUE.
                END IF
            END IF

        END DO

        IF ( EFLAG ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error(s) writing OUT_SRCS files', 2 )
        END IF


C.........  Open and check the IN_EMIS,OUT_EMIS files

        IF ( .NOT.OPEN3( 'IN_EMIS', FSREAD3, PROGNAME ) ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error opening "IN_EMIS"', 2 )
        ELSE IF ( .NOT.DESC3( 'IN_EMIS' ) ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error in  DESC3( IN_EMIS )', 2 )
        END IF

        IF ( NLAYS3D .NE. 1 .OR. NCOLS3D .NE. 1 ) THEN
            CALL M3MESG( 'Unexpected COL, LAY dims for "IN_EMIS"' )
            EFLAG = .TRUE.
        END IF

        IF ( SDATE3D .NE. SDATE ) THEN
            CALL M3MESG( 'Inconsistent SDATE in "IN_EMIS"' )
            EFLAG = .TRUE.
        END IF

        IF ( STIME3D .NE. STIME ) THEN
            CALL M3MESG( 'Inconsistent STIME in "IN_EMIS"' )
            EFLAG = .TRUE.
        END IF

        IF ( NROWS3D .NE. NVECSSRC ) THEN
            CALL M3MESG( 'Inconsistent dim NROWS in "IN_EMIS"' )
            EFLAG = .TRUE.
        END IF

        IF ( EFLAG ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error(s):  inconsistent "IN_EMIS"', 2 )
        END IF

        TSTEP    = TSTEP3D
        NRECS    = MXREC3D
        NROWS3D  = NVECOUT
        NVAREMIS = NVARS3D
        VAREMIS( 1:NVAREMIS ) = VNAME3D( 1:NVAREMIS )
        TYPEMIS( 1:NVAREMIS ) = VTYPE3D( 1:NVAREMIS )

        IF ( .NOT.OPEN3( 'OUT_EMIS', FSUNKN3, PROGNAME ) ) THEN
            CALL M3EXIT( PROGNAME, 0,0, 'Error(s) opening OUT_EMIS file', 2 )
        END IF
        CALL M3MESG( BLANK )


C.........  Write the OUT_EMIS

        JDATE = SDATE
        JTIME = STIME

        DO M = 1, NRECS

            DO V = 1, NVAREMIS

                IF ( TYPEMIS( V ) .EQ. M3INT ) THEN
                    IF ( .NOT.WRINT_AGGGRPS( 'IN_EMIS', 'OUT_EMIS',
     &                                VAREMIS(V), JDATE, JTIME ) ) THEN
                        EFLAG = .TRUE.
                    END IF
                ELSE IF ( TYPEMIS( V ) .EQ. M3REAL ) THEN
                    IF ( .NOT.WRREAL_AGGGRPS( 'IN_EMIS', 'OUT_EMIS',
     &                                 VAREMIS(V), JDATE, JTIME ) ) THEN
                        EFLAG = .TRUE.
                    END IF
                ELSE IF ( TYPEMIS( V ) .EQ. M3DBLE ) THEN
                    IF ( .NOT.WRDBLE_AGGGRPS( 'IN_EMIS', 'OUT_EMIS',
     &                                 VAREMIS(V), JDATE, JTIME ) ) THEN
                        EFLAG = .TRUE.
                    END IF
                END IF

            END DO

            CALL NEXTIME( JDATE, JTIME, TSTEP )

        END DO

        
        IF ( EFLAG ) THEN
            MESG = 'Error(s) writing OUT_EMIS files'
        ELSE
            MESG = 'Successful completion'
        END IF

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )


      CONTAINS    !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


            SUBROUTINE WRRPT_GRPS( )

            INTEGER         I, RDEV

            RDEV = GETEFILE( 'OUT_REPORT', .FALSE., .TRUE., PROGNAME )
            IF ( INRPRT .LT. 0 ) THEN
                CALL M3EXIT( PROGNAME, 0,0,  'Could not open "OUT_REPORT"', 2 )
            END IF

            IF ( CATEGORY .EQ. 'A' .OR. CATEGORY .EQ. 'M' ) THEN

                DO I = 1, NSRCALL
                    WRITE( RDEV, '(I8,1X,A,1X,A,1X,I8,1X,I8)' )
     &                  I, CFIPIN(I), CSCCIN(I), IGROUP(I), SGROUP(I)
                END DO

            ELSE IF ( CATEGORY .EQ. 'B' ) THEN

                DO I = 1, NSRCALL
                    WRITE( RDEV, '(A,1X,I8,1X,I8)' ) CFIPIN(I), IGROUP(I), SGROUP(I)
                END DO

            ELSE IF ( CATEGORY .EQ. 'P' ) THEN

                DO I = 1, NSRCALL
                    WRITE( RDEV, '(I8,1X,A,1X,A,1X,A,1X,A,1X,I8,1X,I8)' )
     &                  I, CFIPIN(I), CSCCIN(I), CPLTIN(I), CPNTIN(I), IGROUP(I), SGROUP(I)
                END DO

            ELSE
                CALL M3EXIT( PROGNAME, 0,0, 'Unrecognized category '//CATEGORY, 2 )
            END IF

            RETURN

            END SUBROUTINE WRRPT_GRPS


        !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


            LOGICAL FUNCTION WRINT_AGGGRPS( IFILE, OFILE, VNAME, JDATE, JTIME )

            CHARACTER(LEN=*), INTENT( IN ) :: IFILE
            CHARACTER(LEN=*), INTENT( IN ) :: OFILE
            CHARACTER(LEN=*), INTENT( IN ) :: VNAME
            INTEGER         , INTENT( IN ) :: JDATE, JTIME

            INTEGER     IBUF( NVECSSRC )
            INTEGER     OBUF( NVECOUT )

            INTEGER     K, L, N, SUM

            CHARACTER(NAMLEN3), PARAMETER :: NONAGGVARS( 7 ) =
     &       (/ 'ISTACK  ', 'IGROUP  ', 'ROW     ', 'COL     ',
     &          'IFIP    ', 'LMAJOR  ', 'LPING   '   /)

            IF ( .NOT.READ3( IFILE, VNAME, 1,JDATE,JTIME, IBUF ) ) THEN
                WRINT_AGGGRPS = .FALSE.
                RETURN
            END IF

            !!  "Select first available data" for NONAGGVARS variables

            IF ( INDEX1( VNAME, 7, NONAGGVARS ) .GT. 0 ) THEN

                IF( VNAME == 'IGROUP  ' ) IBUF = GRPOUT

                K = 1
                DO N = 1, NVECOUT
                    OBUF(N) = IBUF( VECDEX(K) )
                    K = K + VECCNT(N)
                END DO

            ELSE    !!  "Compute Sum"  for  other variables

                K = 0
                DO N = 1, NVECOUT
                    SUM = 0.0d0
                    DO L = 1 , VECCNT(N)
                        K   = K   + 1
                        SUM = SUM + IBUF( VECDEX(K) )
                    END DO
                    OBUF(N) = SUM
                END DO

            END IF

            WRINT_AGGGRPS = WRITE3( OFILE, VNAME,JDATE,JTIME, OBUF )
            RETURN

            END FUNCTION WRINT_AGGGRPS


        !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


            LOGICAL FUNCTION WRREAL_AGGGRPS( IFILE, OFILE, VNAME, JDATE, JTIME )

            CHARACTER(LEN=*), INTENT( IN ) :: IFILE
            CHARACTER(LEN=*), INTENT( IN ) :: OFILE
            CHARACTER(LEN=*), INTENT( IN ) :: VNAME
            INTEGER         , INTENT( IN ) :: JDATE, JTIME

            REAL        IBUF( NVECSSRC )
            REAL        OBUF( NVECOUT )
            REAL*8      SUM             !! high-precision accumulator to minimize round-off error
            INTEGER     K, L, N

            CHARACTER(NAMLEN3), PARAMETER :: NONAGGVARS( 4 ) =
     &       (/ 'LATITUDE  ', 'LONGITUDE ', 'XLOCA     ', 'YLOCA     '   /)

            CHARACTER(NAMLEN3), PARAMETER :: MEANVARS( 5 ) =
     &       (/ 'STKDM ', 'STKHT ', 'STKTK ', 'STKVE ', 'STKFLW' /)

            IF ( .NOT.READ3( IFILE, VNAME, 1,JDATE,JTIME, IBUF ) ) THEN
                WRREAL_AGGGRPS = .FALSE.
                RETURN
            END IF

            !!  "Select first available data" for NONAGGVARS variables

            IF ( INDEX1( VNAME, 4, NONAGGVARS ) .GT. 0 ) THEN

                K = 1
                DO N = 1, NVECOUT
                    OBUF(N) = IBUF( VECDEX(K) )
                    K = K + VECCNT(N)
                END DO

            ELSE IF ( INDEX1( VNAME, 5, MEANVARS ) .GT. 0 ) THEN    !!  "Compute Mean" for MEANVARS variables

                K = 0
                DO N = 1, NVECOUT
                    SUM = 0.0d0
                    DO L = 1 , VECCNT(N)
                        K   = K   + 1
                        SUM = SUM + IBUF( VECDEX(K) )
                    END DO
                    OBUF(N) = SUM / DBLE( VECCNT(N) )
                END DO

            ELSE            !!  "Compute Sum" for other variables

                K = 0
                DO N = 1, NVECOUT
                    SUM = 0.0d0
                    DO L = 1 , VECCNT(N)
                        K   = K   + 1
                        SUM = SUM + IBUF( VECDEX(K) )
                    END DO
                    OBUF(N) = SUM
                END DO

            END IF

            WRREAL_AGGGRPS = WRITE3( OFILE, VNAME,JDATE,JTIME, OBUF )
            RETURN

            END FUNCTION WRREAL_AGGGRPS


        !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


            LOGICAL FUNCTION WRDBLE_AGGGRPS( IFILE, OFILE, VNAME, JDATE, JTIME )

            CHARACTER(LEN=*), INTENT( IN ) :: IFILE
            CHARACTER(LEN=*), INTENT( IN ) :: OFILE
            CHARACTER(LEN=*), INTENT( IN ) :: VNAME
            INTEGER         , INTENT( IN ) :: JDATE, JTIME

            REAL*8      IBUF( NVECSSRC )
            REAL*8      OBUF( NVECOUT )
            REAL*8      SUM
            INTEGER     K, L, N

            IF ( .NOT.READ3( IFILE, VNAME, 1,JDATE,JTIME, IBUF ) ) THEN
                WRDBLE_AGGGRPS = .FALSE.
                RETURN
            END IF

            !!  "Compute Sum"

            DO N = 1, NVECOUT
                SUM = 0.0d0
                DO L = 1 , VECCNT(N)
                    K   = K   + 1
                    SUM = SUM + IBUF( VECDEX(K) )
                END DO
                OBUF(N) = SUM
            END DO

            WRDBLE_AGGGRPS = WRITE3( OFILE, VNAME,JDATE,JTIME, OBUF )
            RETURN

            END FUNCTION WRDBLE_AGGGRPS


      END PROGRAM REGROUP

