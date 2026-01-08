
        SUBROUTINE GENRPRT( FDEV, RCNT, ADEV, MDEV, ENAME, TNAME,
     &                      LNAME, OUTFMT, SMAT, ZEROFLAG, EFLAG )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      The GENRPRT routine is reponsible for generating the columnar 
C      contents of the report.  It will be potentially called multiple 
C      times in one program run to generate many different types of reports.
C
C  PRECONDITIONS REQUIRED:
C      From previous subroutines, we will have a list of records that have 
C      been selected, and this list will include the source IDs (passed
C      through MODREPBN
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C   Created 7/2000 by M Houyoux
C
C   Version 9/2014 by C Coats:  use parallel-mode sparse binning matrices 
C   MODREPBN:<NBINS,ISRCB,GFACB>.  Array TMPBIN reorganization to match
C   Fortran storage order.
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
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
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: CSOURC, POLVAL

C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: QAFMTL3, RPT_, SDATE, STIME, RPTNSTEP,
     &                      AFLAG, ASCREC, NSTEPS, EMLAYS, TSTEP,
     &                      ALLRPT, ALLOUTHR, UCNVFAC, DLFLAG,
     &                      LOC_BEGP, LOC_ENDP

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: NSVARS, NOUTBINS, NOUTREC, OUTSFAC,
     &                      NBINS, BINDATA, ISRCB, ISPRO, GFACB,
     &                      TODOUT, TOSOUT, SPCTOTPR, SPCTOINV,
     &                      INVIDX, TPRIDX, INVTOPRJ, INVTOCMU,
     &                      SPCIDX, OUTSRC, OUTBIN, OUTGFAC,
     &                      BINPOPDIV

C.........  This module contains the temporal profile tables
        USE MODTMPRL, ONLY: NTPDAT, TPNAME

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: ACUMATX, PRMAT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, NIPPA, EAREAD, EANAM, MXCHRS, NCHARS

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
C        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
C        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C...........   EXTERNAL FUNCTIONS
C       CHARACTER(2) CRLF
        INTEGER      MULTUNIT
C       INTEGER      SECSDIFF

C        EXTERNAL     CRLF, MULTUNIT, SECSDIFF
        EXTERNAL     MULTUNIT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT(IN   ) :: FDEV    ! output file unit number
        INTEGER     , INTENT(IN   ) :: RCNT    ! report number
        INTEGER     , INTENT(IN   ) :: ADEV    ! unit no. ASCII elevated file
        INTEGER     , INTENT(IN   ) :: MDEV    ! unit no. source mapping file
        CHARACTER(*), INTENT(IN   ) :: ENAME   ! inventory file name
        CHARACTER(*), INTENT(IN   ) :: TNAME   ! hourly data file name
        CHARACTER(*), INTENT(IN   ) :: LNAME   ! layer fractions file name
        CHARACTER(QAFMTL3),
     &                INTENT(IN   ) :: OUTFMT  ! output record format
        REAL        , INTENT(IN   ) :: SMAT( NSRC, NSVARS ) ! mole spc matrix
        LOGICAL     , INTENT(INOUT) :: ZEROFLAG! true: report zero values
        LOGICAL     , INTENT(INOUT) :: EFLAG   ! true: error occured

C...........   Arrays for source characteristics output formatting
        CHARACTER(300) CHARS ( MXCHRS ) !  source fields for output
        LOGICAL, ALLOCATABLE, SAVE :: LF ( : ) ! true if column should be output

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE, SAVE :: SIDX( : ) ! spc/dat idx for incl pol/act
        REAL   , ALLOCATABLE       :: LFRAC1L( : ) ! layer fractions
        REAL   , ALLOCATABLE       :: TMPBIN( :,:,:,: ) ! array layered data per each time step

C...........   Other local variables
        INTEGER         E, H, I, J, K, L, M, N, P, NC, S, T, V   ! counters and indices
        INTEGER         IS, ID, IE, IP, SE, SP, SS

        INTEGER         IOS               ! i/o status
        INTEGER         JDATE             ! Julian date
        INTEGER         JTIME             ! time (HHMMSS)
        INTEGER      :: KM    = 1         ! index to mult control matrix
        INTEGER      :: KP    = 1         ! index to projection matrix
        INTEGER         LOUT              ! number of output layers
        INTEGER         NDATA             ! number of data columns
        INTEGER         NV                ! number data or spc variables
        INTEGER         SRCNO             ! source no. from ASCII elevated file

        REAL            EMISVAL           ! emissions values from ASCII elevated file
        REAL            F
        REAL*8          BSUM

        LOGICAL      :: FIRSTIME = .TRUE.  ! true: first time routine called
        LOGICAL      :: FIRSTOUT = .TRUE.  ! true: first time routine called
        LOGICAL      :: SFLAG    = .FALSE. ! true: speciation applies to rpt

        CHARACTER(10)         POL         ! species from ASCII elevated file
        CHARACTER(256)        MESG        !  message buffer
        CHARACTER(300)        LINE        !  tmp line buffer
        CHARACTER(IOVLEN3) :: VBUF        !  tmp variable name

        CHARACTER(16), PARAMETER :: PROGNAME = 'GENRPRT' ! program name

C***********************************************************************
C   begin body of subroutine GENRPRT

C.........  First time routine is call
        IF( FIRSTIME ) THEN

C.............  Allocate memory for flagging output non-speciated data
            ALLOCATE( SIDX( NIPPA + NTPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SIDX', PROGNAME )
            
            FIRSTIME = .FALSE.
        END IF

C.............  Allocate memory for LF if not available already
        IF( RPT_%SRCMAP ) THEN
            IF( ALLOCATED( LF ) ) DEALLOCATE( LF )
            ALLOCATE( LF( MXCHRS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LF', PROGNAME )
            LF( 1:NCHARS ) = .TRUE.
            WRITE( MDEV,'(A)' ) 'GROUPID, SRCID, FIPS, FAC_ID, UNIT_ID, REL_POINTID, PROC_ID'
            ZEROFLAG = .TRUE.    ! for a proper match, need to output all sources
        END IF

C.........  Report-specific local settings
        NDATA = ALLRPT( RCNT )%NUMDATA
        RPT_  = ALLRPT( RCNT )     

        SFLAG = ( RPT_%USESLMAT .OR. RPT_%USESSMAT )

C.........  Allocate local memory for reading input data

        N  = NIPPA
        IF( RPT_%USEHOUR ) N = NTPDAT

        ALLOCATE( POLVAL( NSRC, N ),
     &           LFRAC1L( NSRC ),
     &           BINDATA( NOUTBINS, NDATA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVAL...BINDATA', PROGNAME )
        BINDATA = 0.0           !  array

        IF( DLFLAG ) THEN
            ALLOCATE( TMPBIN( NOUTBINS, NDATA, EMLAYS, RPTNSTEP )
     &                       ,STAT=IOS )
            CALL CHECKMEM( IOS, 'TMPBIN', PROGNAME )
            TMPBIN = 0.0    ! array
        END IF

C.........  Set variable loop maxmimum based on speciation status
        NV = NIPPA + NTPDAT
        IF( SFLAG ) NV = NSVARS

C.........  Initialize status of output non-speciated data for this report
        SIDX = 0    ! array

C.........  Loop through time steps
        FIRSTOUT = .TRUE.
        JDATE = SDATE
        JTIME = STIME
        DO T = 1, RPTNSTEP
C.............  Set hour index
            H =  1 + MOD( JTIME / 10000 , 24 )

C...........  Read hourly emissions, if needed
C..............  From temporal file
            IF( RPT_%USEHOUR .AND. .NOT. AFLAG ) THEN
                DO V = 1, NTPDAT

                    VBUF = TPNAME( V )
                    IF( .NOT. READSET( TNAME, VBUF, ALLAYS3, ALLFILES, 
     &                               JDATE, JTIME, POLVAL(1,V) ) ) THEN
                        MESG = 'Could not read "' // TRIM( VBUF ) //
     &                         '" from '// TNAME
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF
                END DO

C.............  From ASCII elevated file
            ELSE IF( RPT_%BYHOUR .AND. AFLAG ) THEN
                DO V = 1, NIPPA
                    DO S = 1, NSRC

                        VBUF = EANAM( V )
                        READ( ADEV, 93010 ) SRCNO, POL, EMISVAL
                        ASCREC = ASCREC + 1

                        IF( SRCNO .NE. S ) THEN
                            POLVAL( S, V ) = 0.
                            BACKSPACE( ADEV )
                            CYCLE

                        ELSE
                            POLVAL( S, V ) = EMISVAL
        
                        END IF

                        IF( POL .NE. VBUF ) THEN
                            WRITE( MESG, '(A,I5)' )
     &                      'Reading in pollutant "' //
     &                      TRIM( VBUF ) // '", but found ' //
     &                      'pollutant "' // TRIM( POL ) //
     &                      '" at line ', ASCREC
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        END IF

                    END DO
                END DO

                IF( T .NE. RPTNSTEP ) THEN
                    DO I = 1, 12
                        ASCREC = ASCREC + 1
                        READ( ADEV, '(A)' ) LINE
                    END DO
                END IF

C...........  Otherwise, read inventory emissions
            ELSE IF( .NOT. RPT_%USEHOUR .AND. .NOT. AFLAG ) THEN
                CALL RDMAPPOL( NSRC, NIPPA, 1, EAREAD, POLVAL )

            ELSE IF( .NOT. RPT_%USEHOUR .AND. AFLAG ) THEN
                POLVAL = 0.
                DO I = 1, NSTEPS
                    DO V = 1, NIPPA
                        DO S = 1, NSRC

                          VBUF = EANAM( V )
                          READ( ADEV, 93010 ) SRCNO, POL, EMISVAL
                          ASCREC = ASCREC + 1

                          IF( SRCNO .NE. S ) THEN
                            BACKSPACE( ADEV )
                            CYCLE

                          ELSE
                            POLVAL( S, V ) = POLVAL( S, V ) +
     &                                       EMISVAL

                          END IF

                          IF( POL .NE. VBUF ) THEN
                            WRITE( MESG, '(A,I5)' )
     &                      'Reading in pollutant "' //
     &                      TRIM( VBUF ) // '", but found ' //
     &                      'pollutant "' // TRIM( POL ) //
     &                      '" at line ', ASCREC
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                          END IF

                        END DO
                    END DO

                    IF( I .NE. NSTEPS ) THEN
                        DO J = 1, 12
                          ASCREC = ASCREC + 1
                          READ( ADEV, '(A)' ) LINE
                        END DO
                    END IF
                END DO

            END IF

C.............  Loop over layers (EMLAYS will be 1 by default)
            LOUT = 1
            IF( RPT_%BYLAYER ) LOUT = EMLAYS
            DO L = 1, LOUT

C.................  If needed for this report, read layer fractions for current
C                   layer, otherwise set to 1.
                IF( RPT_%BYLAYER ) THEN

                    IF( .NOT. READ3( LNAME, 'LFRAC', L,
     &                         JDATE, JTIME, LFRAC1L ) ) THEN
                        MESG = 'Could not read "LFRAC" from '// LNAME
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF

                ELSE
                    LFRAC1L = 1.   ! array
                    
                END IF

C.................  Loop through input data (NV=NIPPA+NTPDAT or NSVARS) and 
C                   sum to bins within list of output records
                DO V = 1, NV

C.....................  Set index to data arrays based on speciation status
                    E = V
                    IF( RPT_%USEHOUR .AND. SFLAG ) THEN
                        E = SPCTOTPR( V )

                    ELSE IF( SFLAG ) THEN
                        E = SPCTOINV( V )

                    END IF

C.....................  Set index from global to actually input pol/act/etype
                    J = INVIDX( E )
                    IF( RPT_%USEHOUR ) J = TPRIDX( E )

C.....................  Skip variable if it is not used for any reports
                    IF( J .EQ. 0 ) CYCLE

C..................  Determine index to projection matrix
                    KP = 1
                    IF( TODOUT( E,RCNT )%PRYN ) 
     &                  KP = MAX( 1, INVTOPRJ( E ) + 1 )

C..................  Determine index to mult control matrix
C..................  Note that first columns is an array of ones
                    KM = 1
                    IF( TODOUT( E,RCNT )%CUYN ) 
     &                  KM = MAX( 1, INVTOCMU( E ) + 1 )

C NOTE: Insert here for reactivity controls.  More formula changes will be needed in
C   N: the formulas below (may be a good idea to apply in separate section?)

C.....................  If speciation, apply speciation factors to appropriate
C                       pollutant and emission types.

                    IF( TODOUT( E,RCNT )%SPCYN ) THEN

C.........................  If current speciation variable used for this report

                        IF( TOSOUT( V,RCNT )%AGG .GT. 0 ) THEN

C.............................  Set index from global to actually input spc vars
C.............................  Set BINDATA-indices

                            K  = SPCIDX( V )

                            IS = TOSOUT(V,RCNT)%SPC
                            IE = TOSOUT(V,RCNT)%ETPSPC
                            IP = TOSOUT(V,RCNT)%PRCSPC
                            SE = TOSOUT(V,RCNT)%SUMETP
                            SP = TOSOUT(V,RCNT)%SUMPOL
                            SS = TOSOUT(V,RCNT)%SUMSPC

C.............................  Sum gridded output records into temporary bins
C.............................  Gridding factor has normalization by cell area

                            IF( RPT_%USEGMAT ) THEN

C.............................  Sum gridded output records into temporary bins
C.............................  Gridding factor has normalization by cell area

!$OMP                           PARALLEL DO
!$OMP&                           DEFAULT( NONE ),
!$OMP&                            SHARED( NOUTBINS, NBINS, ISRCB, GFACB,
!$OMP&                                    POLVAL, SMAT, LFRAC1L, PRMAT,
!$OMP&                                    ACUMATX, BINPOPDIV, BINDATA,
!$OMP&                                    J, K, KP, KM, IS, IE, IP, SE, SP, SS ),
!$OMP&                           PRIVATE( N, BSUM, M, S, F )

                                DO N = 1, NOUTBINS

                                    BSUM = 0.0D0
                                    DO M = NBINS( N-1 )+1, NBINS( N )
                                        S = ISRCB( M )
                                        F = GFACB( M )
                                        BSUM = BSUM + F * POLVAL ( S,J ) *
     &                                                    SMAT   ( S,K ) *
     &                                                    LFRAC1L( S )   *
     &                                                    PRMAT  ( S,KP) *
     &                                                    ACUMATX( S,KM)
                                        IF( RPT_%SRCMAP ) THEN
                                            CALL PARSCSRC( CSOURC( S ), MXCHRS, LOC_BEGP,
     &                                                     LOC_ENDP, LF, NC, CHARS )
                                            IF( FIRSTOUT ) THEN
                                                WRITE( MDEV,94020 ) N, S, ( TRIM(CHARS(K)) ,K=1,NC )
                                            END IF
                                        END IF

                                    END DO

                                    !!  .....  Add temporary bins values to output columns
                                    BSUM = BSUM * BINPOPDIV( N )
                                    IF ( IS .GT. 0 ) BINDATA( N,IS ) = BINDATA( N,IS ) + BSUM
                                    IF ( IE .GT. 0 ) BINDATA( N,IE ) = BINDATA( N,IE ) + BSUM
                                    IF ( IP .GT. 0 ) BINDATA( N,IP ) = BINDATA( N,IP ) + BSUM
                                    IF ( SE .GT. 0 ) BINDATA( N,SE ) = BINDATA( N,SE ) + BSUM
                                    IF ( SP .GT. 0 ) BINDATA( N,SP ) = BINDATA( N,SP ) + BSUM
                                    IF ( SS .GT. 0 ) BINDATA( N,SS ) = BINDATA( N,SS ) + BSUM

                                END DO
                                FIRSTOUT = .FALSE.

                            ELSE      !  else not usegmat

C.............................  Sum non-gridded output records into tmp bins

!$OMP                           PARALLEL DO
!$OMP&                           DEFAULT( NONE ),
!$OMP&                            SHARED( NOUTBINS, NBINS, ISRCB, ISPRO,
!$OMP&                                    POLVAL, SMAT, OUTSFAC, LFRAC1L, PRMAT,
!$OMP&                                    ACUMATX, BINPOPDIV, BINDATA,
!$OMP&                                    J, K, KP, KM, IS, IE, IP, SE, SP, SS ),
!$OMP&                           PRIVATE( N, BSUM, M, S, P )

                                DO N = 1, NOUTBINS

                                    BSUM = 0.0D0
                                    DO M = NBINS( N-1 )+1, NBINS( N )
                                        S = ISRCB( M )
                                        P = ISPRO( M )
                                        BSUM = BSUM + POLVAL ( S,J ) *
     &                                                SMAT   ( S,K ) *
     &                                                LFRAC1L( S )   *
     &                                                PRMAT  ( S,KP) *
     &                                                ACUMATX( S,KM) *
     &                                                OUTSFAC( P )
                                        IF( RPT_%SRCMAP ) THEN
                                            CALL PARSCSRC( CSOURC( S ), MXCHRS, LOC_BEGP,
     &                                                     LOC_ENDP, LF, NC, CHARS )
                                            IF( FIRSTOUT ) THEN
                                                WRITE( MDEV,94020 ) N, S, ( TRIM(CHARS(K)) ,K=1,NC )
                                            END IF
                                        END IF

                                    END DO

                                    !!  .....  Add temporary bins values to output columns
                                    BSUM = BSUM * BINPOPDIV( N )
                                    IF ( IS .GT. 0 ) BINDATA( N,IS ) = BINDATA( N,IS ) + BSUM
                                    IF ( IE .GT. 0 ) BINDATA( N,IE ) = BINDATA( N,IE ) + BSUM
                                    IF ( IP .GT. 0 ) BINDATA( N,IP ) = BINDATA( N,IP ) + BSUM
                                    IF ( SE .GT. 0 ) BINDATA( N,SE ) = BINDATA( N,SE ) + BSUM
                                    IF ( SP .GT. 0 ) BINDATA( N,SP ) = BINDATA( N,SP ) + BSUM
                                    IF ( SS .GT. 0 ) BINDATA( N,SS ) = BINDATA( N,SS ) + BSUM

                                END DO
                                FIRSTOUT = .FALSE.

                            END IF      !  if usegmat, or not

                        END IF

                    END IF       ! end if speciation

C.....................  If used for this report, transfer emission values  
C                       without speciation to temporary bin array
C.....................  Make sure that this data record has not already been 
C                       added to the output columns, as could have happened when
C                       loop is over species (NV=NSVARS)
                    IF( TODOUT( E,RCNT )%AGG .GT. 0 .AND.
     &                ( SIDX( E ) .EQ. 0 .OR. SIDX( E ) .EQ. V ) )THEN

C.........................  Flag data value as already having been added to output

                        SIDX( E ) = V
                        IE = TODOUT( E,RCNT )%ETP
                        ID = TODOUT( E,RCNT )%DAT

                        IF( RPT_%USEGMAT ) THEN

C.........................  Sum gridded output records into temporary bins
C..........................  Gridding factor has normalization by cell area

!$OMP                       PARALLEL DO
!$OMP&                       DEFAULT( NONE ),
!$OMP&                        SHARED( NOUTBINS, NBINS, ISRCB, GFACB,
!$OMP&                                POLVAL, LFRAC1L, PRMAT, ACUMATX,
!$OMP&                                BINDATA, BINPOPDIV, J, KP, KM, ID, IE ),
!$OMP&                       PRIVATE( N, BSUM, M, S, F )

                            DO N = 1, NOUTBINS
                                BSUM = 0.0D0
                                DO M = NBINS( N-1 )+1, NBINS( N )
                                    S = ISRCB( M )
                                    F = GFACB( M )
                                    BSUM = BSUM + F * POLVAL ( S,J ) *
     &                                                LFRAC1L( S )   *
     &                                                PRMAT  ( S,KP) *
     &                                                ACUMATX( S,KM) 
                                    IF( RPT_%SRCMAP ) THEN
                                        CALL PARSCSRC( CSOURC( S ), MXCHRS, LOC_BEGP,
     &                                                 LOC_ENDP, LF, NC, CHARS )
                                        IF( FIRSTOUT ) THEN
                                            WRITE( MDEV,94020 ) N, S, ( TRIM(CHARS(K)) ,K=1,NC )
                                        END IF
                                    END IF
                                END DO

                                !!  .....  Add temporary bins values to output columns
                                BSUM = BSUM * BINPOPDIV( N )
                                IF ( ID .GT. 0 ) BINDATA( N,ID ) = BINDATA( N,ID ) + BSUM
                                IF ( IE .GT. 0 ) BINDATA( N,IE ) = BINDATA( N,IE ) + BSUM

                            END DO
                            FIRSTOUT = .FALSE.

                        ELSE      !  else not usegmat

C.........................  Sum non-gridded output records into temporary bins

!$OMP                       PARALLEL DO
!$OMP&                       DEFAULT( NONE ),
!$OMP&                        SHARED( NOUTBINS, NBINS, ISRCB, OUTSFAC,
!$OMP&                                POLVAL, LFRAC1L, PRMAT, ACUMATX,
!$OMP&                                BINDATA, BINPOPDIV, J, KP, KM, ID, IE ),
!$OMP&                       PRIVATE( N, BSUM, M, S, P )

                            DO N = 1, NOUTBINS
                                BSUM = 0.0D0
                                DO M = NBINS( N-1 )+1, NBINS( N )
                                    S = ISRCB( M )
                                    P = ISPRO( M )
                                    BSUM = BSUM + POLVAL ( S,J ) *
     &                                            LFRAC1L( S )   *
     &                                            PRMAT  ( S,KP) *
     &                                            ACUMATX( S,KM) *
     &                                            OUTSFAC( P )
                                    IF( RPT_%SRCMAP ) THEN
                                        CALL PARSCSRC( CSOURC( S ), MXCHRS, LOC_BEGP,
     &                                                 LOC_ENDP, LF, NC, CHARS )
                                        IF( FIRSTOUT ) THEN
                                            WRITE( MDEV,94020 ) N, S, ( TRIM(CHARS(K)) ,K=1,NC )
                                        END IF
                                    END IF

                                END DO

                                !!  .....  Add temporary bins values to output columns
                                BSUM = BSUM * BINPOPDIV( N )
                                IF ( ID .GT. 0 ) BINDATA( N,ID ) = BINDATA( N,ID ) + BSUM
                                IF ( IE .GT. 0 ) BINDATA( N,IE ) = BINDATA( N,IE ) + BSUM

                            END DO
                            FIRSTOUT = .FALSE.

                        END IF      !  if usegmat, or not

                    END IF  ! End if current pollutant

                END DO      ! End loop on data or speciation variables

C.................  If this is an output hour when it is not BY LAYER
                IF( ALLOUTHR( H,RCNT ) ) THEN

C.....................  Convert units of output data
                    DO J = 1, NDATA
                        BINDATA( :,J ) = BINDATA( :,J ) * UCNVFAC( J )

C.......................... Store tmp bindata for summing daily layered emission later
                        IF( DLFLAG ) TMPBIN( :,J,L,T ) = BINDATA( :,J )
                    END DO

C..................... Skip writing hourly emission when need daily total layered emissions
                    IF( DLFLAG ) THEN
                        BINDATA = 0
                        CYCLE
                    END IF

C.....................  Write emission totals
                    CALL WRREPOUT( FDEV, RCNT, NDATA, JDATE, JTIME, 
     &                             L, RPT_%DELIM, OUTFMT, ZEROFLAG, 
     &                             EFLAG )
C.....................  Reinitialize sum array
                    BINDATA = 0  ! array

                END IF

            END DO          ! End loop on layers

C.............  If error occured, end writing of report
            IF( EFLAG ) THEN
                WRITE( MESG,94010 )
     &                 'WARNING: Incomplete writing for report', RCNT, 
     &                 'because error(s) occurred.'
                CALL M3WARN( PROGNAME, JDATE, JTIME, MESG )
                RETURN
            END IF

C.............  Increment time step
            CALL NEXTIME( JDATE, JTIME, TSTEP )

        END DO    ! End loop over time steps

C.........  Summing hourly layered for total daily laytered emission when DLFLAG is set to TRUE
C           Skipped hourly output at line 464 and sum hourly emissions for daily total
        IF( DLFLAG ) THEN

            DO L = 1, LOUT
                BINDATA = 0   ! array (summing layered hourly emissions)

                DO T = 1, RPTNSTEP-1

                    DO J = 1, NDATA

                        BINDATA( :,J ) = BINDATA( :,J ) 
     &                                   + TMPBIN( :,J,L,T )
                    END DO

                END DO
C.................  Write daily layered emission totals
                CALL WRREPOUT( FDEV, RCNT, NDATA, JDATE, JTIME, L,
     &                      RPT_%DELIM, OUTFMT, ZEROFLAG, EFLAG )

            END DO

        END IF

C.........  Deallocate routine-specific memory
        IF( DLFLAG ) DEALLOCATE( TMPBIN )
        DEALLOCATE( POLVAL, LFRAC1L, BINDATA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( I10, A10, F10.3 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )
94020   FORMAT( 2I8, 7( 1X, A ) )

        END SUBROUTINE GENRPRT

