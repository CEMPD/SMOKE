
        SUBROUTINE GENRPRT( FDEV, RCNT, HWID, ADEV, ENAME, TNAME,
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
C     Created 7/2000 by M Houyoux
C
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
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: POLVAL

C.........  This module contains Smkreport-specific settings
        USE MODREPRT, ONLY: QAFMTL3, RPT_, SDATE, STIME, RPTNSTEP,
     &                      AFLAG, ASCREC, NSTEPS, EMLAYS, TSTEP,
     &                      ALLRPT, ALLOUTHR, UCNVFAC

C.........  This module contains report arrays for each output bin
        USE MODREPBN, ONLY: NSVARS, NOUTBINS, NOUTREC, BINDATA,
     &                      TODOUT, TOSOUT, SPCTOTPR, SPCTOINV,
     &                      INVIDX, TPRIDX, INVTOPRJ, INVTOCMU,
     &                      SPCIDX, OUTSRC, OUTBIN, OUTGFAC,
     &                      BINPOPDIV

C.........  This module contains the temporal profile tables
        USE MODTMPRL, ONLY: NTPDAT, TPNAME

C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: ACUMATX, PRMAT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NSRC, NIPPA, EAREAD, EANAM

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
        INCLUDE 'IODECL3.EXT'   !  I/O API function declarations
        INCLUDE 'FDESC3.EXT'    !  I/O API file description data structures.
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

C...........   EXTERNAL FUNCTIONS
        CHARACTER*2 CRLF
        INTEGER     MULTUNIT
        INTEGER     SECSDIFF

        EXTERNAL    CRLF, MULTUNIT, SECSDIFF

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV    ! output file unit number
        INTEGER     , INTENT (IN) :: RCNT    ! report number
        INTEGER     , INTENT (IN) :: HWID    ! header width
        INTEGER     , INTENT (IN) :: ADEV    ! unit no. ASCII elevated file
        CHARACTER(*), INTENT (IN) :: ENAME   ! inventory file name
        CHARACTER(*), INTENT (IN) :: TNAME   ! hourly data file name
        CHARACTER(*), INTENT (IN) :: LNAME   ! layer fractions file name
        CHARACTER(LEN=QAFMTL3),
     &                INTENT (IN) :: OUTFMT  ! output record format
        REAL        , INTENT (IN) :: SMAT( NSRC, NSVARS ) ! mole spc matrix
        LOGICAL     , INTENT (IN) :: ZEROFLAG! true: report zero values
        LOGICAL     , INTENT(OUT) :: EFLAG   ! true: error occured

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE, SAVE :: SIDX( : ) ! spc/dat idx for incl pol/act

        REAL, ALLOCATABLE, SAVE :: BINARR( : )  ! helping sum to bins of data

        REAL, ALLOCATABLE, SAVE :: LFRAC1L( : ) ! layer fractions

C...........   Other local variables
        INTEGER          E, H, I, J, K, L, N, S, T, V   ! counters and indices

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

        LOGICAL      :: FIRSTIME = .TRUE.  ! true: first time routine called
        LOGICAL      :: SFLAG    = .FALSE. ! true: speciation applies to rpt

        CHARACTER*10              POL         ! species from ASCII elevated file
        CHARACTER*16           :: RNAME = 'IOAPI_DAT' ! logical name for reading pols
        CHARACTER*256             MESG        !  message buffer
        CHARACTER*300             LINE        !  tmp line buffer
        CHARACTER(LEN=IOVLEN3) :: VBUF        !  tmp variable name

        CHARACTER*16 :: PROGNAME = 'GENRPRT' ! program name

C***********************************************************************
C   begin body of subroutine GENRPRT

C.........  First time routine is call
        IF( FIRSTIME ) THEN

C.............  Allocate memory for flagging output non-speciated data
            ALLOCATE( SIDX( NIPPA + NTPDAT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SIDX', PROGNAME )
            
            FIRSTIME = .FALSE.

        END IF

C.........  Report-specific local settings
        NDATA = ALLRPT( RCNT )%NUMDATA
        RPT_ = ALLRPT( RCNT )     

        SFLAG = ( RPT_%USESLMAT .OR. RPT_%USESSMAT )

C.........  Allocate local memory for reading input data
        N  = NIPPA
        IF( RPT_%USEHOUR ) N = NTPDAT

        ALLOCATE( POLVAL( NSRC, N ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVAL', PROGNAME )
        ALLOCATE( LFRAC1L( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LFRAC1L', PROGNAME )

C.........  Allocate local memory for bin helper summing array
        ALLOCATE( BINARR( NOUTBINS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BINARR', PROGNAME )

C.........  Set variable loop maxmimum based on speciation status
        NV = NIPPA + NTPDAT
        IF( SFLAG ) NV = NSVARS

C.........  Initialize status of output non-speciated data for this report
        SIDX = 0    ! array

C.........  Loop through time steps
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

C.............................  Initialize temporary bin sum array
                            BINARR = 0   ! array

C.............................  Set index from global to actually input spc vars
                            K = SPCIDX( V )

C.............................  Sum gridded output records into temporary bins
C.............................  Gridding factor has normalization by cell area
                            IF( RPT_%USEGMAT ) THEN
                                DO I = 1, NOUTREC
                                    S = OUTSRC( I )
                                    N = OUTBIN( I )
                                    BINARR( N ) = BINARR ( N ) + 
     &                                            OUTGFAC( I )   *
     &                                            POLVAL ( S,J ) * 
     &                                            SMAT   ( S,K ) *
     &                                            LFRAC1L( S )   *
     &                                            PRMAT  ( S,KP) *
     &                                            ACUMATX( S,KM) *
     &                                          BINPOPDIV( N )
                                END DO

C.............................  Sum non-gridded output records into tmp bins
                            ELSE
                                DO I = 1, NOUTREC
                                    S = OUTSRC( I )
                                    N = OUTBIN( I )
                                    BINARR( N ) = BINARR ( N ) + 
     &                                            POLVAL ( S,J ) * 
     &                                            SMAT   ( S,K ) *
     &                                            LFRAC1L( S )   *
     &                                            PRMAT  ( S,KP) *
     &                                            ACUMATX( S,KM) *
     &                                          BINPOPDIV( N )
                                END DO

                            END IF

C.............................  Add temporary bins values to output columns
                            CALL UPDATE_OUTCOL( TOSOUT(V,RCNT)%SPC )
                            CALL UPDATE_OUTCOL( TOSOUT(V,RCNT)%ETPSPC )
                            CALL UPDATE_OUTCOL( TOSOUT(V,RCNT)%PRCSPC )
                            CALL UPDATE_OUTCOL( TOSOUT(V,RCNT)%SUMETP )
                            CALL UPDATE_OUTCOL( TOSOUT(V,RCNT)%SUMPOL )

                        END IF

                    END IF       ! end if speciation

C.....................  If used for this report, transfer emission values  
C                       without speciation to temporary bin array
C.....................  Make sure that this data record has not already been 
C                       added to the output columns, as could have happened when
C                       loop is over species (NV=NSVARS)
                    IF( TODOUT( E,RCNT )%AGG .GT. 0 .AND.
     &                ( SIDX( E ) .EQ. 0 .OR. SIDX( E ) .EQ. V ) )THEN

C.........................  Flag data value as already having been added to
C                           output
                        SIDX( E ) = V

C.........................  Initialize temporary bin sum array
                        BINARR = 0   ! array

C.........................  Sum gridded output records into temporary bins
C..........................  Gridding factor has normalization by cell area
                        IF( RPT_%USEGMAT ) THEN
                            DO I = 1, NOUTREC
                                S = OUTSRC( I )
                                N = OUTBIN( I )
                                BINARR( N ) = BINARR ( N ) + 
     &                                        OUTGFAC( I )   *
     &                                        POLVAL ( S,J ) *
     &                                        LFRAC1L( S )   *
     &                                        PRMAT  ( S,KP) *
     &                                        ACUMATX( S,KM) *
     &                                      BINPOPDIV( N )
                            END DO

C.........................  Sum non-gridded output records into temporary bins
                        ELSE
                            DO I = 1, NOUTREC
                                S = OUTSRC( I )
                                N = OUTBIN( I )
                                BINARR( N ) = BINARR ( N ) + 
     &                                        POLVAL ( S,J ) *
     &                                        LFRAC1L( S )   *
     &                                        PRMAT  ( S,KP) *
     &                                        ACUMATX( S,KM) *
     &                                      BINPOPDIV( N )
                            END DO

                        END IF

C.........................  Add temporary bins values to output columns
                        CALL UPDATE_OUTCOL( TODOUT( E,RCNT )%ETP )
                        CALL UPDATE_OUTCOL( TODOUT( E,RCNT )%DAT )

                    END IF  ! End if current pollutant

                END DO      ! End loop on data or speciation variables

C.................  If this is an output hour...
                IF( ALLOUTHR( H,RCNT ) ) THEN

C.....................  Convert units of output data
                    DO J = 1, NDATA
                        BINDATA( :,J ) = BINDATA( :,J ) * UCNVFAC( J )
                    END DO

C.....................  Write emission totals
                    CALL WRREPOUT( FDEV, RCNT, NDATA, JDATE, JTIME, 
     &                             L,  RPT_%DELIM, OUTFMT, ZEROFLAG, 
     &                             EFLAG )

C.....................  Reinitialize sum array
                    BINDATA = 0  ! array

                END IF

            END DO          ! End loop on layers

C.............  If error occured, end writing of report
            IF( EFLAG ) THEN
                WRITE( MESG,94010 ) 'WARNING: Incomplete writing for '//
     &                 'report', RCNT, 'because error(s) occurred.'
                CALL M3MSG2( MESG )
                RETURN
            END IF

C.............  Increment time step
            CALL NEXTIME( JDATE, JTIME, TSTEP )

        END DO    ! End loop over time steps

C.........  Write line to separate reports from each other and from metadata
        WRITE( FDEV, '(/,A,/)' ) REPEAT( '#', HWID )

C.........  Deallocate routine-specific memory
        DEALLOCATE( POLVAL, LFRAC1L, BINARR )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93010   FORMAT( I10, A10, F10.3 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************
 
        CONTAINS
 
C.............  This internal function updates the output bin columns with
C               the available array of data
            SUBROUTINE UPDATE_OUTCOL( OUTCOL )

C.............  Subprogram arguments
            INTEGER, INTENT (IN) :: OUTCOL

C----------------------------------------------------------------------

            IF( OUTCOL .GT. 0 ) THEN

                BINDATA( :,OUTCOL ) = BINDATA( :,OUTCOL ) + BINARR   ! array

            END IF

            RETURN
 
            END SUBROUTINE UPDATE_OUTCOL

        END SUBROUTINE GENRPRT

