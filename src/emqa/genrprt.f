
        SUBROUTINE GENRPRT( FDEV, RCNT, ENAME, TNAME, OUTFMT, 
     &                      SMAT, EFLAG )

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
C COPYRIGHT (C) 2000, MCNC--North Carolina Supercomputing Center
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
C***********************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC

C.........  This module contains Smkreport-specific settings
        USE MODREPRT

C.........  This module contains report arrays for each output bin
        USE MODREPBN

C.........  This module contains the temporal profile tables
        USE MODTMPRL

C.........  This module contains arrays for plume-in-grid and major sources
c       USE MODELEV

C.........  This module contains the lists of unique source characteristics
c       USE MODLISTS

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS
        CHARACTER*2 CRLF
        INTEGER     MULTUNIT
        INTEGER     SECSDIFF

        EXTERNAL    CRLF, MULTUNIT, SECSDIFF

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV    ! output file unit number
        INTEGER     , INTENT (IN) :: RCNT    ! report number
        CHARACTER(*), INTENT (IN) :: ENAME   ! inventory file name
        CHARACTER(*), INTENT (IN) :: TNAME   ! hourly data file name
        CHARACTER(LEN=QAFMTL3),
     &                INTENT (IN) :: OUTFMT  ! output record format
        REAL        , INTENT (IN) :: SMAT( NSRC, NSVARS ) ! mole spc matrix
        LOGICAL     , INTENT(OUT) :: EFLAG   ! true: error occured

C...........   Local allocatable arrays
        INTEGER, ALLOCATABLE, SAVE :: SIDX( : )   ! spc/dat idx for incl pol/act

        REAL, ALLOCATABLE, SAVE :: BINARR( : ) ! helping sum to bins of data

C...........   Other local variables
        INTEGER          E, H, I, J, K, N, S, T, V   ! counters and indices

        INTEGER         IOS               ! i/o status
        INTEGER         JDATE             ! Julian date
        INTEGER         JTIME             ! time (HHMMSS)
        INTEGER         NDATA             ! number of data columns
        INTEGER         NDIN              ! no. data variables to read in
        INTEGER         NV                ! number data or spc variables

        LOGICAL      :: FIRSTIME = .TRUE.  ! true: first time routine called
        LOGICAL      :: SFLAG    = .FALSE. ! true: speciation applies to rpt

        CHARACTER(LEN=IOVLEN3), SAVE :: BNAME    ! file name buffer

        CHARACTER*300          MESG        !  message buffer

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

        BNAME = ENAME
        NDIN  = NIPPA
        IF( RPT_%USEHOUR ) THEN
            BNAME = TNAME
            NDIN  = NTPDAT
        END IF

C.........  Allocate local memory for reading input data
        ALLOCATE( POLVAL( NSRC, NDIN ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVAL', PROGNAME )

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

C.............  Read emissions for whole inventory for possible pollutants
C               and activities for current report
            CALL RDINVPOL( BNAME, NSRC, NDIN, JDATE, JTIME, 
     &                     RDNAMES( 1, RCNT ), POLVAL, IOS )

            IF( IOS .GT. 0 ) THEN
                MESG = 'Problem reading inventory data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Loop through input data (NV=NIPPA+NTPDAT or NSVARS) and sum to 
C               bins within list of output records
            DO V = 1, NV

C.................  Set index to data arrays based on speciation status
                E = V
                IF( SFLAG ) E = SPCTODAT( V )

C.................  Set index from global to actually input pol/act/etype
                J = INVIDX( E )
                IF( RPT_%USEHOUR ) J = TPRIDX( E )

C.................  If speciation, apply speciation factors to appropriate
C                   pollutant and emission types.
                IF( TODOUT( E,RCNT )%SPCYN ) THEN

C.....................  If current speciation variable used for this report
                    IF( TOSOUT( V,RCNT )%AGG .GT. 0 ) THEN

C.........................  Initialize temporary bin sum array
                        BINARR = 0   ! array

C.........................  Set index from global to actually input spc vars
                        K = SPCIDX( V )

C.........................  Sum gridded output records into temporary bins
                        IF( RPT_%USEGMAT ) THEN
                            DO I = 1, NOUTREC
                                S = OUTSRC( I )
                                N = OUTBIN( I )
                                BINARR( N ) = BINARR ( N ) + 
     &                                        OUTGFAC( I )   *
     &                                        POLVAL ( S,J ) * 
     &                                        SMAT   ( S,K )
                            END DO

C.........................  Sum non-gridded output records into temporary bins
                        ELSE
                            DO I = 1, NOUTREC
                                S = OUTSRC( I )
                                N = OUTBIN( I )
                                BINARR( N ) = BINARR ( N ) + 
     &                                        POLVAL ( S,J ) * 
     &                                        SMAT   ( S,K )
                            END DO

                        END IF

C.........................  Add temporary bins values to output columns
                        CALL UPDATE_OUTCOL( TOSOUT( V,RCNT )%SPC )
                        CALL UPDATE_OUTCOL( TOSOUT( V,RCNT )%ETPSPC )
                        CALL UPDATE_OUTCOL( TOSOUT( V,RCNT )%PRCSPC )
                        CALL UPDATE_OUTCOL( TOSOUT( V,RCNT )%SUMETP )
                        CALL UPDATE_OUTCOL( TOSOUT( V,RCNT )%SUMPOL )

                    END IF

                END IF

C.................  If used for this report, transfer emission values without 
C                   speciation to temporary bin array
C.................  Make sure that this data record has not already been added
C                   to the output columns, as could have happened when
C                   loop is over species (NV=NSVARS)
                IF( TODOUT( E,RCNT )%AGG .GT. 0 .AND.
     &            ( SIDX( E ) .EQ. 0 .OR. SIDX( E ) .EQ. V ) )THEN

C.....................  Flag data value as already having been added to output
                    SIDX( E ) = V

C.....................  Initialize temporary bin sum array
                    BINARR = 0   ! array

C.....................  Sum gridded output records into temporary bins
                    IF( RPT_%USEGMAT ) THEN
                        DO I = 1, NOUTREC
                            S = OUTSRC( I )
                            N = OUTBIN( I )
                            BINARR( N ) = BINARR ( N ) + 
     &                                    OUTGFAC( I ) *
     &                                    POLVAL( S,J )
                        END DO

C.....................  Sum non-gridded output records into temporary bins
                    ELSE
                        DO I = 1, NOUTREC
                            S = OUTSRC( I )
                            N = OUTBIN( I )
                            BINARR( N ) = BINARR ( N ) + 
     &                                    POLVAL( S,J )
                        END DO
 
                    END IF

C.....................  Add temporary bins values to output columns
                    CALL UPDATE_OUTCOL( TODOUT( E,RCNT )%ETP )
                    CALL UPDATE_OUTCOL( TODOUT( E,RCNT )%DAT )

                END IF  ! End if current pollutant

            END DO      ! End loop on data or speciation variables

C.............  If this is an output hour...
            IF( ALLOUTHR( H,RCNT ) ) THEN

C.................  Convert units of output data
                DO J = 1, NDATA
                    BINDATA( :,J ) = BINDATA( :,J ) * UCNVFAC( J )
                END DO

C.................  Write emission totals
                CALL WRREPOUT( FDEV, RCNT, NDATA, JDATE, JTIME, 
     &                         RPT_%DELIM, OUTFMT, EFLAG )

C.................  Reinitialize sum array
                BINDATA = 0  ! array

            END IF

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

C.........  Deallocate routine-specific memory
        DEALLOCATE( POLVAL, BINARR )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

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

