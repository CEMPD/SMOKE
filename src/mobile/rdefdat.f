
        SUBROUTINE RDEFDAT( FDEV )
   
C***********************************************************************
C  subroutine RDEFDAT body starts at line < >
C
C  DESCRIPTION:
C      Pre-process emission factors inputs file.  This involves identifying
C      the model(s) being used, flaging the types, flagging combination
C      emission factors, and allocating and storing various arrays.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION HISTORY:
C
C***************************************************************************
C 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C 
C COPYRIGHT (C) 1999, MCNC--North Carolina Supercomputing Center
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

C...........   MODULES for public variables
C.........  This module contains emission factor tables and related
        USE MODEMFAC

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'M5CNST3.EXT'   !  Mobile5a/b constants

C...........   EXTERNAL FUNCTIONS:
        CHARACTER*2  CRLF
        EXTERNAL     CRLF

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! file unit number

C...........   Local alloctable arrays...

C...........   Temporary array for reading in combination PSIs
        INTEGER, ALLOCATABLE :: CPSIS( : )  ! PSIs
        REAL   , ALLOCATABLE :: CFACS( : )  ! factors

C...........   Commons linked with MOBILE subroutine
        CHARACTER*20    M5VRSION ! name of version of MOBILE used
        COMMON / MULTISC / M5VRSION 

C...........   Other arrays
        INTEGER         M5PSI( MXM5SCEN )  ! PSIs codes for pure PSIs

C...........   Local variables
        INTEGER         I, J, K, L1, L2   ! counters and indices

        INTEGER         CINDX       ! index to combination PSIs table
        INTEGER         CNTCOMBO    ! tmp count of 
        INTEGER         CNTPSI      ! tmp count of PSIs in scenario group
        INTEGER         ICHK        ! counter check value
        INTEGER         IOS         ! i/o status
        INTEGER         IREC        ! record counter
        INTEGER         PSI         ! tmp PSI
        INTEGER         PPSI        ! tmp PSI from previous iteration

        LOGICAL      :: EFLAG    = .FALSE.  ! true: error found
        LOGICAL      :: MOB5AON  = .FALSE.  ! true: Mobile5a packet encountered
        LOGICAL      :: MOB5BON  = .FALSE.  ! true: Mobile5b packet encountered
        LOGICAL      :: SKIPLINE = .FALSE.  ! true: skip current line in file

        CHARACTER*300   LINE     ! line from file
        CHARACTER*300   MESG     ! message buffer

        CHARACTER*16 :: PROGNAME = 'RDEFDAT' ! program name

C***********************************************************************
C   begin body of subroutine RDEFDAT

C.........  Loop through emission factor data input file, to determine sizes
C           for allocating memory.  Count the total number of PSIs and the
C           number of combination PSIs
        IREC     = 0
        NPSIDAT  = 0      ! from MODEMFAC
        NCMBOPSI = 0      ! from MODEMFAC
        MXNCOMBO = 0      ! from MODEMFAC
        MOB5AON  = .FALSE.
        MOB5BON  = .FALSE.
        M5VRSION = ' '
        DO

            READ( FDEV, 93000, END=55, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check read status
            SKIPLINE = .FALSE.
            CALL CHECK_READ_STATUS
            IF( SKIPLINE ) CYCLE

C.............  Check for type of line, and determine count of PSIs
            CALL CHECK_EFDAT_LINE( CNTPSI, CNTCOMBO )

C.............  Skip the current line
            IF( SKIPLINE ) CYCLE

C.............  Add to the number of PSIs that are in the file
            NPSIDAT  = NPSIDAT + CNTPSI

C.............  Add to the number of combination PSIs that are in the file
            NCMBOPSI = NCMBOPSI + CNTCOMBO

C.............  Determine maximum number of contributing PSIs
            MXNCOMBO = MAX( MXNCOMBO, CNTCOMBO )

C.............  If a combination factor was read in, skip the next line too
            IF( CNTCOMBO .GT. 0 ) THEN
                READ( FDEV, *, END=999 )
                IREC = IREC + 1
            END IF

        END DO

55      CONTINUE   ! exit from read loop

C.........  Allocate memory for various arrays for storing the information
C           from this file
        ALLOCATE( PSIDAT( NPSIDAT ),STAT=IOS )     ! sorted
        CALL CHECKMEM( IOS, 'PSIDAT', PROGNAME )
        ALLOCATE( PSIDATA( NPSIDAT ),STAT=IOS )    ! unsorted
        CALL CHECKMEM( IOS, 'PSIDATA', PROGNAME )
        ALLOCATE( PDATINDX( NPSIDAT ),STAT=IOS )
        CALL CHECKMEM( IOS, 'PDATINDX', PROGNAME )
        ALLOCATE( PDATTYPE( NPSIDAT ),STAT=IOS )
        CALL CHECKMEM( IOS, 'PDATTYPE', PROGNAME )
        ALLOCATE( PDATPNTR( NPSIDAT ),STAT=IOS )
        CALL CHECKMEM( IOS, 'PDATPNTR', PROGNAME )
        ALLOCATE( PDATLINE( NPSIDAT ),STAT=IOS )
        CALL CHECKMEM( IOS, 'PDATLINE', PROGNAME )
        ALLOCATE( PDATROOT( NPSIDAT ),STAT=IOS )
        CALL CHECKMEM( IOS, 'PDATROOT', PROGNAME )
        ALLOCATE( PDATMCNT( NPSIDAT ),STAT=IOS )
        CALL CHECKMEM( IOS, 'PDATMCNT', PROGNAME )

C.........  Allocate memory for combination PSIs table
        ALLOCATE( CFACS( MXNCOMBO ),STAT=IOS )
        CALL CHECKMEM( IOS, 'CFACS', PROGNAME )
        ALLOCATE( CPSIS( MXNCOMBO ),STAT=IOS )
        CALL CHECKMEM( IOS, 'CPSIS', PROGNAME )
        ALLOCATE( CMBOPSI( MXNCOMBO, NCMBOPSI  ),STAT=IOS )
        CALL CHECKMEM( IOS, 'CMBOPSI', PROGNAME )
        ALLOCATE( CMBOFAC( MXNCOMBO, NCMBOPSI ),STAT=IOS )
        CALL CHECKMEM( IOS, 'CMBOFAC', PROGNAME )

C.........  Rewind file to read in more information
        REWIND( FDEV )

C.........  Loop through file again and store information for each PSI
        I       = 0
        ICHK    = 0
        IREC    = 0
        CINDX   = 0 
        MOB5AON  = .FALSE.
        MOB5BON  = .FALSE.
        M5VRSION= ' '
        DO

C.............  Read line as character string, which can be the
C.............  first line of 2-line records, or can have a MOBILE5* 
C.............  packet in it.
            READ( FDEV, 93000, END=99, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check read status
            SKIPLINE = .FALSE.
            CALL CHECK_READ_STATUS
            IF( SKIPLINE ) CYCLE

C.............  Check for type of line, determine count of PSIs
            CALL CHECK_EFDAT_LINE( CNTPSI, CNTCOMBO )

C.............  Skip the current line if not the start of a PSI
            IF( SKIPLINE ) CYCLE

C.............  Process for MOBILE inputs
            IF( MOB5AON .OR. MOB5BON ) THEN

C.................  Read the PSIs for current group of Mobile inputs.
C.................  Note that the check routine makes sure that CNTPSI is
C                   not too large
                READ( LINE( L1:L2 ), *, IOSTAT=IOS ) 
     &                CNTPSI, ( M5PSI( K ), K = 1, CNTPSI )

                DO K = 2, CNTPSI

C.....................  Make sure PSIs for a give scenario are sequential
                    IF( M5PSI( K ) .NE. M5PSI( K-1 ) + 1 ) THEN

                         EFLAG = .TRUE.
                         WRITE( MESG, 94010 )
     &                          'PSIs for MOBILE5 packet are not ' //
     &                          'sequential at line', IREC, 
     &                          CRLF() // BLANK10 // 
     &                          'in emission factors data file'

                        CALL M3MESG( MESG )
                        CYCLE   ! to head of read loop

                    END IF
                END DO

            END IF

C.............  Increment counter check (so have a value for output message)
            ICHK = I + CNTPSI

C.............  Go to next line if the memory allocate was somehow insufficient
            IF( ICHK .GT. NPSIDAT ) CYCLE

C.............  For MOBILE5* input type, collect all of the PSI records
            IF( CNTCOMBO .EQ. 0 ) THEN

                DO K = 1, CNTPSI

                    I             = I + 1
                    PSIDATA ( I ) = M5PSI( K )
                    PDATINDX( I ) = I
                    PDATTYPE( I ) = CNTCOMBO
                    PDATPNTR( I ) = K
                    PDATLINE( I ) = IREC
                    PDATROOT( I ) = M5PSI( 1 )
                    PDATMCNT( I ) = CNTPSI

                END DO

C.................  For combo PSIs, increment and store properties
            ELSEIF( CNTCOMBO .GT. 0 ) THEN

                I             = I + CNTPSI
                CINDX         = CINDX + 1  ! indx to contributing EFs for combos

                PSIDATA ( I ) = PSI
                PDATINDX( I ) = I
                PDATTYPE( I ) = CNTCOMBO
                PDATPNTR( I ) = CINDX
                PDATLINE( I ) = IREC
                PDATROOT( I ) = PSI
                PDATMCNT( I ) = 1

C.................  Make sure that combination EF doesn't use more than
C                   the legal number of other PSIs
                IF( CNTCOMBO .GT. MXNCOMBO ) THEN

                    EFLAG = .TRUE.

                    WRITE( MESG,94010 ) 'INTERNAL ERROR: ' //
     &                     'dimension for storing combo PSIs was',
     &                     MXNCOMBO, CRLF() // BLANK10 //
     &                     'but actually needed', CNTCOMBO, 'at line',
     &                     IREC
                    CALL M3MESG( MESG )
                    CYCLE

C.................  No overflow, set combo EF properties...
                ELSE IF( CINDX .LE. NCMBOPSI ) THEN

C.....................  Read contributing PSIs for combination PSI
                    READ( FDEV, *, IOSTAT=IOS ) 
     &                  ( CFACS( K ), CPSIS( K ), K = 1, CNTCOMBO )
                    IREC = IREC + 1

C.....................  Check read status
                    CALL CHECK_READ_STATUS
                    IF( SKIPLINE ) CYCLE

C.....................  Store combination PSI information
                    CMBOPSI( 1:CNTCOMBO, CINDX ) = CPSIS( 1:CNTCOMBO ) 
                    CMBOFAC( 1:CNTCOMBO, CINDX ) = CFACS( 1:CNTCOMBO )

                END IF  ! end of storage of combination PSI factors and PSIs

            END IF

        END DO      

99      CONTINUE   ! exit from read loop

C.........  Check total number of PSIs from this file
        IF( ICHK .GT. NPSIDAT ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'PSI data from emission factor data file was',
     &             NPSIDAT, CRLF() // BLANK10 // 
     &             'but actually needed', ICHK
            CALL M3MSG2( MESG )

        END IF

C.........  Check total number of combination PSIs from this file
        IF( CINDX .GT. NCMBOPSI ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //
     &             'combination emission factor data table was',
     &             NCMBOPSI, CRLF() // BLANK10 // 
     &             'but actually needed', CINDX
            CALL M3MSG2( MESG )

        END IF

C.........  Sort PSI list from emission factor data file
        CALL SORTI1( NPSIDAT, PDATINDX, PSIDATA )

C.........  Create sorted list for searching purposes only
C.........  Also check to ensure no duplicate PSIs
        PPSI = -9        
        DO I = 1, NPSIDAT

           J   = PDATINDX( I )
           PSI = PSIDATA ( J )

           IF( PSI .EQ. PPSI ) THEN
               EFLAG = .TRUE.
               WRITE( MESG,94010 ) 'ERROR: Duplicate PSI code', PSI,
     &                'found at line', PDATLINE( J ), 'in emission ' //
     &                'factors data file.'
               CALL M3MSG2( MESG )

           ELSE
               PSIDAT( I ) = PSI

           END IF

           PPSI = PSI

        END DO

C.........  Deallocate local arrays
        DEALLOCATE( CPSIS, CFACS )

C.........  Abort if an error occured with the file
        IF( EFLAG ) THEN
            MESG = 'Problem with emission factors data file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

999     WRITE( MESG,94010 ) 'Unexpected end-of-file at line', IREC,
     &                      'of actual emission factor data file.'

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C..............  This internal subprogram checks the read status and reports
C                and error if there was one.  It sets EFLAG and skips to 
C                next record

            SUBROUTINE CHECK_READ_STATUS

C.............................................................................

            IF( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &                 'I/O error ', IOS, 'reading emission factor ' //
     &                 'data file at line', IREC
                CALL M3MESG( MESG )
                SKIPLINE = .TRUE.

            END IF

C--------------  SUBPROGRAM FORMAT  STATEMENTS   --------------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_READ_STATUS

C.............................................................................
C.............................................................................

C.............  This internal subprogram checks for the MOBILE packets and
C               sets the count of PSIs and the number for counting the
C               contributing PSIs to a combo PSI (or 0 for a MOBILE input)
            SUBROUTINE CHECK_EFDAT_LINE( CNTPSI, CNTCOMBO )

C.............  Subprogram arguments
            INTEGER, INTENT (IN OUT) :: CNTPSI
            INTEGER, INTENT (IN OUT) :: CNTCOMBO

C.............  Local subprogram variables
            INTEGER  MA, MB

C.............................................................................

C.............  Initialize setting for new line
            SKIPLINE = .FALSE.
            CNTPSI   = 0
            CNTCOMBO = 0

C.............  Check for MOBILE5 packets
            MA = INDEX( LINE, 'MOBILE5A' )
            MB = INDEX( LINE, 'MOBILE5B' )
            L2 = LEN_TRIM( LINE )

C.............  Process for Mobile model inputs
            IF( MA .GT. 0 .OR. MB .GT. 0 ) THEN

C.................  Mobile5a-specific input info
        	IF( MA .GT. 0 ) THEN
                    MOB5AON  = .TRUE.
                    M5VRSION = 'MOBILE5a'
                    L1 = MA + 8

C.....................  Error if combining different versions of Mobile model
                    IF( MOB5BON ) THEN
                	EFLAG = .TRUE.
                	WRITE( MESG,94010 )
     &                    'MOBILE5A packet found, but MOBILE5B ' //
     &                    'used previously in emission factors' //
     &                    CRLF()// BLANK10// 'data file at line', IREC
                	CALL M3MESG( MESG )
                    END IF

C.................  Mobile5b-specific input info
        	ELSE IF( MB .GT. 0 ) THEN
                    MOB5BON = .TRUE.
                    M5VRSION = 'MOBILE5b'
                    L1 = MB + 8

C.....................  Error if combining different versions of Mobile model
                    IF( MOB5AON ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )
     &                    'MOBILE5B packet found, but MOBILE5A ' //
     &                    'used previously in emission factors' //
     &                    CRLF()// BLANK10// 'data file at line', IREC
                        CALL M3MESG( MESG )
                    END IF
                END IF

C.................  Store count of PSIs on this line
                READ( LINE( L1:L2 ), *, IOSTAT=IOS ) CNTPSI

                IF( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                'Bad MOBILE5* packet in emission factors data ' //
     &                'file at line', IREC
                    CALL M3MESG( MESG )
                    SKIPLINE = .TRUE. 

                ELSE IF( CNTPSI .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )
     &                'PSI count is zero in emission factors data ' //
     &                'file at line', IREC
                    CALL M3MESG( MESG )
                    SKIPLINE = .TRUE. 

                ELSE IF( CNTPSI .GT. MXM5SCEN ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )
     &                'Max number of PSIs per scenario is', MXM5SCEN,
     &                'but', CNTPSI, 'requested in emission' //
     &                CRLF()// BLANK10// 'factors data file at line', 
     &                IREC
                    CALL M3MESG( MESG )
                    SKIPLINE = .TRUE. 

                END IF

C.............  For combination formats...
C.............  These need to be in the file before any MOBILE input sections
            ELSE IF( .NOT. MOB5AON .AND. .NOT. MOB5BON ) THEN

                CNTPSI = 1
                READ( LINE, *, IOSTAT=IOS ) PSI, CNTCOMBO

                CALL CHECK_READ_STATUS

C.............  Otherwise, skip the current line, which is the meat of the 
C               MOBILE input data
            ELSE
                SKIPLINE = .TRUE. 

            END IF

            RETURN

C--------------  SUBPROGRAM FORMAT  STATEMENTS   --------------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_EFDAT_LINE

        END SUBROUTINE RDEFDAT
