
        SUBROUTINE RDEFTMPR( FDEV, READALL )
   
C***********************************************************************
C  subroutine RDEFTMPR body starts at line < >
C
C  DESCRIPTION:
C      Read emission-factor/temperature file from the preprocessing program
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
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C 
C See file COPYRIGHT for conditions of use.
C 
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C
C smoke@emc.mcnc.org
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C****************************************************************************
 
C...........   MODULES for public variables
C...........   This module contains emission factor tables and related
        USE MODEMFAC

C...........   This module contains the information about the source category
        USE MODINFO

C...........   This module is the derived meteorology data for emission factors
        USE MODMET

        IMPLICIT NONE

C...........   INCLUDES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'FLTERR.EXT'    !  error filter statement functions

C...........   EXTERNAL FUNCTIONS
        CHARACTER*2   CRLF
        INTEGER       FIND1
        INTEGER       INDEX1
        INTEGER       GETFLINE
        INTEGER       STR2INT
        REAL          STR2REAL

        EXTERNAL      CRLF, FIND1, INDEX1, GETFLINE, STR2INT,
     &                STR2REAL

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV    ! Min/max-T/EF file 
        LOGICAL     , INTENT (IN) :: READALL ! true: read data; false: PSILIST

C...........   Local parameters
        INTEGER, PARAMETER :: MXCOL = 5   ! maximum number of input columns

C...........   Local variables
        INTEGER         I, J1, J2, K, K1, K2, N

        INTEGER         IOS    ! i/o status
        INTEGER         IREC   ! record counter
        INTEGER         LPSI   ! previous iteration PSI
        INTEGER, SAVE:: MXMMTEF! max no. min/max tmpr/EF combinations
        INTEGER         PSI    ! tmp parameter scheme index
        INTEGER         PSI1   ! tmp parameter scheme index
        INTEGER         PSI2   ! tmp parameter scheme index
        INTEGER         TMMI   ! min/max temperature index

        REAL            CHKTMN ! reference min tmpr from master list
        REAL            CHKTMX ! reference max tmpr from master list
        REAL            TMAX   ! tmp maximum temperature
        REAL            TMIN   ! tmp minimum temperature

        LOGICAL      :: EFLAG    = .FALSE.  ! true: error found
        LOGICAL, SAVE:: FIRSTIME = .TRUE.   ! true: first time routine called
        LOGICAL      :: LMNFLG   = .FALSE.  ! true: min tmpr not okay
        LOGICAL      :: LMXFLG   = .FALSE.  ! true: max tmpr not okay
        LOGICAL      :: SKIPLINE = .FALSE.  ! true: skip line during read
        LOGICAL      :: WFLAG    = .TRUE.   ! true: write out warning

        CHARACTER*20    SEGMENT( MXCOL ) ! for reading min/max-T/EF file

        CHARACTER*300          LINE      ! message buffer
        CHARACTER*300          MESG      ! message buffer

        CHARACTER(LEN=IOVLEN3) ACT       ! activity
        CHARACTER(LEN=IOVLEN3) LACT      ! previous iteration activity

        CHARACTER*16 :: PROGNAME = 'RDEFTMPR' ! program name

C***********************************************************************
C   begin body of subroutine RDEFTMPR

C.........  Determine maximum memory usage for arrays
        IF( FIRSTIME ) THEN
            MXMMTEF = GETFLINE( FDEV, 'EF-ref/temperature file' )
            FIRSTIME = .FALSE.
        END IF

C.........   When whole file is not to be read...
        IF( .NOT. READALL ) THEN

C.............  Allocate memory for the number of PSIs per activity and for
C               the PSIs.  The number used for the PSIs is way too big, but
C               much simpler so that the file doesn't need to be read again.
            ALLOCATE( NPSI( NIACT ),STAT=IOS )
            CALL CHECKMEM( IOS, 'NPSI', PROGNAME )
            ALLOCATE( PSILIST( MXMMTEF, NIACT ),STAT=IOS )
            CALL CHECKMEM( IOS, 'PSILIST', PROGNAME )

C.............  Determine unique PSIs for each activity from this file to create
C               PSILIST array
            I     = 0
            IREC  = 0
            LPSI  = -9
            LACT  = ' '
            N = 0
            DO 

        	READ( FDEV, 93000, END=55, IOSTAT=IOS ) LINE
        	IREC = IREC + 1

C.................  Check status of read and skip line if bad status
        	CALL CHECK_READ_STATUS
        	IF( SKIPLINE ) CYCLE

C.................  Parse line into sections
        	CALL PARSLINE( LINE, MXCOL, SEGMENT )

C.................  Extract activity from record and search for it in main list
        	ACT = SEGMENT( 5 )
        	K   = INDEX1( ACT, NIACT, ACTVTY )

C.................  Skip over any entries that are not for the activity of 
C                   interest
        	IF( K .LE. 0 ) CYCLE

C.................  Extract PSI from line
        	PSI  = STR2INT ( SEGMENT( 1 ) )

        	IF( ACT .NE. LACT ) THEN
                    N = N + 1
                    NPSI( N ) = 0
        	END IF

        	IF( PSI .NE. LPSI ) THEN
                    NPSI( N ) = NPSI( N ) + 1
                    PSILIST( NPSI( N ), N ) = PSI
                END IF

                LACT = ACT
                LPSI = PSI

            END DO
55          CONTINUE  ! Exit from endless loop

            REWIND( FDEV )

            RETURN   ! That is all if not reading all data

        END IF

C.........  This section for full read (READALL=TRUE) ...

C.........  Allocate memory for the EF-ref/temperature data arrays
        ALLOCATE( MMTEFPSI( MXMMTEF ),STAT=IOS )
        CALL CHECKMEM( IOS, 'MMTEFPSI', PROGNAME )
        ALLOCATE( MMTEFIDX( MXMMTEF ),STAT=IOS )
        CALL CHECKMEM( IOS, 'MMTEFIDX', PROGNAME )
        ALLOCATE( MMTEFPTR( NPSIALL ),STAT=IOS )
        CALL CHECKMEM( IOS, 'MMTEFPTR', PROGNAME )
        ALLOCATE( MMTEFNUM( NPSIALL ),STAT=IOS )
        CALL CHECKMEM( IOS, 'MMTEFNUM', PROGNAME )
        ALLOCATE( MMTEFEND( NPSIALL ),STAT=IOS )
        CALL CHECKMEM( IOS, 'MMTEFEND', PROGNAME )
        MMTEFNUM = 0    ! array
        MMTEFEND = 0    ! array

C.........  Loop through file until we are out of lines
        I     = 0
        IREC  = 0
        DO 

            READ( FDEV, 93000, END=101, IOSTAT=IOS ) LINE
            IREC = IREC + 1

C.............  Check status of read and skip line if bad status
            CALL CHECK_READ_STATUS
            IF( SKIPLINE ) CYCLE

C.............  Parse line into sections
            CALL PARSLINE( LINE, MXCOL, SEGMENT )

C.............  Extract activity from record and search for it in main list
            ACT = SEGMENT( 5 )
            K   = INDEX1( ACT, NIACT, ACTVTY )

C.............  Skip over any entries that are not for the activity of interest
            IF( K .LE. 0 ) CYCLE

C.............  Extract other entries from file
            PSI  = STR2INT ( SEGMENT( 1 ) )
            TMIN = STR2REAL( SEGMENT( 2 ) )
            TMAX = STR2REAL( SEGMENT( 3 ) )
            TMMI = STR2INT ( SEGMENT( 4 ) )

C.............  Search for PSI is in master list
            N = FIND1( PSI, NPSIALL, PSIALL )

C.............  If PSI is not found, then report and set as an error
            IF( N .LE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: PSI', PSI, 
     &                 'in EF-ref/temperature file is not in list of'//
     &                  CRLF() // BLANK10 // 'valid PSIs.  Recreate '//
     &                 'EF-ref/temperature file.'
                CALL M3MESG( MESG )
                CYCLE

            END IF

C.............  Check min/max temperature combo against list of valid ones
C               to ensure using the same min/max bounds
            LMNFLG = .TRUE.
            LMXFLG = .TRUE.
            IF( TMMI .LE. NVLDTMM ) THEN
                CHKTMN = VLDTMIN( TMMI )
                CHKTMX = VLDTMAX( TMMI )                
                LMNFLG = FLTERR( TMIN, CHKTMN )
                LMXFLG = FLTERR( TMAX, CHKTMX )
            END IF

C.............  If min/max temperatures inconsistent, then report and set as 
C               error
            IF( LMNFLG .OR. LMXFLG ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94100 ) 'ERROR: min/max temperatures', TMIN,
     &                 '/', TMAX, 'in EF-ref/temperature file '//
     &                 CRLF() // BLANK10 // 'are inconsistent with '//
     &                 'acceptable values set with SMK_MINT_MIN,' //
     &                 CRLF() // BLANK10 // ' SMK_MAXT_MIN, etc.'
                CALL M3MSG2( MESG )

            END IF

            I = I + 1

            IF( I .LE. MXMMTEF ) THEN

                MMTEFPSI( I ) = PSI
                MMTEFIDX( I ) = TMMI

            END IF

        END DO

101     CONTINUE   ! exit from read loop
        NMMTEF = I 

        IF( NMMTEF .GT. MXMMTEF ) THEN

            EFLAG = .TRUE.
            WRITE( MESG, 94010 )
     &             'INTERNAL ERROR: dimension for storing ' //
     &             'EF-ref/temperature table was', MXMMTEF, 
     &             CRLF() // BLANK10 // 'but actually needed', NMMTEF           
            CALL M3MSG2( MESG )

        END IF

        IF( EFLAG ) THEN
            MESG = 'ERROR reading min/max temperatures with EFs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Note that it is okay to not have minimum and maximum temperatures
C           for a PSI, because perhaps this PSI is for sources that are outside
C           the grid.  In this circumstance, however, need to make sure that 
C           the PSI won't be updated!

C.........  Loop through PSIs to initialize pointer to the position in the 
C           temperature index array of the first occurence of each PSI
C.........  NOTE - FIND1 won't work because it is not guaranteed to find the
C           first position
        DO I = 1, NPSIALL

            PSI = PSIALL( I )

            DO K = 1, NMMTEF
                IF( PSI .EQ. MMTEFPSI( K ) ) EXIT
            END DO

C.............  Duplicate the behavior of FIND1 for PSI not found
            IF( K .GT. NMMTEF ) K = -1

C.............  Store position or -1 to indicate position or that PSI is not
C               in the grid or that it was not in MPLIST file because it
C               is a contributing PSI
            MMTEFPTR( I ) = K

        END DO

C.........  Set the number of temperature indices per PSI.  Loop through
C           the PSIs from the EF/temperature file though, because other
C           PSIs will not have temperatures associated with them.
        DO N = 1, NIACT
            DO I = 1, NPSI( N ) - 1

                PSI1 = PSILIST( I  ,N )
                PSI2 = PSILIST( I+1,N )
                J1 = FIND1( PSI1, NPSIALL, PSIALL )
                J2 = FIND1( PSI2, NPSIALL, PSIALL )

                IF( J1 .LE. 0 .OR. J2 .LE. 0 ) THEN
                    MESG = 'INTERNAL ERROR 1'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                K1 = MMTEFPTR( J1 )
                K2 = MMTEFPTR( J2 )
                IF( K1 .LE. 0 .OR. K2 .LE. 0 ) THEN
                    MESG = 'INTERNAL ERROR 2'
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                MMTEFNUM( J1 ) = K2 - K1
                MMTEFEND( J1 ) = K2 - 1

            END DO

            PSI1 = PSILIST( NPSI( N ),N )
            J1 = FIND1( PSI1, NPSIALL, PSIALL )
            K1 = MMTEFPTR( J1 )
            K2 = NMMTEF
            MMTEFNUM( J1 ) = K2 - K1 + 1
            MMTEFEND( J1 ) = K2 

        END DO

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I9, :, 1X ) )

94100   FORMAT( 10( A, :, F9.2, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C..............  This internal subprogram checks the read status and reports
C                and error if there was one.  It sets EFLAG and skips to 
C                next record

            SUBROUTINE CHECK_READ_STATUS

C.............................................................................

            SKIPLINE = .FALSE.
            IF( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'I/O Error ', IOS,
     &              'reading EF_REF/TEMPERATURE file at line', IREC
                CALL M3MESG( MESG )
                SKIPLINE = .TRUE.

            END IF

C--------------  SUBPROGRAM FORMAT  STATEMENTS   --------------------------

C...........   Internal buffering formats............ 94xxx

94010       FORMAT( 10( A, :, I8, :, 1X ) )

            END SUBROUTINE CHECK_READ_STATUS

        END SUBROUTINE RDEFTMPR
