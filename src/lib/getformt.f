
        INTEGER FUNCTION GETFORMT( FDEV )

C***********************************************************************
C  function body starts at line 109
C
C  DESCRIPTION:
C      This function returns the format of the inventory file, assuming that
C      the source category is already known.
C
C  PRECONDITIONS REQUIRED:
C      Point source inventory file opened on unit FDEV
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C      Functions: I/O API functions, CHKINT, CHKREAL
C
C  REVISION  HISTORY:
C      Created 11/98 by M. Houyoux
C
C**************************************************************************
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
C***************************************************************************

C...........   MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  i/o api parameters

C...........   EXTERNAL FUNCTIONS:
        LOGICAL         CHKINT
        LOGICAL         CHKREAL
        INTEGER         INDEX1
        INTEGER         JUNIT
        INTEGER         LBLANK
        INTEGER         STR2INT
        REAL            STR2REAL

        EXTERNAL        CHKINT, CHKREAL, INDEX1, JUNIT, LBLANK, STR2INT, 
     &                  STR2REAL

C...........   SUBROUTINE ARGUMENTS
        INTEGER         FDEV   ! unit number of input point sources file

C...........   Other local variables
        INTEGER         CSS, I, K, L, L1  ! counters and indices

        INTEGER         AT     ! tmp mobile area type
        INTEGER         CYID   ! tmp county ID
        INTEGER         DAY    ! tmp day number
        INTEGER         IDEV   ! tmp unit number
        INTEGER         IDIU   ! tmp diurnal temporal prof code
        INTEGER         IOS    ! tmp i/o status
        INTEGER         IWEK   ! tmp weekly temporal prof code
        INTEGER         IREC   ! line counter
        INTEGER         FIP    ! tmp FIPS code
        INTEGER         FT     ! tmp mobile facility type
        INTEGER         MONTH  ! tmp month number
        INTEGER         SCC    ! tmp SCC 
        INTEGER         SIC    ! tmp SIC 
        INTEGER         SDATE  ! tmp starting date
        INTEGER         STID   ! tmp state ID
        INTEGER         UZONE  ! tmp UTM zone
        INTEGER         VMT    ! tmp VMT
        INTEGER         YY     ! tmp 2-digit year number

        REAL            EMIS   ! tmp emissions
        REAL            STKD   ! tmp stack diameter
        REAL            STKH   ! tmp stack height
        REAL            STKT   ! tmp stack temperature
        REAL            STKV   ! tmp velocity
        REAL            XLOC   ! tmp x location
        REAL            YLOC   ! tmp y location

        LOGICAL         CFLAG  ! for flagging alphabetic characters

        CHARACTER*5     CPOL 
        CHARACTER*300   MESG 
        CHARACTER*300   NAMTMP 
        CHARACTER*600   LINE

        CHARACTER*2  :: EPSTYPES( 9 ) = 
     &                ( / 'AC', 'AA', 'AB', 'AE', 'AF', 'GB',
     &                    'GN', 'LB', 'LN' / )

        CHARACTER*2  :: EMSTYPES( 3 ) =
     &                ( / 'AA', 'AD', 'DS' / )

        CHARACTER*16 :: PROGNAME = 'GETFORMT'  ! program name

C***********************************************************************
C   begin body of function GETFORMT

        GETFORMT = IMISS3

C.........  Loop through lines of file
        IREC = 0
        DO

            READ( FDEV, 93000, END=201, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 )
     &                 'Error', IOS,  'reading inventory file ' //
     &                 'at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C............. Determine if the file is list format
            L1 = INDEX( LINE, 'INVYEAR' )

C.............  Found INVYEAR packet, so list format
            IF( L1 .GT. 0 ) THEN
                GETFORMT = LSTFMT
                EXIT           ! To rewind and return
            END IF
 
C.............  Determine if the file contains header lines
            L1 = INDEX( LINE, CINVHDR )
            IF( L1 .EQ. 1 ) THEN

C................  Check if format is provided as a header entry
                L = INDEX( LINE, 'IDA' )
                IF( L .GT. 0 ) THEN
                    GETFORMT = IDAFMT
                    EXIT ! To end inner loop
                END IF

                L = INDEX( LINE, 'EMS-95' )
                IF( L .GT. 0 ) THEN
                    GETFORMT = EMSFMT
                    EXIT ! To end inner loop
                END IF

                L = INDEX( LINE, 'EPS2' )
                IF( L .GT. 0 ) THEN
                    GETFORMT = EPSFMT
                    EXIT ! To end inner loop
                END IF

C................  Make sure it is a point source file and not another type
C                  of IDA file. Approach is to read additional lines until
C                  there is not a '#' and check positions 61-69.
                DO 
                    READ( FDEV, 93000, END=201, IOSTAT=IOS ) LINE

                    IF ( IOS .GT. 0 ) THEN
                        WRITE( MESG, 94010 )
     &                         'Error', IOS, 'reading inventory file '//
     &                         'at line', IREC
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ENDIF

                    L1 = INDEX( LINE, CINVHDR )
                    L  = LEN_TRIM( LINE )

C.....................  Skip remaining header lines
                    IF( L1 .EQ. 1 ) THEN
                        CYCLE

C.....................  Skip blank lines
                    ELSE IF( L .EQ. 0 ) THEN
                        CYCLE

                    ELSE IF( CATEGORY .EQ. 'AREA' .AND.
     &                      L .GT. 16 ) THEN            ! Not a perfect test
                            
                        GETFORMT = IDAFMT
                        EXIT ! To end inner loop

                    ELSE IF( CATEGORY .EQ. 'MOBILE' .AND.
     &                      L .GT. 16 ) THEN            ! Not a perfect test
                            
                        GETFORMT = IDAFMT
                        EXIT ! To end inner loop

                    ELSE IF( CATEGORY .EQ. 'POINT' .AND.
     &                      L .GT. 250 ) THEN           ! Not a perfect test
                            
                        GETFORMT = IDAFMT
                        EXIT ! To end inner loop

                    ELSE
                        MESG = 'Input format error or IDA-formatted ' //
     &                         'file is not a ' // 
     &                         CATDESC // ' source file!'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                END DO

                EXIT ! To rewind and return

            END IF
            
C............. Determine if the file is EPS format.  Approach is to scan the
C              correct part of the field for emissions type, FIPS code, and
C              beginning period interval; and do a range checks
            IF( CATEGORY .EQ. 'POINT' ) THEN  
                I = INDEX1( LINE( 9:10 ), 9, EPSTYPES )
            ELSE
                I = INDEX1( LINE( 7:8 ), 9, EPSTYPES )
            END IF

            IF( I .GT. 0 ) THEN

                IF( CATEGORY .EQ. 'POINT' ) THEN  
                    FIP = STR2INT( LINE( 12:16 ) )  ! returns < 0 for bad value
                ELSE
                    FIP = STR2INT( LINE(  9:13 ) )
                END IF

                IF( FIP .GT. 1000 ) THEN

                    IF( CATEGORY .EQ. 'POINT' ) THEN  
                        SDATE = STR2INT( LINE( 60:67 ) )
                    ELSE
                        SDATE = STR2INT( LINE( 40:47 ) )
                    END IF

                    IF( SDATE .GT. 9999999 ) THEN

                        YY = STR2INT( LINE( 60:61 ) )

                        IF( YY .GT. 18 ) THEN
                            GETFORMT = EPSFMT
                            EXIT           ! To rewind and return
                        END IF
                    END IF
                END IF
            END IF
 
C............. Determine if the file is EMS-95 format.  Needs to return EMSFMT
C              for all types of EMS-95 formats.  Approach is to scan the
C              correct part of the fields for FIPS code for all files and do
C              a range check.  Then, we try various other fields for
C              different EMS-95 formats.

            STID = 0
            CYID = 0
            IF( CHKINT( LINE( 1:2 ) ) )
     &          STID = STR2INT( LINE( 1:2 ) )    ! returns < 0 for bad value
 
            IF( CHKINT( LINE( 3:5 ) ) )
     &          CYID = STR2INT( LINE( 3:5 ) )    ! returns < 0 for bad value

            IF( CATEGORY .EQ. 'POINT' ) THEN
 
                IF( STID .GT. 0 .AND. CYID .GT. 0 ) THEN

C.....................  Try device file
                    SIC = -99999
                    IF( CHKINT( LINE( 45:48 ) ) )
     &                  SIC = STR2INT( LINE( 45:48 ) )

C.....................  Check for valid SIC values from 0111 to 9999 or missing
C                       set to 0000 or -9
                    IF( SIC .GT. 111 .OR. 
     &                  SIC .EQ. 0   .OR. 
     &                  SIC .EQ. -9       ) THEN 

                        IWEK = -9
                        IF( CHKINT( LINE( 123:124 ) ) )
     &                      IWEK = STR2INT( LINE( 123:124 ) )

                        IDIU = -9
                        IF( CHKINT( LINE( 121:122 ) ) )
     &                      IDIU = STR2INT( LINE( 121:122 ) )

                        IF( IWEK .GT. -9 .AND. IDIU .GT. -9 ) THEN
                            GETFORMT = EMSFMT
                            EXIT           ! To rewind and return
                        ENDIF

                    ENDIF

C.....................  Try emission file
                    CSS  = LBLANK ( LINE( 57:61 ) )
                    CPOL = LINE   ( MIN(CSS+57,61):61 )  ! char string
                    L = LEN_TRIM( CPOL )

                    CFLAG = .FALSE.    ! Does the field for pollutant name have chars?
                    DO I = 1, L
                        IF( .NOT. CFLAG .AND. 
     &                      CPOL( I:I ) .GT. '9' ) CFLAG = .TRUE.
                    ENDDO
                
                    IF( CFLAG ) THEN

C.........................  Check that emission type is valid for EMS-95 files
                        I = INDEX1( LINE( 114:115 ), 3, EMSTYPES )

                        IF( I .GT. 0 ) THEN
                            EMIS = 0.
                            IF( CHKREAL( LINE( 88:100 ) ) )
     &                          EMIS = STR2REAL( LINE( 88:100 ) )

C.............................  Check for non-integer value or zero
                            IF( EMIS .NE. FLOAT( INT( EMIS ) ) .OR.
     &                          EMIS .EQ. 0. )  THEN
                                GETFORMT = EMSFMT
                                EXIT           ! To rewind and return
                            ENDIF
                        ENDIF
   
                    ENDIF

C.....................  Try facility file
                    IF( CHKINT( LINE( 43:44 ) ) )
     &                  UZONE = STR2INT( LINE( 43:44 ) )

                    IF( UZONE .GT. 0 .AND. UZONE .LE. 60 ) THEN 

                        XLOC = 0.
                        IF( CHKREAL( LINE( 25:33 ) ) ) 
     &                       XLOC = STR2REAL( LINE( 25:33 ) )

                        YLOC = 0.
                        IF( CHKREAL( LINE( 34:42 ) ) ) 
     &                      YLOC = STR2REAL( LINE( 34:42 ) )
    
                        IF( XLOC .GT. 0. .AND. YLOC .GT. 0. ) THEN
                            GETFORMT = EMSFMT
                            EXIT           ! To rewind and return
                        ENDIF
   
                    ENDIF

C.....................  Try stack file
                    STKD = -9.
                    IF( CHKREAL( LINE( 33:40 ) ) ) 
     &                  STKD = STR2REAL( LINE( 33:40 ) )

                    STKH = -9.
                    IF( CHKREAL( LINE( 41:47 ) ) ) 
     &                  STKH = STR2REAL( LINE( 41:47 ) )

                    STKT = -9.
                    IF( CHKREAL( LINE( 48:54 ) ) ) 
     &                  STKT = STR2REAL( LINE( 48:54 ) ) 

                    STKV = -9.
                    IF( CHKREAL( LINE( 55:61 ) ) ) 
     &                  STKV = STR2REAL( LINE( 55:61 ) )

                    IF( STKD .GT. -9. .AND. STKH .GT. -9. .AND.
     &                  STKT .GT. -9. .AND. STKV .GT. -9.       ) THEN
                        GETFORMT = EMSFMT
                        EXIT           ! To rewind and return
                    ENDIF

C.....................  Try process file
                    SCC = 0
                    IF( CHKINT( LINE( 57:64 ) ) ) 
     &                  SCC = STR2INT( LINE( 57:64 ) )

                    IF( SCC .GT. 9999999 ) THEN 
                        GETFORMT = EMSFMT
                        EXIT           ! To rewind and return
                    ENDIF

C.....................  Try day-specific emissions file
                    MONTH = 0
                    IF( CHKINT( LINE( 62:63 ) ) )
     &                  MONTH = STR2INT( LINE( 62:63 ) )

                    DAY   = 0
                    IF( CHKINT( LINE( 65:66 ) ) )
     &                  DAY   = STR2INT( LINE( 65:66 ) )

                    YY    = 0
                    IF( CHKINT( LINE( 68:69 ) ) )
     &                  YY    = STR2INT( LINE( 68:69 ) )

                    IF( MONTH .GT. 0 .AND. MONTH .LE. 12 .AND.
     &                  DAY   .GT. 0 .AND. DAY   .LE. 31 .AND.
     &                  YY    .GT. 18 ) THEN
                        GETFORMT = EMSFMT
                        EXIT           ! To rewind and return
                    END IF

                END IF

            ELSE IF ( CATEGORY .EQ. 'AREA' ) THEN

                CPOL = ADJUSTL( LINE( 21:25 ) )
                L = LEN_TRIM( CPOL ) ! Width of SCC field

                CFLAG = .FALSE.    ! Does the field for pollutant name have chars?
                DO I = 1, L
                    IF( .NOT. CFLAG .AND. 
     &                  CPOL( I:I ) .GT. '9' ) CFLAG = .TRUE.
                END DO
                
                IF( CFLAG ) THEN

                    I = INDEX1( LINE( 95:96 ), 3, EMSTYPES )

                    IF( I .GT. 0 ) THEN
                        EMIS = 0.
                        IF( CHKREAL( LINE( 52:64 ) ) )
     &                      EMIS = STR2REAL( LINE( 52:64 ) )

                        IF( EMIS .NE. FLOAT( INT( EMIS ) ) )  THEN
                            GETFORMT = EMSFMT
                            EXIT           ! To rewind and return
                        END IF
                    END IF
   
                END IF

            ELSE IF ( CATEGORY .EQ. 'MOBILE' ) THEN
                AT = -1
                IF( CHKINT( LINE( 6:6 ) ) ) AT = STR2INT( LINE( 6:6 ) )

                FT = -1
                IF( CHKINT( LINE( 7:10 ) ) ) 
     &              FT = STR2INT( LINE( 7:10 ) )

                VMT = -1.
                IF( CHKREAL( LINE( 11:18 ) ) ) 
     &              VMT = STR2REAL( LINE( 11:18 ) ) 

                IF( AT .GE. 0 .AND. FT .GE. 0 .AND. VMT .GE. 0. ) THEN
                    GETFORMT = EMSFMT
                    EXIT           ! To rewind and return
                ENDIF

            END IF

C.............  Attempt to open LINE as file name
            NAMTMP = LINE( 1:INDEX( LINE,' ' ) )
            IDEV = JUNIT()
            OPEN( IDEV, ERR=99, FILE=NAMTMP, STATUS='OLD' )

C.............  If get here, then file has been opened, so this must be a
C               list file
            GETFORMT = LSTFMT
            CLOSE( IDEV )

            EXIT              ! To rewind and return

C.............  Unsuccessful opening of LINE as file name
99          CONTINUE
            CLOSE( IDEV )

        END DO     ! To head of read loop

201     CONTINUE

999     REWIND( FDEV )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END FUNCTION GETFORMT

