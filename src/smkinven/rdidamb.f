
       SUBROUTINE RDIDAMB ( FDEV, NRAWIN, MXIPOL, INVPNAM, NRAWOUT, 
     &                      EFLAG, NDROP, VDROP )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      
C 
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 1998, MCNC--North Carolina Supercomputing Center
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

C...........   Modules for public variables
C.........  This module contains the inventory arrays
        USE MODSOURC

        IMPLICIT NONE

C...........   Include files

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters
	
C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         INDEX1
        INTEGER         STR2INT
        REAL            STR2REAL
        CHARACTER*2     CRLF

        EXTERNAL        CRLF, INDEX1, PADZERO, STR2INT, STR2REAL

C...........   Subroutine arguments. Note that number and amount of dropped
C              VMT is initialied in calling program.

        INTEGER     , INTENT (IN) :: FDEV        ! file unit number
        INTEGER     , INTENT (IN) :: NRAWIN      ! toal raw record count
        INTEGER     , INTENT (IN) :: MXIPOL      ! max no of inventory pols
        CHARACTER(*), INTENT (IN) :: INVPNAM( MXIPOL ) ! inv pol names + 'VMT'
        INTEGER     , INTENT(OUT) :: NRAWOUT     ! valid raw record count
        INTEGER     , INTENT(OUT) :: NDROP       ! number of records dropped
        REAL        , INTENT(OUT) :: VDROP       ! sum of VMT dropped
        LOGICAL     , INTENT(OUT) :: EFLAG       ! error flag

C...........   Local variables

        INTEGER       I              ! indices

        INTEGER       CYID           ! tmp county code
        INTEGER       FIP            ! tmp country/state/county
        INTEGER       ICC            ! tmp country code
        INTEGER       IDROP          ! local dropped raw record counter
        INTEGER       IOS            ! I/O status
        INTEGER       IRAWOUT        ! valid raw record counter
        INTEGER       IREC           ! record counter
        INTEGER       IYEAR          ! data year
        INTEGER       SCCLEN         ! length of SCC string 
        INTEGER       ISPEED         ! tmp speed, integer
        INTEGER       STID           ! tmp state code
        INTEGER       TPF            ! tmp temporal adjustments setting

        REAL          VMTR           ! tmp vehicle miles traveled (10^6 VMT/yr)

        LOGICAL       INVALID        ! tmp error flag for current record
        LOGICAL       INVCNTRY       ! tmp error flag for country name
        LOGICAL       INVYEAR        ! tmp error flag for year

        CHARACTER*2   SPEEDCHR       ! tmp speed, character
        CHARACTER*10  CNTRY          ! tmp country name
        CHARACTER*35  LINE           ! read buffer for a line
        CHARACTER*300 MESG           ! message buffer

        CHARACTER(LEN=POLLEN3) CCOD  ! character pollutant index to INVPNAM
        CHARACTER(LEN=FIPLEN3) CFIP  ! character FIPS code
        CHARACTER(LEN=LNKLEN3) CLNK  ! tmp link ID 
        CHARACTER(LEN=SCCLEN3) TSCC  ! tmp source classification code
        CHARACTER(LEN=VTPLEN3) VTYPE ! tmp vehicle type

        CHARACTER*16 :: PROGNAME = 'RDIDAMB'   ! program name

C***********************************************************************
C   Begin body of subroutine RDIDAMB

C.........  Find 'VMT' in "pollutants" names list
        I =  INDEX1( 'VMT', MXIPOL, INVPNAM )

C.........  Internal error if it has not been put in the pollutant names array
        IF( I .LE. 0 ) THEN

            MESG = 'INTERNAL ERROR: "VMT" has not been added to ' //
     &             'pollutant names list.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C.........  Store for use later
        ELSE
            WRITE( CCOD,94125 ) I

        END IF

C.........  Initialize before loop
        IREC     = 0
	IDROP    = 0
        IRAWOUT  = 0
	EFLAG    = .FALSE.
        INVCNTRY = .FALSE.
        TPF      = MTPRFAC * WTPRFAC
        CLNK     = ' '

C.........  Make sure the file is at the beginning
        REWIND( FDEV )

C.........  Loop through lines of mobile inventory file
        DO

            READ( FDEV, 93000, END=12, IOSTAT=IOS ) LINE
	
            INVALID = .FALSE.
            IREC = IREC + 1
 
            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &              'I/O error', IOS, 'reading VMT inventory '//
     &              'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
	
            IF ( LINE .EQ. ' ' ) CYCLE      ! skip if line is blank
	
            IF ( LINE(1:1) .EQ. '#' ) THEN  ! determine if line is a header

	        IF ( LINE(2:8) .EQ. 'COUNTRY' ) THEN  ! read in 'country-name'
                    CNTRY   = LINE(9:35)
                    CNTRY   = ADJUSTL( CNTRY )
                    ICC     = INDEX1( CNTRY, MXCNTRY3, CNTRYNM3 )
	     
                    IF ( ICC .EQ. 0 ) THEN
                        WRITE( MESG, 94010 ) 'Invalid country name ' //
     &	                       'encountered on line ', IREC
                        CALL M3MESG( MESG )
                        INVCNTRY = .TRUE.
                    ELSE
                        ICC      = CNTRYCD3( ICC )
                        INVCNTRY = .FALSE.
                    END IF 
	   
                ELSEIF ( LINE(2:5) .EQ. 'YEAR' ) THEN ! read in 'datayear'
                    IYEAR = STR2INT( LINE(6:35) )
                    IF ( IYEAR .LT. 1971 ) THEN
                        WRITE( MESG, 94010 )
     &	                   'Invalid year encountered on line ', IREC
                        CALL M3MESG( MESG )
                        INVYEAR = .TRUE.
                    ELSE
                        INVYEAR = .FALSE.
                    END IF

                ELSE
                    CYCLE   ! cycle if header line does not contain 'country-
                            !    name' or 'datayear' information
                END IF      ! on type of header

c           IF ( ICC .EQ. 0 .OR. IYEAR .LT. 1971 ) THEN ! set error flag to
c                INVCORY = .TRUE.   ! true if either ICC or IYEAR is invalid
c           ELSE
c                INVCORY = .FALSE.
c           END IF
	   
            ELSE   ! Continue for not header

                IF ( INVCNTRY .OR. INVYEAR ) THEN   ! if invalid country code or
                    EFLAG = .TRUE.       ! invalid year has been encountered, 
                    INVALID = .TRUE.     ! set EFLAG=.TRUE and INVALID = .TRUE.
                    WRITE( MESG, 94010 ) ! for all associated records
     &                  'ERROR: Record being dropped due to invalid ' //
     &                  'country code ' // CRLF() // BLANK5 //
     &                  'and/or invalid inventory year'
                    CALL M3MESG( MESG )       
                END IF

C.................    Fill temporary fields using current record line

                STID  = STR2INT ( LINE(1:2) )
	        CYID  = STR2INT ( LINE(3:5) )
                TSCC  = ADJUSTL ( LINE(6:15) )
                VMTR  = STR2REAL ( LINE(16:28) )
                VTYPE = ADJUSTL ( LINE(29:33) )
                SPEEDCHR = ADJUSTL ( LINE(34:35) )

C.................    Determine if vehicle type is valid
	
                IF ( VTYPE .EQ. ' ' ) THEN   ! determine if vehicle type is blank
                     EFLAG = .TRUE.
	             INVALID = .TRUE.
                     WRITE( MESG, 94010 )
     &	             'Invalid (blank) vehicle type encountered on ' //
     &               'line', IREC
                     CALL M3MESG( MESG )
                END IF
            
C.................    Determine if speed is valid

                IF ( SPEEDCHR .EQ. ' ' ) THEN   ! determine if speed is blank
                     EFLAG = .TRUE.
	             INVALID = .TRUE.
                     WRITE( MESG, 94010 )
     &	             'Invalid (blank) speed encountered on line ', IREC
                     CALL M3MESG( MESG )
                END IF

                ISPEED = STR2INT( SPEEDCHR )

                IF ( ISPEED .EQ. 0 ) THEN   ! determine if speed=0
                     EFLAG = .TRUE.
	             INVALID = .TRUE.
                     WRITE( MESG, 94010 )
     &	             'Invalid (zero) speed encountered on line ', IREC
                     CALL M3MESG( MESG )
                END IF

C.................    Determine if SCC is valid

                SCCLEN = LEN_TRIM( TSCC )

                IF ( SCCLEN .NE. 10 ) THEN   ! check if SCC is 10 char. wide
                    EFLAG = .TRUE.
	            INVALID = .TRUE.
                    WRITE( MESG, 94010 )
     &	                   'SCC not 10 characters wide on line ', IREC
                    CALL M3MESG( MESG )
                END IF

C.................  Make adjustments to pad with zeros, if needed
                FIP = ICC*100000+ STID*1000+ CYID
                WRITE( CFIP,94120 ) FIP
                CALL PADZERO( TSCC )

                IF ( INVALID ) THEN
	
                    IDROP = IDROP + 1
                    VDROP = VDROP + VMTR
	     
                ELSE
	
                    IRAWOUT = IRAWOUT + 1

                    IF ( IRAWOUT .LE. NRAWIN ) THEN

                        IFIPA  ( IRAWOUT ) = FIP
                        CSCCA  ( IRAWOUT ) = TSCC
                        CVTYPEA( IRAWOUT ) = VTYPE
                        CLINKA ( IRAWOUT ) = CLNK
                        SPEEDA ( IRAWOUT ) = REAL( ISPEED )
                        TPFLGA ( IRAWOUT ) = TPF
                        INVYRA ( IRAWOUT ) = IYEAR
                        XLOC1A ( IRAWOUT ) = AMISS3
                        YLOC1A ( IRAWOUT ) = AMISS3
                        XLOC2A ( IRAWOUT ) = AMISS3
                        YLOC2A ( IRAWOUT ) = AMISS3
                        POLVLA ( IRAWOUT,1 ) = VMTR * 1000000

                        CALL BLDCSRC( CFIP, TSCC, VTYPE, CLNK, ' ', ' ', 
     &                                ' ', CCOD, CSOURCA( IRAWOUT ) )

                    END IF
	     
                END IF   ! if valid record or not
	
            END IF       ! if header or not

        END DO           ! end of loop for reading input file
              
12      CONTINUE         ! end of read on input file

        NRAWOUT = IRAWOUT
	NDROP   = NDROP + IDROP

        IF( NRAWOUT .GT. NRAWIN ) THEN  ! Check for memory overflow

            WRITE( MESG, 94010 )
     &        'INTERNAL ERROR: Number of valid raw records ' //
     &        'encountered: ', NRAWOUT, CRLF() // BLANK5 //
     &        'Maximum number of raw records allowed: ', NRAWIN

            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )      

        END IF
	
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94120   FORMAT( I6.6 )

94125   FORMAT( I5 )

        END SUBROUTINE RDIDAMB
