
       SUBROUTINE RDIDAMB ( FDEV, NRAWIN, NRAWOUT, EFLAG, NDROP, VDROP )

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
C File: @(#)rdsprof.f	1.1
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

        USE MODSOURC            ! correct?

        IMPLICIT NONE

C...........   Include files

        INCLUDE 'EMCNST3.EXT'   ! valid country codes
	
C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         INDEX1
        INTEGER         STR2INT
        REAL            STR2REAL
	INTEGER         TRIMLEN
        CHARACTER*2     CRLF

        EXTERNAL        CRLF, INDEX1, PADZERO, STR2INT, STR2REAL, TRIMLEN

C...........   Subroutine arguments

        INTEGER     , INTENT    (IN) :: FDEV        ! file unit number
        INTEGER     , INTENT    (IN) :: NRAWIN      ! toal raw record count
        INTEGER     , INTENT   (OUT) :: NRAWOUT     ! valid raw record count
        INTEGER     , INTENT   (OUT) :: NDROP       ! number of records dropped
        REAL        , INTENT   (OUT) :: VDROP       ! sum of VMT dropped
        LOGICAL     , INTENT   (OUT) :: EFLAG       ! error flag

C...........   Local variables

        INTEGER       CYID           ! tmp county code
        INTEGER       ICC            ! tmp country code
        INTEGER       IDROP          ! dropped raw record counter
        INTEGER       IOS            ! I/O status
        INTEGER       IRAWOUT        ! valid raw record counter
        INTEGER       IREC           ! record counter
        INTEGER       IYEAR          ! data year
        INTEGER       SPEED          ! tmp speed, integer
        INTEGER       STID           ! tmp state code
        REAL          VMTR           ! tmp vehicle miles traveled
        LOGICAL       INVALID        ! tmp error flag for current record
        CHARACTER*10  CNTRY          ! tmp country name
        CHARACTER*35  LINE           ! read buffer for a line
        CHARACTER*300 MESG           ! text for M3EXIT()
        CHARACTER*16 :: PROGNAME = 'RDIDAMB'   ! program name
        CHARACTER*10  SCC            ! tmp AMS source classification code
        CHARACTER*2   SPEEDCHR       ! tmp speed, character
        CHARACTER*5   VTYPE          ! tmp vehicle type


C***********************************************************************
C   Begin body of subroutine RDIDAMB

        IREC    = 0
	IDROP   = 0
        IRAWOUT = 0
        VDROP   = 0
	EFLAG   = .FALSE.

        DO

        READ( FDEV, 93000, END=12, IOSTAT=IOS ) LINE
	
        INVALID = .FALSE.
        IREC = IREC + 1
	
        IF ( LINE .EQ. ' ' ) CYCLE      ! skip if line is blank
	
        IF ( LINE(1:1) .EQ. '#' ) THEN           ! determine if line is a header
	
	   IF ( LINE(2:8) .EQ. 'COUNTRY' ) THEN  ! read in 'country-name'
                CNTRY   = LINE(9:35)
                CNTRY   = ADJUSTL( CNTRY )
                ICC     = INDEX1( CNTRY, MXCNTRY3, CNTRYNM3 )
	     
                IF ( ICC .EQ. 0 ) THEN
                     WRITE( MESG, 94010 )
     &	             'Invalid country name encountered on line ', IREC
                     CALL M3MESG( MESG )
                     INVALID = .TRUE.
                END IF 
	   
           ELSEIF ( LINE(2:5) .EQ. 'YEAR' ) THEN ! read in 'datayear'
                IYEAR = STR2INT( LINE(6:35) )
                IF ( IYEAR .LT. 1971 ) THEN
                     WRITE( MESG, 94010 )
     &	             'Invalid year encountered on line ', IREC
                     CALL M3MESG( MESG )
                     INVALID = .TRUE.
                END IF

           ELSE
                CYCLE   ! cycle if header line does not contain 'country-name'
                        ! or 'datayear' information
           END IF
	   
        ELSE

C...........    Fill temporary fields using current record line

        STID  = STR2INT ( LINE(1:2) )
	CYID  = STR2INT ( LINE(3:5) )
        SCC   = ADJUSTR ( LINE(6:15) )
        VMTR  = STR2REAL ( LINE(16:28) )
        VTYPE = ADJUSTL ( LINE(29:33) )
        SPEEDCHR = ADJUSTL ( LINE(34:35) )

C...........    Determine if vehicle type is valid
	
        IF ( VTYPE .EQ. ' ' ) THEN   ! determine if vehicle type is blank
             EFLAG = .TRUE.
	     INVALID = .TRUE.
             WRITE( MESG, 94010 )
     &	     'Invalid vehicle type encountered on line ', IREC
             CALL M3MESG( MESG )
        END IF

C...........    Determine if speed is valid

        IF ( SPEEDCHR .EQ. ' ' ) THEN   ! determine if speed is blank
             EFLAG = .TRUE.
	     INVALID = .TRUE.
             WRITE( MESG, 94010 )
     &	     'Invalid (blank) speed encountered on line ', IREC
             CALL M3MESG( MESG )
        END IF

        SPEED = STR2INT( SPEEDCHR )

        IF ( SPEED .EQ. 0 ) THEN   ! determine if speed=0
             EFLAG = .TRUE.
	     INVALID = .TRUE.
             WRITE( MESG, 94010 )
     &	     'Invalid (zero) speed encountered on line ', IREC
             CALL M3MESG( MESG )
        END IF

C...........    Determine if SCC is valid


        IF ( SCC .EQ. ' ' ) THEN   ! determine if SCC is 10 char. wide
             EFLAG = .TRUE.
	     INVALID = .TRUE.
             WRITE( MESG, 94010 )
     &	     'SCC not 10 characters wide on line ', IREC
             CALL M3MESG( MESG )
        END IF

        IF ( INVALID ) THEN
	
             IDROP = IDROP + 1
             VDROP = VDROP + VMTR
	     
        ELSE
	
             IRAWOUT = IRAWOUT + 1

             IF ( IRAWOUT .LE. NRAWIN ) THEN

             IFIPA  ( IRAWOUT )   = ICC*100000 + STID*1000 + CYID
             CSCCA  ( IRAWOUT )   = SCC
             VMTA   ( IRAWOUT )   = VMTR
             CVTYPEA( IRAWOUT )   = VTYPE
             ISPEEDA( IRAWOUT )   = SPEED
             INVYRA ( IRAWOUT )   = IYEAR
	     
             CALL PADZERO( CSCCA ( IRAWOUT ) )

             END IF
	     
        END IF
	
	END IF

        END DO
            
12      CONTINUE    ! end of read on input file

        NRAWOUT = IRAWOUT
	NDROP   = IDROP

        IF( NRAWOUT .GT. NRAWIN ) THEN  ! Check for memory overflow

            WRITE( MESG, 94010 )
     &        'INTERNAL ERROR: Number of valid raw records ' //
     &        'encountered: ', NRAWOUT, CRLF() // BLANK5 //
     &        'Maximum number of raw records allowed: ', NRAWIN

            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )      

        ENDIF
	
C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END
