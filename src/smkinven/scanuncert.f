
        SUBROUTINE SCANUNCERT( FDEV )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C    The SCANUNCERT routine will be used to scan the uncertainty file and make
C    global program settings.  It determines the number of each packet type,
C    the number of packet entries and the maximum number of packet entries.
C
C  PRECONDITIONS REQUIRED:
C    Uncertainty file is opened
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C    Subroutines:  I/O API subroutines, RMCOMMNT, UPCASE, PARSLINE
C                     PRCLINUC
C    Functions:  I/O API functions, CRLF, BLKORCMT, GETFLINE, GETNLIST
C
C  REVISION  HISTORY:
C    Created 9/2001 by A. Holland
C
C***********************************************************************
C  
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C  
C COPYRIGHT (C) 2001, MCNC--North Carolina Supercomputing Center
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

C.........  This module contains uncertainty-specific settings
        USE MODUNCERT


        IMPLICIT NONE
	
C...........  INCLUDES
    	INCLUDE 'EMCNST3.EXT'	


C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER*2   CRLF
        LOGICAL       BLKORCMT
        INTEGER       GETFLINE
        INTEGER       GETNLIST

        EXTERNAL   CRLF, BLKORCMT, GETFLINE, GETNLIST

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT (IN) :: FDEV   ! report configuration file unit no.

C...........   Line parsing array
        CHARACTER*300 SEGMENT( 100 )
 
C...........   Other local variables
        INTEGER         I, L    ! counters and indices

        INTEGER          ICNT    !  no. non-blank/non-comment lines
        INTEGER          IOS     !  i/o status
        INTEGER          IREC    !  record counter
        INTEGER       :: N = 1   !  no. segments in line

        LOGICAL       :: EFLAG   = .FALSE. !  true: error found

        CHARACTER*300    BUFFER            !  line work buffer
        CHARACTER*300    LINE              !  line input buffer
        CHARACTER*300    MESG              !  message buffer

        CHARACTER*16 :: PROGNAME = 'SCANUNCERT' ! program name

C***********************************************************************
C   begin body of subroutine SCANUNCERT

        MESG = 'Scanning uncertainty file...'
        CALL M3MSG2( MESG )         
       
C.........  Get the number of lines in the file
        NLINE_UC = GETFLINE( FDEV, 'Uncertainty file' )

C.........  Read lines of file to extract settings
        IREC = 0
        DO I = 1, NLINE_UC

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 
     &              'I/O error', IOS, 
     &              'reading uncertainty file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

C.............  Skip blank lines and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Screen for appended comments and remove them
            CALL RMCOMMNT( '##', LINE )

C.............  Left-justify and convert line to upper case
            LINE = ADJUSTL( LINE )
            BUFFER = LINE
            CALL UPCASE( BUFFER )

C.............  Initialize segment from previous iteration
            SEGMENT( 1:N ) = ' '

C.............  Parse line into segments
            L = LEN_TRIM( BUFFER )
            N = GETNLIST( L, BUFFER )
            CALL PARSLINE( BUFFER, N, SEGMENT )

C.............  Interpret line of code.  Set global variables in MODUNCERT.
            CALL PRCLINUC( IREC, N, LINE, SEGMENT )
	    
C.............  Get number of factor assignment entries
	    IF( INFAPKT ) THEN
	        IF( .NOT. FASTART ) FPKTENT = FAENTN
	    END IF	
	    
C.............  Maximum number of emission factor values and probabilities
	    MXEMPDAT  = MAX( MXEMPDAT, EMPENTN, 1 )
	    
C.............  Maximum number of parameters
	    MXPARDAT  = MAX( MXPARDAT, PAR_%NUMP )	    	        

C.............  Keep count of non-blank/non-comment lines
            ICNT = ICNT + 1

        END DO    ! End read loop

C.........  Rewind file
        REWIND( FDEV )

C.........  Store the number of packets of each type
        NFPCKT  = PKTCOUNT( FA_IDX )    ! Number of factor assignment pkts.
        NEPCKT  = PKTCOUNT( EMP_IDX )   ! Number of empirical pkts.
	NPPCKT  = PKTCOUNT( PAR_IDX )   ! Number of parametric pkts.
        

C.........  If there were no non-blank and non-comment lines, error
        IF( ICNT .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: No non-blank or non-comments lines found ' //
     &             'in uncertainty file.'
            CALL M3MSG2( MESG )
        END IF

C.........  If there was an error reading the file
        IF( UC_ERROR ) EFLAG = .TRUE. 


C.........  If there was any error, exit 
        IF( EFLAG ) THEN
             MESG = 'Problem reading uncertainty file.'
             CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

C.........  Error message for reaching the end of file too soon
999     WRITE( MESG,94010 )
     &         'End of file reached unexpectedly at line', IREC, CRLF()
     &         //BLANK10 //'Check format of uncertainty file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

        END SUBROUTINE SCANUNCERT
