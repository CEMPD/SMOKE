
        SUBROUTINE  CYLST( UNIT , FUNIT, HEADER, MAXROW, NROW, NCOL, 
     &                        FMT, NAMES, CYIDS, RVAL           )

C***********************************************************************
C  subroutine body starts at line  44
C
C  DESCRIPTION:
C     Writes out a report of values to opened file UNIT, from array 
C     RVAL for each county using the format in FMT and the HEADER 
C             
C  PRECONDITIONS REQUIRED:  
C     UNIT must be a valid and opened unit number
C
C  SUBROUTINES AND FUNCTIONS CALLED: 
C   
C     STR2INT, TRIMLEN, LBLANK
C
C  REVISION  HISTORY:
C       prototype 8/98 by J Vukovich
C
C***********************************************************************
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
C****************************************************************************

      IMPLICIT NONE

      INCLUDE 'PARMS3.EXT'
      INCLUDE 'IODECL3.EXT'
      INCLUDE 'FDESC3.EXT'

      INTEGER MAXST, MAXCO

      PARAMETER (MAXST = 60 , MAXCO = 3500 )

C...........   ARGUMENTS and their descriptions:
        
        INTEGER         UNIT           ! length of input vector
        INTEGER         FUNIT          ! state-county names unit
        CHARACTER*256   HEADER         ! header for top of output file
        CHARACTER*256   MESG           ! warning message
        INTEGER         MAXROW         ! max length of output vector
        INTEGER         NROW           ! length of output vector
        INTEGER         NCOL           ! length of input vector
        CHARACTER*(*)   FMT            ! format to use for out values
        CHARACTER*(*)   NAMES ( NCOL ) ! names for columns
        INTEGER         CYIDS ( NROW ) ! length of output vector
        REAL            RVAL  ( MAXROW, NCOL )   !  reported values

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER         STR2INT
        INTEGER         TRIMLEN
        INTEGER         LBLANK
        EXTERNAL        LBLANK, STR2INT, TRIMLEN
        
C...........   LOCAL VARIABLES and their descriptions:
        
        REAL            TOTAL( 1000 )           ! columns totals

        INTEGER         I, M, C, R                   ! counters
        INTEGER         KTMP, KST, IOLD, ICHECK      ! check vars 
        INTEGER         JCHECK                       ! check vars
        INTEGER         L2, LE, LF, LG, LI, LP       ! format choices
        INTEGER         LENLAB                 ! length of column label fmt
        INTEGER         TOT1, TOT2             ! 1st and 2nd part of totals fmt
        INTEGER         WIDTH                 ! of format
        INTEGER         TOT_ST, TOT_CO  ! Total number of states and counties
                                        ! in the FIPS file.

        CHARACTER*1     TOTR                   ! temporary totals format type
        CHARACTER*50    LABELFMT               ! output columns lamel format
        CHARACTER*50    OUTFMT                 ! output table format
        CHARACTER*50    TOTFMT                 ! output totals format
        LOGICAL         INTFLAG                ! true if integer FMT on input

        INTEGER      ISTFIP  ( MAXST )  !  All state FIPS IDs
        CHARACTER*20 STATES  ( MAXST )  !  All state names
        INTEGER      ICOFIP  ( MAXCO )  !  All county FIPS IDs
        CHARACTER*20 COUNTIES( MAXCO )  !  All county names

        INTEGER         FIPNUM  !  Buffer for reading the FIPS code
        CHARACTER*20    FIPNAME !  Buffer for reading the state/county name
        CHARACTER*10    BUF1,BUF2  ! Buffer for writing error messages
        INTEGER         IREAD1           ! Indicates if first function call
C
C.....  save these data so will only have to read once

        DATA          IREAD1 / 1 /
        SAVE          IREAD1, ISTFIP, STATES, ICOFIP, COUNTIES
        SAVE          TOT_ST, TOT_CO

C***********************************************************************
C   begin body of subroutine CYLST

C.........  Perform input-format-dependent changes

        INTFLAG = .FALSE.
        LENLAB  = 0

        L2 = TRIMLEN( FMT )     ! Length of format statement
        LI = INDEX( FMT, 'I' )  ! Location of I (for integer format)
        LF = INDEX( FMT, 'F' )  ! Location of F (for float   format)
        LE = INDEX( FMT, 'E' )  ! Location of E (for expontial format)
        LG = INDEX( FMT, 'G' )  ! Location of G (for general format)
        LP = INDEX( FMT, '.' )  ! Location of . ( for F, E, G format)

        IF( LP .EQ. 0 ) LP = L2 ! If no '.', then use end of format

        IF( LI .GT. 0 ) THEN

            INTFLAG = .TRUE.

            LENLAB = STR2INT( FMT( LI+1:L2 ) )

        ELSEIF( LF .GT. 0 ) THEN

            LENLAB = STR2INT( FMT( LF+1:LP-1 ) )
            TOTR   = FMT( LF:LF )
            TOT1   = LENLAB
            TOT2   = MAX( 0, STR2INT( FMT( LP+1:L2 ) ) - 1 )
           
        ELSEIF( LE .GT. 0 ) THEN

            LENLAB = STR2INT( FMT( LE+1:LP-1 ) )
            TOTR   = FMT( LE:LE )
            TOT1   = LENLAB
            TOT2   = STR2INT( FMT( LP+1:L2 ) )

        ELSEIF( LG .GT. 0 ) THEN

            LENLAB = STR2INT( FMT( LG+1:LP-1 ) )
            TOTR   = FMT( LG:LG )
            TOT1   = LENLAB
            TOT2   = MAX( 0, STR2INT( FMT( LP+1:L2 ) ) - 1 )

        ENDIF
       
        IF( LENLAB .EQ. 0 ) THEN
            CALL M3MSG2( 'Bad format in STATELST, using F5.0' )
            LENLAB = 5
            FMT    = 'F5.0'
            TOT1   = LENLAB
            TOT2   = MAX( 0, STR2INT( FMT( LP+1:L2 ) ) - 1 )
        ENDIF

C.........  Determine width
        WIDTH = 16 + LENLAB * NCOL

C.........  Create format statement for column labels

        WRITE( LABELFMT,94010 ) '( A' , LENLAB,  ', 1X, $ )'

C.........  Create format statement for output data rows

        WRITE( OUTFMT,94020 ) NCOL, FMT( 1:L2 ) 

C.........  Create format statement for output totals

        WRITE( TOTFMT,94010 ) '( I5.5, 2X, ', NCOL, '(' //
     &                        TOTR, TOT1, '.', TOT2, ', 1X ), 2X, A )'

C.........  Read county and state names once from FUNIT
C.........  IREAD1 is used so names will only have to be
C.........  read once.
 
        IF ( IREAD1 .EQ. 1 ) THEN

30       CONTINUE ! Begin "infinite" loop; until EOF.

         READ ( FUNIT, 93070, END = 40 ) FIPNUM, FIPNAME
         IF ( FIPNUM .NE. 0 ) THEN  ! If FIPNUM=0, we just read a blank line
                IF ( MOD ( FIPNUM, 1000 ) .EQ. 0 ) THEN
                    TOT_ST = TOT_ST + 1
                    IF ( TOT_ST .GT. MAXST ) THEN
                        WRITE ( BUF1, 600 ) TOT_ST
                        WRITE ( BUF2, 600 ) MAXST
                        CALL M3EXIT( 'CYLST', 0, 0,
     &                      'Number of states in STCYNAM (' //
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXST=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) //
     &                      ')', 2 )
                    END IF

                    ISTFIP ( TOT_ST ) = FIPNUM / 1000
                    STATES ( TOT_ST ) = FIPNAME
                ELSE
                    TOT_CO = TOT_CO + 1
                    IF ( TOT_CO .GT. MAXCO ) THEN
                        WRITE ( BUF1, 600 ) TOT_CO
                        WRITE ( BUF2, 600 ) MAXCO
                        CALL M3EXIT( 'CYLST', 0, 0,
     &                      'Number of counties in STCYNAM (' //
     &                      BUF1 ( LBLANK(BUF1):TRIMLEN(BUF1) ) //
     &                      ') exceeds program limits (MAXCO=' //
     &                      BUF2 ( LBLANK(BUF2):TRIMLEN(BUF2) ) //
     &                      ')', 2 )
                    END IF

                    ICOFIP ( TOT_CO ) = FIPNUM
                    COUNTIES ( TOT_CO ) = FIPNAME
                END IF
         END IF
         GO TO 30
40       CONTINUE
         IREAD1 = 2

        ENDIF
       
C.........  Write header

        WRITE( UNIT, 93030 ) HEADER( 1:TRIMLEN( HEADER ) )

C.........  Write column labels
        WRITE( UNIT, 93000 ) 'FIPS      '
        WRITE( UNIT, 93040 ) 'CODE      '

        DO 66 C = 1, NCOL
            WRITE( UNIT, LABELFMT ) NAMES ( C )
66      CONTINUE

        WRITE( UNIT, 93000 ) 'COUNTY NAME'

C.........  Initialize totals
        DO 88 C = 1, NCOL
            TOTAL( C ) =  0.
88      CONTINUE 

        IOLD = -99

C.........  Write data

        IF( INTFLAG ) THEN

            DO 100 R = 1, NROW

C.........  jcheck used to see if FIPS code is in state and
C.........  and county name file (FUNIT)
      
                JCHECK = 0

                DO 150 I = 1, TOT_CO
                
                 IF( ICOFIP( I ) .EQ. CYIDS( R ) ) THEN

                   JCHECK = 1

                   ICHECK = CYIDS( R ) / 1000

C..........  If the state id changes then write totals of previous state
C..........  and initialize appropriate variables for the new state  

	           IF ( ICHECK .NE. IOLD ) THEN
                      KST = -99
                      DO M = 1, TOT_ST
                        IF ( IOLD .EQ. ISTFIP( M ) ) KST = M
	              ENDDO

                      IF (IOLD .GT. 0 ) THEN

                       WRITE( UNIT, 93050 )

                       IF ( KST .GT. 0 ) THEN
	                 KTMP = IOLD * 1000	
                         WRITE( UNIT,TOTFMT ) KTMP, ( INT( TOTAL( C ) ),
     &                                    C = 1, NCOL   ), STATES( KST )

                       ELSE

                         WRITE( UNIT,TOTFMT ) 0, ( INT( TOTAL( C ) ),
     &                                     C = 1, NCOL   ),'Total'
		     	
	               ENDIF

                       WRITE( UNIT, 93000 ) ' '

C.........  Write column labels for new state

                       WRITE( UNIT, 93000 ) 'FIPS      '
                       WRITE( UNIT, 93040 ) 'CODE      '

                       DO 166 C = 1, NCOL
                         WRITE( UNIT, LABELFMT ) NAMES ( C )
166                    CONTINUE

                       WRITE( UNIT, 93000 ) 'COUNTY NAME'
                  
                      ENDIF
 
C.........   Initialize variables for new state
 
                      DO  188 C = 1, NCOL
                        TOTAL( C ) =  0.
188                   CONTINUE

                      IOLD = ICHECK

                   ENDIF

C.........   Sum emissions for state
  
                   DO 99 C = 1, NCOL
                      TOTAL( C ) = TOTAL( C ) + RVAL( R,C )
99                 CONTINUE

C.........   Write out emissions for specific county

                   WRITE( UNIT,OUTFMT ) CYIDS( R ),
     &                      ( INT( RVAL( R,C ) ), C = 1, NCOL ),
     &                      COUNTIES( I )

 
                 ENDIF

150             CONTINUE

C..........  FIPS code was NOT found in state and county name file
C..........  write warning

            IF ( JCHECK .EQ. 0 ) THEN

                WRITE( MESG , 95010 ) CYIDS( R )  
                CALL M3WARN( 'CYLST', 0, 0, MESG )

            ENDIF 
              

100         CONTINUE

        ELSE 

            DO 200 R = 1, NROW

C.........  jcheck used to see if FIPS code is in state and
C.........  and county name file (FUNIT)

                JCHECK = 0

                DO 250 I = 1, TOT_CO

                 IF( ICOFIP( I ) .EQ. CYIDS( R ) ) THEN
 
                   JCHECK = 1

                   ICHECK = CYIDS( R ) / 1000

C..........  If the state id changes then write totals of previous state
C..........  and initialize appropriate variables for the new state

                   IF ( ICHECK .NE. IOLD ) THEN
                      KST = -99
                      DO M = 1, TOT_ST
                        IF ( IOLD .EQ. ISTFIP( M ) ) KST = M
                      ENDDO

                      IF (IOLD .GT. 0 ) THEN

                       WRITE( UNIT, 93050 )
                      
                       IF ( KST .GT. 0 ) THEN
                         KTMP = IOLD * 1000
                         WRITE( UNIT,TOTFMT ) KTMP, (  TOTAL( C ) ,
     &                               C = 1, NCOL   ), STATES( KST )

                       ELSE

                         WRITE( UNIT,TOTFMT ) 0, (  TOTAL( C ) ,
     &                                     C = 1, NCOL   ),'Total'
 
                       ENDIF

                       WRITE( UNIT, 93000 ) ' '

C.........  Write column labels for new state

                       WRITE( UNIT, 93000 ) 'FIPS      '
                       WRITE( UNIT, 93040 ) 'CODE      '

                       DO 266 C = 1, NCOL
                         WRITE( UNIT, LABELFMT ) NAMES ( C )
266                    CONTINUE

                       WRITE( UNIT, 93000 ) 'COUNTY NAME'

                      ENDIF

C..........  Initialize variables for new state

                      DO  288 C = 1, NCOL
                        TOTAL( C ) =  0.
288                   CONTINUE

                      IOLD = ICHECK

                   ENDIF

C..........   Sum emissions for state

                   DO 299 C = 1, NCOL
                      TOTAL( C ) = TOTAL( C ) + RVAL( R,C )
299                CONTINUE

C..........   Write out emissions for county

                   WRITE( UNIT,OUTFMT ) CYIDS( R ),
     &                      ((  RVAL( R,C ) ), C = 1, NCOL ),
     &                      COUNTIES( I )


                 ENDIF

250             CONTINUE

C..........  FIPS code was NOT found in state and county name file
C..........  write warning

            IF ( JCHECK .EQ. 0 ) THEN

                WRITE( MESG , 95010 ) CYIDS( R )
                CALL M3WARN( 'CYLST', 0, 0, MESG )

            ENDIF


200         CONTINUE

        ENDIF

C..........  Write out the last state's totals

        WRITE( UNIT, 93050 )
        KTMP = IOLD * 1000

C.......... find state name

        KST = -99
        DO M = 1, TOT_ST
           IF ( IOLD .EQ. ISTFIP( M ) ) KST = M
        ENDDO

        IF( INTFLAG ) THEN
            
           IF ( KST .GT. 0 ) THEN

             WRITE( UNIT,TOTFMT ) KTMP, ( INT( TOTAL( C ) ),
     &                             C = 1, NCOL   ), STATES( KST )

           ELSE

             WRITE( UNIT,TOTFMT ) 0, ( INT( TOTAL( C ) ),
     &                             C = 1, NCOL   ),'Total'

           ENDIF
 
        ELSE 

           IF ( KST .GT. 0 ) THEN

             WRITE( UNIT,TOTFMT ) KTMP, (  TOTAL( C ) ,
     &                             C = 1, NCOL   ), STATES( KST )

           ELSE

             WRITE( UNIT,TOTFMT ) 0, (  TOTAL( C ) ,
     &                             C = 1, NCOL   ),'Total'

           ENDIF

        ENDIF

        WRITE( UNIT,93000 ) ' '
        WRITE( UNIT,93000 ) ' '

        RETURN

C******************  FORMAT  STATEMENTS   ******************************
 
C...........   Formatted file I/O formats............ 93xxx

600     FORMAT ( I10 )
 
93000   FORMAT( A )
 
93030   FORMAT( A, //, 22('-'), ' By County ', 22('-') )
 
93040   FORMAT( A, $ )

C		04/21/99 wts
C93050   FORMAT( 1X, <WIDTH>('-') )
93050   FORMAT( 1X, 80('-') )

93070   FORMAT ( 2X, I5, 1X, A )

C...........   Internal buffering formats............ 94xxx
 
94010   FORMAT( 10( A, :, I3, :, 2X ) )
 
94020   FORMAT( '( I5.5, 2X, ', I3, '(', A, ', 1X ), 2X, A )' )

95000   FORMAT( 2x, A )

95010   FORMAT( 'County name record missing for FIPS = ',I5 ,
     &          ' in STCYNAM file.  This county will',
     &          ' not be written out to the report file!!' )
 
C******************  FORMAT  STATEMENTS   ******************************

        END

