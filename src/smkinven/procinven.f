
        SUBROUTINE PROCINVEN( NRAWBP, MXIPOL, FILFMT, PRATIO, INVSTAT )

C***********************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine 

C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 4/99 by M. Houyoux
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

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC 

C.........  This module contains the information about the source category
        USE MODINFO

        IMPLICIT NONE

C...........   INCLUDES

        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions
        CHARACTER*2     CRLF
        LOGICAL         ENVYN
        INTEGER         STR2INT

        EXTERNAL        CRLF, ENVYN, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NRAWBP            ! no.raw records x pol
        INTEGER     , INTENT (IN) :: MXIPOL            ! max no. inv pols
        INTEGER     , INTENT (IN) :: FILFMT            ! input file(s) fmt code
        REAL        , INTENT (IN) :: PRATIO            ! position ratio
        INTEGER     , INTENT(OUT) :: INVSTAT( MXIPOL ) ! (0=not in inventory)

C...........   Other local variables

        INTEGER         I, J, K, LK, LS, L2, S    ! counter and indices

        INTEGER         IOS         !  i/o status
        INTEGER         LCAT   ! length of CATEGORY string
        INTEGER         PIPCOD      !  IPOSCOD of previous iteration of loop
 
        REAL            EMISI       !  inverse emissions value
        REAL            EMISN       !  new emissions value
        REAL            EMISO       !  old emissions value

        LOGICAL         DFLAG             ! true: if should error on duplicates
        LOGICAL      :: EFLAG  = .FALSE.  ! true: error occured

        CHARACTER*5     TPOLPOS     !  Temporary pollutant position
        CHARACTER*300   BUFFER      !  input file line buffer
        CHARACTER*300   MESG        ! message buffer 

        CHARACTER(LEN=ALLLEN3)  LSRCCHR     !  previous CSOURC
        CHARACTER(LEN=ALLLEN3)  TSRCCHR     !  tmporary CSOURC

        CHARACTER*16 :: PROGNAME = 'PROCINVEN' ! program name

C***********************************************************************
C   begin body of subroutine PROCINVEN

C.........  Get settings from the environment
        DFLAG = ENVYN( 'RAW_DUP_CHECK',
     &                 'Check for duplicate species-records',
     &                 .FALSE., IOS )

C.........  Loop through sources X pollutants to determine source IDs and check
C           for duplicates. Also keep a count of the total unique key
C           combinations (CSOURCA without the pollutant position)
C.........  NOTE: The last part of the CSOURCA string is the integer position 
C           of the pollutant for that record in the INVPNAM pollutant array 
        LSRCCHR = EMCMISS3
        LK = IMISS3
        S = 0
        DO I = 1, NRAWBP
            
            J  = INDEXA( I )

            TSRCCHR = CSOURCA( J )(       1:SRCLEN3 ) ! Source characteristics
            TPOLPOS = CSOURCA( J )( POLPOS3:ALLLEN3 ) ! Pos of pollutant (ASCII)

C.............  Update pointer for list of actual pollutants
            K = STR2INT( TPOLPOS )  ! Convert pollutant code to integer
            INVSTAT( K ) = 1
            IPOSCOD( I ) = K
           
C.............  Increment source count by comparing this iteration to previous
            IF( TSRCCHR .NE. LSRCCHR ) THEN
                S = S + 1
                LSRCCHR = TSRCCHR

C.............  Give message of duplicates are not permitted in inventory
C.............  This IF also implies TSRCCHR = LSRCCHR
            ELSE IF( K .EQ. LK ) THEN

                CALL FMTCSRC( TSRCCHR, 8, BUFFER, L2 )

                IF ( DFLAG ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Duplicate emissions found for' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                ELSE
                    MESG = 'WARNING: Duplicate emissions found for' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                END IF

                CALL M3MESG( MESG )

            END IF

            LK = K  ! Store pollutant index for comparison in next iteration

C.............  Assign source ID (to use as an index) for all inv X pol
            SRCIDA( I ) = S

        END DO  ! On sources x pollutants

C.........  Set NSRC in for module MODINFO
        NSRC = S

        IF( EFLAG ) THEN
           MESG = 'Error in raw inventory file(s)'         
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Allocate memory for SMOKE inventory arrays (NOT sources X pollutnts)

        CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .FALSE., NSRC, 
     &               NRAWBP, NPPOL )

C.........  Loop through sources x pollutants to store sorted arrays for output
C           to I/O API file. Use PRATIO to determine the appropriate position 
C           in source-based (non-pollutant) arrays.
        LS = IMISS3
        DO I = 1, NRAWBP

            J = INDEXA( I )
            S = SRCIDA( I )
            K = INT( ( REAL( J ) - 0.5 ) * PRATIO ) + 1

            IF( S .NE. LS ) THEN
                LS  = S
                IFIP  ( S )  = IFIPA  ( K )
                ISIC  ( S )  = ISICA  ( K )
                IORIS ( S )  = IORISA ( K )
                IDIU  ( S )  = IDIUA  ( K )
                IWEK  ( S )  = IWEKA  ( K )
                TPFLAG( S )  = TPFLGA ( K )
                INVYR ( S )  = INVYRA ( K )
                XLOCA ( S )  = XLOCAA ( K )
                YLOCA ( S )  = YLOCAA ( K )
                STKHT ( S )  = STKHTA ( K )
                STKDM ( S )  = STKDMA ( K )
                STKTK ( S )  = STKTKA ( K )
                STKVE ( S )  = STKVEA ( K )
                CSCC  ( S )  = CSCCA  ( K )
                CBLRID( S )  = CBLRIDA( K )
                CPDESC( S )  = CPDESCA( K )

                CSOURC( S )  = CSOURCA( J )( 1:SRCLEN3 )
            ENDIF

        ENDDO

C.........  Deallocate local memory for per-source unsorted arrays
        CALL SRCMEM( CATEGORY, 'UNSORTED', .FALSE., .FALSE., 
     &               1, 1, 1 )

C.........  Allocate memory for aggregating any duplicate pol-specific data
C.........  Note that the POLVAL array that contains the pollutant-specific
C           data is dimensioned to output only one pollutant at a time. This is 
C           because we may need to be able to handle many pollutants in the 
C           future, and the memory requirements would be prohibitive if all of
C           the memory were allocated at the same time.           

        CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .TRUE., NSRC, 
     &               NRAWBP, NPPOL )

        POLVAL = AMISS3

C.........  Store pollutant-specific data in sorted order.  For EPS and EMS-95
C           formats, ensure that any duplicates are aggregated.
C.........  Aggregate duplicate pollutant-specific data (not possible 
C           for IDA format)
C.........  NOTE: We have already checked to ensure that if there are duplicate
C           emissions, they are allowed
C.........  NOTE: Pollutants are stored in output order because they've been
C           previously sorted in part based on their position in the master
C           array of output pollutants
        IF( FILFMT .EQ. IDAFMT ) THEN
            DO I = 1, NRAWBP

                J = INDEXA( I )

                DO K = 1, NPPOL
                    POLVAL( I, K ) = POLVLA( J, K )
                END DO
            END DO

C.............  Set the number of pollutants per source to the constant value
            NPCNT = NRAWBP / NSRC   ! NPCNT is array

        ELSE IF( FILFMT .EQ. EPSFMT .OR. FILFMT .EQ. EMSFMT ) THEN
        
            NPCNT = 0  ! initialize pollutant count array

            K = 0
            PIPCOD = IMISS3  ! Previous iteration IPOSCOD 
            LS = IMISS3      ! Previous iteration S
            DO I = 1, NRAWBP

                J = INDEXA( I )
                S = SRCIDA( I )

                IF( S .NE. LS .OR. IPOSCOD( I ) .NE. PIPCOD ) THEN

C.....................  Sum up the number of pollutants by source, but do this
C                       here only, because this part of the IF statement is for
C                       new pollutants
                    NPCNT( S ) = NPCNT( S ) + 1

                    K = K + 1

                    POLVAL( K, NEM ) = POLVLA( J, NEM )
                    POLVAL( K, NCE ) = POLVLA( J, NCE )
                    POLVAL( K, NRE ) = POLVLA( J, NRE )

                    PIPCOD = IPOSCOD( I ) 

C.................  If the existing value is defined, sum with new emissions
C                   and use weighted average for control factors
C.................  No need to change NPCNT because it is already 1 for all
                ELSE

                    EMISN = POLVLA( J, NEM )
                    EMISO = POLVAL( K, NEM )
                    POLVAL( K, NEM ) = EMISO + EMISN

                    EMISI = 1. / POLVAL( K, NEM )   ! Compute inverse only once
                    POLVAL( K,NCE ) = ( POLVAL( K,NCE )*EMISO + 
     &                                  POLVLA( J,NCE )*EMISN  ) * EMISI
                    POLVAL( K,NRE ) = ( POLVAL( K,NRE )*EMISO + 
     &                                  POLVLA( J,NRE )*EMISN  ) * EMISI

                END IF

                LS = S

            END DO

        END IF

C.........  Deallocate memory for unsorted pollutant arrays
        CALL SRCMEM( CATEGORY, 'UNSORTED', .FALSE., .TRUE., 1, 1, 1 )
                
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE PROCINVEN
