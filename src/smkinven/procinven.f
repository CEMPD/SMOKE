
        SUBROUTINE PROCINVEN( NRAWBP, MXIDAT, FILFMT, PRATIO, INVSTAT )

C***********************************************************************
C  subroutine body starts at line 106
C
C  DESCRIPTION:
C      This subroutine 
C      Many places in the in-line documentation refers to pollutants, but
C      means pollutants or activity data
C
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
        INTEGER         INDEX1
        LOGICAL         ENVYN
        INTEGER         STR2INT

        EXTERNAL        CRLF, INDEX1, ENVYN, STR2INT

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NRAWBP            ! no.raw recs x pol/act
        INTEGER     , INTENT (IN) :: MXIDAT            ! max no. inv pols/actvty
        INTEGER     , INTENT (IN) :: FILFMT            ! input file(s) fmt code
        REAL        , INTENT (IN) :: PRATIO            ! position ratio
        INTEGER     , INTENT(IN OUT) :: INVSTAT( MXIDAT ) ! (<0 actv, >0 pol)

C...........   Variables dimensioned by subroutine arguments
        INTEGER         TMPSTAT( MXIDAT ) ! tmp data status

C...........   Other local variables

        INTEGER         I, J, K, LK, LS, L2, S    ! counter and indices

        INTEGER         IDUP        !  no. dulicate records
        INTEGER         IOS         !  i/o status
        INTEGER         LCAT        !  length of CATEGORY string
        INTEGER         PE, PS      !  pollutant postn end and start in CSOURCA 
        INTEGER         PIPCOD      !  IPOSCOD of previous iteration of loop
        INTEGER         SLEN        !  length of source 

        REAL            EMISI       !  inverse emissions value
        REAL            EMISN       !  new emissions value
        REAL            EMISO       !  old emissions value
        REAL            RIMISS3     !  real typed integer missing value

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

C.........  Initlialze temporary data status
        TMPSTAT = 0  ! array

C.........  Loop through sources X pollutants to determine source IDs and check
C           for duplicates. Also keep a count of the total unique key
C           combinations (CSOURCA without the pollutant position)
C.........  NOTE the last part of the CSOURCA string is the integer position 
C           of the pollutant for that record in the INVPNAM pollutant array 
        LSRCCHR = EMCMISS3
        LK = IMISS3
        S = 0
        SLEN  = SC_ENDP( MXCHRS )
        PS    = SC_BEGP( MXCHRS + 1 )
        PE    = SC_ENDP( MXCHRS + 1 )
        IDUP  = 0
        DO I = 1, NRAWBP
            
            J  = INDEXA( I )

            TSRCCHR = CSOURCA( J )(  1:SLEN ) ! Source characteristics
            TPOLPOS = CSOURCA( J )( PS:PE   ) ! Pos of pollutant (ASCII)

C.............  Update pointer for list of actual pollutants & activities
            K = STR2INT( TPOLPOS )  ! Convert pol/activity position to integer
            
            TMPSTAT( K ) = 2
            IPOSCOD( I ) = K
           
C.............  Increment source count by comparing this iteration to previous
            IF( TSRCCHR .NE. LSRCCHR ) THEN
                S = S + 1
                LSRCCHR = TSRCCHR

C.............  Give message of duplicates are not permitted in inventory
C.............  This IF also implies TSRCCHR = LSRCCHR
            ELSE IF( K .EQ. LK ) THEN

                CALL FMTCSRC( TSRCCHR, NCHARS, BUFFER, L2 )

                IF ( DFLAG ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Duplicate records found for' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                ELSE
                    MESG = 'WARNING: Duplicate records found for' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                END IF

                CALL M3MESG( MESG )

                IDUP = IDUP + 1

            END IF

            LK = K  ! Store pol/activity index for comparison in next iteration

C.............  Assign source ID (to use as an index) for all inv X pol/act
            SRCIDA( I ) = S

        END DO  ! On sources x pollutants/activities

C.........  Set NSRC in for module MODINFO
        NSRC = S

C.........  Report the number of records that were duplicated
        IF( IDUP .NE. 0 ) THEN
            WRITE( MESG,94010 ) 'NOTE: The number of duplicate ' //
     &             'records was', IDUP, '.' 

            IF( .NOT. DFLAG ) THEN
                MESG = MESG( 1:LEN_TRIM( MESG ) ) // CRLF()// BLANK10//
     &                 'The inventory data were summed for ' //
     &                 'these sources.'
            END IF

            CALL M3MSG2( MESG )

        END IF

        IF( EFLAG ) THEN
           MESG = 'Error in raw inventory file(s)'         
           CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

C.........  Use sign of INVSTAT and value of TMPSTAT to set type (pol/act) and
C           indicator of whether it's present or not
        DO I = 1, MXIDAT
            INVSTAT( I ) = INVSTAT( I ) * TMPSTAT( I )
        END DO

C.........  Allocate memory for SMOKE inventory arrays (NOT sources X pollutnts)

        CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .FALSE., NSRC, 
     &               NRAWBP, NPPOL )

C.........  Loop through sources x pollutants/activities to store sorted arrays
C           for output to I/O API file. Use PRATIO to determine the appropriate 
C           position in source-based (non-pollutant, non-activity) arrays.
C.........  Keep case statement outside the loops to speed processing
        LS = IMISS3
        SELECT CASE ( CATEGORY )
        CASE( 'AREA' ) 

             DO I = 1, NRAWBP

                J = INDEXA( I )
                S = SRCIDA( I )
                K = INT( ( REAL( J ) - 0.5 ) * PRATIO ) + 1

                IF( S .NE. LS ) THEN
                    LS  = S
                    IFIP  ( S ) = IFIPA  ( K )
                    TPFLAG( S ) = TPFLGA ( K )
                    INVYR ( S ) = INVYRA ( K )
                    CSCC  ( S ) = CSCCA  ( K )

                    CSOURC( S ) = CSOURCA( J )( 1:SRCLEN3 )
                END IF

            END DO

       CASE( 'MOBILE' )

            DO I = 1, NRAWBP

                J = INDEXA( I )
                S = SRCIDA( I )
                K = INT( ( REAL( J ) - 0.5 ) * PRATIO ) + 1

                IF( S .NE. LS ) THEN
                    LS  = S
                    IFIP  ( S ) = IFIPA  ( K )
                    IRCLAS( S ) = IRCLASA( K )
                    IVTYPE( S ) = IVTYPEA( K ) 
                    TPFLAG( S ) = TPFLGA ( K )
                    INVYR ( S ) = INVYRA ( K )
                    CSCC  ( S ) = CSCCA  ( K )
                    CLINK ( S ) = CLINKA ( K )
                    CVTYPE( S ) = CVTYPEA( K )
                    XLOC1 ( S ) = XLOC1A ( K )
                    YLOC1 ( S ) = YLOC1A ( K )
                    XLOC2 ( S ) = XLOC2A ( K )
                    YLOC2 ( S ) = YLOC2A ( K )
                    SPEED ( S ) = SPEEDA ( K )

                    CSOURC( S ) = CSOURCA( J )( 1:SRCLEN3 )
                END IF

            END DO

        CASE( 'POINT' )

            DO I = 1, NRAWBP

                J = INDEXA( I )
                S = SRCIDA( I )
                K = INT( ( REAL( J ) - 0.5 ) * PRATIO ) + 1

                IF( S .NE. LS ) THEN
                    LS  = S
                    IFIP  ( S )  = IFIPA  ( K )
                    ISIC  ( S )  = ISICA  ( K )
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
                    CORIS ( S )  = CORISA ( K )
                    CBLRID( S )  = CBLRIDA( K )
                    CPDESC( S )  = CPDESCA( K )

                    CSOURC( S )  = CSOURCA( J )( 1:SRCLEN3 )
                END IF

            END DO

        END SELECT

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

C.........  Initialize pollutant/activity-specific values.  Inititalize integer
C           values with the real version of the missing integer flag, since
C           these are stored as reals until output
        POLVAL  = BADVAL3          ! array
        RIMISS3 = REAL( IMISS3 )
        IF( NC1 .GT. 0 ) POLVAL( :,NC1 ) = RIMISS3 ! array
        IF( NC2 .GT. 0 ) POLVAL( :,NC2 ) = RIMISS3 ! array

C.........  Initialize pollutant count per source array
        NPCNT = 0  ! array

C.........  Store pollutant/activity-specific data in sorted order.  For EPS
C           and EMS-95 formats, ensure that any duplicates are aggregated.
C.........  Aggregate duplicate pollutant/activity-specific data (not possible 
C           for IDA format)
C.........  Note that we have already checked to ensure that if there are
C           duplicate data, they are allowed
C.........  Note that pollutants & activities  are stored in output order
C           because they've been previously sorted in part based on their
C           position in the master array of output pollutants/activities
        IF( FILFMT   .EQ. IDAFMT         .OR.
     &    ( FILFMT   .EQ. EMSFMT   .AND.
     &      CATEGORY .EQ. 'MOBILE'      )     ) THEN
            DO I = 1, NRAWBP

                J = INDEXA( I )
                S = SRCIDA( I )

                DO K = 1, NPPOL
                    POLVAL( I, K ) = POLVLA( J, K )
                END DO

                NPCNT( S ) = NPCNT( S ) + 1

            END DO

C.........  For non-IDA formats for area and point sources...
C.........  NOTE- this section assumes that the ozone-day emissions and emission
C           factors are not assigned for this format
        ELSE IF( FILFMT .EQ. EPSFMT .OR. FILFMT .EQ. EMSFMT ) THEN
        
            K = 0
            PIPCOD = IMISS3  ! Previous iteration IPOSCOD 
            LS = IMISS3      ! Previous iteration S
            DO I = 1, NRAWBP

                J = INDEXA( I )
                S = SRCIDA( I )

                IF( S .NE. LS .OR. IPOSCOD( I ) .NE. PIPCOD ) THEN

C.....................  Sum up the number of pollutants/activities by source,
C                       but do this here only, because this part of the IF
C                       statement is for new pollutants
                    NPCNT( S ) = NPCNT( S ) + 1

                    K = K + 1

                    POLVAL( K, NEM ) = POLVLA( J, NEM )
                    POLVAL( K, NCE ) = POLVLA( J, NCE )
                    POLVAL( K, NRE ) = POLVLA( J, NRE )
                    IF( NRP .GT. 0 ) POLVAL( K, NRP ) = POLVLA( J, NRP )

C.....................  Update IPOSCOD (the position of the pol/act in the 
C                       master list) to align properly with summed emissions.
C.....................  NOTE- K is always <= to I, so there is no conflict for
C                       reassigning IPOSCOD.
                    IPOSCOD( K ) = IPOSCOD( I )
                    PIPCOD       = IPOSCOD( I ) 

C.................  If the existing value is defined, sum with new emissions
C                   or activity and use weighted average for control factors
                ELSE

                    EMISN = POLVLA( J, NEM )
                    EMISO = POLVAL( K, NEM )
                    POLVAL( K, NEM ) = EMISO + EMISN

                    EMISI = 1. / POLVAL( K, NEM )   ! Compute inverse only once
                    POLVAL( K,NCE ) = ( POLVAL( K,NCE )*EMISO + 
     &                                  POLVLA( J,NCE )*EMISN  ) * EMISI
                    POLVAL( K,NRE ) = ( POLVAL( K,NRE )*EMISO + 
     &                                  POLVLA( J,NRE )*EMISN  ) * EMISI
                    IF( NRP .GT. 0 ) 
     &              POLVAL( K,NRP ) = ( POLVAL( K,NRP )*EMISO + 
     &                                  POLVLA( J,NRP )*EMISN  ) * EMISI

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
