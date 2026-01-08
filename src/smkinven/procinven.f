
        SUBROUTINE PROCINVEN( NRAWBP, UDEV, YDEV, CDEV, LDEV )

C**************************************************************************
C  subroutine body starts at line 114
C
C  DESCRIPTION:
C      This subroutine sorts and stores the inventory data including source
C      specific characteristics and emissions data.
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
C**************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C 
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C 
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$ 
C
C       Updated with USE M3UTILIO by Huy Tran UNC-IE on 2026-01
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE M3UTILIO

        USE MODSOURC, ONLY: INRECA, POLVLA, INDEXA, IPOSCODA, ICASCODA, 
     &                      CSOURC, NPCNT, POLVAL, IPOSCOD         

C.........  This module contains the lists of unique inventory information
        USE MODLISTS, ONLY: MXIDAT, INVSTAT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NEM, NDY, NEF, NCE, NRE, NRP, 
     &                     NPPOL, NSRC, NC1, NC2, NCHARS 

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
C        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

C...........   EXTERNAL FUNCTIONS and their descriptions
C       CHARACTER(2)    CRLF
C       INTEGER         ENVINT
C       LOGICAL         ENVYN

C       EXTERNAL        CRLF, ENVINT, ENVYN

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (IN) :: NRAWBP   ! no. raw recs x pol/act
        INTEGER , INTENT (IN) :: UDEV     ! unit no. for non-HAP exclusions
        INTEGER , INTENT (IN) :: YDEV     ! unit no. for area-to-point
        INTEGER , INTENT (IN) :: CDEV     ! SCC descriptions unit no.
        INTEGER , INTENT (IN) :: LDEV     ! log file unit no.

C...........   Variables dimensioned by subroutine arguments
        INTEGER         TMPSTAT( MXIDAT ) ! tmp data status
c NOTE: GET RID OF THIS ALLOCATION METHOD!

C...........   Other local variables
        INTEGER         I, J, K, LS, L2, S    ! counter and indices

        INTEGER         CASNUM      !  current CAS number from ICASCODA
        INTEGER      :: IDUP = 0    !  no. dulicate records
        INTEGER         IOS         !  i/o status
        INTEGER         MXERR       !  max no. errors
        INTEGER         MXWARN      !  max no. warnings
        INTEGER      :: NERR  = 0   !  no. errors
        INTEGER      :: NWARN = 0   !  no. warnings
        INTEGER         PCAS        !  previous CAS number
        INTEGER         PIPCOD      !  previous pollutant code
        INTEGER         POLCOD      !  current pollutant code from IPOSCODA

        REAL            EMISI       !  inverse emissions value
        REAL            EMISN       !  new emissions value
        REAL            EMISN_DY    !  new average day emissions value
        REAL            EMISO       !  old emissions value
        REAL            EMISO_DY    !  old average day emissions value
        REAL            EMIST       !  total old and new emissions
        REAL            RIMISS3     !  real typed integer missing value

        LOGICAL         ACTFLAG           ! true: current pollutant is activity
        LOGICAL         DFLAG             ! true: if should error on duplicates
        LOGICAL      :: NEGOK = .FALSE.   ! true: negative values are okay
        LOGICAL      :: EFLAG  = .FALSE.  ! true: error occured
        LOGICAL      :: CE_100_FLAG  = .FALSE. ! true: control eff of 100 found
        LOGICAL      :: RE_ZERO_FLAG = .FALSE. ! true: rule effective of 0 found
        LOGICAL      :: RP_ZERO_FLAG = .FALSE. ! true: rule penetration of 0 found

        CHARACTER(256)     BUFFER      !  input file line buffer
        CHARACTER(256)     MESG        !  message buffer 

        CHARACTER(16) :: PROGNAME = 'PROCINVEN' ! program name

C***********************************************************************
C   begin body of subroutine PROCINVEN

C.........  Get settings from the environment
        DFLAG = ENVYN( 'RAW_DUP_CHECK',
     &                 'Check for duplicate species-records',
     &                 .FALSE., IOS )
        NEGOK = ENVYN( 'ALLOW_NEGATIVE',
     &                 'Allow negative output data',
     &                 .FALSE., IOS )

        MXERR  = ENVINT( ERRSET  , ' ', 100, I )
        MXWARN = ENVINT( WARNSET , ' ', 100, I )

C.........  Allocate memory for sorted inventory data
        ALLOCATE( NPCNT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NPCNT', PROGNAME )
        ALLOCATE( POLVAL( NRAWBP,NPPOL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'POLVAL', PROGNAME )
        ALLOCATE( IPOSCOD( NRAWBP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IPOSCOD', PROGNAME )

C.........  Initialize pollutant/activity-specific values.  
C.........  Initialize annual and average day values with 0.
C.........  Inititalize integer
C           values with the real version of the missing integer flag, since
C           these are stored as reals until output
        IF( CATEGORY /= 'MOBILE' ) THEN
            POLVAL( :,1:2 )    = 0.               ! array
            POLVAL( :,3:NPPOL) = BADVAL3          ! array
        END IF

        RIMISS3 = REAL( IMISS3 )        
        IF( NC1 .GT. 0 ) POLVAL( :,NC1 ) = RIMISS3 ! array
        IF( NC2 .GT. 0 ) POLVAL( :,NC2 ) = RIMISS3 ! array

C.........  Initialize pollutant count per source array
        NPCNT = 0  ! array

C.........  Initialize temporary data status
        TMPSTAT = 0  ! array

C.........  Store pollutant/activity-specific data in sorted order. Ensure that
C           any duplicates are aggregated.
C.........  Give warnings or errors when duplicates are encountered
C.........  Note that pollutants & activities  are stored in output order
C           because they've been previously sorted in part based on their
C           position in the master array of output pollutants/activities
        K = 0
        PIPCOD = IMISS3   ! Previous iteration pollutant code
        PCAS   = IMISS3   ! Previous CAS number
        LS     = IMISS3   ! Previous iteration S
        DO I = 1, NRAWBP

            J = INDEXA( I )
            S = INRECA( J )
            
            POLCOD = IPOSCODA( J )
            CASNUM = ICASCODA( J )

C.............  Update pointer for list of actual pollutants and activities
            TMPSTAT( POLCOD ) = 2

C.............  If current source, pollutant, and CAS number match previous,
C               print duplicate error or warning message            
            IF( S      == LS     .AND. 
     &          POLCOD == PIPCOD .AND. 
     &          CASNUM == PCAS         ) THEN
                
                CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                
                IF( DFLAG .AND. NERR <= MXERR ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Duplicate records found for' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )
                    NERR = NERR + 1
                
                ELSE IF( NWARN <= MXWARN )THEN
                    MESG = 'WARNING: Duplicate records found for' //
     &                     CRLF() // BLANK5 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                END IF
                
                IDUP = IDUP + 1
            
            END IF

C.............  Skip rest of loop if duplicates are not allowed
            IF( EFLAG ) THEN
                LS     = S
                PIPCOD = POLCOD
                PCAS   = CASNUM
                CYCLE
            END IF

C.............  Check if current pollutant is an activity
            ACTFLAG = .FALSE.
            IF( INVSTAT( POLCOD ) < 0 ) ACTFLAG = .TRUE.
            
C.............  Reset emissions values to zero, if it's negative
            IF ( POLVLA( J, NEM ) < 0 .AND.
     &           POLVLA( J, NEM ) > AMISS3 ) THEN

    	        IF( .NOT. NEGOK ) POLVLA( J, NEM ) = 0.

                IF ( NWARN < MXWARN .AND. POLCOD /= PIPCOD ) THEN
                    CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )

C.....................  If negative are not okay
                    IF( NEGOK ) THEN
                        IF( ACTFLAG ) THEN
                            MESG = 'WARNING: Negative inventory data'//
     &                             'retained for:' // CRLF() //
     &                             BLANK5 // BUFFER( 1:L2 )
                        ELSE
                            MESG = 'WARNING: Negative annual data ' //
     &                             'retained for:' //
     &                              CRLF() // BLANK5 // BUFFER( 1:L2 )
                        END IF

C.....................  If negative values are not okay
                    ELSE
                        IF( ACTFLAG ) THEN
                            MESG = 'WARNING: Negative inventory data'//
     &                             'reset to zero for:' // CRLF() //
     &                             BLANK5 // BUFFER( 1:L2 )
                        ELSE
                            MESG = 'WARNING: Negative annual data ' //
     &                             'reset to zero for:' //
     &                              CRLF() // BLANK5 // BUFFER( 1:L2 )
                        END IF
                    END IF

                    CALL M3MESG( MESG )
                    NWARN = NWARN + 1
                END IF
            END IF

            IF( .NOT. ACTFLAG ) THEN                
                IF ( POLVLA( J, NDY ) < 0 .AND.
     &               POLVLA( J, NDY ) > AMISS3 ) THEN

                    IF( .NOT. NEGOK ) POLVLA( J, NDY ) = 0.
    
                    IF ( NWARN < MXWARN .AND. POLCOD /= PIPCOD ) THEN
                        CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )

C.........................  If negative emissions are okay
                        IF( NEGOK ) THEN
                            MESG = 'WARNING: Negative average day ' //
     &                             'data retained for:' //
     &                             CRLF() // BLANK5 // BUFFER( 1:L2 )

C.........................  If negative emissions are not okay
                        ELSE
                            MESG = 'WARNING: Negative average day ' //
     &                             'data reset to zero for:' //
     &                             CRLF() // BLANK5 // BUFFER( 1:L2 )
                        END IF  

                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF
                END IF
            END IF

C.............  Convert 0-100 based values to 0-1 based values.
C.............  Check control efficiency, rule effectiveness, and rule 
C               penetration and if missing, set to default value.
C.............  CE default = 0., RP default = 1., RE default = 1.
C.............  Control efficiency
            IF ( NCE > 0 ) THEN                
                IF( POLVLA( J, NCE ) < 0. ) THEN
                    POLVLA( J, NCE ) = 0.
                ELSE IF( POLVLA( J, NCE ) == 100. ) THEN
                    POLVLA( J, NCE ) = 0.
                    CE_100_FLAG = .TRUE.
                ELSE
                    POLVLA( J, NCE ) = POLVLA( J, NCE ) / 100.
                END IF
            END IF

C.............  Rule effectiveness
            IF ( NRE > 0 ) THEN
                IF( POLVLA( J, NRE ) < 0. ) THEN
                    POLVLA( J, NRE ) = 1.
                ELSE IF( POLVLA( J, NRE ) == 0. ) THEN
                    POLVLA( J, NRE ) = 1.
                    RE_ZERO_FLAG = .TRUE.
                ELSE
                    POLVLA( J, NRE ) = POLVLA( J, NRE ) / 100.
                END IF
            END IF

C.............  Rule penetration
            IF ( NRP > 0 ) THEN
                IF( POLVLA( J, NRP ) < 0. ) THEN
                    POLVLA( J, NRP ) = 1.
                ELSE IF( POLVLA( J, NRP ) == 0. ) THEN
                    POLVLA( J, NRP ) = 1.
                    RP_ZERO_FLAG = .TRUE.
                ELSE
                    POLVLA( J, NRP ) = POLVLA( J, NRP ) / 100.
                END IF
            END IF

C.............  For a new source or a new pollutant code...
            IF( S /= LS .OR. POLCOD /= PIPCOD ) THEN

C.................  Sum up the number of pollutants/activities by source,
C                   but do this here only, because this part of the IF
C                   statement is for new pollutants
                NPCNT( S ) = NPCNT( S ) + 1
                K = K + 1

                POLVAL( K, NEM ) = POLVLA( J, NEM )
                
                IF( .NOT. ACTFLAG ) THEN
                    POLVAL( K, NDY ) = POLVLA( J, NDY )
                ELSE
                    POLVAL( K, NDY ) = BADVAL3
                END IF
                    
                IF( NCE > 0 ) POLVAL( K, NCE ) = POLVLA( J, NCE )
                IF( NRE > 0 ) POLVAL( K, NRE ) = POLVLA( J, NRE )
                IF( NEF > 0 ) POLVAL( K, NEF ) = POLVLA( J, NEF )
                IF( NRP > 0 ) POLVAL( K, NRP ) = POLVLA( J, NRP )

C.................  Store position of the pol/act in the 
C                   master list in sorted array.
                IPOSCOD( K ) = POLCOD

C.............  If the existing value is defined, sum with new emissions
C               or activity and use weighted average for control factors
            ELSE

                EMISN    = 0.
                EMISO    = 0.
                EMISN_DY = 0.
                EMISO_DY = 0.

                IF( POLVAL( K, NEM ) >= 0. ) THEN
                    EMISN = POLVLA( J, NEM )
                    EMISO = POLVAL( K, NEM )
                    POLVAL( K, NEM ) = EMISO + EMISN
                END IF

                IF( .NOT. ACTFLAG ) THEN
                    IF( POLVAL( K, NDY ) >= 0. ) THEN
                        EMISN_DY = POLVLA( J, NDY )
                        EMISO_DY = POLVAL( K, NDY )
                        POLVAL( K, NDY ) = EMISO_DY + EMISN_DY
    
C.........................  Use average day emissions for weighting if 
C                           annual emissions are not available.
                        IF( EMISN == 0. ) EMISN = EMISN_DY
                        IF( EMISO == 0. ) EMISO = EMISO_DY
                    END IF
                END IF

C.................  Compute inverse only once
                EMIST = EMISN + EMISO
                IF( EMIST > 0. ) THEN
                    EMISI = 1. / EMIST

C.................  Continue in loop if zero emissions 
                ELSE
                    CYCLE

                END IF

C.................  Weight the control efficiency, rule effectiveness, and 
C                   rule penetration based on the emission values
                IF ( NCE > 0 ) 
     &          POLVAL( K,NCE ) = ( POLVAL( K,NCE )*EMISO + 
     &                              POLVLA( J,NCE )*EMISN  ) * EMISI
                IF ( NRE > 0 ) 
     &          POLVAL( K,NRE ) = ( POLVAL( K,NRE )*EMISO + 
     &                              POLVLA( J,NRE )*EMISN  ) * EMISI
                IF( NRP > 0 ) 
     &          POLVAL( K,NRP ) = ( POLVAL( K,NRP )*EMISO + 
     &                              POLVLA( J,NRP )*EMISN  ) * EMISI
                IF( NEF > 0 ) THEN
                    IF ( POLVAL( K,NEF ) > 0 ) 
     &              POLVAL( K,NEF ) = ( POLVAL( K,NEF )*EMISO + 
     &                                  POLVLA( J,NEF )*EMISN  ) * EMISI
                END IF

            END IF

            LS     = S
            PIPCOD = POLCOD
            PCAS   = CASNUM

        END DO
        
C.........  Report the number of records that were duplicated
        IF( IDUP /= 0 ) THEN
            WRITE( MESG,94010 ) 'NOTE: The number of duplicate ' //
     &             'records was', IDUP, '.'
                   
            IF( .NOT. DFLAG ) THEN
                MESG = TRIM( MESG ) // CRLF() // BLANK10 //
     &                 'The inventory data were summed for ' //
     &                 'these sources.'
            END IF
            
            CALL M3MSG2( MESG )
        END IF

        IF( EFLAG ) THEN
            MESG = 'Duplicates found in raw inventory file(s)'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        
C.........  Print warnings about changed control efficiency, rule
C           effectiveness, and rule penetration
        IF( CE_100_FLAG ) THEN
            MESG = 'WARNING: Some base year control efficiency values'//
     &             ' that were input as 100%' // CRLF() // BLANK10 //
     &             'were reset to 0%. A control efficiency of 100% ' //
     &             'would mean a control' // CRLF() // BLANK10 // 
     &             'is 100% effective and emissions should be zero.'
            CALL M3MSG2( MESG )
        END IF

        IF( RE_ZERO_FLAG ) THEN
            MESG = 'WARNING: Some base year rule effectiveness values'//
     &             ' that were input as 0%' // CRLF() // BLANK10 //
     &             'were reset to 100%. A rule effectiveness of 0% ' //
     &             'would mean a control' // CRLF() //BLANK10 // 
     &             'works 0% of the time, and this does not make sense.'
            CALL M3MSG2( MESG )
        END IF

        IF( RP_ZERO_FLAG ) THEN
            MESG = 'WARNING: Some base year rule penetration values'//
     &             ' that were input as 0%' // CRLF() // BLANK10 //
     &             'were reset to 100%. A rule penetration of 0% ' //
     &             'would mean a control' // CRLF() //BLANK10 // 
     &             'applies to none of the sources, and this does ' //
     &             'not make sense.'
            CALL M3MSG2( MESG )
        END IF

C.........  Deallocate memory for unsorted pollutant arrays
        DEALLOCATE( POLVLA, IPOSCODA, INDEXA, 
     &              INRECA, ICASCODA )

C.........  Use sign of INVSTAT and value of TMPSTAT to set type (pol/act) and
C           indicator of whether it's present or not
        DO I = 1, MXIDAT
            INVSTAT( I ) = INVSTAT( I ) * TMPSTAT( I )
        END DO

C.........  Call adjustment routine to create area-to-point sources and
C           read in nonhap exclusion file
        CALL ADJUSTINV( K, UDEV, YDEV, CDEV, LDEV ) 
                
        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )
 
        END SUBROUTINE PROCINVEN
