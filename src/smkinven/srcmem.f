
        SUBROUTINE SRCMEM( CATEGORY, SORTTYPE, AFLAG, PFLAG, NDIM1, 
     &                     NDIM2, NDIM3 )

C***************************************************************************
C  subroutine body starts at line 
C
C  DESCRIPTION:
C      This subroutine allocates memory for the inventory sorted and 
C      unsorted arrays. IT deallocates memory, when necessary.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 10/98 by M. Houyoux
C
C***************************************************************************
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
C***************************************************************************

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: INDEXA, TPFLGA, INVYRA, CSCCA,
     &                      IPOSCOD, CSOURCA, POLVLA, IRCLASA, IVTYPEA,
     &                      CVTYPEA, CLINKA, XLOC1A, YLOC1A, XLOC2A,
     &                      YLOC2A, IDIUA, IWEKA, XLOCAA, YLOCAA,
     &                      STKHTA, STKDMA, STKTKA, STKVEA, CORISA,
     &                      CBLRIDA, CPDESCA, CIFIP, TPFLAG, INVYR,
     &                      TZONES, NPCNT, CSCC, CSOURC, POLVAL, XLOCA,
     &                      YLOCA, CELLID, IRCLAS, IVTYPE, CVTYPE,
     &                      XLOC1, YLOC1, XLOC2, YLOC2, CISIC, IDIU,
     &                      STKHT, STKDM, STKTK, STKVE, CORIS, CBLRID,
     &                      CPDESC, SRCIDA, CLINK, IWEK

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CATEGORY  ! source category
        CHARACTER(*), INTENT (IN) :: SORTTYPE  ! sorted or unsorted
        LOGICAL     , INTENT (IN) :: AFLAG     ! true: allocate
        LOGICAL     , INTENT (IN) :: PFLAG     ! true: act on pollutant-spec
        INTEGER     , INTENT (IN) :: NDIM1     ! dim for non-pol-spec arrays
        INTEGER     , INTENT (IN) :: NDIM2     ! dim 1 for pol-spec arrays
        INTEGER     , INTENT (IN) :: NDIM3     ! dim 2 for pol-spec arrays

C...........   Other local variables
        INTEGER         IOS    ! memory allocation status

        LOGICAL         UFLAG  ! true: allocate non-pollutant variables

        CHARACTER(300)  MESG

        CHARACTER(16) :: PROGNAME =  'SRCMEM' ! program name

C***********************************************************************
C   begin body of subroutine SRCMEM

        UFLAG = ( AFLAG .AND. .NOT. PFLAG )

        SELECT CASE ( SORTTYPE )

C......... Unsorted...
        CASE( 'UNSORTED' )

C.............  Allocate variables irrespective of PFLAG
            IF( UFLAG .AND. .NOT. ASSOCIATED( INDEXA ) ) THEN
                ALLOCATE( INDEXA( NDIM2 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )
            END IF

C.............  Allocate for any source category
            IF( UFLAG .AND. .NOT. ASSOCIATED( TPFLGA ) ) THEN
                ALLOCATE( TPFLGA( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'TPFLGA', PROGNAME )
            END IF

            IF( UFLAG .AND. .NOT. ASSOCIATED( INVYRA ) ) THEN
                ALLOCATE( INVYRA( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVYRA', PROGNAME )
            END IF
 
            IF( UFLAG .AND. .NOT. ASSOCIATED( CSCCA ) ) THEN
                ALLOCATE( CSCCA( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSCCA', PROGNAME )
            END IF

            IF( UFLAG .AND. .NOT. ASSOCIATED( SRCIDA ) ) THEN
                ALLOCATE( SRCIDA( NDIM2 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'SRCIDA', PROGNAME )
            END IF
 
            IF( UFLAG .AND. .NOT. ASSOCIATED( IPOSCOD ) ) THEN
                ALLOCATE( IPOSCOD( NDIM2 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IPOSCOD', PROGNAME )
            END IF

            IF( UFLAG .AND. .NOT. ASSOCIATED( CSOURCA ) ) THEN
                ALLOCATE( CSOURCA( NDIM2 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSOURCA', PROGNAME )
            END IF

            IF( PFLAG .AND. AFLAG 
     &                .AND. .NOT. ASSOCIATED( POLVLA ) ) THEN
                ALLOCATE( POLVLA( NDIM2, NDIM3 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'POLVLA', PROGNAME )
            END IF

C.............  Deallocate for any source category
            IF( .NOT. AFLAG .AND. .NOT. PFLAG ) THEN
                IF( ASSOCIATED( TPFLGA )  ) DEALLOCATE( TPFLGA )
                IF( ASSOCIATED( INVYRA )  ) DEALLOCATE( INVYRA )
                IF( ASSOCIATED( CSCCA )   ) DEALLOCATE( CSCCA )
                IF( ASSOCIATED( CSOURCA ) ) DEALLOCATE( CSOURCA )

            ELSE IF( .NOT. AFLAG ) THEN
                 IF( ASSOCIATED( INDEXA  ) ) DEALLOCATE( INDEXA )
                 IF( ASSOCIATED( SRCIDA  ) ) DEALLOCATE( SRCIDA )
                 IF( ASSOCIATED( POLVLA  ) ) DEALLOCATE( POLVLA )

            END IF               

C.............  Allocate specifically based on source category
            SELECT CASE( CATEGORY )
            CASE( 'AREA' )
            CASE( 'MOBILE' )
 
                IF( UFLAG .AND. .NOT. ALLOCATED( IRCLASA ) ) THEN
                    ALLOCATE( IRCLASA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IRCLASA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( IVTYPEA ) ) THEN
                    ALLOCATE( IVTYPEA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IVTYPEA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( CVTYPEA ) ) THEN
                    ALLOCATE( CVTYPEA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CVTYPEA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( CLINKA ) ) THEN
                    ALLOCATE( CLINKA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CLINKA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( XLOC1A ) ) THEN
                    ALLOCATE( XLOC1A( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'XLOC1A', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( YLOC1A ) ) THEN
                    ALLOCATE( YLOC1A( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'YLOC1A', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( XLOC2A ) ) THEN
                    ALLOCATE( XLOC2A( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'XLOC2A', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( YLOC2A ) ) THEN
                    ALLOCATE( YLOC2A( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'YLOC2A', PROGNAME )
                END IF
 
                IF( .NOT. PFLAG .AND. 
     &              .NOT. AFLAG .AND. ALLOCATED( CVTYPEA ) ) 
     &              DEALLOCATE( IRCLASA, IVTYPEA, CVTYPEA, CLINKA, 
     &                          XLOC1A, YLOC1A, XLOC2A, YLOC2A )

            CASE( 'POINT' )
 
                IF( UFLAG .AND. .NOT. ALLOCATED( IDIUA ) ) THEN
                    ALLOCATE( IDIUA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IDIUA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( IWEKA ) ) THEN
                    ALLOCATE( IWEKA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IWEKA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ASSOCIATED( XLOCAA ) ) THEN
                    ALLOCATE( XLOCAA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'XLOCAA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ASSOCIATED( YLOCAA ) ) THEN
                    ALLOCATE( YLOCAA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'YLOCAA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( STKHTA ) ) THEN
                    ALLOCATE( STKHTA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'STKHTA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( STKDMA ) ) THEN
                    ALLOCATE( STKDMA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'STKDMA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( STKTKA) ) THEN
                    ALLOCATE( STKTKA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'STKTKA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( STKVEA ) ) THEN
                    ALLOCATE( STKVEA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'STKVEA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( CORISA ) ) THEN
                    ALLOCATE( CORISA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CORISA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( CBLRIDA ) ) THEN
                    ALLOCATE( CBLRIDA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CBLRIDA', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( CPDESCA ) ) THEN
                    ALLOCATE( CPDESCA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CPDESCA', PROGNAME )
                END IF

                IF( .NOT. PFLAG .AND. 
     &              .NOT. AFLAG .AND. ALLOCATED( IDIUA ) ) 
     &              DEALLOCATE( IDIUA, IWEKA, XLOCAA, YLOCAA,  
     &                          STKHTA, STKDMA, STKTKA, STKVEA,  
     &                          CORISA, CBLRIDA, CPDESCA )

            CASE DEFAULT
                MESG = 'INTERNAL ERROR: Do not know about source ' //
     &                 'category ' // CATEGORY // ' in program ' //
     &                 PROGNAME
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END SELECT  ! select category

C.........  Sorted ...
        CASE( 'SORTED' )

            IF( UFLAG .AND. .NOT. ASSOCIATED( CIFIP ) ) THEN
                ALLOCATE( CIFIP( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CIFIP', PROGNAME )
            END IF

            IF( UFLAG .AND. .NOT. ASSOCIATED( TPFLAG ) ) THEN
                ALLOCATE( TPFLAG( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'TPFLAG', PROGNAME )
            END IF

            IF( UFLAG .AND. .NOT. ASSOCIATED( INVYR ) ) THEN
                ALLOCATE( INVYR( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'INVYR', PROGNAME )
            END IF

            IF( UFLAG .AND. .NOT. ALLOCATED( TZONES ) ) THEN
                ALLOCATE( TZONES( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'TZONES', PROGNAME )
            END IF

            IF( UFLAG .AND. .NOT. ASSOCIATED( NPCNT ) ) THEN
                ALLOCATE( NPCNT( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'NPCNT', PROGNAME )
            END IF

            IF( UFLAG .AND. .NOT. ASSOCIATED( CSCC ) ) THEN
                ALLOCATE( CSCC( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSCC', PROGNAME )  
            END IF

            IF( UFLAG .AND. .NOT. ASSOCIATED( CSOURC ) ) THEN
                ALLOCATE( CSOURC( NDIM1 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'CSOURC', PROGNAME )
            END IF

            IF( UFLAG .AND. .NOT. ASSOCIATED( IPOSCOD ) ) THEN  ! make sure
                ALLOCATE( IPOSCOD( NDIM2 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'IPOSCOD', PROGNAME )
            END IF

            IF( PFLAG .AND. AFLAG .AND. .NOT. ASSOCIATED(POLVAL) ) THEN
                ALLOCATE( POLVAL( NDIM2, NDIM3 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'POLVAL', PROGNAME )
            END IF

C.............  Deallocate for any source category
C.............  NOTE - do not deallocate CSOURC, CSCC, or IFIP, because they 
C               may be needed for reading day- and hour-specific data
            IF( .NOT. AFLAG .AND. .NOT. PFLAG ) THEN
                IF( ASSOCIATED( TPFLAG ) ) DEALLOCATE( TPFLAG )
                IF( ASSOCIATED( INVYR )  ) DEALLOCATE( INVYR )

            ELSE IF( .NOT. AFLAG ) THEN
                IF( ASSOCIATED( IPOSCOD ) ) DEALLOCATE( IPOSCOD )
                IF( ASSOCIATED( POLVAL  ) ) DEALLOCATE( POLVAL )
                IF( ASSOCIATED( NPCNT   ) ) DEALLOCATE( NPCNT )

            END IF               

            SELECT CASE( CATEGORY )
            CASE( 'AREA' )
            
                IF( UFLAG .AND. .NOT. ALLOCATED( XLOCA ) ) THEN
                    ALLOCATE( XLOCA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( YLOCA ) ) THEN
                    ALLOCATE( YLOCA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )
                END IF
                
                IF( UFLAG .AND. .NOT. ALLOCATED( CELLID ) ) THEN
                    ALLOCATE( CELLID( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CELLID', PROGNAME )
                END IF

            CASE( 'MOBILE' )
 
                IF( UFLAG .AND. .NOT. ALLOCATED( IRCLAS ) ) THEN
                    ALLOCATE( IRCLAS( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IRCLAS', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( IVTYPE ) ) THEN
                    ALLOCATE( IVTYPE( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IVTYPE', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( CVTYPE ) ) THEN
                    ALLOCATE( CVTYPE( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CVTYPE', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( CLINK ) ) THEN
                    ALLOCATE( CLINK( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CLINK', PROGNAME )
                END IF
  
                IF( UFLAG .AND. .NOT. ALLOCATED( XLOC1 ) ) THEN
                    ALLOCATE( XLOC1( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'XLOC1', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( YLOC1 ) ) THEN
                    ALLOCATE( YLOC1( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'YLOC1', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( XLOC2 ) ) THEN
                    ALLOCATE( XLOC2( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'XLOC2', PROGNAME )
                END IF
 
                IF( UFLAG .AND. .NOT. ALLOCATED( YLOC2 ) ) THEN
                    ALLOCATE( YLOC2( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'YLOC2', PROGNAME )
                END IF
 
                IF( .NOT. PFLAG .AND. 
     &              .NOT. AFLAG .AND. ALLOCATED( CVTYPE ) ) 
     &              DEALLOCATE( IRCLAS, IVTYPE, CVTYPE, CLINK, 
     &                          XLOC1, YLOC1, XLOC2, YLOC2 )

            CASE( 'POINT' )
 
                IF( UFLAG .AND. .NOT. ASSOCIATED( CISIC ) ) THEN
                    ALLOCATE( CISIC( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CISIC', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( IDIU ) ) THEN
                    ALLOCATE( IDIU( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IDIU', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( IWEK ) ) THEN
                    ALLOCATE( IWEK( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'IWEK', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( XLOCA ) ) THEN
                    ALLOCATE( XLOCA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'XLOCA', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( YLOCA ) ) THEN
                    ALLOCATE( YLOCA( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'YLOCA', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( STKHT ) ) THEN
                    ALLOCATE( STKHT( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'STKHT', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( STKDM ) ) THEN
                    ALLOCATE( STKDM( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'STKDM', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( STKTK ) ) THEN
                    ALLOCATE( STKTK( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'STKTK', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( STKVE ) ) THEN
                    ALLOCATE( STKVE( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'STKVE', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( CORIS ) ) THEN
                    ALLOCATE( CORIS( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CORIS', PROGNAME )
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( CBLRID ) ) THEN
                    ALLOCATE( CBLRID( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CBLRID', PROGNAME )  
                END IF

                IF( UFLAG .AND. .NOT. ALLOCATED( CPDESC ) ) THEN
                    ALLOCATE( CPDESC( NDIM1 ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'CPDESC', PROGNAME )  
                END IF

C.................  Do not deallocate plant description in case needed for
C                   point sources report.
                IF( .NOT. PFLAG .AND. 
     &              .NOT. AFLAG .AND. ASSOCIATED( CISIC ) ) 
     &              DEALLOCATE( CISIC, XLOCA, YLOCA, 
     &                          STKHT, STKDM, STKTK, STKVE )

            END SELECT

        CASE DEFAULT

            MESG = 'INTERNAL ERROR: Do not know about sorting type ' //
     &             SORTTYPE // ' in program ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END SELECT      ! sorted or unsorted

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE SRCMEM
