      MODULE LABELS
C
C     This module defines the data structure that is used to hold a 
C     list of the labels that are allowed to appear on the labeled 
C     input records in a Mobile6 input file.
C
C  Changes : (Last change first)
C
C  02 Sep 09 @EPA-bag bug22 Added ROADFE; allows ROADFE to be external data
C            file label as an alternate to MPG ESTIMATES
C  02 Jul 02 AIR Task 36: Added RELATIVE HUMIDITY and BAROMETRIC PRES, NLABS=197.
C  20 Mar 02 AIR Task 32: Added NO HIGH EMITTERS, OBD EXH RESPONSE,
C            OBD EVAP RESPONSE, removed OBD RESPONSE; NLABS=195
C  02 May 02 @EPA-bag Added HOT SOAK LQD LEAK; NLABS=193
C  06 Nov 01 AIR Task 27: Added MPG ESTIMATES (192)
C  31 Aug 01 AIR Task 18: Added ADDITONAL HAPS.
C  30 Aug 01 AIR Task 22: Added TOXICS SPREADSHEET (174)
C  03 Aug 01 AIR Tasks 15-21 Removed XXX command for END OF RUN (190)
C  02 Jul 01 AIR Task 15: Added Air Toxic OXYGENATE and RVP WAIVER commands
C            increased count to 190.
C  23 May 01 AIR Task 17: Added PM Speed Correction and Particulate EF,
C           increased count to 188.
C  10 May 01 AIR Task 11/5: Added PM Idle, increased count to 186.
C  30 Apr 01 AIR Task 12: Added Tier2 PM label, increased count to 185.
C  30 Jan 01 AIR Task 11: New labels for particulates and air toxics.
C  19 Jul 01 @EPA-djb Wish235, Add Average Speed command (173).
C  16 Apr 01 @EPA-EG increased NLABS to 172 for post 1996 I&M label
C  19 Jan 01 @EPA-bag increased NLABS to 171 for 2007 rule for heavy duty vehicles
C  20 Nov 00 @EPA-EG  New labels for clean air act
C  24 Oct 00 @EPA-BG  New labels for database age, model year and hour selection; NLABS
C            increased to 169
C  19 Oct 00 @EPA-djb Task W2, Increased NLABS to 166 for AGGREGATED OUTPUT label.
C  25 Oct 00 @EPA-djb Bug292, Added clarification of SECTION description.
C  27 Oct 00 @EPA-BG  BUG293 OBD REPAIRED EFS put back; NLABS increased to 165
C  06 Oct 00 @EPA-BG  BUG155&137 removed NEW EVAP TEST, IDLE EFS, INDUSTRY AVERAGE FUEL,
C            STARTS AS GM/ST, FEDERAL TEST PROC, NO RUNNING LOSS CF, OBD PERCENT EQUIPPED;
c            reduced NLABS to 164
C  15 Sep 00 AIR Task 08: Added EVAP I/M FRAC BEN. NLAB and GLL REDUCTION
C            NLABS from 173 to 175
C  15 Sep 00 @EPA-elg Added new label for T2 EVAP Phase-in inputs.
C  08 Aug 00 @EPA-elg Added new label to shut down evap calculations,
C            NLABS from 171 to 172.
C  11 Jul 00 AIR Task 06: Added NO EVP SAWTOOTH, FIXED P/P RATES,
C            HOT SOAK ACTIVITY, and DIURN SOAK ACTIVITY.
C  18 Jul 00 @EPA-bag BUG#54  Removed AMBIENT TEMPERATURE label; Changed
C                NLABS from 168 to 167
C  21 Jul 00 @EPA-RAG running loss; Reduced NLABS to 168
C  03 Jul 00 EPA@bag Task 04: Changed NLABS from 170 to 169 for AIR
C                task04
C  14 Jun 00 AIR Task 03: Added NO WEATHERED RVP and NO COLDCO STD.
C            Removed PHASE 2 CO STANDARDS, COLD CO CERTIFICATION DATA,
C            and NO CO OFFSET.
C   8 Jun 00 @EPA-ddj Bug120, Decreased NLABS to 174 (removed labels
C            associated with spreadsheet output).
C  23 May 00 @EPA-djb Added START DIST (180)
C  05 MAY 00 @EPA-hmc Increased NLABS to 179 for Daily Database Output
C   2 May 00 @EPA-ddj BUG7, decreased NLABS to 177 to account for removal
C            of SCENARIO TITLE label.
C  28 Apr 00 @EPA-elg Increased NLABS to 178 for SFTP.
C  17 Apr 00 @EPA-elg Increased NLABS to 176 for NGV.
C  10 Apr 00 AIR Task 02: Removed DATABASE POLLUTANST and abbreviations,
C            added POLLUTANTS, removed CALCULATE HC, CO, NOX ONLY
C  06 Apr 00 AIR Task 02: Added NO DESC OUTPUT and EXPAND EXHAUST.
C  20 Feb 00 @DynTel-bag 1-005 Increased NLABS to 177 for liquid leaks.
C  13 Jan 00 @DynTel-ddj 1-013  Removed ALHFLG CORRECTION FACTOR label
C            and added AC labels (increased NLABS to 174).
C  30 Dec 99 @DynTel-ddj 1-021 Increased NLABS to 171.
C  29 Nov 99 @DynTel-bag 1-006  Decreased NLABS to 169; hourly evaporative
C            emissions sent to output in CCEVRT.FOR
C   01 Nov 99 @Dyntel-nh 1-004 Increased NLABS to 170.
C   22 Sep 99 @DynTel-mla 2-684  Changed NLABS to 166 for new input and
C            database
C   7 Sep 99 @DynTel-ddj 2-710  Changed NLABS to 140 to accomodate
C            new FTP input options.
C  29 Jul 99 @DynTel-MLA 2-000  Removed module procedures CHKLAB and 
C            INILAB. They both are now simple external procedures, and
C            not module procedures.
C
C  Notes:
C
C     The derived type LABELDEF defines one variant form of an input
C     label. The structure contains the full and abbreviated forms of
C     the label text. It also contains a logical field, NEWLAB, which
C     is set to .TRUE. if there are no variant forms of the label. If
C     there are variant forms of a label, they should all be grouped
C     together in the LABELDEF list. (Values are assigned to LABELDEF
C     in the INILAB subroutine below.) The NEWLAB field should be set to
C     .TRUE. in the first record in the group and to .FALSE. in all the
C     remaining records in the group. The LABID field is set by the code.
C     (The INILAB subroutine uses the NEWLAB field to assign unique ID
C     numbers to different labels, but all of the variants of the same
C     label will be assigned the same ID number. The initial value of
C     the LABID field is not important.)
C
C     The SECTION field should be set to:
C
C       0 if the label can appear in the header portion of an input file.
C       1 if the label can appear only in the run portion of an input file.
C       2 if the label can appear either in the run or scenario portion
C         of an input file.
C       3 if the label can appear only in the scenario portion of an input file.
C
C     NOTE: The run data section is also called the "Global" data or
C           the "One-Time" data section.
C
C     The LABLIST array lists all of the input labels that can be used
C     in Mobile6 input files.
C
      IMPLICIT NONE
C
      TYPE :: LABELDEF
         CHARACTER(LEN=40) :: FULLAB
         CHARACTER(LEN=39) :: MINLAB
         LOGICAL           :: NEWLAB
         INTEGER           :: LABID
         INTEGER           :: SECTION
      END TYPE LABELDEF
C
      INTEGER, PARAMETER  :: NLABS = 198
C
      TYPE(LABELDEF), DIMENSION(NLABS) :: LABLIST
C
      END MODULE LABELS
