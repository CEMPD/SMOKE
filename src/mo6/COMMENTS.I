C  May 4, 2000-Subscripts were added to "comments.i" file and related to 
C              the routines with which they are associated.  Further 
C              explanation and elaboration of subscripts is necessary to 
C              make the file complete. 
C 
C 
C  MOBILE6 (xx-xxxember-98) is derived from MOBILE5b (14-September-96). 
C 
C  MOBILE6 is the sixth version of the FORTRAN program implementation of the 
C  Mobile Source Emission Model. It updates and replaces MOBILE5.  The model 
C  estimates HC, CO and NOx exhaust and HC evaporative emission factors by 
C  motor vehicles.  For information on using the program and on the 
C  differences between MOBILE5 and MOBILE6, consult the User's Guide to 
C  MOBILE5, and the User's Guide to MOBILE6.  Otherwise questions should be 
C  directed to: 
C 
C         Office of Mobile Sources 
C         Assessment and Modeling Division 
C         2565 Plymouth Road 
C         Ann Arbor, Mich. 48105 
C 
C 
C  MOBILE6/Standards     8/05/98 
C 
C 
C  Source Code Documentation 
C  ------------------------- 
C 
C  The source code itself is commented to assist the users who require a 
C  more detailed understanding of the program than provided in the Guides. 
C  This section provides the following information: 
C 
C  1) Execution summary 
C  2) Source code structure 
C  3) Source code comments structure 
C     a) Subroutines and functions header format 
C     b) BLOCK DATA subprograms header format 
C     c) In-line comments 
C  4) Subscript dictionary 
C  5) Subroutine / function parameter dictionary 
C 
C 
C  1) Execution summary 
C 
C  To run the program, the logical I/O units need to be attached to the 
C  appropriate value with respect to which platform the user/programmer is 
C  utilizing. 
C 
C    Platform                                      IBM   Macintosh 
C    -------------------------------------------------------------------- 
C    Input file          IOUGENerate                5       9 
C    Output file         IOUREPort                  6       6 
C    I/M credit file     IOUIMDeck                  4       4 
C    Error file          IOUERRor                   6       6 
C    Screen              IOUASK                     0       9 
C 
C  These I/O units need to be set in the routines DRIVER. Scan for 
C  'Change for different computers'. 
C 
C  Logical I/O Unit 3 is attached to the I/M credit file Tech12.D in all 
C  platforms. 
C 
C  Refer to the User's Guide to MOBILE for more detail on setting up runs 
C  and check your operating system documentation for the system commands for 
C  compiling and executing the program.  There is discussion of running the 
C  program on the IBM PC, Macintosh, and Mainframe in the User's Guide. 
C 
C 
C 
C  2) Source code structure 
C 
C     This section describes the layout of the code plus the numbering, 
C     naming and programming conventions and guidelines. 
C 
C     The source code for the subprograms is grouped as follows: 
C 
C        MAIN (DRIVER) 
C        MOBILE subroutine 
C        QUITER subroutine 
C        input subroutines 
C        computation subroutines & functions 
C        pointer functions 
C        output subroutines 
C        BLOCK DATA subprograms 
C 
C     Since MOBILE was written as a callable subroutine for release MOBILE4.1, 
C     the program has been susceptible to batch mode errors.  A batch mode error 
C     occurs when a variable is altered in one run and used again in a 
C     subsequent run without reinitialization.  Any variable that is not 
C     static must be reinitialized prior to the next run.  Otherwise, the 
C     variable holds the incorrect value. 
C 
C     Block data is suitable for holding values of static variables.  However, 
C     if the data need to be reinitialized, block data is sometimes unnecessary 
C     and redundant.  It may be more appropriate and efficient to have the 
C     variables reinitialized by a calling module at the beginning of the 
C     program or just before the variable is used rather than using block data. 
C 
C     Statement number ranges are defined for usage as follows: 
C 
C          1 -   9  -  not used 
C         10 -  98  -  branching 
C               99  -  branch to last RETURN / RETURN1 / STOP before END 
C        100 - 199  -  READ formats 
C        200 - 299  -  WRITE prompt and output formats 
C        300 - 499  -  WRITE error and warning messages 
C        500+       -  not used 
C 
C     Statement numbers are local.  For example, within each subprogram, 
C     branch numbers start at 10, READ formats start at 100, etc.  When ever 
C     possible, consecutive statement numbers should be incremented by 10, 
C     5 if necessary. 
C 
C     Typed entities should be explicitly declared.  Beginning with MOBILE6, 
C     DO NOT use implicit typing.  This applies to variables, functions, etc. 
C 
C     There should be only one RETURN statement per routine.  Unless there is 
C     a need for modified RETURNs (eg. RETURN1). 
C 
C     Functions should have their value set in one place only; before the last 
C     RETURN statement. 
C 
C     GOTO statements should be avoided when possible.  A GOTO 99 where 99 is a 
C     RETURN statement at the end of the program is preferable to multiple 
C     RETURNS. 
C 
C     Flags should adhere to the value convention of 1=NO and 2 or greater=YES. 
C 
C     IF THEN ELSE statements should be used instead of multiple FORTRAN 66 
C     type IF statements when testing the value of one token. 
C 
C     DO loops should end with a CONTINUE or END DO statement. 
C 
C     Indenting IF THEN ELSE statements DO loops etc. is encouraged for easier 
C     comprehension. 
C 
C     Given FORTRAN 77's 6 character maximum on subprogram name lengths, useful 
C     naming conventions were hard to come by.  Beginning with MOBILE6, which 
C     will use more Fortran 90 features, the 6 character limit has been expanded 
C     to 31.  However, given the DOS restriction of 8 characters for a filename, 
C     one must restrict subprogram and common block names to 8 characters. 
C     There was also the limitation posed by not wanting to change arbitrarily 
C     all the subprogram names of the predecessor MOBILE programs.  Hence, the 
C     only conventions are: 
C 
C     a) Input driver sections all have the suffix 'SEC'. 
C     b) Input routines should use the prefix 'GET'. 
C     c) Pointer functions all begin with 'I' (to get default INTEGER typing). 
C     d) Pointer functions all contain/end with 'PT' or 'PTR'. 
C     e) Output subroutines all have the prefix 'OUT'. 
C 
C     Otherwise, new subprogram names were picked to indicate the purpose of 
C     the code. 
C 
C  Note: some features during the course of the model development were 
C  subsequently removed or "turned off".  The latter case would occur when 
C  the prospect of restoring the feature for in-house MOBILE4 or for MOBILE5 
C  seemed likely.  Commented out code to this end is done with a 'CC'. 
C 
C 
C  3) Source code comments structure 
C 
C     Comments should be professional.  They should be meaningful to someone 
C     who has never read the code before.  Pseudo computer language should be 
C     avoided.  Comments should be made in English using upper and lower case. 
C 
C     Comments that are in the body of the code should be consistently placed, 
C     before the line of code that they refer to. 
C 
C     Each subprogram contains a documentation header, set up in a standard 
C     format as described below.  All variables and arrays used in the code 
C     are defined in the comments.  The subroutine / function headers define 
C     local parameters and local variables / arrays.  The BLOCK DATA headers 
C     describe the COMMON storage areas.  Array subscripts and subprogram 
C     parameters are specified in dictionaries following this section.  The 
C     intent is to allow the user to scan the source code on-line for 
C     definitions, rather than have to leaf through hard-copy documentation 
C     such as the Programmer's Guide to MOBILE4. 
C 
C     For readability all comments should be made to fit within 80 columns. 
C 
C     Changes should be signed and dated in the module in which the 
C     change occurred. The format for signing and dating is: 
C     Day-Month-Year @ DynTel-fml (first middle last initial) Request # 2-nnn. 
C     This should be at the top of the "Changes:" section with the newest 
C     entries at the very top of the section. 
C 
C     a) Subroutines and functions header format 
C 
C        <subprogram name> 
C        C 
C        C  <short description of purpose> 
C        C 
C        C  Called by <subprogram name(s)> 
C        C 
C        C  Calls <subprogram name(s)> 
C        C 
C        C  Changes: 
C        C  <date and signature stamp of latest change ie. Day-Month-Year @ 
C        C  DynTel-fml (first middle last initial) Request # 2-nnn> 
C        C 
C        C  Input on call: 
C        C 
C        C    constant(s): <names> 
C        C    agrument list: <names> 
C        C    common block(s): 
C        C    /<cb name>/ <names> 
C        C         .         . 
C        C         .         . 
C        C 
C        C  Output on return: 
C        C 
C        C    argument list: <names> 
C        C    common block(s): 
C        C    /<cb name>/ <names> 
C        C         .         . 
C        C         .         . 
C        C 
C        C  Local array subscripts: 
C        C 
C        C  <name>(<#>,...)  -  <name>(<name>,...) 
C        C         .                    . 
C        C         .                    . 
C        C 
C        C  Local variable / array dictionary: 
C        C 
C        C   Name   Type              Description 
C        C  ------  ----  -------------------------------------------------- 
C        C  <name> <type>         <short description> 
C        C    .      .                     . 
C        C    .      .                     . 
C        C 
C        C  Notes: 
C        C 
C        C  <general notes on subprogram> 
C        C 
C        C 
C            IMPLICIT NONE 
C            <code> 
C 
C        Notes on format: 
C 
C        A local variable is one not passed in via the argument list 
C        or a common block.  It can be used in the argument list of a 
C        subprogram call made by the defining subprogram. 
C 
C        Also, if a section of the header is not needed, its format is 
C        omitted.  For example, if no subprograms are called, then the 
C        'C  Calls <names>' line is not included.  "Notes" are always included. 
C 
C        Finally, 'Output on return' includes any variable / array whose 
C        contents can be modified by the subprogram.  This usually will not 
C        include data that can be changed by subprograms called by the 
C        subprogram, since the comment headers on those subprograms will 
C        indicate what is being passed back up. 
C 
C     b) BLOCK DATA subprogram header format 
C 
C              BLOCK DATA 
C        C 
C        C  BLOCK DATA subprogram <#> 
C        C 
C        C  <short description of data initialized in this BLOCK DATA> 
C        C 
C        C  Changes: 
C        C  <date and signature stamp of latest change ie. Day-Month-Year @ 
C        C  DynTel-fml (first middle last initial) Request # 2-nnn> 
C        C 
C        C  Common block array subscripts: 
C        C 
C        C  <name>(<#>,...)  -  <name>(<name>,...) 
C        C         .                    . 
C        C         .                    . 
C        C 
C        C  Common block dictionary: 
C        C 
C        C   Name   Type              Description 
C        C  ------  ----  -------------------------------------------------- 
C        C  /<cb name>/: 
C        C  <name> <type>         <short description> 
C        C    .      .                     . 
C        C    .      .                     . 
C        C 
C        C  Local array subscripts: 
C        C 
C        C  <name>(<#>,...)  -  <name>(<name>,...) 
C        C         .                    . 
C        C         .                    . 
C        C 
C        C  Local array dictionary: 
C        C 
C        C   Name   Type              Description 
C        C  ------  ----  -------------------------------------------------- 
C        C  <name> <type>         <short description> 
C        C    .      .                     . 
C        C    .      .                     . 
C        C 
C        C  Notes: 
C        C 
C        C 
C            IMPLICIT NONE 
C 
C            COMMON /<name 1>/ ... 
C                 . 
C                 . 
C                 . 
C        C 
C              DIMENSION (or type) <array 2>(subscript(s)) 
C                 . 
C                 . 
C        C 
C              EQUIVALENCE (<array 1, subscripted>,<array 2, subscripted>) 
C                 . 
C                 . 
C                 . 
C        C 
C        C  /<name 1>/: <any comments on <name 1>'s labels, purpose or data 
C        C               or on the structure of the initializing DATA 
C        C               statement(s)> 
C        C 
C              DATA <var/arr name 1>/<data>... 
C                 . 
C                 . 
C                 . 
C                 .../ 
C        C 
C        C  /<name 2>/: ... 
C                 . 
C                 . 
C                 . 
C              END 
C 
C        Notes on format: 
C 
C        Local arrays are defined and then equivalenced to common block arrays 
C        that are too large be initialized in 1 DATA statement.  Too large = 
C        would require greater than 39 continuation lines, the Fortran 90 limit. 
C 
C        The local arrays dictionary is often compressed by providing a 
C        generic definition of the conventions used to construct and name the 
C        set of local arrays to be equivalenced to a given common array. 
C 
C        "/<name #>/:" prefix is included only if a large number of common 
C        blocks are being initialized in a single BLOCK DATA (see BD 23). 
C 
C 
C     c) INCLUDE file header format 
C 
C        Include files were introduced to the model after the release of 
C        MOBILE5a.  Include files should hold the description of the common 
C        block variables rather than block data. 
C 
C        C 
C        C  The data is initialized in BD<name or #> 
C        C 
C        C  Notes/Changes: 
C        C 
C 
C 
C     d) In-line statements 
C 
C        In-line comments are intended to describe the algorithm being used 
C        by the code.  Where MOBILE4 is identical to MOBILE3 in operation, 
C        the MOBILE3 comments have usually been left intact.  Where the 
C        method has been updated or the whole section / subprogram is new, 
C        the in-line comments describe the MOBILE4 methodology.   In-line com- 
C        ments should be preceded by a line with just a 'C' and followed by a 
C        line with just a 'C'. 
C 
C        Particular attention is paid to the use of indices and pointers. 
C        Some of the ambiguities of MOBILE3, such as the use of IDX vs JDX, 
C        have been cleared up.  The usage of new keys and pointers is specified 
C        in detail where they occur. 
C 
C 
C  4) Subscript dictionary 
C 
C     Subscript variables used to index arrays in MOBILE4 are globally 
C     defined, not in the sense that the values are available to all 
C     subprograms, but rather that the same label is used for the same 
C     subscript function throughout the code.  IP is the 3 pollutants index 
C     for all of MOBILE4, although its instances are usually as a local 
C     variable.  As array subscripts, all these variables must be of type 
C     INTEGER. The default INTEGER typing syntax rule is employed by having 
C     all subscripts begin with a letter from the range I through N. 
C 
C      Name    Range               Description 
C     ------  -------  ---------------------------------------------------- 
C     I         1-3    (ossall)  
C     IAD       1-5    index into ID dimension of ATR for ADJFAC equation 
C 
C     IADJ      1-2    index into CSIZE adjustment vector ADJFAC 
C 
C     IAGE1     1-20   (rdt12c) 
C 
C     IAIM      0-MAXYRS(initim) 
C 
C     IATP      1-2    (iniflg) 
C 
C     IAY       1-2    anti-tampering program start year flag: 
C                        1 = in year range up thru start year 
C                        2 = in year range after start year 
C 
C     IB        1-3    bag number: 
C                        1 = cold start 
C                        2 = hot stabilized 
C                        3 = hot start 
C 
C     IBER      1-3    terms of user supplied new emission rate equation: 
C                        1 = zero mile level (intercept) 
C                        2 = deterioration rate (slope) 
C                        3 = deterioration rate (slope), LDGV 81+ HC/CO 50K+ 
C 
C     IBIN      1-12   (gett2frc) 
C 
C     IC        1-7    overlap effects categories computed direct from base 
C                      tamper rates; each case is the set intersection of 
C                      the indicated disablement types: 
C                        1 = AorE * CATS 
C                        2 = AorE        * NCKS 
C                        3 = AorE               * TNKS 
C                        4 = AIRS * CATS * NCKS 
C                        5 = AIRS * CATS        * TNKS 
C                        6 =        CATS * NCKS 
C                        7 =        CATS        * TNKS 
C                      AorE = AIRS if HC or CO 
C                           = EGRS if NOx 
C 
C               8-12   no overlap effects categories computed from the 7 
C                      overlap cases: 
C                        8 = AIRS 
C                        9 = CATS 
C                       10 = NCKS 
C                       11 = TNKS 
C                       12 = <no tampering effect> 
C 
C     ICAT      1-3    catalyst technology configurations: 
C                        1 = ox catalyst / no catalyst 
C                        2 = 3 way catalyst 
C                        3 = 3 way + ox catalyst 
C 
C     ICH       1-20   character string array position index: 
C                        1 = 1st 4/8 characters (CHARACTER*4/8) 
C                        2 = 2nd 4/8 characters 
C                                . 
C                                . 
C                      also serves to subscript any other CHARACTER type array. 
C 
C     ICOEF2    1-6    coefficient number of 2nd degree polynomial 
C                      (0th (intercept) to 2nd decree term coefficients) 
C 
C     ICOEF5    1-6    coefficient number of 5th degree polynomial 
C                      (0th (intercept) to 5th decree term coefficients) 
C 
C     ICOL      1-3    index into EFF's column dimension (BD 24 defines abbrs): 
C                      1 = PREV CAT ; 2 = PREV ; 3 = SUBS 
C 
C     ICOLD     1-2    index into cold CO standard case of MYCOL & CCSTD 
C 
C     ICS       1-NT   secondary index: selects pointers to category 
C                      sizes to be summed up to get SUB (in C8TO11) and 
C                      DIV (in A8TO11) 
C 
C     IC11      1-11   positive tampering effects categories 1-11; 12 is the 
C                      no effects case. 
C 
C     IC1       1-11   (initap) 
C 
C 
C     IC2       1-11   (initap)  
C     IC7       1-7    overlap tampering effects categories 1-7 
C 
C     ID        1-9    disablement types (1-5,9 for bag, 6-8 for evap): 
C                        1 = AIRS = air pump disabled 
C                        2 = CATS = catalyst removed 
C                        3 = NCKS = misfueled - filler neck 
C                        4 = TNKS = misfueled - tank (all but filler neck) 
C                        5 = EGRS = EGR disabled 
C                        6 = EVAP = evaporative system tampered - cannister only 
C                        7 = PCVS = PCV disabled 
C                        8 = CAPS = evaporative system tampered - cannister and 
C                                   gas (fuel inlet) cap 
C                        9 = MISF = misfueled - all sources = NCKS + TNKS 
C 
C     IDBAG     1-5    bag exhaust disablement types (see ID) 
C 
C     IDPB      1-3    misfueling related disablement types: 
C                        1 = NCKS = misfueled - filler neck 
C                        2 = TNKS = misfueled - tank (all but filler neck) 
C                        3 = MISF = misfueled - all sources = NCKS + TNKS 
C 
C     IDU       1-5    index for Diurnal (DU) indicating full, 8AM-11AM, 
C                      10AM-3PM, 8AM-2PM, and multiple 
C 
C     IDUSER    1-8    ID cases & order of cases presented to user in TAMFLG 
C                      (TAMZML & TAMDR) & ATPFLG (DISTYP) options.  IDUSER 
C                      differs as follows: 4 = MISF.  GETTAM derives TNKS = 
C                      MISF - NCKS for TAMDR & TAMZML and stores the read and 
C                      computed values in ID order. ATPAER does not need TNKS 
C                      & uses DISTYP in the read in IDUSER order. 
C 
C     IDX       1-25   model year (my) window index: 
C                        1 = 25+ my age group cell => tamper effect is 
C                            always 0.0 => this case skipped in TAMPER (only) 
C                        2 = 24th my cell 
C                        3 = 23rd my cell 
C                          . 
C                          . 
C                        . 
C                       24 = 2nd my cell 
C                       25 = 1st my cell = calendar year (ICY) 
C 
C     IECT      1-2     (soakfcn) 
C 
C     IEG       1-7    tampering effects group: 
C                        1 = air pump disabled, air pump only equipped 
C                        2 = air pump disabled, air pump/catalyst equipped 
C                        3 = catalyst removed 
C                        4 = misfueled catalyst 
C                        5 = EGR disabled 
C                        6 = EGR disabled / catalyst removed 
C                        7 = EGR disabled / misfueled 
C 
C     IEK       1-3    catalyst emission effects groups: IEK = IEG - 1 
C 
C     IELST     4-7    last non-zero tampering effects group size, given IPG 
C 
C     IEQU      1-20   pointer to idle equation to be zeroed out due to entry 
C                      of new exhaust ef equation by user 
C 
C     IERNEW    1-20   pointer to location in default ef rates & model years 
C                      arrays where new ef equation's parameters are to be 
C                      stored (IERNEW = 0 => no space for new equation) 
C 
C     IETP      1-7    (getetp) 
C 
C     IETPY     1-NETPYR(etpcal) 
C 
C     IEVEQ     1-9    HS and DU equation term coefficient - see definitions 
C                      for HSEQ in BD 09 and DUEQ in BD 11 
C 
C     IEVP      1-4    index for gram evap rates EVP, GREVP and VGREVP: 
C                        1=hot soak 
C                        2=in-use diurnal 
C                        3=multiple diurnal 
C                        4=crankcase 
C 
C     IEVX      1-32   model year (my) window index for evaporative data input, 
C                      i.e. to the IEVXth myg's rate(s) for a given N<> set 
C 
C     IE1ST     1-3    first non-zero tampering effects group size, given IPG 
C 
C     IFAC      1-100  index into new (user supplied) emission rates arrays; 
C                      the user can enter up to 100 replacement cases. 
C 
C     IFDS      1-3    fuel delivery system type: 
C                        1 = carbureted      2 = fuel injected    OR 
C                        1 = CARbureted      2 = TBI      3 = PFI 
C 
C     IFG       1-13   flag names for prompts (listed by common block used): 
C                         /FLAGS1/      /FLAGS2/       /FLAGS3/      /FLAGS4/ 
C                        1 = TAMFLG    4 = MYMRFG     8 = ATPFLG   12 = PRTFLG 
C                        2 = RNLFLG    5 = NEWFLG     9 =          13 = 
C                        3 = VMFLAG    6 = IMFLAG    10 = RLFLAG   14 = NMHFLG 
C                                      7 = ALHFLG    11 = OUTFMT   15 = HCFLAG 
C                      The 14th flag (or 0 = PROMPT) is read without a prompt. 
C 
C     IFL       1-2    myg bounds:  1 = First year    2 = Last year 
C 
C     IFUEL     1-3    fuel type:   1 = Gasoline, 2 = Ether Blend, 
C                                   3 = Alcohol Blend 
C     IFUEL_NUM1(?)    (sets the lower limit to the IFUEL do loop                                         index,allows.57"d:) 	 
C     IFUEL_NUM2(?)    (sets the upper limit to the IFUEL do loop                          			    index,allows.57"d:)  
C     IFUEL_SYS 1-7    (dirnl_hr)  
C 
C     IG        1-MAXOXY(hcoxy) 
C 
C     IG<i>     1-20   model year group (myg) pointer: 
C                      <i> = null string,C,CAT,D,E,F,L,LM,S,T,U,W,1,2,3,4 
C                      IG links a data array to its corresponding myg index 
C                      array. These data arrays all have my as a subscript 
C                      plus have cases of contiguous subsets of years with 
C                      the same data point.  Core space is then saved by 
C                      storing each point only once per my subset and then 
C                      constructing an index array to enable the program to 
C                      find the location IG<i> of the data point for a given 
C                      my. 
C 
C     IGCSF     1-3    model year group index in the array CSF1ST: 
C                        1 = pre-1981 model years 
C                        2 = 1981-1983 model years 
C                        3 = 1984+ model years 
C 
C     IGD       1-3    MYG index into 3rd dimension of TAMZML & TAMDR 
C                        1 = pre-1981  2 = 1981-83  3 = 1984+ for LDGV/T 
C                        1 = 1981+ for HDGV 
C     IGER      1-20   myg pointer into basic emission rate arrays from IERPTR 
C 
C     IGFID     1-MAXFID(nmogfid) 
C 
C     IGIDL     1-15   myg pointer into idle emission rate arrays from IDLPTR 
C 
C     IG50      1-12   myg pointer into the 50K+ deterioration rate arrays 
C                      ERB50K & ERU50K (the arrays differ in lookup algorithm) 
C 
C     IGTL      1-4    running loss HC: emission rate MYG pointer 
C 
C     IH        1-3    non-exhaust HC emission case: 
C                        1 = hot soak evaporative 
C                        2 = diurnal evaporative 
C                        3 = crankcase 
C 
C     IHG       1-2    non-exhaust HC technolo00gy / tampering case: 
C                        1 = evaporative 
C                        2 = crankcase 
C 
C     IHDSAL    1-2    HDGV sales fraction case 
C                        1 = HDGV weighing < 14000 lbs 
C                        2 = HDGV weighing >= 14000 lbs 
C     IHOUR     1-24   (hccalx)  
C 
C     IHRD      1-6    evaporative temperature types: 
C                        1 = hot soak 
C                        2 = running loss 
C                        3 = diurnal, 8 AM 
C                        4 = diurnal, 10-11 AM 
C                        5 = diurnal, 2 PM 
C                        6 = resting loss 
C 
C     IHR       1-14    (rnloss) 
C 
C     IHR_SOAK  1-72    (dirnl_hr) 
C 
C     IK        1-2     (hccalx)  
C 
C     IL        2-4     (?)  
C 
C     ILU       1-2    bounds of a value range (ex.: tags of nonzero 
C                      category sizes): 
C                        1 = lower bound 
C                        2 = upper bound 
C 
C     IM        1-2    I/M program flag: 
C                        1 = I/M program not in effect for MY x IV x ICY case 
C                        2 = I/M program is in effect for MY x IV x ICY case 
C 
C     IMIN      1-60   (Hotsok) 
C     
C     IMINPUT   2-6    (TESTIM)  
C 
C     IMPGM     MAXIMPGM-1(improg)  
C 
C     IMPNUM    1-2    I/M program number in effect 
C                        1 = first I/M program is in effect. 
C                        2 = second I/M program is in effect. 
C 
C     IMDL      1-2    coefficients for trips per day model equation: 
C                        1 = intercept 
C                        2 = slope 
C     IMLOG_NP  RSD and COIM corrections do not apply unless IMLOG_NP>0 
C 
C      
C     IMSKIP    1-400 (rdt12c)   
C 
C     IMVG      1-4   (reinit) 
C 
C     INH       1-2   (hrloop) 
C 
C     INO       1-2    index into NOYES character string vector: 
C                        1 = No     2 = Yes 
C 
C     IOPT      1-300 (EVWTPRE)  
 
C 
C     IP        1-3    pollutants: 
C                        1 = HC  = hydrocarbon 
C                        2 = CO  = carbon monoxide 
C                        3 = NOx = oxides of nitrogen 
C 
C     IPAR      1-5    parameters identifying where and when to apply user 
C                      supplied new emission rates: 
C                        1 = region 
C                        2 = vehicle type     4 = first my of range covered 
C                        3 = pollutant        5 = last my of range covered 
C 
C     IPARTL    1-3    index for Partial Diurnal (PARD) indicating 
C                      8AM-11AM, 10AM-3PM and 8AM-2PM        
C 
C     IPF       1-2    (iniflg) 
C 
C     IPG       1-2    pollutants grouped for tampering effects & SK2 scfs: 
C                        1 = HC and CO 
C                        2 = NOX 
C 
C     IPH       2-3    (Local) 
C 
C     IPIM      1-3    pollutants affected by I/M programs: 
C                        1 = HC 
C                        2 = CO 
C                        3 = NOx 
C 
C     IPIV      1-3    (Getcol) 
C 
C     IPOXY     1-2    (hcoxy) 
C 
C     IPPS      (?)    (DO loop index for the three pass/fail pressure and purge                          tests)                         
C     IPPS_NUM  (?)    (Upper limit of IPPS DOloopfor the three pass/fail                          pressure and purge tests)  
C 
C 
C     IP1       1-3    lower bound of output IP loop (first pollutant to 
C                      have results printed) 
C 
C     IP2       1-3    upper bound of output IP loop (last pollutant to 
C                      have results printed) 
C 
C     IP50      1-2    pollutant types affected by 1981+ LDGV 50K+ change in 
C                      deterioration rate 
C 
C     IQG       1-6    emission control equipment group; 2 more groups 
C                      computed from the stored 6 groups data: 
C                        1 = air pump only 
C                        2 = air pump / catalyst 
C                        3 = catalyst only 
C                        4 = EGR only 
C                        5 = EGR / 3-way catalyst 
C                        6 = no EGR / 3-way catalyst 
C 
C     IR        1-2    region code: 
C                      IR has same value range as IR00EJN,  
C                      but is not necessarily  
C                      the scenario field's value 
C                        1 = low altitude 
C                        2 = high altitude 
C 
C     IREC      1-999  input record number, on a given READ 
C 
C     IRNG      1-7    (Enforc)  
C 
C     IROAD     1-5    (Hrloop) 
C 
C     IROW      1-15   index into EFF's row dimension = the 8 ID cases split 
C                      into 15 subcases (see BD 24 definition of EFF) 
C 
C     IRT       1-4    temperature cases in RVP cf calculations, the values 
C                        referenced depending on the (sub)array.  See /RVPEX1/. 
C 
C     IRQ1      1-3    RVP cf equation parameters set 1 (RVP > 9.0 adj.): 
C                        RADJCF=(b1+b2*RVP)/b3 
C 
C     IRQ2      1-2    RVP cf equation param0eters set 2 (T interpolation) 
C                          1 = initial temperature (T range lb) 
C                          2 = temperature difference (interval length) 
C 
C     IRQ3      1-3    RVP & high T combined cf formula coefficients: 
C                        TCF(IB)=EXP(a1*(RVP-9)+a2*(T-75)+a3*(RVP-9)*(T-75)) 
C 
C     IRVP      1-2/4  types of non-certification RVP: 1 = Base, 2 = In-use 
C                                             3 = Uncontrolled refueling loss 
C                                             4 = Hot soak (weathered) 
C 
      IS        1-4    overlap tampering effects category sizes 8-11 
C                      ( IS = IC - 7 ); used in A8TO11. 
C 
C     ISC       1-5    underscore string repetition counter - for example, 
C                      ISC = 3 => repeat the string 3 times 
C     ISEG      JAN_SEP, OCT_DEC (improg) 
C 
C     ISOAK     1-20   (Dirnltem) 
C 
C     ISP       1-21   my group pointer into speed correction equation array; 
C                      found by ISPPTR; indexes SK1 & SCUNA1. 
C     ISR       1-2    (Bigin2)                                                                                                                  C       
C     ISTD      1-10   (hrloop) 
C 
C                      1994 and newer model years are split into 11 emission standard 
C                      groupings.  These groupings are: 
C 
C                      ISTD = 1 : Tier 0 
C                             2 : Intermediate Tier 1 
C                             3 : Tier 1 
C                             4 : Tier 2 
C                             5 : Intermediate Transitional Low Emission Vehicle 
C                             6 : Transitional Low Emission Vehicle 
C                             7 : Intermediate Low Emission Vehicle 
C                             8 : Low Emission Vehicle (LEV) 
C                             9 : Transitional Ultra Low Emission Vehicle 
C                            10 : Ultra Low Emission Vehicle (ULEV) 
C                            11 : Zero Emission Vehicle (ZEV) 
C 
C 
C     ISTR      1-5    (rdt12c) 
C 
C     IS2a      1-2/3  IS2a      1-2/3 number of coefficients in a new scf                         equation: 
C     IS2C = 2 (A/X+B);  IS2N = 3 (A+B*X+C*X**2) 
C 
C     IT        1-9    pointer to temperature correction equation 
C                      coefficients; found by ITCPTR; indexes TTA/4/7COF 
C                      arrays. 
C 
C     ITECH     1-2    (rdt12c) 
C 
C     ITECH<a>  1-15   emission control technology class wrt I/M impact; 
C                      f(MY,IVLD) as determined by IMPTR; indexes I/M credits 
C                      arrays CR12HC/CO & CRD4VA/B. ITECHA/B/C are in BD8&34. 
C 
C     ITEM      1-2    NAMTEM index: 1 = 'mini', 2 = 'maxi' prefixes 
C 
C     ITER      1-4    unspecified string iterations counter; same usage as 
C                      ISC, only with any string 
C 
C     ITL       1-13   running loss HC: defined trip lengths for 3 cycle types: 
C                          1-6 - NYCC      7-12 - LA-4      13 - HFET 
C 
C     ITLR      1-4    running loss HC: RVPs at which measured rl HC available 
C 
C     ITLT      1-4    running loss HC: temperatures at which measured rl HC 
C                                       emission rates are available 
C     ITN       1-NTF  (initim) 
C 
C     ITY       1-3    type tampering rates wrt I/MxIGD table heading substring 
C                      pointer 
C 
C     IUDI      1-5    index into Uncontrolled Diurnal rate UDI: 
C                        1 = standard FTP conditions 
C                        2 = in-use conditions 
C                        3 = 8AM-11AM conditions 
C                        4 = 10AM-3PM conditions 
C                        5 = 8AM-2PM conditions 
C 
C     IV        1-8    vehicle types evaluated by MOBILE4: 
C                        1 = LDGV  = light duty (ld) gasoline fueled vehicle 
C                        2 = LDGT1 = ld gasoline fueled truck, <= 6000 lbs. 
C                        3 = LDGT2 = ld gasoline fueled truck, 6001-8500 lbs. 
C                        4 = HDGV  = heavy duty gasoline fueled vehicle 
C                        5 = LDDV  = light duty diesel fueled vehicle 
C                        6 = LDDT  = light duty diesel fueled truck 
C                        7 = HDDV  = heavy duty diesel fueled vehicle 
C                        8 = MC    = motorcycle 
C 
C     IV<i>     1-6    vehicle type subgroups; index factor arrays; 
C                      <i> = LM, TEMP, A, B or LDG, defined as: 
C                         IVLM        IVTEMP      IVA         IVB       IV        C                      correction LDG 
C                        1 = LDGV    1 = LDGV    1 = LDGV    1 = HDGV  1 = LDGV 
C                        2 = LDGT1   2 = LDGT1   2 = LDGT1   2 = LDDV  2 = LDGT1 
C                        3 = LDGT2   3 = LDGT2   3 = LDGT2   3 = LDDT  3 = LDGT2 
C                        4 = LDDV    4 = HDGV    4 = MC      4 = HDDV 
C                        5 = LDDT    5 = MC 
C                        6 = MC 
C                      indexes: 
C                        MYGSP1      MYGTCF      SCADJ1      GPBSCO    TTFC 
C                        BFRCOF                  ISPGRP 
C                                                ALHRET 
C                      See subroutines BIGCFX & BIGALH and function ISPPTR. 
C     IVC       1-16   (Vehcount) 
C 
C     IVCUT     1-4    ( 
C 
C     IVD       1-30   (Vehcount) 
C 
C     IVEFF     1-2    vehicle class groupings for Stage II vrs percent 
C                      efficency in reducing HC refueling losses: 
C                        1 = LDGV, LDGT1 & LDGT2    2 = HDGV 
C 
C     IVEXP     1-4    (Reinit) 
C 
C     IVF1      1-6    vehicle class types for MOBILE1 numeric output 
C                      format (OUTFMT = 1): 
C                        1 = LDV      3 = LDT2     5 = HDD 
C                        2 = LDT1     4 = HDG      6 = MC 
C 
C     IVF2      1-9    vehicle class types for MOBILE2 numeric output 
C                      format (OUTFMT = 2): 
C                        1 = LDGV     4 = LDGT     7 = LDDT 
C                        2 = LDGT1    5 = HDGV     8 = HDDV 
C                        3 = LDGT2    6 = LDDV     9 = MC 
C 
C     IVF3      1-10   vehicle class types for MOBILE4 descriptive output 
C                      format (OUTFMT = 3): same as IVF2 plus 10 = All Vehicles 
C 
C     IVIM      1-4    vehicle class types that can be included in MOBILE4 
C                      inspection / maintenance (I/M) programs: 
C                        1 = LDGV    2 = LDGT1    3 = LDGT2    4 = HDGV 
C 
C     IVHD      1-2    vehicle class types for heavy duty only (used for 
C                      units conversion of user entered hd zml & dr): 
C                        1 = HDGV     2 = HDDV 
C     IVLDG     1-5    (Getcol) 
C 
C 
C     IVP       1-3    vehicle types whose ef's can be modified by correction 
C                      factors for extra load, trailer towing &/or CO offset: 
C                        1 = LDGV     2 = LDGT1    3 = LDGT2 
C 
C     IVREV     1-4    vehicle types covered by the revised idle EF algorithm: 
C                        1 = LDGV     2 = LDGT1    3 = LDGT2    4 = HDGV (HC&CO) 
C                      IVI = IVREV, except 3 & 4 -> 3 (HDGV uses LDGT2 data) 
C 
C     IVSCAN    1-MAXVEH(bigmet) 
C 
C     IVTAM     1-4    vehicle types impacted by tampering: 
C                        1 = LDGV     2 = LDGT1    3 = LDGT2    4 = HDGV 
C     IVTEMP    1-14   (iniscn) 
C 
C 
C     IV9       1-9    vehicle type categories for output arrays EFEVAP, 
C                      EFIDLE, & EFFTP.  First 8 are identical to IV's. 
C                        9 = LDGT 
C                          = LDGT1 + LDGT2, each weighted by its vmt share 
C 
C     IWK       IWK_DAY-IWK_END(iniscn) 
C 
C     IY        IVAGEIM (ISEG), JDX-1,7(improg) 
C 
C     IYR       9-MAXLYR (no_tier2) 
C 
C     I2ND      1-4    secondary index (pointer into a pointer vector) 
C 
C     J         1-3    (reinit)  
C     J1        1-25   (reinit)  
C     JCS       1-7    pointer to category size to be summed into DIV 
C                      (in A8TO11) or (* TGSUSE) into EGS (in SETEGS) 
C 
C     JDX       25-1   model year window index = age of vehicle relative to 
C                      ICY + 1. In detail: 
C                       25 = 25+ my age group cell => tamper effect is 
C                            always 0.0 => this case skipped in TAMPER (only) 
C                       24 = 24th my cell 
C                       23 = 23rd my cell 
C                          . 
C                          . 
C                . 
C                        2 = 2nd my cell 
C                        1 = 1st my cell = calendar year (ICY) 
C 
C     JG        1-20   model year group pointer = IG or 21 - IG, depending on 
C                      whether looking for match or unused cell, respectively, 
C                      in MYGERB.  Also used (as is KG) to distinguish (from IG) 
C                      myg dimension size differences in BD comments. 
C 
C     JPGD      1-5    pointer to base tamper rate to be used in 
C                      computing CSIZE(1-3,IAY).  Depends on IPG. 
C 
C     JSTD      1-6    emissions standard model year group for gas fueled IV. 
C                      The group key is the same as ISTD, but JSTD is set 
C                      internally, instead of being passed in via EVMAIN. 
C                      (The ISTD definition is in the Parameter Dictionary.) 
C 
C     JV        1-2    vapor pressure calculation loop & vectors index (CALUDI) 
C     K         1-24   (reinit) 
C     K1        1-2    (reinit)  
C     KEYEQ     1-3,   (iniscn) 
C 
C     KJDX      KAGE-1-KAGE+1(EVWTENH)  
C 
C 
C     KSTD      1-2    mapping of JSTD used to index R2ALT: 
C                          KSTD = 1, unless JSTD = 4 -> KSTD = 2 
C 
C     KDX       1-21   model year window index = ICY - 1979, unless result is 
C                      < 1 (=> = 1) or > 21 ( => = 21); used by HDDMYM 
C 
C     LAY       1-2    same as IAY, but set for DISCAL call only. 
C 
C     LIM       1-2    = IM, except = 1 when there is an I/M kink, but it is 
C                        after LDXSY & the PREV segment is being evaluated: 
C 
C                      IM   IMKINK   LKINK   LIM   Comments 
C                      --   ------   -----   ---   -------------------------- 
C                       1      1       1      1    P & S: no I/M 
C                       2      1       1      2    P & S: all I/M 
C                       2      2       1      1    P: no I/M; S: kinked I/M 
C                       2      2       2      2    P: kinked I/M; S: all I/M 
C 
C     LIMIT     1-2    index for range check on a variable's value assignment: 
C                      1 = lower bound    2 = upper bound 
C 
C     LINES     1-9    number of lines in user entered data set title 
C 
C     LJDX      25-1   model year relative to LAPSY = JDX - LDXSY + 1; 
C                      used with CUMMIL to get mileage of "previous" segment 
C 
C     LKINK     1-2    same as IMKINK, except IMDXSY must be within (JDX,LDXSY) 
C 
C     M         1-MO   (IMPROG_NP) 
C 
C     MY        NEWMYF-NEWMYL(Hdcgpm) 
C 
C     MYR       1-25   (Getnvt) 
C 
C     NCKIAY    1-3    index into ECK = ICOL index 
C 
C     NDATA     1-6    (Getipr)  
C 
C  5) Subroutine / function parameter dictionary 
C 
C     Variables /arrays occurring in 1 or more of the subprogram argument 
C     lists are defined in this section, unless the variable is a subscript, 
C     in which case look it up in 4)'s subscript dictionary.  The motivation 
C     for defining all the parameter lists in one place is to have only one 
C     instance of a definition in the code plus to encourage only one use 
C     for a mnemonic in the code plus the combined lists contain only 18 
C     names after subscripts have been removed.  Whether scanning on-line or 
C     off a listing, the lookup will be quick. 
C 
C     Parameter array subscripts 
C     -------------------------- 
C     CSFUSE(7)         -  CSFUSE ( IC7 ) 
C     CSF1ST(7,2,3)     -  CSF1ST ( IC7, IPG, IGCSF ) 
C     S2LEFT(4)         -  S2LEFT ( IVTAM ) 
C 
C      Name   Type              Description 
C     ------  ----  ------------------------------------------------------- 
C     CSFUSE   R    category size factors used * ATR to calculate CSIZE 
C     CSF1ST   R    category size factor for 1st calculation of CSIZE 
C     FLGERR   I    number of flag values found to be out of range 
C     ICASE    I    flags the my groups array to be searched. 
C     ICY      I    calendar year for which emission factors are calculated 
C     INERR    I    cumulative count of data entry errors 
C     IOUOUT   I    unit number on which title is printed 
C     IVALUE   I    entered integer value causing error / warning message 
C     IY       I    year to test and adjust to 4-digit year 
C     IY1      I    lower 4-digit bound for IY 
C     IY2      I    upper 4-digit bound for IY 
C     KEYER    I    key to which intercept and slope arrays pair (triple, if 
C                   LDGV 81+ HC or CO) to use in the ER equation calculation 
C                   in BEF's current IP x IV x IR case: 
C                   1 = (ERBZML, ERBDR & ERB50K)  2 = (ERUZML, ERUDR & ERU50K) 
C     MESSAG   I    message code sent by calling subprogram 
C     MY       I    model year 
C     NERR     I    number of error if IY is out of bounds 
C     NREC     I    number of records to be skipped on next READ 
C     RVALUE   R    entered real value causing error / warning message 
C     S2LEFT   R    See local dictionary definition in REFUEL header. 
C     VMLDGT   R    sum of vmt for LDGT1 and LDGT2 (= X23 in MOBILE2) 
C     VMTAGE   R    vehicle cumulative mileage 
C 
C  Evap Study II Parameter List Dictionary for CCEVRT block of routines: 
C 
C      Name    Range               Description 
C     ------  -------  ---------------------------------------------------- 
C     EVRATE     R     dependent variable value from evaporative component 
C                      equation(s) calculation: 
C                        KEYEVP     EVRATE 
C                          1    => hot soak rate 
C                          2    => diurnal rate 
C 
C     ISTD      1-7    emissions standard model year group for gas fueled IV: 
C                        IVGAS    value    myg 
C                        -----    -----    -------------------------------- 
C     (the name of this parameter has changed) 
C                          1        1      pre-1971 
C                                   2      1971 
C                                   3      1972-1977 
C                                   4      1978-1980 
C                                   5      1981+ carbureted (KEYEQ=1 or 2) 
C                                   6      1981+ TBI KEYEQ=1 or 2 
C                                   7      1981+ PFI KEYEQ=1 or 2 
C                          2        1      pre-1971 
C                                   2      1971 
C                                   3      1972-1977 
C                                   4      1978-1980 
C                                   5      1981+ carbureted (KEYEQ=1 or 2) 
C                                   6      1981+ TBI KEYEQ=1 or 2 
C                                   7      1981+ PFI KEYEQ=1 or 2 
C                          3        1      pre-1979 
C                                   2      1979-1980 
C                                   3      1981+ carbureted (KEYEQ=1 or 2) 
C                                   4      1981+ TBI KEYEQ=1 or 2 
C                                   5      1981+ PFI KEYEQ=1 or 2 
C                          4        1      pre-1985 
C                                   2      1985+ (KEYEQ=1 or 2) 
C                          5        1      pre-1978 
C                                   2      1978-1979 
C                                   3      1980-81 
C                                   4      1982-84 
C                                   5      1985+ 
C 
C     ITRSTL    1-2    technology type for resting losses calculation 
C                        1=fuel injected with open bottom canister 
C                        2=fuel injected with closed canister 
C                       (currently carbureted=fuel injected, open) 
C 
C     IVGAS     1-15   gasoline fueled vehicle types : 
C                       1 = LDGV    2 = LDGT1   3 = LDGT2  4 = LDGT3  5 = LDGT4 
C                       6 = HDGV2B  7 = HDGV3   8 = HDGV4  9 = HDGV5 10 = HDGV6 
C                      11 = HDGV7  12 = HDGV8A 13 = HDGV8 14 = MC    15 = GAS     C                           BUS 
C                      IVGAS = IV, except MC and GAS BUS 
C 
C     IVTL      1-4    running loss HC: vehicle classes with rl HC emissions 
C 
C     KEYEQ     1-3    HS/DU equation group selection key: 
C                        1 = passing 
C                        2 = failed purge 
C                        3 = failed pressure 
C 
C     KEYEVP    1-2    Evaporative component calculation selection key 
C                        1 = hot soak  2 = diurnal 
C 
C     RIVAL      R     independent variable value in evaporative component 
C                      equation(s) calculation: 
C                        KEYEVP                  RIVAL 
C                          1    => hot soak => fuel RVP 
C                          2    => diurnal  => diurnal Index 
C 
C  *************************************************************************** 
C 
C     IC, IDX, IFDS, IG, IM, IP, IP1, IP2, IR, IV and JSTD are parameters 
C     defined in the subscript dictionary. 
C 
C  Notes: 
C 
C  **** NOT IN USE **** = the subroutine/function is not called by the current 
C  version of MOBILE4.  Release MOBILE4 (final version), for example, does not 
C  use WTLRNL or OUTTAM. 
C 
C 
