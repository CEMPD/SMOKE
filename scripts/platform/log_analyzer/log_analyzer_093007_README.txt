Usage: log_analyzer.py [options] <log_files>

Options:
 -h, --help            show this help message and exit
 -k KNOWNMSGFILE       known message file name
 -f OUTFILE            output file to write summary data
 -d INDIR              logs directory
 -l LEVEL              level of analysis: 1, 2, 3.  Default is level 1.
 --sort=SORTTYPE       sort output by column: priority, known, count,
                       message, file.  Default is priority.
 --delimiter=DELIMITER
                       delimiter to separate columns.  Default is ','
 --list_unknowns       List all unknown messages.  If not set, will sum all
                       the unknown messages and list as 1 line in output
 --unknown_priority=UNKNOWNPRIORITY
                       Priority of unknown messages. Default is 0
 --warnings_only       Search only for warning messages.  Default is False,
                       search for both warnings and errors.
 --errors_only         Search only for error messages.  Default is False,
                       search for both warnings and errors.
 -e EXITPRIORITY       If priority is less than or equal to the exit
                       priority, then log analyzer will exit with a non zero
                       exit status.  If this isn't set, default, it will
                       return exit status equal to 0, independent of what
                       type of messages are found




Some examples:

1.  Level 3 analysis searching for all files in a directory. Note: a
directory search does an internal 'find' command, which allows this to
deal w/ very large # of log files

$ log_analyzer.py -k known_messages.list  -l 3  -f ~/tmp/test_log5.csv
 -d ~/emf/smoke_optim/logs

2.  Level 2 analysis searching for specific files and sorting based on
   message count.  Include a separate line for each unknown message:

$ log_analyzer.py -k known_messages.list  -l 2 --sort count -f ~/tmp/test_log4.csv ~/smoke_emf/2002/smoke/intermed/v3/2002ac/alm/logs/*.log --list_unknowns
 3.  Level 1 analysis, return nonzero exit if find messages of priority
   less than or equal to 1:

$ log_analyzer.py -k known_messages.list  -e 1 -f ~/tmp/test_log3.csv ~/tmp/logs/*.log


Known message types:   I included my working version of known messages file.  The
 priorities are probably not right.  Notice, that some special
 characters need to be "escaped" to prevent the python regular
 expression analyzer from mis-interpreting them.  The key wild card I
 used was ".*" . I found that I it was easier not to try and preserve
 the number of spaces around a variable or number.  See the python
 regular expression syntax if you want to be more specific in your
 searches:  http://www.python.org/doc/current/lib/re-syntax.html

Some example outputs:
 test_analysis_lev1.csv -- level 1, sort by priority

 test_analysis_lev1_count.csv -- lev 1, sort by message count

 test_analysis_lev1_missing.csv -- lev 1, summed missing

 test_analysis_lev1_listmissing.csv -- lev 1, list all missing

 test_analysis_lev2.csv -- lev 2, sort by count

 test_analysis_lev3.csv  -- lev 3, sort by priority, 1 log file



Some timing stats on running the log analyzer:


1.  On a "typical" individual log file
lev 1 ~ .5 sec
lev 2 ~ .5 sec
lev 3 ~ .5 sec


2.  Directory with 74 log files in it (my alm dir)
lev 1 ~ 1.7 sec
lev 2 ~ 1.8
lev 3 ~ 1.6

3.  EPA's 2020 logs directory (all logs minus 12 largest -
   smkinven_ptipm_<mon>_2020.log) - 13,604 log files
   Note: I separated the internal grep from the classification (internal grep
   ~11.5 minuts)
additional time over the grep:
lev 1 -  30 sec
lev 2 - 600 sec
lev 3 -  30 sec
