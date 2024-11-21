#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


###################################################################
try:
    from builtins import map
except ImportError:
    from __builtin__ import map
try:
    from builtins import filter
except ImportError:
    from __builtin__ import filter

"""
Log analyzer for SMOKE logs.

Parses and characterizes warnings and error messages from the SMOKE
logs.


Inital version: Alexis Zubrow, IE UNC-Chapel Hill, September 2007 
28 Mar 2013: Moves from os.popen to subprocess.Popen 
"""
###################################################################

import re, subprocess
import string, os


### Some functions used for sorting keys

def getPriority(msgLog):
    """
    get the priority from a msgLog
    """
    return msgLog.priority

def getFileName(msgLog):
    """
    get the fileName from a msgLog
    """
    return msgLog.fileName

def getMsgTypeStr(msgLog):
    """
    get the message type string from a msgLog
    """
    return msgLog.msgTypeStr

def getMsgTypeInd(msgLog):
    """
    get the message type index from a msgLog
    """
    return msgLog.msgTypeInd

def getMsg(msgLog):
    """
    get the message from a msgLog
    """
    return msgLog.msg

## Example of how to use in sort (python 2.4 on):
## msgLogLst.sort(key=getMsgTypeInd)

## series of functions to test equality
## uses a non-local variable, bad idea but
## helps w/ maps and filters
def testInd(msg):
    return msg.msgTypeInd == cmpVal

def testFile(msg):
    return msg.fileName == cmpVal

## test not a new line
def notNewLine(a):
    return a != "\n"

## convert a series of elements to list
def lineLst2strLst(aLst):
    return ["%s" %elem for elem in aLst]

## convert a message to level 3 line list
def msg2lev3(msg):
    lineLst = [msg.priority, "'" + msg.msgTypeStr + "'", msg.line, \
               "'" + msg.msg + "'", "%s\n" %msg.fileName]
    return lineLst



#####################
## Class definitions
#####################


class MessageLog(object):
    """
    Message log class.  Organizes all the info for a particular
    message.  A message is one warning or error line from the SMOKE
    logs. Contains information on message, message type, line number
    and file.
    """

    def __init__(self, msgLine):
        """
        Input:

          msgLine  - one warning or error line.

        Note:
          1.  Format expected for msgLine is 'filename: line #: message'
          
        """
        self.msg = None
        self.fileName = None
        self.line = None
        self.msgTypeStr = None
        self.msgTypeInd = None
        self.known = False
        self.priority = None

        ## parse the message line
        self.parseLine(msgLine)

    def parseLine(self, msgLine):
        """
        Input:

          msgLine  - one warning or error line.

        Note:
          1.  Format expected for msgLine is 'filename: line #: message'
          
        """

        ## split the line based on ":"
        ## element 0 = file, element 1 = line #, element 2: = message
        tmpLst = msgLine.split(":")
        self.fileName = tmpLst[0]
        self.line = int(tmpLst[1])
        
        ## for the message [2:], recreate w/ the ':' and left strip
        ## any spaces and remove any newlines
        self.msg = ':'.join(tmpLst[2:]).strip()


    def checkType(self, msgType):
        """
        Checks if message matches the message type.  If the msgType is
        known, it does a regular expression search.  If the message is
        unkown it tests if the message and message type strings are
        equal.

        Inputs:
          msgType - message type object
          
        
        Returns
           True - matches
           False - does not match
           
        """

        matchObj = None

        ## for known message types, use regular expression for search
        if msgType.known:
            matchObj = msgType.typeRegexp.search(self.msg)
            if matchObj is not None:
                ## matched
                self.known = True
                self.priority = msgType.priority
                self.msgTypeStr = msgType.msgTypeStr
                return True

            ## did not match
            return False

        else:
            ## unknown message, therefore test w/ string equivalence
            if self.msg == msgType.msgTypeStr:
                ## matched
                self.known = False
                self.priority = msgType.priority
                self.msgTypeStr = msgType.msgTypeStr
                return True

            ## did not match
            return False

    def selectType(self, msgTypeLst, lastInd = 0):
        """
        This selects the matching message type from the list of
        message types.  If no message type matches the message, it
        returns None

        Inputs:
           msgTypeLst - list of current message types (known and unknown)

           lastInd - the index in the msgTypeLst that the previous message matched.

        Returns:
           ind - index of matching msgType in msgTypeLst. None if no
                 message types match this message
           
        Note:

          1.  The last index is used to speed up the search through
              msgTypeLst.  The idea is that messages often come in
              groups.  If that last message checked was of this
              message type, than likely this message is of this
              message type.  In other words, begin search at the last
              message type index and search through the end of list.
              If there still isn't a match, start search from
              beginning of msgTypeLst and search through last index.

          2.  There is probably a more efficient way to do this.
          
        """

        cnt = 0
        ind = 0
        
        ## start the search at the last index for a matching message
        for i,msgType in enumerate(msgTypeLst[lastInd:]):
            cnt += 1
            if self.checkType(msgType):
                ind = i + lastInd
                self.msgTypeInd = ind
                return ind

        ## Didn't find the message from last index on, start search
        ## from beginning of list through last index
        for i,msgType in enumerate(msgTypeLst[:lastInd]):
            cnt += 1
            if self.checkType(msgType):
                ind = i
                self.msgTypeInd = ind
                return ind

        ## Doesn't match any of the known message types
        return None
    
##         ## message doesn't match any of current message types create a
##         ## new unknown message type and append it to message type
##         ## list. Use checkType to update the data
##         newMsgType = MessageType(self.msg, False)
##         if not self.checkType(newMsgType): raise LookupError("Message does not match self"
##         msgTypeLst.append(newMsgType)
##         ind = len(msgTypeLst) - 1
##         self.msgTypeInd = ind 
##         print("Search lastInd = %d, index = %d, search cnt = %d, unknown"  %(lastInd, ind, cnt)
##         return ind

            
            
########################################


class MessageType(object):
    """
    Message type class.  A message type is a generalized warning or
    error message.  It consists of message type string and a regular
    expression.  A message type is used to group similar messages
    together from the SMOKE logs.

    There are two ways that message types are treated in the
    analysis. Known message types are used as regular expressions to
    match from the beginning of the message.  Unknown message types
    are used as strings to test equivalency with the message.
    """

    def __init__(self, typeLine, known = True, priority = 1):
        """
        Input:

          typeLine  - one message type line.

          known - is it a known message, ie. from the known list file

          priority - the priority of the message (1 highest)

        Note:
          1.  Format expected for typeLine is "message type string", priority

          2.  If known is false, the message type string will be
              recorded w.o converting to a regular expression, with
              the priority from the parameter.  If known is true, the
              priority is taken from the typeLine and priority
              parameter is ignored.
          
        """

        self.msgTypeStr = None
        self.typeRegexp = None
        self.priority = None
        self.known = known
        
        if known:
            ## parse the type line
            self.parseTypeLine(typeLine)

        else:
            ## unknown message so record the whole line as the
            ## msg type string
            self.priority = priority
            self.msgTypeStr = typeLine

    def parseTypeLine(self, typeLine):
        """
        Input:

          typeLine  - one message type line.


        Note:
        
          1.  Format expected for typeLine is 'priority, "message type
              string"'.  This is comma separated.

          2.  This creates a regular expression from the message type
              string.  It is used as a search from the beginning of
              the string.
          
        """


        ## separate the priority from the message type string

        ## split line on ","
        ## element 0 = priority, element [1:] is message type
        tmpLst = typeLine.split(",")
        self.priority = int(tmpLst[0])

        ## recreate message type [1:] by readding w/ the ',' and
        ## left strip any spaces and right strip any spaces or new lines
        self.msgTypeStr = ','.join(tmpLst[1:]).strip()

        ## create a regular expression for this message type
        self.createRegExp()
            
            

    def createRegExp(self):
        """
        Creates a regular expression for this message type string.
        Assumes that the user has added the appropriate wild cards
        into the message type string.
        """

        ## compile the regular expression so the message type string
        ## must match the beginning of the message
        self.typeRegexp = re.compile("^" + self.msgTypeStr, re.MULTILINE)

## --------------------
def pipethrough(cmd, lineLst=[]):

    """
        Allows one to send a command to the shell and get back the
        stdout.  This is a very general program w/ applications
        outside of this projects.

        Based on 'Python in a Nutshell' by Alex Martelli, pg 290.

        input:
          cmd -   a command to run, simply a string.

          lineLst - not sure??? I leave it blank and it seems to work.

        Output:
          resultLst -  the stdout and stderror.  Each line is the element of the
                       list.

        Example:
            buff = pipethrough('ls -l *.py')

            buff will contain the result of:
               $ ls -l *.py
            With each line of stdout as a seperate element of the
            buff list.

    """
    ## get standard in and standard out/error
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True, close_fds=True)
    p.stdin.writelines(lineLst)
    p.stdin.close()
    resultLst = p.stdout.readlines()
    resultLst = [x.decode("utf-8") for x in resultLst] # converts bytes to strings
    p.stdout.close()
    return resultLst


##-----------------------------------------------------------------------##
## Actual script


if __name__ == "__main__":


    from operator import itemgetter
#    from string import atoi # not in python 3.x
    import sys
    from optparse import OptionParser

    ## exit status if returns an error b.c. of "-e" option
    errorStatus = 10

    ## Command line options
    usage = "usage: %prog [options] <log_files>"
    parser = OptionParser(usage=usage)

    parser.add_option("-k", dest="knownMsgFile", default=None,
                      help="known message file name")

    parser.add_option("-f", dest="outFile", default=None,
                      help="output file to write summary data")
    parser.add_option("-d", dest="inDir", default=None,
                      help="logs directory")
    parser.add_option("-l", dest="level", default=1, type="int",
                      help="level of analysis: 1, 2, 3.  Default is level 1.")
    parser.add_option("--sort", dest="sortType", default="priority",
                      help="sort output by column: priority, known, count, "
                      "message, file.  Default is priority.")
    parser.add_option("--delimiter", dest="delimiter", default=",",
                      help="delimiter to separate columns.  Default is ','")
    parser.add_option("--list_unknowns", dest="listUnknown", default=False,
                      action="store_true",
                      help="List all unknown messages.  If not set, will sum "
                      "all the unknown messages and list as 1 line in output")
    parser.add_option("--unknown_priority", dest="unknownPriority", default=0,
                      type="int", help="Priority of unknown messages. "
                      "Default is 0")
    parser.add_option("--warnings_only", dest="onlyWarnings", default=False,
                      action="store_true",
                      help="Search only for warning messages.  Default "
                      "is False, search for both warnings and errors.")
    parser.add_option("--errors_only", dest="onlyErrors", default=False,
                      action="store_true",
                      help="Search only for error messages.  Default "
                      "is False, search for both warnings and errors.")
    parser.add_option("-e", dest="exitPriority", default=None,
                      type="int", help="If priority is less than or equal to "
                      "the exit priority, then log analyzer will exit with "
                      "a non zero exit status.  If this isn't set, default, it will "
                      "return exit status equal to 0, independent of what type "
                      "of messages are found")

    
    (options, args) = parser.parse_args()

    logFilesLst = None
    if len(args) > 0:
        logFilesLst = args
        ## expand user
        logFilesLst = list(map(os.path.expanduser, logFilesLst))
    
    ## get options
    knownMsgFile = options.knownMsgFile
    outFile = options.outFile
    inDir = options.inDir
    level = options.level
    sortType = options.sortType
    delimiter = options.delimiter
    listUnknown = options.listUnknown
    unknownPriority = options.unknownPriority
    onlyWarnings = options.onlyWarnings
    onlyErrors = options.onlyErrors
    exitPriority = options.exitPriority

    ## Some info to stdout
    print("log analyzer")
    ## if log files list and inDir are None, then raise an error
    if (logFilesLst is None) and (inDir is None):
        raise LookupError("Must either define an input logs directory or pass log files as arguments")
        
    print("Getting message data (might take some time).... ")

    ## if using a logs directory, do a find exec grep -- prevents too
    ## many files
    if inDir is not None:
        ## check if real directory
        inDir = os.path.expanduser(inDir)
        if not os.path.isdir(inDir):
            raise IOError("The logs directory was not found: %s" %inDir)

        ## find exec grep to avoid too many files
        progStrWarning = "find %s -exec grep -Hani warning {} \; " %inDir
        progStrError = "find %s -exec grep -Hani error {} \; " %inDir
    else:
        ## check if files exist
        logFilesExist = list(filter(os.path.isfile,logFilesLst))
        if (len(logFilesExist) < len(logFilesLst)):
            ## some files don't exist
            raise IOError("One or more of the input logs doesn't exist")

        ## make all files into one long string and define a simple grep
        logFiles = ''.join(logFilesLst)  
        progStrWarning = "grep -Hani warning %s" %logFiles
        progStrError = "grep -Hani error %s" %logFiles

    ## grep the messages
    if onlyWarnings:
        messagesGrep = pipethrough(progStrWarning)
    elif onlyErrors:
        messagesGrep = pipethrough(progStrError)
    else:
        ## grep for both warnings and errors
        messagesW = pipethrough(progStrWarning)
        messagesE = pipethrough(progStrError)
        messagesGrep = messagesW + messagesE
    if len(messagesGrep) == 0:
        print("No warning or error messages. Exiting")
        sys.exit(0)

    ## Use the known messages file to get messagetypes
    if knownMsgFile is not None:
        g = open(knownMsgFile, "rt")  ## open file converting new lines to unix type

        ## Read the known message list into a buffer
        ## strip out lines that are just newline
        buff = g.readlines()    
        buffKnLst = list(filter(notNewLine, buff))
        g.close()

        ## get a list of message types
        knTypLst = list(map(MessageType, buffKnLst))
    else:
        ## no known messages defined
        knTypLst = []
        
    ## Get a list of messageLog 
    msgLst = list(map(MessageLog, messagesGrep))
    print("Finished getting data")

    ## store messages into a known and unknown list
    knMsgLst = []
    unknMsgLst = []

    lastInd = 0

    ## loop over each message and type match it, if it is not found
    ## append to the unknown list.  Use that index to speed up the
    ## search
    print("Classifying message types...")
    fileLst = []
    for msg in msgLst:
        newInd = msg.selectType(knTypLst, lastInd)
        if newInd is None:
        ## unkown message
            msg.priority = unknownPriority
            unknMsgLst.append(msg)
        else:
        ## known message, store and record the index for next pass
            knMsgLst.append(msg)
            lastInd = newInd

    ## final list of summarised data to print(out
    outLst = []

    ## Some diagnostics:
    print("Total number of known messages: %d" %len(knMsgLst))
    print("Total number of unknown messages: %d" %len(unknMsgLst))
    
    if (level == 1):
        print("Level 1 analysis...")
        ## loop over message type
        for ind, msgType in enumerate(knTypLst):

            ## filter the known message list, use non-local variable (cmpVal)
            ## to define the index to test against
            cmpVal = ind
            subMsgLst = list(filter(testInd, knMsgLst))

            numMsg = len(subMsgLst)
            if (numMsg > 0):
            ## there are messages of this type,
                ## write a new list to output list
                ## append a newline character to the end of filename,
                lineLst = [msgType.priority, msgType.known, numMsg,\
                            "'" + msgType.msgTypeStr + "'",\
                           "%s\n" %subMsgLst[0].fileName]
                outLst.append(lineLst)

    elif (level == 2):
        print("Level 2 analysis...")
        ## filter by msg type and file -- slow
        
        ## loop over message type, same as level 1
        for ind, msgType in enumerate(knTypLst):

            ## filter the known message list, use non-local variable (cmpVal)
            ## to define the index to test against
            cmpVal = ind
            subMsgLst = list(filter(testInd, knMsgLst))

            numMsg = len(subMsgLst)
            if (numMsg > 0):
            ## there are messages of this type,
                ## loop over unique files 
                fileLst = list(map(getFileName, subMsgLst))
                uniqueFileLst = list(set(fileLst))

                ## loop over each file and create a subsubLst of msgtype by file
                for uniqueFile in uniqueFileLst:
                    ## filter the known message list, use non-local variable (cmpVal)
                    ## to define the file name to test against
                    cmpVal = uniqueFile
                    subsubMsgLst = list(filter(testFile, subMsgLst))
                    numMsg = len(subsubMsgLst)

                    ## there are messages in this file by definition
                    ## write a new list to output list
                    ## append a newline character to the end of filename,
                    lineLst = [msgType.priority, msgType.known, numMsg,\
                                "'" + msgType.msgTypeStr + "'",\
                               "%s\n" %uniqueFile]
                    outLst.append(lineLst)


    else:
        ## level 3 (do later)
        pass

    
    numUnknMsg = len(unknMsgLst)
    if (numUnknMsg > 0):
        if listUnknown:
            ## add all the unknown messages to the list
            ## use the message as message type
            for msg in unknMsgLst:
                lineLst = [msg.priority, msg.known, 1, "'" + msg.msg + "'", \
                           "%s\n" %msg.fileName]
                outLst.append(lineLst)
        else:
            ## just sum up # of messages
            firstFile = unknMsgLst[0].fileName
            lineLst = [unknownPriority, False, numUnknMsg , \
                       "all unknown message types", "%s\n" %firstFile]
            outLst.append(lineLst)

    ## sort outLst based on sort type
    if (sortType == "priority"):
        item = 0
        reverse = False
    elif (sortType == "known"):
        item = 1
        reverse = False
    elif (sortType == "count"):
        item = 2
        reverse = True  ## sort in reverse order (largest # first)
    elif (sortType == "message"):
        item = 3
        reverse=False
    elif (sortType == "file"):
        item = 4
        reverse=False
    else:
        raise LookupError("Unknown sort type: %s" %sortType)

    if level == 3:
        ## level 3
        ## List all messages, no unknown messages included
        print("Level 3 analysis...")

        ## for each line get a list: priority, msgType, line, msg, file
        outLst = list(map(msg2lev3,knMsgLst))
        
        ## level 3 sort outLst based on sort type
        if (sortType == "priority"):
            item = 0
            reverse = False
        elif (sortType == "known"):
            ## all are known so just sort by priority
            item = 0
            reverse = False
        elif (sortType == "count"):
            ## no count so just sort by priority
            item = 0
            reverse = False
        elif (sortType == "message"):
            item = 1
            reverse=False
        elif (sortType == "file"):
            item = 4
            reverse=False
        else:
            raise LookupError("Unknown sort type: %s" %sortType)


    ## sort by that item # of the sub list
    outLst.sort(key=itemgetter(item), reverse=reverse)

    ## convert each elem of outLst to list of strings
    outLst = list(map(lineLst2strLst, outLst))
    ## insert header in front of sorted list
    ## level 1 analysis
    if (level == 1):
        ## output of the program, headers of the file
        headerLst = ["priority", "known", "count", "message type", "first file\n"]
    elif (level == 2):
        headerLst = ["priority", "known", "count", "message type", "file\n"]
    else:
        ## level 3
        headerLst = ["priority", "message type", "line", "full message", "file\n"]

    ## Take each line list and join into a string using the delimiter
    outStrLst = [delimiter.join(elem) for elem in outLst]
    headerStr = delimiter.join(headerLst)

    print("Finished classifying message types")

    if outFile is not None:
        ## write to out file
        k = open(os.path.expanduser(options.outFile),"w")
        k.write(headerStr)
        k.writelines(outStrLst)
        k.close()
    else:
        ## print to stdout
        print("")
        print(headerStr.rstrip())
        for lineStr in outStrLst:
            print(lineStr.rstrip())  ## get rid of new line char

    ## test priority for non-zero exit
    if exitPriority is not None:
        print("Testing for exit priority <= ", exitPriority)
        ## resort the outLst based on priority
        outLst.sort(key=itemgetter(0))

        ## check if the first element is less than or equal to the exit
        ## exit priority.  IF so, it is an error
        if int(outLst[0][0]) <= exitPriority:
##             raise LookupError(\
##                   "At least one message's priority is less than or equal to the exit priority: %d" \
##                   %exitPriority
            print("ERROR: At least one message's priority is less than or equal to the exit priority: %d" %exitPriority)
            sys.exit(errorStatus)
        else:
            print("All message priorities > ", exitPriority)
            sys.exit(0)
    else:
        sys.exit(0)
