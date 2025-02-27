#!/work/EMIS/python/miniforge3/envs/smoke_python/bin/python
# The "new2" program for approach 2.  
# Creates csv annual sums from temperature adjusted and unadjusted onroad netCDFs.
# 10/21/08 James Beidler <beidler.james@epa.gov>
# Updated 4/16/09
# Update 7/13/09 by C. Allen - changed definition of dayPath so that this script works for different years
# Updated 12/23/09 by J. Beidler - Fixed input paths to work with unadjusted summaries.
# Updated 1/8/10 by J. Beidler - Changed parameter file loading.  Added calendar based dates.
# Update 25 Mar 2013 by J. Beidler - code overhaul for using new libraries and consolidating into one script
# Update 11 June 2013 by J. Beidler - Modified to work without SCCs, giving just county/state reports for faster processing

from __future__ import division
from __future__ import print_function
from past.builtins import execfile
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import numpy as np
from scipy.io.netcdf import *
import sys, os, csv, calendar 

def checkEV(evName):
	"""
	Checks if an environment variable is set.  If not, exits.  If it is, returns the variable.
	Takes the name of the environment variable.
	"""
	try: 
		var = os.environ[evName]
	except:
		print("ERROR: Environment variable '%s' is not defined." %evName)
		sys.exit(1)
	else: 
		return var

baseYear = checkEV('BASE_YEAR') # Modeling year, used in definition of dayPath
scriptdir = checkEV('SCRIPTS')
parameter_file = '%s/annual_report/parameter_file_cmaq_cb05.txt' %scriptdir

# Path and filename of molecular weight conversion files
mwDict = {'cmaq_cb05_soa': parameter_file, 'cmaq_cb05': parameter_file, 'cmaq_cb05_tx': parameter_file, 'saprc07t': '/garnet/work/bte/pyQA/parameter_file_saprc07t.txt'} 
mwDefault = 'cmaq_cb05'  #  Default molecular weight file 
# End User Config

def openFile(fileName, accessType = 'r'):
	"""
	Tests to see if a file is available for access.  If not it returns an error and exits.
	If it is available it returns an open file object.
	"""
	try: file = open(fileName, accessType)
	except:
		print("ERROR: %s not available for access." %fileName)
		sys.exit(1)
	else: return file

def fixLen(x, y = 2):
	x = str(x)
	while len(x) < y:
		x = '0' + x
	return x 

class createReport(object):
	"""
	"""

	def __init__(self, Mon = ""):
		self.Mon = Mon

	def nameInFile(self, inDate):
		"""
		Set the infile name based on the SMOKE conventions.
		""" 
		inFileName = 'emis_mole_' + sector + '_adj_%s_' %inDate + grid + '_' + spec + '_' + case + '.ncf'
#		inPath = os.path.join(imPath, sector + '_adj')
		inPath = '%s_adj' %premerged
		return os.path.join(inPath, inFileName)

	def moleToMassFactor(self, speciesName):
		"""
		Get the moles to mass conversion factor for the specified species.
		"""
		if speciesName == 'NAPHTH_72': speciesName = 'NAPHTHALENE'
		if '_72' in speciesName:
			tmp = speciesName.split('_72')
			speciesName = tmp[0]
		if speciesName in molecDct: 
			factor = molecDct[speciesName]
		else: 
			factor = 1
		return factor

	def sumSpeciesDay(self, species):
		"""
		Sum a species over a day
		"""
		dataIn = species[:]
		speciesDaySum = np.zeros([species.shape[2], species.shape[3]], '>f')
		for hour in range(species.shape[0] - 1):
			speciesDaySum += dataIn[hour][0][:][:]
		return speciesDaySum

	def convValue(self, monSum, speciesName):
		"""
		Convert monthly value from moles per second to tons per hour.
		"""
		monSum = monSum * self.moleToMassFactor(speciesName)   # Convert moles per second to grams per second
		monSum = monSum * (0.00000110231131)  # Convert grams per second to tons per second
		monSum = monSum * 3600   # Convert tons per second to tons per hour
		return monSum  

	def __apportion(self, specSum, specDict):
		"""
		Apportion the cell values to FIPS based on cell-county inven report ratios
		"""

		outDict = {}
		for fips in list(specDict.keys()):
			outDict[fips] = 0

			for cell in list(specDict[fips].keys()):
				col = int(cell.split(',')[0]) - 1
				row = int(cell.split(',')[1]) - 1
				outDict[fips] += specSum[row,col] * specDict[fips][cell]
		return outDict

	def sumSpeciesMonth(self, facs):
		"""
		Sum a species over a month.
		"""
		mon = int(self.Mon) 
		inDate = Year+self.Mon+'01'

		factorDict = facs.countyFac

		# Open the infile for reading
		inFileName = self.nameInFile(inDate)
		try:
			inFile = netcdf_file(inFileName, 'r')
		except:
			sys.exit('ERROR: Could not open %s' %inFileName)

		# Get the list of variables from the in file
		variableNames = list(inFile.variables.keys())
		vNames = [species for species in variableNames if species != 'TFLAG']
		speciesMonSum = dict(list(zip(vNames, [{} for x in range(len(vNames))])))

		cpl = 0
		cMax = float(len(vNames))
		sys.stdout.write("\r%.2f%%" %(old_div(float(cpl),cMax)))
		for speciesName in vNames:
			for day in range(1, calendar.monthrange(int(Year), mon)[1] + 1):
				Day = fixLen(day)

				# Set in date and out date and open the in and out files.
				inDate = Year+self.Mon+Day

				# Open the infile for reading
				inFileName = self.nameInFile(inDate)
#				print("In File: " + inFileName)
				try:
					inFile = netcdf_file(inFileName, 'r')
				except:
					sys.exit('ERROR: Could not open %s' %inFileName)


				species = inFile.variables[speciesName]
				# Sum a species over a day
				speciesDaySum = self.sumSpeciesDay(species)

				# Calculate and populate monthly sums
				if day == 1:
					speciesMonSum[speciesName] = speciesDaySum
				else:
					speciesMonSum[speciesName] += speciesDaySum

			try:
				speciesMonSum[speciesName] = self.convValue(speciesMonSum[speciesName], speciesName)
			except:
				sys.exit('%s' %list(speciesMonSum[speciesName].keys()))
			speciesMonSum[speciesName] = self.__apportion(speciesMonSum[speciesName], factorDict[speciesName])

			cpl += 1
			sys.stdout.flush()
			sys.stdout.write("\r%.2f%%" %((old_div(float(cpl),cMax))*100.))

		sys.stdout.write("\n")
		return speciesMonSum

class repFac(object):
	"""
	Create a class for output dictionary
	"""
	def __init__(self):
		self.countyDescDict = {}
		self.stateDescDict = {}
		self.speciesList = []
		self.repFileName = self.__nameRepFile()
		self.repDict = self.readRepFile()
		self.countyFac = self.calcCountyFac()

	def __nameRepFile(self):
		"""
		Set the infile name for the REP file based on the SMOKE conventions.
		Christine Allen updated 6 Dec 2018 for othptdust, which is a monthly sector, and as such there is no annual report.
		  Since this report is only used for spatial allocation, use September report to represent an "average" month.
		Updated again 31 Mar 2023 to cover new sector names canada_afdust and canada_ptdust.
		22 Dec 2023 Christine update: Canada_afdust is monthly only in 2020.
		"""
		if sector == 'othptdust' or sector == 'canada_ptdust' or (sector == 'canada_afdust' and baseYear == '2020'):
			inFileName = 'rep_' + sector + '_sep_' + case + '_invgrid_cell_county_' + grid + '.txt' 
		else:
			inFileName = 'rep_' + sector + '_' + case + '_invgrid_cell_county_' + grid + '.txt' 
		inPath = os.path.join(home, case + '/reports/inv')
		return os.path.join(inPath, inFileName)

	def readRepFile(self):
		"""
		Read the rep file and put into a list.
		"""
		repDict = {}
		try:
			print('Loading %s' %self.repFileName)
			repFile = open(self.repFileName,'r')
		except:
			sys.exit('ERROR: Could not open %s' %self.repFileName)

		for line in csv.reader(repFile, delimiter='|'):
			line = [col.strip() for col in line]
			# On the line that starts with #LABEL...
			if '#Label' in line[0]:
				# Added by C. Allen, 1/30/13, to support Smkreport v3.1
				# If pollutant name contains "I-", drop the "I-"
				for i,speciesName in enumerate(line):
					if 'I-' in speciesName[:3]: 
						line[i] = speciesName[2:]
						print("NOTE: Smkreport v3.1 detected, using %s as %s" %(speciesName, line[i]))
				# Create the specieslist starting at the 8th col 
				self.createSpeciesList(line[7:])
				# Create the reference dictionary of the columns
				colDict = dict(list(zip(line, list(range(len(line))))))

			# If the line is a comment, the skip it 
			if line[0][0] == "#": continue

			col = line[colDict['X cell']]
			row = line[colDict['Y cell']]
			cell = col+','+row
			fips = line[colDict['Region']][1:]
			# scc = line[colDict['SCC Tier 3']] # Commented out because newer reports have full SCCs (C. Allen, 11 Aug 2014)
			state = line[colDict['State']]
			county = line[colDict['County']]

			for speciesName in self.speciesList:	
				value = float(line[colDict[speciesName]])
				if fips not in repDict:
					repDict[fips] = { cell: { speciesName: value } } 
				elif cell not in repDict[fips]: 
					repDict[fips][cell] = { speciesName: value } 
				elif speciesName not in repDict[fips][cell]: 
					repDict[fips][cell][speciesName] = value
				else:
					repDict[fips][cell][speciesName] += value
	 
			if fips not in self.countyDescDict: 
					self.countyDescDict[fips] = {'st': state, 'cty': county}
			if fips[:8] not in self.stateDescDict: 
					self.stateDescDict[fips[:8]] = {'st': state}
		return repDict

	def createSpeciesList(self, speciesLine):
		"""
		Create a list of species from the report file
		"""
		for speciesName in speciesLine:
			if 'S-' in speciesName[:3]: 
				continue
			# Skip PM10 and PM2_5
#			if 'PM10' in speciesName: 
#				continue
#			if 'PM2_5' in speciesName: 
#				continue
			if 'XPMC' in speciesName:
				continue
			else:
				self.speciesList.append(speciesName)

	def calcCountyFac(self):
		"""
		Calculate the county allocation factors.
		"""
		print("Calculating the county factors")

		crSum = {}
		fipsCRSum = {}
		for fips in self.countyDescDict:
			fipsCRSum[fips] = {} 
			for cell in self.repDict[fips]:
				if cell not in crSum: 
					crSum[cell] = {} 
				if cell not in fipsCRSum[fips]: 
					fipsCRSum[fips][cell] = {}
				for speciesName in self.repDict[fips][cell]:
					if speciesName not in crSum[cell]: 
						crSum[cell][speciesName] = 0
					if speciesName not in fipsCRSum[fips][cell]: 
						fipsCRSum[fips][cell][speciesName] = 0 
					crSum[cell][speciesName] += self.repDict[fips][cell][speciesName]
					fipsCRSum[fips][cell][speciesName] += self.repDict[fips][cell][speciesName]
		factorDict = {}
		for fips in self.countyDescDict:
			for cell in crSum:
				if cell in fipsCRSum[fips]:
					for speciesName in self.speciesList:
						if crSum[cell][speciesName] != 0:
							ratio = (old_div(fipsCRSum[fips][cell][speciesName], crSum[cell][speciesName]))
						else:
							ratio = 0 

						if speciesName not in factorDict:
							factorDict[speciesName] = {fips: {cell: ratio}}
						elif fips not in factorDict[speciesName]:
							factorDict[speciesName][fips] = {cell: ratio}
						else:
							factorDict[speciesName][fips][cell] = ratio
		return factorDict

class repOut(object):
	def __init__(self):
		self.outDict = {} 

	def __loadOut(self, areaType, dateType, speciesList):
		"""
		Set the outfile name.
		"""
		outFileName = 'rep_' + areaType + '_' + dateType + '_adj_' + sector + '_' + case + '_' + grid + '.txt'
		byType = 'emissions'
		outPath = os.path.join(home, case + '/reports/programs')
		outFileName = os.path.join(outPath, outFileName)

		if outFileName not in self.outDict:
			outFile = openFile(outFileName, 'w')
			outFile.write( '# %s adjusted %s %s %s for months %s-%s\n' %(sector, areaType, dateType, byType, monList[0], monList[ len(monList) - 1 ]) )
			print('Writing %s %s adjusted %s to: %s' %(areaType, dateType, byType, outFileName))

			if dateType == 'monthly': 
				headerLine = '# Month,FIPS,State Name,'
			else: 
				headerLine = '# FIPS,State Name,'
			if areaType == 'county': 
				headerLine = headerLine + 'County Name,'
			outFile.write(headerLine + ','.join(speciesList) + '\n')
			self.outDict[outFileName] = outFile
		return self.outDict[outFileName]

	def getState(self, repSum):
		"""
		Roll up the county data to state
		"""
		stateRep = {}
		for species in list(repSum.keys()):
			if species not in stateRep:
				stateRep[species] = {}
			for fips in list(repSum[species].keys()):
				stfips = fips[:8]
				if stfips not in stateRep[species]:
					stateRep[species][stfips] = repSum[species][fips]
				else:
					stateRep[species][stfips] += repSum[species][fips]
		return stateRep

	def writeOut(self, repSum, facs, Mon):
		"""
		Output to the afdust reports with SCC and without SCC
		Automatically determine whether to write monthly, annual, county, or state
		"""
		speciesList = sorted(repSum.keys())
		if Mon != 'annual':
			dateType = 'monthly'
			lineMon = [Mon,]
		else:
			dateType = timeLabel #'annual'
			lineMon = []

		for area in ['county', 'state']:
			if area == 'county':
				fipsDescDict = facs.countyDescDict
			elif area == 'state':
				repSum = self.getState(repSum)
				fipsDescDict = facs.stateDescDict
			
			fipsList = sorted(repSum[speciesList[0]].keys())

			for fips in fipsList:
				noSccOut = self.__loadOut(area, dateType, speciesList)
				noSccLine = lineMon + [fips, fipsDescDict[fips]['st']]

				if area == 'county':
					noSccLine.append(fipsDescDict[fips]['cty'])

				fipsSum = dict(list(zip(speciesList, [0 for x in speciesList])))

				for species in speciesList:
					try:
						val = repSum[species][fips]
					except:
						val = 0

					fipsSum[species] += val

				noSccLine += [str(fipsSum[species]) for species in speciesList]
				noSccOut.write('%s\n' %','.join(noSccLine))

	def emfWrite(self, repSum, facs):
		"""
		Write the county and state annual reports in _emf.csv format
		"""
		speciesList = sorted(repSum.keys())

		for area in ['county', 'state']:
			outFileName = '%s_%s_%s_adj_%s_%s_%s_emf.csv' %(timeLabel, case, sector, grid, spec, area)
			outPath = os.path.join(home, case + '/reports/annual_report')
			outFileName = os.path.join(outPath, outFileName)
			print('Writing EMF %s report to %s' %(timeLabel, outFileName))
			outFile = open(outFileName,'w')

			if area == 'county':
				fipsDescDict = facs.countyDescDict
				head = ['FIPS','State','County']
			else:
				repSum = self.getState(repSum)			
				fipsDescDict = facs.stateDescDict
				head = ['State',]

			head += ['Sector','Species','ann_emis']
			outFile.write('%s\n' %','.join(head))

			fipsList = sorted(repSum[speciesList[0]].keys())

			for fips in fipsList:

				if area == 'county':
					fipsLine = [fixLen(fips,6), fipsDescDict[fips]['st'], fipsDescDict[fips]['cty'], sector + '_adj']
				else:
					fipsLine = [fipsDescDict[fips]['st'], sector + '_adj']

				fipsSum = dict(list(zip(speciesList, [0 for x in speciesList]))) 
				for species in speciesList:
					try:
						val = repSum[species][fips]
					except:
						val = 0
					outLine = fipsLine + [species, str(val)]
					outFile.write('%s\n' %','.join(outLine))

# Check environment variables
Year = checkEV('BASE_YEAR')
grid = checkEV('GRID')
sector = checkEV('SECTOR')
spec = checkEV('EMF_SPC')
case = checkEV('CASE')
imPath = checkEV('IMD_ROOT')
premerged = checkEV('PREMERGED')
home = checkEV('PROJECT_ROOT')

# Load the months list
monList = checkEV('MONTHS_LIST') # Months to run as reported by the set_months script
monList = sorted([ int(mon.strip()) for mon in monList.split(' ') ])

# Detect whether script is being run annually or for only part of the year, and then set output label accordingly
# (added by C.Allen, 20 Sep 2018)
#print("first " + fixLen(monList[0]) + "  last " + fixLen(monList[-1]))
if fixLen(monList[0]) == "01" and fixLen(monList[-1]) == "12":
	timeLabel = 'annual'
elif fixLen(monList[0]) == fixLen(monList[-1]):
	timeLabel = 'month' + fixLen(monList[0])
else:
	timeLabel = 'month' + fixLen(monList[0]) + '-' + fixLen(monList[-1])
#print(timeLabel)

# Find the species used in this case and load the molecular weight conversion dictionary
if spec in mwDict:
	mwSpec = spec.lower()
else:
	mwSpec = mwDefault
print("Parameter file: " + mwDict[mwSpec])
execfile(mwDict[mwSpec])  # Load up the molecular weight conversion dictionary from the external file

######

# Main loop
facs = repFac()
reps = repOut()

annualCountyDict = {}
annualStateDict = {}	
annSum = {}

for mon in monList:
	Mon = fixLen(mon)
	print("Processing month: " + Mon)
	reportMonth = createReport(Mon) 
	monSum = reportMonth.sumSpeciesMonth(facs)

	# Write monthly county + state
	reps.writeOut(monSum, facs, Mon)

	# Roll up to annual
	for species in list(monSum.keys()):
		if species not in annSum:
			annSum[species] = {}
		for fips in list(monSum[species].keys()):
			if fips not in annSum[species]:
				annSum[species][fips] = monSum[species][fips]
			else:
				annSum[species][fips] += monSum[species][fips]

# Write annual county + state
reps.writeOut(annSum, facs, 'annual')
reps.emfWrite(annSum, facs)
