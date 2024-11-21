from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
import sys

def percentError(speciesList, fileName):
	'''
	Calculate the percentError between two NCF files.
	'''
	sys.exit('PE not yet implemented')
	if len(fileName.split(',')) != 2: \
		sys.exit('ERROR: You must specify two input filenames using the -f argument.')
	file1 = dataFile(fileName.split(',')[0])
	file2 = dataFile(fileName.split(',')[1])

	print('Running percent error calculation...')
	outDict = {}  # Create the output dictionary that will be of shape { row: { col: { speciesname: val ... } ... } ... }

	#  Run percent error a special way to get correct state and county percentages
	if region == 'state' or region == 'county':
		if not grid: 
			sys.exit('No grid specified.  Needed to write state or county based csv.')
		ratioTable = parseRatio(region)
		outFile = open(outFileName, 'w')
		outFile.write('fips,' + ','.join(speciesName for speciesName in speciesList) + '\n')
		fipsList = sorted(ratioTable.keys())
		for fips in fipsList:
			outLine = fips
			lineDict = {}
			for speciesName in speciesList:
				if verbosity: 
					print('Getting percent error for species: %s' %speciesName)
				if speciesName not in lineDict:
					lineDict[speciesName] = { 'a1': 0, 'a2': 0}
				array1 = speciesArray(file1.sumVal(speciesName, layer, allHours), speciesName)
				array2 = speciesArray(file2.sumVal(speciesName, layer, allHours), speciesName)
				for cell in ratioTable[fips]:
					col = int(cell.split(',')[0])
					row = int(cell.split(',')[1])
					if col == 1: 
						print('Col')
					if row == 1: 
						print('Row')
					if row not in list(range(array1().shape[0])) or col not in list(range(array1().shape[1])): 
						continue  # Skip ratios that are outside the domain
					a1val = array1()[row][col] * ratioTable[fips][cell]
					a2val = array2()[row][col] * ratioTable[fips][cell]
					lineDict[speciesName]['a1'] = lineDict[speciesName]['a1'] + a1val
					lineDict[speciesName]['a2'] = lineDict[speciesName]['a2'] + a2val
				if lineDict[speciesName]['a1'] == 0:
					lineDict[speciesName]['pe'] = 0
				else:
					lineDict[speciesName]['pe'] = round(((old_div((lineDict[speciesName]['a2'] - lineDict[speciesName]['a1']), lineDict[speciesName]['a1'])) * 100), 2)
				if tons: 
					outLine = '%s,%s' %(outLine, moles2tons(lineDict[speciesName]['pe'], speciesName))
				else: 
					outLine = '%s,%s' %(outLine, lineDict[speciesName]['pe'])
			outFile.write(outLine + '\n')
		sys.exit(0)
	# Run percent error the standard way for gridded percentages 
	else:
		for speciesName in speciesList:
			if verbosity: 
				print('Getting percent error for species: %s' %speciesName)
			array1 = speciesArray(file1.sumVal(speciesName, layer, allHours), speciesName)
			array2 = file2.sumVal(speciesName, layer, allHours)
			PE = array1.pctErr(array2)

			if tons: 
				outDict[speciesName] = moles2tons(PE, speciesName)
			else:
				outDict[speciesName] = PE 

		return outDict

