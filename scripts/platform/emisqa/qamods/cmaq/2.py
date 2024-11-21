import numpy as np
import scipy.io.netcdf as ncf
import gzip, random
from sys import exit

class NCFFile(object):
	"""
	Instance of a NCFFile.
	"""
	def __init__(self, inFileName, verbosity = False, zipDict = {}):
		self.inFileName = inFileName
		self.verbosity = verbosity
		self.NCF = self.openNCF(zipDict)
		self.sdate = getattr(self.NCF, 'SDATE')
		self.speciesList = self.NCF.variables.keys()

	def __str__(self):
		return self.inFileName

	def __call__(self):
		return self.NCF

	def getSpecies(self, speciesName, grid = '', gridDesc = '', ignoreSpec = False, inln = False, stacks = ''):
		if speciesName not in self.speciesList:
			if ignoreSpec == True:
				print 'WARNING: The species %s does not exist in the file %s.' %(speciesName, self.inFileName)
				specShape = self.NCF.variables[self.speciesList[0]].shape
				species = np.zeros(specShape, '>f')
			else:
				exit('ERROR: The species %s does not exist in the file %s.' %(speciesName, self.inFileName))
		else:
			species = self.NCF.variables[speciesName]

		dataIn = species[:]

		# Place inline data into a 2D grid
		if inln == True:
			dataIn = self.gridInln(dataIn, stacks, grid, gridDesc)

		return dataIn

	def openNCF(self, zipDict = {}):
		'''
		Opens the netCDF input file and returns an open file object.
		'''
		if self.verbosity == True: 
			print 'Opening file for reading: %s' %self.inFileName
		try: 
			fileIn = ncf.netcdf_file(self.inFileName, 'r')
		except:
			print 'WARNING: %s not available for access.  Attempting to open zipped version.' %self.inFileName
			if self.verbosity == True: 
				print 'Opening file for reading: %s.gz' %self.inFileName
			# Try to read a zipped version of the input file
			try: 
				zipIn = gzip.open(self.inFileName+'.gz','rb')
			except: 
				exit('ERROR: %s.gz not available for access.' %self.inFileName)
			else:
				if self.inFileName in zipDict:
					# Check if the file has already been unzipped
					print 'Previously unzipped version found.'
					fileIn = ncf.netcdf_file(zipDict[self.inFileName], 'r')
					return fileIn
				tmpFileName = '%s/pyqa.ncf-%s' %(tmpDir, str(random.randint(100, 9999999)))
				tmpFile = open(tmpFileName, 'wb')
				for data in dataBlocks(zipIn):
					tmpFile.write(data)
				zipIn.close()
				tmpFile.close()
				zipDict[self.inFileName] = tmpFileName  # Store the unzipped temporary file name in the zip dictionary
				fileIn = ncf.netcdf_file(tmpFileName, 'r')
				return fileIn
		else: 
			return fileIn

	def gridInln(self, dataIn, stacks, grid, gridDesc):
		"""
		Process the input species and adjust based on the ratio table
		"""
		gridInfo = getGrid(grid, gridDesc)
	 
		dataOut = np.zeros([dataIn.shape[0],1,gridInfo['rows'],gridInfo['cols']], '>f4')
 
		for stack in range(stacks.stkNum):
			row = stacks.ptXref[stack]['row'] - 1
			col = stacks.ptXref[stack]['col'] - 1
 
			if row not in range(dataOut.shape[2]) or col not in range(dataOut.shape[3]):
				print 'stack: %s at col: %s row: %s outside of bounds' %(stack + 1, col + 1, row + 1)
				continue
 
			try:
				dataOut[:,0,row,col] = dataIn[:,stack,0] + dataOut[:,0,row,col]
			except:
				exit('ERROR: Inline to grid problem at: ROW %s COL %s STACK %s' %(row,col,stack))
	 
		return dataOut

	def closeFile(self):
		'''
		Closes the open file
		'''
		self.NCF.close()

