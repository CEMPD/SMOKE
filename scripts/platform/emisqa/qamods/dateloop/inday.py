from datetime import timedelta, date, datetime
import time 
import os.path

class InDay(object):
	"""
	Current date and object for looping.
	Takes a Gregorian date (YYYYMMDD) to run as input.

	Operations-
	obj.iterday(): Move to the next Gregorian date.
	obj: returns current Gregorian date when called as a string 
	"""

	def __init__(self, gsdate, rep_days, run_days, smkdates_path):
		self.today = date(int(gsdate[:4]), int(gsdate[4:6]), int(gsdate[6:8]))
		self.y = self.today.year
		self.m = self.today.month
		self.first_day = self.today.strftime('%Y%m%d')
		self.rep_days = rep_days
		self.smkdates_path = smkdates_path

		self.last_day = self.today + timedelta(run_days - 1)
		self.last_day = self.last_day.strftime('%Y%m%d')

		if self.rep_days: 
			self.date_dict = self._parse_smkdates()

	def __str__(self):
		"""
		When called as a string returns the representative day
		"""
		current_date = self.today.strftime('%Y%m%d')
		if self.rep_days:
			return "%s" %self.date_dict[current_date]['rep']
		else:
			return current_date

	def _parse_smkdates(self):
		"""
		Parse in the SMOKE dates file for this month.  Creating a dictionary containing representative day information.
		"""
		infile_name = os.path.join(self.smkdates_path, str(self.y), 'smk_merge_dates_%s%0.2d.txt' %(self.y, self.m))
		in_file = open(infile_name)

		for line in in_file:
			row = [record.strip().upper() for record in line.split(',')]
			if 'Date' in line:
				try: 
					col_num = row.index(self.rep_days)
				except ValueError: 
					raise ValueError('Representative type %s not found in SMOKE merge dates file.' %self.rep_days)
				else: 
					date_dict = {}
					rep_dict = {}
				continue

			in_date = row[0]
			rep_date = int(row[col_num])

			# Read only the days that fall in the date range
			if (int(in_date) >= int(self.first_day)) and (int(in_date) <= int(self.last_day)):
				if rep_date not in rep_dict:
					rep_dict[rep_date] = {'first': in_date, 'mult': 1}
				else:
					rep_dict[rep_date]['mult'] = rep_dict[rep_date]['mult'] + 1

				date_dict[in_date] = {'rep': rep_date, 'mult': 0}

		for rep_date in rep_dict:
			first_date = rep_dict[rep_date]['first']
			mult = rep_dict[rep_date]['mult']
			date_dict[first_date]['mult'] = mult

		return date_dict

	def current_date(self):
		"""
		Return the non-representative current date.
		"""
		return self.today.strftime('%Y%m%d')

	def current_mult(self):
		"""
		Return the current representative day multiplier for the day.
		"""
		if self.rep_days:
			mult = self.date_dict[self.today.strftime('%Y%m%d')]['mult'] 
		else:
			mult = 1

		return mult
	
	def iterday(self):
		"""
		Advance to the next Gregorian day.
		"""
		self.today = self.today + timedelta(1)
		self.y = self.today.year
		old_month = self.m
		self.m = self.today.month
		if self.rep_days:
			if old_month != self.m:
				self.date_dict = self._parse_smkdates()

