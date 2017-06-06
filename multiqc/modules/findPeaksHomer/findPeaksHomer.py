#!/usr/bin/env python

""" MultiQC module to parse output from findPeaks (HOMER) """

from __future__ import print_function
from collections import OrderedDict
import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table, beeswarm

# Initialize the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
	""" findPeaks module """

	def __init__(self):

		# Initialize the parent object
		super(MultiqcModule, self).__init__(name='findPeaksHomer', anchor='findPeaksHomer', href="http://homer.ucsd.edu/homer/ngs/peaks.html", info="finds enriched peaks, regions, and transcripts")

		# Find and load findPeaks reports
		self.findPeaksHomer_data = dict()
		for f in self.find_log_files('findPeaksHomer', filehandles=True):
			self.parse_findPeaks(f)

		# Filter to strip out ignored sample names
		self.findPeaksHomer_data = self.ignore_samples(self.findPeaksHomer_data)

		if len(self.findPeaksHomer_data) == 0:
			log.debug("Could not find any reports in {}".format(config.analysis_dir))
			raise UserWarning

		log.info("Found {} reports".format(len(self.findPeaksHomer_data)))

		# Write parsed report data to a file
		self.write_data_file(self.findPeaksHomer_data, 'multiqc_findPeaksHomer')

		# Basic Stats Table
		self.findPeaks_general_stats_table()

		# Plot
		self.add_section( plot = self.findPeaks_plot() )


	def parse_findPeaks(self, f):

		# Sample name
		name = re.search(r'(.+)\.peaks1kb',f['s_name'])
		if name:
			s_name = name.group(1)

		for l in f['f']:

			# Total peaks
			total = re.search(r'# total peaks = (\d+)',l)
			if total:
				total_peaks = float(total.group(1))

			# IP efficiency
			efficiency = re.search(r'# Approximate IP efficiency = (\d+\.\d+)%',l)
			if efficiency:
				ip_efficiency = float(efficiency.group(1))

			# Expected tags per peak
			tags = re.search(r'# expected tags per peak = (\d+\.\d+)',l)
			if tags:
				expected_tags = float(tags.group(1))

		if s_name in self.findPeaksHomer_data:
			log.debug("Duplicate sampel name found! Overwriting: {}".format(s_name))
		self.add_data_source(f, s_name)
		self.findPeaksHomer_data[s_name] = {
			'total_peaks': total_peaks,
			'ip_efficiency': ip_efficiency,
			'expected_tags': expected_tags
		}


	def findPeaks_general_stats_table(self):
		""" Take the parsed stats and add it to the basic stats table at the top of the report """

		headers = OrderedDict()
		headers['total_peaks'] = {
			'title': 'Total Peaks',
			'description': 'Total peaks'
		}
		headers['ip_efficiency'] = {
			'title': 'Approximate IP efficiency',
			'description': 'Approximate IP efficiency',
			'suffix': '%'
		}
		headers['expected_tags'] = {
			'title': 'Expected tags per peak',
			'description': 'Expected tags per peak'
		}
		self.general_stats_addcols(self.findPeaksHomer_data, headers)


	def findPeaks_plot(self):
		""" Plot the statistics """		
		
		headers = OrderedDict()
		headers['total_peaks'] = {
			'title': 'Total Peaks',
			'description': 'Total peaks'
		}
		headers['ip_efficiency'] = {
			'title': 'Approximate IP efficiency',
			'description': 'Approximate IP efficiency',
			'suffix': '%'
		}
		headers['expected_tags'] = {
			'title': 'Expected tags per peak',
			'description': 'Expected tags per peak'
		}

		return table.plot(self.findPeaksHomer_data, headers)	
