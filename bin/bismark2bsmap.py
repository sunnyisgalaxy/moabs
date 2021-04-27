import os
import gzip
import pysam
import argparse
import sys
import logging

DEBUG=False
NotDEBUG=not DEBUG
parser = argparse.ArgumentParser(description="convert the bam file from the bismark output format to bsmap format with ZS tag",
							 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', action='store', nargs='?', help='bam file from bismark mapping', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="bam file of bsmap format", required=NotDEBUG)

args = parser.parse_args()
if DEBUG: 
	args.input = ''
	args.output = ''
	

logger = logging.getLogger('bismark2bsmap')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

ibam = args.input
obam = args.output
#temp = obam+'.temp'

save = pysam.set_verbosity(0)
with pysam.AlignmentFile(ibam,'rb') as f:

	wf = pysam.AlignmentFile(obam, "wb", template=f)
	logger.info("Processing reads ...")

	for line in f:
		#flag = line.flag
		#mismatch = line.get_tag('NM');
		XR = line.get_tag('XR');
		XG = line.get_tag('XG');
		
		if XR == 'GA' and XG == 'GA':
			line.set_tag('ZS','--')
		elif XR == 'CT' and XG == 'GA':
			line.set_tag('ZS','-+')
		elif XR == 'CT' and XG == 'CT':
			line.set_tag('ZS','++')
		elif XR == 'GA' and XG == 'CT':
			line.set_tag('ZS','+-')	
		else:
			logger.info("unexpected tags in bismark bam")

		wf.write(line)

	logger.info("Adding tag ZS is done.")

	# 
	# try:
		# wf.header['HD']['SO']
	# except KeyError:
		# logger.info('sorting bam file')
		# pysam.sort(temp, temp + '.sorted')
		# os.remove(temp)
		# os.rename(temp + 'sorted.bam', obam)
		# pysam.index(obam)
	# else:
		# if test_head.header['HD']['SO'] == 'coordinate':
			# pass
		# else:
			# logger.info('sorting bam file')
			# pysam.sort(temp, temp + '.sorted')
			# os.remove(temp)
			# os.rename(temp + 'sorted.bam', obam)
		# pysam.index(obam)

	logger.info("bismark2bsmap is done. Maybe you could sort the bam file if the original bam was not sorted!")

pysam.set_verbosity(save)

