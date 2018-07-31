# embl to fasta
import re
import sys

if len(sys.argv) != 3:
	print('Usage: embl2fasta.py <embl_input_file> <fasta_output_name>')
	sys.exit()

IDset = False
inFile = sys.argv[1]	# 'example.embl'
outFile = sys.argv[2]	# 'example.fa'

output = open(outFile,'w')

with open(inFile) as f:
	for line in f:
		if line.startswith('ID') and not IDset:
			output.write('>'+line.lstrip('ID').split(';')[0].lstrip(' ')+'\n')
			IDset = True
		if line.startswith(' '):
			output.write(re.sub('[ 0-9\n]','', line).upper()+'\n')

output.close()
print('Done!')
