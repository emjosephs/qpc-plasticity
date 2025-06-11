#read through the .hdf file and get out 500,000 random snps for making a kinship matrix

import h5py, numpy
import sys

def main():

	#read in file
	#f = h5py.File('1001_SNP_MATRIX/imputed_snps_binary.hdf5','r')
	f = h5py.File(_args.input,'r')

	#randomly pick n sites
	positions = f['positions'][:]
	n = _args.num
	#mysamp = numpy.random.choice(positions, n, replace=False)
	mysamp = numpy.random.choice(range(0, len(positions)-1), n, replace=False)
	#mysnps = [f['snps'][x] for x in mysamp]

	#write out
	myOut = open(_args.output, 'w')
	myOut.write('position	'+'	'.join(map(str, f['accessions'][:].tolist())))
	for x in mysamp:
		mysnp = f['snps'][x]
		mypos = positions[x]
		myOut.write("\n"+str(mypos)+"	"+'	'.join(map(str, mysnp.tolist())))


import argparse
_args=None
def parseArgs():
	parser=argparse.ArgumentParser(description = "pick snps from hd5 file")
	parser.add_argument('-i', '--input',type=str)
	parser.add_argument('-n', '--num',type=int)
	parser.add_argument('-o', '--output',type=str)

	global _args
	_args = parser.parse_args()
	sys.stderr.write(str(_args)+"\n")	

if __name__ == "__main__":
	parseArgs()
	main()



