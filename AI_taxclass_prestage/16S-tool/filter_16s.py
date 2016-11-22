import os, sys
from Bio import SeqIO
from Bio.Seq import Seq


path = '/'.join(sys.argv[0].split('/')[:-1])
BACTERIA_FWD = os.path.join(path, 'hmm/16s_bact_for3.hmm')
BACTERIA_REV = os.path.join(path, 'hmm/16s_bact_rev3.hmm')
ARCHAEA_FWD = os.path.join(path, 'hmm/16s_arch_for3.hmm')
ARCHAEA_REV = os.path.join(path, 'hmm/16s_arch_rev3.hmm')


def do_search(inp, profile, matches, do_trim=True, reverse=False, cpu=1):
	# Do HMM search
	cmd = 'hmmsearch --cpu {2} --notextw -E 10E-5 --noali {1} {0}'
	cmd = cmd.format(inp, profile, cpu)
	sys.stdout.flush()
	pipe = os.popen(cmd)
	lines = pipe.readlines()
	pipe.close()
	for i in range(len(lines)):
		line = lines[i]
		if line[0:3] == '>> ':
			skip_symbols = ('', '..', '[.', '.]', '[]')
			val = lines[i+3].strip().split(' ')
			val = [v for v in val if (v not in skip_symbols)]
			if int(val[7])-int(val[6]) < 100:
				continue
			pos1 = int(val[8])
			pos2 = int(val[9])
			matches[line.strip()[3:]] = (pos1, pos2, reverse)

def extract_sequences(matches, inp, do_trim=True):		
	# Extract matching sequences from input fasta file
	fasta_file = open(inp, 'r')
	matching_headers = matches.keys()
	for rec in SeqIO.parse(fasta_file, 'fasta'):
		if rec.description in matching_headers:
			matching_headers.remove(rec.description)
			reverse = matches[rec.description][2]
			print '>'+rec.description
			if not do_trim:
				if reverse:
					print rec.seq.reverse_complement().tostring()
				else:
					print rec.seq.tostring()
			elif do_trim:
				pos1 = matches[rec.description][0]
				pos2 = matches[rec.description][1]
				if reverse:
					print rec.seq.reverse_complement().tostring()[pos1-1:pos2]
				else:
					print rec.seq.tostring()[pos1-1:pos2]


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser()
	parser.add_argument('input', help="Input FASTA file", type=str)
	parser.add_argument('-a', help="Only screen for archaeal 16S rRNA", action="store_true")
	parser.add_argument('-b', help="Only screen for bacterial 16S rRNA", action="store_true")
	parser.add_argument('-t', help="Trim sequences to regions matching 16S rRNA", action="store_true")
	parser.add_argument('-n', help="Number of CPUs to use", type=int, default=1)
	parser.add_argument('-o', help="Output file name", nargs=1, type=str)
	args = parser.parse_args()
	
	if args.o:
		f = open(args.o[0], 'w')
		sys.stdout = f
		
	matches = {}
	
	if args.a and not args.b:
		do_search(args.input, ARCHAEA_FWD, matches, cpu=args.n)
		do_search(args.input, ARCHAEA_REV, matches, reverse=True, cpu=args.n)
	elif args.b and not args.a:
		do_search(args.input, BACTERIA_FWD, matches, cpu=args.n)
		do_search(args.input, BACTERIA_REV, matches, reverse=True, cpu=args.n)
	else:
		do_search(args.input, ARCHAEA_FWD, matches, cpu=args.n)
		do_search(args.input, ARCHAEA_REV, matches, reverse=True, cpu=args.n)
		do_search(args.input, BACTERIA_FWD, matches, cpu=args.n)
		do_search(args.input, BACTERIA_REV, matches, reverse=True, cpu=args.n)
	
	extract_sequences(matches, args.input, args.t)
	
	if args.o:
		sys.stdout.flush()
		f.close()






