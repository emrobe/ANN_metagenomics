import subprocess
import argparse
import urlparse
import shlex

#Needs to be set manually
python_path = 'python'
blastn_path = 'blastn'
silvamod = 'silvamod.fasta'

#Get task id which references current dataset in table.tsv
parser = argparse.ArgumentParser()
parser.add_argument('task', help="Task ID for the current job (required)", type=int)
parser.add_argument('-n', help="Number of cores to use. (Default =1)", type=int, default=1)
args = parser.parse_args()


# Parse overview
overview = open("target_samples/metadata.csv")
data = []
for line in overview:
	entry = line.split("\t")
	if entry[0].startswith(str(args.task)):
		data = entry
		print data
		break

# Download dataset
urlprefix = 'http://api.metagenomics.anl.gov/1/download/'
dataset_id = data[1]
rrna_out = dataset_id+'16s_out'
blast_out = dataset_id+'blast_out'
lca_out = dataset_id+'lca_out'
urlsuffix = '?file=150.1'
tmp = urlparse.urljoin(urlprefix, dataset_id)
url = urlparse.urljoin(tmp, urlsuffix)

# Produce classification.
#wget_command = ['wget', '-U', 'USER', '-O', str(dataset_id), '-q',str(url)]
rrna_pred_command = [python_path, '16S-tool/filter_16s_new.py', '-t', '-n', str(args.n), '-o', rrna_out, str(dataset_id)]
megablast_command = [blastn_path, '-query', rrna_out, '-task', 'megablast', '-db', silvamod, '-out', blast_out, '-num_threads', str(args.n), '-outfmt', '5', '-num_alignments', '100']
lca_command = [python_path, 'LCAClassifier/classify.py', blast_out, '-1', lca_out, '-2', '2', '-3', '3', '-p']
#subprocess.call(wget_command)
subprocess.call(rrna_pred_command)
subprocess.call(megablast_command)
subprocess.call(lca_command)

#Clean up, move
clean = ['rm', rrna_out, blast_out, '3', '2']
result_path = './target_samples/results/'+dataset_id
move = ['mv', lca_out, result_path]
subprocess.call(clean)
subprocess.call(move)

