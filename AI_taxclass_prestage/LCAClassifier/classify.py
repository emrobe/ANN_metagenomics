'''
Created on Sep 22, 2011

Classifier based on assignment to LCA (Lowest Common Ancestor), using BLAST
homology search results to a reference database (currently bit-score).

@author: andersl
'''
import sys
from optparse import OptionParser

from Bio import SeqIO
from Bio.Blast import NCBIXML

from taxa import Node, Read, Tree
from config import config

BSR_DEFAULT = 0.98
MS_DEFAULT = 155

sfLimits = {Tree.SPECIES: .99, # 97
            Tree.GENUS: .97, # 95
            Tree.FAMILY: .95, # 90
            Tree.ORDER: .90, # 85
            Tree.CLASS: .85, # 80
            Tree.PHYLUM: .80   # 75
            }


class LCAClassifier(Tree):
    """A classifier instance inherits from the standard Tree and assigns
    reads to it representing the sequence reads of the classified dataset"""

    def __init__(self, name, minFilter=True, fastafile=None, qualfile=None):
        Tree.__init__(self, Node(name="root", nodeID=1), name=name)
        self.bsr = BSR_DEFAULT
        self.ms = MS_DEFAULT
        self.minFilter = minFilter
        self.noHits = Node("No hits", parent=self.root)
        self.addNode(self.noHits)
        self.seqs = {}
        self.qual = {}

        if fastafile:
            fstream = open(fastafile, 'r')
            for seq_record in SeqIO.parse(fstream, "fasta"):
                self.seqs[seq_record.id] = seq_record.seq
                
        if qualfile:
            qstream = open(qualfile, 'r')
            for q_record in SeqIO.parse(qstream, "qual"):
                self.qual[q_record.id] = q_record

    def assign(self, records, dataset=None, verbose=False):
        """Accepts a of biopython blast iterator and carries out LCA
        assignments to a given dataset"""

        for record in records:
            qName = record.query.split(" ")[0]
            if (record.alignments and
                record.alignments[0].hsps[0].bits >= self.ms):
                best_hsp = record.alignments[0].hsps[0]
                topScore = best_hsp.bits
                if self.seqs and qName in self.seqs.keys():
                    qSeq = self.seqs[qName]
                else:
                    qSeq = str(best_hsp.query).replace("-", "")
                hitname = record.alignments[0].hit_def.split()[0]
                node = self.getNode(hitname)
                if not node:
                    sys.stderr.write("Best-scoring node %s not found!\n" %
                                     hitname)
                    sys.stderr.write("Cannot assign read %s\n" % qName)
                else:
                    parents = node.getPhylogeny()[1:]

                    # Iterate through rest of hits until falling below treshold
                    for a in record.alignments[1:]:
                        if a.hsps[0].bits < float(topScore) * self.bsr:
                            break
                        hitname = a.hit_def.split()[0]
                        n = self.getNode(hitname)

                        if not n:
                            sys.stderr.write("Node " + hitname +
                                             " not found! Ignoring.\n")
                        else:
                            p = n.parent
                            # iterate through parents until found in the
                            # parents list
                            while p not in parents:
                                p = p.parent
                            parents = parents[parents.index(p):]

                    # Take a look at similarity, print info if verbose and
                    # kick up if filter
                    hsp_sim = (float(best_hsp.identities) /
                               float(best_hsp.align_length))
                    if verbose and hsp_sim >= .99:
                        print ("Read %s in %s is %s percent similar to %s" %
                               (qName, dataset, hsp_sim * 100,
                                record.alignments[0].hit_def))

                    if self.minFilter:
                        maxRankLimit = Tree.SPECIES
                        maxRank = maxRankLimit
                        d = maxRankLimit
                        ranks = sfLimits.keys()
                        ranks.sort()
                        ranks.reverse()
                        for rank in ranks:
                            if hsp_sim < sfLimits[rank]:
                                maxRank = rank - 1
                            else:
                                break

                        while (maxRank < Tree.SPECIES and
                               maxRank < parents[0].getHighestRank()):
                            d = min(parents[0].getHighestRank(), maxRankLimit)
                            if verbose:
                                print ("Read %s in %s cannot be assigned to "
                                       "rank %s (similarity=%s)" %
                                       (qName, dataset, Tree.depths[d],
                                        hsp_sim))
                            parents = parents[1:]

                        if d < maxRankLimit:
                            novelName = ("Unknown %s %s" %
                                         (parents[0].name, Tree.depths[d]))
                            nn = self.getNode(novelName)
                            if nn:
                                novelNode = nn
                            else:
                                depth = parents[0].getHighestRank() + 1
                                novelNode = Node(novelName, parent=parents[0],
                                                 depth=depth)
                                self.addNode(novelNode)
                            parents = [novelNode] + parents

                    # Handle assignment
                    read = Read(qName, seq=qSeq)
                    parents[0].assignRead(read, dataset=dataset, primary=True,
                                          recursive=True)
            else:
                #No hits
                if self.seqs and qName in self.seqs.keys():
                    qSeq = self.seqs[qName]
                elif record.alignments:
                    qSeq = record.alignments[0].hsps[0].query.replace("-", "")
                else:
                    qSeq = None
                nhr = Read(name=qName, seq=qSeq)
                self.noHits.assignRead(nhr, dataset=dataset, primary=True)
                self.root.assignRead(nhr, dataset=dataset, primary=False)

    def setBitscoreRange(self, percent):
        self.bsr = 1 - float(percent) / 100

    def setMinScore(self, minScore):
        self.ms = minScore
        
    def printAssignmentsRDPQual(self, node, dataset=None, printFile=None,
                                 newTabStyle=False):
        assignments = node.getAssignment(dataset)
        if assignments and assignments.primReads:

            for r in assignments.primReads:
                toPrint = (">%s\t%s\n" %
                           (r.name,
                            node.getPhylogenyRDPStyle(
                                      root=False, newTabStyle=newTabStyle)))
                for line in self.qual[str(r.name)].format("qual").split("\n")[1:]:
                    toPrint+=(line + " ")
                if printFile:
                        printFile.write(toPrint + "\n")
                else:
                    print toPrint
        if assignments:
            for child in node.children:
                self.printAssignmentsRDPQual(node=child, dataset=dataset, 
                                              printFile=printFile,
                                              newTabStyle=newTabStyle)


def main():

    parser = OptionParser()

    parser.add_option("-d", "--dbname",
                      dest="dbname",
                      type="string",
                      default="silvamod",
                      help="the taxonomy to be used (defaut = silvamod)")

    parser.add_option("-r", "--range",
                      dest="bitScoreRange",
                      type="int",
                      default=2,
                      help=("bitscore-range (the range of blast hits to find "
                            "LCA of, given in percent drop from highest "
                            "score; default = 2)"))

    parser.add_option("-s", "--minscore",
                      dest="minScore",
                      type="int",
                      default=MS_DEFAULT,
                      help=("minimum bit-score given as an integer; "
                            "default = 155"))

    parser.add_option("-n", "--normalisetobase",
                      action="store_true", dest="normBase",
                      default=False,
                      help=("Normalise domain level and lower to those "
                            "classified at base level"))

    parser.add_option("-f", "--nofilter",
                      action="store_false", dest="minFilter",
                      default=True,
                      help=("Deactivate minimum identity filter preventing "
                            "classification to higher ranks when a minimum "
                            "rank-identity is not met (3% for species, 5% for "
                            "genera, 10% for family"))

    parser.add_option("-u", "--includeunknown",
                      action="store_true", dest="outputNovel",
                      default=False,
                      help=('When the minimum percent identity is is not met, '
                            'sequences are classified as "unknown". By '
                            'default this is not hidden from the composition '
                            'table.'))

    parser.add_option("-o", "--overview",
                      action="store_true", dest="altCompo",
                      default=False,
                      help="Print tcomposition overview in alternative \
                            parser-friendly tab-separated text format.")

    parser.add_option("-p", "--rdp",
                      action="store_true", dest="rdpOut",
                      default=False,
                      help="Print Classifications in RDP Classifier format")

    parser.add_option("-a", "--fasta",
                      action="store_true", dest="faOut",
                      default=False,
                      help="Write fasta file with classified sequences")

    parser.add_option("-i", "--fastain",
                      dest="fastafile",
                      type="str",
                      default=None,
                      help=("Output sequences from submitted FASTA file with "
                            "classification instead of aligned ones from "
                            "blast result"))
    
    parser.add_option("-q", "--qualin",
                      dest="qualfile",
                      type="str",
                      default=None,
                      help=("Output sequences from submitted FASTA quality " 
                      "file with annotation"))

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      default=False,
                      help=("Verbose mode (output information about "
                            "interesting assignments"))

    parser.add_option('-c', '--config',
                      dest='config',
                      default=None,
                      help=('Configuration file.'))
                   
    parser.add_option("-1", 
                      dest="outp_composition",
                      type="str",
                      default=None,
                      help=("."))
    parser.add_option("-2", 
                      dest="outp_tree",
                      type="str",
                      default=None,
                      help=("."))
    parser.add_option("-3", 
                      dest="outp_assignments",
                      type="str",
                      default=None,
                      help=("."))

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("No blast-file specified!")

    if options.config is not None:
        config.configure(options.config)

    #Initiate from tree
    lca = LCAClassifier("reftree", minFilter=options.minFilter,
                        fastafile=options.fastafile, qualfile=options.qualfile)

    mapFile = ("%s/%s.map" %
               (config.DATABASES[options.dbname], options.dbname))
    treFile = ("%s/%s.tre" %
               (config.DATABASES[options.dbname], options.dbname))

    lca.initFrom(mapFile, treFile) #, options.synFile)
    lca.setBitscoreRange(options.bitScoreRange)
    lca.setMinScore(options.minScore)
    datasets = []

    for a in args:
        name = a.replace(".xml", "").replace(".XML", "")
        datasets.append(name)
        results = open(a)

        #Parse blast (only XML)
        records = NCBIXML.parse(results)

        #Classify!
        lca.assign(records, dataset=name, verbose=options.verbose)
        results.close()

    #Crop trees from those nodes with no reads asssigned in ANY LCAClassifier
    lca.pruneUnassigned()

    #Write results
    if options.altCompo:
        altFile = open("All_Composition.txt", 'w')
        lca.printCompositionAlternative(altFile, datasets, options.outputNovel)
        altFile.close()

    for name in datasets:
        #tabFile = open(name + "_Composition.txt", 'w')
        tabFile = open(options.outp_composition, 'w')
        if options.rdpOut:
            rdpFile = open(options.outp_assignments, 'w')
        treeFile = open(options.outp_tree, 'w')
        if options.faOut:
            faFile = open('/dev/null', 'w')
        if options.qualfile:
            qualFile= open("/dev/null", 'w')
        if options.rdpOut:
            lca.printAssignmentsRDPStyle(name, rdpFile)
        if options.faOut:
            lca.printAssignmentsRDPFasta(name, faFile)
        if options.qualfile:
            lca.printAssignmentsRDPQual(node=lca.root, dataset=name, printFile=qualFile)
        lca.printAsTree(popDataset=name, showLeaves=True, printFile=treeFile)

        for level in Tree.depths:
            tabFile.write("Assingments at %s level\n\n" %
                          Tree.depths[level])
            lca.printPopulationsAtDepth(level, outFile=tabFile, dataset=name,
                                        outputNovel=options.outputNovel,
                                        normaliseToBase=options.normBase)
            tabFile.write("\n")
        tabFile.close()

        if options.rdpOut:
            rdpFile.close()
        treeFile.close()
        if options.faOut:
            faFile.close()


if __name__ == "__main__":

    main()
