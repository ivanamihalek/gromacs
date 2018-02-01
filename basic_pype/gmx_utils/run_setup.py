import argparse
import sys

#########################################
def parse_commandline():
    if False and len(sys.argv)<=1:
        print "Usage:"
        print "%s -h" % sys.argv[0]
        print "or"
        print "%s -p <pdb_root> [-l <ligand list>] [--remd] [--min] [--posres]" % sys.argv[0]
        exit(1)

    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pdb', metavar='pdbfile', type=str,  default="none",
                        dest="pdb", help='input pdb file (peptide, ions, and ligands)')
    parser.add_argument('-l','--ligands', metavar='ligand', type=str, nargs='+', default=[],
                        dest="ligands", help='list of ligand names')
    parser.add_argument('-m','--min', action='store_true', default=False,
                        dest="minimization", help='energy minimization run (defaults to False)')
    parser.add_argument('-r','--remd', action='store_true', default=False,
                        dest="remd", help='REMD run (defaults to False)')
    parser.add_argument('-s','--posres', action='store_true', default=False,
                        dest="posres", help='position restrained run (defaults to False)')

    run_options = parser.parse_args(sys.argv[1:])

    return run_options
