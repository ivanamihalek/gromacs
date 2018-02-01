#!/usr/bin/python -u

from gmx_utils import run_setup

#########################################
def main():
    run_options = run_setup.parse_commandline()
    print run_options.pdb
    print run_options.ligands
    print run_options.minimization
    return True


#########################################
if __name__ == '__main__':
    main()
