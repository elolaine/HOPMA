# -*- coding: utf-8 -*-

# Copyright (c) 2020: Elodie Laine
# This code is part of the gemme package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

import sys
import os
import argparse
import re
import subprocess
import math
from hopmaAnal import *


def check_argument_groups(parser, arg_dict, group, argument):
    """
    Check for use of arguments.
    Raise an error if the parser uses an argument that belongs to an argument
    group (i.e. mode) but the group flag is not used or if the argument is
    required in the argument group and it's not used.
    Notes:
    - argument and group should be strings and should start with - or --
    - group should be a flag with action 'store_true'
    - argument shouldn't be a flag, and should be None if it isn't used.
    >>> import argparse
    >>> parser = argparse.ArgumentParser()
    >>> parser.add_argument('--phylo', action='store_true')
    _StoreTrueAction(...)
    >>> parser.add_argument('--inseq')
    _StoreAction(...)
    >>> args = parser.parse_args(['--phylo', '--inseq', 'not_none'])
    >>> check_argument_groups(parser, vars(args), '--phylo', '--inseq', True)
    """
    c = 0
    for myArg in argument:
    	arg_name = myArg.replace("-", "")
    	if arg_dict[arg_name]!='':
    		c = c+1
    group_name = group.replace("-", "")
    if arg_dict[group_name]=="input":
    	if c!=1:
    		parser.error("HOPMA requires " + str(argument) +
                         " if " + group + " is set to input.")
    else:
        if c>0:
            parser.error("HOPMA requires " + group +
                         " to be set to input if " + str(argument) + " is used.")
    return None

def parse_command_line():
    """
    Parse command line.
    It uses argparse to parse phylosofs' command line arguments and returns the
    argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="HOPMA",
        description="""
        HOPMA (proteiN cOloRed MAps) 
        is a tool to build protein elastic networks with enhanced dynamical potential
        """,
        epilog="""
        If you use it, please cite:
        Laine E, Grudinin S. HOPMA: Boosting protein structures dynamical potential with colored contact maps. submitted.
        """,
    )

    parser.add_argument(
    	'input', metavar='pdbFile', type=str,
        help='input 3D structure file in PDB format. Please provide the name of the file without the .pdb extension',
    )

    parser.add_argument(
        '-c', '--chain',
        help='name of the protein chain to be analysed ("A" by default)',
        default='A'
    )

    retMet_args = parser.add_argument_group('Detection of the contact patches', """
        Arguments used for detecting the contact patches to be removed .
        """)
    
    retMet_args.add_argument(
        '-s', '--size',
        help='minimum size allowed for a patch',
        default='625'
    )

    retMet_args.add_argument(
        '-k', '--padding',
        help='padding for the smoothing window',
        default='2'
    )

    retMet_args.add_argument(
        '-p', '--propensity',
        help='cutoff value for the contact propensity matrix. By default, the cutoff values 0.1 and 0.001 are combined',
        default='combi'
    )

    retMet_args.add_argument(
        '-dmin', '--dmin',
        help='minimal distance bound',
        default='3'
    )
        
    retMet_args.add_argument(
        '-dmax', '--dmax',
        help='maximal distance bound',
        default='15'
    )

    args = parser.parse_args()

    arg_dict = vars(args)

    #check_argument_groups(parser, arg_dict, '--retrievingMethod', ['--blastFile','--fastaFile'])
 
    # Check flag arguments
    if args.input is None:
    	parser.error("HOPMA requires an input protein 3D structure.")

    return args


def doit(pdb_code,chain,size,k,cutoff,dmin,dmax):
    """
    doit(args.input,args.chain,args.size,args.padding,args.propensity,args.dmin,args.dmax)
    """
    preprocessInput(pdb_code, chain)
 #   print "computing maps..."
    createColoredMaps(pdb_code, chain, k, cutoff, dmin, dmax)
    extractExcludedContacts(pdb_code, chain, k, size, cutoff)
 #   print "done"


def main():
    args = parse_command_line()
    doit(args.input,args.chain,args.size,args.padding,args.propensity,args.dmin,args.dmax)

if (__name__ == '__main__'):
    main()

