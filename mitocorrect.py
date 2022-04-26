#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Thomas J. Creedy
"""

# Imports
import argparse
import textwrap as _textwrap
import mitocorrect_modules as mcm
import shutil
import multiprocessing
import functools

# Class definitions


class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width,
                                                 initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

# Function definitions


def required_multiple(multiple):
    class RequiredMultiple(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not len(values) % multiple == 0:
                msg = 'argument "{f}" requires a multiple of {multiple} values'
                msg = msg.format(f=self.dest, multiple=multiple)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)

    return RequiredMultiple

# Argument parser


parser = argparse.ArgumentParser(description=
            """Mitocorrect corrects annotations on a mitochondrial genome, supplied as one or more
            genbank-format flat file(s). Annotations are corrected according to the parameters 
            defined in a specifications table, supplied as a tab-separated values file. As part of
            the process of assessing potential annotations, mitocorrect aligns potential gene 
            sequences with a reference (aka profile) alignment for a given gene. These should be 
            constructed prior to running and paths given in another tsv file. MAFFT should be 
            installed and on the PATH available to mitocorrect; biopython should also be installed. 
            |n
            Mitocorrect is multithreaded so outputs relatively little intermediate information. 
            Detailed information can be viewed in the log file during running and/or in the results
            file upon completion. The corrected annotations are output in one or several genbank-
            format flat file(s).
            |n
            See the readme at https://github.com/tjcreedy/mitocorrect for more information 
            """,
                                 formatter_class=MultilineFormatter)

parser._optionals.title = "arguments"

# Required
parser.add_argument('-s', '--specifications', type = str, metavar = 'path', required = True,
                    help = "path to a specifications file in tab-separated variables format")
parser.add_argument('-g', '--genbank', type = str, nargs = '+', required = True, metavar = 'path',
                    help = "path(s) to genbank-format files containing annotated mitogenomes to "
                           "correct")
parser.add_argument('-a', '--alignmentpaths', type = str, required = True, metavar = 'path',
                    help = "path to a two column tab-separated variables file giving the name and"
                           "path to a profile alignment for each gene")
parser.add_argument('-b', '--translationtable', type = int, required = True, metavar = 'n',
                    help = "the appropriate genetic code translation table number (see "
                           "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi")

# Optional
parser.add_argument('-t', '--threads', type = int, default = 1, metavar = 'n',
                    help = "number of simultaneous threads to use for correcting mitogenomes "
                           "(default: 1).")
parser.add_argument('-c', '--alignmenttype', type = str, default = 'nt', choices = ['aa', 'nt'],
                    help = "the type of sequence data used for the supplied profile alignments ( "
                           "default: nt)")
parser.add_argument('-o', '--outputdirectory', type = str, default = 'mitocorrect_output',
                    metavar = 'path',
                    help = "path to a directory (created if needed) to which outputs will be "
                           "written (default: a directory called mitocorrect_output)")
parser.add_argument('-l', '--logfile', type = str, default = 'mitocorrect.log', metavar = 'name',
                    help = "name of a file to which log messages will be written, saved to the "
                           "output directory (default: mitocorrect.log)")
parser.add_argument('-r', '--detailedresults', default=False, action='store_true',
                    help = "write a detailed table of individual "
                           "potential annotation statistics and scores to filtering_results.tsv "
                           "in the output directory (default: not written)")
parser.add_argument('-1', '--onefile', type = str, metavar = 'name',
                    help = "if a file name is supplied, write all corrected sequences to a single "
                           "genbank-format file with this name in the output directory (default: "
                           "each input file written to a separate output file)")
parser.add_argument('-k', '--keepalignments', default = False, action = 'store_true',
                    help = "retain all alignments of individual potential annotation to profiles "
                           "(default: individual alignments cleaned up after each assessement)")
parser.add_argument('-n', '--namevariants', type = str, metavar = 'path',
                    help = "path to a file specifying variant gene names to add to those retrieved"
                           "from the list at https://github.com/tjcreedy/genenames (default: no "
                           "additional variants)")  # Additional name variants file, optional
parser.add_argument('-m', '--maxinternalstops', type = int, default = 0, metavar = 'n',
                    help = "number of stops permitted in the translation of annotations (default: "
                           "0")
parser.add_argument('-e', '--framefree', default = False, action  = 'store_true',
                    help = "allow codon search strings to not start in frame with one another; "
                           "not generally recommended (default: codon search strings start in "
                           "frame)")
parser.add_argument('-p', '--potentialfeatures', default = False, action = 'store_true',
                    help = "add all potential annotations as features to the output genbank-format"
                           " file(s); not generally recommended (default: only selected "
                           "annotations recorded)")
parser.add_argument('-f', '--fullassessment', action = 'store_true', default = False,
                    help = "assess the alignment of all potential annotations, even if their "
                           "overlap score means they cannot possibly be selected; not generally "
                           "recommended (default: only align potential annotations that have good "
                           "enough overlap scores to be contenders for selection)")

# Main
if __name__ == "__main__":

    args = parser.parse_args()

    # arglist = ("-s /home/thomas/programming/bioinformatics/mitocorrect/specifications/mitocorrect_specifications_coleoptera_2022-02-24.tsv "
    #            "-g SRAAredux_2022-04-15/annotated_contigs.gb "
    #            "-l testlog.txt "
    #            "-a aaalignfile.tsv "
    #            "-o testout/ "
    #            "-t 2 -b 5 -c aa -r -m 0 -f").split()
    # import os
    # os.chdir('/home/thomas/work/iBioGen_postdoc/MMGdatabase/')
    # args = parser.parse_args(arglist)

    # Parse the arguments into the main utility variables
    utilityvars = mcm.initialise(args)

    # Initialise the queue manager and pool
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(args.threads + 4)

    # Start the writers first in their own threads
    # TODO: if outputing filtering results, add file with selected score of each
    # gene present, plus total average score
    writers, watchers = mcm.start_writers(pool, manager, args)

    # Do the work
    seqrecordgen = mcm.get_seqrecords(args.genbank, args.onefile)
    issues = pool.map(functools.partial(mcm.process_seqrecord, args, utilityvars, writers),
                      seqrecordgen)

    # Delete temporary alignment directory
    if not args.keepalignments:
        shutil.rmtree(utilityvars[4])

    # Close down the writers and the pools
    for w in writers:
        if w is not None:
            w.put(None)
    pool.close()
    pool.join()

    # Process the issues and write to terminal
    mcm.process_issues(issues)

    # Get any errors in any of the writing processes
    for w in watchers:
        w.get()

    exit()
