"""Translate DNA fastas into all six frames concatonated into single protein seq for each DNA record.

Usage:
    translate_fastas.py DNA_PATH OUTPUT_PATH

Arguments:
    DNA_PATH        path to DNA fasta file
    OUTPUT_PATH     path for fasta to be created

Options:
  -h --help     Show this screen.
  --version     Show version.
"""

from docopt import docopt

import pyfaidx

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


__version__ = "0.0.1"


# # params
#
# # input
# split_fastas = snakemake.input.split_fastas
#
# # output
# translated_fastas = snakemake.output.translated_fastas



def translate_6frames_and_concat(seq):
    """Return all 6 translated frames of `seq` concatonated by 50 Xs.

    Expects input to be a pyfaidx obj.
    Returns a biopython.Seq object.
    """
    seq_str = seq[:].seq

    seq_str = pad_to_multiple_of_3(seq_str)

    dna = Seq(seq_str, generic_dna)
    dna_rc = dna.reverse_complement()

    Xs = "X"*50

    prot_fwd = dna.translate(stop_symbol="X") + Xs + (dna[1:] + 'N').translate(stop_symbol="X") + Xs + (dna[2:] + 'NN').translate(stop_symbol="X")

    prot_rvc = dna_rc.translate(stop_symbol="X") + Xs + (dna_rc[1:] + 'N').translate(stop_symbol="X") + Xs + (dna_rc[2:] + 'NN').translate(stop_symbol="X")

    return prot_fwd + Xs + prot_rvc

def pad_to_multiple_of_3(seq_str):
    """Return `seq_str` with enough Ns added to end to be a multiple of three."""
    olen = len(seq_str)
    left_over = olen % 3
    to_add = 3 - left_over

    if to_add == 3:
        to_add = 0

    return seq_str + 'N'*to_add


def do_translation(dna_path, output_path):
    """Translate DNA fastas into all six frames concatonated into single protein seq for each DNA record."""
    dfasta = pyfaidx.Fasta(dna_path)

    with open(output_path,'w') as out:
        for name in dfasta.keys():
            peptid = translate_6frames_and_concat(seq=dfasta[name])
            rec_str = ">{name}\n{seq_str}\n"
            out.write(rec_str.format(name=name,seq_str=str(peptid)))

if __name__ == "__main__":
    args = docopt(__doc__, version=__version__)

    do_translation(dna_path=args["DNA_PATH"], output_path=args["OUTPUT_PATH"])
