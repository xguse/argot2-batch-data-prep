"""Functions used in this pipeline."""

import pandas as pd
import numpy as np

import subprocess as sp

def get_total_record_count(fasta_file):
    """Count number of headers found in fasta file."""
    headers = sp.run(['grep', '>', fasta_file], stdout=sp.PIPE).stdout.split(b'\n')
    
    return len(headers)


def get_number_of_groups(total, group_size):
    """Return number of groups to split the original fasta file into."""
    return int(total/group_size)+1
