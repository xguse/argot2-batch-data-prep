"""Split up main fasta file into many to facilitate parallelized DB searches."""

import pyfaidx

# params
fasta_header_delim = snakemake.params.fasta_header_delim

# input
fasta_path = snakemake.input.fasta_path

#output
split_fastas = snakemake.output.split_fastas



def write_group_to_file(name_group, recs, path):
    """Write the formated fasta recs for each rec in `name_group` to file with path `path`."""
    with open(path, 'w') as out:
        for name in name_group:
            out.write(rec_to_string(recs, name))

def rec_to_string(recs, name):
    """Format a rec into correct string format for writing to file."""
    rec_name = name
    rec_seq = recs[name][:].seq
    return ">{rec_name}\n{rec_seq}\n".format(rec_name=rec_name,rec_seq=rec_seq)


# Load the main fasta file
orig_records = pyfaidx.Fasta(fasta_path, key_function=lambda key: key.split(fasta_header_delim)[0])


# set important grouping variables
seq_names = list(orig_records.keys())

num_of_seqs = len(seq_names)
num_of_groups = len(split_fastas)

num_per_group = int(num_of_seqs/num_of_groups)+1


# collect seq_names into their groups
name_groups = [seq_names[i:i+num_per_group] for i in range(0, num_of_seqs, num_per_group)]

if len(split_fastas) != len(name_groups):
    raise ValueError("len(split_fastas) != len(name_groups)")

# write out the groups to their split_fastas files
for group, out_path in zip(name_groups, split_fastas):
    write_group_to_file(name_group=group, recs=orig_records, path=out_path)
