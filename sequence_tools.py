import re
from Bio.Seq import Seq
from Bio import SeqIO
import tempfile
import os
import csv
#=====================================================================

def validate_seq(sequence):
    """Validate DNA sequence without header."""
    sequence = sequence.strip()
    sequence = sequence.replace(" ", "")
    sequence.upper()
    regex = re.compile('^[ACTGNRYSWKMBDHVEFILPQSXZ]*$', re.I)
    if regex.search(sequence) is not None:
        return True
    else:
        return False


def validate_single_fasta_seq(sequence):
    """Validate sequence in single Fasta format."""
    sequence = sequence.replace(" ", "")
    sequence.upper()
    regex = re.compile('>\S*\n[ACTGNRYSWKMBDHVEFILPQSXZ]', re.MULTILINE)
    if regex.search(sequence) is not None:
        return True
    else:
        return False


def validate_fasta_seq(sequence):
    """Validate sequence in multiple Fasta format."""
    sequence = sequence.replace(" ", "")
    sequence.upper()
    regex = re.compile('>\S*\n[ACTGNRYSWKMBDHVEFILPQSXZ]*', re.MULTILINE)
    if regex.search(sequence) is not None:
        return True
    else:
        return False


def validate_csv_seq(sequence):
    """Validate sequence in csv format. Only comma or semicolon separated """
    if sequence.find(',') != -1 or sequence.find(';') != -1:
        return True
    else:
        return False


def reverse_complement(sequence_file_location, rc):
    """Reverse_complement a sequence."""
    f_fasta = open(sequence_file_location, "r")
    temp_seq_file = tempfile.mkstemp()
    f_seq = open(temp_seq_file[1] + '.txt', 'w')
    for seq_record in SeqIO.parse(f_fasta, "fasta"):
        id = seq_record.id
        seq = str(seq_record.seq)
        if rc:
            sequence_reverse = str(Seq(seq).reverse_complement())
        else:
            sequence_reverse = seq[::-1]
        f_seq.write('>' + str(id) + '\n')
        f_seq.write(sequence_reverse)
        f_seq.write('\n')
        sequence_file_location = temp_seq_file[1] + '.txt'
    f_fasta.close()
    f_seq.close()
    return sequence_file_location


def id_extract(ids, db_location, database_file):
    """Finds substring in fasta id and extract sequences."""
    # For exact searching, need to be adapted for ignored letter
    # SeqIO.write((seq for seq in seqiter if seq.id in ids), temp_seq_file[1]+'.txt', "fasta")

    os.chdir(db_location)
    temp_seq_file = tempfile.mkstemp()
    f_results = open(temp_seq_file[1] + '.txt', 'w')

    seqiter = SeqIO.parse(open(database_file), 'fasta')

    for id_to_find in ids:
        for sequence in seqiter:
            if sequence.id.find(id_to_find) is not -1:
                f_results.write('>' + sequence.id + '\n' + str(sequence.seq) + '\n')
    f_results.close()

    return temp_seq_file[1] + '.txt'


def set_blast_settings(blast_type):
    # Settings:
    # evalue, min_percent, min_length, wordsize, strand, processor, genetic_code, matrix, gap_open, gap_ext, ungapped

    if blast_type == 'blastn':
        settings = ['1e-10', '80', '100', '11', 0, '1', 'enabled', 'enabled', '5', '2', 2]
    elif blast_type == 'blastp':
        settings = ['1e-10', 'enabled', '100', '3', 'enabled', '1', 'enabled', 4, '9', '1', 'enabled']
    elif blast_type == 'blastx':
        settings = ['1e-10', 'enabled', '100', '3', 0, '1', 0, 4, '11', '1', 'enabled']
    elif blast_type == 'tblastn':
        settings = ['1e-10', 'enabled', '100', '3', 'enabled', '1', 'enabled', 4, '11', '1', 'enabled']
    elif blast_type == 'tblastx':
        settings = ['1e-10', 'enabled', '100', '3', 0, '1', 0, 4, 'enabled', 'enabled', 'enabled']

    return settings


def add_header(in_file, file_type):
    """Adds a header a first line of outfile and return new file."""

    if file_type == 'bowtie':
        header = "Read name\t" + "Reference strand\t" + "Name of reference sequence\t" \
                 + "Position alignment occurs\t" + "Read sequence\t" + "Read qualities\t" \
                 + "Ceiling\t" + "Mismatch descriptors\n"
    else:
        header = ''

    # Temp file for final results including header
    temp_out = tempfile.mkstemp()
    f_in = open(in_file, 'r')
    results = f_in.read()
    f_out = open(temp_out[1] + '.txt', 'w')
    f_out.write(header)
    f_out.write(results)

    f_in.close()
    f_out.close()
    return temp_out[1] + '.txt'


def convert_seq_fomats(in_file, conversion_type):
    """Converts sequence file format from fasta to excel and vice versa. Return converted file."""
    temp_out = tempfile.mkstemp()
    if conversion_type == "fasta_to_excel":
        seq_iter = SeqIO.parse(open(in_file), 'fasta')
        writer = csv.writer(open(temp_out[1] + ".csv", "wb"))
        for seq in seq_iter:
            writer.writerow(["Sequence ID", "Sequence"])
            writer.writerow([seq.id, str(seq.seq)])
        return temp_out[1] + ".csv"
    else:
        try:
            f_in = open(in_file, "rb")
            # Check weather csv is semicolon or comma separated.
            dialect = csv.Sniffer().sniff(f_in.read(1024), delimiters=";,")
            f_in.seek(0)
            reader = csv.reader(f_in, dialect)
            f_out = open(temp_out[1] + '.txt', 'w')
            for row in reader:
                print row
                f_out.write('>' + row[0:2][0] + '\n' + row[0:2][1] + '\n')
            f_out.close()
            return temp_out[1] + '.txt'
        except IndexError:
            return None
