import re
import sys

#to run this code use python zee_a2.py testfas.fas enzymes.txt
def run(seq_file, enzyme_file):
	"""Runs our code

	Arguments/Args:
		seq_file: FASTA file
		enzyme_file: Enzyme file
	"""
	sequence_read_lines = read_files(seq_file, enzyme_file)[0]
	recog_site_read = read_files(seq_file, enzyme_file)[1]
	recog_site_read_lines = read_files(seq_file, enzyme_file)[2]
	print("Restriction enzyme analysis of sequence from file " + seq_file)
	print("Cutting with enzymes found in file " + enzyme_file)
	create_report(sequence_read_lines, recog_site_read_lines)


def read_files(seq, recog_site):
        """Reads provided files

        Args:
                seq: sequence/FASTA file
                recog_site: enzyme file

        Returns:
                f_readlines: list of lines from first FASTA file
                e_readlines: list of lines from enzyme file
                e_read: string of lines from enzyme file

        """
        with open(seq, "r") as f, open(recog_site, "r") as e:
                f_readlines = f.readlines()
                e_read = e.read()
                e.seek(0) #returns seek to 0 because it moves to end of file after each read
                e_readlines = e.readlines()

        return f_readlines, e_read, e_readlines


def get_cutting_sites(seq, cut_enzyme):
	"""Get the cutting site position for the provided sequence.

	Args:
		seq: The provided sequence.
		cut_enzyme: Where to cut the sequence.

	Returns:
		sites: list of the sites to cut.
	"""
	before_cut = cut_enzyme.split("^")[0]

	sites = []
	for site in re.finditer(cut_enzyme.replace("^", ""), seq):
		sites.append(site.start() + len(before_cut))
	return sites


def cut_sequence(sequence, cutting_sites):
	"""Cuts the provided sequence according to cutting sites

	Args:
		sequence: The provided sequence
		cutting_sites: Enzyme cutting sites of the sequence

	Returns:
		framents_list: A list of fragments after cutting the sequence

	"""
	fragments_list = []
	previous_site = 0
	for site in range(0, len(cutting_sites)):
		fragment = sequence[previous_site:cutting_sites[site]]
		fragments_list.append(fragment)
		previous_site = cutting_sites[site]
	fragments_list.append(sequence[previous_site:len(sequence)])

	return fragments_list


def split_sequence(sequence, val):
	"""Splits the sequence according to provided value

	Args:
		sequence: The sequence to be splitted.
		val: Value of where to split the sequence.
	"""
	splitted_sequence = []
	for i in range(0, len(sequence), val):
		splitted_sequence.append(sequence[i:i+val])
	return " ".join(splitted_sequence), splitted_sequence


def create_report(sequence_read_lines, recog_site_read_lines):
	"""Creates the report

	Args:
		sequence_read_lines: list of lines from the sequence file
		recog_site_read_lines: list of lines from the enzyme file
	"""
	line_length = 63 #an arbitrary amount of lines added
	full_sequence = sequence_read_lines[1].strip()
	total_sequence_length = len(full_sequence)
	print("-"*line_length)
	print("Sequence name: " + sequence_read_lines[0].strip().split(">")[1])
	print("Sequence is " + str(total_sequence_length) + " bases long")

	for line in recog_site_read_lines:
		print("-"*line_length)
		enzyme_name = line.strip().split(";")[0] 
		cut_enzyme = line.strip().split(";")[1] 
		cutting_sites = get_cutting_sites(full_sequence, cut_enzyme)

		if not cutting_sites:
			print("There are no sites for " + enzyme_name + "\n")
			continue #exits the forloop if there are no cutting sites but proceeds if there are
		print("There are " + str(len(cutting_sites)) + " cutting sites for " + enzyme_name + ", cutting at " + cut_enzyme)
		fragments = cut_sequence(full_sequence, cutting_sites)
		print("There are " + str(len(fragments)) + " fragments:\n")

		LENGTH = "Length:" #not really needed but helps to format like the example
		fragment_start = 1
		for fragment in fragments:
			print(LENGTH + str(len(fragment)))

			if len(fragment) > 60: #splits frgaments by 60 if longer than 60 bases
				fragments = split_sequence(fragment, 60)[1]
			else:
				fragments = [fragment]

			for split_fragment in fragments:
				split_bases = split_sequence(split_fragment, 10)[0] #grouping bases of 10
				print(str(fragment_start) + " "*(len(LENGTH) - len(str(fragment_start))) + (split_bases))
				fragment_start += len(split_fragment)
		print("")


fasta_f = sys.argv[1]
enzyme_f = sys.argv[2]
run(fasta_f, enzyme_f)
