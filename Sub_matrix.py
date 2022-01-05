import os
import numpy as np
import pandas as pd
from Bio import AlignIO
import math


class SubstitutionMatrix:

    @staticmethod
    def count_aa_in_seq(sequence: str):
        """
        Function to count amino acids in a sequence
        :param sequence:
        :return:

        Just a function to try out
        """
        aa_count = {"A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "Q": 0, "E": 0, "G": 0, "H": 0, "I": 0, "L": 0, "K": 0,
                    "M": 0, "F": 0, "P": 0, "S": 0,
                    "T": 0, "W": 0, "Y": 0, "V": 0, "Total": 0}

        for pos in sequence:
            aa_count[pos] += 1
            aa_count["Total"] += 1

        return aa_count

    @staticmethod
    def count_frequencies_in_file(frequency_table: pd.DataFrame, path: str):
        """
        Function to count frequencies of change in aminoacids.
        Compares all pairs of possible sequences in a file.
        Sums into a pd table of amino acid change.

        Parameters
        ----------
        path : str
            Path to the fasta file we want to read in.
        frequency_table : pd.DataFrame

        Returns
        -------
        A pandas dataframe with the frequency table

        """
        assert path.endswith('.fa') or path.endswith('.fasta'), 'Cannot read a non Fasta file.'

        with open(path, 'r') as file:
            # use of Biopython's implementation, because it is more efficient
            alignment = AlignIO.read(file, 'fasta')

            alignment_length = len(alignment)

            # iterate over each line of the alignment object and add it to our list
            for record_i in range(alignment_length):

                for record_j in range(alignment_length):

                    if record_i != record_j:

                        sequence_length = len(alignment[record_i].seq)

                        for pos_i in range(sequence_length):

                            letter_first_seq = str(alignment[record_i].seq[pos_i])
                            letter_second_seq = str(alignment[record_j].seq[pos_i])

                            frequency_table[letter_first_seq][letter_second_seq] += 1

        return frequency_table
    @staticmethod
    def read_directory(path: str):
        """
        Function to read in multiple preprocessed fasta files in a directory and calculate the frequencies.
        Internally just reads in every file individually and calculates the frequencies.
        Other file types that are in the directory are skipped.

        Parameters
        ----------
        path : str
            Path to a directory that contains the fasta files we want to read in. It can include other files.

        Returns
        -------
        A pandas dataframe with the frequency table
        """

        assert os.path.isdir(path), 'Path is not a directory'  # assure that path is a directory

        # list of aminoacids
        aa_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

        # generate frequencies table to sum
        frequency_table = pd.DataFrame(np.zeros((len(aa_list), len(aa_list))), columns=aa_dict, index=aa_dict)


        # iterate over all files in the directory and calculate the added frequencies
        for filename in os.listdir(path):

            # if the file is not a fasta file, skip it
            if not (filename.endswith('.fa') or filename.endswith('.fasta')):
                continue

            frequency_table = SubstitutionMatrix.count_frequencies_in_file(frequency_table, f'{path}/{filename}')


        return frequency_table

    @staticmethod
    def sub_matrix(path: str):
        aa_dict = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                   "M", "F", "P", "S",
                   "T", "W", "Y", "V"]
        table = SubstitutionMatrix.read_directory(path)
        # Calculate the probability of occurence for each amino acid pair
        all_pairs = np.triu(table).sum()
        table = table / all_pairs

        # Calculate the probability of occurence for each amino acid
        table_pi = {"A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "Q": 0, "E": 0, "G": 0, "H": 0, "I": 0, "L": 0, "K": 0,
                    "M": 0, "F": 0, "P": 0, "S": 0,
                    "T": 0, "W": 0, "Y": 0, "V": 0}

        table["sum"] = table.sum(axis=1)

        for i in aa_dict:
            table["sum"][i] = table["sum"][i] - table[i][i]

        for i in aa_dict:
            table_pi[i] = table[i][i] + (1 / 2) * table["sum"][i]

        # Calculate the expected probability
        exp_table = pd.DataFrame(np.zeros((len(aa_dict), len(aa_dict))), columns=aa_dict, index=aa_dict)
        for i in aa_dict:
            for j in aa_dict:
                if i == j:
                    exp_table[i][j] = table_pi[i] * table_pi[j]
                else:
                    exp_table[i][j] = 2 * table_pi[i] * table_pi[j]

        # Calculate the log ratio

        sub_matrix = pd.DataFrame(np.zeros((len(aa_dict), len(aa_dict))), columns=aa_dict, index=aa_dict)
        for i in aa_dict:
            for j in aa_dict:
                sub_matrix[i][j] = int(math.log((table[i][j] / exp_table[i][j]), 2))

        return sub_matrix



if __name__ == '__main__':

    """
    Test 1
    """

    #print(SubstitutionMatrix.count_aa_in_seq("AGGGAAAYYY"))

    """
    Test 2: Testing SubstitutionMatrix.count_frequencies_in_file with "_7_processed_0.fasta"
    """

    aa_dict = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
               "M", "F", "P", "S",
               "T", "W", "Y", "V"]

    # frequency_table = pd.DataFrame(np.zeros((len(aa_dict), len(aa_dict))), columns=aa_dict, index=aa_dict)
    #
    # path = "processed/Proteobacteria_priest_2021_0_processed_7.fasta"
    #
    # out_matrix = SubstitutionMatrix.count_frequencies_in_file(frequency_table, path)

    #print(out_matrix)

    """
    Test 3: Testing SubstitutionMatrix.read_directory with "processed" directory
    """

    path = "processed/"

    #table= SubstitutionMatrix.read_directory(path)

    #print(table)


    """
    Calculate the substitution matrix
    """
    sub_matrix = SubstitutionMatrix.sub_matrix(path)
    print(sub_matrix)