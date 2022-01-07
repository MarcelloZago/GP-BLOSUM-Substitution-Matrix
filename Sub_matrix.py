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
        Function to count frequencies of change in amino acids.
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
    def read_directory(path: str, identity_value: int):
        """
        Function to read in multiple preprocessed fasta files in a directory and calculate the frequencies.
        Internally just reads in every file individually and calculates the frequencies.
        Other file types that are in the directory are skipped.

        Parameters
        ----------
        path : str
            Path to a directory that contains the fasta files we want to read in. It can include other files.
        identity_value : int
            The identity threshold for the use of the alignments for the calculation.


        Returns
        -------
        A pandas dataframe with the frequency table
        """

        assert os.path.isdir(path), 'Path is not a directory'  # assure that path is a directory

        # list of amino acids
        aa_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

        # generate frequencies table to sum
        frequency_table = pd.DataFrame(np.zeros((len(aa_list), len(aa_list))), columns=aa_dict, index=aa_dict)

        # iterate over all files in the directory and calculate the added frequencies
        for filename in os.listdir(path):

            # if the file is not a fasta file, skip it
            if not (filename.endswith('.fa') or filename.endswith('.fasta')):
                continue

            # check if the identity can be used, if not skip the file
            if SubstitutionMatrix.__is_usable_identity(identity_value, filename):
                frequency_table = SubstitutionMatrix.count_frequencies_in_file(frequency_table, f'{path}/{filename}')

        return frequency_table

    @staticmethod
    def __is_usable_identity(wanted_identity: int, filename: str) -> bool:
        """
        Private function that checks whether the given file can be used with the given identity threshold.

        Parameters
        ----------
        wanted_identity : int
            The threshold value for the matrix we want to calculate.
        filename : str
            The file name of the file that is currently checked.
        Returns
        -------
        True if the identity of the given alignment is smaller than the given identity threshold.
        """

        split_filename = filename.split('_')
        parsed_identity = int(split_filename[-3])  # the actual identity of the file

        return wanted_identity >= parsed_identity

    @staticmethod
    def sub_matrix(path: str, identity_value: int):
        aa_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                   "M", "F", "P", "S",
                   "T", "W", "Y", "V"]
        table = SubstitutionMatrix.read_directory(path, identity_value=identity_value)
        # Calculate the probability of occurrence for each amino acid pair
        all_pairs = np.triu(table).sum()
        table = table / all_pairs  # q

        # Calculate the probability of occurrence for each amino acid
        table_pi = {"A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "Q": 0, "E": 0, "G": 0, "H": 0, "I": 0, "L": 0, "K": 0,
                    "M": 0, "F": 0, "P": 0, "S": 0,
                    "T": 0, "W": 0, "Y": 0, "V": 0}

        table["sum"] = table.sum(axis=1)

        for i in aa_list:
            table["sum"][i] = table["sum"][i] - table[i][i]

        for i in aa_list:
            table_pi[i] = table[i][i] + (1 / 2) * table["sum"][i]  # p

        # Calculate the expected probability
        exp_table = pd.DataFrame(np.zeros((len(aa_list), len(aa_list))), columns=aa_list, index=aa_list)
        for i in aa_list:
            for j in aa_list:
                if i == j:
                    exp_table[i][j] = table_pi[i] * table_pi[j]
                else:
                    exp_table[i][j] = 2 * table_pi[i] * table_pi[j]

        # Calculate the log ratio

        sub_matrix = pd.DataFrame(np.zeros((len(aa_list), len(aa_list))), columns=aa_list, index=aa_list)

        for i in aa_list:
            for j in aa_list:
                sub_matrix[i][j] = round(2 * math.log((table[i][j] / exp_table[i][j]), 2))

        return sub_matrix


if __name__ == '__main__':

    """
    Test 1
    """

    # print(SubstitutionMatrix.count_aa_in_seq("AGGGAAAYYY"))

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

    # print(out_matrix)

    """
    Test 3: Testing SubstitutionMatrix.read_directory with "processed" directory
    """

    path = "processed\\"

    # table= SubstitutionMatrix.read_directory(path)

    # print(table)

    """
    Calculate the substitution matrix
    """
    identity_threshold = 62

    sub_matrix = SubstitutionMatrix.sub_matrix(path, identity_value=identity_threshold)
    print(sub_matrix.to_string)

    sub_matrix.to_csv(f'matrices\\our_matrix_{identity_threshold}.csv')
