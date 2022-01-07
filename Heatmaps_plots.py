import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple


def get_matrices(blosum_path: str, our_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Function that reads in the files and returns the matrices as numpy arrays. 
    
    Parameters
    ----------  
    blosum_path: str
        Path to the blosum matrix, that will be returned.
    our_path: str
        Path to our matrix file, that will be returned.

    Returns
    -------
    A Tuple of numpy arrays: list of the one letter code of all amino acids, the blosum matrix values, the values of our
    matrix.
    """
    
    # load in the whole file
    table = np.genfromtxt(blosum_path, delimiter=',', dtype=str)

    amino_acids = table[0]  # extract the list of amino acids
    amino_acids = np.array([aa.strip() for aa in amino_acids])

    blosum_62 = table[1:].astype(int)  # extract the actual values

    # load in our matrix
    table = np.genfromtxt(our_path, delimiter=',', dtype=int)

    our_matrix = table[1:, 1:]  # drop the amino acid labels

    return amino_acids, blosum_62, our_matrix


def plot_direct_heatmaps(blosum_path: str, our_path: str) -> None:
    """
    Function that generates a plot of two heatmaps (one for each matrix) for direct comparison of both.

    Parameters
    ----------
    blosum_path: str
        Path to the Blosum matrix that will be compared.
    our_path: str
        Path to our matrix that will be compared to.

    Returns
    -------
    None
    """
    amino_acids, blosum_62, our_matrix = get_matrices(blosum_path, our_path)

    fig, (ax1, ax2) = plt.subplots(1, 2)

    min_val = min(np.amin(blosum_62), np.amin(our_matrix))
    max_val = max(np.amax(blosum_62), np.amax(our_matrix))

    # generate the heatmap
    hm1 = ax1.imshow(blosum_62, vmin=min_val, vmax=max_val)
    hm2 = ax2.imshow(our_matrix, vmin=min_val, vmax=max_val)

    # set the tick range
    ax1.set_xticks(list(range(20)))
    ax1.set_yticks(list(range(20)))
    ax2.set_xticks(list(range(20)))
    ax2.set_yticks(list(range(20)))

    # label the axes with the amino acids
    ax1.set_xticklabels(amino_acids)
    ax1.set_yticklabels(amino_acids)
    ax2.set_xticklabels(amino_acids)
    ax2.set_yticklabels(amino_acids)

    # set plot names
    ax1.set_title('Blosum62')
    ax2.set_title('Our matrix')

    # add a colorbar
    fig.colorbar(hm1, ax=ax1)
    fig.colorbar(hm2, ax=ax2)

    plt.show()


def main() -> None:
    blosum = 'matrices\\blosum62.csv'
    our = 'matrices\\our_matrix_2.csv'
    plot_direct_heatmaps(blosum, our)
    # plot_difference_plot(blosum, our)


if __name__ == '__main__':
    main()
