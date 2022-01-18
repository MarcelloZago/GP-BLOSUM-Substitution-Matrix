import matplotlib.cm
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


def plot_direct_heatmaps(blosum_path: str, our_path: str, identity: int) -> None:
    """
    Function that generates a plot of two heatmaps (one for each matrix) for direct comparison of both.

    Parameters
    ----------
    blosum_path: str
        Path to the Blosum matrix that will be compared.
    our_path: str
        Path to our matrix that will be compared to.
    identity: int
        Identity that was used for the matrices. It is just used for the individual Heatmap titles.

    Returns
    -------
    None
    """
    amino_acids, blosum_62, our_matrix = get_matrices(blosum_path, our_path)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    min_val = min(np.amin(blosum_62), np.amin(our_matrix))  # minimal value for both matrices
    max_val = max(np.amax(blosum_62), np.amax(our_matrix))  # maximal value for both matrices

    max_max = max(abs(min_val), abs(max_val))  # absolute max, for symmetrical colormap

    # generate the heatmap
    cmap = matplotlib.cm.get_cmap('Spectral_r')
    hm1 = ax1.imshow(blosum_62, vmin=-max_max, vmax=max_max, cmap=cmap)
    hm2 = ax2.imshow(our_matrix, vmin=-max_max, vmax=max_max, cmap=cmap)

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
    ax1.set_title(f'Blosum{identity}')
    ax2.set_title(f'Our matrix ({identity})')

    # add a colorbar
    fig.colorbar(hm1, ax=ax1)
    fig.colorbar(hm2, ax=ax2)

    plt.show()


def plot_difference_plot(blosum_path: str, our_path: str) -> None:
    """
    Function that generates a heatmap of the absolute difference (blosum - our_matrix) of the two matrices.

    Parameters
    ----------
    blosum_path: str
        Path for the Blosum matrix.
    our_path: str
        Path for our matrix.

    Returns
    -------
    None
    """
    amino_acids, blosum, our_matrix = get_matrices(blosum_path, our_path)

    difference = np.abs(blosum - our_matrix)

    plt.imshow(difference)

    plt.colorbar()

    plt.show()


def main() -> None:

    for identity in [62, 90]:
        blosum = f'matrices\\blosum{identity}.csv'
        our = f'matrices\\our_matrix_{identity}.csv'

        plot_direct_heatmaps(blosum, our, identity)

    # Not used
    # plot_difference_plot(blosum, our)


if __name__ == '__main__':
    main()
