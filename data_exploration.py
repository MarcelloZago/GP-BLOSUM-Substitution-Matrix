import os

from Bio import AlignIO


def print_all_alignments(path: str) -> None:
    print('Lengths of all ebola sequences: \n')

    for filename in os.listdir(path):
        assert filename.endswith('.fa')

        with open(f'{path}\\{filename}', 'r') as file:
            alignment = AlignIO.read(file, 'fasta')

            print(f'{filename:<12}: {alignment.get_alignment_length():<4} aa, {len(alignment)} sequences')


def print_al():
    with open('Data\\L.fa') as file:
        al = AlignIO.read(file, 'fasta')
        i = 0

        for read in al:
            i += 1
            print(read.seq)
            if i == 100:
                break

    return

def main() -> None:
    print_all_alignments('Data')
    #print_al()


if __name__ == '__main__':
    main()
