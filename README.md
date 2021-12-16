# Group Project BLOSUM Substitution Matrix
Github for the Group Project about BLOSUM Substitution Matrices


## Coding outline

This section is a collection of all functions and tasks that need to be implemented. Please add your name in brackets so that the others know that the respective function/ task is already worked on.

- Deletion of gaps in the alignment (Preprocessing) [Marcello]
- Calculation of identicality []
- Counting of frequencies []
- Calculation of the actual log-odd-scores (the matrix) []

## Project structure and naming convention

We will save the preprocessed alignments in a directory ```processed```. Here all the preprocessed ```.fasta``` files can be stored locally. The identicality of the alignment, when computed will be appended as a prefix to the file name in following manner:

```__[ident]__[filename].fasta ```

It should be sufficient to store the indenticality as an integer percentage values (in between 0-100), but that can be discussed.

Some examples:
- ```__99__e_coli_genome.fasta ```
- ```__62__shortreadfile.fa ```
- ```__1__averyobscurealignment.fasta ```
