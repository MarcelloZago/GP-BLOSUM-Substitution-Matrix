# Group Project BLOSUM Substitution Matrix
Github for the Group Project about BLOSUM Substitution Matrices


## Project structure

The unprocessed data that is used to calculate the matrix needs to be stored in a directory called ```Data```. Here can subdirectories of aquired data be stored. Each subdirectory needs to be processed manually. All processed files will be stored in the ```processed``` directory.

The ```hand-in``` directory contains all additional data for the hand-in of the project.

**Example project structure:**
```
.
├── Data
│   ├── subdirectory_1
│   │     ├── alignment_1.fasta
|   |     ...
│   │     └── alignment_k.fasta
│   ├── subdirectory_2
│   │     ├── alignment_1.fasta
|   |     ...
│   │     └── alignment_k.fasta
|   ...
│   └── subdirectory_m
│         ├── alignment_1.fasta
|         ...
│         └── alignment_k.fasta
|   
├── hand-in
|   ...
|   
├── processed
│   ├── subdirectory_1_63_processed_1.fasta
│   ├── subdirectory_1_67_processed_2.fasta
|   ...
│   └── subdirectory_m_99_processed_k.fasta
|
├── Fasta.py
├── Heatmap_plots.py
├── README.md
└── Sub_matrix.py
```


## Naming convention of processed files

We will save the preprocessed alignments in the ```processed``` directory. Here all the preprocessed ```.fasta``` files can be stored locally. The calculated identities are "stored" in the filename for easy handling for the future calculation of the matrices. Because of the directory wise processing, the subdirectory name is added to the filename as well as a simple id for this subdirectory. Here is a depiction of the naming convention:


```[subdirectory name]_[identity]_[filename]_[id].fasta ```


## Short Tutorial

### Process new data
First, a new directory containing all ```.fasta``` files needs to be added to the ```Data``` directory. These files can now be processed by adding following line to the main function of ```Fasta.py```:

```AlignmentFasta.preprocess_directory('Data/new_data_directory')```

After running the main function, all files are automatically processed, the identies calculated and saved to the ```processed``` directory.

### Calculate a matrix with all processed files
Change the ```identity_threshold``` to a wanted percentage $0-100\%$ in the main function of ```Sub_matrix.py```. Afterwards, the file just needs to be runned and the new matrix will be calculated and saved as a ```.csv``` file.

