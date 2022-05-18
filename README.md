#kmer-map

This program can only be used on Linux systems

### Installation

To install Noise kmermap using the source:  
1. Dowload:
```bash  
     git clone --https://github.com/xianjia10/kmer-map.git
```  

2. Compile:  
```bash  
    cd kmer-map/src  
    make  
```
3. You may wish to add the path to your copy of the kmermap
directory to your PATH variable.


### Prerequisites

* gcc or similar C compiler and linker


### Usage Overview

The simplest use, searching reads.fa for the repeated motif GGCAT:

1.build
```bash 
kmermap [-t <threads>] [-o <out_path>] build [-q <query>] [-r <reffile>] <kmerfile>
```

2.generate kmer map
```bash
kmermap genmap <r_pos> <q_pos> <paf>
```

### Usage Details

```bash  
    -t, --threads         Number of threads
    -o, --out-put         Diractory of output
    build:
        -q, --query       Query/Reads fasta file
        -r, --reference   Reference fasta file
        <kmerfile>        Kmer file with creating by jellyfish

    genmap:
        <r_pos>           Ref pos file built with build
        <q_pos>           Query pos file built with build
        <paf>             Alignment file

    -h, --help            Show this page
    -V, --version         Version
```
### Output file description
The *.sort.paf is the result of *.paf sorted by reference name and query name.
The *.sort.kmermap file contains the final result.