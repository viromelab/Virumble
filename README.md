# Virumble

Reference-free assembly of viral sequences using a compression-based methods.

## Install Virumble ##

Then, clone the repository and compile the code using the following instructions:
```
git clone https://github.com/viromelab/Virumble
cd Virumble/
make clean
make
```

### PARAMETERS ###

To see the possible reconstruction options type
<pre>
./virumble -h
</pre>

This will print the following options:
```

Usage: ./virumble -f <forward.fastq> -r <reverse.fastq> -c <context_size> [options]

Options:
  -h, --help              Show this help
  -f, --forward FILE      Forward FASTQ file (required)
  -r, --reverse FILE      Reverse FASTQ file (required)
  -o, --output FILE       Output TSV file (default: output.tsv)
  -t, --threads N         Number of threads (default: 1)
  -c, --context N         Context size (kmer length, >=2)
  -s, --substitutions N   Max substitutions for merging (default: 3)
  -e, --overlap N         Minimum overlap for merging contigs (default: 30)
  -v, --verbose           Verbose output

```


### CITATION ###

On using this software/method please cite:

* Pending

### ISSUES ###

For any issue let us know at [issues link](https://github.com/viromelab/Virumble/issues).

### LICENSE ###

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>
