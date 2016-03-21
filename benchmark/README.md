## Benchmark for HGT detection software

### Scripts details

1. *launch_software.sh*
   launch all tested software

2. *reformat_input.py*
   reformats input data (ex. species and gene trees in Newick format)
   to input supported by each software (ex. Nexus format)

3. *parse_output.py*
   parses the output of each software to report the number of HGTs detected

4. *compute_accuracy.py*
   uses the output of ```launch_software.sh``` to compute the number of
   true positive, false positive and false negative HGTs, as well as
   precision, recall and F-score for each software

5. *test_1_simulate_genomes.sh*
   simulates various data sets using the [ALF simulator](http://mbe.oxfordjournals.org/content/early/2011/12/07/molbev.msr268)

6. *run_*`.sh*
   contain commands specific to running each HGT detection tool

### Running benchmark

See [INSTALL.md](https://github.com/biocore/WGS-HGT/blob/master/benchmark/INSTALL.md)
for the dependencies list.

### Additional information

This benchmark was setup to allow easy integration of other
software. The user must include a system call for the tool
in *launch_software.sh* and input/output parsing scripts
(if applicable) to *reformat_input.py* and *parse_output.py*.