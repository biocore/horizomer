### Benchmark for HGT detection software

## Scripts details

	* `launch_software.sh` launch all tested software
	* `reformat_input.py` reformats input data (ex. species
		and gene trees in Newick format) to input supported
		by each software (ex. Nexus format)
	* `parse_output.py` parses the output of each software
		to report the number of HGTs detected
	* `compute_accuracy.py` uses the output of `launch_software.sh`
		to compute the number of true positive, false positive
		and false negative HGTs, as well as precision, recall
		and F-score for each software
	* `test_1_simulate_genomes.sh` simulates various data sets
		using the [ALF simulator](http://mbe.oxfordjournals.org/content/early/2011/12/07/molbev.msr268)

## Running benchmark

See [INSTALL.md](https://github.com/biocore/WGS-HGT/blob/master/benchmark/INSTALL.md)
for the dependencies list.

## Example

```
bash launch_software.sh /path/to/output/folder /path/to/WGS-HGT/benchmark /path/to/WGS-HGT/benchmark/test_benchmark/RealTree.nwk None None None None /path/to/WGS-HGT/benchmark/test_benchmark/GeneTrees/ /path/to/WGS-HGT/benchmark/test_benchmark/MSA/ /path/to/phylonet/install/dir /path/to/jane4/install/dir /path/to/trex/install/dir false
```

## Additional information

This benchmark was setup to allow easy integration of other
software. The user must include a system call for the tool
in `launch_software.sh` and input/output parsing scripts
(if applicable) to `reformat_input.py` and `parse_output.py`.