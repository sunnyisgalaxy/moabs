# MOABS v1.3.8.2 @ 20200108

## Two newly added features

1. bseqc2

`bseqc2` performs an automatic QC for a sorted BAM file. The BAM file can be
single-end mapping or paired-end mapping. Failed samples will be skipped for
downstream analysis. This feature is suitable for a quick QC.

2. pefilter

`pefilter` aims to filter invalid mappings in a paired-end BAM file. Mapping
errors may introduce false estimation of methylation levels. True paired-end
mappings adhere to the library construction protocol, i.e., traditional library
and Pico library. This new feature will detect the library protocol and rule
out invalid PE mappings.

# MOABS v1.3.7.8 @ 20191031

## Feature improvement for `--withVariance` in mcomp

- To speedup merge by skipping existing chromosome files

# MOABS v1.3.7.7 @ 20191029

## Changes since v1.3.4

1. The source codes are organized using Autotools.

2. Embeded Boost library has been moved from the package, and the system Boost will be linked instead.

3. Feature fix for `--withVariance` in mcomp.

    3.1 The two-way pipe decouples the Beta-binomial fitting from the main program.

    3.2 Implement the scalabe running by chromosomes.

    3.3 Speedup for consensus replicates without variances.

