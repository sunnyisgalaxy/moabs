# MOABS v1.3.9.5 @ 20200626

Feature release:

1. To support chromosome split when generating statistics file in `MCALL`.

2. To Refine the error message of short query length in `BSeQC2`.

# MOABS v1.3.9.4 @ 20200603

Minor release: add notes about usage of `-i` and `--def`.

# MOABS v1.3.9.3 @ 20200526

Compiling fixes on Mac OS X.

1. To rename version to version.txt under gzstream, to avoid name clash with STL by clang.

2. To suppress the diagnosis on narrow conversions in BSMAP using `-Wno-narrowing` under OSX. Otherwise, clang fails the compilation due to narrow conversions.

# MOABS v1.3.9.2 @ 20200524

A feature improvement to speedup numerical integral for big counts in pdiff() using Gauss-Kronrod quadrature.

# MOABS v1.3.9.1 @ 20200517

A new feature for BSeQC2 and PEFilter: add a traditional construction protocol by Nextera transposase-based approach.

# MOABS v1.3.9.0 @ 20200506

New feature: add a second quality control for the PEFilter BAM.

# MOABS v1.3.8.9 @ 20200505

Bug fix: enable quality control for single-end mapping in BSeQC2.

# MOABS v1.3.8.8 @ 20200502

Feature release and bug fixes.

1. Support a reference FASTA header with additional descriptions in MCALL and BSeQC2.

2. Remove `Rplots.pdf` in BSeQC methyaltion bias plots.

3. Refine the output header in statistics file of the `-f|--predefinedFeature` option in MCOMP.

# MOABS v1.3.8.7 @ 20200207

Feature release and bug fixes.

1. Multiple `-A` is implemented to BSMAP in MOABS configuration file.

2. Bug fix for the sorted FASTQ files to BSMAP.

3. Bug fix for creating directory in BSeQC2.

# MOABS v1.3.8.5 @ 20200122

Feature releases and minor bug fixes.

1. Add LiBis to support low input bisufilte experiments.

2. Add an option `-U` to `BSMAP` to turn off the automatic BAM sort.

3. Bug fix for `MCALL` with out-of-range local seq.

4. Fix warning message report for `bseqc2`.

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

