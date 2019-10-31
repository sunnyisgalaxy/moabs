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

