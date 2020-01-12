# MOABS: MOdel based Analysis of Bisulfite Sequencing data
A comprehensive, accurate and efficient solution for analysis of large scale base-resolution DNA methylation data, bisulfite sequencing or single molecule direct sequencing.

## MOABS at Bioconda

MOABS has been deployed in https://anaconda.org/bioconda/moabs. To install MOABS from Bioconda, the following channels should be added. Namely,

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install moabs
```

## MOABS at UseGalaxy.eu

MOABS has been depolyed in https://usegalaxy.eu/root?tool_id=moabs for public usage.

## To use LiBis

LiBis aims to rescue unmapped reads in low input bisulfite experiments. See more at https://github.com/Dangertrip/LiBis

In order to invoke LiBis as a module in MOABS, LiBis should be installed separately. I.e.,

```
conda install libis
```

## Abstract

MOABS seamlessly integrates alignment, methylation calling, identification of hypomethylation for one sample and differential methylation for multiple samples, and other downstream analysis.

## Comprehensive Document

The PDF version is available at

http://dldcc-web.brc.bcm.edu/lilab/deqiangs/moabs/moabs-v1.2.2.pdf

The online HTML version created from 'latex2html' command is available at

http://dldcc-web.brc.bcm.edu/lilab/deqiangs/moabs/moabs.html

## Download

The source code, prebuilt binary on x86_64 Linux system, and test data

https://s3.amazonaws.com/deqiangsun/software/moabs-v1.3.0.src.x86_64_Linux.data.tar.gz

If you download the source from github: The data files in folder bin/ have been deleted in github repository due to size limit. Please make sure you download them and put them in bin/ folder.

## Installation of prebuilt binaries

Download moabs-v1.3.0.src.x86_64_Linux.data.tar.gz to /your/path/, Add the bin/ to your $PATH variable by command export PATH=/your/path/moabs-v1.3.0.src.x86_64_Linux.data/bin/:$PATH, Then use the prebuilt executables in bin/ ! Note the moabs path need be in front of $PATH because there is a system program named mcomp too.

Please do not move or copy mcomp to a different location because it need to read the database files in the same directory.

## Installation from source codes

You need install the Perl module Config::Simple for the pipeline. This module is not required if you use individual modules.

```
cd moabs-v1.3.0.src.x86_64_Linux.data
make
make install 
```

The make command will compile source files and generate system dependent dynamic binary executables. The make install command will overwrite the static binaries in moabs-v1.3.0.src.x86_64_Linux.data/bin/ with the dynamic ones.

You need also export PATH=/your/path/moabs-v1.3.0.src.x86_64_Linux.data/bin/:$PATH to use the command without involving full path.

If you encounter make errors, please contact us.

## Fast Manual

One may go to the test/ and type moabs --cf mytestrun.cfg to get started if you have installed the Perl Config::Simple module. In addition, you may also directly use the binary for each individual function in the moabs-v1.2.3/bin/. In short, one may simply finish the whole processing of bisulfite data for two conditions by typing

```
moabs --cf my_research_config_file
```

or

```
moabs -i wt_r1.fq -i wt_r2.fq -i ko_r1.fq -i ko_r2.fq
```

Done!

An example use of individual module mcomp setting replicate variance containing biological variance

```
cd test/
../bin/mcomp -r wt_r1.bam.G.bed,wt_r2.bam.G.bed -r ko_r1.bam.G.bed,ko_r2.bam.G.bed -m wildtype -m knockout -c comp.wiVar.txt --withVariance 1 -p 4 
```

## Cite Our Paper

Deqiang Sun, Yuanxin Xi, Benjamin Rodriguez, Hyun Jung Park, Pan Tong, Mira Meong, Margaret A Goodell and Wei Li. MOABS: model based analysis of bisulfite sequencing data. Genome Biology, 15 (2014)

Link:http://genomebiology.com/2014/15/2/R38/abstract

## Contact

The package is developed by Deqiang Sun. Please post any questions, suggestions or problems to the [MOABS Discussion google group](https://groups.google.com/d/forum/moabs_msuite) or send email to Deqiang Sun at moabs_msuite@googlegroups.com.

You are welcome to subscribe to the [moabs MOABS Discussion google group](https://groups.google.com/d/forum/moabs_msuite) for updates.

## Related projects

mSuite: [DNA Methylation analysis in one Suite](https://code.google.com/p/msuite/), which packages ExactNumCI, MOABS, MEPS and other tools.

MEPS: MEthylation Pipeline and Service on cloud

