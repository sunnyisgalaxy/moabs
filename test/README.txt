1.
Please edit mytestrun.cfg to modify the reference path for mm9.fa. File mm9.chrom.sizes must be in the same directory.
If you do not have them at hand, download it at https://s3.amazonaws.com/deqiangsun/software/moabs/ref/mm9.fa.gz. || mm9.chrom.sizes.

2.
Please enter command
../bin/moabs -cf mytestrun.cfg
to see how it works.

3.
To see how mcall and mcomp work, please enter command
rm *.bam.G.bed comp.*
and go to 2.

4.
To start a new run starting from fastq files, please enter comand
rm *.bam *.bam.G.bed comp.*
and go to 2.

5. 
I have put results from different runs into several folders.
Please see command.log and mytestrun.cfg to see what modules are used.

6.
One example to directly use the individual command is like
../bin/mcomp -r wt_r1.bam.G.bed,wt_r2.bam.G.bed -r ko_r1.bam.G.bed,ko_r2.bam.G.bed -m wildtype -m knockout -c comp.wiVar.txt --withVariance 1 -p 4 --reference ~/moabs/ref/mm9.fa
where --withVariance 1 tells that there is biological variance 
while --withVariance 0 tells that there is no biolgical variance but only technical variance (technical >> biological).
