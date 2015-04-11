#!/bin/sh 

samtools view -h $1 | perl -F'\t' -lane 'if($_ =~ /^\@/){print;next;} $_ =~ s/XS\:Z\:/ZS\:Z\:/; $_ =~ s/read1_F/\+\+/; $_ =~ s/read1_R/\-\+/; $_ =~ s/read2_F/\-\-/; $_ =~ s/read2_R/\+\-/; print;' | samtools view -bS - > $1.ok.bam
