#!/bin/sh 

outbam=${1%.*}.bam

if [ ! -f $1 ]; then
	echo "$1 does not exist."
	exit 1
fi

echo "Sorting BAM ..."
mv $1 $1.unsorted.bam
samtools sort $1.unsorted.bam ${outbam%.*}
if [ $? -ne 0 ]; then
	echo "BAM file sorting not sucessful."
	echo "$outbam is in unsorted BAM format".
	mv $1.unsorted.bam $outbam
	exit 1
fi
rm $1.unsorted.bam

echo "Indexing BAM ..."
samtools index $outbam
if [ $? -ne 0 ]; then
	echo "BAM file indexing not sucessful."
	exit 1
fi
exit 0
