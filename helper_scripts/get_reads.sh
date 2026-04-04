#!/bin/bash
# Collect reads for a particular sample
# Usage: get_reads.sh <sample ID>
# Example: get_reads.sh HG00099

SAMPLE_ID=$1

PROJ_DIR=/private/groups/patenlab/fokamoto/centrolign
BED_DIR=/private/groups/patenlab/mira/centrolign/batch_submissions/centrolign/release2_QC_v2/per_smp_asat_beds

SAMPLE_TMP=$TMPDIR/faith_$SAMPLE_ID
mkdir -p $SAMPLE_TMP

reads=`grep "^$SAMPLE_ID," "$PROJ_DIR/to_align/aws_file_locations.csv" | cut -f3 -d ","`
echo "Downloading reads for $SAMPLE_ID from AWS"
aws s3 --no-sign-request cp "$reads" $SAMPLE_TMP/full.bam &> /dev/null
echo "Download complete"

for hap_num in 1 2; do
    hap_id="${SAMPLE_ID}.${hap_num}"

    # Subset all-chr centromere coordinate file to have only our chosen ones
    fgrep -e chr4 -e chr6 -e chr9 -e chr10 -e chr11 -e chr12 -e chr17 \
        ${BED_DIR}/${hap_id}_asat_arrays.bed > $SAMPLE_TMP/chosen_chroms.bed
    if [ -s $SAMPLE_TMP/chosen_chroms.bed ]; then
        echo "BED coordinates found for $hap_id"
        cat $SAMPLE_TMP/chosen_chroms.bed
    else
        echo "No BED coordinates found for $hap_id"
        continue
    fi

    # Subset to all chrs first, so we only have to use the giant file once
    samtools view -@32 -L $SAMPLE_TMP/chosen_chroms.bed -b $SAMPLE_TMP/full.bam > $SAMPLE_TMP/centromeres.bam
    samtools index $SAMPLE_TMP/centromeres.bam

    cat $SAMPLE_TMP/chosen_chroms.bed | while read line; do
        echo "$line"
        old_path_name=`echo $line | cut -f1 -d " "`
        start=`echo $line | cut -f2 -d " "`
        end=`echo $line | cut -f3 -d " "`
        chrom=`echo $line | cut -f4 -d " "`
        new_path_name="${SAMPLE_ID}#${hap_num}#${hap_id}"

        echo "samtools view -@32 -h $SAMPLE_TMP/centromeres.bam "$old_path_name" > $SAMPLE_TMP/chrom.sam"
        samtools view -@32 -h $SAMPLE_TMP/centromeres.bam "$old_path_name" > $SAMPLE_TMP/chrom.sam
        echo "Subset $hap_id reads for $chrom"

        samtools view $SAMPLE_TMP/chrom.sam | sed "s/$old_path_name/$new_path_name/" > $SAMPLE_TMP/no_header.sam
        samtools view -H $SAMPLE_TMP/chrom.sam | sed "s/$old_path_name/$new_path_name/" > $SAMPLE_TMP/header.txt
        # Filter for reads which appear within the BED file's boundaries
        # while also editing the coordinates to be graph-friendly
        ./helper_scripts/edit_sam.py --start "$start" --end "$end" \
            $SAMPLE_TMP/no_header.sam > $SAMPLE_TMP/edited.sam
        cat $SAMPLE_TMP/header.txt $SAMPLE_TMP/edited.sam > $SAMPLE_TMP/combined.sam
        
        # Single-hap graph is convenient for injection & simulation
        single_hap_graph=$PROJ_DIR/graph/haploid/${chrom}.${hap_id}
        ./helper_scripts/create_single_path_ref.sh $PROJ_DIR/graph/unsampled/${chrom}.pg \
            "${new_path_name}#0" "$single_hap_graph"
        read_prefix=$PROJ_DIR/to_align/${chrom}.${hap_id}
        # Get truth positions
        vg inject -x ${single_hap_graph}.gbz $SAMPLE_TMP/combined.sam > ${read_prefix}.real.gam
        vg filter --tsv-out "name;nodes" ${read_prefix}.real.gam > ${read_prefix}.real.tsv
        # Convert to FASTQ
        vg view --fastq-out ${read_prefix}.real.gam > ${read_prefix}.real.fastq

        echo "Produced real reads"

        # Generate similar simulated reads
        num_reads=`grep -c "^@" ${read_prefix}.real.fastq`
        vg sim -x ${single_hap_graph}.gbz --num-reads "$num_reads" \
            --random-seed 42 --threads 20 --fastq ${read_prefix}.real.fastq \
            --align-out --use-average-length > ${read_prefix}.sim.gam
        vg filter --tsv-out "name;nodes" ${read_prefix}.sim.gam > ${read_prefix}.sim.tsv
        # Convert to FASTQ
        vg view --fastq-out ${read_prefix}.sim.gam > ${read_prefix}.sim.fastq

        echo "Produced simulated reads"

        # No longer needed, strictly speaking
        rm ${read_prefix}.real.gam ${read_prefix}.sim.gam
    done
done

rm -rf $SAMPLE_TMP