process  BWA_INDEX {
    
    storeDir "reference"

    input:
    path(reference)
    path(refsubsets)
    val(cramreference)

    output:
    path("reference.*")
    path("S*.list")

    script:
    """
    if [ ! -e "${cramreference}" ]; then 
      samtools faidx ${reference}
      cp -L ${reference.name} ${cramreference}
      cp -L ${reference.name}.fai ${cramreference}.fai
      cp -L ${reference.name}.gzi ${cramreference}.gzi
    fi

    zcat ${reference} > reference.fa
    samtools faidx reference.fa
    bwa-mem2 index reference.fa
    samtools dict reference.fa > reference.dict

    tar xzf ${refsubsets}
    """
}

process ALIGN_BWA {

    tag "${meta.id}"

    input:
    tuple val(meta),path(fastqs)
    path(refindex)

    output:
    tuple val(meta), path("aligned.bam"), path("aligned.bam.bai")

    script:
    def refname = "${refindex[0].simpleName}.fa" 
    def bwa2_threads  = task.cpus - 2 
    def readgroup= "@RG\\tID:${meta.id}\\tLB:Lib\\tSM:${meta.id}\\tPL:Plfm"
    
    """
    set -euo pipefail

    fastq1="\$(echo ${fastqs} | xargs -n1 | grep -F -e .1.fq.gz -e _1.fq.gz -e .1.fastq.gz -e _1.fastq.gz -e R1.fastq.gz -e _R1_001.fastq.gz | sort | tr '\n' ' ' )"
    fastq2="\$(echo ${fastqs} | xargs -n1 | grep -F -e .2.fq.gz -e _2.fq.gz -e .2.fastq.gz -e _2.fastq.gz -e R2.fastq.gz -e _R2_001.fastq.gz | sort | tr '\n' ' ' )"
    echo \${fastq1}
    echo \${fastq2}
    
    bwa-mem2 mem -R \"${readgroup}\" -t ${bwa2_threads} -K 100000000 -Y ${refname}  <(cat \${fastq1}) <(cat \${fastq2}) \
    | samtools view -h - | samtools fixmate -m -  mapped.bam

    echo ${fastqs} | xargs realpath | xargs rm 
 
    if ${params.filterbam}; then
        samtools sort -@ ${bwa2_threads} mapped.bam -o sorted.bam
        rm mapped.bam 
        samtools view -h sorted.bam | filter_bam.sh | samtools view -b -o aligned.bam -
	    rm sorted.bam
    else
        samtools sort -@ ${bwa2_threads} mapped.bam -o aligned.bam
        rm mapped.bam 
    fi

    samtools index -@ ${bwa2_threads} aligned.bam
    """
}

process REALIGN_INDELS {

    tag "${meta.id}"

    input:
    tuple val(meta), path(alignedbam), path(alignedbai), path(subset)
    path(refindex)
    
    output:
    tuple val(meta), path("*_markdup.bam"), path("*_markdup.bam.bai"), val("${subset.simpleName}")

    script:
    def refname = "${refindex[0].simpleName}.fa"
    def outname = "${meta.id}_${subset.simpleName}_markdup.bam"

    """    
    gatk3_4G \
    -T RealignerTargetCreator \
    -R ${refname} \
    -I ${alignedbam} \
    -L ${subset} \
    -nt ${task.cpus} \
    -o targets.intervals
    
    gatk3_4G \
    -T IndelRealigner \
    -R ${refname} \
    -I ${alignedbam} \
    -L ${subset} \
    -targetIntervals targets.intervals \
    -o realigned.bam

    samtools markdup realigned.bam "${outname}"
    samtools index "${outname}"   
    rm realigned.bam
    """
}

process MERGE_BAMS {

    tag "${meta.id}"
    
    input:
    tuple val(meta), path(bams), path(bais)
    
    output:
    tuple val(meta), path("${meta.id}_markdup.bam"), path("${meta.id}_markdup.bam.bai")

    script:
    """
    samtools cat -o "${meta.id}_markdup.bam" ${bams}
    samtools index "${meta.id}_markdup.bam"

    echo ${bams} | xargs realpath | xargs rm
    echo ${bams} | xargs realpath | sed 's/bam/bam.bai/g' | xargs rm
    """
}


process DELETE_BAMS {
    
    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(bam2), path(bai2)
    
    script:
    """
    echo ${bam} | xargs realpath | xargs rm
    echo ${bai} | xargs realpath | xargs rm
    """

} 

process RUN_HAPLOTYPECALLER {

    tag "${meta.id}"
    
    input:
    tuple val(meta), path(alignedbam), path(alignedbai), path(subset)
    path(refindex)

    output:
    tuple val(meta), path("*.gvcf.gz"), path("*.gvcf.gz.tbi"), val("${subset.simpleName}")

    script:
    def refname = "${refindex[0].simpleName}.fa"
    def outname = "${meta.id}_${subset.simpleName}.gvcf.gz"

    """
    gatk4 --java-options '-Xmx4g' HaplotypeCaller \
    -R ${refname} \
    -I "${alignedbam}" \
    -L ${subset} \
    -O "${outname}" \
    -ERC GVCF
    """
}

process MERGE_GVCFS {

    tag "${meta.id}"
    
    input:
    tuple val(meta), path(gvcfs), path(tbis)
    
    output:
    tuple val(meta), path("${meta.id}_merged.gvcf.gz"), path("${meta.id}_merged.gvcf.gz.tbi")

    script:
    """
    bcftools concat -Oz -o "${meta.id}_merged.gvcf.gz" ${gvcfs}
    tabix "${meta.id}_merged.gvcf.gz" 

    echo ${gvcfs} | xargs realpath | xargs rm
    echo ${gvcfs} | xargs realpath | sed 's/gz/gz.tbi/g' | xargs rm
    """
}

process CONVERT_CRAM_AND_PUBLISH {

    tag "${meta.id}"
   
    publishDir "${file(meta.bam).parent.toUriString()}", mode: 'copy'

    input:
    tuple val(meta), path(alignedbam), path(alignedbai)
    
    output:
    tuple val(meta), path("${file(meta.bam).name}"), path("${file(meta.bam).name}.crai")

    script:
    """    
    samtools view -@ ${task.cpus} -C ${alignedbam} -o ${file(meta.bam).name} -T ${params.cramreference}
    samtools index ${file(meta.bam).name}
    """
}

process PUBLISH_GVCF {

    tag "${meta.id}"
   
    publishDir "${file(meta.gvcf).parent.toUriString()}", mode: 'copy'

    input:
    tuple val(meta), path(gvcf), path(tbi)
    
    output:
    tuple val(meta), path("${file(meta.gvcf).name}"), path("${file(meta.gvcf).name}.tbi")

    script:
    def outname = file(meta.gvcf).name
    def s3_dir  = file(meta.gvcf).parent.toUriString()
    """
    cat > s3cfg << EOF
    [default]
    access_key = ${System.getenv('NF_AWS_ACCESS_KEY')}
    secret_key = ${System.getenv('NF_AWS_SECRET_KEY')}
    host_base  = a3s.fi
    host_bucket = %(bucket)s.a3s.fi
    use_https  = True
    EOF

    s3cmd --config s3cfg put -F ${gvcf}   ${s3_dir}/${outname}
    s3cmd --config s3cfg put -F ${tbi}    ${s3_dir}/${outname}.tbi

    s3cmd --config s3cfg info ${s3_dir}/${outname}     > /dev/null || { echo "ERROR: gvcf upload failed"; exit 1; }
    s3cmd --config s3cfg info ${s3_dir}/${outname}.tbi > /dev/null || { echo "ERROR: tbi upload failed"; exit 1; }

    mv ${gvcf} ${outname}
    mv ${tbi} ${outname}.tbi
    """
}

process DELETE_FINALS {
    
    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(cram), path(crai)
    tuple val(meta2), path(gvcf), path(tbi)
    tuple val(meta3), path(gvcf2), path(tbi2)
    
    script:
    """
    echo ${bam} | xargs realpath | xargs rm
    echo ${bai} | xargs realpath | xargs rm
    echo ${cram} | xargs realpath | xargs rm
    echo ${crai} | xargs realpath | xargs rm
    echo ${gvcf} | xargs realpath | xargs rm
    echo ${tbi} | xargs realpath | xargs rm
    echo ${gvcf2} | xargs realpath | xargs rm
    echo ${tbi2} | xargs realpath | xargs rm
    """

} 

/////////////////////////////////////////////////////7

process ALIGN_MINIMAP2 {

    tag "${meta.id}"

    input:
    tuple val(meta), path(fastqs)
    path(refindex)

    output:
    tuple val(meta),
          path("${meta.id}_mm2aligned.bam"),
          path("${meta.id}_mm2aligned.bam.bai")

    script:
    def refname = "${refindex[0].simpleName}.fa"
    def mm2_threads  = task.cpus - 2
    def sort_threads = 2

    def files = fastqs.collect { "\"$it\"" }.join(' ')
    def bam_count = fastqs.count { it.name.endsWith('.bam') }
    def fq_count  = fastqs.size() - bam_count

    """
    set -euo pipefail

    BAM_COUNT=${bam_count}
    FQ_COUNT=${fq_count}

    if [[ \$BAM_COUNT -gt 0 && \$FQ_COUNT -gt 0 ]]; then
        echo "ERROR: Mixed BAM and FASTQ inputs are not supported"
        exit 1
    fi

    if [[ \$BAM_COUNT -gt 0 ]]; then
        echo "Input detected as BAM"

        samtools merge -u - ${files} | \
        samtools fastq - | \
        minimap2 -t ${mm2_threads} -ax map-hifi -Y --eqx ${refname} - | \
        samtools sort -@ ${sort_threads} -o ${meta.id}_mm2aligned.bam

    else
        echo "Input detected as FASTQ"

        minimap2 -t ${mm2_threads} -ax map-hifi -Y --eqx ${refname} \
          ${files} | \
        samtools sort -@ ${sort_threads} -o ${meta.id}_mm2aligned.bam
    fi

    samtools index ${meta.id}_mm2aligned.bam
    """
}

process RUN_DEEPVARIANT {
 
    tag "${meta.id}"

    input:
    tuple val(meta), path(alignedbam), path(alignedbai)
    path(refindex)

    output:
    tuple val(meta), path("${meta.id}_deepv.gvcf.gz"), path("${meta.id}_deepv.gvcf.gz.tbi"), path("${meta.id}_deepv.vcf.gz"), path("${meta.id}_deepv.vcf.gz.tbi")

    script:
    def refname = "${refindex[0].simpleName}.fa"

    """
    mkdir -p \$(echo \$TMPDIR | sed 's#/tmp#temp#')    
    singularity run \
        -B temp:/tmp -B /scratch:/scratch -B /usr/lib/locale/:/usr/lib/locale/ \
        /projappl/project_2006976/singularity/deepvariant_1.10.0-beta.sif  /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${refname} --reads ${alignedbam} --output_gvcf "${meta.id}_deepv.gvcf.gz" --output_vcf "${meta.id}_deepv.vcf.gz" \
        --num_shards ${task.cpus} --intermediate_results_dir /tmp

    tabix -f "${meta.id}_deepv.gvcf.gz"
    tabix -f "${meta.id}_deepv.vcf.gz"
    """
}

process PUBLISH_GVCF_VCF {

    tag "${meta.id}"
   
    publishDir "${file(meta.gvcf).parent.toUriString()}", mode: 'copy'

    input:
    tuple val(meta), path(gvcf), path(gtbi), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("${file(meta.gvcf).name}"), path("${file(meta.gvcf).name}.tbi"), path("${file(meta.gvcf).simpleName}.vcf.gz"), path("${file(meta.gvcf).simpleName}.vcf.gz.tbi")

    script:
    """    
    mv ${gvcf} ${file(meta.gvcf).name}
    mv ${gtbi} "${file(meta.gvcf).name}.tbi"
    mv ${vcf} "${file(meta.gvcf).simpleName}.vcf.gz"
    mv ${tbi} "${file(meta.gvcf).simpleName}.vcf.gz.tbi"
    """
}

process CONVERT_BAM {

    tag "${meta.id}"
   
    input:
    tuple val(meta), path(alignedcram), path(alignedcrai)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai")

    script:
    """    
    samtools view -@ ${task.cpus} ${alignedcram} -o "${file(meta.bam).simpleName}.bam" -b
    samtools index "${file(meta.bam).simpleName}.bam"
    """
}
