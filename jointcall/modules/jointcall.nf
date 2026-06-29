process FETCH_REFERENCE {
    
    storeDir "reference"

    input:
    path(reference)
    
    output:
    path("reference.*")

    script:
    """
    zcat ${reference} > reference.fa
    samtools faidx reference.fa
    samtools dict reference.fa > reference.dict
    """
}

process FETCH_SUBSETS {
    
    input:
    path(refsubsets)
    
    output:
    path("S*.list")

    script:
    """
    tar xzf ${refsubsets}
    """
}

process STAGE_FILES {
    stageInMode "copy"

    input:
    tuple val(sample), path(gvcfgz), path(gvcftbi)

    output:
    tuple val(sample), path(gvcfgz), path(gvcftbi)

    script:
    """
    # No-op
    """
}

// An attempt to write the output only to TMPDIR. Works but ay not scale. Given up without proper testing for other ideas.
//
process COMBINE_GVCFS_TEMP {

    input:
    tuple val(sample), path(gvcfgz), path(gvcftbi), path(subset)
    path(refindex)

    output:
    tuple path("${subset.simpleName}.gvcf.gz"), path("${subset.simpleName}.gvcf.gz.tbi"), path(subset)

    script:
    def refname = "${refindex[0].simpleName}.fa"
    def outgvcf = "${subset.simpleName}.gvcf.gz"
    def java_mem = task.memory.toMega() - 1024
    """
    temp_dir=\$(mktemp -d -p "$LOCAL_SCRATCH")
    echo "\$temp_dir"
    df -h "\$temp_dir"

    files=\$(echo ${gvcfgz} | xargs | sed 's/ / -V /g')
    
    gatk4 --java-options "-Xmx${java_mem}m" CombineGVCFs \
    -R ${refname} \
    -V \${files} \
    -L ${subset} \
    -O "\$temp_dir/${outgvcf}"

    ln -s  "\$temp_dir/${outgvcf}" "${outgvcf}" 
    ln -s  "\$temp_dir/${outgvcf}.tbi" "${outgvcf}.tbi" 
    """
} 

process COMBINE_GVCFS {

    publishDir "temp_gvcf", mode: 'copy'

    input:
    tuple val(sample), path(gvcfgz), path(gvcftbi), path(subset)
    path(refindex)

    output:
    tuple path("${subset.simpleName}.gvcf.gz"), path("${subset.simpleName}.gvcf.gz.tbi"), path(subset)

    script:
    def refname = "${refindex[0].simpleName}.fa"
    def outgvcf = "${subset.simpleName}.gvcf.gz"
    def java_mem = task.memory.toMega() - 1024
    """
    files=\$(echo ${gvcfgz} | xargs | sed 's/ / -V /g')
    
    gatk4 --java-options "-Xmx${java_mem}m" CombineGVCFs \
    -R ${refname} \
    -V \${files} \
    -L ${subset} \
    -O ${outgvcf}
    """
} 

process CALL_GVCFS {

    publishDir "temp_vcf", mode: 'copy'

    input:
    tuple path(gvcfgz), path(gvcftbi), path(subset)
    path(refindex)

    output:
    tuple path("${subset.simpleName}.vcf.gz"), path("${subset.simpleName}.vcf.gz.tbi"), val("${subset.simpleName}")

    script:
    def refname = "${refindex[0].simpleName}.fa"
    def outvcf = "${subset.simpleName}.vcf.gz"
    def java_mem = task.memory.toMega() - 1024
    """
    gatk4 --java-options "-Xmx${java_mem}m" GenotypeGVCFs \
    -R ${refname} \
    -V ${gvcfgz} \
    -L ${subset} \
    -O ${outvcf} 
    """
} 


process COMBINE_CALL_GVCFS {

    input:
    tuple val(sample), path(gvcfgz), path(gvcftbi), path(subset)
    path(refindex)

    output:
    tuple path("${subset.simpleName}.vcf.gz"), path("${subset.simpleName}.vcf.gz.tbi"), val("${subset.simpleName}")

    script:
    def refname = "${refindex[0].simpleName}.fa"
    def outgvcf = "${subset.simpleName}.gvcf.gz"
    def outvcf = "${subset.simpleName}.vcf.gz"
    def java_mem = task.memory.toMega() - 1024
    """
    files=\$(echo ${gvcfgz} | xargs | sed 's/ / -V /g')
    
    gatk4 --java-options "-Xmx${java_mem}m" CombineGVCFs \
    -R ${refname} \
    -V \${files} \
    -L ${subset} \
    -O ${outgvcf} 
    
    gatk4 --java-options "-Xmx${java_mem}m" GenotypeGVCFs \
    -R ${refname} \
    -V ${outgvcf} \
    -L ${subset} \
    -O ${outvcf} 
    """
} 

process IMPORT_CALL_GVCFS {

    input:
    tuple val(sample), path(gvcfgz), path(gvcftbi), path(subset)
    path(refindex)

    output:
    tuple path("${subset.simpleName}.vcf.gz"), path("${subset.simpleName}.vcf.gz.tbi"), val("${subset.simpleName}")

    script:
    def refname = "${refindex[0].simpleName}.fa"
    def outvcf = "${subset.simpleName}.vcf.gz"
 
    def entries = [sample, gvcfgz].transpose().collect { s, g -> "${s}\t${g.simpleName}.g.vcf.gz" }.join("\\n")
    def rename1 = gvcfgz.collect { g -> "mv ${g.name} ${g.simpleName}.g.vcf.gz"}.join("\\n")
    def rename2 = gvcftbi.collect { g -> "mv ${g.name} ${g.simpleName}.g.vcf.gz.tbi"}.join("\\n")
    
    def java_mem = task.memory.toMega() - 1024
    """
    printf "${entries}\\n" > sample_map.txt
    printf "${rename2}\\n" > rename.sh
    source rename.sh
    printf "${rename1}\\n" > rename.sh
    source rename.sh

    gatk4 --java-options "-Xmx${java_mem}m" GenomicsDBImport \
       --genomicsdb-shared-posixfs-optimizations true \
       --batch-size 50 \
       --bypass-feature-reader \
       --genomicsdb-workspace-path my_database \
       --sample-name-map sample_map.txt \
       --tmp-dir $TMPDIR \
       --reader-threads ${task.cpus} \
       -L "${subset}" 

    gatk4 --java-options "-Xmx${java_mem}m" GenotypeGVCFs \
       -R ${refname} \
       -V gendb://my_database \
       -L ${subset} \
       -O ${outvcf} 
    """

}

process MERGE_VCFS {
    input:
    tuple path(vcfs), path(tbis)

    output:
    tuple path("merged.vcf.gz"), path("merged.vcf.gz.tbi")

    script:
    """
    bcftools concat -Oz -o merged.vcf.gz ${vcfs}
    tabix merged.vcf.gz 
    """
}

process PUBLISH_VCF {

    publishDir "results", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    val outname
    
    output:
    tuple path("${outname}"), path("${outname}.tbi")

    script:
    """
    mv ${vcf} ${outname}
    mv ${tbi} ${outname}.tbi
    """
}
