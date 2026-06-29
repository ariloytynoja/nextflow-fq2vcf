#!/usr/bin/env nextflow

import groovy.yaml.YamlSlurper

params.samples         = 'samples.yaml'
params.reference       = 's3://processed_tuukka/superXY/reference/Saimaa01_superXY.fa.gz'
params.refsubsets      = 's3://processed_tuukka/superXY/reference/ref_subsets.tgz'
params.cramreference   = '/scratch/project_2009641/references/Saimaa01_superXY.fa.gz'
params.filterbam       = true

include { BWA_INDEX; ALIGN_BWA; REALIGN_INDELS; MERGE_BAMS; RUN_HAPLOTYPECALLER; MERGE_GVCFS; CONVERT_CRAM_AND_PUBLISH; PUBLISH_GVCF } from './modules/preprocess.nf'
include { ALIGN_MINIMAP2; RUN_DEEPVARIANT; CONVERT_CRAM_AND_PUBLISH as CONVERT_CRAM_AND_PUBLISH_HIFI; PUBLISH_GVCF_VCF } from './modules/preprocess.nf'
include { DELETE_BAMS; DELETE_FINALS; CONVERT_BAM } from './modules/preprocess.nf'

workflow {

    // PREPARATIONS
    // Download and bwa indexing of genome, download of genome subset files, and setup of local cram reference
    //  
    (refindex,subsets) = BWA_INDEX(params.reference,params.refsubsets,params.cramreference)
    subsets_count_ch = subsets.map { it.size() }

    // sample data from yaml file, allows mixing Illumina and Pacbio (but not yet)
    //
    inputs = new YamlSlurper().parse(params.samples as File)

    samples = Channel.fromList(inputs['samples']).
        map { s ->
            tuple([ id: s.sample_id, bam: s.bam, gvcf: s.gvcf, type: s.type ], s.fastq )
        }
        .branch { s ->
            sread_fq:  s[0].type == "SR" && s[1] != null
            sread_bam: s[0].type == "SR" && s[1] == null && s[0].bam != null
            hifi:  s[0].type == "HIFI"
        }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    // SHORT-READs 
    // mapping, realigning and merging of split analysis
    //
    bwaBams = ALIGN_BWA(samples.sread_fq,refindex)
    
    bamSubs = bwaBams.combine(subsets.flatten())
    realBams = REALIGN_INDELS(bamSubs,refindex)

    groupedBams = realBams
        .combine(subsets_count_ch) 
        .map { meta, bam, bai, subset, total -> 
            [ groupKey(meta.id.toString(), total.toInteger()), bam, bai, subset, meta ] 
        }
        .groupTuple()
        .map { gkey, bams, bais, subsets, metas ->
            def sorted_indices = new ArrayList(subsets.indices).sort { subsets[it] }
            def meta = metas[0] 
            return [ meta, sorted_indices.collect { bams[it] }, sorted_indices.collect { bais[it] } ]            
        }
    
    DELETE_BAMS(bwaBams,groupedBams)

    finalBams = MERGE_BAMS(groupedBams)
    

    existingCrams = samples.sread_bam.map { meta, fastq -> [ meta, meta.bam, "${meta.bam}.crai" ] }
    existingBams = CONVERT_BAM(existingCrams)

    allSRBams = finalBams.mix(existingBams)

    // variant calling with GATK/HC
    //
    finalSubs = allSRBams.combine(subsets.flatten())
    hcGvcfs = RUN_HAPLOTYPECALLER(finalSubs,refindex)
    
    groupedGvcfs = hcGvcfs
        .combine(subsets_count_ch)
        .map { meta, gvcf, tbi, subset, total -> 
            [ groupKey(meta.id.toString(), total.toInteger()), gvcf, tbi, subset, meta ] 
        }
        .groupTuple()
        .map { gkey, gvcfs, tbis, subsets, metas ->
            def sorted_indices = new ArrayList(subsets.indices).sort { subsets[it] }
            def meta = metas[0] 
            return [ meta, sorted_indices.collect { gvcfs[it] }, sorted_indices.collect { tbis[it] } ]            
        }
    finalGvcfs = MERGE_GVCFS(groupedGvcfs)
    
    // publishing of final files
    //
    publishedCrams = CONVERT_CRAM_AND_PUBLISH(finalBams)
    publishedGvcfs = PUBLISH_GVCF(finalGvcfs)
    
    DELETE_FINALS(finalBams,publishedCrams,finalGvcfs,publishedGvcfs)


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    // HIFI READs 
    //
    finalBams  = ALIGN_MINIMAP2(samples.hifi,refindex)
    finalGvcfs  = RUN_DEEPVARIANT(finalBams,refindex)


    // publishing of final files
    //
    publishedCrams = CONVERT_CRAM_AND_PUBLISH_HIFI(finalBams)
    publishedGvcfs = PUBLISH_GVCF_VCF(finalGvcfs)

 }
