#!/usr/bin/env nextflow

import groovy.yaml.YamlSlurper

params.samples         = 'samples.yaml'
params.reference       = 's3://processed_ribbon/reference/Ribbon_PR_18_6_2024.fa.gz'
params.refsubsets      = 's3://processed_ribbon/reference/ref_subsets.tgz'
params.vcfout          = "seals.vcf.gz"

// Many of these processes aren't scaling up to hundreds of samples
// 
include { FETCH_REFERENCE; FETCH_SUBSETS; STAGE_FILES; COMBINE_GVCFS; CALL_GVCFS; COMBINE_CALL_GVCFS; IMPORT_CALL_GVCFS; MERGE_VCFS; PUBLISH_VCF } from './modules/jointcall.nf'

workflow {

    // PREPARATIONS
    // Download genome, download of genome subset files
    //  
    refindex = FETCH_REFERENCE(params.reference)

    subsets = FETCH_SUBSETS(params.refsubsets)
    subsets_count_ch = subsets.map { it.size() }

    // sample data from yaml file, allows mixing Illumina and Pacbio (but not yet)
    //
    inputs = new YamlSlurper().parse(params.samples as File)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    gvcfs_ch = Channel.fromList(inputs['samples'])
        .map { s -> [ s.sample_id, file(s.gvcfgz), file(s.gvcftbi) ] }
        .collect()  
        .map { pairs ->
            def names = pairs.collate(3).collect { it[0] }  // extract sample IDs
            def gvcfs = pairs.collate(3).collect { it[1] }  // extract gvcf files
            def tbis  = pairs.collate(3).collect { it[2] }  // extract tbi files
            tuple(names, gvcfs, tbis)
        }

    // n option to stage *LOCAL* files in TMPDIR. Should not be done for s3 files as they are staged automatically.
    /* 
    gvcfs_ch = STAGE_FILES(gvcfs_ch) //*/

    gvcfSubs = gvcfs_ch.combine(subsets.flatten())

    // GenomicsDBImport and GenotypeGVCFs in one process
    /* vcfSubs = IMPORT_CALL_GVCFS(gvcfSubs,refindex)  // crashes, too slow? */
    
    // CombineGVCFs and GenotypeGVCFs in one process
    /* vcfSubs = COMBINE_CALL_GVCFS(gvcfSubs,refindex) // too slow */
    
    // CombineGVCFs and then storing combined data locally
    combGvcfs = COMBINE_GVCFS(gvcfSubs,refindex)
    combGvcfs | view

    // GenotypeGVCFs and then storing alled data locally
    vcfSubs = CALL_GVCFS(combGvcfs,refindex)
    vcfSubs | view

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    groupedVcfs = vcfSubs
        .combine(subsets_count_ch)
        .map { vcf, tbi, subset, total ->
            [ groupKey("all_vcfs", total.toInteger()), vcf, tbi, subset ]
        }
        .groupTuple()
        .map { gkey, vcfs, tbis, subsets ->
            def sorted_indices = new ArrayList(subsets.indices).sort { subsets[it] }
            [
                sorted_indices.collect { vcfs[it] },
                sorted_indices.collect { tbis[it] }
            ]
        }


    finalVcfs = MERGE_VCFS(groupedVcfs)

    // publishing of final files
    //
    publishedVcfs = PUBLISH_VCF(finalVcfs,params.vcfout)
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
        

 }
