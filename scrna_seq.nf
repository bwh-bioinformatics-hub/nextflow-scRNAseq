#!/usr/bin/ nextflow

nextflow.enable.dsl=1

/*
========================================================================================
   QC Report Generator Nextflow Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------
*/


// Replace this with the path to a directory containing raw fastq files
params.fastqs_dir = '/mnt/data0/projects/donglab/EWA_Ruifeng2023/data/fastqs/'
fastq_path = params.fastqs_dir


Channel
    .fromPath( params.samples_csv )
    .splitCsv( header: true, sep: ',' )
    .map { row ->  row.sample_id }
    .set { sample_id_ch }

(sample_id,samples,sample) = sample_id_ch.into(3)


println """\
         RNA Seq - N F   P I P E L I N E
         ===================================
         Experiment       		     : ${params.experiment_id}
         Samplesheet        		   : ${params.in_key}
         CellRangersOuts Directory : ${params.cellrangers_outdir}
         QC Report input directory : ${params.qc_in_dir}
         QC Report Output directory: ${params.qc_output}
         """
         .stripIndent()

process scrublet {

  publishDir (
        path: "${params.outdir}/scrubletdir",
        mode: 'copy',
        overwrite: 'true',
  )	 
    input:
    each samples 
    path(outs_dir) from cellrangers_outs.collect()            
    output: 
    
    path("${samples}_scrublet.{logic,score}") into scrublet_out
    
    script:

    """
	python3 ${baseDir}/scripts/scrublet_multi.py ${params.cellrangers_dir} ${params.scrublet_SUFFIX} ${samples}
    """

}

process add_meta {

  publishDir (
        path: "${params.outdir}",
        mode: 'copy',
        overwrite: 'true',
  )	
        
    input:
    each sample 
    path(ranger) from cellranger.collect()          
    output: 
    
    	path("${sample}.h5") into meta_added
    
    script:

    """
    Rscript path/to/rna_seq_pipeline_bwh/tenx_metadata_rna_adder.r \
  -i ${params.cellrangers_dir}/${sample}/outs/filtered_feature_bc_matrix.h5   \
  -l ${params.cellrangers_dir}/${sample}/outs/molecule_info.h5 \
  -s ${params.cellrangers_dir}/${sample}/outs/metrics_summary.csv \
  -k ${params.in_key} \
  -j ${sample} \
  """

}

process QC_Report {

    publishDir params.qc_output 			// output dir
	
    input:
    
    file(in_h5) from meta_added.collect()
    file(scrub_dir) from scrublet_out.collect()
    path(cellrangerouts) from counts.collect()
    output: 
    	
    script:
    """
    Rscript /home/acicalo/bioinf_hub/pipeline_dir/HUB_rachel_2023/qcreporter/qc_batch_summary.r \
    	-e  ${params.experiment_id} \
    	-m  'scrna' \
    	-i  ${params.qc_in_dir} \
    	-z  ${cellrangerouts} \
    	-u  ${scrub_dir} \
    	-f  ${params.refdir} \
    	-k  ${params.in_key}   \
    	-d  ${params.qc_output} \
    	-o  ${params.qc_output}/${params.experiment_id}_rnaseq_sample_report.html 
  """
}

