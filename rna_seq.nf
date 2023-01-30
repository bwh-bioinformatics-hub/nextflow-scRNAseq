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

println """\
         RNA Seq - N F   P I P E L I N E
         ===================================
         Experiment       		   : ${params.experiment_id}
         Samplesheet        		   : ${params.in_key}
         Metadata input directory          : ${params.data_dir}
         QC Report input directory 	   : ${params.qc_in_dir}
         QC Report Output directory 	   : ${params.qc_output}
         """
         .stripIndent()







sample_id = params.sample_id


process add_meta {

	publishDir params.outdir, mode: params.publish_dir_mode

        
    input:
   
    	each sample_id 

               
    output: 
    
    	path("${sample_id}.h5") into meta_added
    
    script:

    """
    Rscript /home/acicalo/harvard/scRNA_scATAC_pipeline_dev/development_scripts/developmental_software/rna_seq_pipeline_bwh/tenx_metadata_rna_adder.r \
  -i ${params.data_dir}/${sample_id}/outs/filtered_feature_bc_matrix.h5   \
  -l ${params.data_dir}/${sample_id}/outs/molecule_info.h5 \
  -s ${params.data_dir}/${sample_id}/outs/metrics_summary.csv \
  -k ${params.in_key} \
  -j ${sample_id} \
  """

}

process QC_Report {

    publishDir params.qc_output 			// output dir
	
    input:
    
    	path in_h5 from meta_added.collect()


               
    output: 
    	
    script:
    """
    Rscript /home/acicalo/harvard/software/rna_atac_seq_qc_reporter/qc_batch_summary.r \
    	-e  ${params.experiment_id} \
    	-m  'scrna' \
    	-i  ${params.outdir} \
    	-k  ${params.in_key}   \
    	-d  ${params.qc_output} \
    	-o  ${params.qc_output}/${params.experiment_id}_rnaseq_sample_report.html 
  """

}


