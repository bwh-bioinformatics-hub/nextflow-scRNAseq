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
sample_id = params.sample_id

println """\
         RNA Seq - N F   P I P E L I N E
         ===================================
         Experiment       		   : ${params.experiment_id}
         Samplesheet        		   : ${params.in_key}
         QC Report input directory 	   : ${params.qc_in_dir}
         QC Report Output directory 	   : ${params.qc_output}
         """
         .stripIndent()

process scrublet {
	
	publishDir params.scrublet_OUTS, mode: params.publish_dir_mode
	
	input:
	
	path params.scrublet_PREFIX into PREFIX
	path params.scrublet_SUFFIX into SUFFIX
	each sample_id 
	
	output:
	
	script:
	"""
	
	python3 scrublet_multi.py ${PREFIX} ${SUFFIX} ${params.scrublet_OUTS} ${sample_id}
	
	"""
}











process add_meta {

	publishDir params.outdir, mode: params.publish_dir_mode

        
    input:
   
    	each sample_id 

               
    output: 
    
    	path("${sample_id}.h5") into meta_added
    
    script:

    """
    Rscript /home/acicalo/software/rna_seq_pipeline_bwh/tenx_metadata_rna_adder.r \
  -i ${params.cellrangersout}/${sample_id}/outs/filtered_feature_bc_matrix.h5   \
  -l ${params.cellrangersout}/${sample_id}/outs/molecule_info.h5 \
  -s ${params.cellrangersout}/${sample_id}/outs/metrics_summary.csv \
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
    Rscript /home/acicalo/software/qcreporter/qc_batch_summary_rachel_2022.r \
    	-e  ${params.experiment_id} \
    	-m  'scrna' \
    	-i  ${params.qc_in_dir} \
    	-z  ${params.cellrangersout} \
    	-u  ${params.scrublet_OUTS} \
    	-f  ${params.refdir} \
    	-k  ${params.in_key}   \
    	-d  ${params.qc_output} \
    	-o  ${params.qc_output}/${params.experiment_id}_rnaseq_sample_report.html 
  """

}


