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
         CellRangersOuts Directory         : ${params.cellrangers_outdir}
         QC Report input directory 	   : ${params.qc_in_dir}
         QC Report Output directory 	   : ${params.qc_output}
         """
         .stripIndent()


process cellranger_count {

  when:
  params.fastq_path

  input:
  ID from params.ID
  ref from params.refdir
  fastq_path from params.fastq_path

  output:
  file("cellrangersouts") into cellrangers_outs
  script: 
  """
  mkdir -p ${baseDir}/results/cellrangersouts/
  cd ${baseDir}/results/cellrangersouts/

  ${params.cellranger_software_path}/cellranger count --id=$ID \
                   --transcriptome=$ref \
                   --fastqs=$fastq_path \
                   --sample=$ID \
                   --localcores=30 \
                   --localmem=64
  """

}

process scrublet {
	
	publishDir params.scrublet_OUT, mode: params.publish_dir_mode
	
	input:

	    each sample_id 
	output:

	script:
	"""
	
	python3 ${baseDir}/scrublet_multi.py ${params.cellrangers_outdir} ${params.scrublet_SUFFIX} ${params.scrublet_OUT} ${sample_id}
	
	"""
}

Channel
    .fromPath( params.samples_csv )
    .splitCsv( header: true, sep: ',' )
    .map { row ->  row.sample_id }
    .set { sample_id }

process add_meta {

	publishDir params.outdir, mode: params.publish_dir_mode

        
    input:
    	each sample_id 

               
    output: 
    
    	path("${sample_id}.h5") into meta_added
    
    script:

    """
    Rscript /home/acicalo/software/rna_seq_pipeline_bwh/tenx_metadata_rna_adder.r \
  -i ${params.cellrangers_outdir}/${sample_id}/outs/filtered_feature_bc_matrix.h5   \
  -l ${params.cellrangers_outdir}/${sample_id}/outs/molecule_info.h5 \
  -s ${params.cellrangers_outdir}/${sample_id}/outs/metrics_summary.csv \
  -k ${params.in_key} \
  -j ${sample_id} \
  """

}


process QC_Report {

    publishDir params.qc_output 			// output dir
	
    input:
    
    in_h5 from meta_added.collect()
    output: 
    	
    script:
    """
    Rscript /home/acicalo/software/qcreporter/qc_batch_summary.r \
    	-e  ${params.experiment_id} \
    	-m  'scrna' \
    	-i  ${params.qc_in_dir} \
    	-z  ${params.cellrangersout} \
    	-u  ${params.scrublet_dir} \
    	-f  ${params.refdir} \
    	-k  ${params.in_key}   \
    	-d  ${params.qc_output} \
    	-o  ${params.qc_output}/${params.experiment_id}_rnaseq_sample_report.html 
  """

}
