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


Channel
    .fromPath( params.samples_csv )
    .splitCsv( header: true, sep: ',' )
    .map { row ->  tuple(row.sample_id, file(row.fastq_path)) }
    .set { sample_run_ch }

(cellranger_input, meta_input,scrublet_input) = sample_run_ch.into(3)


Channel
    .fromPath( params.samples_csv )
    .splitCsv( header: true, sep: ',' )
    .map { row ->  row.sample_id }
    .set { sample_id }

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

/* set ref directory e.g. (mouse,human,fly) */
ref = params.refdir

process cellranger_count {

  publishDir (
  path: "${params.outdir}/cellrangersouts",
  mode: 'copy',
  overwrite: 'true',
    )

  input:
  ref 
  tuple val(sample_id), file(fastq_path) from cellranger_input

  output:
  file("cellrangersouts") into cellrangers_outs,outs_dir
  script: 
  """
  ${params.cellranger_software_path}/cellranger count --id=$sample_id \
                   --transcriptome=$ref \
                   --fastqs=$fastq_path \
                   --sample=$sample_id \
                   --localcores=30 \
                   --localmem=64
  """

}

process scrublet {
	
  publishDir (
        path: "${params.outdir}",
        mode: 'copy',
        overwrite: 'true',
  )		
	input:

  each sample_id 
	file(cellranger) from cellrangers_outs.collect()

  output:

  file("scrubletdir") into scrubletdir
	script:
	"""
  mkdir -p scrubletdir
	python3 ${baseDir}/scripts/scrublet_multi.py ${cellranger} ${params.scrublet_SUFFIX} scrubletdir/ ${sample_id}
	"""
}

process add_meta {

  publishDir (
        path: "${params.outdir}",
        mode: 'copy',
        overwrite: 'true',
  )	
        
    input:
    each sample_id 
    file(outs) from outs_dir.collect()
               
    output: 
    
    	path("${sample_id}.h5") into meta_added
    
    script:

    """
    Rscript /mnt/data0/projects/biohub/software/rna_seq_pipeline_bwh/tenx_metadata_rna_adder.r \
  -i ${outs}/${sample_id}/outs/filtered_feature_bc_matrix.h5   \
  -l ${outs}/${sample_id}/outs/molecule_info.h5 \
  -s ${outss}/${sample_id}/outs/metrics_summary.csv \
  -k ${params.in_key} \
  -j ${sample_id} \
  """

}


process QC_Report {

    publishDir params.qc_output 			// output dir
	
    input:
    
    file(in_h5) from meta_added.collect()
    output: 
    	
    script:
    """
    Rscript /mnt/data0/projects/biohub/software/qc_batch_summary.r \
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
