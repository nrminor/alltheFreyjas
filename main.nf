#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info 	"""
			Do all the things, Freyja!
			(version 0.1.0)
			==========================

			I/O
			------
			bam directory           : ${params.bam_dir}
			reference fasta         : ${params.ref_fasta}
			results directory       : ${params.results}

			Which Freyja subcommands to run
			-------------------------------
			freyja demix            : ${params.demix}
			freyja get-lineage-def  : ${params.get_lineage_def}
			freyja boot             : ${params.boot}
			freyja aggregate        : ${params.aggregate}
			freyja plot             : ${params.plot}
			freyja dash             : ${params.dash}
			freyja relgrowthrate    : ${params.relgrowthrate}
			freyja extract          : ${params.extract}
			freyja filter           : ${params.filter}
			freyja covariants       : ${params.covariants}
			freyja plot-covariants  : ${params.plot_covariants}

			"""
			.stripIndent()

// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	// input channels
	ch_bams = Channel
		.fromPath( "${params.bam_dir}/*.bam" )
		.map { bam -> tuple( file(bam).getSimpleName(), file(bam) ) }

	ch_ref = Channel
		.fromPath( params.ref_fasta )

	// Workflow steps
	FREYJA_VARIANTS (
		ch_bams,
		ch_ref
	)

	FREYJA_DEMIX (
		FREYJA_VARIANTS.out
	)

	// FREYJA_GET_LINEAGE_DEF ( )

	// FREYJA_BOOT ( )

	FREYJA_AGGREGATE (
		FREYJA_DEMIX.out
		.map { sample_id, tsv -> file(tsv) }
		.collect( )
	)

	FREYJA_PLOT (
		FREYJA_AGGREGATE.out
	)

	// FREYJA_DASH ( )

	// FREYJA_RELGROWTHRATE ( )

	// FREYJA_EXTRACT ( )

	// FREYJA_FILTER ( )

	// FREYJA_COVARIANTS ( )

	// FREYJA_PLOT_COVARIANTS ( )


}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Additional parameters that are derived from parameters set in nextflow.config
params.variants_result = params.results + "/variants"
params.demix_result = params.results + "/demix"
params.get_lineage_def_result = params.results + "/get_lineage_def"
params.boot_result = params.results + "/boot"
params.aggregate_result = params.results + "/aggregate"
params.plot_result = params.results + "/plot"
params.dash_result = params.results + "/dash"
params.relgrowthrate_result = params.results + "/relgrowthrate"
params.extract_result = params.results + "/demextractix"
params.filter_result = params.results + "/filter"
params.covariants_result = params.results + "/covariants"
params.plot_covariants_result = params.results + "/plot_covariants"

// --------------------------------------------------------------- //




// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process FREYJA_VARIANTS {

	/*
	*/

	tag "${sample_id}"
	publishDir params.variants_result, mode: 'copy'

	input:
	tuple val(sample_id), path(bam)
	each path(ref)

	output:
	tuple val(sample_id), path("${sample_id}_variants.tsv"), path("${sample_id}.depth")

	script:
	"""
	freyja variants \
	${bam} \
	--variants ${sample_id}_variants.tsv \
	--depths ${sample_id}.depth \
	--ref ${ref}
	"""

}

process FREYJA_DEMIX {

	/*
	*/

	tag "${sample_id}"
	publishDir params.demix_result, mode: 'copy'

	input:
	tuple val(sample_id), path(variants), path(depths)

	output:
	tuple val(sample_id), path("${sample_id}.demix.tsv")

	when:
	params.demix == true

	script:
	"""
	freyja demix \
	${variants} \
	${depths} \
	--output ${sample_id}.demix.tsv \
	--confirmedonly
	"""

}

process FREYJA_GET_LINEAGE_DEF {

	/*
	*/

	tag "${sample_id}"
	publishDir params.get_lineage_def_result, mode: 'copy'

	input:

	output:

	when:
	params.get_lineage_def == true

	script:
	"""
	"""

}

process FREYJA_BOOT {

	/*
	*/

	tag "${sample_id}"
	publishDir params.boot_result, mode: 'copy'

	input:

	output:

	when:
	params.boot == true

	script:
	"""
	"""

}

process FREYJA_AGGREGATE {

	/*
	*/

	publishDir params.aggregate_result, mode: 'copy'

	input:
	path demix_files

	output:
	path "aggregated.tsv"

	when:
	params.aggregate == true

	script:
	"""
	freyja aggregate ./ --output aggregated.tsv --ext ".demix.tsv"
	"""

}

process FREYJA_PLOT {

	/*
	*/

	publishDir params.plot_result, mode: 'copy'

	input:
	path aggregated_tsv

	output:
	path "aggregated_plot.pdf"

	when:
	params.plot == true

	script:
	"""
	freyja plot ${aggregated_tsv} --output aggregated_plot.pdf
	"""

}

process FREYJA_DASH {

	/*
	*/

	tag "${sample_id}"
	publishDir params.dash_result, mode: 'copy'

	input:

	output:

	when:
	params.dash == true

	script:
	"""
	"""

}

process FREYJA_RELGROWTHRATE {

	/*
	*/

	tag "${sample_id}"
	publishDir params.relgrowthrate_result, mode: 'copy'

	input:

	output:

	when:
	params.relgrowthrate == true

	script:
	"""
	"""

}

process FREYJA_EXTRACT {

	/*
	*/

	tag "${sample_id}"
	publishDir params.extract_result, mode: 'copy'

	input:

	output:

	when:
	params.extract == true

	script:
	"""
	"""

}

process FREYJA_FILTER {

	/*
	*/

	tag "${sample_id}"
	publishDir params.filter_result, mode: 'copy'

	input:

	output:

	when:
	params.filter == true

	script:
	"""
	"""

}

process FREYJA_COVARIANTS {

	/*
	*/

	tag "${sample_id}"
	publishDir params.covariants_result, mode: 'copy'

	input:

	output:

	when:
	params.covariants == true

	script:
	"""
	"""

}

process FREYJA_PLOT_COVARIANTS {

	/*
	*/

	tag "${sample_id}"
	publishDir params.plot_covariants_result, mode: 'copy'

	input:

	output:

	when:
	params.plot_covariants == true

	script:
	"""
	"""

}

// --------------------------------------------------------------- //