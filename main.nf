#!/usr/bin/env nextflow
/*
========================================================================================
                         hybrid-assembly
========================================================================================
 hybrid-assembly Analysis Pipeline. Started 2018-03-20.
 #### Homepage / Documentation
 https://github.com/kevinmenden/hybrid-assembly
 #### Authors
 Kevin Menden kevinmenden <kevin.menden@t-online.de> - https://github.com/kevinmenden>
 Tristan Kast tristankast <tristan.kast@gmx.de> - https://github.com/TrisKast
----------------------------------------------------------------------------------------
*/
//

def helpMessage() {
    log.info"""
    =========================================
     hybrid-assembly //v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run tristankast/hybrid_assembly --shortReads '*_R{1,2}.fastq.gz' --longReads nano_reads.fastq.gz --assembler masurca --masurca_genomsize 100200000  -profile galaxy

    Mandatory arguments:
      --assembler                   The assembler pipeline to choose. One of 'spades' | 'canu' | 'masurca'
      --shortReads                  The paired short reads
      --longReads                   The long reads
      -profile                      Hardware config to use.

    References                      If you want to use a reference genome
      --fasta                       Path to Fasta reference

    Options:
      --lr_type                     Long read technology. One of 'nanopore' | 'pacbio' . Default: 'nanopore'

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.fasta = false
params.shortReads = ""
params.longReads = ""
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.assembler = "masurca"
params.genomeSize = 0

multiqc_config = file(params.multiqc_config)

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
    if( params.genome ) log.info "Genome and reference specified. Reference will be used!"
} else {
    fasta = file("placeholder") // create placeholder
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


/*
 * Create a channel for input short read files
 */
Channel
    .fromFilePairs( params.shortReads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.shortReads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!" }
    .into { short_reads_qc; short_reads_assembly; short_reads_correction;sv_detection_sniffles_short}

///*
// * Create a channel for input long read files
// */
Channel
        .fromPath( params.longReads )
        .ifEmpty { exit 1, "Cannot find any long reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .into { long_reads_qc; long_reads_assembly; long_reads_scaffolding; sv_detection_sniffles_long }
        
///*
// * Create a channel for reference fasta file
// */
Channel
        .fromPath( params.fasta )
        .ifEmpty { exit 1, "Cannot find any reference file matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .into { spades_assembly; quast_reference_before_polishing; quast_reference_after_polishing; sv_reference_sniffles; sv_reference_assemblytics }



// Header log info
log.info "========================================="
log.info " hybrid-assembly"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Short Reads']  = params.shortReads
summary['Long Reads']   = params.longReads
summary['Fasta Ref']    = params.fasta
summary['Assembler']    = params.assembler
summary['Genome Size (Masurca)']    = params.masurca_genomesize
summary['Genome Size (Canu)']    = params.genomeSize
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container']    = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    spades.py --version > v_spades.txt
    canu --version > v_canu.txt
    quast --version > v_quast.txt
    minimap2 --version > v_minimap.txt
    pilon --version > v_pilon.txt
    nucmer -v > v_nucmer.txt
    ngmlr --version v_ngmlr.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/**
 * STEP 1.1 QC for short reads
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from short_reads_qc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/**
 * STEP 1.2 QC for long reads
 */
process nanoqc {
    tag "${lreads.baseName}"
    publishDir "${params.outdir}/nanoqc", mode: 'copy'

    input:
    file lreads from long_reads_qc

    output:
    file "*" into nanoqc_results

    script:
    ftype = (lreads.extension == "fasta" || lreads.extension == "fa") ? "--fasta" : "--fastq"
    """
    source activate nanoqc-env
    NanoPlot $ftype $lreads
    source deactivate nanoqc-env
    """
}

/**
 * STEP 2 Pre-processing
 */
 
/* process prinseq {
 
    publishDir "${params.outdir}/prinseq", mode: 'copy'
    
    input:
    file lreads from long_reads_filtering
    
    output:
    file "prinseq_good.fasta" into filtered_longreads
    file "*" into prinseq_results
    
    script:
    """
    prinseq-lite.pl -fasta $lreads -min_len 500 -out_good prinseq_good

    """
 
 }
 filtered_longreads.into{ long_reads_assembly; 
}
 */
 

/**
 * STEP 3 Assembly
 */

/**
 * SPAdes assembly workflow
 * 1. Hybrid assembly with SPAdes
 * 2. quast for assesment
 */
if (params.assembler == 'spades') {

    // Create assembly with SPAdes
    process spades {
        tag "$name"
        publishDir "${params.outdir}/spades", mode: 'copy'

        input:
        file fasta from spades_assembly
        set val(name), file(sreads) from short_reads_assembly
        file lreads from long_reads_assembly
        
        output:
        file "scaffolds.fasta" into assembly_results
        file "contigs.fasta" into assembly_results_contigs
        file "*" into spades_results

        script:
        ref_genome = params.fasta ? "--trusted-contigs $fasta" : ''
        lr = (params.lr_type == 'nanopore') ? '--nanopore' : '--pacbio'
        kmers = params.kmers
        """
        spades.py -o "spades_results" -t ${task.cpus} \\
        -m $params.mem_spades \\
        -1 ${sreads[0]} -2 ${sreads[1]} \\
        $lr $lreads $ref_genome \\
        -k $kmers
        mv spades_results/scaffolds.fasta scaffolds.fasta
        mv spades_results/contigs.fasta contigs.fasta
        """

    }
    assembly_results.into{ assembly_mapping; assembly_pilon; quast_wo_pilon }

}

/**
 * Canu assembly workflow
 * 1. assembly with Canu
 * 2. map short reads with minimap2
 * 3. polish assembly with pilon
 * 4. quast for assesment
 */
if (params.assembler == 'canu') {
    if (params.genomeSize == 0){
        log.error "No genome size specified. Necessary for Canu assembly workflow"
        exit 1
    }
    // Create assembly with Canu
    process canu {
        tag "${lreads.baseName}"
        publishDir "${params.outdir}/canu", mode: 'copy'

        input:
        file lreads from long_reads_assembly
        
        output:
        file "*contigs.fasta" into assembly_results
        file "*" into canu_results

        script:
        """
        canu \\
        -p ${lreads.baseName} genomeSize=$params.genomeSize -nanopore-raw $lreads gnuplotTested=true \\
        correctedErrorRate=$params.correctedErrorRate \\
        rawErrorRate=$params.rawErrorRate \\
        minReadLength=$params.minReadLength \\
        minOverlapLength=$params.minOverlapLength
        """
    }
    assembly_results.into{ assembly_mapping; assembly_pilon; quast_wo_pilon }

}

/**
 * MaSuRCA assembly workflow
 */
if (params.assembler == 'masurca') {
    // Generate MaSuRCA config file and run assembler
    process masurca {
        tag "$name"
        publishDir "${params.outdir}/masurca", mode: 'copy'

        input:
        set val(name), file(sreads) from short_reads_assembly
        file lreads from long_reads_assembly

        output:
        file "masurca_config.txt" into masurca_config_file
        file "final.genome.scf.fasta" into assembly_results
        
        
        script:
        cg = params.close_gaps ? "--close_gaps" : ""
        hc = params.high_cov ? "--high_cov" : ""
        """
        masurca_config.py \\
        --sr1 ${sreads[0]} --sr2 ${sreads[1]} \\
        --isize $params.insert_size --stdev $params.insert_stdv \\
        --lr $lreads --lr_type $params.lr_type \\
        --genome_size $params.masurca_genomesize \\
        $cg $hc -p ${task.cpus}

        masurca masurca_config.txt

        ./assemble.sh

        mv CA.mr*/final.genome.scf.fasta final.genome.scf.fasta
        """
    }
    assembly_results.into{ assembly_mapping; assembly_pilon; quast_without_polishing }
    
}

/**
 * STEP 4 Assembly Evaluation before polishing
 */

 // Assess assembly without polishing with quast
    process quast_before_polishing{
        publishDir "${params.outdir}/quast_plain_assembly", mode: 'copy'

        input:
        file fasta from quast_reference_before_polishing
        file scaffolds from quast_without_polishing
        

        output:
        file "*" into quast_results_without_polishing

        script:
        """
        quast $scaffolds -R $fasta --large --threads 20
        """
}


/**
 * STEP 5 Polishing
 */

  // Map short reads to assembly with minimap2
  process minimap2_assembly_polishing {
      tag "${sreads[0].baseName}"
      publishDir "${params.outdir}/minimap2", mode: 'copy'

      input:
      file assembly from assembly_mapping
      set val(name), file(sreads) from short_reads_correction

      output:
      file "*" into minimap_alignment_results
      file "*.sorted.bam" into short_reads_mapped_bam

      script:
      """
      minimap2 -t 15 -ax sr $assembly ${sreads[0]} ${sreads[1]} > sreads_assembly_aln.sam
      samtools view -h -b sreads_assembly_aln.sam > sreads_assembly_aln.bam
      samtools sort sreads_assembly_aln.bam > sreads_assembly_aln.sorted.bam
      """
  }


  // Polish assembly with pilon
  process pilon {
       publishDir "${params.outdir}/pilon", mode: 'copy'

       input:
       file sr_bam from short_reads_mapped_bam
       file assembly from assembly_pilon

       output:
       file "*" into pilon_results
       file "pilon.fasta" into pilon_scaffold

       script:
       """
       samtools index $sr_bam

       java -Xmx80G -jar /opt/conda/pkgs/pilon-1.22-1/share/pilon-1.22-1/pilon-1.22.jar --threads 15 --genome $assembly --bam $sr_bam
       """

  }
  pilon_scaffold.into{ quast_assembly_after_polishing; sv_detection_sniffles_assembly; sv_detection_assemblytics_assembly }


/**
 * STEP 6 Assembly Evaluation after polishing
 */


    // Assess assembly with quast
    process quast_after_polishing{
        publishDir "${params.outdir}/quast_after_polishing", mode: 'copy'

        input:
        file fasta from quast_reference_after_polishing
        file scaffolds from quast_assembly_after_polishing

        output:
        file "*" into quast_results

        script:
        """
        quast $scaffolds -R $fasta --large --threads 20
        """
}
 
 
 /**
 * STEP 7 SV-Detection
 */
 
 /**
 * STEP 7.1 Mapping-based approach
 */
 
 process minimap2_SV_calling{
      publishDir "${params.outdir}/ngmlr", mode: 'copy'
      
      input:
      file fasta from sv_reference_sniffles
      file lr from sv_detection_sniffles_long
      set val(name), file(sreads) from sv_detection_sniffles_short
      
      output:
      file "aln_long_sorted.bam" into sv_bam
      
      
      script:
      """
      minimap2 -ax map-ont $fasta $lr --MD > aln_long.sam
      samtools view -Sb aln_long.sam > aln_long.bam
      samtools sort aln_long.bam > aln_long_sorted.bam
      """
 }
 
 process sniffles{
        publishDir "${params.outdir}/sniffles", mode: 'copy'
        
        input: 
        file sorted_long from sv_bam

        
        output: 
        file "sniffles_only_long_reads.vcf" into sniffles_only_long_reads_vcf
        file "*" into sniffles_results
        
        
        script:
        """
        sniffles -m $sorted_long -v sniffles_only_long_reads.vcf
        """
 
 }
 
  /**
 * STEP 7.2 Assembly-based approach
 */
 
 process nucmer {
     publishDir "${params.outdir}/nucmer", mode: 'copy'
      
     input:
     file ref from sv_reference_assemblytics
     file assembly from sv_detection_assemblytics_assembly
      
     output:
     file "out.delta" into delta_file
     file "*" into nucmer_results
 
     script:
     """
     nucmer $ref $assembly
     """
 }
 
 process assemblytics {
      publishDir "${params.outdir}/assemblytics", mode: 'copy'
 
    input:
    file delta from delta_file
    
    output:
    file "*" into assemblytics_results
 
    script:
    """
    source activate python2_7-env
    
    Assemblytics $delta assemblytics 50 /Assemblytics/
    
    source deactivate python2_7-env
    """
 
 }

 
 
/*
 * Step 8 MultiQC
 * collect the results
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file ('quast_results/*') from quast_results

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {
    log.info "[hybrid-assembly] Pipeline Complete"

}
