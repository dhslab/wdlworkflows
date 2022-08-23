version 1.0

workflow Chipseq {
  input {
    String Name
    String ID
    String SM
    String LB
    String PU
    String PL = "illumina"
    String OutputDir
    String JobGroup
    String Queue = "timley"
    
    File? Fastq1
    File? Fastq2
    File? AlignedReads
    String? AlignedReadsReference

    String Reference    
    String ReferenceIndex
    String Dictionary
    String Blacklist
    String Annotation

  }
  
  if (defined(AlignedReads)) {
    call revert_aligned_reads {
      input:
        Bam = AlignedReads,
        Reference = select_first([AlignedReadsReference, Reference]),
        queue = Queue,
        jobGroup = JobGroup
    }
  }
  
  call trim_fastqs {
    input:
      F1 = select_first([revert_aligned_reads.read1, Fastq1]),
      F2 = select_first([revert_aligned_reads.read2, Fastq2]),
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call align_and_tag {
    input:
      F1 = trim_fastqs.read1,
      F2 = trim_fastqs.read2,
      refFasta = Reference,
      label = Name,
      readGroup = '@RG' + '\\tID:' + ID + '\\tLB:' + LB + '\\tSM:' + SM + '\\tPU:' + PU + '\\tPL:' + PL,
      queue = Queue,
      jobGroup = JobGroup
  }

  call collect_insert_metrics {
    input:
      in = align_and_tag.bamfile,
      label = Name,
      ref = Reference,
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call fingerprint {
    input:
      in = align_and_tag.bamfile,
      label = Name,
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call flagstat {
    input:
      in = align_and_tag.bamfile,
      label = Name,
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call gene_coverage {
    input:
      BW = make_bw.bigwig_file,
      GTF = Annotation,
      Name = Name,
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call make_bw {
    input:
      in = align_and_tag.bamfile,
      label = Name,
      queue = Queue,
      jobGroup = JobGroup,
      Blacklist = Blacklist
  }
  
  call gather_result as gather_files {
    input:
      dir = OutputDir,
      files = [trim_fastqs.report, align_and_tag.bamfile,align_and_tag.bamindex, collect_insert_metrics.isOut, collect_insert_metrics.isPDF, fingerprint.fingerprint_table, fingerprint.fingerprint_plot, gene_coverage.profile, gene_coverage.plot, flagstat.fsOut, make_bw.bigwig_file],
      queue = Queue,
      jobGroup = JobGroup
  }
  
  call remove_files as rm_files {
    input:
      files = select_all([revert_aligned_reads.read1, revert_aligned_reads.read2, trim_fastqs.read1, trim_fastqs.read2]),
      order_by = align_and_tag.bamfile,
      queue = Queue,
      jobGroup = JobGroup
  } 
}

task make_bw {
  input {
    String in
    String label
    String queue
    String jobGroup
    String Blacklist
    Int? genome_size
  }


  output {
    File bigwig_file = "~{label}.bw"
  }
  command <<<

    export PYTHONPATH=/opt/conda/lib/python3.6/site-packages/ && \
    /opt/conda/bin/bamCoverage --bam ~{in} -o "~{label}.bw" --effectiveGenomeSize ~{default="2451960000"  genome_size} --normalizeUsing RPGC --ignoreDuplicates -bl ~{Blacklist} --binSize 10 --minMappingQuality 1 --extendReads -p 4 -ignore X Y MT
  
  >>>
  runtime {
    queue: queue
    docker_image: "dhspence/docker-deeptools"
    cpu: "4"
    job_group: jobGroup
    memory: "32 G"
  }

}

task fingerprint {
  input {
    String in
    String label
    String queue
    String jobGroup
  }


  output {
    File fingerprint_table = "~{label}_fingerprint.tab"
    File fingerprint_plot = "~{label}_fingerprint.png"
  }
  command <<<

    export PYTHONPATH=/opt/conda/lib/python3.6/site-packages/ && \	
    /opt/conda/bin/plotFingerprint -b ~{in} --minMappingQuality 10 --skipZeros -T "~{label}_fingerprint" --plotFile "~{label}_fingerprint.png" \
    --outRawCounts "~{label}_fingerprint.tab" -p 2
  
  >>>
  runtime {
    queue: queue
    docker_image: "dhspence/docker-deeptools"
    cpu: "2"
    job_group: jobGroup
    memory: "16 G"
  }

}

task align_and_tag {
  input {
    String F1
    String F2
    String refFasta
    String label
    String readGroup
    String queue
    String jobGroup
  }

  command <<<

    set -eo pipefail && \
    /usr/local/bin/bwa-mem2 mem -K 100000000 -t 8 -Y -R '~{readGroup}' ~{refFasta} ~{F1} ~{F2} | \
    /usr/local/bin/samblaster --addMateTags | samtools view -Sb - > out.bam && \
    /usr/local/bin/samtools sort -@ 4 --write-index -T bamsort -O BAM -o "~{label}.bam##idx##~{label}.bam.bai" out.bam && \
    rm -f out.bam
  
  >>>
  runtime {
    queue: queue
    docker_image: "henrycwong/chipseq:tagged-alignment"
    cpu: "8"
    job_group: jobGroup
    memory: "64 G"
  }
  output {
    File bamfile = "~{label}.bam"
    File bamindex = "~{label}.bam.bai"	
  }
}

task convert_to_cram {
  input {
    String refFasta
    String bam
    String label
    String queue
    String jobGroup
  }
  
  output {
    File cramfile = "~{label}.cram"
    File cramindex = "~{label}.cram.crai"
  }
  command <<<

    /usr/local/bin/samtools view -C -T ~{refFasta} ~{bam} > "~{label}.cram"; /usr/local/bin/samtools index "~{label}.cram"
  
  >>>
  runtime {
    queue: queue
    docker_image: "registry.gsc.wustl.edu/genome/samtools-1.3.1-2:2"
    cpu: "1"
    job_group: jobGroup
    memory: "8 G"
  }

}

task exit_early {
  input {
    String message
    String queue
    String jobGroup
    String docker
  }

  command <<<

    echo ~{message}
  
  >>>
  runtime {
    queue: queue
    docker_image: docker
    cpu: "1"
    job_group: jobGroup
    memory: "1 G"
  }

}

task revert_aligned_reads {
  input {
    String? Bam
    String Reference
    String queue
    String jobGroup
  }


  output {
    File read1 = "R1.fastq.gz"
    File read2 = "R2.fastq.gz"
    Array[String] readgroup = read_lines("read_group.txt")
  }
  command <<<

    set -eo pipefail && \
    /usr/local/bin/samtools view -H ~{Bam} | /usr/bin/awk '/^@RG/' | sed "s/\"//g" | sed 's/\t/\\t/g' > "read_group.txt" && \
    (set -eo pipefail && /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar RevertSam INPUT=~{Bam} OUTPUT=/dev/stdout SORT_ORDER=queryname VALIDATION_STRINGENCY=SILENT | \
    /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar SamToFastq I=/dev/stdin F="R1.fastq.gz" F2="R2.fastq.gz" NON_PF=true VALIDATION_STRINGENCY=SILENT)
  
  >>>
  runtime {
    queue: queue
    docker_image: "dhspence/docker-trimgalore"
    cpu: "1"
    job_group: jobGroup
    memory: "32 G"
  }

}

task flagstat {
  input {
    String in
    String label
    String queue
    String jobGroup
  }


  output {
    File fsOut = "~{label}.flagstat.txt"
  }
  command <<<

    /usr/local/bin/samtools flagstat ~{in} > "~{label}.flagstat.txt"
  
  >>>
  runtime {
    queue: queue
    docker_image: "registry.gsc.wustl.edu/genome/tagged-alignment:2"
    cpu: "1"
    job_group: jobGroup
    memory: "10 G"
  }

}

task gene_coverage {
  input {
    String BW
    String GTF
    String Name
    String queue
    String jobGroup
  }


  output {
    File profile = "~{Name}.genes.tab.gz"
    File plot = "~{Name}.genes.png"
  }
  command <<<

    /opt/conda/bin/computeMatrix scale-regions -S ~{BW} -R ~{GTF} -p 4 -b 2000 -m 2000 -a 2000 -bs 100 -o "~{Name}.genes.tab.gz" && \
    /opt/conda/bin/plotProfile -m "~{Name}.genes.tab.gz" -o "~{Name}.genes.png" && \
    echo -e "# title: 'Normalized coverage across genes'\n# description: 'This plot shows normalized coverage across gene bodies and 2kb upstream and downstream of gene starts and ends, respectively'\n# section: 'Custom Data File'\n# format: 'tsv'\n# plot_type: 'linegraph'\n# pconfig:\n#    id: 'gene_coverage_graph'\n#    categories: True\n#    ylab: 'Normalized coverage'\n#    xlab: 'Relative gene position'" > "~{Name}.gene_coverage_mqc.tsv" && \
    /usr/bin/Rscript --quiet -e 'require(rjson); args <- commandArgs(T); pars <- fromJSON(gsub("@","",readLines(args[1],n=1))); dat <- colMeans(read.table(args[1],sep="\t",skip=1)[,-c(1:6)],na.rm=T); labs <- c(seq(-pars$upstream+pars[[4]],-pars[[4]],by=pars[[4]]),paste0(seq(0,100,length=pars[[3]]/pars[[4]]+1),"%"),paste0("+",seq(pars[[4]],pars$downstream,by=pars[[4]]))); write.table(cbind(labs,dat),file=args[2],quote=F,sep="\t",append=T,col.names=F,row.names=F)' "~{Name}.genes.tab.gz" "~{Name}.gene_coverage_mqc.tsv"
  
  >>>
  runtime {
    queue: queue
    docker_image: "dhspence/docker-deeptools"
    cpu: "4"
    job_group: jobGroup
    memory: "8 G"
  }

}

task collect_insert_metrics {
  input {
    String in
    String label
    String ref
    String queue
    String jobGroup
  }


  output {
    File isOut = "~{label}.insert_size_summary.txt"
    File isPDF = "~{label}.insert_size.pdf"
  }
  command <<<

    /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar CollectInsertSizeMetrics INPUT=~{in} OUTPUT="~{label}.insert_size_summary.txt" \
    HISTOGRAM_FILE="~{label}.insert_size.pdf" ASSUME_SORTED=true
  
  >>>
  runtime {
    queue: queue
    docker_image: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
    cpu: "2"
    job_group: jobGroup
    memory: "10 G"
  }

}

task gather_result {
  input {
    String dir
    Array[String] files
    String queue
    String jobGroup
  }


  output {
    String out = stdout()
  }
  command <<<

    /bin/mv -t ~{dir}/ ~{sep=" "  files}
  
  >>>
  runtime {
    docker_image: "ubuntu:xenial"
    queue: queue
    job_group: jobGroup
    memory: "2 G"
  }

}

task trim_fastqs {
  input {
    String F1
    String F2
    String queue
    String jobGroup
  }


  output {
    File report = glob("*.html")[0]
    File read1 = glob("*R1*trim.fastq.gz")[0]
    File read2 = glob("*R2*trim.fastq.gz")[0]
  }
  command <<<

    fastp -w 8 -i ~{F1} -o $(basename ~{F1})".trim.fastq.gz" -I ~{F2} -O $(basename ~{F2})".trim.fastq.gz"
  
  >>>
  runtime {
    queue: queue
    docker_image: "dhspence/docker-fastp"
    cpu: "8"
    job_group: jobGroup
    memory: "32 G"
  }

}

task remove_file {
  input {
    String file
    String order_by
    String queue
    String jobGroup
  }


  output {
    String out = stdout()
  }
  command <<<

    /bin/rm ~{file}
  
  >>>
  runtime {
    docker_image: "ubuntu:xenial"
    memory: "2 G"
    queue: queue
    job_group: jobGroup
  }

}

task remove_files {
  input {
    Array[String] files
    String order_by
    String queue
    String jobGroup
  }


  output {
    String out = stdout()
  }
  command <<<

    /bin/rm ~{sep=" "  files}
  
  >>>
  runtime {
    docker_image: "ubuntu:xenial"
    memory: "2 G"
    queue: queue
    job_group: jobGroup
  }

}

