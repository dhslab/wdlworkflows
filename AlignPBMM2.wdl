version 1.0

workflow AlignPBMM2 {

  input {
    File SampleSheet
    String Name
    String OutputDir

    String Reference = "/storage1/fs1/dspencer/Active/spencerlab/refdata/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.mmi"

    String JobGroup = "/dspencer/adhoc"
    String Queue = "dspencer"
    String Docker = "dhspence/docker-pacbio:latest"
  }

  # sample sheet should be like: ccs_read_path name
  Array[Array[String]] inputData = read_tsv(SampleSheet)
  
  scatter(samples in inputData){
    call pbmm2_align_sort {
      input: InBam=samples[0],
      Name=samples[1],
      Reference=Reference,
      jobGroup=JobGroup,
      queue=Queue,
      docker="dhspence/docker-pacbio:latest"
    }
  }

  call merge_bams {
    input: Bams=pbmm2_align_sort.bam,
    Name=Name,
    jobGroup=JobGroup,
    queue=Queue,
    docker="dhspence/docker-baseimage:031822"
  }

  call gather_files {
    input: OutputFiles=[merge_bams.bam,merge_bams.bai],
    	   OutputDir=OutputDir,
	   jobGroup=JobGroup,
      	   queue=Queue,
  }
  
}



task pbmm2_align_sort {

  input {
    String InBam
    String Name
    String Reference

    String queue
    String docker
    String jobGroup

  }
  
  command {
    /opt/conda/bin/pbmm2 align ~{Reference} ~{InBam} "~{Name}.aligned.bam" --preset CCS -j 8 --sort --log-level INFO
  }

  runtime {
      docker_image: docker
      cpu: "8"
      memory: "64 G"
      queue: queue
      job_group: jobGroup
  }

  output {
      File bam = "~{Name}.aligned.bam"
      File bai = "~{Name}.aligned.bam.bai"
  }
}


task merge_bams {
  input {
    Array[File] Bams
    String Name
    String queue
    String jobGroup
    String docker
  }
  
  command {
    samtools merge -@ 8 --write-index -o "~{Name}.5mc_aligned.bam" ${sep=' ' Bams}
  }
  
  runtime {
    docker_image: docker
    cpu: "8"
    memory: "64 G"
    queue: queue
    job_group: jobGroup
  }
  output {
    File bam = "~{Name}.5mc_aligned.bam"
    File bai = "~{Name}.5mc_aligned.bam.bai"
  }
}

task gather_files {
  input {
    Array[String] OutputFiles
    String OutputDir
    String queue
    String jobGroup
  }


  output {
    String out = stdout()
  }
  command <<<

    /bin/mv -f -t ~{OutputDir} ~{sep=" "  OutputFiles}
  
  >>>
  runtime {
    queue: queue
    docker_image: "ubuntu:xenial"
    cpu: "1"
    job_group: jobGroup
    memory: "10 G"
  }

}