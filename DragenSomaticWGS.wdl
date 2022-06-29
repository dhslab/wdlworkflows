version 1.0

workflow DragenSomaticWGS {

  input {
    File TumorFastq1
    File TumorFastq2
    File Fastq1
    File Fastq2
    String Name
    String ID
    String SM
    String LB
    String PU
    String TumorID
    String TumorSM
    String TumorLB
    String TumorPU
    String OutputDir

    String DBSNP
    String DOCM
    String Reference

    String JobGroup
    String DragenQueue = "duncavagee"
    String DragenDocker = "docker.io/seqfu/centos7-dragen-3.10.4:latest"
  }

  call dragen_align {
      input: TumorFQ1=TumorFastq1,TumorFQ2=TumorFastq2,
      NormalFQ1=Fastq1,NormalFQ2=Fastq2,
      Name=Name,
      RGID=ID,
      RGSM=SM,
      RGLB=LB,
      RGPU=PU,
      TumorRGID=TumorID,
      TumorRGSM=TumorSM,
      TumorRGLB=TumorLB,
      TumorRGPU=TumorPU,
      Reference=Reference,
      DBSNP=DBSNP,
      DOCM=DOCM,
      jobGroup=JobGroup,
      queue=DragenQueue,
      docker=DragenDocker,
      OutputDir=OutputDir
    }

}



task dragen_align {

  input {
    File TumorFQ1
    File TumorFQ2  
    File NormalFQ1
    File NormalFQ2
  
    String Name
    String RGID
    String RGSM
    String RGLB
    String RGPU
    String TumorRGID
    String TumorRGSM
    String TumorRGLB
    String TumorRGPU
    
    String LocalAlignDir = "/staging/tmp/" + Name

    String OutputDir
  
    String Reference
    String? DBSNP
    String? DOCM
    String queue
    String docker
    String jobGroup

  }
  command {
      if [ ! -d "~{LocalAlignDir}" ]; then
        /bin/mkdir ~{LocalAlignDir}
      fi

      /opt/edico/bin/dragen -r ~{Reference} \
      -1 ~{NormalFQ1} \
      -2 ~{NormalFQ2} \
      --tumor-fastq1 ~{TumorFQ1} \
      --tumor-fastq2 ~{TumorFQ2} \
      --RGID ~{RGID} --RGSM ~{RGSM} --RGLB ~{RGLB} --RGPU ~{RGPU} \
      --RGID-tumor ~{TumorRGID} --RGSM-tumor ~{TumorRGSM} --RGLB-tumor ~{TumorRGLB} --RGPU-tumor ~{TumorRGPU} \
      --read-trimmers adapter \
      --trim-adapter-read1 /staging/runs/dev/adapter1.fa \
      --trim-adapter-read2 /staging/runs/dev/adapter2.fa \
      --enable-map-align true \
      --enable-map-align-output true \
      --enable-bam-indexing true \
      --enable-duplicate-marking true \
      --qc-coverage-ignore-overlaps=true \
      --gc-metrics-enable=true \
      --enable-variant-caller true \
      --vc-combine-phased-variants-distance 3 \
      ${"--dbsnp " + DBSNP} \
      ${"--vc-somatic-hotspots " + DOCM} \
      --enable-sv true \
      --sv-output-contigs true \
      --sv-hyper-sensitivity true \
      --sv-enable-liquid-tumor-mode true \
      --enable-cnv true \
      --cnv-use-somatic-vc-baf true \
      --output-format CRAM \
      --output-directory ~{LocalAlignDir}/ \
      --output-file-prefix ~{Name} && \      
      /bin/mv ~{LocalAlignDir}/* ~{OutputDir} && rm -Rf ~{LocalAlignDir}
  }

  runtime {
      docker_image: docker
      cpu: "20"
      memory: "200 G"
      queue: queue
      job_group: jobGroup
  }

  output {
      String outdir = "~{OutputDir}"
  }
}


task annotate_variants {
  input {
    String Vcf
    String refFasta
    String Vepcache
    String CustomAnnotationVcf
    String CustomAnnotationIndex
    String CustomAnnotationParameters
    String OutputDir
    String Name
    String queue
    String jobGroup
    String? tmp
    String docker
  }
  
  command {
    set -eo pipefail && \
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
    --format vcf --vcf --fasta ${refFasta} --hgvs --symbol --term SO --per_gene -o ${Name}.annotated.vcf \
    -i ~{Vcf} --custom ~{CustomAnnotationVcf},~{CustomAnnotationParameters} --offline --cache --max_af --dir ~{Vepcache} && \
    /opt/htslib/bin/bgzip -c ~{Name}.annotated.vcf > ~{OutputDir}/~{Name}.annotated.vcf.gz && \
    /usr/bin/tabix -p vcf ~{OutputDir}/~{Name}.annotated.vcf.gz && \
  }
  
  runtime {
    docker_image: docker
    cpu: "1"
    memory: "32 G"
    queue: queue
    job_group: jobGroup
  }
  output {
    File annotated_vcf = "~{OutputDir}/~{Name}.annotated.vcf.gz"
    File annotated_vcf_index = "~{OutputDir}/~{Name}.annotated.vcf.gz.tbi"
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