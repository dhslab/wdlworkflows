version 1.0

struct FastqPair {
	String name
    File read_1
    File read_2
}

workflow micro_c{
	input{
		Array[FastqPair] fastq
		File reference_index
		String outdir
		String jobGroup = "/spencerlab/micro-c"
		String docker = "atex91/micro-c"
		Boolean run_lib_complexity
		Boolean run_cleanup
	}

	call create_index{
		input: idx_tar = reference_index,
		jobGroup = jobGroup,
		docker = docker
	}

	scatter(fastq_pair in fastq){

		call makedir{
			input: name = fastq_pair.name,
			outdir = outdir,
			jobGroup = jobGroup
		}

		call align{
			input: idx_tar = reference_index,
			fastq_pair = fastq_pair,
			jobGroup = jobGroup,
			docker = docker
		}

		call parse{
			input: ref_genome = create_index.out_genome,
			aligned_bam = align.out_aligned,
			jobGroup = jobGroup,
			docker = docker
		}

		call split{
			input: deduped = parse.out_parsed,
			jobGroup = jobGroup,
			docker = docker,
			name = fastq_pair.name
		}

		call samtools_sort{
			input: unsorted_bam = split.out_bam,
			jobGroup = jobGroup,
			docker = docker,
			name = fastq_pair.name
		}

		call bam_to_bed{
			input: sorted_bam = samtools_sort.out_bam,
			jobGroup = jobGroup,
			docker = docker
		}

		call remove_files as rm_early{
			input: files = select_all([align.out_aligned, parse.out_parsed]),
			order_by = split.out_mapped,
			jobGroup = jobGroup
		}

		call remove_files as rm_late{
			input: files = select_all([split.out_bam]),
			order_by = samtools_sort.out_index,
			jobGroup = jobGroup
		}
	
		call calc_stats{
			input: stats = parse.stats,
			jobGroup = jobGroup,
			docker = docker,
			name = fastq_pair.name
		}

		if (run_lib_complexity){
			call lib_complexity{
				input: input_bed = bam_to_bed.out_bed,
				jobGroup = jobGroup,
				docker = docker,
				name = fastq_pair.name
			}
		}

		call create_hic{
			input: ref_genome = create_index.out_genome,
			mapped_pairs = split.out_mapped,
			jobGroup = jobGroup,
			name = fastq_pair.name
		}

		call create_cooler{
			input: ref_genome = create_index.out_genome,
			mapped_pairs = split.out_mapped,
			jobGroup = jobGroup,
			docker = docker,
			name = fastq_pair.name
		}

		if (run_lib_complexity){
			call fetch_output as with_complexity{
				input: dir = select_first(read_lines(makedir.out)),
				files = select_all([samtools_sort.out_bam, samtools_sort.out_index, split.out_mapped, calc_stats.qc_stats, lib_complexity.out_complexity, create_hic.hic, create_cooler.cooler]),
				jobGroup = jobGroup
			}
		}

		if (!run_lib_complexity){
			call fetch_output as without_complexity{
				input: dir = select_first(read_lines(makedir.out)),
				files = select_all([samtools_sort.out_bam, samtools_sort.out_index, split.out_mapped, calc_stats.qc_stats, create_hic.hic, create_cooler.cooler]),
				jobGroup = jobGroup,
			}
		}

		if (run_cleanup){
			call remove_files as cleanup{
				input: files = select_all([bam_to_bed.out_bed, parse.stats, create_cooler.fixed_res_cooler]),
				order_by = select_first([with_complexity.out, without_complexity.out]),
				jobGroup = jobGroup
			}
		}
	}

	output{
		Array[File] bam = samtools_sort.out_bam
		Array[File] index = samtools_sort.out_index
		Array[File] pairs = split.out_mapped
		File genome = create_index.out_genome
		Array[File] stats = calc_stats.qc_stats
		Array[File?] complexity = lib_complexity.out_complexity
		Array[File] hic = create_hic.hic
		Array[File] cooler = create_cooler.cooler
	}
	
}

task create_index{
	input{
		File idx_tar
		String jobGroup
		String docker
	}

	command <<<
		mkdir reference
		cd reference && tar -xvf ~{idx_tar}
		index_folder=$(ls)
		reference_fasta=$(ls | head -1)
		reference_folder=$(pwd)
		reference_index_path=$reference_folder/$reference_fasta
		cd ..
		samtools faidx $reference_index_path
		cut -f1,2 $reference_index_path.fai > hg38.genome
	>>>

	output{
		File out_genome = "hg38.genome"
	}

	runtime{
		cpu: "1"
		memory: "4 G"
		job_group: jobGroup
		docker_image: docker
	}
}

task align{
	input{
		File idx_tar
		FastqPair fastq_pair
		String jobGroup
		String docker
	}

	command <<<
		mkdir reference
		cd reference && tar -xvf ~{idx_tar}
		index_folder=$(ls)
		reference_fasta=$(ls | head -1)
		reference_folder=$(pwd)
		reference_index_path=$reference_folder/$reference_fasta
		cd ..
		bwa mem -5SP -T0 -t16 $reference_index_path ~{fastq_pair.read_1} ~{fastq_pair.read_2} | \
		samtools view -hbS > aligned.bam
	>>>

	output{
		File out_aligned = "aligned.bam"
	}

	runtime{
		cpu: "16"
		memory: "64 G"
		job_group: jobGroup
		docker_image: docker
	}
}

task parse{
	input{
		String ref_genome
		String aligned_bam
		String jobGroup
		String docker
	}

	command <<<
		mkdir temp/
		cd temp/
		tmpdir=$(pwd)
		cd ..
		pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 \
		--nproc-out 8 --chroms-path ~{ref_genome} ~{aligned_bam} | pairtools sort --nproc 8 --tmpdir=$tmpdir | \
		pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt \
		--output-dups dedup.pairsam.gz --output dedup.pairsam.gz
		rm -r $tmpdir
	>>>

	output{
		File out_parsed = "dedup.pairsam.gz"
		File stats = "stats.txt"
	}

	runtime{
		cpu: "8"
		memory: "32 G"
		job_group: jobGroup
		docker_image: docker
	}
}

task split{
	input{
		String deduped
		String jobGroup
		String docker
		String name
	}

	command <<<
		pairtools split --nproc-in 8 --nproc-out 8 --output-pairs ~{name}.mapped.pairs \
		--output-sam unsorted.bam ~{deduped}
	>>>

	output{
		File out_mapped = "~{name}.mapped.pairs"
		File out_bam = "unsorted.bam"
	}

	runtime{
		cpu: "8"
		memory: "16 G"
		job_group: jobGroup
		docker_image: docker
	}
}

task samtools_sort{
	input{
		String unsorted_bam
		String jobGroup
		String docker
		String name
	}

	command <<<
		mkdir temp/
		cd temp/
		tmpdir=$(pwd)
		cd ..
		touch $tmpdir/tempfile.bam
		/opt/samtools/bin/samtools sort -@16 -T $tmpdir/tempfile.bam -o ~{name}.mapped.PT.bam ~{unsorted_bam}
		rm -r $tmpdir
		/opt/samtools/bin/samtools index ~{name}.mapped.PT.bam
	>>>

	output{
		File out_bam = "~{name}.mapped.PT.bam"
		File out_index = "~{name}.mapped.PT.bam.bai"
	}

	runtime{
		cpu: "8"
		memory: "16 G"
		job_group: jobGroup
		docker_image: docker
	}
}

task bam_to_bed{
	input{
		String sorted_bam
		String jobGroup
		String docker
	}

	command <<<
		bamToBed -i ~{sorted_bam} > mapped.PT.bed
	>>>

	output{
		File out_bed = "mapped.PT.bed"
	}

	runtime{
		cpu: "8"
		memory: "16 G"
		job_group: jobGroup
		docker_image: docker
	}
}

task remove_files{
	input{
		Array[String] files
		String order_by
		String jobGroup
	}

	command <<<
		/bin/rm -f ~{sep=" " files}
	>>>

	output{
		String out = stdout()
	}

	runtime{
		cpu: "1"
		memory: "4 G"
		job_group: jobGroup
		docker_image: "ubuntu:xenial"
	}
}

task calc_stats{
	input{
		String stats
		String jobGroup
		String docker
		String name
	}

	command <<<
		python3 /usr/bin/get_qc.py -p ~{stats} > ~{name}.qc_stats.txt
	>>>

	output{
		File qc_stats = "~{name}.qc_stats.txt"
	}

	runtime{
		cpu: "1"
		memory: "16 G"
		job_group: jobGroup

		docker_image: docker
	}
}

task lib_complexity{
	input{
		String input_bed
		String jobGroup
		String docker
		String name
	}

	command <<<
		/opt/preseq/bin/preseq lc_extrap ~{input_bed} -pe -extrap 2100000000 -step 100000000 -seg_len 1000000000 -output ~{name}.out.preseq
	>>>

	output{
		File out_complexity = "~{name}.out.preseq"
	}

	runtime{
		cpu: "4"
		memory: "16 G"
		job_group: jobGroup
		docker_image: docker
	}
}

task create_hic{
	input{
		String ref_genome
		String mapped_pairs
		String jobGroup
		String name
	}

	command <<<
		java -Xmx64g  -Djava.awt.headless=true -jar /opt/scripts/common/juicer_tools.jar pre --threads 16 ~{mapped_pairs} ~{name}.contact_map.hic ~{ref_genome}
	>>>

	output{
		File hic = "~{name}.contact_map.hic"
	}

	runtime{
		cpu: "8"
		memory: "64 G"
		job_group: jobGroup
		docker_image: "henrycwong/pipeline:encodedcc"
	}
}

task create_cooler{
	input{
		String ref_genome
		String mapped_pairs
		String jobGroup
		String docker
		String name
	}

	command <<<
		cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ~{ref_genome}:1000 ~{mapped_pairs} cooler.1000.cool 
		cooler zoomify cooler.1000.cool -o ~{name}.contact_map.mcool
	>>>

	output{
		File cooler = "~{name}.contact_map.mcool"
		File fixed_res_cooler = "cooler.1000.cool"
	}

	runtime{
		cpu: "8"
		memory: "64 G"
		job_group: jobGroup
		docker_image: docker
	}
}

task makedir{
	input{
		String name
		String outdir
		String jobGroup
	}

	command <<<
		mkdir -p ~{outdir}/~{name}
		echo ~{outdir}/~{name} > dirname.txt
	>>>

	output{
		File out = "dirname.txt"
	}

	runtime{
		cpu: "1"
		memory: "2 G"
		job_group: jobGroup
		docker_image: "ubuntu:xenial"
	}
}

task fetch_output{
	input{
		String dir
		Array[String] files
		String jobGroup
	}

	command <<<
		/bin/mv -t ~{dir}/ ~{sep=" " files}
		echo "~{sep=" " files} moved" > outfile.txt
	>>>

	output{
		File out = "outfile.txt"
	}

	runtime{
		cpu: "1"
		memory: "4 G"
	    job_group: jobGroup
	    docker_image: "ubuntu:xenial"
	}
}