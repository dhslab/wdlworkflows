#!/usr/local/bin/python

import sys, os, getopt, json, pandas as pd

def arg_parse():
	r1 = None #read 1 filepath. For specifying individual samples
	r2 = None #read 2 filepath. For specifying individual samples
	sample_names = None
	ref = "/storage1/fs1/dspencer/Active/spencerlab/refdata/hg38/all_sequences.tar.gz" #reference tarball
	samplemap = None #samplemap.csv path.
	inputs = None #directory containing fastq pairs. 
	wdl = "/storage1/fs1/dspencer/Active/spencerlab/alex/micro_c/pipeline/micro_c.wdl" #path to the .wdl
	outdir = None #directory to output files to
	clean = False #Optional. Use --clean to delete cromwell executions folder for each sample ran
	bulk = False #Optional. Use --bulk to run all samples in a directory when no Samplemap.csv is present
	complexity = False #Optional. Use --complexity to run the library complexity estimation task.
	docker = "atex91/micro-c" #sets the virtual environment to run the majority of the workflow in.
	jobgroup = "/spencerlab/micro-c-tasks" #sets the jobGroup parameter for LSF backends.
	config = "/storage1/fs1/dspencer/Active/spencerlab/alex/micro_c/pipeline/application.conf" #directory to LSF backend config file
	norun = False #Boolean. If True, no command will be run, but the input JSON will be printed
	groups = 1 #The number of groups to split your samples into. Defaults to running all samples at once (groups=1)

	argv = sys.argv[1:]

	try:
		opts, args = getopt.getopt(argv, "w:o:s:i:r:n:d:j:c:g:", ["bulk", "r1=", "r2=", "clean", "complexity", "norun"])
	except getopt.GetoptError as err:
		print(err)
		opts = []

	for opt, arg in opts:
		if opt in ["--bulk"]:
			bulk = True
		if opt in ["--r1"]:
			r1 = arg
		if opt in ["--r2"]:
			r2 = arg
		if opt in ["--clean"]:
			clean = True
		if opt in ["--complexity"]:
			complexity = True
		if opt in ["--norun"]:
			norun = True
		if opt in ["-w"]:
			wdl = arg
		if opt in ["-o"]:
			outdir = arg
		if opt in ["-s"]:
			samplemap = arg
		if opt in ["-i"]:
			inputs = arg
		if opt in ["-r"]:
			ref = arg
		if opt in ["-n"]:
			sample_names = [arg]
		if opt in ["-d"]:
			docker = arg
		if opt in ["-j"]:
			jobgroup = arg
		if opt in ["-c"]:
			config = arg
		if opt in ["-g"]:
			groups = int(arg)

	return r1, r2, samplemap, wdl, outdir, clean, inputs, ref, sample_names, bulk, complexity, docker, jobgroup, config, norun, groups

def parse_samplemap(samplemap, inputs):
	#Parses Samplemap.csv information into the appropriate format for use with the WDL
	samplemap = pd.read_csv(samplemap)
	samplemap = samplemap.sort_values("File Name")
	#samplemap['Sample Name'] = samplemap['Sample Name'] + '_' + samplemap['Flow Cell ID']
	ones = samplemap[samplemap['File Name'].str.contains('R1')]['File Name']
	twos = samplemap[samplemap['File Name'].str.contains('R2')]['File Name']
	sample_names = list(samplemap['Sample Name'].unique())
	fastq_pairs = [{"name": name, "read_1": f"{inputs}/{read_1}", "read_2": f"{inputs}/{read_2}"} for name, read_1, read_2 in zip(sample_names, ones, twos)]
	return fastq_pairs

def parse_directory(inputs):
	#Parses the data directory and generates a list of fastq pairs for use with the WDL. Also generates sample names for gathering output
	files = [file for file in os.listdir(inputs) if '.fastq' in file]
	ones = [file for file in files if '_R1_' in file]
	twos = [file for file in files if '_R2_' in file]
	sample_names = [fname.split('_R1_')[0] for fname in sorted(ones)]
	fastq_pairs = [{"name": name, "read_1": f"{inputs}/{read_1}", "read_2": f"{inputs}/{read_2}"} for name, read_1, read_2 in zip(sample_names, sorted(ones), sorted(twos))]
	return fastq_pairs

def get_output(outdir, sample_names):
	#Retrieves the output specified by the wdl workflow.
	with open(f"{outdir}/outfile.json", "r") as infile:
		my_json = json.load(infile)
	for key in my_json['outputs']:
		if type(my_json['outputs'][key]) is list:
			for path, sample in zip(my_json['outputs'][key], sample_names):
				if os.path.exists(f"{outdir}/{sample}"):
					os.system(f"mv {path} {outdir}/{sample}")
				else:
					os.mkdir(f"{outdir}/{sample}")
					os.system(f"mv {path} {outdir}/{sample}")
		else:
			for sample in sample_names:
				if os.path.exists(f"{outdir}/{sample}"):
					os.system(f"mv {my_json['outputs'][key]} {outdir}/{sample}")
				else:
					os.mkdir(f"{outdir}/{sample}")
					os.system(f"mv {my_json['outputs'][key]} {outdir}/{sample}")
			os.system(f"rm {my_json['outputs'][key]}")

def cleanup(outdir):
	#Removes the cromwell-executions directory for the sample.
	with open(f"{outdir}/outfile.json", "r") as infile:
		my_json = json.load(infile)
	executions = my_json['workflowRoot']
	os.system(f"rm -r {executions}")

def list_split(lst, num_groups):
	#Splits list lst into num_groups groups. Use num_groups=1 to run all samples at once.
	k, m = divmod(len(lst), num_groups)
	return (lst[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(num_groups))

def main():
	r1, r2, samplemap, wdl, outdir, clean, inputs, ref, sample_names, bulk, complexity, docker, jobgroup, config, norun, groups = arg_parse()
	
	if samplemap:
		fastq_pairs = parse_samplemap(samplemap, inputs)
	elif bulk:
		fastq_pairs = parse_directory(inputs)
	else:
		fastq_pairs = [{"name": sample_names[0], "read_1": r1, "read_2": r2}]

	for group_num, fastq_group in enumerate(list_split(fastq_pairs, groups)):
		input_json = {}
		input_json["micro_c.fastq"] = fastq_group
		input_json["micro_c.reference_index"] = ref
		input_json["micro_c.run_lib_complexity"] = complexity
		input_json["micro_c.run_cleanup"] = clean
		input_json["micro_c.outdir"] = f"{outdir}/group{group_num}"
		input_json["micro_c.docker"] = docker
		input_json["micro_c.jobGroup"] = jobgroup
		
		if norun:
			print(json.dumps(input_json, indent=4, sort_keys=True))
		else:
			os.system(f"mkdir -p {outdir}/group{group_num}")

			with open(f"{outdir}/group{group_num}/inputs.json", "w") as outfile:
				json.dump(input_json, outfile)

			os.system(f"LSF_DOCKER_VOLUMES=\"/storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active /scratch1/fs1/dspencer:/scratch1/fs1/dspencer $HOME:$HOME\" bsub -oo {outdir}/group{group_num}/out.log -eo {outdir}/group{group_num}/err.log -q dspencer -G compute-oncology -g /spencerlab/micro-c -R \"span[hosts=1]\" -a \"docker(henrycwong/cromwell)\" /usr/bin/java -Dconfig.file={config} -jar /cromwell/cromwell.jar run -t wdl -m {outdir}/group{group_num}/outfile.json -i {outdir}/group{group_num}/inputs.json {wdl}")

if __name__ == "__main__":
	main()