from __future__ import division
import os, sys, argparse, json, subprocess, re, glob
import pandas as pd
import pysam

def getch():
    """Read single character from standard input without echo."""
    import sys, tty, termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch

def yes_or_no(question):
    c = ""
    print(question + " (Y/N/C): ", end = "", flush = True,file=sys.stderr)
    while c not in ("y", "n","c"):
        c = getch().lower()
    if (c == "c"):
        print(file=sys.stderr)
        exit(1)

    if (c=='y'):
        print(file=sys.stderr)
    
    return c == 'y'
    
# static variables
conf = '/storage1/fs1/dspencer/Active/spencerlab/workflow/dspencer.default.cromwell.conf'
womtool='/storage1/fs1/dspencer/Active/spencerlab/workflow/womtool-71.jar'


#
# Note that all WDLs launched with this script need the following input arguments:
#
# Name, OutputDir, ID, SM, LB, PU, Reference, and either Fastq1 and Fastq2 or AlignedReads
#
# (other arguments are of course ok, but these are mandatory for the script to work


queue = 'dspencer'
group = '/dspencer/adhoc'

# arguments
parser = argparse.ArgumentParser(description='WDL launcher')
parser.add_argument("-w","--wdl",help="WDL file")
parser.add_argument("-i","--inputs",help="WDL file")
parser.add_argument("-n", "--name",help="Name of this process")
parser.add_argument("-t", "--testrun",help="Parse inputs and prepare to run, but dont make any output files or launch",action='store_true')
parser.add_argument('wdlvars',nargs='*',help="WDL inputs, e.g.: Key1=Value1 Key2=Value2. Can also accept these on stdin by entering '-'")

# parse arguments
args = parser.parse_args()

norun = False
if args.testrun:
    norun = True

if args.name is None:
    print('ERROR: Run name is required',file=sys.stderr)
    print()
    parser.print_help()
    sys.exit()
    
name = args.name

# if no workflow supplied
if args.wdl is None:
    parser.print_help()
    sys.exit()
    
wdl = str(args.wdl)


if os.path.exists(wdl) is False and os.path.exists('/storage1/fs1/dspencer/Active/spencerlab/workflow/' + wdl) is False:   
    print('ERROR: wdl file' + wdl + " does not exist",file=sys.stderr)
    print()
    parser.print_help()
    sys.exit()

if os.path.exists('/storage1/fs1/dspencer/Active/spencerlab/workflow/' + wdl) is True:
    wdl = '/storage1/fs1/dspencer/Active/spencerlab/workflow/' + wdl

    
# this is the input from the wdl. used to check that all the required ones are there.
wdlinput = json.loads(subprocess.run(['java','-jar',womtool,'inputs',wdl], stdout=subprocess.PIPE).stdout)
wdlname = [*wdlinput][0].split(".")[0]

if args.inputs:
    f = open(args.inputs)
    otherinput = json.load(f)
    for k in otherinput:
        if k in wdlinput.keys():
            wdlinput[k] = otherinput[k]
    f.close()
    
# get provided inputs
myargs = {}
for onearg in list(args.wdlvars):
    (k, v) = onearg.split("=")
    
    # add wdl name prefix, if necessary
    if re.search(wdlname,k) is None:
        k = wdlname + '.' + k

    myargs[k] = v

# iterate over all expected inputs
myinputs = {}
requiredinputs = []

for k in wdlinput:

    # if var was provided as input
    if k in myargs.keys():

        # get absolute path for files
        if re.search('File',str(wdlinput[k])) is not None:
            if os.path.exists(os.path.realpath(myargs[k])):
                myinputs[k] = os.path.realpath(myargs[k])
            else:
                print('ERROR: WDL input ' + k + ' is not a valid file',file=sys.stderr)

        else:
            myinputs[k] = myargs[k]

    elif re.search('optional',str(wdlinput[k])) is None and k not in myinputs.keys() and k not in myargs and wdlinput[k] == "":
            requiredinputs = requiredinputs + [k]

if len(requiredinputs) > 1:
    print('ERROR: WDL missing required inputs:\n' + "\n".join(requiredinputs),file=sys.stderr)
    print("\nAll inputs for WDL " + wdlname + ": \n",file=sys.stderr)    
    print(*[f"{': '.join(map(str,v))}" for i,v in enumerate(list(wdlinput.items()))], sep='\n',file=sys.stderr)
    sys.exit()

print("Launching " + wdlname + "...",file=sys.stderr)

with open(os.path.join(name + "." + wdlname + ".json"),'w') as jf:
    json.dump(myinputs, jf)
    jf.close()

err = os.path.join(name + "." + wdlname + ".err")
log = os.path.join(name + "." + wdlname + ".log")
inputs = os.path.join(name + "." + wdlname + ".json")

if norun is not True:
    os.system(f"LSF_DOCKER_VOLUMES=\"/scratch1/fs1/gtac-mgi/CLE:/scratch1/fs1/gtac-mgi/CLE /storage1/fs1/gtac-mgi/Active/CLE:/storage1/fs1/gtac-mgi/Active/CLE /storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active /storage1/fs1/timley/Active:/storage1/fs1/timley/Active /storage1/fs1/duncavagee/Active:/storage1/fs1/duncavagee/Active /storage1/fs1/danielclink/Active:/storage1/fs1/danielclink/Active /scratch1/fs1/dspencer:/scratch1/fs1/dspencer /scratch1/fs1/duncavagee:/scratch1/fs1/duncavagee /scratch1/fs1/danielclink:/scratch1/fs1/danielclink $HOME:$HOME\" bsub -eo {err} -oo {log} -q {queue} -G compute-dspencer -g {group} -R \"select[mem>=32000] span[hosts=1] rusage[mem=32000]\" -M 32000000 -a \"docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-20)\" /usr/bin/java -Dconfig.file={conf} -jar /opt/cromwell.jar run -t wdl -i {inputs} {wdl}")

print(f"LSF_DOCKER_VOLUMES=\"/scratch1/fs1/gtac-mgi/CLE:/scratch1/fs1/gtac-mgi/CLE /storage1/fs1/gtac-mgi/Active/CLE:/storage1/fs1/gtac-mgi/Active/CLE /storage1/fs1/dspencer/Active:/storage1/fs1/dspencer/Active /storage1/fs1/timley/Active:/storage1/fs1/timley/Active /storage1/fs1/duncavagee/Active:/storage1/fs1/duncavagee/Active /storage1/fs1/danielclink/Active:/storage1/fs1/danielclink/Active /scratch1/fs1/dspencer:/scratch1/fs1/dspencer /scratch1/fs1/duncavagee:/scratch1/fs1/duncavagee /scratch1/fs1/danielclink:/scratch1/fs1/danielclink $HOME:$HOME\" bsub -eo {err} -oo {log} -q {queue} -G compute-dspencer -g {group} -R \"select[mem>=32000] span[hosts=1] rusage[mem=32000]\" -M 32000000 -a \"docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-20)\" /usr/bin/java -Dconfig.file={conf} -jar /opt/cromwell.jar run -t wdl -i {inputs} {wdl}\n",file=sys.stderr)
