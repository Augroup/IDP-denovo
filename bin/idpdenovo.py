#!/usr/bin/env python
import argparse, sys, os
from shutil import rmtree, copy
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

def main(args):
  # Back up environment
  _environ = dict(os.environ)

  curr_dir = os.getcwd()
  sys.stderr.write("starting in directory: "+curr_dir+"\n")
  exec_dir = os.path.dirname(os.path.realpath(__file__))
  sys.stderr.write("executing in directory: "+exec_dir+"\n")
  args.tempdir = os.path.abspath(args.tempdir)
  sys.stderr.write("working directory: "+args.tempdir+"\n")

  sys.stderr.write("building the configuration file\n")
  build_cfg(args)

  os.chdir(args.tempdir) # move into wthe working directory

  sys.stderr.write("storing environment variables temporarily\n")
  os.environ['PATH'] = exec_dir+'/../utilities' +':'+ os.environ['PATH']
  os.environ['PATH'] = exec_dir+'/../plugins/BIN' +':'+ os.environ['PATH']
  os.environ['PATH'] = exec_dir+'/../plugins/BIN/bin' +':'+ os.environ['PATH']
  os.environ['PATH'] = exec_dir+'/../plugins/IDP_0.1.9/bin' +':'+ os.environ['PATH']

  if args.test:
    # Run the test command
    cmd = exec_dir+'/pre-test ./config_file'
    sys.stderr.write("executing: "+cmd+"\n")
    p = Popen(cmd.split(),stderr=sys.stderr,stdout=sys.stdout)
    p.communicate()
    return
  else:
    # Run the actual command
    cmd = exec_dir+'/../bin/idpdenovo_core ./config_file'
    sys.stderr.write("executing: "+cmd+"\n")
    p = Popen(cmd.split(),stderr=sys.stderr)
    p.communicate()
  
  os.chdir(curr_dir) # move into the working directory we originally were in
  if not os.path.exists(args.output):
    os.makedirs(args.output)

  sys.stderr.write("Saving outputs\n")
  save_results(args)


  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)
  # Restore environment
  sys.stderr.write("restoring environment variables\n")
  os.environ.clear()
  os.environ.update(_environ)

def save_results(args):
  copy(args.tempdir+'/idpdenovo_output/combine_seq',args.output+'/')
  copy(args.tempdir+'/idpdenovo_output/seq_cluster',args.output+'/')
  copy(args.tempdir+'/idpdenovo_output/report_seq',args.output+'/')
  copy(args.tempdir+'/idpdenovo_output/report.gpd',args.output+'/')
  copy(args.tempdir+'/idpdenovo_output/confirmed_gap',args.output+'/')
  copy(args.tempdir+'/idpdenovo_output/lr_quantify',args.output+'/')
  copy(args.tempdir+'/idpdenovo_output/sr_quantify',args.output+'/')
  copy(args.tempdir+'/idpdenovo_output/lr_input',args.output+'/')

def build_cfg(args):
  of = open(args.tempdir+'/config_file','w')
  of.write('scaffold='+os.path.abspath(args.SR_scaffold)+"\n")
  of.write('lr='+os.path.abspath(args.long_reads)+"\n")
  of.write('leftSr='+os.path.abspath(args.SR_left)+"\n")
  of.write('rightSr='+os.path.abspath(args.SR_right)+"\n")
  of.write('kmer='+str(args.k_mer_length)+"\n")
  of.write('kCutoff='+str(args.k_mer_frequency)+"\n")
  of.write('nThreads='+str(args.threads)+"\n")
  of.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="IDP-denovo",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('SR_scaffold',help="Short read scaffold",type=lambda x: check_file(parser,x))
  parser.add_argument('long_reads',help="Long reads in FASTA format",type=lambda x: check_file(parser,x))
  parser.add_argument('SR_left',help="1s Left mate short reads reads in FASTA format",type=lambda x: check_file(parser,x))
  parser.add_argument('SR_right',help="2s Right mate short reads reads in FASTA format",type=lambda x: check_file(parser,x))
  parser.add_argument('-k','--k_mer_length',default=15,type=int,help="k-mer length used in clustering of unaligned long reads")
  parser.add_argument('-f','--k_mer_frequency',default=0.05,type=float,help="k-mer frequency cutoff")
  parser.add_argument('-o','--output',required=True,help="REQUIRED Output directory")
  parser.add_argument('-t','--threads',type=int,default=1,help="INT number of threads to run.")
  # Hidden arguments for internal testing
  parser.add_argument('--test',action='store_true',help=argparse.SUPPRESS)

  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()

  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

def check_file(parser,fname):
  if not os.path.exists(fname):
    parser.error("file "+fname+" not found\n")
  else:
    return fname

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="idpdenovouser.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return 

if __name__=="__main__":
  #do our inputs
  args = do_inputs()
  main(args)
