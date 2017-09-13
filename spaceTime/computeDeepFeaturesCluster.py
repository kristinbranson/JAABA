#!/usr/bin/python
import os,sys, getopt
import math
import pwd

def main(argv):
  if len(argv)!=5:
    print "Usage: python computeDeepFeaturesCluster.py expdir stationary nframes outdir recompute"
    return;
  mcrpath = '/groups/branson/bransonlab/mayank/MCR/v85'
  moviename = argv[0] + '/movie.ufmf'
  trxname = argv[0] + '/registered_trx.mat'
  stationary = argv[1]
  outdir = argv[3]
  blocksize = 500
  nframes = int(argv[2])
  nblocks = int(math.ceil(float(nframes)/blocksize))
  expdir = argv[0]
  recompute = int(argv[4])
  username = pwd.getpwuid( os.getuid() ).pw_name
  if expdir[-1]=='/':
    expdir = expdir[0:-1]

  (temp,expname) = os.path.split(expdir)

  for curb in range(1,nblocks+1):
    curname = expname + '_' + str(curb)
    curoutname = outdir + '/' + expname
    curjob = 'export MCR_CACHE_ROOT="/scratch/' +  username + '/mcr_cache_root.' + curname + '";\n'
    curjob = curjob + './run_computeFeaturesCompiled.sh' + ' ' + mcrpath + ' ' + moviename + ' '+ trxname + ' ' + stationary + ' deep-sup ' + str(curb) +' ' + str(blocksize) +' ' + curoutname + ';\n'
    outfile = curoutname + '_' + str(curb) + '_log.txt'
    shfile = curoutname +  '_' + str(curb) + '_script.sh'
    outmatfile = curoutname + '_' + str(curb) + '.mat'
    if (not recompute) & os.path.isfile(outmatfile):
      continue
    f = open(shfile,'w')
    f.write(curjob)
    f.close()
    cmd = "qsub -pe batch 1 -N deepFtrs" + curname + " -j y -o " + outfile +  " -b y -cwd -V 'bash " +  shfile +  " > " +  outfile +  ".out '"
    print cmd
    retval = os.system(cmd)


#    curjob = curjob + "bash /groups/branson/home/kabram/JAABA/perframe/run_StrawmenInfluenceCluster_par.sh /groups/branson/bransonlab/projects/olympiad/MCR/v717 " + curstr1
#    cwd = "/groups/branson/home/kabram/JAABA/perframe"
#    outfile = cwd + "/loaded_influence/outfiles/" + curstr + ".log"
#    shfile = cwd + '/loaded_influence/scripts/' + curstr + ".sh"
#    f = open(shfile,'w')
#    f.write(curjob)
#    f.close()
#    cmd = "qsub -pe batch 1 -l new=true -N strawmen" + curstr + " -j y -o " + outfile +  " -b y -cwd -V 'bash " +  shfile +  " > " +  outfile +  ".out '"
#    print cmd
#    retval = os.system(cmd)


if __name__ == "__main__":
  main(sys.argv[1:])

