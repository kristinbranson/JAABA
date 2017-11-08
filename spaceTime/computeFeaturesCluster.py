#!/usr/bin/python
import os,sys, getopt
import math
import pwd
import time

def main(argv):
  if len(argv)!=6:
    print "Usage: python computeFeaturesCluster.py expdir stationary nframes outdir recompute method"
    print('Method is one of deep-sup or hs-sup')

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
  method = argv[5]
  username = pwd.getpwuid( os.getuid() ).pw_name
  if expdir[-1]=='/':
    expdir = expdir[0:-1]

  (temp,expname) = os.path.split(expdir)
  if not os.path.isdir(outdir):
	  os.makedirs(outdir)

  for curb in range(1,nblocks+1):
    curname = expname + '_' + str(curb)
    curoutname = outdir + '/' + expname
    curjob = 'export MCR_CACHE_ROOT="/scratch/' +  username + '/mcr_cache_root.' + curname + '";\n'
    dirname = os.path.dirname(os.path.realpath(__file__))
    fullscriptname = os.path.join(dirname,'run_computeFeaturesCompiled.sh')
    curjob = curjob + fullscriptname + ' ' + mcrpath + ' ' + moviename + ' '+ trxname + ' ' + stationary + ' ' + method + ' '+  str(curb) +' ' + str(blocksize) +' ' + curoutname + ';\n'
    outfile = curoutname + '_' + str(curb) + '_log.txt'
    shfile = curoutname +  '_' + str(curb) + '_script.sh'
    outmatfile = curoutname + '_' + str(curb) + '.mat'
    if (not recompute) & os.path.isfile(outmatfile):
      continue
    f = open(shfile,'w')
    f.write(curjob)
    f.close()
    timeoutpath = os.path.join(dirname,'mytimeout')
    cmd = "qsub -pe batch 4 -N deepFtrs" + curname + " -j y -o " + outfile +  ''' -b y -cwd -V ' ''' + timeoutpath + ''' -m 30000000 "bash ''' +  shfile +  " > " +  outfile +  '''.out" ' '''
    print cmd
    time.sleep(0.5)
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

