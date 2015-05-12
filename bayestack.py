#!/usr/bin/env python

"""
Execute this code as

mpirun -np 4 ./bayestack.py bayestack_settings.py

"""

import os,sys
import time,shutil,glob
import importlib
import pymultinest
from bayestackClasses import countModel
from utils import touch,remark,remarks,dump_variable_values
from mpi4py import MPI
import dill
MPI._p_pickle.dumps = dill.dumps
MPI._p_pickle.loads = dill.loads

__name_cached=__name__
if __name__=='__main__':
    param_file=sys.argv[-1]
    settingsf=param_file.split('.')[-2]
    set_module=importlib.import_module(settingsf)
    globals().update(set_module.__dict__)
    __name__=__name_cached

#-------------------------------------------------------------------------------

def main():

    settingsf=param_file.split('.')[-2]
    expt=countModel(modelFamily,nlaws,settingsf,dataset,binStyle,floatNoise)

    # Set up MPI
    world=MPI.COMM_WORLD
    rank=world.rank
    size=world.size
    master = rank==0

    if master:
        set_module=importlib.import_module(settingsf)
        globals().update(set_module.__dict__)

    note='MPI processors checked in: rank/size = (%i/%i)' % (rank,size)
    print note

    if master:
        try:
            os.mkdir(outdir)
            # Fix permissions
            os.chmod(outdir,0755)
        except OSError:
            pass

        logf=os.path.join(outdir,logfile)
        if master and os.path.exists(logf): os.remove(logf)
        log=open(logf,'w')
        remark(log,note)

    # Wait here after check-in...
    world.Barrier()
    if master: print 'All %i processors checked in...' % size

    # Broadcast global settings variables
    if master:
        set_dict = set_module.__dict__
    else:
        set_dict = None
    set_dict = world.bcast(set_dict,root=0)

    if not master:
        globals().update(set_dict)
        #print globals()

    # Wait here after broadcast...
    world.Barrier()
    if master: print 'All %i processors received OK...\n' % size

    # Write settings variables to file
    if master:
        variablesf=os.path.join(outdir,variablesfile)
        dump_variable_values(set_module,variablesf,verbose=False)

        startTime = time.strftime('%X %x %Z')
        shutil.copy(param_file,os.path.join(outdir,'bayestack_settings.py'))
        shutil.copy(datafile,outdir)
        notes=['Time now is %s' % startTime,\
               'Settings file: %s' % param_file,\
               'Data file: %s' % datafile]
        remarks(log,notes)

        # This is to allow import of settings from outdir
        # i.e. from outdir import * [or whatever]
        init_file='__init__.py'
        initf=os.path.join(outdir,init_file)
        touch(initf)

        notes=['Bins taken from %s' % datafile,\
               '# Bin occupancies [i uJy uJy field^-1]:']
        remarks(log,notes)
        for ibin in xrange(expt.nbins-1):
            try:
                line='%i %f %f %f'%(ibin+1,expt.bins[ibin],expt.bins[ibin+1],expt.data[ibin])
            except IndexError:
                print "Probably your binstyle doesn't match the datafile bins"
                sys.exit(0)
            remark(log,line)


    # Run MultiNest
    if master: t0 = time.time()
    try:
        pymultinest.run(expt.loglike,expt.logprior,expt.nparams,\
                    resume=RESUME,verbose=True,\
                    multimodal=multimodal,max_modes=max_modes,write_output=True,\
                    n_live_points=n_live_points,\
                    evidence_tolerance=evidence_tolerance,\
                    # mode_tolerance=-1e90 bugfix for earlier versions
                    # of PyMultiNest
                    mode_tolerance=-1e90,seed=SEED_SAMP,max_iter=max_iter,\
                    importance_nested_sampling=do_INS,\
                    outputfiles_basename=os.path.join(outdir,outstem),\
        # NB MPI is already init'ed by mpi4py (crashes otherwise)
                    init_MPI=False)
    except:
        return 1

    if master:
        stopTime=time.strftime('%X %x %Z')
        t1 = time.time()
        dt=t1-t0

        # Touch the output dir so Dropbox picks it up
        touch(outdir)
        
        notes=['Time then was %s' % startTime,\
               'Time now is %s' % stopTime,\
               'Execution took %6.4f sec (~ %i min) with %i cores' % \
                                         (dt,int(round(dt/60.0)),size),\
               'Arguments: %s' % ' '.join(sys.argv),\
               'INS   = %s' % do_INS,\
               'nlive = %i' % n_live_points,\
               'Run comment: %s' % comment,\
               'Now execute:',\
               '\n./plot.py %s' % outdir,\
               'and\n./reconstruct.py %s' % outdir]

        remarks(log,notes)
        log.close()

        print 'Parameters were:',expt.parameters
        
        # Copy the stats file so it's legible on my iPhone, Google, email etc.
        stats_dotdat= '%(od)s/%(os)sstats.dat' % {'od':outdir,'os':outstem}
        stats_dottxt= '%(od)s/%(os)sstats.txt' % {'od':outdir,'os':outstem}
        shutil.copy(stats_dotdat,stats_dottxt)

        # Now make all the files world readable
        globlist=glob.glob(os.path.join(outdir,'*'))
        print 'Globlist:'
        print globlist
        [os.chmod(f,644) for f in globlist]
        print 'End Globlist'

        print 'Run finished.'
        
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

