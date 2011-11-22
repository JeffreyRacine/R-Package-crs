# This R profile can be used when a cluster does not allow spawning or a job 
# scheduler is required to launch any parallel jobs. Saving this file as 
# .Rprofile in the working directory or root directory. For unix platform, run
# mpirexec -n [cpu numbers] R --no-save -q
# For windows platform with mpich2, use mpiexec wrapper and specify a working 
# directory where .Rprofile is inside.
# Cannot be used as Rprofile.site because it will not work

# Following system libraries are not loaded automatically. So manual loads are 
# needed.
library(utils)
library(stats)
library(datasets)
library(grDevices)
library(graphics)
library(methods)

if (!invisible(library(npRmpi,logical.return = TRUE))){
    warning("npRmpi cannot be loaded")
    q(save = "no")
}

options(error=quote(assign(".mpi.err", FALSE, env = .GlobalEnv)))

if (mpi.comm.size(0) > 1)
    invisible(mpi.comm.dup(0,1))

if (mpi.comm.rank(0) >0){
    #sys.load.image(".RData",TRUE)
    options(echo=FALSE)
    .comm <- 1
    mpi.barrier(0)
    repeat 
	try(eval(mpi.bcast.cmd(rank=0,comm=.comm)),TRUE)
	#try(eval(mpi.bcast.cmd(rank=0,comm=.comm),env=sys.parent()),TRUE)
    #mpi.barrier(.comm)
    if (is.loaded("mpi_comm_disconnect")) 
        mpi.comm.disconnect(.comm)
    else mpi.comm.free(.comm)
    mpi.quit()
}
    
if (mpi.comm.rank(0)==0) {
    #options(echo=TRUE)
    mpi.barrier(0)
    if(mpi.comm.size(0) > 1)
        slave.hostinfo(1)
}

.Last <- function(){
    if (is.loaded("mpi_initialize")){
        if (mpi.comm.size(1) > 1){
            print("Please use mpi.close.Rslaves() to close slaves")
            mpi.close.Rslaves(comm=1)
    	}
    }
    print("Please use mpi.quit() to quit R")
    mpi.quit()
}
