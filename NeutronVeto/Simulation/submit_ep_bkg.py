#!/apps/bin/python

from subprocess import Popen, PIPE
from numpy import linspace

for m in range(1,100): #100
    command="""#!/bin/sh
#SBATCH --job-name=ep_bkg_{0}
#SBATCH --output=/farm_out/erins/ep_bkg_{0}_output.out
#SBATCH --error=/farm_out/erins/ep_bkg_{0}_error.err
#SBATCH --partition=production
#SBATCH --account=clas
#SBATCH --mem-per-cpu=2048
#SBATCH -t 1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=farm1968,farm160133
bash ./add_ep_bkg.sh {0}
""".format(str(m))
    p=Popen(args=["sbatch"],stdin=PIPE);
    p.communicate(command)
