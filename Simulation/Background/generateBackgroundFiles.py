import os
import glob
import subprocess
import random


class slurmWrapper:
    def __init__(self):
        self.header = """sbatch << EOF
#!/bin/bash
#SBATCH --cpus-per-task=1\n
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000                                                                                            
#SBATCH --account=clas12                                                                                             
#SBATCH --job-name=mc_rgm_gcf                                                                                              
#SBATCH --partition=production                                                               
#SBATCH --time=20:00:00                                                                                               
#SBATCH --output=/farm_out/%u/%x-%j-%N.out                                                                           
#SBATCH --error=//farm_out/%u/%x-%j-%N.err                                                                           
        """

    def addCommand(self,command):
        self.header += "\n" + command
    
    def printCommand(self):
        print(self.header)

    def submit(self):
        self.header += "EOF"
        print("Submitting %s" % self.printCommand())
        os.system(self.header)


def find_sub_dirs(path):
    paths = []
    directories = glob.glob(path)
    print("Looking in path %s" % path)
    directories = [int(sub[len(sub) -5 :]) for sub in directories if sub[len(sub) -5 :].isdigit() ]
    return directories

os.system('source ../../environment.csh')

random.seed(12345)

decoded_path='/cache/clas12/rg-m/production/decoded/v10.0.1/'
output_path='/volatile/clas12/rg-m/bkgfiles/'
target_events = 3000000 # target number of events
trigger_bit = 39 #random trigger bit FC

energy_setting = {"6gev":5986, "4gev":4029,"2gev":2070}
min_currents   = {"6gev":20, "4gev":-1,"2gev":-1}
#min_currents   = {"6gev":10, "4gev":5,"2gev":0.5}
mag_setting    = {"6gev":"tor-1.00_sol-1.00","4gev":"tor-1.00_sol-1.00","2gev":"tor+0.50_sol-1.00"}
estimates      = {"6gev":5000, "4gev":2000,"2gev":1000}

#print(glob.glob("/cache/clas12/rg-m/production/pass1/**/**/*", recursive = True))
cooking_paths = glob.glob("/cache/clas12/rg-m/production/pass1/**/**/")
print(cooking_paths)
for cooking_path in cooking_paths:
    cooking_path = cooking_path + "dst/recon/*"

    energy_mev = [val for key,val in energy_setting.items() if key in cooking_path][0]
    setting  = [val for key,val in mag_setting.items() if key in cooking_path][0]
    energy   = [key for key,val in mag_setting.items() if key in cooking_path][0]
    target   = cooking_path.split("/")[cooking_path.split("/").index(str(energy)) + 1]
    min_current = [val for key,val in min_currents.items() if key in cooking_path][0]
    avg_bkg     = [val for key,val in estimates.items() if key in cooking_path][0]
    

    if energy != "6gev":
        continue

    print("")
    print("Beam energy is %s Target %s Magnetic Setting %s Min Current %f" % (energy,target,setting,min_current))
    
    output_dir = "%srgm_fall2021/%s/%s_%sMeV/current/" % (output_path,setting,target,str(energy_mev))
    
    print("New DIRECTORY %s " % output_dir)
    if not(os.path.isdir(output_dir)):
        command = "mkdir -p " + output_dir
        os.system(command)
        
        
    # get cooked runs from directory 
    cooked_runs = find_sub_dirs(cooking_path)
    #print(cooked_runs)
    paths = [decoded_path + '0' + str(run) for run in cooked_runs]
    
    files = [file for file in [glob.glob(path + '/*.hipo') for path in paths]]
    # make into one big list of files
    files = sum(files,[])
    
    est_total = len(files)*avg_bkg #estimated total number of events
    sample_size = int(target_events/avg_bkg) #number of files we want to process
    
    print("Sampling %i files to meet %i events" % (sample_size,target_events))
    print("Estimated total BKG events in full data set %i" % est_total)
    print("")
    
    #randomly sample files
    #files = random.sample(files,2) #for testing 
    if sample_size < len(files):
        files = random.sample(files,sample_size) 

    print(len(files))
    print(sample_size)

    for file in files:
        slurm = slurmWrapper()
        slurm.addCommand("source ../environment_gemc.sh")
        
#        print("Processing %s" % file)
#        print(file.split("/"))
#        print(file.split("/")[-1])
        outputfile_name = output_dir + "bkg_" + file.split("/")[-1]
#        print("Output files %s" % outputfile_name)
        
        command = "trigger-filter -b %i -c %f %s -o %s \n" % (trigger_bit,min_current,file,outputfile_name)
        #    slurm.addCommand(command)
        
        #    slurm.submit()


