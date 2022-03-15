# Module import
import os
      
class bsubmissions:
    """ A class for submitting jobs into the lsf 10 cluster dtu uses. All variables need to be filled otherwise a job will fail.
    We set jobName, which appears in the console of the jobs
    coresAsked -> How many cores we want to use
    scriptName -> What we will call the script we want to submit
    mem -> The memory we wants
    maxMem -> Maximum memory to use before we let program stop
    wallTime -> Max time out program should run for
    """
    def __init__(self, jobName, coresAsked, scriptName, numHosts, email, mem, maxMem, wallTime):
        # Set variables
        self.currentFolder = os.getcwd()
        self.jobName = jobName
        self.coresAsked = coresAsked
        self.scriptName = scriptName
        self.numHosts = numHosts
        self.email = email
        self.mem = mem
        self.maxMem = maxMem
        self.wallTime = wallTime
    
    # make error paths
    def make_subfolders(self):
        os.system('mkdir errors')
        os.system('mkdir outputs')
        os.system('mkdir data')


    def make_bsub_file(self):
        """ 
        A bsub file maker, ready to go, fills in all the necesisties and creates two subfolders for errors and output files. 
        """
        # Making subfolders for error
        self.make_subfolders()

        errorPath = self.currentFolder + '/errors/'
        outputPath = self.currentFolder + '/outputs/'

        bsubString = f"""
        #!/bin/sh                                                                          
        ### General options                                                                           
        ### -- specify queue --                                                                         
        #BSUB -q hpc                                                                                   
        ### -- set the job Name --                                                                     
        #BSUB -J {self.jobName}                                                                              
        ### -- ask for number of cores (default: 1) --                                                 
        #BSUB -n {self.coresAsked}                                                                                     
        ### -- specify that the cores must be on the same host --                                                            
        #BSUB -R "span[hosts={self.numHosts}]"                                                                               
        ### -- specify that we need 2GB of memory per core/slot --                                                                                  
        #BSUB -R "rusage[mem={self.mem}]"                                                                                 
        ### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot --                                                                             
        #BSUB -M {self.maxMem}                           
        ### -- set walltime limit: hh:mm --                            
        #BSUB -W {self.wallTime}                            
        ### -- set the email address --                                                       
        # please uncomment the following line and put in your e-mail address,                           
        # if you want to receive e-mail notifications on a non-default address                           
        ###BSUB -u {self.email}                                                    
        ### -- send notification at start --                                                       
        #BSUB -B                                                                                                            
        ### -- send notification at completion --                                                                                  
        #BSUB -N                                                                                                             
        ### -- Specify the output and error file. %J is the job-id -- 
        ### -- -o and -e mean append, -oo and -eo mean overwrite -- 
        #BSUB -o {outputPath}Output_%J.out 
        #BSUB -e {errorPath}Error_%J.err 


        # here follow the commands you want to execute
        cd /zhome/06/4/136859/bachelor
        source venv/bin/activate

        cd {self.currentFolder}
        mpiexec -n {self.coresAsked} python3 {self.scriptName}
        """
        shFile = open(f'{self.jobName}.sh','w')
        shFile.write(bsubString)


    def delete_bsub_file(self):
        os.remove(f'{self.jobName}.sh')


    def do_submission(self):
        os.system(f'bsub < {self.jobName}.sh')

def make_list_name(theList):
    theString = ""
    for i in theList:
        theString = theString + str(i)
    return theString