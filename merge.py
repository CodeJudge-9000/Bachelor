import pandas as pd
import os
from os import listdir
from os.path import isfile, join

def mergeTDs(mergeName = "defaultMerge", subFolder = ""):
    # Get current working directory
    cwd = os.getcwd()
    
    # Go to subfolder relative path, and update current working directory
    os.chdir(f"{cwd}{subFolder}")
    cwd = f"{cwd}{subFolder}"
    
    # Find all files in folder
    allFiles = [f for f in os.listdir(cwd) if isfile(join(cwd, f))]
    
    # Remove all filenames that are not .csv files
    allFiles = [f for f in allFiles if os.path.splitext(f)[1] == ".csv"]
    
    # For each file load using pandas, and concatenate everything into one object
    pandaList = []
    for f in allFiles:
        tempPanda = pd.read_csv(f)
        pandaList.append(tempPanda)
    
    pandaResult = pd.concat(pandaList)
    
    # Now save this to a single file
    pandaResult.to_csv(f"{mergeName}.csv")
    
    return

