from trace import Trace
import pandas as pd
import os
from os import listdir
from os.path import isfile, join
from tabulate import tabulate
from math import isnan 
from ase.io import read
from FingerPrints import *
from atomic_annihilator import remove_atoms
from ase.io import Trajectory

def get_raw_data_frame(directory = './data') -> pd.DataFrame:
    """ get_raw_data_frame - goes through the specified directory and collects the raw data into a singular dataframe
    It goes for csv files and collect them to a dataFrame.

    Parameters:
    -----------
     - `directory` [str]: The directory where we look for the csv files we want to recombine

    Returns:
    --------
     - `df` [pd.DataFrame]: A complete dataframe
    """

    # Going down to the folder and extract all filenames
    os.chdir(directory)
    mypath = os.getcwd()
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and ".csv" in f]

    # Reading all the files and add them to a dataframe array
    dataArray = []
    for file in onlyfiles:
        df = pd.read_csv(file)
        dataArray.append(df)
    os.chdir("..")

    #refining the output
    dataArray2 = pd.concat(dataArray)
    dataArray2 = dataArray2.iloc[: , 1:]
    df = dataArray2
    df = df.reset_index()
    return df


def fill_missing(df) -> pd.DataFrame:
    """ fill_missing - A function which takes a dataframe and cleans old datapoints such that missing new_id's after removals are set if they where nan
    fixes so older simple fingers goes from length to 4 to all 5
    fixes so all instances will have a SCF

    Parameters:
    -----------
     - `df` [pd.Dataframe]: the dataframe we want to clean
    
    Returns:
    --------
     - `df` [pd.Dataframe]: the cleaned dataframe

    """

    #Here we loop through all fingerprints and clean the simple fingerprints which was 4 long, but needed to be 5
    for g in range(0,len(df)):
        h = df.loc[g].finger
        res = h.strip('][').split(',')
        listen = []
        for a in res:
            listen.append(int(a))

        while len(listen) < 5:
            listen.append(0)


        first = listen[:1] 
        rest = listen[1:]
        rest.sort(reverse = True)
        listen = first + rest

        df = df.replace(h, str(listen))
    
    #Here we loop through all atoms and makes sure that we replace new id, for them which was a NAN
    for g in range(0,len(df)):
        h = df.loc[g].new_id
        if isnan(h) and  df.loc[g].startsystemname == "eq_struct_6x6.traj":
            system = read("eq_struct_6x6.traj")
            ogId = df.loc[g].original_id
            listt = df.loc[g].removedList
            res = listt.strip('][').split(',')
            listen = []
            for a in res:
                listen.append(int(a))

            idNew = remove_atoms(system = system, atomIndex = ogId, atomRemoveIndex = listen, relax = False)
            df.new_id[g] = idNew
        

    #A lot of the data, did not have SCF, so this section is to add it to the table for the points which are missing
    cumList = []
    lens = []
    for g in range(0,len(df)):
        newId = df.loc[g].new_id
        newId = int(newId)
        filename = df.loc[g].trajFile
        if(filename[-6] == "_"):
            filename = filename[:-6] + ".traj"
        os.chdir("./trajec")
        T = Trajectory(filename)
        lens.append(len(T))
        T = T[0]
        cumFinger = ShortCommonFinger(system = T, index = newId)
        cumList.append(cumFinger)
        os.chdir("..")
        df['shortFinger'] = cumList
        df['Frames'] = lens
                
        return df

def removed_categorize(df) -> pd.DataFrame:
    """ removed_categorize - A function that takes a pandas dataframe, with our model and checks which of the neighbohrs we have removed was on the same side 
    or the different side

    Parameters:
    -----------
     - `df` [pd.DataFrame]: The data under consideration

    Returns:
    --------
     - `df` [pd.DataFrame]: with the added [D,S,D,...] vector
    """

    # A loop for getting the coordinate of all the removed atoms, so we can check them later on
    zcoordMovedAtom = []
    zzListen = []
    for g in range(0,len(df)):
            system = read("eq_struct_6x6.traj")

            ogId = df.loc[g].original_id
            zcoordMovedAtom.append(system[ogId].position[2])

            listt = df.loc[g].removedList
            res = listt.strip('][').split(',')
            listen = []

            for a in res:
                if(a == ""):
                    listen.append([])
                else:
                    listen.append(int(a))

            zListen = []
            for d in listen:
                if d == []:
                    zListen.append([])
                else:
                    zListen.append(system[d].position[2])

            zzListen.append(zListen)   

    # Section for assigning if the removed are the same or different side
    categoryVector = []
    for i, zvec in enumerate(zzListen):
        categoryList = []

        for zpos in zvec:
            if(zpos == []):
                categoryList.append(" ")
            elif(zpos > 8.5 and zcoordMovedAtom[i] > 8.5): ## Hvis begge er på toppen
                categoryList.append("S")
            elif(zpos < 8.5 and zcoordMovedAtom[i] < 8.5): ## Hvis begge er i bunden
                categoryList.append("S")
            elif(zpos > 8.5 and zcoordMovedAtom[i] < 8.5): ## Hvis den ene er top og den anden er bund
                categoryList.append("D")
            elif(zpos < 8.5 and zcoordMovedAtom[i] > 8.5): ## Hvis den ene er bund og den anden er top
                categoryList.append("D")

        categoryVector.append(categoryList)    

    df['neigbohrs'] = categoryVector
    return df


def finger_print_search(df, fingerPrint = []) -> pd.DataFrame:
    """ finger_print_search - scrolls through all the instances in the pandas dataFrame 

    Parameters:
    -----------
     - `df` [pd.DataFrame]: the dataframe we want to search for
     - `fingerPrint` [list]: the fingerprint we want to search for
    
    Returns:
    --------
     - `df2` [pd.DataFrame]: A dataframe with the entries that match the fingerPrint
    """

    # Getting the indeces where the shortFinger match
    ins = []
    for i in range(0,len(df)):
        ins.append(df['shortFinger'][i] == fingerPrint)


    df2 = df[ins]
    try:
        df2 = df2.drop(['level_0'], axis = 1)
    except:
        pass
    df2 = df2.reset_index() 
    return df2


# Similar system, similar fingerprint  ## SORTING SO WE ONLY See the 6x6 results
def remove_small_system(df2) -> pd.DataFrame:
    """ remove_small_system - takes out dataframe and sort out all the data which involves 5x5 struct because they are unreliable

    Parameters:
    -----------
     - `df2` [pd.DataFrame]: The dataframe we want to sort in

    Returns:
    --------
     - `df3` [pd.Dataframe]: The outsorted dataframe
    """

    df3 = df2[df2['startsystemname'] == 'eq_struct_6x6.traj']

    # Refreshing indexing
    try:
        df3 = df3.drop(['level_0'], axis = 1)
    except:
        pass
    df3 = df3.reset_index()
    return df3


def get_trajectories_from_table(df, trajectoryFolder = "./trajec") -> list(Trajectory):
    """ get_trajectories_from_table - takes a dataframe and goes into the subfolder ./trajec and extract the trajectories

    Parameters:
    -----------
     - `df` [pd.Dataframe]:
     - `trajectoryFolder` [str]: the subfolder where the trajectory is placed

    Returns:
    --------
     - `Trajs` [list(Trajectory)]: A list with alle the trajectories, so we can have easy acces to things we want to look up
    """
    # go into the trajectoryfolder and extract all the filename with their trajectories
    os.chdir(trajectoryFolder)
    Trajs = []
    for g in range(0,len(df)):
        filename = df.loc[g].trajFile
        # if the filename from the csv file is od from the trajectory filename, we remove _ add a .traj
        if(filename[-6] == "_"):
            filename = filename[:-6] + ".traj"
        
        T = Trajectory(filename)
        Trajs.append(T)
    #Go back so we do not have to deal with that
    os.chdir("..")
    
    return Trajs


def sort_out_id(df, index) -> pd.DataFrame:
    """ sort_out_id - Just a regular sorter, where we look at the original id and sort out and also refreshes the table
     - `df` [dataframe]: dataframe we want to sort in
     - `index` [int]: the originale index of the atom we want to sort out
    """

    viewTable = df[df['original_id'] == index].sort_values('removedList', ascending=False)
    viewTable = viewTable.drop(['level_0'], axis = 1)
    viewTable = viewTable.reset_index()
    return viewTable

def nn_search(df, neigbohrsCategory = [])  -> pd.DataFrame:
    """ nn_search - takes a dataframe and neigbohrsCategorizaiton to sort out for example so we only get ['S'] or ['S','D']    
    """
    # Here we just finds the neigbohrs which match the category we are interested in
    ins = []
    for i in range(0,len(df)):
        ins.append( df['neigbohrs'][i] == neigbohrsCategory)
    df2 = df[ins]

    # if dataframe is not clean in indexing we clean it
    try:
        df2 = df2.drop(['level_0'], axis = 1)
    except:
        pass
    df2 = df2.reset_index() 

    return df2

def get_same_siders(df)  -> pd.DataFrame:
    """ get_same_siders - takes a dataframe and looks so we get out instances from the table where ['S'] and ['S', 'S'] and ['S', 'S', 'S'] instances

    Input and Output:
    -----------------
     - `df` [pd.Dataframe]: the table we want to sort in and the table where we have sorted.
    """

    # Get 3 nn_searces and conjugate them
    df1 =  nn_search(df, ['S'])
    df2 =  nn_search(df, ['S','S']) 
    df3 =  nn_search(df, ['S','S','S'])
    dff = pd.concat([df1, df2, df3])

    # Clean the indexing of the new dataframe
    try:
        dff = dff.drop(['level_0'], axis = 1)
    except:
        pass
    dff = dff.reset_index()
    return dff
    