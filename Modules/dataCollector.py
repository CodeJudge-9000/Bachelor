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

def get_raw_data_frame(directory = './data'):
    theCwd = os.getcwd()
    os.chdir(directory)

    mypath = os.getcwd()
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and ".csv" in f]

    dataArray = []
    for file in onlyfiles:
        df = pd.read_csv(file)
        dataArray.append(df)

    dataArray2 = pd.concat(dataArray)
    os.chdir("..")
    dataArray2 = dataArray2.iloc[: , 1:]
    #dataArray2 = dataArray2.reset_index()


    df = dataArray2
    df = df.reset_index()
    return(df)


def fill_missing(df):
    for g in range(0,len(df)):
        h = df.loc[g].finger
        res = h.strip('][').split(',')
        listen = []
        for a in res:
            listen.append(int(a))

        while len(listen) < 5:
            listen.append(0)


        first = listen[0:1] 
        rest = listen[1:5]
        rest.sort(reverse = True)
        listen = first + rest

        df = df.replace(h, str(listen))
        
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
                #df = df.replace(h, idNew)
                df.new_id[g] = idNew

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

                #df = df.replace(h, idNew)
                df.new_id[g] = idNew

        cumList = []
        lens = []
        g = 0
        # Tilføj så vi også kan se z-coordinaten af den fjernede
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
                
        return(df)
    
def removed_categorize(df):
        # Tilføj så vi ser om de fjernede er same side eller opposite size
    zcoordMovedAtom = []
    zzListen = []

    for g in range(0,len(df)): #len(df)
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

    # Section for seeing if the removed are the same or different side
    categoryVector = []
    for i, zvec in enumerate(zzListen):
        categoryList = []

        for zpos in zvec:
            if(zpos == []):
                categoryList.append(" ")
            elif(zpos > 8.5 and zcoordMovedAtom[i] > 8.5): ## Hvis begge er på toppen
                categoryList.append("S")
            elif(zpos < 8.5 and zcoordMovedAtom[i] < 8.5): ## Hvis begge er på toppen
                categoryList.append("S")
            elif(zpos > 8.5 and zcoordMovedAtom[i] < 8.5): ## Hvis den ene er top og den anden er bund
                categoryList.append("D")
            elif(zpos < 8.5 and zcoordMovedAtom[i] > 8.5): ## Hvis den ene er bund og den anden er top
                categoryList.append("D")


        categoryVector.append(categoryList)    

    ## Inversion section
    #### For dem med original id zcoordinat over 8.5 Skal der ikke ske noget, for dem hvor den er under, skal de invertes.   
    df['neigbohrs'] = categoryVector
    return(df)


def finger_print_search(df, fingerPrint = []):
    ## Kigger i denne section på fingerprints som er lignene
    ins = []
    for i in range(0,len(df)):
        ins.append( df['shortFinger'][i] == fingerPrint)

    df2 = df[ins]
    return(df2)


# Similar system, similar fingerprint  ## SORTING SO WE ONLY See the 6x6 results
def remove_small_system(df2):
    df3 = df2[df2['startsystemname'] == 'eq_struct_6x6.traj']
    df3
    df3 = df3.reset_index()
    return(df3)



def get_trajectories_from_table(df):
    # LAD OS EXTRAHERER VORES TING FRA TABELLEN
    Trajs = []
    for g in range(0,len(df)):
        filename = df.loc[g].trajFile
        if(filename[-6] == "_"):
            filename = filename[:-6] + ".traj"
        os.chdir("./trajec")
        T = Trajectory(filename)
        Trajs.append(T)
        os.chdir("..")