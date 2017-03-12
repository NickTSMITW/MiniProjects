import matplotlib.pyplot as plt
import numpy as np
import os

path = "/home/cdt1602/Documents/MiniProjects/Preconditioners/SCF_script"

i = 0 #dummy

dirListing = os.listdir(path)
file_list = [] #array of strings of filenames

for item in dirListing:
    if ".cell" in item:
        file_list.append(item) #adds the filenames ending in .cell to the array file_list
        i = i + 1


for j in range(0,i):

    file_list[j] = file_list[j].replace(' ', '')[:-5] #gets rid of the .cell in the filename string

    #os.system("rm " + file_list[j] + ".castep") #removes current filename.castep 

    #runs current filename.castep
    #os.system("mpirun -np 4 /home/cdt1602/Documents/CASTEP_source_edit/CASTEP-16.11/obj/linux_x86_64_gfortran4.8/castep.mpi " + file_list[j])


#this loop uses string "Final energy" in the .castep to locate the linenumber of the final SCF calcultion
#then prints the first 3 characters of the line of the final SCF cycle -- which is the number of SCF cycles
#to reach convergence 
for j in range(0,i):
    with open(file_list[j] + ".castep") as dummy_file:
        for line_num, line in enumerate(dummy_file, 1):
            if "Final energy, E" in line:
                line_num = line_num - 4
                file = open(file_list[j] + ".castep")
                SCF_cycle_num = file.readlines()
                SCF_cycle_num[line_num] = SCF_cycle_num[line_num].strip()[:3]

                SCF_file = open("SCF_cycles.dat","a")
                SCF_file.write(SCF_cycle_num[line_num] + "\n")

