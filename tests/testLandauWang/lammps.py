#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 08:40:19 2024

@author: hj-home

This script reads the output from oILAB and gives out the energy of the system

"""
import os
import subprocess
import numpy as np
import numpy.linalg as la
import sys

"""
def read_oILAB_output(path):
    # Reads data from oILAB putput file 
    # Return a list consisting of : [Array of atom positions (type, x,y,z,radius)
    #                               Array of box limits (3x3)
    #                               Array of box origin (1x3)
    #                               Boundary conditions
    #                               Other details in xyz file]
    sys.stdout = open(os.devnull, "w")
    sys.stderr = open(os.devnull, "w")

    atoms = []
    message = "Reading " + path
    print(message)
    with open(path,'r') as file:
        for line in file:
            fields = line.split(' ')
            if len(fields) == 1:
                number_atoms = int(fields[0])
            elif len(fields) >= 10 and fields[0].split("\"")[0] == "Lattice=":
                details = fields
                count_lattice = pbc_count = flag_pbc = flag_origin = lattice_flag = 0
                pbc = []
                origin = []
                lattice = []
                for j in range(len(details)):
                    if details[j] != "":
                        if count_lattice == 0 and details[j].split("\"")[0] == "Lattice=":
                            #lattice.append(float(details[j].split("\"")[1]))
                            #count_lattice += 1 
                            lattice_flag = 1
                        if count_lattice<9 and lattice_flag == 1 and details[j].split("\"")[0] != "Lattice=":
                            lattice.append(float(details[j].strip("\"")))
                            count_lattice += 1
                            
                        if len(details[j])>30:
                            properties = details[j]
                            
                        if details[j].split("\"")[0] == "PBC=" and pbc_count < 3:
                            pbc.append(details[j][-1])
                            pbc_count += 1
                            flag_pbc =1
                        if flag_pbc == 1 and pbc_count <3 and details[j].split("\"")[0] != "PBC=":
                            pbc.append(details[j].strip("\""))
                            pbc_count += 1
                            if pbc_count >2:
                                flag_pbc = 0
                            
                        if details[j].split("\"")[0] == "origin=":
                            #origin.append(float(details[j].split("\"")[1]))
                            flag_origin = 1
                        if flag_origin == 1 and details[j].split("\"")[0] != "origin=":
                            origin.append(float(details[j].split("\"")[0]))
                        
            if len(fields) >4 and fields[0].split("\"")[0] != "Lattice=":
                row = []
                for j in range(len(fields)):
                    if fields[j] != "":
                        row.append(float(fields[j]))
                atoms.append(row)
    return [np.array(atoms),np.reshape(lattice,(3,3)),np.array(origin),pbc,properties]
"""

def read_oILAB_output(path):
    # Reads data from oILAB putput file
    # Return a list consisting of : [Array of atom positions (type, x,y,z,radius)
    #                               Array of box limits (3x3)
    #                               Array of box origin (1x3)
    #                               Boundary conditions
    #                               Other details in xyz file]
    sys.stdout = open(os.devnull, "w")
    sys.stderr = open(os.devnull, "w")
    atoms = []
    message = "Reading " + path
    print(message)
    with open(path,'r') as file:
        for line in file:
            fields = line.split(' ')
            if len(fields) == 1:
                number_atoms = int(fields[0])
            elif len(fields) >= 10 and fields[0].split("\"")[0] == "Lattice=":
                details = fields
                count_lattice = pbc_count = flag_pbc = flag_origin = lattice_flag = 0
                pbc = []
                origin = []
                lattice = []
                lattice_formatting_flag = 0
                for j in range(len(details)):
                    if details[j] != "":
                        # Diagnostic section to check if the Lattice information is formatted the way we want
                        if details[j].split("\"")[0] == "Lattice=" and details[j].split("\"")[1] == "":
                            lattice_formatting_flag = 1
                        if count_lattice == 0 and details[j].split("\"")[0] == "Lattice=":
                            if lattice_formatting_flag == 0:
                                lattice.append(float(details[j].split("\"")[1]))
                                count_lattice += 1
                            lattice_flag = 1
                        if count_lattice<9 and lattice_flag == 1 and details[j].split("\"")[0] != "Lattice=":
                            lattice.append(float(details[j].strip("\"")))
                            count_lattice += 1
                        if len(details[j])>30:
                            properties = details[j]
                        if details[j].split("\"")[0] == "PBC=" and pbc_count < 3:
                            pbc.append(details[j][-1])
                            pbc_count += 1
                            flag_pbc =1
                        if flag_pbc == 1 and pbc_count <3 and details[j].split("\"")[0] != "PBC=":
                            pbc.append(details[j].strip("\""))
                            pbc_count += 1
                            if pbc_count >2:
                                flag_pbc = 0
                        if details[j].split("\"")[0] == "origin=":
                            #origin.append(float(details[j].split("\"")[1]))
                            flag_origin = 1
                        if flag_origin == 1 and details[j].split("\"")[0] != "origin=":
                            origin.append(float(details[j].split("\"")[0]))
            if len(fields) >4 and fields[0].split("\"")[0] != "Lattice=":
                row = []
                for j in range(len(fields)):
                    if fields[j] != "":
                        row.append(float(fields[j]))
                atoms.append(row)
    return [np.array(atoms),np.reshape(lattice,(3,3)),np.array(origin),pbc,properties]


def align_box_lammps(box1):
    # For LAMMPS input , the first dimension vector needs to be aligned to [100], the second to [010] and the third to [001]
    # This function find the appropriate rotations and changes the box
    lammps_basis = np.eye(3)
    lattice_basis = np.zeros((3,3))
    for i in range(len(lattice_basis)):
        lattice_basis[:,i] = box1[i,:]/la.norm(box1[i,:])
    R = np.matmul(lammps_basis,lattice_basis.T)
    
    return R

def write_lammps_datafile(outfolder,file,box,atom_data,atom_types):

    file = outfolder + file
    message = "Writing " + file
    print(message)
    f = open(file,"w")
    f.write("#LAMMPS data file via write_data, contains minimized configurations selected, units = metal\n")
    f.write("\n")
    f.write("%d atoms\n"%len(atom_data))
    f.write("%d atom types\n"%atom_types)
    f.write("\n")
    f.write("%f %f xlo xhi\n"%(box[0,0],box[0,1]))
    f.write("%f %f ylo yhi\n"%(box[1,0],box[1,1]))
    f.write("%f %f zlo zhi\n"%(box[2,0],box[2,1]))
    f.write("0 0 0 xy xz yz\n")
    f.write("\n")
    f.write("Atoms # atomic\n")
    f.write("\n")
    for i in range(len(atom_data)):
        f.write("%d %d %f %f %f 0 0 0\n"%(atom_data[i,0],atom_data[i,1],atom_data[i,2],atom_data[i,3],atom_data[i,4]))

    return file


def write_lammps_input_script(folder,file,infile,outfile,lattice_constant,cohesive_energy,gb_thickness_parameter,potential_file_path,output_dump_file):
    
    
    f = open(folder+file,"w")
    
    f.write("#Find energy of a config and write in %s file\n"%(outfile))
    f.write("\n")
    f.write("variable a equal 0\n")
    f.write("variable a0 equal %f\n"%lattice_constant)
    f.write("variable coh equal %f\n"%cohesive_energy)
    f.write("variable potential_path string %s\n"%potential_file_path)
    f.write("clear\n")
    f.write("units metal\n")
    f.write("boundary p p p\n")
    f.write("atom_style atomic\n")
    f.write("neighbor 1.0 bin\n")
    f.write("neigh_modify every 1 delay 2 check yes\n")
    f.write("read_data %s\n"%infile)
    f.write("pair_style      eam/alloy\n")
    f.write("pair_coeff      * * ${potential_path} Cu Cu\n")
    f.write("change_box      all x scale ${a0} y scale ${a0} z scale ${a0} remap units box\n")
    f.write("variable        xlobox equal xlo-10\n")
    f.write("variable        xhibox equal xhi-10\n")
    f.write("delete_atoms overlap 1e-2 all all\n")
    f.write("variable area equal ly*lz\n")
    f.write("variable        xlogb equal xlo+%d\n"%gb_thickness_parameter)
    f.write("variable        xhigb equal xhi-%d\n"%gb_thickness_parameter)
    f.write("region         GB      block ${xlogb} ${xhigb} INF INF INF INF side in units box\n")
    f.write("group           GB region GB\n")
    f.write("compute         peratom GB pe/atom\n")
    f.write("compute         pe GB reduce sum c_peratom\n")
    f.write("variable        peGB equal c_pe\n")
    f.write("variable        atomsGB equal count(GB)\n")
    f.write("\n")
    f.write("thermo 1000\n")
    f.write("compute 1 all ke/atom\n")
    f.write("compute cna all cna/atom 3.08133\n")
    f.write("compute csys all centro/atom  fcc\n")
    f.write("compute 3 all pe/atom\n")
    f.write("compute 4 all stress/atom NULL pair\n")
    f.write("timestep        0.001\n")
    f.write("thermo_style custom step temp ke pe etotal press pxx pyy pzz pxy pxz pyz ly lx lz xy xz yz c_pe v_atomsGB\n")
    f.write("dump                    OUT0 all custom 10 %s id type x y z fx fy fz c_3 c_1 vx vy vz c_4[1] c_4[2] c_4[3] c_4[4] c_4[5] c_4[6]\n"%output_dump_file)
    f.write("run                     0\n")
    f.write("variable        GBene equal (${peGB}-${coh}*${atomsGB})/${area}*16.02\n")
    f.write("print \"a = ${a} energy = ${peGB} numAtoms = ${atomsGB} GBene = ${GBene} area = ${area}\" file %s\n"%outfile)
    f.write("\n")
    f.write("\n")
    
    f.close()

    
    return f

def read_python_outfile(path):
    # Reading the outfile that this script created containing the energy of the configuration under consideration
    data_python = []
    message = "Reading " + path
    print(message)
    with open(path,'r') as file:
        for line in file:
            fields = line.split(' ')
            for i in range(len(fields)):
                if fields[i]=="a":
                    state_id = float(fields[i+2])
                elif fields[i] == "energy":
                    total_energy = float(fields[i+2])
                elif fields[i] == "numAtoms":
                    natoms_gb = float(fields[i+2])
                elif fields[i] == "GBene":
                    gb_energy = float(fields[i+2])
                elif fields[i] == "area":
                    area = float(fields[i+2].strip("\n"))
                    break
            gb_density = natoms_gb/area
            data_python.append([state_id,area,gb_energy,gb_density])
    return data_python

def energy():
    '''
    input_folder = "/Users/hj-home/Desktop/Research/testGbMesoStates/sigma3_itb/sigma3/"
    out_folder = "/Users/hj-home/Desktop/Research/testGbMesoStates/sigma3_itb/sigma3_python/"
    potential_file_path = "/Users/hj-home/Desktop/Research/testGbMesoStates/"
    '''
    input_folder = os.getcwd()+"/"
    # Folder where you want the output to be stored (and it serves as the folder where the temporary files are stored)
    out_folder = os.getcwd()+"/"
    # Location of potential file to be used
    potential_file_path = os.getcwd()+"/Cu_mishin1.eam.alloy"
    # Lammps object location
    lammps_loc = "/Users/Nikhil/Documents/Academic/Software/lammps-15May15/src/"
    # Lattice constant of the material you want to simulate
    lattice_constant = 3.615
    # Cohesive energy corresponsing to the potential file
    cohesive_energy = -3.54
    # Location of mpi object
    mpi_loc = "/opt/homebrew/bin/"

    lattice_constant = 3.615
    cohesive_energy = -3.54

    # Read data
    file = "temp_reference1.txt"
    data = read_oILAB_output(input_folder+file)
    box = data[1]
    atoms = data[0]
    origin = data[2]
    print(box)
    # Find rotation and new box
    R = align_box_lammps(box)
    new_box = np.zeros((3,3))
    for i in range(len(new_box)):
        new_box[i,:] = np.matmul(R,box[i,:])
    new_atoms = np.zeros((len(atoms),5))
    for i in range(len(atoms)):
        new_atoms[i,2:5] = np.matmul(R,atoms[i,1:4])
        new_atoms[i,0] = i+1
        new_atoms[i,1] = atoms[i,0]
    
    
    nbox = np.zeros((3,2))
    for i in range(len(nbox)):
        if i == 0:
            nbox[i,0] = -new_box[i,i]/2
            nbox[i,1] = new_box[i,i]/2
        elif i == 1:
            nbox[i,0] = origin[i]
            nbox[i,1] = new_box[i,i]
        else:
            nbox[i,0] = 0
            nbox[i,1] = new_box[i,i]
    
    # Write data
    os.makedirs(out_folder,exist_ok = True)
    file = "data.lammps_input"
    f =  write_lammps_datafile(out_folder, file, nbox, new_atoms, 2)
    
    # Write lammps input script
    lammps_file = "in.find_energy"
    outfile = out_folder + "lmp_mesostate_energies.txt"
    gb_thickness_parameter = 0.75*(new_box[0,0]*lattice_constant/2)
    output_dump_file = out_folder+"dump.lammpsConfigs"
    lf = write_lammps_input_script(out_folder,lammps_file,f,outfile,lattice_constant,cohesive_energy,gb_thickness_parameter,potential_file_path,output_dump_file)
    
    # Run the lammps script
    multicore = False # Change this to run the lammps script on single or multiple cores
    lammps_object_path = lammps_loc
    if multicore == False:
        command = lammps_object_path+"lmp_serial -in " + out_folder+lammps_file
    else:
        command = mpi_loc+"mpirun -np 6 "+lammps_loc+"lmp_mpi -in " + out_folder+lammps_file
    subprocess.run([command],shell=True,stdout = subprocess.DEVNULL) 
    #subprocess.run([command],shell=True)
    
    # Read energy
    data_energy = read_python_outfile(outfile)
    
    return data_energy[0][2], data_energy[0][3]
print(energy())


