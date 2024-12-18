//
// Created by Nikhil Chandra Admal on 11/19/24.
//

#ifndef OILAB_LAMMPS_H
#define OILAB_LAMMPS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <Eigen/Eigen>
#include <iomanip>

std::tuple<Eigen::MatrixXd,
           Eigen::Matrix3d,
           Eigen::Vector3d> read_oILAB_output(const std::string& path)
{
    Eigen::MatrixXd atoms;
    Eigen::Matrix3d box;
    box.setZero();
    Eigen::Vector<double,9> boxVectorized;
    boxVectorized.setZero();
    Eigen::Vector3d origin;
    origin.setZero();

    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error opening file for reading oilab config file: " << path << std::endl;
    }

    std::string line;
    int number_atoms = 0;

    int lineCount= 0;
    int atomCount= 0;
    int boxCount= -1;
    int originCount= -1;

    while (std::getline(file, line)) {
        lineCount++;
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;

        while (ss >> field) {
            fields.push_back(field);
        }

        if (lineCount == 1) {
            number_atoms = std::stoi(fields[0]);
            atoms.resize(number_atoms,5);
            atoms.setZero();
        }
        else if (lineCount == 2)
        {
            for (size_t i = 0; i < fields.size(); ++i) {
                if (fields[i].find("Lattice=") != std::string::npos) {
                    boxCount++;
                    continue;
                }
                if (boxCount>=0 && boxCount<9) {
                    boxVectorized(boxCount) = std::stod(fields[i]);
                    boxCount++;
                    continue;
                }
                if (fields[i].find("origin=") != std::string::npos)
                {
                    originCount++;
                    continue;
                }
                if (originCount>-1 && originCount<3) {
                    origin(originCount) = std::stod(fields[i]);
                    originCount++;
                    continue;
                }
            }
        }
        else if (lineCount > 2) {
            // Reading atom data (type, x, y, z, radius)
            Eigen::VectorXd atom_data(5);
            int coordCount= 0;
            for (const auto &field : fields) {
                atom_data(coordCount)= std::stod(field);
                coordCount++;
            }
            atoms.row(atomCount)= atom_data;
            atomCount++;
        }
    }
    box= boxVectorized.reshaped(3,3);

    return {atoms, box, origin};

}

/*-----------------------*/
//std::string write_lammps_datafile(const std::string &filename, const std::vector<std::vector<double>> &box, const std::vector<std::vector<double>> &atom_data, int atom_types) {
void write_lammps_datafile(const std::string &filename,
                           const std::vector<std::vector<double>> &box,
                           const Eigen::MatrixXd& atom_data,
                           int atom_types)
{
    std::ofstream file(filename);
    if (!file.is_open())
        std::cerr << "Error opening file for writing lammps configuration file: " << filename << std::endl;

    file << "# LAMMPS data file via write_data\n";
    file << "\n";
    file << atom_data.rows() << " atoms\n";
    file << atom_types << " atom types\n";
    file << "\n";
    file << std::setprecision(8) << box[0][0] << " " << std::setprecision(8) << box[0][1] << " xlo xhi\n";
    file << std::setprecision(8) << box[1][0] << " " << std::setprecision(8) << box[1][1] << " ylo yhi\n";
    file << std::setprecision(8) << box[2][0] << " " << std::setprecision(8) << box[2][1] << " zlo zhi\n";
    file << "0 0 0 xy xz yz\n";
    file << "\n";
    file << "Atoms # atomic\n";
    file << "\n";
    for (size_t i = 0; i < atom_data.rows(); ++i) {
        //file << atom_data(i,0) << " " << atom_data(i,1) << " " << atom_data(i,2) << " " << atom_data(i,3) << " " << atom_data(i,4) << " 0 0 0\n";
        file << std::setprecision(8) << atom_data.row(i) << " 0 0 0\n";
    }
}

/*-----------------------*/
void write_lammps_input_script(const std::string &filename, const std::string &infile, const std::string &outfile, double cohesive_energy, double gb_thickness_parameter, const std::string &potential_file_path, const std::string &output_dump_file) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for writing lammps input script: " << filename << std::endl;
        return;
    }

    file << "# Find energy of a config and write in " << outfile << " file\n";
    file << "\n";
    file << "variable potential_path string " << potential_file_path << "\n";
    file << "clear\n";
    file << "units metal\n";
    file << "boundary p p p\n";
    file << "atom_style atomic\n";
    file << "neighbor 1.0 bin\n";
    file << "neigh_modify every 1 delay 2 check yes\n";
    file << "read_data " << infile << "\n";
    file << "pair_style      eam/alloy\n";
    file << "pair_coeff      * * ${potential_path} Cu Cu\n";
    file << "delete_atoms overlap 1e-2 all all\n";
    file << "variable area equal ly*lz\n";
    file << "variable        xlogb equal xlo+" << std::setprecision(8) << gb_thickness_parameter << "\n";
    file << "variable        xhigb equal xhi-" << std::setprecision(8) << gb_thickness_parameter << "\n";
    file << "variable        xlobulk equal xlo+" << std::setprecision(8) << gb_thickness_parameter << "\n";
    file << "variable        xhibulk equal xlo+" << std::setprecision(8) << 1.5*gb_thickness_parameter << "\n";
    file << "region         GB      block ${xlogb} ${xhigb} INF INF INF INF side in units box\n";
    file << "region         BULK    block ${xlobulk} ${xhibulk} INF INF INF INF side in units box\n";
    file << "group           GB region GB\n";
    file << "group           BULK region BULK\n";
    file << "compute         peratom GB pe/atom\n";
    file << "compute         peratombulk BULK pe/atom\n";
    file << "compute         pe GB reduce sum c_peratom\n";
    file << "compute         pebulk GB reduce sum c_peratombulk\n";
    file << "variable        peGB equal c_pe\n";
    file << "variable        peBULK equal c_pebulk\n";
    file << "variable        atomsGB equal count(GB)\n";
    file << "variable        atomsBULK equal count(BULK)\n";
    file << "thermo 1000\n";
    file << "compute 1 all ke/atom\n";
    file << "compute cna all cna/atom 3.08133\n";
    file << "compute csys all centro/atom  fcc\n";
    file << "compute 3 all pe/atom\n";
    file << "compute 4 all stress/atom NULL pair\n";
    file << "timestep        0.001\n";
    file << "thermo_style custom step temp ke pe etotal press pxx pyy pzz pxy pxz pyz ly lx lz xy xz yz c_pe v_atomsGB v_peBULK v_atomsBULK\n";
    file << "dump                    OUT0 all custom 10 " << output_dump_file << " id type x y z fx fy fz c_3 c_1 vx vy vz c_4[1] c_4[2] c_4[3] c_4[4] c_4[5] c_4[6]\n";
    file << "run                     0\n";
    file << "variable        coh equal (${peBULK}/${atomsBULK})\n";
    file << "variable        GBene equal (${peGB}-${coh}*${atomsGB})\n";
    file << "print \"coh = ${coh} energy = ${peGB} numAtoms = ${atomsGB} GBene = ${GBene} area = ${area}\" file " << outfile << "\n";
    file << "\n";
}

/*-----------------------*/
std::vector<std::vector<double>> read_python_outfile(const std::string &path) {
    std::vector<std::vector<double>> data;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error opening file for reading lammps output: " << path << std::endl;
        return data;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        double state_id, area, total_energy, gb_energy, gb_density;

        while (ss >> field)
        {
            if (field == "coh") {
                ss >> field;
                ss >> state_id;
            }
            else if (field == "energy") {
                ss >> field;
                ss >> total_energy;
            }
            else if (field == "GBene") {
                ss >> field;
                ss >> gb_energy;
            }
            else if (field == "numAtoms") {
                ss >> field;
                ss >> gb_density;
            }
            else if (field == "area")
            {
                ss >> field;
                ss >> area;
                break;
            }
        }

        //gb_density = gb_density / area;
        data.push_back({state_id, area, gb_energy, gb_density});
    }

    return data;
}


std::pair<double, double> energy(const std::string& lammpsLocation,
                                 const std::string& oilabConfigFile,
                                 const std::string& potentialFile)
{
    // Write data
    std::string threadNumber= std::to_string(omp_get_thread_num());
    std::string lammpsInputFile= "in"+  threadNumber +".find_energy";
    std::string lammpsDataFile= "data" + threadNumber + ".lammps_input";
    std::string lammpsDumpFile= "dump" + threadNumber + ".lammpsConfigs";
    std::string outfile = "lmp_mesostate_energies" + threadNumber + ".txt";

    //double cohesive_energy = -3.539999890784346;
    double cohesive_energy = -3.54022;

    // Read data
    auto [atoms, box, origin] = read_oILAB_output(oilabConfigFile);

    // Find rotation and new box
    Eigen::Matrix3d R= box.transpose();
    for(int i=0; i<3; ++i)
        R.row(i).normalize();
    Eigen::Matrix3d new_box= R*box;

    // New atoms
    Eigen::MatrixXd new_atoms(atoms.rows(),5);
    for (size_t i = 0; i < atoms.rows(); ++i) {
        new_atoms(i,Eigen::seq(2,4))= R*(atoms(i,Eigen::seq(1,3)).transpose());
        /*
        new_atoms[i][2] = new_atoms[i][2] + R[0][0]*atoms[i][1] + R[0][1]*atoms[i][2] + R[0][2]*atoms[i][3];
        new_atoms[i][3] = new_atoms[i][3] + R[1][0]*atoms[i][1] + R[1][1]*atoms[i][2] + R[1][2]*atoms[i][3];
        new_atoms[i][4] = new_atoms[i][4] + R[2][0]*atoms[i][1] + R[2][1]*atoms[i][2] + R[2][2]*atoms[i][3];
         */
        new_atoms(i,0)= i + 1;
        new_atoms(i,1) = atoms(i,0);
    }

    // Determine nbox values
    std::vector<std::vector<double>> nbox(3, std::vector<double>(2, 0.0));
    for (size_t i = 0; i < 3; ++i) {
        if (i == 0) {
            nbox[i][0] = -new_box(i,i) / 2;
            nbox[i][1] = new_box(i,i) / 2;
        } else if (i == 1) {
            nbox[i][0] = origin(i);
            nbox[i][1] = new_box(i,i);
        } else {
            nbox[i][0] = 0;
            nbox[i][1] = new_box(i,i);
        }
    }

    //double gb_thickness_parameter = 0.75 * (new_box(0,0) * 0.5);
    double gb_thickness_parameter= 6;

    // Write files
    write_lammps_datafile(lammpsDataFile, nbox, new_atoms, 2);
    write_lammps_input_script(lammpsInputFile, lammpsDataFile, outfile, cohesive_energy, gb_thickness_parameter, potentialFile, lammpsDumpFile);

    // Run the LAMMPS script
    std::string command = lammpsLocation +" -in " + lammpsInputFile + " > /dev/null 2>&1";
    std::system(command.c_str());

    // Read energy
    auto data_energy = read_python_outfile(outfile);

    return {data_energy[0][3], data_energy[0][2]};
}

#endif //OILAB_LAMMPS_H
