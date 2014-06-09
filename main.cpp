#include "pdbhelper.h"

#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <memory>
#include <math.h>
#include <algorithm>
#include <typeinfo>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <sstream>

std::vector<PDBAtom> extractAtomsFromXYZ(std::string filename) {

	std::vector<PDBAtom> atoms;
	std::ifstream xyz_file(filename.c_str());
	if (xyz_file.is_open() && xyz_file.good()) {

		// Read for 2 data lines.
		std::string comment_line;
		int number_of_atoms;
		xyz_file >> number_of_atoms;
		std::getline(xyz_file, comment_line);


		std::string atomic_symbol, x, y, z;

		for (int i = 0; i < number_of_atoms; ++i) {
			
			xyz_file >> atomic_symbol >> x >> y >> z;
			atoms.push_back(PDBAtom(atomic_symbol, atof(x.c_str()), atof(y.c_str()), atof(z.c_str())));
		}

	}
	xyz_file.close();

	return atoms;
}

int main(int argc, char **argv) {
	
	int c;
	std::string amber_file_s, pdb_file_s, xyz_file_s, output_pdb_file_s;
	while ( (c = getopt(argc, argv, "a:p:x:o:")) != -1 ) {

		switch(c) {
			case 'a':
				amber_file_s.assign(optarg);
				break;
			case 'p':
				pdb_file_s.assign(optarg);
				break;
			case 'x':
				xyz_file_s.assign(optarg);
				break;
			case 'o':
				output_pdb_file_s.assign(optarg);
				break;
			default:
				exit(EXIT_FAILURE);
		}
	}

	if (amber_file_s.size() == 0) {
		std::cout << "Amber topology file required. Use -a <file>." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (pdb_file_s.size() == 0) {
		std::cout << "PDB file required. Use -p <file>." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (xyz_file_s.size() == 0) {
		std::cout << "Input XYZ file is required. Use -x <command>." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (output_pdb_file_s.size() == 0) {
		std::cout << "Output pdb file is required. Use -o <command>." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::vector<PDBAtom> new_atoms = extractAtomsFromXYZ(xyz_file_s);

	PDBHelper pdb_helper(pdb_file_s, amber_file_s, output_pdb_file_s, std::vector<std::string>());
	pdb_helper.updateAllAtomsFromXYZ(new_atoms);

	// std::cout << output_pdb_file_s << std::endl;
	pdb_helper.writePDBFile(output_pdb_file_s);
	// pdb_helper.writePDBFile("output.pdb");

	// pdb_helper->updateAllNonEXAFSAtomsFromXYZ(dcd_helper.getXYZAtFrame(0));
	// pdb_helper->writePDBFile();
	
	// std::cout << "Energy: " << vmd_helper->calculateEnergy() << std::endl;
/*
	Genfig fitness_config(fitness_file);
	Genfig ga_config(ga_file);

	if (ga_config.hasKey("seed")) {
		
		seed = ga_config.getString("seed");
		srand(ga_config.getInt("seed"));
		std::cout << "Seed: " << ga_config.getInt("seed") << std::endl;
	} else if (seed.size() != 0) {

		srand(atoi(seed.c_str()));
		std::cout << "Seed: " << atoi(seed.c_str()) << std::endl;
	} else {

		time_t initial_seed = time(NULL);
		srand(initial_seed);

		std::stringstream ss;
		ss << initial_seed;
		seed = ss.str();
		
		std::cout << "Seed: " << initial_seed << std::endl;
	}

	// Check if scratch directory exists
	struct stat sb;
	if (stat(fitness_config.getString("folder-name").c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
		std::cout << "Scratch directory is missing." << std::endl;
		return 0;
	}

	// Create directory for scratch stuff.
	std::string temp_folder = fitness_config.getString("folder-name") + "/" + seed;

	if (mkdir(temp_folder.c_str(), 0755) != 0){
		do {
			temp_folder = temp_folder + "d";
		} while(mkdir(temp_folder.c_str(), 0755) != 0);
	}

	std::string temp_pdb = "temp_pdb.pdb";
	PDBHelper* pdb_helper = new PDBHelper(fitness_config.getString("pdb-file"), fitness_config.getString("amber-topology-file"), temp_folder + "/" + temp_pdb, fitness_config.getStringList("exafs-atoms"));

	IFEFFITHelper* ifeffit_helper = new IFEFFITHelper(temp_folder, pdb_helper->getAllEXAFSAtoms(), fitness_config.getString("target-atom"), fitness_config.getString("experimental-exafs"), fitness_config.getDouble("x-min"), fitness_config.getDouble("x-max"), fitness_config.getString("feff"), fitness_config.getString("ifeffit"));
	// IFEFFITHelper* ifeffit_helper = new IFEFFITHelper(fitness_config.getString("folder-name"), pdb_helper->getAllEXAFSAtoms(), fitness_config.getString("target-atom"), fitness_config.getString("experimental-exafs"), fitness_config.getDouble("x-min"), fitness_config.getDouble("x-max"), fitness_config.getString("feff"), fitness_config.getString("ifeffit"));

	VMDHelper* vmd_helper = new VMDHelper(temp_folder, temp_pdb, fitness_config.getString("amber-topology-file"), fitness_config.getString("namd2-path"), fitness_config.getString("vmd-path"));
	EXAFSEvaluator* exafs_evaluator = new EXAFSEvaluator(ifeffit_helper, pdb_helper, vmd_helper);

	if (ga_config.getString("eval-type").compare("solo") == 0) {
		std::cout << "Solo" << std::endl;
		std::cout << "Initial pdb file" << std::endl;
		double initial_rmsd = ifeffit_helper->run(pdb_helper->getEXAFSAtoms(), true);
		std::cout << "RMSD = " << initial_rmsd << std::endl;

		// DCDHelper dcd_helper = DCDHelper(ga_config.getString("dcd-file"));
		// pdb_helper->updateAllNonEXAFSAtomsFromXYZ(dcd_helper.getXYZAtFrame(0));
		pdb_helper->writePDBFile();
		
		std::cout << "Energy: " << vmd_helper->calculateEnergy() << std::endl;

	} else if (ga_config.getString("eval-type").compare("dir_potential") == 0) {
		std::cout << "Directory Potential" << std::endl;
		std::cout << "Reading directory: " << ga_config.getString("pot-dir") << std::endl;

		std::vector<std::string> files;
		getdir(ga_config.getString("pot-dir"),files);

		DCDHelper dcd_helper = DCDHelper(ga_config.getString("dcd-file"));

		for (std::vector<std::string>::iterator file = files.begin(); file != files.end(); ++file) {
			
			std::string pdb_file = ga_config.getString("pot-dir") + "/" + *file;
			std::cout << "Analyzing :" << pdb_file << std::endl;

			PDBHelper temp_pdb_helper = PDBHelper(pdb_file, fitness_config.getString("amber-topology-file"), temp_folder + "/" + temp_pdb, fitness_config.getStringList("exafs-atoms"));
			temp_pdb_helper.updateAllNonEXAFSAtomsFromXYZ(dcd_helper.getXYZAtFrame(0));
			temp_pdb_helper.writePDBFile();
			std::cout << "Energy: " << vmd_helper->calculateEnergy() << std::endl;
		}

	} else if (ga_config.getString("eval-type").compare("ga") == 0) {
		std::cout << "GA" << std::endl;

		EXAFSGA ga(exafs_evaluator, ga_config.getDouble("mutation"), ga_config.getDouble("crossover"), ga_config.getBool("elitism"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::cout << "Getting initial population." << std::endl;

		std::vector< std::vector<PDBAtom> > initial_population;
		switch (ga_config.getInt("population-type")) {
			case 0: {
				std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZs(ga_config.getString("dcd-file"), ga_config.getInt("population-size"));
				for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

					pdb_helper->updateEXAFSAtomsFromXYZ(*i);
					initial_population.push_back( pdb_helper->getEXAFSAtoms() );
				}
				break;
			}
			case 1: {
				std::vector< std::vector<PDBAtom> > initial_exafs_population = random_population(pdb_helper->getEXAFSAtoms(), ga_config.getInt("population-size"), 0.05);
				for (std::vector< std::vector<PDBAtom> >::iterator i = initial_exafs_population.begin(); i != initial_exafs_population.end(); ++i) {
					initial_population.push_back( *i );
				}
				break;
			}
		}

		std::cout << "GA: Begin" << std::endl;

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			initial_populations.push_back(initial_population);
		}

		ga.begin(initial_populations);

	} else if (ga_config.getString("eval-type").compare("ga_recenter") == 0) {
		std::cout << "GA Recentering" << std::endl;

		EXAFSGA ga(exafs_evaluator, ga_config.getDouble("mutation"), ga_config.getDouble("crossover"), ga_config.getBool("elitism"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::cout << "Getting initial population." << std::endl;

		std::vector< std::vector<PDBAtom> > initial_population;
		switch (ga_config.getInt("population-type")) {
			case 0: {
				std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZs(ga_config.getString("dcd-file"), ga_config.getInt("recentering-population"));
				for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

					pdb_helper->updateEXAFSAtomsFromXYZ(*i);
					initial_population.push_back( pdb_helper->getEXAFSAtoms() );
				}
				break;
			}
			case 1: {
				std::vector< std::vector<PDBAtom> > initial_exafs_population = random_population(pdb_helper->getEXAFSAtoms(), ga_config.getInt("population-size"), 0.05);
				for (std::vector< std::vector<PDBAtom> >::iterator i = initial_exafs_population.begin(); i != initial_exafs_population.end(); ++i) {
					initial_population.push_back( *i );
				}
				break;
			}
		}

		std::cout << "GA: Begin" << std::endl;

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			initial_populations.push_back(initial_population);
		}

		ga.begin_recentering(initial_populations, ga_config.getInt("population-size"), ga_config.getDouble("convergence-rate"), ga_config.getInt("recentering"));

	} else if (ga_config.getString("eval-type").compare("xyz") == 0) {
		std::cout << "XYZ" << std::endl;

		std::vector<std::string> xyz_files = ga_config.getStringList("xyz-files");

		for (std::vector<std::string>::iterator file = xyz_files.begin(); file != xyz_files.end(); ++file) {
			
			std::cout << "File: " << *file << std::endl;
			std::vector<PDBAtom> atoms = extractAtomsFromXYZ(*file);
			IFEFFITHelper ifeffit_helper = IFEFFITHelper(fitness_config.getString("folder-name"), atoms, fitness_config.getString("target-atom"), fitness_config.getString("experimental-exafs"), fitness_config.getDouble("x-min"), fitness_config.getDouble("x-max"), fitness_config.getString("feff"), fitness_config.getString("ifeffit"));

			double rmsd = ifeffit_helper.run(atoms, true);
			std::cout << "RMSD = " << rmsd << std::endl;

			std::vector< std::pair<double, double> > target_exafs = ifeffit_helper.getTargetEXAFS();
			std::vector< std::pair<double, double> > calculated_exafs = ifeffit_helper.getEXAFSData();

			std::ofstream exafs_output((*file + ".csv").c_str());
			for (int i = 0; i < (int)target_exafs.size() && i < (int)calculated_exafs.size() - 2; ++i) {
				
				exafs_output << target_exafs[i].first << "," << target_exafs[i].second << "," << calculated_exafs[i+1].second << std::endl;
			}
			exafs_output.close();

		}
	} else if (ga_config.getString("eval-type").compare("indexes") == 0) {

		std::cout << "Indexes" << std::endl;

		std::vector<int> indexes = dcdGetIndexLessThan(ga_config.getString("index-file"), ga_config.getDouble("index-max"));

	} else if (ga_config.getString("eval-type").compare("index_ga") == 0) {

		std::cout << "Index_GA" << std::endl;

		std::vector<int> indexes = dcdGetIndexLessThan(ga_config.getString("index-file"), ga_config.getDouble("index-max"));
		std::vector< std::vector<PDBAtom> > initial_population;

		EXAFSGA ga(exafs_evaluator, ga_config.getDouble("mutation"), ga_config.getDouble("crossover"), ga_config.getBool("elitism"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZsByIndex(ga_config.getString("dcd-file"), indexes);
		std::cout << "DCD Population = " << initial_dcd_population.size() << std::endl;
		for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

			pdb_helper->updateEXAFSAtomsFromXYZ(*i);
			initial_population.push_back( pdb_helper->getEXAFSAtoms() );
		}

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			
			// Get subset of population.
			std::random_shuffle(initial_population.begin(), initial_population.end());
			std::vector< std::vector<PDBAtom> > subset_pop = std::vector< std::vector<PDBAtom> >(initial_population.begin(), initial_population.begin()+ga_config.getInt("population-size"));

			initial_populations.push_back(subset_pop);
		}

		DCDHelper dcd_helper = DCDHelper(ga_config.getString("dcd-file"));
		pdb_helper->updateAllNonEXAFSAtomsFromXYZ(dcd_helper.getXYZAtFrame(0));

		std::cout << "GA: Begin" << std::endl;
		ga.begin(initial_populations);

	} else if (ga_config.getString("eval-type").compare("index_ga_recenter") == 0) {

		std::cout << "Index_GA_Recenter" << std::endl;

		std::vector<int> indexes = dcdGetIndexLessThan(ga_config.getString("index-file"), ga_config.getDouble("index-max"));
		std::vector< std::vector<PDBAtom> > initial_population;

		EXAFSGA ga(exafs_evaluator, ga_config.getDouble("mutation"), ga_config.getDouble("crossover"), ga_config.getBool("elitism"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZsByIndex(ga_config.getString("dcd-file"), indexes);
		std::cout << "DCD Population = " << initial_dcd_population.size() << std::endl;
		for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

			pdb_helper->updateEXAFSAtomsFromXYZ(*i);
			initial_population.push_back( pdb_helper->getEXAFSAtoms() );
		}

		std::cout << "GA: Begin" << std::endl;

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			
			// Get subset of population.
			std::random_shuffle(initial_population.begin(), initial_population.end());
			std::vector< std::vector<PDBAtom> > subset_pop = std::vector< std::vector<PDBAtom> >(initial_population.begin(), initial_population.begin()+ga_config.getInt("recentering-population"));

			initial_populations.push_back(subset_pop);
		}
		
		ga.begin_recentering(initial_populations, ga_config.getInt("population-size"), ga_config.getDouble("convergence-rate"), ga_config.getInt("recentering"));

	} else if (ga_config.getString("eval-type").compare("index_de") == 0) {

		std::cout << "Index_DE" << std::endl;

		std::vector<int> indexes = dcdGetIndexLessThan(ga_config.getString("index-file"), ga_config.getDouble("index-max"));
		std::vector< std::vector<PDBAtom> > initial_population;

		EXAFSDE de(exafs_evaluator, ga_config.getDouble("f"), ga_config.getDouble("cr"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZsByIndex(ga_config.getString("dcd-file"), indexes);
		std::cout << "DCD Population = " << initial_dcd_population.size() << std::endl;
		for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

			pdb_helper->updateEXAFSAtomsFromXYZ(*i);
			initial_population.push_back( pdb_helper->getEXAFSAtoms() );
		}

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			
			// Get subset of population.
			std::random_shuffle(initial_population.begin(), initial_population.end());
			std::vector< std::vector<PDBAtom> > subset_pop = std::vector< std::vector<PDBAtom> >(initial_population.begin(), initial_population.begin()+ga_config.getInt("population-size"));

			initial_populations.push_back(subset_pop);
		}

		DCDHelper dcd_helper = DCDHelper(ga_config.getString("dcd-file"));
		pdb_helper->updateAllAtomsFromXYZ(dcd_helper.getXYZAtFrame(0));
		
		std::cout << "DE: Begin" << std::endl;
		de.begin(initial_populations);

	} else if (ga_config.getString("eval-type").compare("index_de_recenter") == 0) {

		std::cout << "Index_DE_Recenter" << std::endl;

		std::vector<int> indexes = dcdGetIndexLessThan(ga_config.getString("index-file"), ga_config.getDouble("index-max"));
		std::vector< std::vector<PDBAtom> > initial_population;

		EXAFSDE de(exafs_evaluator, ga_config.getDouble("f"), ga_config.getDouble("cr"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZsByIndex(ga_config.getString("dcd-file"), indexes);
		std::cout << "DCD Population = " << initial_dcd_population.size() << std::endl;
		for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

			pdb_helper->updateEXAFSAtomsFromXYZ(*i);
			initial_population.push_back( pdb_helper->getEXAFSAtoms() );
		}

		std::cout << "DE: Begin" << std::endl;

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			
			// Get subset of population.
			std::random_shuffle(initial_population.begin(), initial_population.end());
			std::vector< std::vector<PDBAtom> > subset_pop = std::vector< std::vector<PDBAtom> >(initial_population.begin(), initial_population.begin()+ga_config.getInt("recentering-population"));

			initial_populations.push_back(subset_pop);
		}
		
		de.begin_recentering(initial_populations, ga_config.getInt("population-size"), ga_config.getDouble("keep-percentage"), ga_config.getInt("recentering"));

	} else if (ga_config.getString("eval-type").compare("index_pso") == 0) {

		std::cout << "Index_PSO" << std::endl;

		std::vector<int> indexes = dcdGetIndexLessThan(ga_config.getString("index-file"), ga_config.getDouble("index-max"));
		std::vector< std::vector<PDBAtom> > initial_population;

		EXAFSPSO pso(exafs_evaluator, ga_config.getDouble("inertia"), ga_config.getDouble("social"), ga_config.getDouble("cognitive"), ga_config.getDouble("velocity-range"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZsByIndex(ga_config.getString("dcd-file"), indexes);
		std::cout << "DCD Population = " << initial_dcd_population.size() << std::endl;
		for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

			pdb_helper->updateEXAFSAtomsFromXYZ(*i);
			initial_population.push_back( pdb_helper->getEXAFSAtoms() );
		}

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			
			// Get subset of population.
			std::random_shuffle(initial_population.begin(), initial_population.end());
			std::vector< std::vector<PDBAtom> > subset_pop = std::vector< std::vector<PDBAtom> >(initial_population.begin(), initial_population.begin()+ga_config.getInt("population-size"));

			initial_populations.push_back(subset_pop);
		}

		DCDHelper dcd_helper = DCDHelper(ga_config.getString("dcd-file"));
		pdb_helper->updateAllAtomsFromXYZ(dcd_helper.getXYZAtFrame(0));

		std::cout << "PSO: Begin" << std::endl;
		pso.begin(initial_populations);

	} else if (ga_config.getString("eval-type").compare("full-dcd-potential") == 0) {
		std::cout << "Full DCD Potential" << std::endl;

		DCDHelper dcd_helper = DCDHelper(ga_config.getString("dcd-file"));

		// Get all the xyz entries.
		// std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZs(ga_config.getString("dcd-file"));

		std::ofstream dcd_results((std::string("dcd_potential.csv")).c_str());
		// for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {
		for (int i = 0; i < dcd_helper.numberOfFrames(); ++i) {
			
			pdb_helper->updateAllAtomsFromXYZ(dcd_helper.getXYZAtFrame(i));
			pdb_helper->writePDBFile();

			double potential_energy = vmd_helper->calculateEnergy();
			std::cout << "Energy: " << i << " - " << potential_energy << std::endl;

			dcd_results << potential_energy << std::endl;
		}
		// 	pdb_helper->updateAllAtomsFromXYZ(*i);
			
		// 	pdb_helper->writePDBFile();
		// 	double potential_energy = vmd_helper->calculateEnergy();
		// 	std::cout << "Energy: " << potential_energy << std::endl;

		// 	dcd_results << potential_energy << std::endl;
		// }
		dcd_results.close();
	} else {
		std::cout << "Other" << std::endl;

		// Get all the xyz entries.
		std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZs(ga_config.getString("dcd-file"));

		std::ofstream dcd_results((std::string("dcd_results.csv")).c_str());
		for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {
			pdb_helper->updateEXAFSAtomsFromXYZ(*i);
			double initial_rmsd = ifeffit_helper->run(pdb_helper->getEXAFSAtoms(), true);
			
			std::cout << initial_rmsd << std::endl;
			dcd_results << initial_rmsd << std::endl;
		}
		dcd_results.close();
	}

	// system(("rm -rf " + temp_folder).c_str());
*/

	return 0;
}
