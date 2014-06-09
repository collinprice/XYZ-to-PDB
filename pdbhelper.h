#pragma once
#include "pdbatom.h"

#include <string>
#include <vector>

/*
	PDBHelper

	Used to read/write a PDB file.

	output_pdb_file - Default output file name.
	target_atoms - Used in EXAFS. This list is used to filter only allowed atom symbols into the EXAFS atoms.
*/

class PDBHelper {

public:

	std::string output_pdb_file;
	std::vector<std::string> target_atoms;

	PDBHelper(std::string pdb_file, std::string amber_topology_file, std::string output_pdb_file, std::vector<std::string> target_atoms);

	int numberOfAtoms();

	/*
		Returns all the atoms that are involved in EXAFS. EXAFS atoms are distinguished by have an occupancy of 1. Non-EXAFS atoms have an occupancy of 0.
	*/
	std::vector<PDBAtom> getEXAFSAtoms();

	std::vector<PDBAtom> getAllEXAFSAtoms();

	void updateEXAFSAtoms(std::vector<PDBAtom> atoms);

	// void updateEXAFSAtomsFromXYZ(std::string filename);

	/*
		Update atoms based on an list of atoms. XYZ files only contain atomic symbol and 3d coordinates. Indexes cannot be trusted.
	*/
	void updateEXAFSAtomsFromXYZ(std::vector<PDBAtom> atoms);

	void updateAllAtomsFromXYZ(std::vector<PDBAtom> atoms);
	void updateAllNonEXAFSAtomsFromXYZ(std::vector<PDBAtom> atoms);

	// void updateAtomsFromList(std::vector<PDBAtom> atoms);
	void writePDBFile(std::string filename);
	void writePDBFile();
	
private:

	/*
		This function checks the occupancy and atom symbol to determine if the atom is an EXAFS atom.
	*/
	bool validEXAFSAtom(PDBAtom atom);
	bool validEXAFSAtomicSymbol(std::string atom);
	bool validEXAFSIndex(int index);
	std::string atomAtIndex(int index);

};