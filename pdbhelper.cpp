#include "pdbhelper.h"

extern "C" {
#include "file_read_write.h"
}

#include <iostream>
#include <fstream>
#include <cstdlib>

#define C_TEXT( text ) ((char*)std::string( text ).c_str())

PDBHelper::PDBHelper(std::string pdb_file, std::string amber_topology_file, std::string output_pdb_file, std::vector<std::string> target_atoms) {

	read_pdb(C_TEXT(pdb_file));
	read_parm7(C_TEXT(amber_topology_file)); // Restore reliable atom names
	mass_to_element();

	this->output_pdb_file = output_pdb_file;
	this->target_atoms = target_atoms;
}

int PDBHelper::numberOfAtoms() {
	return N;
}

std::string PDBHelper::atomAtIndex(int index) {
	if (Element[index][1] == '\0') {
		return std::string(1, Element[index][0]);
	} else {
		return std::string(1, Element[index][0]) + std::string(1, Element[index][1]);
	}
}

std::vector<PDBAtom> PDBHelper::getEXAFSAtoms() {

	std::vector<PDBAtom> exafs_atoms;

	for (int i = 0; i < N; ++i) {
		if (this->validEXAFSIndex(i) && this->validEXAFSAtomicSymbol(this->atomAtIndex(i))) {
			exafs_atoms.push_back(PDBAtom(this->atomAtIndex(i), i, X[i], Y[i], Z[i]));
		}
	}

	return exafs_atoms;
}

std::vector<PDBAtom> PDBHelper::getAllEXAFSAtoms() {

	std::vector<PDBAtom> exafs_atoms;

	for (int i = 0; i < N; ++i) {
		if (this->validEXAFSIndex(i)) {
			exafs_atoms.push_back(PDBAtom(this->atomAtIndex(i), i, X[i], Y[i], Z[i]));
		}
	}

	return exafs_atoms;
}
	
void PDBHelper::updateEXAFSAtoms(std::vector<PDBAtom> atoms) {

	for (std::vector<PDBAtom>::iterator atom = atoms.begin(); atom != atoms.end(); ++atom) {
		
		int index = atom->getIndex();

		if (this->validEXAFSAtom(*atom)) {
			X[index] = atom->x;
			Y[index] = atom->y;
			Z[index] = atom->z;
		}
		
	}
}

void PDBHelper::updateEXAFSAtomsFromXYZ(std::vector<PDBAtom> atoms) {

	if ((int)atoms.size() != N) {
		throw "PDB File does not contain the same amount of atoms as input XYZ file.";
	}

	for (int i = 0; i < N; ++i) {
		
		if (this->validEXAFSIndex(i) && this->validEXAFSAtomicSymbol(this->atomAtIndex(i))) {
			X[i] = atoms[i].x;
			Y[i] = atoms[i].y;
			Z[i] = atoms[i].z;
		}
	}
}

void PDBHelper::updateAllAtomsFromXYZ(std::vector<PDBAtom> atoms) {

	if ((int)atoms.size() != N) {
		throw "PDB File does not contain the same amount of atoms as input XYZ file.";
	}

	for (int i = 0; i < N; ++i) {
		
		X[i] = atoms[i].x;
		Y[i] = atoms[i].y;
		Z[i] = atoms[i].z;
	}
}

void PDBHelper::updateAllNonEXAFSAtomsFromXYZ(std::vector<PDBAtom> atoms) {

	if ((int)atoms.size() != N) {
		throw "PDB File does not contain the same amount of atoms as input XYZ file.";
	}

	for (int i = 0; i < N; ++i) {
		
		if (!this->validEXAFSIndex(i)) {
			X[i] = atoms[i].x;
			Y[i] = atoms[i].y;
			Z[i] = atoms[i].z;
		}
	}
}

void PDBHelper::writePDBFile(std::string filename) {
	write_one_pdb(C_TEXT(filename));
}

void PDBHelper::writePDBFile() {
	this->writePDBFile(this->output_pdb_file);
}

bool PDBHelper::validEXAFSAtom(PDBAtom atom) {

	if (occupancy[atom.getIndex()] == 1) {
		for (std::vector<std::string>::iterator i = this->target_atoms.begin(); i != this->target_atoms.end(); ++i) {
			if ((*i).compare(atom.atomic_symbol) == 0) return true;
		}
	}

	return false;
}

bool PDBHelper::validEXAFSAtomicSymbol(std::string atom) {
	
	for (std::vector<std::string>::iterator i = this->target_atoms.begin(); i != this->target_atoms.end(); ++i) {
		if (*i == atom) return true;
	}
	return false;
}

bool PDBHelper::validEXAFSIndex(int index) {
	return occupancy[index] == 1;
}