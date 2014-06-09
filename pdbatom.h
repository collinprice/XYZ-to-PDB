#pragma once
#include <string>
#include <string.h>

/*
	PDBAtom

	Represents an atom within a pdb file.

	index - The unique atom index within a pdb file. Only way to tell atoms apart besides atomic symbol.
*/

class PDBAtom {

public:
	
	std::string atomic_symbol;
	double x;
	double y;
	double z;

	PDBAtom(double x, double y, double z);
	PDBAtom(int index, double x, double y, double z);
	PDBAtom(std::string atomic_symbol, double x, double y, double z);
	PDBAtom(std::string atomic_symbol, int index, double x, double y, double z);
	int getAtomicNumber();
	int getIndex();

	/*
		Returns the distance between two atoms based on their Euclidean distance.
	*/
	double distance(PDBAtom atom);

	static std::string atomicNumberToSymbol(int atomic_number);
	static int atomicSymbolToNumber(std::string atomic_symbol);

	PDBAtom operator+(const PDBAtom& c) {
		return PDBAtom(this->atomic_symbol, this->index, this->x + c.x, this->y + c.y, this->z + c.z);
	};

	PDBAtom operator-(const PDBAtom& c) {
		return PDBAtom(this->atomic_symbol, this->index, this->x - c.x, this->y - c.y, this->z - c.z);
	};

	PDBAtom operator*(const double& c) {
		return PDBAtom(this->atomic_symbol, this->index, this->x * c, this->y * c, this->z * c);
	};

private:

	int index;
	const static int periodic_table_size = 103;
	const static char* periodic_table[];

};