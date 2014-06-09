#include "pdbatom.h"
#include "math.h"

PDBAtom::PDBAtom(double x, double y, double z) {

	this->index = -1;
	this->x = x;
	this->y = y;
	this->z = z;
}

PDBAtom::PDBAtom(int index, double x, double y, double z) {

	this->index = index;
	this->x = x;
	this->y = y;
	this->z = z;
}

PDBAtom::PDBAtom(std::string atomic_symbol, double x, double y, double z) {
	
	this->atomic_symbol = atomic_symbol;
	this->index = -1;
	this->x = x;
	this->y = y;
	this->z = z;	
}

PDBAtom::PDBAtom(std::string atomic_symbol, int index, double x, double y, double z) {

	this->atomic_symbol = atomic_symbol;
	this->index = index;
	this->x = x;
	this->y = y;
	this->z = z;
}

int PDBAtom::getAtomicNumber() {
	for (int i = 0; i < this->periodic_table_size; ++i) {
		if (strcmp(this->periodic_table[i], this->atomic_symbol.c_str()) == 0) {
			return (i+1);
		}
	}

	return -1;
}

int PDBAtom::getIndex() {
	return this->index;
}

double PDBAtom::distance(PDBAtom atom) {

	return sqrt(pow(atom.x - this->x, 2) + pow(atom.y - this->y, 2) + pow(atom.z - this->z, 2));
}

std::string PDBAtom::atomicNumberToSymbol(int atomic_number) {
	return PDBAtom::periodic_table[atomic_number-1];
}

int PDBAtom::atomicSymbolToNumber(std::string atomic_symbol) {
	for (int i = 0; i < PDBAtom::periodic_table_size; ++i) {
		if (strcmp(PDBAtom::periodic_table[i], atomic_symbol.c_str()) == 0) {
			return (i+1);
		}
	}

	return -1;
}

const char* PDBAtom::periodic_table[] = {
	"H", "He","Li","Be","B", "C", "N", "O", "F", "Ne",
	"Na","Mg","Al","Si","P", "S", "Cl","Ar","K", "Ca",
	"Sc","Ti","V", "Cr","Mn","Fe","Co","Ni","Cu","Zn",
	"Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y", "Zr",
	"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
	"Sb","Te","I", "Xe","Cs","Ba","La","Ce","Pr","Nd",
	"Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
	"Lu","Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg",
	"Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
	"Pa","U", "Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
	"Md","No","Lr"
};