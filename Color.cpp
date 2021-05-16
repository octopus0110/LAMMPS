#include "/mnt/c/Ubuntu/home/kaito/lammps/header/LAMMPS.h"
#include <iostream>

#define DEBUG false

int Num_atoms_in_monomor = 18; // モノマーの中の原子数
map<int, vector<int>> atoms_of_molecule_ID; // 各分子にどの原子が含まれているか

void Color(Data& Data)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  for (int i = 0; i < Data.atoms; ++i)
  {
    if (!Data.Atoms.at(i+1).have_bond) Data.Atoms.at(i+1).molecule_ID = 0;
    else
    {
      int atom_ID, molecule_ID;
      atom_ID = Data.Atoms.at(i+1).atom_ID;
      molecule_ID = Data.Atoms.at(i+1).molecule_ID;
      atoms_of_molecule_ID[molecule_ID].emplace_back(atom_ID);
    }
  }
  for (pair<int, vector<int>> x: atoms_of_molecule_ID)
  {
    int molecule_ID, Num_atoms, PD;
    molecule_ID = x.first;
    Num_atoms = x.second.size();
    PD = Num_atoms / Num_atoms_in_monomor;
    for (int ID: x.second) Data.Atoms.at(ID).molecule_ID = PD;
  }
}

int main(int argc, char** argv)
{
  Data Data;
  if (argc == 1) Data.Read();
  if (argc >= 2) Data.Read(argv[1]);
  if (argc >= 3) Num_atoms_in_monomor = stoi(argv[2]);
  Color(Data);
  Data.Remove_Molecule_ID({0});
  Data.Write();

  return 0;
}