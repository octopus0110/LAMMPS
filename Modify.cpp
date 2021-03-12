#include "/mnt/c/Ubuntu/home/usr/lammps/header/LAMMPS.h"
#include <iostream>
#include <queue>

#define DEBUG false

map<int, vector<int>> atoms_of_molecule_ID; // 各分子にどの原子が含まれているか
double CR;  // Conversion Ratio

// Monomor
double MW_monomor;  // モノマーの分子量
map<int, int> Num_atom_types_in_monomor;  // モノマーの中にどの原子が何個含まれているか
int Num_atoms_in_monomor; // モノマーの中の原子数
int Num_monomor;  // 初めのモノマーの数

// Polymer
int Num_polymer;  // 分子数
map<int, int> Num_atoms_of_molecule_ID; // 各分子に含まれている原子数
map<int, int> PD_of_molecule_ID;  // 各分子の重合度
map<int, int> Num_PD; // 各重合度の分子が何個あるか
int Num_ring_polymer; // 環状分子の数

int lower(const int& a, const int& b)
{
  if (a < b) return a;
  else return b;
}

int higher(const int& a, const int& b)
{
  if (a > b) return a;
  else return b;
}

void get_atoms_of_molecule_ID(const Data& Data)
{
  for (int i = 0; i < Data.atoms; ++i)
  {
    atoms_of_molecule_ID[Data.Atoms.at(i+1).molecule_ID].emplace_back(Data.Atoms.at(i+1).atom_ID);
  }
}

int get_molecule_ID_monomor(const Data& Data)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int molecule_ID_monomor = 1;
  for (int i = 0; i < Data.atoms; ++i)
  {
    if (Data.Atoms.at(i+1).molecule_ID == molecule_ID_monomor)
    {
      if (!Data.Atoms.at(i+1).have_bond)
      {
        i = -1;
        ++molecule_ID_monomor;
      }
      else break;
    }
  }

  return molecule_ID_monomor;
}

void get_monomor(const Data& Data)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  MW_monomor = 0.0;
  Num_atoms_in_monomor = 0;
  Num_monomor = 0;
  int molecule_ID_monomor;
  get_atoms_of_molecule_ID(Data);
  molecule_ID_monomor = get_molecule_ID_monomor(Data);
  for (int i = 0; i < Data.atoms; ++i)
  {
    if (Data.Atoms.at(i+1).molecule_ID == molecule_ID_monomor)
    {
      MW_monomor += Data.Atoms.at(i+1).mass;
      ++Num_atom_types_in_monomor[Data.Atoms.at(i+1).atom_type];
      ++Num_atoms_in_monomor;
    }
  }
  for (int i = 0; i < Data.atoms; ++i)
  {
    if (Data.Atoms.at(i+1).have_bond) ++Num_monomor;
  }
  Num_monomor /= Num_atoms_in_monomor;
}

void get_polymer(const Data& Data)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  for (int i = 0; i < Data.atoms; ++i)
  {
    if (!Data.Atoms.at(i+1).have_bond) continue;
    int ID = Data.Atoms.at(i+1).molecule_ID;
    ++Num_atoms_of_molecule_ID[ID];
  }
  for (int i = 0; i < Num_polymer; ++i)
  {
    int PD;
    if (Num_atoms_of_molecule_ID[i+1] == 0) continue;
    PD = Num_atoms_of_molecule_ID[i+1] / Num_atoms_in_monomor;
    PD_of_molecule_ID[i+1] = PD;
    ++Num_PD[PD];
  }
}

void Write_Data(string filename = "")
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  priority_queue< pair<int, int>, std::vector<pair<int, int>>, std::greater<pair<int, int>> > queue;
  if (filename == "") filename = "Data.data";
  ofstream ofs(filename);
  CR = double(Num_monomor-Num_PD[1]) / double(Num_monomor);
  ofs << "Molecule Weight of Monomor : " << MW_monomor << endl;
  ofs << "Molecule Weight of All     : " << MW_monomor*double(Num_monomor) << endl;
  ofs << "Initial Number of Monomor  : " << Num_monomor << endl;
  ofs << "Number of Polymers         : " << Num_polymer << endl;
  ofs << "Conversion Ratio           : " << CR*100.0 << endl;
  // ofs << "Number of Ring : " << Num_ring_polymer << endl;
  ofs << endl << "<Distribution>" << endl;
  for (pair<int, int> x: Num_PD)
  {
    int PD, Num;
    PD = x.first;
    Num = x.second;
    if (Num > 0) queue.push({PD, Num});
  }
  while(!queue.empty())
  {
    pair<int, int> x = queue.top();
    queue.pop();
    int PD, Num;
    PD = x.first;
    Num = x.second;
    ofs << PD << " " << Num << " " << double(PD)*MW_monomor << endl;
  }
}

void Modify(Data& Data)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  Num_polymer = 0;
  Num_ring_polymer = 0;
  map<int, int> convert_ID;
  map<int, bool> check;

  for (int i = 0; i < Data.atoms; ++i)
  {
    if (!Data.Atoms.at(i+1).have_bond) Data.Atoms.at(i+1).molecule_ID = 0;
  }

  for (int i = 0; i < Data.bonds; ++i)
  {
    if (Data.Bonds.at(i+1).atom[0]->molecule_ID != Data.Bonds.at(i+1).atom[1]->molecule_ID)
    {
      int ID_low, ID_high;
      if (Data.Bonds.at(i+1).atom[0]->molecule_ID < Data.Bonds.at(i+1).atom[1]->molecule_ID)
      {
        ID_low = Data.Bonds.at(i+1).atom[0]->molecule_ID;
        ID_high = Data.Bonds.at(i+1).atom[1]->molecule_ID;
      }
      else
      {
        ID_low = Data.Bonds.at(i+1).atom[1]->molecule_ID;
        ID_high = Data.Bonds.at(i+1).atom[0]->molecule_ID;
      }

      for (int atom_ID: atoms_of_molecule_ID[ID_high])
      {
        Data.Atoms.at(atom_ID).molecule_ID = ID_low;
        atoms_of_molecule_ID[ID_low].emplace_back(atom_ID);
      }
    }
  }

  for (int i = 0; i < Data.atoms; ++i) check[Data.Atoms.at(i+1).molecule_ID] = true;

  int ID_new = 0;
  for (pair<int, int> x: check)
  {
    ++ID_new;
    int ID_old = x.first;
    convert_ID[ID_old] = ID_new;
  }
  Num_polymer = check.size() - 1;
  for (int i = 0; i < Data.atoms; ++i) Data.Atoms.at(i+1).molecule_ID = convert_ID.at(Data.Atoms.at(i+1).molecule_ID);
}

void Run(Data& Data)
{
  get_monomor(Data);
  Modify(Data);
  get_polymer(Data);
  Data.Write();
  Write_Data();
}

int main(int argc, char** argv)
{
  Data Data;
  if (argc == 2) Data.Read(argv[1]);
  else Data.Read();
  Run(Data);

  return 0;
}