#include "/mnt/c/Ubuntu/home/kaito/lammps/header/LAMMPS.h"
#include <iostream>
#include <unordered_map>

#define DEBUG false

#define oxygennumber_of_monomor 2.0

void develop(Data& Data, const int& threshold, const int& atom_type_oxygen)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int molecules = 0;
  unordered_map<int, int> oxygennumber;
  unordered_map<int, bool> count;
  vector<int> remove;
  for (int i = 0; i < Data.atoms; ++i)
  {
    if (Data.Atoms.at(i+1).molecule_ID <= 0) continue;
    count[Data.Atoms.at(i+1).molecule_ID] = true;
    if ((int)Data.Atoms.at(i+1).mass == (int)Data.Masses.at(atom_type_oxygen).mass)
    {
      ++oxygennumber[Data.Atoms.at(i+1).molecule_ID];
    }
  }
  molecules = count.size();
  for (int i = 0; i < molecules; ++i)
  {
    oxygennumber[i+1] = oxygennumber[i+1];
    double polymerization;
    polymerization = oxygennumber.at(i+1) / oxygennumber_of_monomor;
    if (polymerization <= threshold) remove.emplace_back(i+1);
  }

  cout << "-------------------------------------------------------------------------------------------------" << endl;
  cout << "Molecules         : " << molecules << endl;
  cout << "Removed Molecules : " << remove.size() << endl;

  Data.Remove_Molecule_ID(remove);
}

int main(int argc, char** argv)
{
  Data Data;
  int threshold;
  int atom_type_oxygen;
  if (argc >= 4)
  {
    Data.Read(argv[1]);
    threshold = stoi(argv[2]);
    atom_type_oxygen = stoi(argv[3]);
  }
  else
  {
    Data.Read();
    cout << "threshold = ";
    cin >> threshold;
    cout << "The atom type of oxygen : ";
    cin >> atom_type_oxygen;
  }
  develop(Data, threshold, atom_type_oxygen);

  Data.Write("Developed.data");

  return 0;
}