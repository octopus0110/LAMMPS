#include "/mnt/c/Ubuntu/home/usr/lammps/header/LAMMPS.h"
#include <iostream>

#define DEBUG false

map<pair<string, string>, string> Pair_Coeffs;

void update_atom_style(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (D1.atom_style == D2.atom_style) New.atom_style = D1.atom_style;
  else if (D1.atom_style == "full" || D2.atom_style == "full")
  {
    if (D1.atom_style == "atomic" || D2.atom_style == "atomic") New.atom_style = "full";
    else
    {
      cerr << __FUNCTION__ << " : 未実装[" << D1.atom_style << ", " << D2.atom_style << "]" << endl;
      exit(0);
    }
  }
  else
  {
    cerr <<  __FUNCTION__ << " : 未実装[" << D1.atom_style << ", " << D2.atom_style << "]" << endl;
    exit(0);
  }
}

double higher(const double& a, const double& b)
{
  if (a > b) return a;
  else return b;
}

double lower(const double& a, const double& b)
{
  if (a < b) return a;
  else return b;
}

void update_simulation_size(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  New.xlo = lower(D1.xlo, D2.xlo);
  New.xhi = higher(D1.xhi, D2.xhi);
  New.ylo = lower(D1.ylo, D2.ylo);
  New.yhi = higher(D1.yhi, D2.yhi);
  New.zlo = lower(D1.zlo, D2.zlo);
  New.zhi = higher(D1.zhi, D2.zhi);
}

void update_header(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  New.dimension = D1.dimension;
  New.line_header = "LAMMPS Combined";
  update_atom_style(New, D1, D2);
  New.atoms = D1.atoms + D2.atoms;
  New.atom_types = D1.atom_types + D2.atom_types;
  if (New.atom_style == "full")
  {
    New.bonds = D1.bonds + D2.bonds;
    New.bond_types = D1.bond_types + D2.bond_types;
    New.angles = D1.angles + D2.angles;
    New.angle_types = D1.angle_types + D2.angle_types;
    New.dihedrals = D1.dihedrals + D2.dihedrals;
    New.dihedral_types = D1.dihedral_types + D2.dihedral_types;
    New.impropers = D1.impropers + D2.impropers;
    New.improper_types = D1.improper_types + D2.improper_types;
  }
  update_simulation_size(New, D1, D2);
}

void update_Masses(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int type = 0;
  for (int i = 0; i < D1.atom_types; ++i)
  {
    ++type;
    New.Masses[type] = D1.Masses.at(i+1);
    New.Masses[type].atom_type = type;
  }
  for (int i = 0; i < D2.atom_types; ++i)
  {
    ++type;
    New.Masses[type] = D2.Masses.at(i+1);
    New.Masses[type].atom_type = type;
  }
}

void update_pair_coeff_style(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (New.atom_style != "full" && New.atom_style != "atomic")
  {
    cerr << __FUNCTION__ << " : 未実装[" << New.atom_style << "]" << endl;
    exit(0);
  }
  else if (New.atom_style == "atomic") return;
  if (D1.atom_style == "full" && D2.atom_style == "atomic") New.pair_coeff_style = D1.pair_coeff_style;
  else if (D1.atom_style == "atomic" && D2.atom_style == "full") New.pair_coeff_style = D2.pair_coeff_style;
  else if (D1.pair_coeff_style != D2.pair_coeff_style)
  {
    cerr << __FUNCTION__ << " : D1.pair_coeff_style != D2.pair_coeff_style" << endl;
    exit(0);
  }
  else New.pair_coeff_style = D1.pair_coeff_style;
}

void update_Pair_Coeffs(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  update_pair_coeff_style(New, D1, D2);
  if (New.pair_coeff_style == "") return;
  int id = 0;
  if (D1.atom_style == "atomic")
  {
    for (int i = 0; i < D1.atom_types; ++i)
    {
      ++id;
      New.Pair_Coeffs[id].atom_type = id;
      New.Pair_Coeffs[id].coeffs = Pair_Coeffs.at({"Si", New.pair_coeff_style});
    }
  }
  else
  {
    for (int i = 0; i < D1.atom_types; ++i)
    {
      ++id;
      New.Pair_Coeffs[id] = D1.Pair_Coeffs.at(i+1);
      New.Pair_Coeffs[id].atom_type = id;
    }
  }
  if (D2.atom_style == "atomic")
  {
    for (int i = 0; i < D2.atom_types; ++i)
    {
      ++id;
      New.Pair_Coeffs[id].atom_type = id;
      New.Pair_Coeffs[id].coeffs = Pair_Coeffs.at({"Si", New.pair_coeff_style});
    }
  }
  else
  {
    for (int i = 0; i < D2.atom_types; ++i)
    {
      ++id;
      New.Pair_Coeffs[id] = D2.Pair_Coeffs.at(i+1);
      New.Pair_Coeffs[id].atom_type = id;
    }
  }
}

void update_bond_coeff_style(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (New.atom_style != "full" && New.atom_style != "atomic")
  {
    cerr << __FUNCTION__ << " : 未実装[" << New.atom_style << "]" << endl;
    exit(0);
  }
  else if (New.atom_style == "atomic") return;
  if (D1.atom_style == "full" && D2.atom_style == "atomic") New.bond_coeff_style = D1.bond_coeff_style;
  else if (D1.atom_style == "atomic" && D2.atom_style == "full") New.bond_coeff_style = D2.bond_coeff_style;
  else if (D1.bond_coeff_style != D2.bond_coeff_style)
  {
    cerr << __FUNCTION__ << " : D1.bond_coeff_style != D2.bond_coeff_style" << endl;
    exit(0);
  }
  else New.bond_coeff_style = D1.bond_coeff_style;
}

void update_Bond_Coeffs(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  update_bond_coeff_style(New, D1, D2);
  if (New.bond_coeff_style == "") return;
  int id = 0;
  if (D1.bond_coeff_style != "")
  {
    for (int i = 0; i < D1.bond_types; ++i)
    {
      ++id;
      New.Bond_Coeffs[id] = D1.Bond_Coeffs.at(i+1);
      New.Bond_Coeffs[id].bond_type = id;
    }
  }
  if (D2.bond_coeff_style != "")
  {
    for (int i = 0; i < D2.bond_types; ++i)
    {
      ++id;
      New.Bond_Coeffs[id] = D2.Bond_Coeffs.at(i+1);
      New.Bond_Coeffs[id].bond_type = id;
    }
  }
}

void update_angle_coeff_style(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (New.atom_style != "full" && New.atom_style != "atomic")
  {
    cerr << __FUNCTION__ << " : 未実装[" << New.atom_style << "]" << endl;
    exit(0);
  }
  else if (New.atom_style == "atomic") return;
  if (D1.atom_style == "full" && D2.atom_style == "atomic") New.angle_coeff_style = D1.angle_coeff_style;
  else if (D1.atom_style == "atomic" && D2.atom_style == "full") New.angle_coeff_style = D2.angle_coeff_style;
  else if (D1.angle_coeff_style != D2.angle_coeff_style)
  {
    cerr << __FUNCTION__ << " : D1.angle_coeff_style != D2.angle_coeff_style" << endl;
    exit(0);
  }
  else New.angle_coeff_style = D1.angle_coeff_style;
}

void update_Angle_Coeffs(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  update_angle_coeff_style(New, D1, D2);
  if (New.angle_coeff_style == "") return;
  int id = 0;
  if (D1.angle_coeff_style != "")
  {
    for (int i = 0; i < D1.angle_types; ++i)
    {
      ++id;
      New.Angle_Coeffs[id] = D1.Angle_Coeffs.at(i+1);
      New.Angle_Coeffs[id].angle_type = id;
    }
  }
  if (D2.angle_coeff_style != "")
  {
    for (int i = 0; i < D2.angle_types; ++i)
    {
      ++id;
      New.Angle_Coeffs[id] = D2.Angle_Coeffs.at(i+1);
      New.Angle_Coeffs[id].angle_type = id;
    }
  }
}

void update_dihedral_coeff_style(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (New.atom_style != "full" && New.atom_style != "atomic")
  {
    cerr << __FUNCTION__ << " : 未実装[" << New.atom_style << "]" << endl;
    exit(0);
  }
  else if (New.atom_style == "atomic") return;
  if (D1.atom_style == "full" && D2.atom_style == "atomic") New.dihedral_coeff_style = D1.dihedral_coeff_style;
  else if (D1.atom_style == "atomic" && D2.atom_style == "full") New.dihedral_coeff_style = D2.dihedral_coeff_style;
  else if (D1.dihedral_coeff_style != D2.dihedral_coeff_style)
  {
    cerr << __FUNCTION__ << " : D1.dihedral_coeff_style != D2.dihedral_coeff_style" << endl;
    exit(0);
  }
  else New.dihedral_coeff_style = D1.dihedral_coeff_style;
}

void update_Dihedral_Coeffs(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  update_dihedral_coeff_style(New, D1, D2);
  if (New.dihedral_coeff_style == "") return;
  int id = 0;
  if (D1.dihedral_coeff_style != "")
  {
    for (int i = 0; i < D1.dihedral_types; ++i)
    {
      ++id;
      New.Dihedral_Coeffs[id] = D1.Dihedral_Coeffs.at(i+1);
      New.Dihedral_Coeffs[id].dihedral_type = id;
    }
  }
  if (D2.dihedral_coeff_style != "")
  {
    for (int i = 0; i < D2.dihedral_types; ++i)
    {
      ++id;
      New.Dihedral_Coeffs[id] = D2.Dihedral_Coeffs.at(i+1);
      New.Dihedral_Coeffs[id].dihedral_type = id;
    }
  }
}

void update_improper_coeff_style(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (New.atom_style != "full" && New.atom_style != "atomic")
  {
    cerr << __FUNCTION__ << " : 未実装[" << New.atom_style << "]" << endl;
    exit(0);
  }
  else if (New.atom_style == "atomic") return;
  if (D1.atom_style == "full" && D2.atom_style == "atomic") New.improper_coeff_style = D1.improper_coeff_style;
  else if (D1.atom_style == "atomic" && D2.atom_style == "full") New.improper_coeff_style = D2.improper_coeff_style;
  else if (D1.improper_coeff_style != D2.improper_coeff_style)
  {
    cerr << __FUNCTION__ << " : D1.improper_coeff_style != D2.improper_coeff_style" << endl;
    exit(0);
  }
  else New.improper_coeff_style = D1.improper_coeff_style;
}

void update_Improper_Coeffs(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  update_improper_coeff_style(New, D1, D2);
  if (New.improper_coeff_style == "") return;
  int id = 0;
  if (D1.improper_coeff_style != "")
  {
    for (int i = 0; i < D1.improper_types; ++i)
    {
      ++id;
      New.Improper_Coeffs[id] = D1.Improper_Coeffs.at(i+1);
      New.Improper_Coeffs[id].improper_type = id;
    }
  }
  if (D2.improper_coeff_style != "")
  {
    for (int i = 0; i < D2.improper_types; ++i)
    {
      ++id;
      New.Improper_Coeffs[id] = D2.Improper_Coeffs.at(i+1);
      New.Improper_Coeffs[id].improper_type = id;
    }
  }
}

void update_Atoms(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int id = 0;
  for (int i = 0; i < D1.atoms; ++i)
  {
    ++id;
    New.Atoms[id] = D1.Atoms.at(i+1);
    New.Atoms[id].atom_ID = id;
  }
  for (int i = 0; i < D2.atoms; ++i)
  {
    ++id;
    New.Atoms[id] = D2.Atoms.at(i+1);
    New.Atoms[id].atom_ID = id;
    New.Atoms[id].atom_type += D1.atom_types;
  }
}

void update_Bonds(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int atom_num = 2;
  int id = 0;
  for (int i = 0; i < D1.bonds; ++i)
  {
    ++id;
    New.Bonds[id] = D1.Bonds.at(i+1);
    New.Bonds[id].bond_ID = id;
    for (int j = 0; j < atom_num; ++j)
    {
      New.Bonds.at(id).atom[j] = &New.Atoms.at(New.Bonds.at(id).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < D2.bonds; ++i)
  {
    ++id;
    New.Bonds[id] = D2.Bonds.at(i+1);
    New.Bonds[id].bond_ID = id;
    for (int j = 0; j < atom_num; ++j)
    {
      New.Bonds.at(id).atom[j] = &New.Atoms.at(New.Bonds.at(id).atom[j]->atom_ID + D1.atoms);
    }
  }
}

void update_Angles(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int atom_num = 3;
  int id = 0;
  for (int i = 0; i < D1.angles; ++i)
  {
    ++id;
    New.Angles[id] = D1.Angles.at(i+1);
    New.Angles[id].angle_ID = id;
    for (int j = 0; j < atom_num; ++j)
    {
      New.Angles.at(id).atom[j] = &New.Atoms.at(New.Angles.at(id).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < D2.angles; ++i)
  {
    ++id;
    New.Angles[id] = D2.Angles.at(i+1);
    New.Angles[id].angle_ID = id;
    for (int j = 0; j < atom_num; ++j)
    {
      New.Angles.at(id).atom[j] = &New.Atoms.at(New.Angles.at(id).atom[j]->atom_ID + D1.atoms);
    }
  }
}

void update_Dihedrals(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int atom_num = 4;
  int id = 0;
  for (int i = 0; i < D1.dihedrals; ++i)
  {
    ++id;
    New.Dihedrals[id] = D1.Dihedrals.at(i+1);
    New.Dihedrals[id].dihedral_ID = id;
    for (int j = 0; j < atom_num; ++j)
    {
      New.Dihedrals.at(id).atom[j] = &New.Atoms.at(New.Dihedrals.at(id).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < D2.dihedrals; ++i)
  {
    ++id;
    New.Dihedrals[id] = D2.Dihedrals.at(i+1);
    New.Dihedrals[id].dihedral_ID = id;
    for (int j = 0; j < atom_num; ++j)
    {
      New.Dihedrals.at(id).atom[j] = &New.Atoms.at(New.Dihedrals.at(id).atom[j]->atom_ID + D1.atoms);
    }
  }
}

void update_Impropers(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int atom_num = 4;
  int id = 0;
  for (int i = 0; i < D1.impropers; ++i)
  {
    ++id;
    New.Impropers[id] = D1.Impropers.at(i+1);
    New.Impropers[id].improper_ID = id;
    for (int j = 0; j < atom_num; ++j)
    {
      New.Impropers.at(id).atom[j] = &New.Atoms.at(New.Impropers.at(id).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < D2.impropers; ++i)
  {
    ++id;
    New.Impropers[id] = D2.Impropers.at(i+1);
    New.Impropers[id].improper_ID = id;
    for (int j = 0; j < atom_num; ++j)
    {
      New.Impropers.at(id).atom[j] = &New.Atoms.at(New.Impropers.at(id).atom[j]->atom_ID + D1.atoms);
    }
  }
}

void update_body(Data& New, const Data& D1, const Data& D2)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  update_Masses(New, D1, D2);
  update_Pair_Coeffs(New, D1, D2);
  update_Bond_Coeffs(New, D1, D2);
  update_Angle_Coeffs(New, D1, D2);
  update_Dihedral_Coeffs(New, D1, D2);
  update_Improper_Coeffs(New, D1, D2);
  update_Atoms(New, D1, D2);
  update_Bonds(New, D1, D2);
  update_Angles(New, D1, D2);
  update_Dihedrals(New, D1, D2);
  update_Impropers(New, D1, D2);
}

Data Combine(const Data& D1, const Data& D2)
{
  Data New;
  update_header(New, D1, D2);
  update_body(New, D1, D2);
  return New;
}

int main(int argc, char** argv)
{
Pair_Coeffs[{"Si", "lj/cut/coul/long"}] = "0.31 3.80414";

  Data New, D1, D2;
  if (argc == 3)
  {
    D1.Read(argv[1]);
    D2.Read(argv[2]);
  }
  else
  {
    D1.Read();
    D2.Read();
  }
  New = Combine(D1, D2);
  New.Write();
  return 0;
}