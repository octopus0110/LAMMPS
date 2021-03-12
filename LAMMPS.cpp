#include "/mnt/c/Ubuntu/home/usr/lammps/header/LAMMPS.h"
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

#define setpre(n) fixed << setprecision((n))
#define DIGITS 16
#define DEBUG false

Atom::Atom()
{
  // Atom
  atom_ID = 0; molecule_ID = 0; atom_type = 0; etag = 0; spin = 0; spx = 0; spy = 0; spz = 0; template_atom = 0; template_index = 0;
  ccN = 0.0; contact_radius = 0.0; cs_re = 0.0; cs_im = 0.0; cv = 0.0; density = 0.0; diameter = 0.0; espf = 0.0; edpd_temp = 0.0; edpd_cv = 0.0; eradius = 0.0; kernel_radius = 0.0; mass = 0.0; mux = 0.0; muy = 0.0; muz = 0.0; q = 0.0; rho = 0.0; sp = 0.0; theta = 0.0; volume = 0.0; x = 0.0; y = 0.0; z = 0.0; x0 = 0.0; y0 = 0.0; z0 = 0.0;
  ellipsoidflag = false; lineflag = false; triangleflag = false, nx = false, ny = false, nz = false;

  // Velocity
  vx = 0.0; vy = 0.0; vz = 0.0; lx = 0.0; ly = 0.0; lz = 0.0; wx = 0.0; wy = 0.0; wz = 0.0; ervel = 0.0;

  have_bond = false;
}

void Atom::write_Atom(const string& ATOM_STYLE, ofstream& ofs)
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  if (ATOM_STYLE == "full")
  {
    ofs << setpre(DIGITS) << atom_ID << " " << molecule_ID << " " << atom_type << " " << q << " " << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << endl;
  }
  else if (ATOM_STYLE == "atomic")
  {
    ofs << setpre(DIGITS) << atom_ID << " " << atom_type << " " << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << endl;
  }
  else
  {
    cerr << "Atom::write_Atom() : 未実装[" << ATOM_STYLE << "]" << endl;
    exit(0);
  }
}

void Atom::output_Atom(const string& ATOM_STYLE)
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  if (ATOM_STYLE == "full")
  {
    cout << setpre(DIGITS) << atom_ID << " " << molecule_ID << " " << atom_type << " " << q << " " << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << endl;
  }
  else if (ATOM_STYLE == "atomic")
  {
    cout << setpre(DIGITS) << atom_ID << " " << atom_type << " " << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << endl;
  }
  else
  {
    cerr << "Atom::output_Atom() : 未実装[" << ATOM_STYLE << "]" << endl;
    exit(0);
  }
}

void Atom::write_Velocity(const string& ATOM_STYLE, ofstream& ofs)
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  if (ATOM_STYLE == "electron" || ATOM_STYLE == "ellipsoid" ||\
      ATOM_STYLE == "splere"   || ATOM_STYLE == "hybrid")
  {
    cerr << "Atom::write_Velocity() : 未実装[" << ATOM_STYLE << "]" << endl;
    exit(0);
  }
  else
  {
    ofs << setpre(DIGITS) << atom_ID << " " << vx << " " << vy << " " << vz << endl;
  }
}

void Atom::output_Velocity(const string& ATOM_STYLE)
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  if (ATOM_STYLE == "electron" || ATOM_STYLE == "ellipsoid" ||\
      ATOM_STYLE == "splere"   || ATOM_STYLE == "hybrid")
  {
    cerr << "Atom::output_Velocity() : 未実装[" << ATOM_STYLE << "]" << endl;
    exit(0);
  }
  else
  {
    cout << setpre(DIGITS) << atom_ID << " " << vx << " " << vy << " " << vz << endl;
  }
}

Mass::Mass(): atom_type(0), mass(0.0){}

void Mass::write(ofstream& ofs)
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  ofs << setpre(DIGITS) << atom_type << " " << mass << endl;
}

void Mass::output()
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  cout << setpre(DIGITS) << atom_type << " " << mass << endl;
}

Bond::Bond(): atom_num(2), bond_ID(0), bond_type(0), atom(atom_num){}

void Bond::write(ofstream& ofs)
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  ofs << bond_ID << " " << bond_type << " ";
  for (int i = 0; i < atom_num; ++i) ofs << atom[i]->atom_ID << " ";
  ofs << endl;
}

void Bond::output()
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  cout << bond_ID << " " << bond_type << " ";
  for (int i = 0; i < atom_num; ++i) cout << atom[i]->atom_ID << " ";
  cout << endl;
}

Angle::Angle(): atom_num(3), angle_ID(0), angle_type(0), atom(atom_num){}

void Angle::write(ofstream& ofs)
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  ofs << angle_ID << " " << angle_type << " ";
  for (int i = 0; i < atom_num; ++i) ofs << atom[i]->atom_ID << " ";
  ofs << endl;
}

void Angle::output()
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  cout << angle_ID << " " << angle_type << " ";
  for (int i = 0; i < atom_num; ++i) cout << atom[i]->atom_ID << " ";
  cout << endl;
}

Dihedral::Dihedral(): atom_num(4), dihedral_ID(0), dihedral_type(0), atom(atom_num){}

void Dihedral::write(ofstream& ofs)
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  ofs << dihedral_ID << " " << dihedral_type << " ";
  for (int i = 0; i < atom_num; ++i) ofs << atom[i]->atom_ID << " ";
  ofs << endl;
}

void Dihedral::output()
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  cout << dihedral_ID << " " << dihedral_type << " ";
  for (int i = 0; i < atom_num; ++i) cout << atom[i]->atom_ID << " ";
  cout << endl;
}

Improper::Improper(): atom_num(4), improper_ID(0), improper_type(0), atom(atom_num){}

void Improper::write(ofstream& ofs)
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  ofs << improper_ID << " " << improper_type << " ";
  for (int i = 0; i < atom_num; ++i) ofs << atom[i]->atom_ID << " ";
  ofs << endl;
}

void Improper::output()
{
  // if (DEBUG) cout << __FUNCTION__ << " Called" << endl;
  cout << improper_ID << " " << improper_type << " ";
  for (int i = 0; i < atom_num; ++i) cout << atom[i]->atom_ID << " ";
  cout << endl;
}

Pair_Coeff::Pair_Coeff(): atom_type(0), coeffs(""){}

void Pair_Coeff::write(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  ofs << atom_type << " " << coeffs << endl;
}

void Pair_Coeff::output()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << atom_type << " " << coeffs << endl;
}

Bond_Coeff::Bond_Coeff(): bond_type(0), coeffs(""){}

void Bond_Coeff::write(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  ofs << bond_type << " " << coeffs << endl;
}

void Bond_Coeff::output()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << bond_type << " " << coeffs << endl;
}

Angle_Coeff::Angle_Coeff(): angle_type(0), coeffs(""){}

void Angle_Coeff::write(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  ofs << angle_type << " " << coeffs << endl;
}

void Angle_Coeff::output()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << angle_type << " " << coeffs << endl;
}

Dihedral_Coeff::Dihedral_Coeff(): dihedral_type(0), coeffs(""){}

void Dihedral_Coeff::write(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  ofs << dihedral_type << " " << coeffs << endl;
}

void Dihedral_Coeff::output()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << dihedral_type << " " << coeffs << endl;
}

Improper_Coeff::Improper_Coeff(): improper_type(0), coeffs(""){}

void Improper_Coeff::write(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  ofs << improper_type << " " << coeffs << endl;
}

void Improper_Coeff::output()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << improper_type << " " << coeffs << endl;
}

// Constructor
Data::Data(): dimension(0), atom_style("atomic"), pair_coeff_style(""), bond_coeff_style(""), angle_coeff_style(""), dihedral_coeff_style(""), improper_coeff_style(""), filename(""), line_header("")
{
  initialize_parameters();
}

// Initialize
void Data::initialize_parameters()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  atoms = 0; bonds = 0; angles = 0; dihedrals = 0; impropers = 0; atom_types = 0; bond_types = 0; angle_types = 0; dihedral_types = 0; improper_types = 0;
  extra_bond_per_atom = 0; extra_angle_per_atom = 0; extra_dihedral_per_atom = 0; extra_improper_per_atom = 0; extra_special_per_atom = 0;
  ellipsoids = 0; lines = 0; triangles = 0; bodies = 0;
  xlo = -0.5; xhi = 0.5; ylo = -0.5; yhi = 0.5; zlo = -0.5; zhi = 0.5;
  xy = 0.0; xz = 0.0; yz = 0.0;
}

// Get Header
void Data::get_atoms(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  atoms = stoi(tmp);
}

void Data::get_atom_types(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  atom_types = stoi(tmp);
}

void Data::get_bonds(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  bonds = stoi(tmp);
}

void Data::get_bond_types(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  bond_types = stoi(tmp);
}

void Data::get_angles(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  angles = stoi(tmp);
}

void Data::get_angle_types(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  angle_types = stoi(tmp);
}

void Data::get_dihedrals(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  dihedrals = stoi(tmp);
}

void Data::get_dihedral_types(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  dihedral_types = stoi(tmp);
}

void Data::get_impropers(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  impropers = stoi(tmp);
}

void Data::get_improper_types(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp;
  istringstream iss(line);
  iss >> tmp;
  improper_types = stoi(tmp);
}

void Data::get_size_x(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp1, tmp2;
  istringstream iss(line);
  iss >> tmp1 >> tmp2;
  xlo = convert(tmp1);
  xhi = convert(tmp2);
  ++dimension;
}

void Data::get_size_y(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp1, tmp2;
  istringstream iss(line);
  iss >> tmp1 >> tmp2;
  ylo = convert(tmp1);
  yhi = convert(tmp2);
  ++dimension;
}

void Data::get_size_z(const string& line)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string tmp1, tmp2;
  istringstream iss(line);
  iss >> tmp1 >> tmp2;
  zlo = convert(tmp1);
  zhi = convert(tmp2);
  ++dimension;
}

// Get Body
void Data::get_Masses(const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string atom_type; // int
  string mass; // double

  for (int i = start+2; i < start+2+atom_types; ++i)
  {
    istringstream iss(data[i]);
    iss >> atom_type >> mass;
    Masses[stoi(atom_type)].atom_type = stoi(atom_type);
    Masses[stoi(atom_type)].mass = convert(mass);
  }
}

string Data::get_coeffs(const string& line)
{
  int flag = 0;
  string coeffs = "";
  for (char c: line)
  {
    if (flag == 1) coeffs += c;
    else if (c == ' ') ++flag;
  }
  return coeffs;
}

void Data::get_Pair_Coeffs(const string& line, const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string atom_type; //int
  string coeffs;  //string

  int flag = 0;
  pair_coeff_style = "";
  for (char c: line)
  {
    if (c == '#') ++flag;
    else if (flag == 1) ++flag;
    else if (flag == 2) pair_coeff_style += c;
  }
  for (int i = start+2; i < start+2+atom_types; ++i)
  {
    istringstream iss(data[i]);
    iss >> atom_type;
    coeffs = get_coeffs(data[i]);
    Pair_Coeffs[stoi(atom_type)].atom_type = stoi(atom_type);
    Pair_Coeffs[stoi(atom_type)].coeffs = coeffs;
  }
}

void Data::get_Bond_Coeffs(const string& line, const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string bond_type; //int
  string coeffs;  //string

  int flag = 0;
  bond_coeff_style = "";
  for (char c: line)
  {
    if (c == '#') ++flag;
    else if (flag == 1) ++flag;
    else if (flag == 2) bond_coeff_style += c;
  }
  for (int i = start+2; i < start+2+bond_types; ++i)
  {
    istringstream iss(data[i]);
    iss >> bond_type;
    coeffs = get_coeffs(data[i]);
    Bond_Coeffs[stoi(bond_type)].bond_type = stoi(bond_type);
    Bond_Coeffs[stoi(bond_type)].coeffs = coeffs;
  }
}

void Data::get_Angle_Coeffs(const string& line, const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string angle_type; //int
  string coeffs;  //string

  int flag = 0;
  angle_coeff_style = "";
  for (char c: line)
  {
    if (c == '#') ++flag;
    else if (flag == 1) ++flag;
    else if (flag == 2) angle_coeff_style += c;
  }
  for (int i = start+2; i < start+2+angle_types; ++i)
  {
    istringstream iss(data[i]);
    iss >> angle_type;
    coeffs = get_coeffs(data[i]);
    Angle_Coeffs[stoi(angle_type)].angle_type = stoi(angle_type);
    Angle_Coeffs[stoi(angle_type)].coeffs = coeffs;
  }
}

void Data::get_Dihedral_Coeffs(const string& line, const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string dihedral_type; //int
  string coeffs;  //string

  int flag = 0;
  dihedral_coeff_style = "";
  for (char c: line)
  {
    if (c == '#') ++flag;
    else if (flag == 1) ++flag;
    else if (flag == 2) dihedral_coeff_style += c;
  }
  for (int i = start+2; i < start+2+dihedral_types; ++i)
  {
    istringstream iss(data[i]);
    iss >> dihedral_type;
    coeffs = get_coeffs(data[i]);
    Dihedral_Coeffs[stoi(dihedral_type)].dihedral_type = stoi(dihedral_type);
    Dihedral_Coeffs[stoi(dihedral_type)].coeffs = coeffs;
  }
}

void Data::get_Improper_Coeffs(const string& line, const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string improper_type; //int
  string coeffs;  //string

  int flag = 0;
  improper_coeff_style = "";
  for (char c: line)
  {
    if (c == '#') ++flag;
    else if (flag == 1) ++flag;
    else if (flag == 2) improper_coeff_style += c;
  }
  for (int i = start+2; i < start+2+improper_types; ++i)
  {
    istringstream iss(data[i]);
    iss >> improper_type;
    coeffs = get_coeffs(data[i]);
    Improper_Coeffs[stoi(improper_type)].improper_type = stoi(improper_type);
    Improper_Coeffs[stoi(improper_type)].coeffs = coeffs;
  }
}

void Data::get_Atoms(const string& line, const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string atom_ID, molecule_ID, atom_type, etag, spin, spx, spy, spz, template_atom, template_index; //int
  string ccN, contact_radius, cs_re, cs_im, cv, density, diameter, espf, edpd_temp, edpd_cv, eradius, kernel_radius, mass, mux, muy, muz, q, rho, sp, theta, volume, x, y, z, x0, y0, z0; //double
  string ellipsoidflag, lineflag, triangleflag, nx, ny, nz; //bool

  int flag = 0;
  atom_style = "";
  for (char c: line)
  {
    if (c == '#') ++flag;
    else if (flag == 1) ++flag;
    else if (flag == 2) atom_style += c;
  }
  for (int i = start+2; i < start+2+atoms; ++i)
  {
    istringstream iss(data[i]);
    if (atom_style == "atomic")
    {
      iss >> atom_ID >> atom_type >> x >> y >> z >> nx >> ny >> nz;
      Atoms[stoi(atom_ID)].atom_ID = stoi(atom_ID);
      Atoms[stoi(atom_ID)].atom_type = stoi(atom_type);
      Atoms[stoi(atom_ID)].x = convert(x);
      Atoms[stoi(atom_ID)].y = convert(y);
      Atoms[stoi(atom_ID)].z = convert(z);
      Atoms[stoi(atom_ID)].nx = stoi(nx);
      Atoms[stoi(atom_ID)].ny = stoi(ny);
      Atoms[stoi(atom_ID)].nz = stoi(nz);
      Atoms[stoi(atom_ID)].mass = Masses.at(Atoms[stoi(atom_ID)].atom_type).mass;
    }
    else if (atom_style == "full")
    {
      iss >> atom_ID >> molecule_ID >> atom_type >> q >> x >> y >> z >> nx >> ny >> nz;
      Atoms[stoi(atom_ID)].atom_ID = stoi(atom_ID);
      Atoms[stoi(atom_ID)].molecule_ID = stoi(molecule_ID);
      Atoms[stoi(atom_ID)].atom_type = stoi(atom_type);
      Atoms[stoi(atom_ID)].q = convert(q);
      Atoms[stoi(atom_ID)].x = convert(x);
      Atoms[stoi(atom_ID)].y = convert(y);
      Atoms[stoi(atom_ID)].z = convert(z);
      Atoms[stoi(atom_ID)].nx = stoi(nx);
      Atoms[stoi(atom_ID)].ny = stoi(ny);
      Atoms[stoi(atom_ID)].nz = stoi(nz);
      Atoms[stoi(atom_ID)].mass = Masses.at(Atoms[stoi(atom_ID)].atom_type).mass;
    }
    else
    {
      cerr << "Data::get_Atoms() : 未実装[" << atom_style << "]" << endl;
      exit(0);
    }
  }
}

void Data::get_Velocities(const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string atom_ID; //int
  string vx, vy, vz, lx, ly, lz, wx, wy, wz, ervel;  //double

  for (int i = start+2; i < start+2+atoms; ++i)
  {
    istringstream iss(data[i]);
    if (atom_style == "electron" || atom_style == "ellipsoid" ||\
        atom_style == "splere"   || atom_style == "hybrid")
    {
      cerr << "Data::get_Velocities() : 未実装[" << atom_style << "]" << endl;
      exit(0);
    }
    else
    {
      iss >> atom_ID >> vx >> vy >> vz;
      Atoms.at(stoi(atom_ID)).vx = convert(vx);
      Atoms.at(stoi(atom_ID)).vy = convert(vy);
      Atoms.at(stoi(atom_ID)).vz = convert(vz);
    }
  }
}

void Data::get_Bonds(const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int atom_num = 2;
  string bond_ID, bond_type;  // int
  vector<string> atom(atom_num); // int

  for (int i = start+2; i < start+2+bonds; ++i)
  {
    istringstream iss(data[i]);
    if (atom_style == "atomic")
    {
      cerr << endl << "Data::get_Bonds() : " << atom_style << " Does Not Have Bonds" << endl;
      exit(0);
    }
    else if (atom_style == "full")
    {
      iss >> bond_ID >> bond_type;
      for (int j = 0; j < atom_num; ++j) iss >> atom[j];
      Bonds[stoi(bond_ID)].bond_ID = stoi(bond_ID);
      Bonds[stoi(bond_ID)].bond_type = stoi(bond_type);
      for (int j = 0; j < atom_num; ++j)
      {
        Atoms.at(stoi(atom[j])).have_bond = true;
        Bonds.at(stoi(bond_ID)).atom[j] = &Atoms.at(stoi(atom[j]));
      }
    }
    else
    {
      cerr << "Data::get_Bonds() : 未実装[" << atom_style << "]" << endl;
      exit(0);
    }
  }
}

void Data::get_Angles(const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int atom_num = 3;
  string angle_ID, angle_type;  // int
  vector<string> atom(atom_num); // int

  for (int i = start+2; i < start+2+angles; ++i)
  {
    istringstream iss(data[i]);
    if (atom_style == "atomic")
    {
      cerr << endl << "Data::get_Angles() : " << atom_style << " Does Not Have Angles" << endl;
      exit(0);
    }
    else if (atom_style == "full")
    {
      iss >> angle_ID >> angle_type;
      for (int j = 0; j < atom_num; ++j) iss >> atom[j];
      Angles[stoi(angle_ID)].angle_ID = stoi(angle_ID);
      Angles[stoi(angle_ID)].angle_type = stoi(angle_type);
      for (int j = 0; j < atom_num; ++j) Angles.at(stoi(angle_ID)).atom[j] = &Atoms.at(stoi(atom[j]));
    }
    else
    {
      cerr << "Data::get_Angles() : 未実装[" << atom_style << "]" << endl;
      exit(0);
    }
  }
}

void Data::get_Dihedrals(const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int atom_num = 4;
  string dihedral_ID, dihedral_type;  // int
  vector<string> atom(atom_num); // int

  for (int i = start+2; i < start+2+dihedrals; ++i)
  {
    istringstream iss(data[i]);
    if (atom_style == "atomic")
    {
      cerr << endl << "Data::get_Dihedrals() : " << atom_style << " Does Not Have Dihedrals" << endl;
      exit(0);
    }
    else if (atom_style == "full")
    {
      iss >> dihedral_ID >> dihedral_type;
      for (int j = 0; j < atom_num; ++j) iss >> atom[j];
      Dihedrals[stoi(dihedral_ID)].dihedral_ID = stoi(dihedral_ID);
      Dihedrals[stoi(dihedral_ID)].dihedral_type = stoi(dihedral_type);
      for (int j = 0; j < atom_num; ++j) Dihedrals.at(stoi(dihedral_ID)).atom[j] = &Atoms.at(stoi(atom[j]));
    }
    else
    {
      cerr << "Data::get_Dihedrals() : 未実装[" << atom_style << "]" << endl;
      exit(0);
    }
  }
}

void Data::get_Impropers(const int& start)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int atom_num = 4;
  string improper_ID, improper_type;  // int
  vector<string> atom(atom_num); // int

  for (int i = start+2; i < start+2+impropers; ++i)
  {
    istringstream iss(data[i]);
    if (atom_style == "atomic")
    {
      cerr << endl << "Data::get_Impropers() : " << atom_style << " Does Not Have Impropers" << endl;
      exit(0);
    }
    else if (atom_style == "full")
    {
      iss >> improper_ID >> improper_type;
      for (int j = 0; j < atom_num; ++j) iss >> atom[j];
      Impropers[stoi(improper_ID)].improper_ID = stoi(improper_ID);
      Impropers[stoi(improper_ID)].improper_type = stoi(improper_type);
      for (int j = 0; j < atom_num; ++j) Impropers.at(stoi(improper_ID)).atom[j] = &Atoms.at(stoi(atom[j]));
    }
    else
    {
      cerr << "Data::get_Impropers() : 未実装[" << atom_style << "]" << endl;
      exit(0);
    }
  }
}

// Read
string Data::check_extension(const string& FILENAME)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  string ret = "";
  for (char c: FILENAME)
  {
    if (c == '.') return FILENAME;
    else ret += c;
  }
  ret += ".data";
  return ret;
}

void Data::get_filename(const string& FILENAME)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (FILENAME == "")
  {
    cout << "-------------------------------------------------------------------------------------------------" << endl;
    cout << "The name of Data File" << endl;
    cout << ">> ";
    cin >> filename;
    filename = check_extension(filename);
    ifstream ifs(filename);
    bool flg = ifs.fail();
    while (flg)
    {
      cerr << "Failed to open file." << endl;
      cout << "The name of Data File" << endl;
      cout << ">> ";
      cin >> filename;
      filename = check_extension(filename);
      ifstream ifs2(filename);
      flg = ifs2.fail();
    }
  }
  else filename = FILENAME;
}

void Data::get_data(const string& FILENAME)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  get_filename(FILENAME);
  ifstream ifs(filename);
  if (ifs.fail())
  {
    cerr << endl << "Cannot Open File : " << filename << endl;
    exit(0);
  }
  string line;
  while(getline(ifs, line)) data.emplace_back(line);
  line_header = data[0];
}

void Data::read_data()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int line_number = 0;
  for (string line: data)
  {
    if (is_included(line, " atoms")) get_atoms(line);
    else if (is_included(line, " atom types")) get_atom_types(line);
    else if (is_included(line, " bonds")) get_bonds(line);
    else if (is_included(line, " bond types")) get_bond_types(line);
    else if (is_included(line, " angles")) get_angles(line);
    else if (is_included(line, " angle types")) get_angle_types(line);
    else if (is_included(line, " dihedrals")) get_dihedrals(line);
    else if (is_included(line, " dihedral types")) get_dihedral_types(line);
    else if (is_included(line, " impropers")) get_impropers(line);
    else if (is_included(line, " improper types")) get_improper_types(line);
    else if (is_included(line, " xlo xhi")) get_size_x(line);
    else if (is_included(line, " ylo yhi")) get_size_y(line);
    else if (is_included(line, " zlo zhi")) get_size_z(line);
    else if (is_included(line, "Masses")) get_Masses(line_number);
    else if (is_included(line, "Pair Coeffs")) get_Pair_Coeffs(line, line_number);
    else if (is_included(line, "Bond Coeffs")) get_Bond_Coeffs(line, line_number);
    else if (is_included(line, "Angle Coeffs")) get_Angle_Coeffs(line, line_number);
    else if (is_included(line, "Dihedral Coeffs")) get_Dihedral_Coeffs(line, line_number);
    else if (is_included(line, "Improper Coeffs")) get_Improper_Coeffs(line, line_number);
    else if (is_included(line, "Atoms")) get_Atoms(line, line_number);
    else if (is_included(line, "Velocities")) get_Velocities(line_number);
    else if (is_included(line, "Bonds")) get_Bonds(line_number);
    else if (is_included(line, "Angles")) get_Angles(line_number);
    else if (is_included(line, "Dihedrals")) get_Dihedrals(line_number);
    else if (is_included(line, "Impropers")) get_Impropers(line_number);
    ++line_number;
  }
}

void Data::Read(const string& FILENAME)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  get_data(FILENAME);
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  cout << "Read File : " << filename << endl;
  read_data();
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  cout << "atom style : " << atom_style << endl;
  output_header();
}

// Write
void Data::Write(string FILENAME)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (FILENAME == "") FILENAME = "Data.out";
  ofstream ofs(FILENAME);
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  cout << "Write File : " << FILENAME << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  cout << "atom style : " << atom_style << endl;
  write_header(ofs);
  output_header();
  write_body(ofs);
}

void Data::write_header(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  ofs << line_header << endl << endl;
  ofs << atoms << " atoms" << endl;
  ofs << atom_types << " atom_types" << endl;
  if (atom_style == "full")
  {
    ofs << bonds << " bonds" << endl;
    ofs << bond_types << " bond_types" << endl;
    ofs << angles << " angles" << endl;
    ofs << angle_types << " angle_types" << endl;
    ofs << dihedrals << " dihedrals" << endl;
    ofs << dihedral_types << " dihedral_types" << endl;
    ofs << impropers << " impropers" << endl;
    ofs << improper_types << " improper_types" << endl;
  }
  else if (atom_style != "atomic")
  {
    cerr << "Data::wrtie_header() : 未実装[" << atom_style << "]" << endl;
    exit(0);
  }
  ofs << endl;
  if (dimension >= 1) ofs << setpre(DIGITS) << xlo << " " << xhi << " xlo xhi" << endl;
  if (dimension >= 2) ofs << setpre(DIGITS) << ylo << " " << yhi << " ylo yhi" << endl;
  if (dimension >= 3) ofs << setpre(DIGITS) << zlo << " " << zhi << " zlo zhi" << endl;
  ofs << endl;
}

void Data::write_body(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  write_Masses(ofs);
  write_Pair_Coeffs(ofs);
  write_Bond_Coeffs(ofs);
  write_Angle_Coeffs(ofs);
  write_Dihedral_Coeffs(ofs);
  write_Improper_Coeffs(ofs);
  write_Atoms(ofs);
  write_Velocities(ofs);
  write_Bonds(ofs);
  write_Angles(ofs);
  write_Dihedrals(ofs);
  write_Impropers(ofs);
}

void Data::write_Masses(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  ofs << "Masses" << endl << endl;
  for (int i = 0; i < atom_types; ++i) Masses.at(i+1).write(ofs);
  ofs << endl;
}

void Data::write_Pair_Coeffs(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (atom_style == "full" && pair_coeff_style != "")
  {
    ofs << "Pair Coeffs # " << pair_coeff_style << endl << endl;
    for (int i = 0; i < atom_types; ++i) Pair_Coeffs.at(i+1).write(ofs);
  }
  ofs << endl;
}

void Data::write_Bond_Coeffs(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (atom_style == "full" && bond_coeff_style != "")
  {
    ofs << "Bond Coeffs # " << bond_coeff_style << endl << endl;
    for (int i = 0; i < bond_types; ++i) Bond_Coeffs.at(i+1).write(ofs);
  }
  ofs << endl;
}

void Data::write_Angle_Coeffs(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (atom_style == "full" && angle_coeff_style != "")
  {
    ofs << "Angle Coeffs # " << angle_coeff_style << endl << endl;
    for (int i = 0; i < angle_types; ++i) Angle_Coeffs.at(i+1).write(ofs);
  }
  ofs << endl;
}

void Data::write_Dihedral_Coeffs(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (atom_style == "full" && dihedral_coeff_style != "")
  {
    ofs << "Dihedral Coeffs # " << dihedral_coeff_style << endl << endl;
    for (int i = 0; i < dihedral_types; ++i) Dihedral_Coeffs.at(i+1).write(ofs);
  }
  ofs << endl;
}

void Data::write_Improper_Coeffs(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (atom_style == "full" && improper_coeff_style != "")
  {
    ofs << "Improper Coeffs # " << improper_coeff_style << endl << endl;
    for (int i = 0; i < improper_types; ++i) Improper_Coeffs.at(i+1).write(ofs);
  }
  ofs << endl;
}

void Data::write_Atoms(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  ofs << "Atoms # " << atom_style << endl << endl;
  for (int i = 0; i < atoms; ++i) Atoms.at(i+1).write_Atom(atom_style, ofs);
  ofs << endl;
}

void Data::write_Velocities(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  ofs << "Velocities" << endl << endl;
  for (int i = 0; i < atoms; ++i) Atoms.at(i+1).write_Velocity(atom_style, ofs);
  ofs << endl;
}

void Data::write_Bonds(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (atom_style == "full" && bonds > 0)
  {
    ofs << "Bonds" << endl << endl;
    for (int i = 0; i < bonds; ++i) Bonds.at(i+1).write(ofs);
    ofs << endl;
  }
}

void Data::write_Angles(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (atom_style == "full" && angles > 0)
  {
    ofs << "Angles" << endl << endl;
    for (int i = 0; i < angles; ++i) Angles.at(i+1).write(ofs);
    ofs << endl;
  }
}

void Data::write_Dihedrals(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (atom_style == "full" && dihedrals > 0)
  {
    ofs << "Dihedrals" << endl << endl;
    for (int i = 0; i < dihedrals; ++i) Dihedrals.at(i+1).write(ofs);
    ofs << endl;
  }
}

void Data::write_Impropers(ofstream& ofs)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  if (atom_style == "full" && impropers > 0)
  {
    ofs << "Impropers" << endl << endl;
    for (int i = 0; i < impropers; ++i) Impropers.at(i+1).write(ofs);
    ofs << endl;
  }
}

// Output
void Data::Output()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  output_header();
  output_body();
}

void Data::output_header()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  cout << line_header << endl << endl;
  cout << atoms << " atoms" << endl;
  cout << atom_types << " atom_types" << endl;
  if (atom_style == "full")
  {
    cout << bonds << " bonds" << endl;
    cout << bond_types << " bond_types" << endl;
    cout << angles << " angles" << endl;
    cout << angle_types << " angle_types" << endl;
    cout << dihedrals << " dihedrals" << endl;
    cout << dihedral_types << " dihedral_types" << endl;
    cout << impropers << " impropers" << endl;
    cout << improper_types << " improper_types" << endl;
  }
  else if (atom_style != "atomic")
  {
    cerr << "Data::wrtie_header() : 未実装[" << atom_style << "]" << endl;
    exit(0);
  }
  cout << endl;
  if (dimension >= 1) cout << setpre(DIGITS) << xlo << " " << xhi << " xlo xhi" << endl;
  if (dimension >= 2) cout << setpre(DIGITS) << ylo << " " << yhi << " ylo yhi" << endl;
  if (dimension >= 3) cout << setpre(DIGITS) << zlo << " " << zhi << " zlo zhi" << endl;
  cout << endl;
}

void Data::output_body()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  output_Masses();
  output_Pair_Coeffs();
  output_Bond_Coeffs();
  output_Angle_Coeffs();
  output_Dihedral_Coeffs();
  output_Improper_Coeffs();
  output_Atoms();
  output_Velocities();
  output_Bonds();
  output_Angles();
  output_Dihedrals();
  output_Impropers();
}

void Data::output_Masses()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  cout << "Masses" << endl << endl;
  for (int i = 0; i < atom_types; ++i) Masses.at(i+1).output();
  cout << endl;
}

void Data::output_Pair_Coeffs()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  if (atom_style == "full" && pair_coeff_style != "")
  {
    cout << "Pair Coeffs # " << pair_coeff_style << endl << endl;
    for (int i = 0; i < atom_types; ++i) Pair_Coeffs.at(i+1).output();
  }
  cout << endl;
}

void Data::output_Bond_Coeffs()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  if (atom_style == "full" && bond_coeff_style != "")
  {
    cout << "Bond Coeffs # " << bond_coeff_style << endl << endl;
    for (int i = 0; i < bond_types; ++i) Bond_Coeffs.at(i+1).output();
  }
  cout << endl;
}

void Data::output_Angle_Coeffs()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  if (atom_style == "full" && angle_coeff_style != "")
  {
    cout << "Angle Coeffs # " << angle_coeff_style << endl << endl;
    for (int i = 0; i < angle_types; ++i) Angle_Coeffs.at(i+1).output();
  }
  cout << endl;
}

void Data::output_Dihedral_Coeffs()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  if (atom_style == "full" && dihedral_coeff_style != "")
  {
    cout << "Dihedral Coeffs # " << dihedral_coeff_style << endl << endl;
    for (int i = 0; i < dihedral_types; ++i) Dihedral_Coeffs.at(i+1).output();
  }
  cout << endl;
}

void Data::output_Improper_Coeffs()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  if (atom_style == "full" && improper_coeff_style != "")
  {
    cout << "Improper Coeffs # " << improper_coeff_style << endl << endl;
    for (int i = 0; i < improper_types; ++i) Improper_Coeffs.at(i+1).output();
  }
  cout << endl;
}

void Data::output_Atoms()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  cout << "Atoms # " << atom_style << endl << endl;
  for (int i = 0; i < atoms; ++i) Atoms.at(i+1).output_Atom(atom_style);
  cout << endl;
}

void Data::output_Velocities()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  cout << "Velocities" << endl << endl;
  for (int i = 0; i < atoms; ++i) Atoms.at(i+1).output_Velocity(atom_style);
  cout << endl;
}

void Data::output_Bonds()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  if (atom_style == "full" && bonds > 0)
  {
    cout << "Bonds" << endl << endl;
    for (int i = 0; i < bonds; ++i) Bonds.at(i+1).output();
    cout << endl;
  }
}

void Data::output_Angles()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  if (atom_style == "full" && angles > 0)
  {
    cout << "Angles" << endl << endl;
    for (int i = 0; i < angles; ++i) Angles.at(i+1).output();
    cout << endl;
  }
}

void Data::output_Dihedrals()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  if (atom_style == "full" && dihedrals > 0)
  {
    cout << "Dihedrals" << endl << endl;
    for (int i = 0; i < dihedrals; ++i) Dihedrals.at(i+1).output();
    cout << endl;
  }
}

void Data::output_Impropers()
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  cout << "-------------------------------------------------------------------------------------------------" << endl;
  if (atom_style == "full" && impropers > 0)
  {
    cout << "Impropers" << endl << endl;
    for (int i = 0; i < impropers; ++i) Impropers.at(i+1).output();
    cout << endl;
  }
}

// Functions
void Data::Remove_Molecule_ID(const vector<int>& MOLECULE_ID)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int id;
  bool flag;
  int tmp_atoms = atoms, tmp_bonds = bonds, tmp_angles = angles, tmp_dihedrals = dihedrals, tmp_impropers = impropers;
  int atom_num_bond = 2, atom_num_angle = 3, atom_num_dihedral = 4, atom_num_improper = 4;
  map<int, Atom> New_Atoms;
  map<int, Bond> New_Bonds;
  map<int, Angle> New_Angles;
  map<int, Dihedral> New_Dihedrals;
  map<int, Improper> New_Impropers;

  id = 0;
  for (int i = 0; i < tmp_atoms; ++i)
  {
    flag = false;
    for (int ID: MOLECULE_ID)
    {
      if (Atoms.at(i+1).molecule_ID == ID)
      {
        --atoms;
        flag = true;
        break;
      }
    }
    if (flag) continue;
    ++id;
    New_Atoms[id] = Atoms.at(i+1);
    New_Atoms[id].atom_ID = id;
    Atoms.at(i+1) = New_Atoms.at(id);
  }

  id = 0;
  for (int i = 0; i < tmp_bonds; ++i)
  {
    flag = false;
    for (int ID: MOLECULE_ID)
    {
      for (int j = 0; j < atom_num_bond; ++j)
      {
        if (Bonds.at(i+1).atom[j]->molecule_ID == ID)
        {
          flag = true;
          --bonds;
          break;
        }
      }
      if (flag) break;
    }
    if (flag) continue;
    ++id;
    New_Bonds[id] = Bonds.at(i+1);
    New_Bonds[id].bond_ID = id;
    for (int j = 0; j < atom_num_bond; ++j)
    {
      New_Bonds.at(id).atom[j] = &New_Atoms.at(New_Bonds.at(id).atom[j]->atom_ID);
    }
  }

  id = 0;
  for (int i = 0; i < tmp_angles; ++i)
  {
    flag = false;
    for (int ID: MOLECULE_ID)
    {
      for (int j = 0; j < atom_num_angle; ++j)
      {
        if (Angles.at(i+1).atom[j]->molecule_ID == ID)
        {
          flag = true;
          --angles;
          break;
        }
      }
      if (flag) break;
    }
    if (flag) continue;
    ++id;
    New_Angles[id] = Angles.at(i+1);
    New_Angles[id].angle_ID = id;
    for (int j = 0; j < atom_num_angle; ++j)
    {
      New_Angles.at(id).atom[j] = &New_Atoms.at(New_Angles.at(id).atom[j]->atom_ID);
    }
  }

  id = 0;
  for (int i = 0; i < tmp_dihedrals; ++i)
  {
    flag = false;
    for (int ID: MOLECULE_ID)
    {
      for (int j = 0; j < atom_num_dihedral; ++j)
      {
        if (Dihedrals.at(i+1).atom[j]->molecule_ID == ID)
        {
          flag = true;
          --dihedrals;
          break;
        }
      }
      if (flag) break;
    }
    if (flag) continue;
    ++id;
    New_Dihedrals[id] = Dihedrals.at(i+1);
    New_Dihedrals[id].dihedral_ID = id;
    for (int j = 0; j < atom_num_dihedral; ++j)
    {
      New_Dihedrals.at(id).atom[j] = &New_Atoms.at(New_Dihedrals.at(id).atom[j]->atom_ID);
    }
  }

  id = 0;
  for (int i = 0; i < tmp_impropers; ++i)
  {
    flag = false;
    for (int ID: MOLECULE_ID)
    {
      for (int j = 0; j < atom_num_improper; ++j)
      {
        if (Impropers.at(i+1).atom[j]->molecule_ID == ID)
        {
          flag = true;
          --impropers;
          break;
        }
      }
      if (flag) break;
    }
    if (flag) continue;
    ++id;
    New_Impropers[id] = Impropers.at(i+1);
    New_Impropers[id].improper_ID = id;
    for (int j = 0; j < atom_num_improper; ++j)
    {
      New_Impropers.at(id).atom[j] = &New_Atoms.at(New_Impropers.at(id).atom[j]->atom_ID);
    }
  }

  Atoms = New_Atoms;
  Bonds = New_Bonds;
  Angles = New_Angles;
  Dihedrals = New_Dihedrals;
  Impropers = New_Impropers;
  for (int i = 0; i < bonds; ++i)
  {
    for (int j = 0; j < atom_num_bond; ++j)
    {
      Bonds.at(i+1).atom[j] = &Atoms.at(Bonds.at(i+1).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < angles; ++i)
  {
    for (int j = 0; j < atom_num_angle; ++j)
    {
      Angles.at(i+1).atom[j] = &Atoms.at(Angles.at(i+1).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < dihedrals; ++i)
  {
    for (int j = 0; j < atom_num_dihedral; ++j)
    {
      Dihedrals.at(i+1).atom[j] = &Atoms.at(Dihedrals.at(i+1).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < impropers; ++i)
  {
    for (int j = 0; j < atom_num_improper; ++j)
    {
      Impropers.at(i+1).atom[j] = &Atoms.at(Impropers.at(i+1).atom[j]->atom_ID);
    }
  }
}

void Data::Remove_Atom_Type(const vector<int>& ATOM_TYPE)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int id;
  bool flag;
  int tmp_atoms = atoms, tmp_bonds = bonds, tmp_angles = angles, tmp_dihedrals = dihedrals, tmp_impropers = impropers;
  int atom_num_bond = 2, atom_num_angle = 3, atom_num_dihedral = 4, atom_num_improper = 4;
  map<int, Atom> New_Atoms;
  map<int, Bond> New_Bonds;
  map<int, Angle> New_Angles;
  map<int, Dihedral> New_Dihedrals;
  map<int, Improper> New_Impropers;

  id = 0;
  for (int i = 0; i < tmp_atoms; ++i)
  {
    flag = false;
    for (int TYPE: ATOM_TYPE)
    {
      if (Atoms.at(i+1).atom_type == TYPE)
      {
        --atoms;
        flag = true;
        break;
      }
    }
    if (flag) continue;
    ++id;
    New_Atoms[id] = Atoms.at(i+1);
    New_Atoms[id].atom_ID = id;
    Atoms.at(i+1) = New_Atoms.at(id);
  }

  id = 0;
  for (int i = 0; i < tmp_bonds; ++i)
  {
    flag = false;
    for (int TYPE: ATOM_TYPE)
    {
      for (int j = 0; j < atom_num_bond; ++j)
      {
        if (Bonds.at(i+1).atom[j]->atom_type == TYPE)
        {
          flag = true;
          --bonds;
          break;
        }
      }
      if (flag) break;
    }
    if (flag) continue;
    ++id;
    New_Bonds[id] = Bonds.at(i+1);
    New_Bonds[id].bond_ID = id;
    for (int j = 0; j < atom_num_bond; ++j)
    {
      New_Bonds.at(id).atom[j] = &New_Atoms.at(New_Bonds.at(id).atom[j]->atom_ID);
    }
  }

  id = 0;
  for (int i = 0; i < tmp_angles; ++i)
  {
    flag = false;
    for (int TYPE: ATOM_TYPE)
    {
      for (int j = 0; j < atom_num_angle; ++j)
      {
        if (Angles.at(i+1).atom[j]->atom_type == TYPE)
        {
          flag = true;
          --angles;
          break;
        }
      }
      if (flag) break;
    }
    if (flag) continue;
    ++id;
    New_Angles[id] = Angles.at(i+1);
    New_Angles[id].angle_ID = id;
    for (int j = 0; j < atom_num_angle; ++j)
    {
      New_Angles.at(id).atom[j] = &New_Atoms.at(New_Angles.at(id).atom[j]->atom_ID);
    }
  }

  id = 0;
  for (int i = 0; i < tmp_dihedrals; ++i)
  {
    flag = false;
    for (int TYPE: ATOM_TYPE)
    {
      for (int j = 0; j < atom_num_dihedral; ++j)
      {
        if (Dihedrals.at(i+1).atom[j]->atom_type == TYPE)
        {
          flag = true;
          --dihedrals;
          break;
        }
      }
      if (flag) break;
    }
    if (flag) continue;
    ++id;
    New_Dihedrals[id] = Dihedrals.at(i+1);
    New_Dihedrals[id].dihedral_ID = id;
    for (int j = 0; j < atom_num_dihedral; ++j)
    {
      New_Dihedrals.at(id).atom[j] = &New_Atoms.at(New_Dihedrals.at(id).atom[j]->atom_ID);
    }
  }

  id = 0;
  for (int i = 0; i < tmp_impropers; ++i)
  {
    flag = false;
    for (int TYPE: ATOM_TYPE)
    {
      for (int j = 0; j < atom_num_improper; ++j)
      {
        if (Impropers.at(i+1).atom[j]->atom_type == TYPE)
        {
          flag = true;
          --impropers;
          break;
        }
      }
      if (flag) break;
    }
    if (flag) continue;
    ++id;
    New_Impropers[id] = Impropers.at(i+1);
    New_Impropers[id].improper_ID = id;
    for (int j = 0; j < atom_num_improper; ++j)
    {
      New_Impropers.at(id).atom[j] = &New_Atoms.at(New_Impropers.at(id).atom[j]->atom_ID);
    }
  }

  Atoms = New_Atoms;
  Bonds = New_Bonds;
  Angles = New_Angles;
  Dihedrals = New_Dihedrals;
  Impropers = New_Impropers;
  for (int i = 0; i < bonds; ++i)
  {
    for (int j = 0; j < atom_num_bond; ++j)
    {
      Bonds.at(i+1).atom[j] = &Atoms.at(Bonds.at(i+1).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < angles; ++i)
  {
    for (int j = 0; j < atom_num_angle; ++j)
    {
      Angles.at(i+1).atom[j] = &Atoms.at(Angles.at(i+1).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < dihedrals; ++i)
  {
    for (int j = 0; j < atom_num_dihedral; ++j)
    {
      Dihedrals.at(i+1).atom[j] = &Atoms.at(Dihedrals.at(i+1).atom[j]->atom_ID);
    }
  }
  for (int i = 0; i < impropers; ++i)
  {
    for (int j = 0; j < atom_num_improper; ++j)
    {
      Impropers.at(i+1).atom[j] = &Atoms.at(Impropers.at(i+1).atom[j]->atom_ID);
    }
  }
  check_bonds();
}

// Others
bool Data::is_included(const string& S, const string& s)
{
  if (S.find(s) != std::string::npos) return true;
  else return false;
}

string Data::slice(const string& str, int l, int r)
{
  if (r < 0) r = str.size() + r + 1;
  if (l < 0)
  {
    cerr << endl << "Data::slice() : l is smaller than 0" << endl;
    exit(0);
  }
  if (r > int(str.size()))
  {
    cerr << endl << "Data::slice() : r is larger than size" << endl;
    exit(0);
  }
  if (l > r)
  {
    cerr << endl << "Data::slice() : l is larger than r" << endl;
    exit(0);
  }
  string ret;
  for (int i = 0; i < int(str.size()); ++i) if (l <= i && i < r) ret += str[i];

  return ret;
}

string Data::shift_right(string str, const int& n)
{
  int position = 0;
  for (char c: str)
  {
    if (c == '.') break;
    else ++position;
  }
  for (int i = 0; i < n; ++i)
  {
    if (position == int(str.size())-1)
    {
      str.pop_back();
      str += "0.";
    }
    else
    {
      str = slice(str, 0, position) + str[position+1] + "." + slice(str, position+2, -1);
    }
    ++position;
  }

  return str;
}

string Data::shift_left(string str, const int& n)
{
  int position = 0;
  for (char c: str)
  {
    if (c == '.') break;
    else ++position;
  }
  for (int i = 0; i < n; ++i)
  {
    if (position == 0)
    {
      str = ".0" + slice(str, 1, -1);
    }
    else
    {
      str = slice(str, 0, position-1) + "." + str[position-1] + slice(str, position+1, -1);
      --position;
    }
  }

  return str;
}

double Data::convert(const string& str)
{
  if (!is_included(str, "e") && !is_included(str, "E")) return stod(str);
  if (!is_included(str, ".")) return stod(str);
  string base = "", exponential = "", value;
  int e, ct = -1;
  bool flag = false, minus = false;
  for (char c: str)
  {
    ++ct;
    if (ct == 0 && c == '-')
    {
      minus = true;
      continue;
    }
    if (flag) exponential += c;
    else if (c == 'e' || c == 'E') flag = true;
    else base += c;
  }
  e = stoi(exponential);
  if (exponential[0] == '+') value = shift_right(base, e);
  else value = shift_left(base, -e);
  if (minus) value = "-" + value;

  return stod(value);
}

void Data::check_bonds()
{
  for (int i = 0; i < atoms; ++i)
  {
    Atoms.at(i+1).have_bond = false;
  }
  for (int i = 0; i < bonds; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {
      Atoms.at(Bonds.at(i+1).atom[j]->atom_ID).have_bond = true;
    }
  }
}