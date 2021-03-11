// #ifdef FILE_NAME_HPP
// #define FILE_NAME_HPP
#pragma once

#include <vector>
#include <map>
#include <fstream>
using namespace std;

struct Atom
{
  // Atom
  int atom_ID, molecule_ID, atom_type, etag, spin, spx, spy, spz, template_atom, template_index;
  double ccN, contact_radius, cs_re, cs_im, cv, density, diameter, espf, edpd_temp, edpd_cv, eradius, kernel_radius, mass, mux, muy, muz, q, rho, sp, theta, volume, x, y, z, x0, y0, z0;
  bool ellipsoidflag, lineflag, triangleflag, nx, ny, nz;

  // Velocity
  double vx, vy, vz, lx, ly, lz, wx, wy, wz, ervel;

  Atom();

  void write_Atom(const string&, ofstream&);
  void output_Atom(const string&);
  void write_Velocity(const string&, ofstream&);
  void output_Velocity(const string&);
};

struct Mass
{
  int atom_type;
  double mass;

  Mass();

  void write(ofstream&);
  void output();
};

struct Bond
{
  int atom_num;
  int bond_ID, bond_type;
  vector<Atom*> atom;

  Bond();

  void write(ofstream&);
  void output();
};

struct Angle
{
  int atom_num;
  int angle_ID, angle_type;
  vector<Atom*> atom;

  Angle();

  void write(ofstream&);
  void output();
};

struct Dihedral
{
  int atom_num;
  int dihedral_ID, dihedral_type;
  vector<Atom*> atom;

  Dihedral();

  void write(ofstream&);
  void output();
};

struct Improper
{
  int atom_num;
  int improper_ID, improper_type;
  vector<Atom*> atom;

  Improper();

  void write(ofstream&);
  void output();
};

struct Pair_Coeff
{
  int atom_type;
  string coeffs;

  Pair_Coeff();

  void write(ofstream&);
  void output();
};

struct Bond_Coeff
{
  int bond_type;
  string coeffs;

  Bond_Coeff();

  void write(ofstream&);
  void output();
};

struct Angle_Coeff
{
  int angle_type;
  string coeffs;

  Angle_Coeff();

  void write(ofstream&);
  void output();
};

struct Dihedral_Coeff
{
  int dihedral_type;
  string coeffs;

  Dihedral_Coeff();

  void write(ofstream&);
  void output();
};

struct Improper_Coeff
{
  int improper_type;
  string coeffs;

  Improper_Coeff();

  void write(ofstream&);
  void output();
};

struct Data
{
  int dimension;
  string atom_style;
  string pair_coeff_style, bond_coeff_style, angle_coeff_style, dihedral_coeff_style, improper_coeff_style;
  string filename, line_header;
  vector<string> data;

  // Header Keywords
  int atoms, bonds, angles, dihedrals, impropers, atom_types, bond_types, angle_types, dihedral_types, improper_types;
  int extra_bond_per_atom, extra_angle_per_atom, extra_dihedral_per_atom, extra_improper_per_atom, extra_special_per_atom;
  int ellipsoids, lines, triangles, bodies;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double xy, xz, yz;

  // Body Keywords
  map<int, Atom> Atoms;
  map<int, Bond> Bonds;
  map<int, Angle> Angles;
  map<int, Dihedral> Dihedrals;
  map<int, Improper> Impropers;
  map<int, Mass> Masses;
  map<int, Pair_Coeff> Pair_Coeffs;
  map<int, Bond_Coeff> Bond_Coeffs;
  map<int, Angle_Coeff> Angle_Coeffs;
  map<int, Dihedral_Coeff> Dihedral_Coeffs;
  map<int, Improper_Coeff> Improper_Coeffs;

  // Costructor
  Data();

  // Initialize
  void initialize_parameters();

  // Get Header
  void get_atoms(const string&);
  void get_atom_types(const string&);
  void get_bonds(const string&);
  void get_bond_types(const string&);
  void get_angles(const string&);
  void get_angle_types(const string&);
  void get_dihedrals(const string&);
  void get_dihedral_types(const string&);
  void get_impropers(const string&);
  void get_improper_types(const string&);
  void get_size_x(const string&);
  void get_size_y(const string&);
  void get_size_z(const string&);

  // Get Body
  void get_Masses(const int&);
  string get_coeffs(const string&);
  void get_Pair_Coeffs(const string&, const int&);
  void get_Bond_Coeffs(const string&, const int&);
  void get_Angle_Coeffs(const string&, const int&);
  void get_Dihedral_Coeffs(const string&, const int&);
  void get_Improper_Coeffs(const string&, const int&);
  void get_Atoms(const string&, const int&);
  void get_Velocities(const int&);
  void get_Bonds(const int&);
  void get_Angles(const int&);
  void get_Dihedrals(const int&);
  void get_Impropers(const int&);

  // Read
  string check_extension(const string&);
  void get_filename(const string&);
  void get_data(const string&);
  void read_data();
  void Read(const string& FILENAME = "");

  // Write
  void Write(string FILENAME = "");
  void write_header(ofstream&);
  void write_body(ofstream&);
  void write_Masses(ofstream&);
  void write_Pair_Coeffs(ofstream&);
  void write_Bond_Coeffs(ofstream&);
  void write_Angle_Coeffs(ofstream&);
  void write_Dihedral_Coeffs(ofstream&);
  void write_Improper_Coeffs(ofstream&);
  void write_Atoms(ofstream&);
  void write_Velocities(ofstream&);
  void write_Bonds(ofstream&);
  void write_Angles(ofstream&);
  void write_Dihedrals(ofstream&);
  void write_Impropers(ofstream&);

  // Output
  void Output();
  void output_header();
  void output_body();
  void output_Masses();
  void output_Pair_Coeffs();
  void output_Bond_Coeffs();
  void output_Angle_Coeffs();
  void output_Dihedral_Coeffs();
  void output_Improper_Coeffs();
  void output_Atoms();
  void output_Velocities();
  void output_Bonds();
  void output_Angles();
  void output_Dihedrals();
  void output_Impropers();

  // Functions
  void Remove_Molecule_ID(const vector<int>&);
  void Remove_Atom_Type(const vector<int>&);

  // Others
  bool is_included(const string&, const string&);
  string slice(const string&, int, int);
  string shift_right(string, const int&);
  string shift_left(string, const int&);
  double convert(const string&);
};

// #endif