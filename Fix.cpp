#include "/mnt/c/Ubuntu/home/kaito/lammps/header/LAMMPS.h"
#include <iostream>
#include <unordered_map>

#define DEBUG false

template<class MAP, class VALUE> static bool contains_key(MAP m, VALUE v){ return m.find(v) != m.end();}

template<class T>
struct UnionFind
{
  vector<T> parent;
  UnionFind(const T& n): parent(n, -1){}

  // xの親を検索
  T find(const T& x)
  {
    if (parent[x] < 0) return x;
    return parent[x] = find(parent[x]);
  }

  // 同じ集合に属するかどうかをチェック
  bool same_check(const T& x, const T& y)
  {
    return find(x) == find(y);
  }

  // 自分の属する集合のサイズを返す
  T size(const T& x)
  {
    T y = find(x);
    return -parent[y];
  }

  // xとyの属する集合を併合
  bool unite(const T& x, const T& y)
  {
    T x_root = find(x), y_root = find(y);
    if (x_root == y_root) return false;
    else
    {
      if (-parent[x_root] < -parent[y_root]) swap(x_root, y_root);
      parent[x_root] += parent[y_root];
      parent[y_root] = x_root;
      return true;
    }
  }
};

void Fix(Data& Data)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  int molecule_ID = 0;
  unordered_map<int, int> tomolecule_ID;
  UnionFind<int> UF(Data.atoms);
  for (auto x: Data.Bonds)
  {
    int atom_i, atom_j;
    atom_i = x.second.atom[0]->atom_ID-1;
    atom_j = x.second.atom[1]->atom_ID-1;
    UF.unite(atom_i, atom_j);
  }
  for (int i = 0; i < Data.atoms; ++i)
  {
    if (Data.Atoms.at(i+1).molecule_ID <= 0) continue;
    int par = UF.find(i);
    if (par < 0) par = i;
    if (!contains_key(tomolecule_ID, par))
    {
      tomolecule_ID[par] = ++molecule_ID;
    }
    Data.Atoms.at(i+1).molecule_ID = tomolecule_ID.at(par);
  }
}

int main(int argc, char** argv)
{
  Data Data;
  if (argc >= 2)
  {
    Data.Read(argv[1]);
  }
  else
  {
    Data.Read();
  }
  Fix(Data);
  Data.Write("Fixed.data");

  return 0;
}