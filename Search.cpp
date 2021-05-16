#include "/mnt/c/Ubuntu/home/kaito/lammps/header/LAMMPS.h"
#include <iostream>
#include <unordered_map>
#include <queue>
#include <stack>

#define DEBUG false

const int INF = 1e9;

template<class T> using PRIORITY_QUEUE = priority_queue< T, vector<T>, greater<T> >;
template<class T> inline bool chmax(T &a, T b){if (a < b) {a = b; return true;} return false;}
template<class MAP, class VALUE> static bool contains_key(MAP m, VALUE v){ return m.find(v) != m.end();}

int bond_type_mainchain, atom_type_Carbon;
map<pair<int, int>, int> tobondID;

void warshall_floyd(vector<vector<int>>& G)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  /*
  ワ―シャルフロイド法。全点対間の最短距離を求める。計算量は O(V^3)。
  負の閉路がなければ負の辺にも対応可。
  引数のグラフは 0-indexed の隣接行列で与える必要がある。
  G[i][j] : 点iから点jへの辺の重み。辺が存在しない所はINFで初期化しておく。
  最終的にG[i][j]はiからjへの最小コストを格納する。
  中継点をkとしたとき、
  G[i][j] = G[i][k] + G[k][j]
  と表せることから、DPによりすべての点について、その点を中継点としたときのiからjへのパスの組み合わせをすべて計算することで、最終的な最小コストが求まる。
  */
  int V = G.size();
  for (int i = 0; i < V; ++i) G[i][i] = 0;
  for (int k = 0; k < V; ++k) for (int i = 0; i < V; ++i) for (int j = 0; j < V; ++j) G[i][j] = min(G[i][j], G[i][k] + G[k][j]);
}

unordered_map<int, unordered_map<int, vector<int>>> get_graph(const Data& Data)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  unordered_map<int, unordered_map<int, vector<int>>> Graphs;

  for (int idx = 0; idx < Data.bonds; ++idx)
  {
    tobondID[{Data.Bonds.at(idx+1).atom[0]->atom_ID, Data.Bonds.at(idx+1).atom[1]->atom_ID}] = Data.Bonds.at(idx+1).bond_ID;
    tobondID[{Data.Bonds.at(idx+1).atom[1]->atom_ID, Data.Bonds.at(idx+1).atom[0]->atom_ID}] = Data.Bonds.at(idx+1).bond_ID;

    if (Data.Bonds.at(idx+1).bond_type == bond_type_mainchain)
    {
      int molecule_ID;
      pair<int, int> path;
      molecule_ID = Data.Bonds.at(idx+1).atom.at(0)->molecule_ID;
      path = {Data.Bonds.at(idx+1).atom.at(0)->atom_ID, Data.Bonds.at(idx+1).atom.at(1)->atom_ID};
      if (contains_key(Graphs, molecule_ID))
      {
        if (contains_key(Graphs.at(molecule_ID), path.first))
        {
          Graphs.at(molecule_ID).at(path.first).emplace_back(path.second);
        }
        else
        {
          Graphs.at(molecule_ID)[path.first]= {path.second};
        }
        if (contains_key(Graphs.at(molecule_ID), path.second))
        {
          Graphs.at(molecule_ID).at(path.second).emplace_back(path.first);
        }
        else
        {
          Graphs.at(molecule_ID)[path.second]= {path.first};
        }
      }
      else
      {
        Graphs[molecule_ID][path.first] = {path.second};
        Graphs[molecule_ID][path.second] = {path.first};
      }
    }
  }

  return Graphs;
}

stack<int> bfs(const unordered_map<int, vector<int>>& Graph, const int& start, const int& end)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  PRIORITY_QUEUE<int> queue;
  map<int, pair<int, int>> dis;
  bool flag = true;
  stack<int> st;
  queue.push(start);
  dis[start] = {0, start};

  while(!queue.empty() && flag)
  {
    int u = queue.top();
    queue.pop();
    for (auto v: Graph.at(u))
    {
      if (!contains_key(dis, v))
      {
        dis[v] = {dis.at(u).first+1, u};
        if (v == end)
        {
          flag = false;
          break;
        }
        queue.push(v);
      }
    }
  }

  int node = end;
  st.push(node);
  while (node != start)
  {
    node = dis.at(node).second;
    st.push(node);
  }

  return st;
}

stack<int> determine_mainchain(const unordered_map<int, vector<int>>& Graph)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  PRIORITY_QUEUE<int> queue;
  map<int, int> atomIDtoID;
  map<int, int> IDtoatomID;
  int ID = 0;
  vector<vector<int>> NewGraph(Graph.size(), vector<int>(Graph.size(), INF));
  int max_l = 0;

  for (auto x: Graph)
  {
    int u;
    u = x.first;
    queue.push(u);
  }
  while (!queue.empty())
  {
    int atom_ID;
    atom_ID = queue.top();
    queue.pop();
    atomIDtoID[atom_ID] = ID;
    IDtoatomID[ID] = atom_ID;
    ++ID;
  }
  for (auto x: Graph)
  {
    int u;
    u = x.first;
    for (auto v: x.second)
    {
      NewGraph[atomIDtoID.at(u)][atomIDtoID.at(v)] = 1;
    }
  }

  warshall_floyd(NewGraph);
  int start, end;
  for (int i = 0; i < NewGraph.size(); ++i)
  {
    for (int j = 0; j < NewGraph.size(); ++j)
    {
      if (NewGraph[i][j] == INF) continue;
      if (chmax(max_l, NewGraph[i][j]))
      {
        start = IDtoatomID[i];
        end = IDtoatomID[j];
      }
    }
  }

  if (Graph.at(Graph.at(start)[0]).size() > Graph.at(Graph.at(end)[0]).size())
  {
    end = Graph.at(end)[0];
  }
  else
  {
    start = Graph.at(start)[0];
  }

  return bfs(Graph, start, end);
}

void search(Data& Data)
{
  if (DEBUG) cout << "--" << __FUNCTION__ << " Called" << endl;
  unordered_map<int, unordered_map<int, vector<int>>> Graphs;
  stack<int> st;
  Graphs = get_graph(Data);
  for (auto x: Graphs)
  {
    int molecule_ID = x.first;
    st = determine_mainchain(x.second);
    while (st.size() > 2)
    {
      int u, v;
      st.pop();
      u = st.top();
      st.pop();
      v = st.top();
      Data.Bonds.at(tobondID.at({u, v})).bond_type = Data.bond_types+1;
    }
  }
  Bond_Coeff BC;
  BC.bond_type = ++Data.bond_types;
  BC.coeffs = Data.Bond_Coeffs.at(bond_type_mainchain).coeffs;
  Data.Bond_Coeffs[Data.bond_types] = BC;
}

int main(int argc, char** argv)
{
  Data Data;
  if (argc >= 4)
  {
    Data.Read(argv[1]);
    atom_type_Carbon = stoi(argv[2]);
    bond_type_mainchain = stoi(argv[3]);
  }
  else
  {
    Data.Read();
    cout << "The bond type of main chain : ";
    cin >> bond_type_mainchain;
    cout << "The atom type of Carbon : ";
    cin >> atom_type_Carbon;
  }
  search(Data);

  Data.Write("Searched.data");

  return 0;
}