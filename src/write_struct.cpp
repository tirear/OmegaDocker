// OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
// Copyright (C) 2026 Gao Quan-Ze.
// SPDX-License-Identifier: LGPL-2.1-or-later

#include "write_struct.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <set>

#include "log.hpp"

namespace od {

namespace {

const std::set<std::string>& standardResidues() {
  static const std::set<std::string> s{"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                                       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                                       "A",   "C",   "G",   "U",   "I",   "DA",  "DC",  "DG",  "DT",  "DU",
                                       "DI"};
  return s;
}

// Residue fields of an atom; molecules without residue information get one
// residue UNL in chain L (MobM) or S (StaM), sequence number 1.
struct ResFields {
  std::string chain, name, iCode;
  int seq = 1;
  bool standard = false;
};

ResFields resFields(const Molecule& m, int a) {
  const int r = m.atoms[a].residue;
  if (!m.hasResidueInfo || r < 0) return {m.role == Role::MobM ? "L" : "S", "UNL", "", 1, false};
  const ResidueInfo& i = m.residues[r].info;
  return {i.chain, i.name, i.iCode, i.seq, standardResidues().count(i.name) > 0};
}

// Atom names; atoms without names: element symbol + running number per element.
std::vector<std::string> atomNames(const Molecule& m) {
  std::vector<std::string> n(m.atoms.size());
  std::array<int, 128> count{};
  for (std::size_t a = 0; a < m.atoms.size(); ++a) {
    if (!m.atoms[a].name.empty()) {
      n[a] = m.atoms[a].name;
      continue;
    }
    const int z = std::clamp(m.atoms[a].atomicNum, 0, 127);
    n[a] = std::format("{}{}", elementSymbol(z), ++count[z]);
  }
  return n;
}

std::string upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

// PDB columns 13-16: names of 1-letter elements start in column 14.
std::string pdbAtomName(const std::string& name, int atomicNum) {
  if (name.size() >= 4) return name.substr(0, 4);
  if (std::string(elementSymbol(atomicNum)).size() == 1) return std::format(" {:<3}", name);
  return std::format("{:<4}", name);
}

std::string chargeField(int q) {
  if (q == 0) return "  ";
  return std::format("{}{}", std::abs(q), q > 0 ? '+' : '-');
}

std::string cifValue(const std::string& s) {
  if (s.empty()) return "?";
  bool quote = s.find_first_of(" \t\"'") != std::string::npos || std::string("_#$'\"[];").find(s[0]) != std::string::npos;
  if (!quote) return s;
  return s.find('\'') == std::string::npos ? "'" + s + "'" : "\"" + s + "\"";
}

std::string titleLine(const std::string& t) { return t.empty() ? "OmegaDocker" : t; }

void writeFile(const std::string& path, const std::string& text) {
  std::ofstream f(path, std::ios::binary);
  if (!f) fail("cannot write {}", path);
  f << text;
  if (!f) fail("cannot write {}", path);
}

std::vector<std::vector<std::pair<int, const Bond*>>> adjacency(const Molecule& m) {
  std::vector<std::vector<std::pair<int, const Bond*>>> adj(m.atoms.size());
  for (const Bond& b : m.bonds) {
    adj[b.a].push_back({b.b, &b});
    adj[b.b].push_back({b.a, &b});
  }
  return adj;
}

// A bond between two residues other than a standard backbone link (peptide C-N, nucleic O3'-P):
// written as CONECT (PDB) and _struct_conn (mmCIF).
bool extraLink(const Molecule& m, const Bond& b) {
  const int ra = m.atoms[b.a].residue, rb = m.atoms[b.b].residue;
  if (!m.hasResidueInfo || ra < 0 || rb < 0 || ra == rb) return false;
  const auto& sr = standardResidues();
  if (sr.count(m.residues[ra].info.name) && sr.count(m.residues[rb].info.name)) {
    const std::string& x = m.atoms[b.a].name;
    const std::string& y = m.atoms[b.b].name;
    if ((x == "C" && y == "N") || (x == "N" && y == "C") || (x == "O3'" && y == "P") || (x == "P" && y == "O3'"))
      return false;
  }
  return true;
}

}  // namespace

bool fitsPdb(std::span<const StructPart> parts) {
  std::size_t n = 0;
  for (const auto& p : parts) {
    n += p.mol->atoms.size();
    for (std::size_t a = 0; a < p.mol->atoms.size(); ++a) {
      const ResFields r = resFields(*p.mol, static_cast<int>(a));
      if (r.chain.size() > 1 || r.name.size() > 3 || r.iCode.size() > 1 || r.seq > 9999 || r.seq < -999) return false;
    }
  }
  return n <= 99999;
}

std::string pdbText(std::span<const StructPart> parts, const std::string& title) {
  std::string s = std::format("REMARK   1 {}\n", titleLine(title));
  int serial = 0;
  std::string conect;
  for (const auto& p : parts) {
    const Molecule& m = *p.mol;
    const auto names = atomNames(m);
    std::vector<int> ser(m.atoms.size());
    std::vector<char> het(m.atoms.size());
    std::string lastChain;
    ResFields last;
    for (std::size_t a = 0; a < m.atoms.size(); ++a) {
      const Atom& at = m.atoms[a];
      const ResFields r = resFields(m, static_cast<int>(a));
      if (a > 0 && r.chain != lastChain) s += std::format("TER   {:>5}      {:>3} {:1}{:>4}{:1}\n", ++serial, last.name, last.chain, last.seq, last.iCode);
      lastChain = r.chain;
      last = r;
      ser[a] = ++serial;
      het[a] = !r.standard;
      const auto& x = p.xyz[a];
      s += std::format("{:<6}{:>5} {:4}{:1}{:>3} {:1}{:>4}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2}{:2}\n",
                       het[a] ? "HETATM" : "ATOM", ser[a], pdbAtomName(names[a], at.atomicNum), "", r.name, r.chain,
                       r.seq, r.iCode, x[0], x[1], x[2], 1.0, 0.0, upper(elementSymbol(at.atomicNum)),
                       chargeField(at.formalCharge));
    }
    if (!m.atoms.empty())
      s += std::format("TER   {:>5}      {:>3} {:1}{:>4}{:1}\n", ++serial, last.name, last.chain, last.seq, last.iCode);
    // CONECT for bonds with a HETATM atom and for links between residues other than the backbone.
    std::vector<std::vector<int>> nb(m.atoms.size());
    for (const Bond& b : m.bonds)
      if (het[b.a] || het[b.b] || extraLink(m, b)) {
        nb[b.a].push_back(ser[b.b]);
        nb[b.b].push_back(ser[b.a]);
      }
    for (std::size_t a = 0; a < m.atoms.size(); ++a)
      for (std::size_t k = 0; k < nb[a].size(); k += 4) {
        conect += std::format("CONECT{:>5}", ser[a]);
        for (std::size_t j = k; j < nb[a].size() && j < k + 4; ++j) conect += std::format("{:>5}", nb[a][j]);
        conect += "\n";
      }
  }
  return s + conect + "END\n";
}

std::string cifText(std::span<const StructPart> parts, const std::string& title) {
  std::string block = titleLine(title);
  for (char& c : block)
    if (std::isspace(static_cast<unsigned char>(c))) c = '_';
  std::string s = std::format("data_{}\n#\nloop_\n", block);
  for (const char* f : {"group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
                        "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y",
                        "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge", "auth_seq_id", "auth_comp_id",
                        "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"})
    s += std::format("_atom_site.{}\n", f);
  int serial = 0;
  for (const auto& p : parts) {
    const Molecule& m = *p.mol;
    const auto names = atomNames(m);
    for (std::size_t a = 0; a < m.atoms.size(); ++a) {
      const Atom& at = m.atoms[a];
      const ResFields r = resFields(m, static_cast<int>(a));
      const auto& x = p.xyz[a];
      s += std::format("{} {} {} {} . {} {} . {} {} {:.3f} {:.3f} {:.3f} 1.00 0.00 {} {} {} {} {} 1\n",
                       r.standard ? "ATOM" : "HETATM", ++serial, elementSymbol(at.atomicNum), cifValue(names[a]),
                       cifValue(r.name), cifValue(r.chain), r.standard ? std::to_string(r.seq) : std::string("."),
                       r.iCode.empty() ? "?" : cifValue(r.iCode), x[0], x[1], x[2], at.formalCharge, r.seq,
                       cifValue(r.name), cifValue(r.chain), cifValue(names[a]));
    }
  }
  // Bonds within residues that are not standard residues (e.g. UNL without CCD entry), one entry per residue name.
  std::string bonds;
  std::set<std::string> done;
  for (const auto& p : parts) {
    const Molecule& m = *p.mol;
    const auto names = atomNames(m);
    for (const Bond& b : m.bonds) {
      const ResFields ra = resFields(m, b.a), rb = resFields(m, b.b);
      if (ra.standard || ra.name != rb.name || ra.chain != rb.chain || ra.seq != rb.seq || ra.iCode != rb.iCode) continue;
      if (m.hasResidueInfo && m.atoms[b.a].residue != m.atoms[b.b].residue) continue;
      const std::string key = ra.name + " " + names[b.a] + " " + names[b.b];
      if (!done.insert(key).second) continue;
      const char* order = b.order == BondOrder::Double ? "DOUB" : b.order == BondOrder::Triple ? "TRIP" : "SING";
      bonds += std::format("{} {} {} {} {}\n", cifValue(ra.name), cifValue(names[b.a]), cifValue(names[b.b]), order,
                           b.aromatic ? "Y" : "N");
    }
  }
  std::string conn;
  int nConn = 0;
  for (const auto& p : parts) {
    const Molecule& m = *p.mol;
    const auto names = atomNames(m);
    for (const Bond& b : m.bonds) {
      if (!extraLink(m, b)) continue;
      const bool metal = isMetalBond(m, b);
      const bool ss = m.atoms[b.a].atomicNum == 16 && m.atoms[b.b].atomicNum == 16;
      const char* type = metal ? "metalc" : ss ? "disulf" : "covale";
      conn += std::format("{}{} {}", type, ++nConn, type);
      for (int a : {b.a, b.b}) {
        const ResFields r = resFields(m, a);
        conn += std::format(" {} {} {} {} {} {} {} {}", cifValue(r.chain), cifValue(r.name),
                            r.standard ? std::to_string(r.seq) : std::string("."), cifValue(names[a]),
                            cifValue(r.chain), r.seq, cifValue(r.name), r.iCode.empty() ? "?" : cifValue(r.iCode));
      }
      conn += "\n";
    }
  }
  if (!conn.empty()) {
    s += "#\nloop_\n_struct_conn.id\n_struct_conn.conn_type_id\n";
    for (int k = 1; k <= 2; ++k)
      for (const char* f : {"label_asym_id", "label_comp_id", "label_seq_id", "label_atom_id", "auth_asym_id",
                            "auth_seq_id", "auth_comp_id", "PDB_ins_code"})
        s += std::string(std::string(f) == "PDB_ins_code" ? std::format("_struct_conn.pdbx_ptnr{}_", k)
                                                          : std::format("_struct_conn.ptnr{}_", k)) + f + "\n";
    s += conn;
  }
  if (!bonds.empty())
    s += "#\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n"
         "_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n" + bonds;
  return s + "#\n";
}

std::string sdfText(std::span<const StructPart> parts, const std::string& title) {
  std::string s;
  for (const auto& p : parts) {
    const Molecule& m = *p.mol;
    const std::size_t na = m.atoms.size(), nb = m.bonds.size();
    auto order = [](const Bond& b) { return b.order == BondOrder::Metal ? 1 : static_cast<int>(b.order); };
    s += std::format("{}\n  OmegaDocker 3D\n{}\n", m.id, titleLine(title));
    if (na <= 999 && nb <= 999) {
      s += std::format("{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n", na, nb);
      for (std::size_t a = 0; a < na; ++a) {
        const auto& x = p.xyz[a];
        s += std::format("{:10.4f}{:10.4f}{:10.4f} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n", x[0], x[1], x[2],
                         elementSymbol(m.atoms[a].atomicNum));
      }
      for (const Bond& b : m.bonds) s += std::format("{:>3}{:>3}{:>3}  0\n", b.a + 1, b.b + 1, order(b));
      std::vector<std::size_t> charged;
      for (std::size_t a = 0; a < na; ++a)
        if (m.atoms[a].formalCharge != 0) charged.push_back(a);
      for (std::size_t k = 0; k < charged.size(); k += 8) {
        const std::size_t n = std::min<std::size_t>(8, charged.size() - k);
        s += std::format("M  CHG{:>3}", n);
        for (std::size_t j = k; j < k + n; ++j) s += std::format(" {:>3} {:>3}", charged[j] + 1, m.atoms[charged[j]].formalCharge);
        s += "\n";
      }
      s += "M  END\n";
    } else {
      s += "  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\n";
      s += std::format("M  V30 COUNTS {} {} 0 0 0\nM  V30 BEGIN ATOM\n", na, nb);
      for (std::size_t a = 0; a < na; ++a) {
        const auto& x = p.xyz[a];
        s += std::format("M  V30 {} {} {:.4f} {:.4f} {:.4f} 0", a + 1, elementSymbol(m.atoms[a].atomicNum), x[0], x[1], x[2]);
        if (m.atoms[a].formalCharge != 0) s += std::format(" CHG={}", m.atoms[a].formalCharge);
        s += "\n";
      }
      s += "M  V30 END ATOM\nM  V30 BEGIN BOND\n";
      for (std::size_t k = 0; k < nb; ++k)
        s += std::format("M  V30 {} {} {} {}\n", k + 1, order(m.bonds[k]), m.bonds[k].a + 1, m.bonds[k].b + 1);
      s += "M  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    }
    s += "> <atom.dprop.PartialCharge>\n";
    for (std::size_t a = 0; a < na; ++a) s += std::format("{}{:.4f}", a ? " " : "", m.atoms[a].partialCharge);
    s += "\n\n$$$$\n";
  }
  return s;
}

namespace {

using Adjacency = std::vector<std::vector<std::pair<int, const Bond*>>>;

std::string sybylTypeAdj(const Molecule& m, const Adjacency& adj, int a) {
  const Atom& at = m.atoms[a];
  const int z = at.atomicNum;
  int doubles = 0, triples = 0;
  bool aromatic = at.aromatic;
  for (const auto& [o, b] : adj[a]) {
    if (b->order == BondOrder::Double) ++doubles;
    if (b->order == BondOrder::Triple) ++triples;
    aromatic |= b->aromatic;
  }
  const int deg = static_cast<int>(adj[a].size());
  auto terminalO = [&](int c) {  // terminal O neighbours of atom c
    int n = 0;
    for (const auto& [o, b] : adj[c])
      if (m.atoms[o].atomicNum == 8 && adj[o].size() == 1) ++n;
    return n;
  };
  auto hasDoubleTo = [&](int c, int elem) {
    for (const auto& [o, b] : adj[c])
      if (b->order == BondOrder::Double && m.atoms[o].atomicNum == elem) return true;
    return false;
  };
  switch (z) {
    case 1: return "H";
    case 6: {
      if (aromatic) return "C.ar";
      if (triples > 0 || doubles >= 2) return "C.1";
      int nN = 0;
      bool dN = false;
      for (const auto& [o, b] : adj[a])
        if (m.atoms[o].atomicNum == 7) {
          ++nN;
          dN |= b->order == BondOrder::Double;
        }
      if (nN == 3 && dN) return "C.cat";
      return doubles > 0 ? "C.2" : "C.3";
    }
    case 7: {
      if (aromatic) return deg >= 3 && at.formalCharge == 0 ? "N.pl3" : "N.ar";  // pyrrole-type N: N.pl3
      if (triples > 0) return "N.1";
      if (doubles == 0 && deg == 4) return "N.4";
      if (doubles > 0) return "N.2";
      for (const auto& [o, b] : adj[a])
        if (m.atoms[o].atomicNum == 6 && (hasDoubleTo(o, 8) || hasDoubleTo(o, 16))) return "N.am";
      for (const auto& [o, b] : adj[a]) {
        const int oz = m.atoms[o].atomicNum;
        bool sp2 = m.atoms[o].aromatic;
        for (const auto& [p, c] : adj[o]) sp2 |= c->order == BondOrder::Double || c->aromatic;
        if (sp2 && (oz == 6 || oz == 7 || oz == 8 || oz == 16)) return "N.pl3";
      }
      return "N.3";
    }
    case 8: {
      if (deg == 1) {
        const int c = adj[a][0].first;
        const int cz = m.atoms[c].atomicNum;
        if ((cz == 6 || cz == 15) && terminalO(c) >= 2) return "O.co2";
      }
      return doubles > 0 ? "O.2" : "O.3";
    }
    case 16: {
      const int o = terminalO(a);
      if (o >= 2) return "S.O2";
      if (o == 1 && deg >= 3) return "S.O";
      return doubles > 0 ? "S.2" : "S.3";
    }
    case 15: return "P.3";
    default: return elementSymbol(z);
  }
}

}  // namespace

std::string sybylType(const Molecule& m, int a) { return sybylTypeAdj(m, adjacency(m), a); }

std::string mol2Text(std::span<const StructPart> parts, const std::string& title) {
  std::string s;
  for (const auto& p : parts) {
    const Molecule& m = *p.mol;
    const auto names = atomNames(m);
    const Adjacency adj = adjacency(m);
    std::vector<std::string> types(m.atoms.size());
    for (std::size_t a = 0; a < m.atoms.size(); ++a) types[a] = sybylTypeAdj(m, adj, static_cast<int>(a));
    // Substructures: one per residue in order of first appearance.
    std::vector<int> subst(m.atoms.size());
    std::vector<std::pair<std::string, ResFields>> subs;
    std::vector<int> root;
    for (std::size_t a = 0; a < m.atoms.size(); ++a) {
      const ResFields r = resFields(m, static_cast<int>(a));
      const std::string key = std::format("{}{}{}", r.name, r.seq, r.iCode);
      if (subs.empty() || subs.back().first != key || subs.back().second.chain != r.chain) {
        subs.push_back({key, r});
        root.push_back(static_cast<int>(a) + 1);
      }
      subst[a] = static_cast<int>(subs.size());
    }
    s += std::format("# {}\n@<TRIPOS>MOLECULE\n{}\n{} {} {} 0 0\n{}\n{}\n\n", titleLine(title), m.id, m.atoms.size(),
                     m.bonds.size(), subs.size(), m.hasResidueInfo ? "BIOPOLYMER" : "SMALL",
                     m.chargesFromFile ? "USER_CHARGES" : "GASTEIGER");
    s += "@<TRIPOS>ATOM\n";
    for (std::size_t a = 0; a < m.atoms.size(); ++a) {
      const auto& x = p.xyz[a];
      s += std::format("{:>7} {:<6} {:10.4f} {:10.4f} {:10.4f} {:<6} {:>4} {:<8} {:9.4f}\n", a + 1, names[a], x[0], x[1],
                       x[2], types[a], subst[a], subs[subst[a] - 1].first,
                       m.atoms[a].partialCharge);
    }
    s += "@<TRIPOS>BOND\n";
    for (std::size_t k = 0; k < m.bonds.size(); ++k) {
      const Bond& b = m.bonds[k];
      std::string t;
      if (b.aromatic) t = "ar";
      else if (b.order == BondOrder::Metal) t = "1";
      else {
        t = std::to_string(static_cast<int>(b.order));
        if (b.order == BondOrder::Single) {
          const std::string& ta = types[b.a];
          const std::string& tb = types[b.b];
          if ((ta == "N.am" && tb == "C.2") || (ta == "C.2" && tb == "N.am")) t = "am";
        }
      }
      s += std::format("{:>6} {:>5} {:>5} {}\n", k + 1, b.a + 1, b.b + 1, t);
    }
    s += "@<TRIPOS>SUBSTRUCTURE\n";
    for (std::size_t k = 0; k < subs.size(); ++k)
      s += std::format("{:>6} {:<8} {:>6} {} 1 {} {}\n", k + 1, subs[k].first, root[k],
                       m.hasResidueInfo ? "RESIDUE" : "GROUP", subs[k].second.chain.empty() ? "****" : subs[k].second.chain,
                       subs[k].second.name);
  }
  return s;
}

std::vector<std::string> writeStructure(const std::string& base, std::span<const StructPart> parts,
                                        const OutFormats& f, const std::string& title, bool& pdbFallback) {
  std::vector<std::string> files;
  const std::string name = std::filesystem::path(base).filename().string();
  pdbFallback = f.pdb && !fitsPdb(parts);
  if (f.pdb && !pdbFallback) {
    writeFile(base + ".pdb", pdbText(parts, title));
    files.push_back(name + ".pdb");
  }
  if (f.cif || pdbFallback) {
    writeFile(base + ".cif", cifText(parts, title));
    files.push_back(name + ".cif");
  }
  if (f.sdf) {
    writeFile(base + ".sdf", sdfText(parts, title));
    files.push_back(name + ".sdf");
  }
  if (f.mol2) {
    writeFile(base + ".mol2", mol2Text(parts, title));
    files.push_back(name + ".mol2");
  }
  return files;
}

}  // namespace od
