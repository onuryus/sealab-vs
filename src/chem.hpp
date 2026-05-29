#pragma once
#include <memory>
#include <string>

#include "options.hpp"

// Forward declarations — keep RDKit out of the public header (pimpl).
namespace RDKit {
class ROMol;
}
class ExplicitBitVect;

// Per-molecule result. One Record = one CSV row.
struct Record {
  bool valid = false;     // true → passed all gates and is written to CSV
  long line_no = 0;
  std::string smiles;     // canonical SMILES after standardization

  // Physicochemical descriptors
  double mw = 0, logp = 0, mr = 0, tpsa = 0, tpsa_sp = 0, fsp3 = 0, bbb_mpo = 0;
  int hbd = 0, hba = 0, rotb = 0, arom = 0, heavy = 0, natoms = 0, rings = 0,
      carbons = 0, hetero = 0, charge = 0, stereo = 0;

  // QED score (0–1; Bickerton 2012)
  double qed = 0;

  // Rule-set pass/fail flags
  bool lipinski = false, veber = false, ghose = false, egan = false,
       muegge = false, ro3 = false, leadlike = false, reos = false,
       gsk = false, pfizer = false;

  // Alert catalogs (true = CLEAN, i.e. no substructure match)
  bool pains_free = false, brenk_free = false;
  bool lilly_free = true, inpharmatica_free = true;
  bool ames_free = true, aggregators_free = true;
  bool extra_free = true;          // NIH+ZINC+CHEMBL bundled (legacy "extra")

  // BOILED-Egg permeability regions
  bool gia = false, bbb = false;

  // Similarity to target (the screening output)
  double similarity = 0;
};

// Shared, READ-ONLY filter catalogs. Built once, all threads use them via
// const methods (FilterCatalog::hasMatch is thread-safe). Avoids paying
// catalog-load cost per thread.
class SharedCatalogs {
 public:
  explicit SharedCatalogs(bool include_extra);
  ~SharedCatalogs();

  bool has_pains(const RDKit::ROMol& m) const;
  bool has_brenk(const RDKit::ROMol& m) const;
  bool has_lilly(const RDKit::ROMol& m) const;           // CHEMBL_LINT
  bool has_inpharmatica(const RDKit::ROMol& m) const;    // CHEMBL_Inpharmatica
  bool has_ames(const RDKit::ROMol& m) const;            // Benigni-Bossa SMARTS
  bool has_aggregators(const RDKit::ROMol& m) const;     // Shoichet-style SMARTS
  bool has_extra(const RDKit::ROMol& m) const;           // NIH+ZINC+CHEMBL bundled
  bool extra_enabled() const;

  // Number of PAINS + BRENK substructure matches (used as the ALERTS
  // input for QED, as a Bickerton-list-free proxy).
  int alert_count(const RDKit::ROMol& m) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// Per-thread chemistry context. Holds only mutable, non-thread-safe
// objects (RDKit standardizers). Borrows a pointer to the shared catalog set.
class ChemContext {
 public:
  ChemContext(const Options& o, const SharedCatalogs& cats);
  ~ChemContext();
  ChemContext(ChemContext&&) noexcept;

  // Run the full screening pipeline on one .smi line.
  // valid=false → parse failure or filtered out.
  Record process(const std::string& raw, long line_no, const Options& o,
                 const ExplicitBitVect* target_fp) const;

  // Stage-1 fast path: parse + fingerprint + similarity only.
  // No standardization, no descriptors, no rules, no alerts, no egg.
  // Used for the prefilter step of two-stage screening.
  Record process_prefilter(const std::string& raw, long line_no, const Options& o,
                           const ExplicitBitVect* target_fp) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  const SharedCatalogs* cats_;   // non-owning
};

// Free helpers used by main.cpp.
ExplicitBitVect* target_fingerprint(const std::string& smiles,
                                    const std::string& method,
                                    int bits, int radius, std::string& err);
void disable_rdkit_logs();
