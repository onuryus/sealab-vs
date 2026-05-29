#pragma once
#include <set>
#include <string>

// Command-line options. Sensible defaults so a minimal invocation
// (--lib + --target) does something useful out of the box.
struct Options {
  // Required
  std::string lib;                       // --lib    : input .smi file
  std::string target;                    // --target : query SMILES (similarity reference)

  // Output
  std::string out = "results.csv";       // --out

  // Fingerprint
  std::string method = "morgan";         // --method : morgan | aptt | rd
  std::string similarity = "Tanimoto";   // --similarity : Tanimoto|Dice|Sokal|McConnaughey|Kulczynski|Cosine
  int fp_bits = 2048;                    // --fp-bits
  int fp_radius = 2;                     // --fp-radius (Morgan only; others ignore)

  // Drug-likeness rule gates (boolean — all selected rules must pass).
  // Default is empty (no rule applied unless requested).
  std::set<std::string> rules = {"lipinski"};

  // Alert filters (multi-select); each selected catalog must be CLEAN
  // (no substructure match). Options: pains, brenk, lilly, inpharmatica,
  // ames, aggregators, extra (= NIH+ZINC+CHEMBL bundled).
  std::set<std::string> alerts;

  // Continuous-score thresholds
  bool qed_min_set = false;
  double qed_min = 0.0;                  // --qed-min      (Bickerton 2012)
  bool fsp3_min_set = false;
  double fsp3_min = 0.0;                 // --fsp3-min     (Lovering 2009)
  bool stereo_max_set = false;
  int stereo_max = 0;                    // --stereo-max

  // BOILED-Egg permeability region + BBB-MPO threshold
  std::string permeability;              // --permeability : bbb | gia | bbb_and_gia
  bool bbb_threshold_set = false;
  double bbb_threshold = 0.0;            // --bbb-threshold

  // Performance + behavior
  int threads = 0;                       // --threads (0 = all cores)
  bool standardize = true;               // --no-standardize disables full pipeline
  bool tautomer = true;                  // --no-tautomer skips tautomer canonicalization only

  // Two-stage screening for large libraries:
  //   K = 0  → single stage (full pipeline on every molecule)
  //   K > 0  → fast fingerprint-only stage 1 → top-K → full pipeline (stage 2)
  long prefilter_top = 0;                // --prefilter K

  // Back-compat shortcut (equivalent to --alerts extra)
  bool extra_alerts = false;             // --extra-alerts

  // Allowed identifiers (validated by parse_args).
  static const std::set<std::string>& valid_rules();
  static const std::set<std::string>& valid_alerts();
};

bool parse_args(int argc, char** argv, Options& o, std::string& err);
void print_usage(const char* prog);
