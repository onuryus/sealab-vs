#include "chem.hpp"

#include <cmath>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/MolStandardize/Fragment.h>
#include <GraphMol/MolStandardize/Normalize.h>
#include <GraphMol/MolStandardize/Charge.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/AtomPairs.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>
#include <RDGeneral/RDLog.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

#include "egg_data.hpp"

namespace bg = boost::geometry;
using BPoint = bg::model::d2::point_xy<double>;
using BPolygon = bg::model::polygon<BPoint>;
namespace RDD = RDKit::Descriptors;

// ============================================================================
// Anonymous namespace: file-local helpers, constants, and SMARTS lists.
// ============================================================================
namespace {

// ---------------------------------------------------------------------------
// BOILED-Egg polygons (built once at first use; read-only thereafter).
// ---------------------------------------------------------------------------
BPolygon make_polygon(const std::vector<std::pair<double, double>>& pts) {
  BPolygon poly;
  for (const auto& p : pts) bg::append(poly.outer(), BPoint(p.first, p.second));
  bg::correct(poly);  // close + orient
  return poly;
}

const BPolygon& gia_polygon() {
  static const BPolygon p = make_polygon(GIA_COORDS);
  return p;
}
const BPolygon& bbb_polygon() {
  static const BPolygon p = make_polygon(BBB_COORDS);
  return p;
}

// ---------------------------------------------------------------------------
// Ames mutagenicity alerts — subset of Benigni & Bossa 2008.
// 20 SMARTS covering: aromatic amines, nitroaromatics, aziridines, alkyl
// halides, nitrosamines, hydrazines, diazo, epoxides, acrylamides, etc.
// Any match flags a potential mutagenicity concern.
// ---------------------------------------------------------------------------
constexpr const char* const AMES_SMARTS[] = {
    "[NX3;H2,H1;!$(NC=O)]c1ccccc1",                // primary/secondary aromatic amine
    "[N+](=O)([O-])c1ccccc1",                      // nitrobenzene
    "[N+](=O)([O-])[#6]",                          // aliphatic nitro
    "N=N",                                         // azo
    "[NX2]=[NX2+]=[NX1-]",                         // azide
    "N(=O)N",                                      // N-nitroso
    "[NX3](N)c1ccccc1",                            // aryl hydrazine
    "N#N",                                         // diazonium
    "C1OC1",                                       // epoxide
    "C1NC1",                                       // aziridine
    "C1SC1",                                       // thiirane
    "[CX4]([F,Cl,Br,I])[CX4]([F,Cl,Br,I])",        // vicinal dihaloalkane
    "[CH2][F,Cl,Br,I]",                            // alpha-haloalkyl
    "C(=O)C[F,Cl,Br,I]",                           // alpha-halo ketone
    "[CX3](=O)[OX2][NX3]",                         // O-acyl hydroxylamine
    "[NX3]([OH])[#6]",                             // N-hydroxylamine
    "[#6]=[#6][N+](=O)[O-]",                       // alpha,beta-unsaturated nitro
    "C(=O)N(O)[#6]",                               // N-hydroxyamide
    "[CX3]=[NX2][OX2H]",                           // certain oximes
    "OC(=O)C(=O)O",                                // certain oxalates
};

// ---------------------------------------------------------------------------
// Aggregator alerts — small subset of Shoichet/FAF-Drugs aggregator-prone
// scaffolds that frequently produce false positives in HTS assays.
// ---------------------------------------------------------------------------
constexpr const char* const AGGREGATOR_SMARTS[] = {
    "c1ccc(N=Nc2ccccc2)cc1",                       // azobenzene
    "[NX3](c1ccccc1)(c2ccccc2)c3ccccc3",           // triphenylamine
    "c1ccc2c(c1)Sc1ccccc1S2",                      // thianthrene
    "C(=S)N",                                      // generic thioamide
    "[CX3](=O)Oc1ccccc1[OH]",                      // salicylate
    "c1ccc2[nH]c3ccccc3c2c1",                      // carbazole
    "C1=CC=C2C(=C1)C=CC3=CC=CC=C23",               // anthracene-like
    "[CX3](=O)NC1=CC=CC=C1[OH]",                   // hydroxyanilide
    "[#6]=C1C(=O)Nc2ccccc21",                      // oxoindole derivative
    "S(=O)(=O)N=[N,O]",                            // sulfonyl-azo/hydrazide
};

// Compile a SMARTS array to ROMol patterns (alive for the program lifetime).
std::vector<std::unique_ptr<RDKit::ROMol>> compile_smarts(
    const char* const* arr, size_t n) {
  std::vector<std::unique_ptr<RDKit::ROMol>> out;
  out.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    try {
      RDKit::ROMol* m = RDKit::SmartsToMol(arr[i]);
      if (m) out.emplace_back(m);
    } catch (...) {
      // silently skip uncompilable patterns
    }
  }
  return out;
}

bool any_substruct_match(
    const RDKit::ROMol& mol,
    const std::vector<std::unique_ptr<RDKit::ROMol>>& patterns) {
  RDKit::MatchVectType m;
  for (const auto& p : patterns)
    if (RDKit::SubstructMatch(mol, *p, m)) return true;
  return false;
}

// ---------------------------------------------------------------------------
// QED — Quantitative Estimate of Drug-likeness (Bickerton et al. 2012).
//
// Eight properties (MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS) are each
// mapped to [0,1] via an asymmetric double sigmoid (ADS) "desirability"
// function, then combined as a weighted geometric mean.
//
// The ALERTS input uses the count of PAINS + BRENK substructure matches as
// a proxy for the Bickerton-specific alert list (which is not bundled here).
// ---------------------------------------------------------------------------
struct ADSCoefs {
  double a, b, c, d, e, f, dmax;
};

// Bickerton 2012, Table S1
constexpr ADSCoefs Q_MW    = {2.817065973,  392.5754953, 290.7489764, 2.419764353,
                              49.22325677,  65.37051707, 104.9805561};
constexpr ADSCoefs Q_ALOGP = {3.172690585,  137.8624751,   2.534937431, 4.581497897,
                               0.822739154,   0.576295591, 131.3186604};
constexpr ADSCoefs Q_HBA   = {2.948620388,  160.4605972,   3.615294657, 4.435986202,
                               0.290141953,   1.300669958, 148.7763046};
constexpr ADSCoefs Q_HBD   = {1.618662227, 1010.051101,    0.985094388, 1e-9,
                               0.713820843,   0.920922555, 258.1632616};
constexpr ADSCoefs Q_PSA   = {1.876861559,  125.2232657,  62.90773554, 87.83366614,
                              12.01999824,  28.51324732, 104.5686167};
constexpr ADSCoefs Q_ROTB  = {0.01,         272.4121427,   2.558379970, 1.565547684,
                               1.271567166,   2.758063707, 105.4420403};
constexpr ADSCoefs Q_AROM  = {3.217788970,  957.7374108,   2.274627939, 1e-9,
                               1.317690384,   0.375760881, 312.3372610};
constexpr ADSCoefs Q_ALERT = {0.01,        1199.094025,   -0.09002883,  0.185904477,
                               0.875193782,   0.417215146, 417.7253140};

// QED_w (the "mean" weighting from Bickerton et al.)
constexpr double W_MW = 0.66, W_ALOGP = 0.46, W_HBA = 0.05, W_HBD = 0.61;
constexpr double W_PSA = 0.06, W_ROTB = 0.65, W_AROM = 0.48, W_ALERT = 0.95;

double ads(double x, const ADSCoefs& c) {
  const double left  = 1.0 + std::exp(-(x - c.c + c.d / 2.0) / c.e);
  const double right = 1.0 + std::exp(-(x - c.c - c.d / 2.0) / c.f);
  return ((c.a + c.b / left) * (1.0 - 1.0 / right)) / c.dmax;
}

double compute_qed(double mw, double alogp, int hba, int hbd, double psa,
                   int rotb, int arom, int alerts) {
  auto d = [](double x, const ADSCoefs& c) {
    return std::max(ads(x, c), 1e-12);
  };
  const double sum_w = W_MW + W_ALOGP + W_HBA + W_HBD + W_PSA + W_ROTB +
                       W_AROM + W_ALERT;
  const double log_sum =
      W_MW    * std::log(d(mw,                          Q_MW))    +
      W_ALOGP * std::log(d(alogp,                       Q_ALOGP)) +
      W_HBA   * std::log(d(static_cast<double>(hba),    Q_HBA))   +
      W_HBD   * std::log(d(static_cast<double>(hbd),    Q_HBD))   +
      W_PSA   * std::log(d(psa,                         Q_PSA))   +
      W_ROTB  * std::log(d(static_cast<double>(rotb),   Q_ROTB))  +
      W_AROM  * std::log(d(static_cast<double>(arom),   Q_AROM))  +
      W_ALERT * std::log(d(static_cast<double>(alerts), Q_ALERT));
  return std::exp(log_sum / sum_w);
}

// ---------------------------------------------------------------------------
// Stereocenter count: assigned chiral tags + unspecified potential centers.
// ---------------------------------------------------------------------------
int count_stereocenters(const RDKit::ROMol& mol) {
  using namespace RDKit;
  int n = 0;
  for (const auto a : mol.atoms()) {
    if (a->getChiralTag() != Atom::CHI_UNSPECIFIED) ++n;
  }
  std::vector<Chirality::StereoInfo> info =
      Chirality::findPotentialStereo(const_cast<ROMol&>(mol));
  for (const auto& si : info) {
    if (si.type == Chirality::StereoType::Atom_Tetrahedral &&
        si.specified == Chirality::StereoSpecified::Unspecified) {
      ++n;
    }
  }
  return n;
}

std::string trim(const std::string& s) {
  const char* ws = " \t\r\n'()";
  const auto b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
}

// ---------------------------------------------------------------------------
// BBB-MPO-like score (Wager 2010-style multi-parameter optimization).
// Direct port of the original Python implementation; predicts blood-brain
// barrier permeation likelihood from descriptors.
// ---------------------------------------------------------------------------
double bbb_mpo_score(const Record& r) {
  const double h = r.heavy;
  const double mwhbn = (r.mw > 0) ? (r.hba + r.hbd) / std::sqrt(r.mw) : 0.0;
  const double ar = r.arom;
  const double tpsa = r.tpsa;
  const double pKa = 7.4 - 0.1 * r.logp + 0.05 * tpsa;  // rough pKa estimate

  static const double p_arom_lut[5] = {0.336367, 0.816016, 1.0, 0.691115, 0.199399};
  const double p_arom = (ar >= 0 && ar <= 4) ? p_arom_lut[static_cast<int>(ar)] : 0.0;
  const double p_heavy =
      (h > 5 && h <= 45)
          ? (0.0000443 * h * h * h - 0.004556 * h * h + 0.12775 * h - 0.463) / 0.624231
          : 0.0;
  const double p_mwhbn =
      (mwhbn > 0.05 && mwhbn <= 0.45)
          ? (26.733 * mwhbn * mwhbn * mwhbn - 31.495 * mwhbn * mwhbn +
             9.5202 * mwhbn - 0.1358) / 0.72258
          : 0.0;
  const double p_tpsa =
      (tpsa > 0 && tpsa <= 120) ? (-0.0067 * tpsa + 0.9598) / 0.9598 : 0.0;
  const double p_pka =
      (pKa > 3 && pKa <= 11)
          ? (0.00045068 * std::pow(pKa, 4) - 0.016331 * std::pow(pKa, 3) +
             0.18618 * pKa * pKa - 0.71043 * pKa + 0.8579) / 0.597488
          : 0.0;
  return p_arom + p_heavy + 1.5 * p_mwhbn + 2.0 * p_tpsa + 0.5 * p_pka;
}

// ---------------------------------------------------------------------------
// Descriptor + rule computation
// ---------------------------------------------------------------------------
void compute_descriptors(const RDKit::ROMol& mol, Record& r) {
  r.mw = RDD::calcExactMW(mol);
  RDD::calcCrippenDescriptors(mol, r.logp, r.mr);   // WLogP + MR
  r.tpsa = RDD::calcTPSA(mol);                       // default (rules + MPO)
  r.tpsa_sp = RDD::calcTPSA(mol, false, true);       // S&P-augmented (BOILED-Egg)
  r.hbd = RDD::calcNumHBD(mol);
  r.hba = RDD::calcNumHBA(mol);
  r.rotb = RDD::calcNumRotatableBonds(mol);
  r.arom = RDD::calcNumAromaticRings(mol);
  r.heavy = RDD::calcNumHeavyAtoms(mol);
  r.fsp3 = RDD::calcFractionCSP3(mol);
  r.charge = RDKit::MolOps::getFormalCharge(mol);
  r.rings = mol.getRingInfo()->numRings();
  r.stereo = count_stereocenters(mol);

  int natoms = 0, carbons = 0, hetero = 0;
  for (const auto atom : mol.atoms()) {
    const int z = atom->getAtomicNum();
    natoms += 1 + atom->getTotalNumHs();   // H-inclusive atom count (Ghose)
    if (z == 6)
      ++carbons;
    else if (z != 1)
      ++hetero;
  }
  r.natoms = natoms;
  r.carbons = carbons;
  r.hetero = hetero;
}

void compute_rules(Record& r) {
  // Lipinski 1997 — Rule of Five
  r.lipinski = (r.mw <= 500 && r.hbd <= 5 && r.hba <= 10 && r.logp <= 5);
  // Veber 2002
  r.veber = (r.rotb <= 10 && r.tpsa <= 140);
  // Ghose 1999
  r.ghose = (r.mw >= 160 && r.mw <= 480 && r.logp >= -0.4 && r.logp <= 5.6 &&
             r.mr >= 40 && r.mr <= 130 && r.natoms >= 20 && r.natoms <= 70);
  // Egan 2000
  r.egan = (r.logp <= 5.88 && r.tpsa <= 131.6);
  // Muegge 2001
  r.muegge = (r.mw >= 200 && r.mw <= 600 && r.logp >= -2 && r.logp <= 5 &&
              r.tpsa <= 150 && r.rings <= 7 && r.carbons > 4 && r.hetero > 1 &&
              r.rotb <= 15 && r.hba <= 10 && r.hbd <= 5);
  // Rule of Three (Congreve 2003)
  r.ro3 = (r.mw <= 300 && r.logp <= 3 && r.hbd <= 3 && r.hba <= 3 &&
           r.rotb <= 3 && r.tpsa <= 60);
  // Lead-likeness (Teague 1999)
  r.leadlike = (r.mw <= 350 && r.logp <= 3.5);
  // REOS (Walters 2002)
  r.reos = (r.mw >= 200 && r.mw <= 500 && r.logp >= -5 && r.logp <= 5 &&
            r.hbd <= 5 && r.hba <= 10 && r.rotb <= 8 && r.heavy >= 15 &&
            r.heavy <= 50 && r.charge >= -2 && r.charge <= 2);
  // GSK 4/400 (Gleeson 2008) — reduced toxicity risk
  r.gsk = (r.logp < 4.0 && r.mw < 400.0);
  // Pfizer 3/75 (Hughes 2008) — reduced toxicity risk
  r.pfizer = (r.tpsa > 75.0 && r.logp < 3.0);
}

// ---------------------------------------------------------------------------
// Gates — split into three stages, ordered cheapest first.
//
//   (1) rules_pass     : descriptor-based; no substructure matching needed
//   (2) qed_pass       : QED threshold; requires PAINS+BRENK counts already
//   (3) alerts_pass    : substructure matches (lilly, inpharmatica, ames, ...)
//   ((4) egg_pass      : BOILED-Egg region + BBB-MPO)
// ---------------------------------------------------------------------------
bool rules_pass(const Record& r, const Options& o) {
  for (const auto& k : o.rules) {
    if (k == "lipinski" && !r.lipinski) return false;
    if (k == "veber"    && !r.veber)    return false;
    if (k == "ghose"    && !r.ghose)    return false;
    if (k == "egan"     && !r.egan)     return false;
    if (k == "muegge"   && !r.muegge)   return false;
    if (k == "ro3"      && !r.ro3)      return false;
    if (k == "leadlike" && !r.leadlike) return false;
    if (k == "reos"     && !r.reos)     return false;
    if (k == "gsk"      && !r.gsk)      return false;
    if (k == "pfizer"   && !r.pfizer)   return false;
  }
  if (o.fsp3_min_set   && r.fsp3   < o.fsp3_min)   return false;
  if (o.stereo_max_set && r.stereo > o.stereo_max) return false;
  return true;
}

bool qed_pass(const Record& r, const Options& o) {
  return !o.qed_min_set || r.qed >= o.qed_min;
}

bool alerts_pass(const Record& r, const Options& o) {
  if (o.alerts.count("pains")        && !r.pains_free)        return false;
  if (o.alerts.count("brenk")        && !r.brenk_free)        return false;
  if (o.alerts.count("lilly")        && !r.lilly_free)        return false;
  if (o.alerts.count("inpharmatica") && !r.inpharmatica_free) return false;
  if (o.alerts.count("ames")         && !r.ames_free)         return false;
  if (o.alerts.count("aggregators")  && !r.aggregators_free)  return false;
  if (o.alerts.count("extra")        && !r.extra_free)        return false;
  if (o.extra_alerts                 && !r.extra_free)        return false;  // back-compat
  return true;
}

bool egg_pass(const Record& r, const Options& o) {
  if (o.permeability == "bbb"          && !r.bbb)              return false;
  if (o.permeability == "gia"          && !r.gia)              return false;
  if (o.permeability == "bbb_and_gia"  && !(r.bbb && r.gia))   return false;
  if (o.bbb_threshold_set && !(r.bbb_mpo > o.bbb_threshold))   return false;
  return true;
}

// ---------------------------------------------------------------------------
// Fingerprint creation. AtomPair/RDKit generators are deprecated but still
// work; the deprecation warning is suppressed locally.
// ---------------------------------------------------------------------------
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
ExplicitBitVect* make_fingerprint(const RDKit::ROMol& mol, const std::string& y,
                                  int bits, int radius) {
  if (y == "morgan")
    return RDKit::MorganFingerprints::getFingerprintAsBitVect(
        mol, static_cast<unsigned>(radius), static_cast<unsigned>(bits));
  if (y == "aptt")
    return RDKit::AtomPairs::getHashedAtomPairFingerprintAsBitVect(
        mol, static_cast<unsigned>(bits));
  if (y == "rd")
    return RDKit::RDKFingerprintMol(mol, 1, 7, static_cast<unsigned>(bits));
  return nullptr;
}
#pragma GCC diagnostic pop

double similarity_of(const ExplicitBitVect& a, const ExplicitBitVect& b,
                     const std::string& m) {
  if (m == "Dice")         return DiceSimilarity(a, b);
  if (m == "Sokal")        return SokalSimilarity(a, b);
  if (m == "McConnaughey") return McConnaugheySimilarity(a, b);
  if (m == "Kulczynski")   return KulczynskiSimilarity(a, b);
  if (m == "Cosine")       return CosineSimilarity(a, b);
  return TanimotoSimilarity(a, b);   // default
}

}  // anonymous namespace

// ============================================================================
// SharedCatalogs — built once, used by every thread.
// ============================================================================
struct SharedCatalogs::Impl {
  std::unique_ptr<RDKit::FilterCatalog> pains, brenk, lilly, inpharmatica, extra;
  std::vector<std::unique_ptr<RDKit::ROMol>> ames_patterns;
  std::vector<std::unique_ptr<RDKit::ROMol>> aggregator_patterns;
  bool use_extra = false;
};

SharedCatalogs::SharedCatalogs(bool include_extra)
    : impl_(std::make_unique<Impl>()) {
  using FCP = RDKit::FilterCatalogParams;
  { FCP p; p.addCatalog(FCP::PAINS);              impl_->pains       = std::make_unique<RDKit::FilterCatalog>(p); }
  { FCP p; p.addCatalog(FCP::BRENK);              impl_->brenk       = std::make_unique<RDKit::FilterCatalog>(p); }
  { FCP p; p.addCatalog(FCP::CHEMBL_LINT);        impl_->lilly       = std::make_unique<RDKit::FilterCatalog>(p); }
  { FCP p; p.addCatalog(FCP::CHEMBL_Inpharmatica); impl_->inpharmatica = std::make_unique<RDKit::FilterCatalog>(p); }
  impl_->ames_patterns =
      compile_smarts(AMES_SMARTS, sizeof(AMES_SMARTS) / sizeof(*AMES_SMARTS));
  impl_->aggregator_patterns = compile_smarts(
      AGGREGATOR_SMARTS,
      sizeof(AGGREGATOR_SMARTS) / sizeof(*AGGREGATOR_SMARTS));
  if (include_extra) {
    FCP p;
    p.addCatalog(FCP::NIH);
    p.addCatalog(FCP::ZINC);
    p.addCatalog(FCP::CHEMBL);
    impl_->extra = std::make_unique<RDKit::FilterCatalog>(p);
    impl_->use_extra = true;
  }
}

SharedCatalogs::~SharedCatalogs() = default;
bool SharedCatalogs::has_pains       (const RDKit::ROMol& m) const { return impl_->pains       ->hasMatch(m); }
bool SharedCatalogs::has_brenk       (const RDKit::ROMol& m) const { return impl_->brenk       ->hasMatch(m); }
bool SharedCatalogs::has_lilly       (const RDKit::ROMol& m) const { return impl_->lilly       ->hasMatch(m); }
bool SharedCatalogs::has_inpharmatica(const RDKit::ROMol& m) const { return impl_->inpharmatica->hasMatch(m); }
bool SharedCatalogs::has_ames        (const RDKit::ROMol& m) const { return any_substruct_match(m, impl_->ames_patterns); }
bool SharedCatalogs::has_aggregators (const RDKit::ROMol& m) const { return any_substruct_match(m, impl_->aggregator_patterns); }
bool SharedCatalogs::extra_enabled() const { return impl_->use_extra; }
bool SharedCatalogs::has_extra(const RDKit::ROMol& m) const {
  return impl_->use_extra && impl_->extra->hasMatch(m);
}
int SharedCatalogs::alert_count(const RDKit::ROMol& m) const {
  // QED's ALERTS input: number of PAINS + BRENK substructure matches.
  int n = 0;
  n += static_cast<int>(impl_->pains->getMatches(m).size());
  n += static_cast<int>(impl_->brenk->getMatches(m).size());
  return n;
}

// ============================================================================
// ChemContext — per-thread standardizer set + a borrowed catalog reference.
// ============================================================================
struct ChemContext::Impl {
  RDKit::MolStandardize::LargestFragmentChooser lfc;
  RDKit::MolStandardize::Normalizer normalizer;
  RDKit::MolStandardize::Uncharger uncharger;
  RDKit::MolStandardize::TautomerEnumerator taut;
};

ChemContext::ChemContext(const Options&, const SharedCatalogs& cats)
    : impl_(std::make_unique<Impl>()), cats_(&cats) {}

ChemContext::~ChemContext() = default;
ChemContext::ChemContext(ChemContext&&) noexcept = default;

Record ChemContext::process(const std::string& raw, long line_no, const Options& o,
                            const ExplicitBitVect* target_fp) const {
  Record r;
  r.line_no = line_no;

  const std::string cleaned = trim(raw);
  if (cleaned.empty()) return r;
  const std::string smi = cleaned.substr(0, cleaned.find_first_of(" \t"));
  if (smi.empty()) return r;

  std::unique_ptr<RDKit::ROMol> parsed;
  try {
    parsed.reset(RDKit::SmilesToMol(smi));
  } catch (...) {
    parsed.reset();
  }
  if (!parsed) return r;

  // Standardization pipeline: largest fragment → normalize → uncharge →
  // (canonical tautomer). The tautomer step dominates and is opt-out.
  std::unique_ptr<RDKit::ROMol> mol;
  if (o.standardize) {
    try {
      std::unique_ptr<RDKit::ROMol> a(impl_->lfc.choose(*parsed));
      std::unique_ptr<RDKit::ROMol> b(impl_->normalizer.normalize(*a));
      std::unique_ptr<RDKit::ROMol> c(impl_->uncharger.uncharge(*b));
      if (o.tautomer) {
        mol.reset(impl_->taut.canonicalize(*c));
      } else {
        mol = std::move(c);
      }
    } catch (...) {
      mol.reset();
    }
  }
  if (!mol) mol = std::move(parsed);

  compute_descriptors(*mol, r);
  r.smiles = RDKit::MolToSmiles(*mol);
  compute_rules(r);

  // Early-exit ordering: cheapest gate first.
  if (!rules_pass(r, o)) return r;

  // PAINS + BRENK are needed both for the alerts gate and for the QED
  // ALERTS input. Pay this once and cache.
  r.pains_free = !cats_->has_pains(*mol);
  r.brenk_free = !cats_->has_brenk(*mol);
  r.qed = compute_qed(r.mw, r.logp, r.hba, r.hbd, r.tpsa, r.rotb, r.arom,
                      cats_->alert_count(*mol));

  if (!qed_pass(r, o)) return r;

  // Optional alert catalogs are checked only if the user selected them.
  if (o.alerts.count("lilly"))        r.lilly_free        = !cats_->has_lilly(*mol);
  if (o.alerts.count("inpharmatica")) r.inpharmatica_free = !cats_->has_inpharmatica(*mol);
  if (o.alerts.count("ames"))         r.ames_free         = !cats_->has_ames(*mol);
  if (o.alerts.count("aggregators"))  r.aggregators_free  = !cats_->has_aggregators(*mol);
  if (cats_->extra_enabled())         r.extra_free        = !cats_->has_extra(*mol);

  if (!alerts_pass(r, o)) return r;

  // BOILED-Egg region + BBB-MPO
  BPoint pt(r.tpsa_sp, r.logp);
  r.gia = bg::within(pt, gia_polygon());
  r.bbb = bg::within(pt, bbb_polygon());
  r.bbb_mpo = bbb_mpo_score(r);

  if (!egg_pass(r, o)) return r;

  std::unique_ptr<ExplicitBitVect> fp(
      make_fingerprint(*mol, o.method, o.fp_bits, o.fp_radius));
  if (fp && target_fp) r.similarity = similarity_of(*fp, *target_fp, o.similarity);

  r.valid = true;
  return r;
}

// Stage-1 fast path for two-stage screening.
Record ChemContext::process_prefilter(const std::string& raw, long line_no,
                                      const Options& o,
                                      const ExplicitBitVect* target_fp) const {
  Record r;
  r.line_no = line_no;

  const std::string cleaned = trim(raw);
  if (cleaned.empty()) return r;
  const std::string smi = cleaned.substr(0, cleaned.find_first_of(" \t"));
  if (smi.empty()) return r;

  std::unique_ptr<RDKit::ROMol> parsed;
  try {
    parsed.reset(RDKit::SmilesToMol(smi));
  } catch (...) {
    parsed.reset();
  }
  if (!parsed) return r;

  r.smiles = smi;   // raw SMILES — full pipeline will canonicalize in stage 2
  std::unique_ptr<ExplicitBitVect> fp(
      make_fingerprint(*parsed, o.method, o.fp_bits, o.fp_radius));
  if (fp && target_fp) r.similarity = similarity_of(*fp, *target_fp, o.similarity);
  r.valid = true;
  return r;
}

// ============================================================================
// Free helpers
// ============================================================================
ExplicitBitVect* target_fingerprint(const std::string& smiles,
                                    const std::string& method,
                                    int bits, int radius, std::string& err) {
  std::unique_ptr<RDKit::ROMol> mol;
  try {
    mol.reset(RDKit::SmilesToMol(smiles));
  } catch (...) {
    mol.reset();
  }
  if (!mol) {
    err = "invalid target SMILES: " + smiles;
    return nullptr;
  }
  ExplicitBitVect* fp = make_fingerprint(*mol, method, bits, radius);
  if (!fp) err = "invalid fingerprint method: " + method;
  return fp;
}

void disable_rdkit_logs() { boost::logging::disable_logs("rdApp.*"); }
