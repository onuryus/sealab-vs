// SEALab Virtual Screening — main entry point.
//
// Full per-molecule flow (parallelized across cores via OpenMP):
//   parse SMILES → standardize → descriptors → rule gates →
//   PAINS/BRENK (+ QED) → optional alert catalogs → BOILED-Egg + BBB-MPO →
//   fingerprint → similarity to target → CSV row.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <omp.h>
#include <DataStructs/ExplicitBitVect.h>

#include "chem.hpp"
#include "options.hpp"

// ============================================================================
// Allowed rule and alert identifiers (validated at CLI parse time).
// ============================================================================
const std::set<std::string>& Options::valid_rules() {
  static const std::set<std::string> v = {
      "lipinski", "veber", "ghose",    "egan",   "muegge",
      "ro3",      "leadlike", "reos",  "gsk",    "pfizer"};
  return v;
}

const std::set<std::string>& Options::valid_alerts() {
  static const std::set<std::string> v = {
      "pains", "brenk",        "lilly", "inpharmatica",
      "ames",  "aggregators",  "extra"};
  return v;
}

// ============================================================================
// CLI parsing
// ============================================================================
void print_usage(const char* prog) {
  std::cerr <<
      "Usage: " << prog << " --lib <file.smi> --target <SMILES> [options]\n\n"
      "Required:\n"
      "  --lib <path>           Library to screen (.smi)\n"
      "  --target <SMILES>      Query molecule (similarity reference)\n\n"
      "Output:\n"
      "  -o, --out <path>       Output CSV (default: results.csv)\n\n"
      "Fingerprint:\n"
      "  --method <m>           morgan | aptt | rd (default: morgan)\n"
      "  --similarity <m>       Tanimoto | Dice | Sokal | McConnaughey |\n"
      "                         Kulczynski | Cosine (default: Tanimoto)\n"
      "  --fp-bits <N>          Fingerprint bit length (default: 2048)\n"
      "  --fp-radius <N>        Morgan radius (default: 2; ignored for aptt/rd)\n\n"
      "Rule gates (comma list; ALL listed rules must pass):\n"
      "  --rules <list>         lipinski, veber, ghose, egan, muegge, ro3,\n"
      "                         leadlike, reos, gsk (4/400), pfizer (3/75)\n"
      "                         (default: lipinski)\n\n"
      "Continuous-score gates:\n"
      "  --qed-min <x>          QED score >= x (Bickerton 2012; typical 0.5)\n"
      "  --fsp3-min <x>         Fsp3 >= x (Lovering 2009; typical 0.42)\n"
      "  --stereo-max <n>       Stereocenter count <= n (typical 4)\n\n"
      "Alert filters (comma list; each selected catalog must be CLEAN):\n"
      "  --alerts <list>        pains, brenk, lilly, inpharmatica,\n"
      "                         ames, aggregators, extra\n"
      "                         ('extra' bundles NIH+ZINC+CHEMBL)\n\n"
      "Permeability (BOILED-Egg + BBB-MPO):\n"
      "  --permeability <m>     bbb | gia | bbb_and_gia\n"
      "  --bbb-threshold <x>    Keep molecules with BBB-MPO > x\n\n"
      "Behavior:\n"
      "  --no-standardize       Skip the full standardization pipeline\n"
      "  --no-tautomer          Skip only the tautomer canonicalization step\n"
      "  --prefilter <K>        Two-stage screening: fast similarity-only\n"
      "                         pre-filter → top-K → full pipeline\n"
      "  --threads <n>          Thread count (0 = all cores)\n"
      "  -h, --help             Show this help\n";
}

namespace {
std::string lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(c));
  return s;
}

std::vector<std::string> split_csv(const std::string& s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    const auto b = item.find_first_not_of(" \t");
    if (b == std::string::npos) continue;
    const auto e = item.find_last_not_of(" \t");
    out.push_back(item.substr(b, e - b + 1));
  }
  return out;
}
}  // namespace

bool parse_args(int argc, char** argv, Options& o, std::string& err) {
  auto need = [&](int& i) -> const char* {
    if (i + 1 >= argc) {
      err = std::string("missing value for: ") + argv[i];
      return nullptr;
    }
    return argv[++i];
  };
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      print_usage(argv[0]);
      std::exit(0);
    } else if (a == "--lib") {
      auto v = need(i); if (!v) return false; o.lib = v;
    } else if (a == "--target") {
      auto v = need(i); if (!v) return false; o.target = v;
    } else if (a == "-o" || a == "--out") {
      auto v = need(i); if (!v) return false; o.out = v;
    } else if (a == "--method") {
      auto v = need(i); if (!v) return false; o.method = lower(v);
    } else if (a == "--similarity") {
      auto v = need(i); if (!v) return false; o.similarity = v;
    } else if (a == "--rules") {
      auto v = need(i); if (!v) return false;
      o.rules.clear();
      for (auto& k : split_csv(lower(v))) {
        if (!Options::valid_rules().count(k)) {
          err = "unknown rule: " + k; return false;
        }
        o.rules.insert(k);
      }
    } else if (a == "--alerts") {
      auto v = need(i); if (!v) return false;
      o.alerts.clear();
      for (auto& t : split_csv(lower(v))) {
        if (!Options::valid_alerts().count(t)) {
          err = "unknown alert catalog: " + t; return false;
        }
        o.alerts.insert(t);
      }
    } else if (a == "--qed-min") {
      auto v = need(i); if (!v) return false;
      o.qed_min = std::stod(v); o.qed_min_set = true;
    } else if (a == "--fsp3-min") {
      auto v = need(i); if (!v) return false;
      o.fsp3_min = std::stod(v); o.fsp3_min_set = true;
    } else if (a == "--stereo-max") {
      auto v = need(i); if (!v) return false;
      o.stereo_max = std::stoi(v); o.stereo_max_set = true;
    } else if (a == "--permeability") {
      auto v = need(i); if (!v) return false; o.permeability = lower(v);
    } else if (a == "--bbb-threshold") {
      auto v = need(i); if (!v) return false;
      o.bbb_threshold = std::stod(v); o.bbb_threshold_set = true;
    } else if (a == "--threads") {
      auto v = need(i); if (!v) return false; o.threads = std::stoi(v);
    } else if (a == "--extra-alerts") {
      o.extra_alerts = true;
    } else if (a == "--no-standardize") {
      o.standardize = false;
    } else if (a == "--no-tautomer") {
      o.tautomer = false;
    } else if (a == "--fp-bits") {
      auto v = need(i); if (!v) return false;
      o.fp_bits = std::stoi(v);
      if (o.fp_bits < 64) { err = "--fp-bits must be >= 64"; return false; }
    } else if (a == "--fp-radius") {
      auto v = need(i); if (!v) return false;
      o.fp_radius = std::stoi(v);
      if (o.fp_radius < 0) { err = "--fp-radius must be >= 0"; return false; }
    } else if (a == "--prefilter") {
      auto v = need(i); if (!v) return false;
      o.prefilter_top = std::stol(v);
      if (o.prefilter_top < 0) { err = "--prefilter must be >= 0"; return false; }
    } else {
      err = "unknown argument: " + a;
      return false;
    }
  }
  if (o.lib.empty())    { err = "--lib is required";    return false; }
  if (o.target.empty()) { err = "--target is required"; return false; }
  return true;
}

// ============================================================================
// CSV output
// ============================================================================
namespace {
void write_csv(const std::string& path, const std::vector<Record>& hits) {
  std::ofstream out(path);
  out << "rank,smiles,similarity,mw,logp,tpsa,hbd,hba,rotb,arom,fsp3,heavy,"
         "stereo,qed,bbb_mpo,gia,bbb,"
         "pains_free,brenk_free,lilly_free,inpharmatica_free,"
         "ames_free,aggregators_free,extra_free,"
         "lipinski,veber,ghose,egan,muegge,ro3,leadlike,reos,gsk,pfizer\n";
  out.setf(std::ios::fixed);
  long rank = 0;
  for (const auto& r : hits) {
    out.precision(4);
    out << ++rank << ',' << r.smiles << ',' << r.similarity << ',';
    out.precision(2);
    out << r.mw << ',' << r.logp << ',' << r.tpsa << ',' << r.hbd << ','
        << r.hba << ',' << r.rotb << ',' << r.arom << ',' << r.fsp3 << ','
        << r.heavy << ',' << r.stereo << ',';
    out.precision(3);
    out << r.qed << ',';
    out.precision(2);
    out << r.bbb_mpo << ',' << (r.gia ? 1 : 0) << ','
        << (r.bbb ? 1 : 0) << ','
        << (r.pains_free ? 1 : 0) << ',' << (r.brenk_free ? 1 : 0) << ','
        << (r.lilly_free ? 1 : 0) << ',' << (r.inpharmatica_free ? 1 : 0) << ','
        << (r.ames_free ? 1 : 0) << ',' << (r.aggregators_free ? 1 : 0) << ','
        << (r.extra_free ? 1 : 0) << ','
        << (r.lipinski ? 1 : 0) << ',' << (r.veber ? 1 : 0) << ','
        << (r.ghose ? 1 : 0) << ',' << (r.egan ? 1 : 0) << ','
        << (r.muegge ? 1 : 0) << ',' << (r.ro3 ? 1 : 0) << ','
        << (r.leadlike ? 1 : 0) << ',' << (r.reos ? 1 : 0) << ','
        << (r.gsk ? 1 : 0) << ',' << (r.pfizer ? 1 : 0) << '\n';
  }
}
}  // namespace

// ============================================================================
// Splash banner — one of two ASCII screens is shown at startup, picked
// pseudo-randomly so the program has a bit of personality.
// ============================================================================
namespace {

constexpr const char* kBanner1 = R"banner(
           _____ ______          _           _
          / ____|  ____|   /\   | |         | |
         | (___ | |__     /  \  | |     __ _| |__
          \___ \|  __|   / /\ \ | |    / _` | '_ \
          ____) | |____ / ____ \| |___| (_| | |_) |
         |_____/|______/_/    \_\______\__,_|_.__/

   ligand-based virtual screening   ·   C++ / RDKit / OpenMP


 _._     _,-'""`-._                  __________
(,-.`._,'(       |\`-/|             / ___  ___ \
    `-.-' \ )-`( , o o)            / / @ \/ @ \ \
          `-    \`_`"'-            \ \___/\___/ /
                                    \____\/____/

   A Sealab Project
   Created by   : Hacı Aslan Onur İşcil, Ecem Bulut Turhan, Saliha Ece Acuner
)banner";

constexpr const char* kBanner2 = R"banner(
   _      _      _      _      _      _      _      _      _      _
 _( )_  _( )_  _( )_  _( )_  _( )_  _( )_  _( )_  _( )_  _( )_  _( )_
(_ o _)(_ o _)(_ o _)(_ o _)(_ o _)(_ o _)(_ o _)(_ o _)(_ o _)(_ o _)
 (_,_)  (_,_)  (_,_)  (_,_)  (_,_)  (_,_)  (_,_)  (_,_)  (_,_)  (_,_)

           _____ ______          _           _
          / ____|  ____|   /\   | |         | |
         | (___ | |__     /  \  | |     __ _| |__
          \___ \|  __|   / /\ \ | |    / _` | '_ \
          ____) | |____ / ____ \| |___| (_| | |_) |
         |_____/|______/_/    \_\______\__,_|_.__/

           ..
          ( '`<
           )(
    ( ----'  '.       _          _          _          _          _
    (         ;     >(')____,  >(')____,  >(')____,  >(')____,  >(')___,
     (_______,'      (` =~~/    (` =~~/    (` =~~/    (` =~~/    (` =~~/
~^~^~^~^~^~^~^~^~^~^^~^`---'~^~^`---'~^~^`---'~^~^`---'~^~^`---'~^~^~

   A Sealab Project
   Created by   : Hacı Aslan Onur İşcil, Ecem Bulut Turhan, Saliha Ece Acuner
)banner";

void print_intro() {
  const bool pick_two = (std::time(nullptr) % 2) == 0;
  std::fputs(pick_two ? kBanner2 : kBanner1, stderr);
  std::fputc('\n', stderr);
  std::fflush(stderr);
}

}  // namespace

// ============================================================================
// Line count helper (single pass, only counts non-empty lines, so the
// progress percentage matches what the workers actually process).
// ============================================================================
long count_non_empty_lines(const std::string& path) {
  std::ifstream in(path);
  if (!in) return -1;
  long n = 0;
  std::string line;
  while (std::getline(in, line))
    if (!line.empty()) ++n;
  return n;
}

// ============================================================================
// Single-line progress reporter (\r-redrawn). If `expected` is positive,
// shows percentage and ETA; otherwise just an elapsed counter.
// ============================================================================
void print_progress(const char* label, long done, long expected,
                    std::chrono::steady_clock::time_point start,
                    bool final_line) {
  using namespace std::chrono;
  const auto now = steady_clock::now();
  const double elapsed = duration_cast<milliseconds>(now - start).count() / 1000.0;
  const double rate = elapsed > 0 ? done / elapsed : 0;
  if (expected > 0) {
    const double pct = 100.0 * done / expected;
    const long eta = (rate > 0 && expected > done)
                         ? static_cast<long>((expected - done) / rate)
                         : 0;
    std::fprintf(stderr,
                 "\r%s: %ld / %ld  (%5.1f%%)  %6.0f mol/s  ETA %4lds   ",
                 label, done, expected, pct, rate, eta);
  } else {
    std::fprintf(stderr, "\r%s: %ld processed  %6.0f mol/s   ", label, done, rate);
  }
  if (final_line) std::fputc('\n', stderr);
  std::fflush(stderr);
}

// ============================================================================
// Generic streaming driver: reads the .smi file in locked batches and
// dispatches each line through a caller-provided process function.
// process_fn is either ChemContext::process (full pipeline) or
// ChemContext::process_prefilter (stage-1 fast path).
// ============================================================================
template <typename ProcessFn>
void stream_process(const std::string& path, const Options& o,
                    const SharedCatalogs& cats,
                    const ExplicitBitVect* target_fp, ProcessFn process_fn,
                    std::vector<Record>& sink, long& total_lines,
                    long& parsed_ok, const char* label, long expected_total) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "Could not open file: " << path << "\n";
    return;
  }

  total_lines = 0;
  parsed_ok = 0;
  constexpr int BATCH = 256;

  std::atomic<long> done_count{0};
  const auto t_start = std::chrono::steady_clock::now();
  auto t_last_print = t_start;

#pragma omp parallel
  {
    ChemContext ctx(o, cats);
    std::vector<Record> local;
    std::vector<std::pair<long, std::string>> batch;
    batch.reserve(BATCH);
    long local_parsed_ok = 0;
    bool done = false;

    while (!done) {
      batch.clear();
#pragma omp critical(stream_read)
      {
        std::string line;
        while ((int)batch.size() < BATCH && std::getline(in, line)) {
          ++total_lines;
          if (!line.empty())
            batch.emplace_back(total_lines, std::move(line));
        }
        if (batch.empty()) done = true;

        // Progress is printed inside the already-serialized read section,
        // so it does not need an additional lock.
        const auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(
                now - t_last_print).count() >= 300) {
          t_last_print = now;
          print_progress(label, done_count.load(), expected_total, t_start,
                         /*final_line=*/false);
        }
      }
      if (batch.empty()) break;

      for (auto& [ln, line] : batch) {
        Record r = process_fn(ctx, line, ln, o, target_fp);
        if (!r.smiles.empty()) ++local_parsed_ok;
        if (r.valid) local.push_back(std::move(r));
      }
      done_count.fetch_add(static_cast<long>(batch.size()));
    }

#pragma omp critical(stream_merge)
    {
      parsed_ok += local_parsed_ok;
      sink.insert(sink.end(), std::make_move_iterator(local.begin()),
                  std::make_move_iterator(local.end()));
    }
  }

  print_progress(label, done_count.load(), expected_total, t_start,
                 /*final_line=*/true);
}

// ============================================================================
int main(int argc, char** argv) {
  Options o;
  std::string err;
  if (!parse_args(argc, argv, o, err)) {
    std::cerr << "Error: " << err << "\n\n";
    print_usage(argv[0]);
    return 1;
  }

  disable_rdkit_logs();
  if (o.threads > 0) omp_set_num_threads(o.threads);

  print_intro();

  std::cout << "Library:    " << o.lib << "\n";
  std::cout << "Target:     " << o.target << "\n";
  std::cout << "Method:     " << o.method << " / " << o.similarity
            << " (bits=" << o.fp_bits;
  if (o.method == "morgan") std::cout << ", radius=" << o.fp_radius;
  std::cout << ")\n";
  if (o.prefilter_top > 0)
    std::cout << "Two-stage screening: prefilter top-" << o.prefilter_top << "\n";

  std::unique_ptr<ExplicitBitVect> target_fp(target_fingerprint(
      o.target, o.method, o.fp_bits, o.fp_radius, err));
  if (!target_fp) {
    std::cerr << "Error: " << err << "\n";
    return 1;
  }

  const auto tc0 = std::chrono::steady_clock::now();
  const SharedCatalogs cats(o.extra_alerts || o.alerts.count("extra"));
  const auto tc1 = std::chrono::steady_clock::now();
  {
    ChemContext probe(o, cats);
    const auto tc2 = std::chrono::steady_clock::now();
    auto ms = [](auto a, auto b) {
      return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
    };
    std::cerr << "[time] catalog setup (shared, 1x): " << ms(tc0, tc1) << " ms\n";
    std::cerr << "[time] standardizer setup (per thread): " << ms(tc1, tc2)
              << " ms\n";
  }

  long total_lines = 0;
  long parsed_ok = 0;
  std::vector<Record> hits;

  std::cerr << "Counting lines... ";
  std::cerr.flush();
  const long expected = count_non_empty_lines(o.lib);
  if (expected >= 0)
    std::cerr << expected << " molecules\n";
  else
    std::cerr << "(could not open file; progress will lack percentage)\n";

  const auto t0 = std::chrono::steady_clock::now();

  if (o.prefilter_top > 0) {
    // -------- Stage 1: fast similarity-only prefilter --------
    std::vector<Record> stage1;
    long s1_lines = 0, s1_parsed = 0;
    stream_process(
        o.lib, o, cats, target_fp.get(),
        [](ChemContext& ctx, const std::string& line, long ln,
           const Options& opt, const ExplicitBitVect* tfp) {
          return ctx.process_prefilter(line, ln, opt, tfp);
        },
        stage1, s1_lines, s1_parsed, "Prefilter", expected);

    total_lines = s1_lines;

    const long K = std::min<long>(o.prefilter_top, (long)stage1.size());
    if (K < (long)stage1.size()) {
      std::nth_element(stage1.begin(), stage1.begin() + K, stage1.end(),
                       [](const Record& a, const Record& b) {
                         return a.similarity > b.similarity;
                       });
      stage1.resize(K);
    }

    const auto t_s1 = std::chrono::steady_clock::now();
    std::cerr << "[time] Stage 1 (prefilter): "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t_s1 - t0).count()
              << " ms, " << s1_parsed << " parsed, top-" << K << " kept\n";

    // -------- Stage 2: full pipeline on top-K --------
    std::atomic<long> s2_done{0};
    const auto t_s2_start = std::chrono::steady_clock::now();
    auto t_last = t_s2_start;
#pragma omp parallel
    {
      ChemContext ctx(o, cats);
      std::vector<Record> local;
#pragma omp for schedule(dynamic, 16)
      for (long i = 0; i < (long)stage1.size(); ++i) {
        Record r = ctx.process(stage1[i].smiles, stage1[i].line_no, o,
                               target_fp.get());
        if (r.valid) local.push_back(std::move(r));
        const long d = s2_done.fetch_add(1) + 1;
        const auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(now - t_last)
                .count() >= 300) {
#pragma omp critical(stage2_print)
          {
            t_last = now;
            print_progress("Full pipeline", d, (long)stage1.size(), t_s2_start,
                           false);
          }
        }
      }
#pragma omp critical(stage2_merge)
      {
        hits.insert(hits.end(), std::make_move_iterator(local.begin()),
                    std::make_move_iterator(local.end()));
      }
    }
    print_progress("Full pipeline", s2_done.load(), (long)stage1.size(),
                   t_s2_start, true);
    parsed_ok = s1_parsed;
  } else {
    // -------- Classic single-stage: full pipeline streamed --------
    stream_process(
        o.lib, o, cats, target_fp.get(),
        [](ChemContext& ctx, const std::string& line, long ln,
           const Options& opt, const ExplicitBitVect* tfp) {
          return ctx.process(line, ln, opt, tfp);
        },
        hits, total_lines, parsed_ok, "Screening", expected);
  }

  const auto t1 = std::chrono::steady_clock::now();
  const double secs =
      std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() / 1000.0;

  // Dedup by canonical SMILES (keep first occurrence).
  {
    std::unordered_set<std::string> seen;
    std::vector<Record> dedup;
    dedup.reserve(hits.size());
    for (auto& r : hits)
      if (seen.insert(r.smiles).second) dedup.push_back(std::move(r));
    hits = std::move(dedup);
  }

  std::sort(hits.begin(), hits.end(),
            [](const Record& a, const Record& b) {
              return a.similarity > b.similarity;
            });

  write_csv(o.out, hits);

  std::cout << "\n--- Summary ---\n";
  std::cout << "Total lines       : " << total_lines << "\n";
  std::cout << "Parsed molecules  : " << parsed_ok << " / " << total_lines << "\n";
  std::cout << "Passed filters    : " << hits.size() << " (after dedup)\n";
  std::cout << "Time              : " << secs << " s  ("
            << (secs > 0 ? static_cast<long>(parsed_ok / secs) : 0)
            << " mol/s)\n";
  std::cout << "Output            : " << o.out << "\n";
  if (!hits.empty())
    std::cout << "Highest similarity: " << hits.front().similarity << "  ("
              << hits.front().smiles << ")\n";
  return 0;
}
