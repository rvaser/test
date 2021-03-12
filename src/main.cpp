// Copyright (c) 2020 Robert Vaser

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "ram/minimizer_engine.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};

namespace {

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(
    const std::string& path) {
  auto is_suffix = [] (const std::string& str, const std::string& suff) -> bool {  // NOLINT
    return str.size() < suff.size() ? false :
        str.compare(str.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fa")    || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fastq.gz") ||
      is_suffix(path, ".fq")    || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::NucleicAcid>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    return 0;
  }

  auto rp = CreateParser(argv[1]);
  if (rp == nullptr) {
    return 1;
  }
  auto r = rp->Parse(-1);

  auto sp = CreateParser(argv[2]);
  if (sp == nullptr) {
    return 1;
  }

  auto tp = std::make_shared<thread_pool::ThreadPool>();

  auto me = ram::MinimizerEngine(tp);
  me.Minimize(r.begin(), r.end());

  while (true) {
    auto s = sp->Parse(1ULL << 32);
    if (s.empty()) {
      break;
    }

    std::vector<std::future<void>> f;
    for (auto& it : s) {
      f.emplace_back(tp->Submit(
          [&] (decltype(it) jt) -> void {
            if (jt->inflated_len < 1000) {
              jt.reset();
              return;
            }

            auto o = me.Map(jt, false, false);
            if (o.empty()) {
              jt.reset();
              return;
            }

            // remove duplications
            std::sort(o.begin(), o.end(),
                [] (const biosoup::Overlap& lhs,
                    const biosoup::Overlap& rhs) -> bool {
                  return lhs.lhs_end - lhs.lhs_begin >
                         rhs.lhs_end - rhs.lhs_begin;
                });
            while (true) {
              bool ic = 0;
              for (std::uint32_t i = 1; i < o.size(); ++i) {
                for (std::uint32_t j = 0; j < i; ++j) {
                  if (o[j].lhs_begin <= o[i].lhs_begin && o[i].lhs_end <= o[j].lhs_end) {  // NOLINT
                    o.erase(o.begin() + i);
                    ic = 1;
                    break;
                  }
                }
                if (ic) {
                  break;
                }
              }
              if (!ic) {
                break;
              }
            }

            // combine
            std::sort(o.begin(), o.end(),
                [] (const biosoup::Overlap& lhs,
                    const biosoup::Overlap& rhs) -> bool {
                  return lhs.lhs_begin <  rhs.lhs_begin;
                });
            std::string d = "";
            for (const auto& it : o) {
              biosoup::NucleicAcid na(
                  "",
                  r[it.rhs_id]->InflateData(
                      it.rhs_begin,
                      it.rhs_end - it.rhs_begin));
              if (it.strand == 0) {
                na.ReverseAndComplement();
              }
              d += na.InflateData();
              if (d.size() > 0.98 * jt->inflated_len) {
                break;
              }
            }

            jt = std::unique_ptr<biosoup::NucleicAcid>(
                new biosoup::NucleicAcid(jt->name, d));
          },
          std::ref(it)));
    }
    for (const auto& it : f) {
      it.wait();
    }

    for (const auto& it : s) {
      if (it) {
        std::cout << ">" << it->name << std::endl
                  << it->InflateData() << std::endl;
      }
    }
  }

  return 0;
}
