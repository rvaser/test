// Copyright (c) 2020 Robert Vaser

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "ram/minimizer_engine.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace {

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path) {
  auto is_suffix = [] (const std::string& str, const std::string& suff) -> bool {  // NOLINT
    return str.size() < suff.size() ? false :
        str.compare(str.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fa")    || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fastq.gz") ||
      is_suffix(path, ".fq")    || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastqParser>(path);  // NOLINT
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
  auto r = rp->Parse(-1);

  auto sp = CreateParser(argv[2]);
  auto s = sp->Parse(-1);

  auto me = ram::MinimizerEngine();
  me.Minimize(r.begin(), r.end());

  std::vector<double> l;
  for (const auto& it : s) {
    auto o = me.Map(it, false, false);
    if (o.empty()) {
      continue;
    }

    std::sort(o.begin(), o.end(),
        [] (const biosoup::Overlap& lhs, const biosoup::Overlap& rhs) -> bool {
          return lhs.lhs_end - lhs.lhs_begin > rhs.lhs_end - rhs.lhs_begin;
        });
    auto b = o.front();

    if (b.lhs_end - b.lhs_begin / static_cast<double>(it->data.size()) > 0.98) {
      auto n = biosoup::Sequence("", r[b.rhs_id]->data.substr(b.rhs_begin, b.rhs_end - b.rhs_begin));  // NOLINT
      if (b.strand) {
        n.ReverseAndComplement();
      }
      if (!it->quality.empty()) {
        std::cout << "@" << it->name << std::endl
                  << n.data << std::endl
                  << "+" << std::endl
                  << std::string(n.data.size(), '^') << std::endl;
      } else {
        std::cout << ">" << it->name << std::endl
                  << n.data << std::endl;
      }
    }
  }

  return 0;
}
