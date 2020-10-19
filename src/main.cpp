// Copyright (c) 2020 Robert Vaser

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "ram/minimizer_engine.hpp"
#include "thread_pool/thread_pool.hpp"

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

  auto tp =  std::make_shared<thread_pool::ThreadPool>(12);

  auto me = ram::MinimizerEngine(15, 5, 500, 4, 100, 10000, tp);
  me.Minimize(r.begin(), r.end());

  std::vector<double> l;
  std::vector<std::future<void>> f;
  for (auto& it : s) {
    f.emplace_back(tp->Submit(
        [&] (decltype(it) jt) -> void {
          auto o = me.Map(jt, false, false);
          if (!o.empty()) {
            std::sort(o.begin(), o.end(),
                [] (const biosoup::Overlap& lhs, const biosoup::Overlap& rhs) -> bool {  // NOLINT
                  return lhs.lhs_end - lhs.lhs_begin > rhs.lhs_end - rhs.lhs_begin;  // NOLINT
                });
            auto b = o.front();

            if (b.lhs_end - b.lhs_begin / static_cast<double>(jt->data.size()) > 0.98) {  // NOLINT
              jt->data = r[b.rhs_id]->data.substr(b.rhs_begin, b.rhs_end - b.rhs_begin);  // NOLINT
              if (jt->quality.size()) {
                jt->quality = std::string(jt->data.size(), '^');
              }
              if (b.strand) {
                jt->ReverseAndComplement();
              }
            }
          }
        }, std::ref(it)));
  }
  for (std::uint32_t i = 0; i < f.size(); ++i) {
    f[i].wait();
    auto& it = s[i];
    if (!it->quality.empty()) {
      std::cout << "@" << it->name << std::endl
                << it->data << std::endl
                << "+" << std::endl
                << it->quality << std::endl;
    } else {
      std::cout << ">" << it->name << std::endl
                << it->data << std::endl;
    }
  }

  return 0;
}
