/** \file
 *
 * Unit tests for kff.cpp, which provides tools for working with KFF files.
 */

#include "../kff.hpp"
#include "../utility.hpp"

#include "catch.hpp"

#include <algorithm>

namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

std::vector<std::string> kmers {
    "GATTA",
    "CATTAC",
    "GATTACA",
    "ATTACCAT",
};

void encode_decode(const std::string& kmer, const uint8_t* encoding) {
    std::vector<uint8_t> encoded = kff_encode(kmer, encoding);
    std::string decoding = kff_invert(encoding);
    std::string decoded = kff_decode(encoded.data(), kmer.length(), decoding);
    REQUIRE(decoded == kmer);
}

void minimizer_recode(const std::string& kmer, const uint8_t* encoding) {
    gbwtgraph::Key64 encoded = gbwtgraph::Key64::encode(kmer);
    std::vector<uint8_t> recoded = kff_recode(encoded.get_key(), kmer.length(), encoding);
    std::string decoding = kff_invert(encoding);
    std::string decoded = kff_decode(recoded.data(), kmer.length(), decoding);
    REQUIRE(decoded == kmer);

    kff_recoding_t recoding = kff_recoding(encoding);
    gbwtgraph::Key64 round_trip(kff_recode(recoded.data(), kmer.length(), recoding));
    REQUIRE(round_trip == encoded);
}

void recode_trivial(const std::string& kmer) {
    uint8_t encoding[4] = { 0, 1, 2, 3 };
    REQUIRE(kff_is_trivial(encoding));
    std::vector<uint8_t> encoded = kff_encode(kmer, encoding);

    kff_recoding_t recoding = kff_recoding(encoding);
    gbwtgraph::Key64::value_type truth = kff_recode(encoded.data(), kmer.length(), recoding);

    size_t bytes = kff_bytes(kmer.length());
    gbwtgraph::Key64::value_type recoded = kff_recode_trivial(encoded.data(), kmer.length(), bytes);
    REQUIRE(recoded == truth);
}

void rev_comp(const std::string& kmer, const uint8_t* encoding) {
    std::string expected = reverse_complement(kmer);

    {
        std::vector<uint8_t> encoded = kff_encode(kmer, encoding);
        std::vector<uint8_t> complemented = kff_reverse_complement(encoded.data(), kmer.length(), encoding);
        std::string decoding = kff_invert(encoding);
        std::string kff_decoded = kff_decode(complemented.data(), kmer.length(), decoding);
        REQUIRE(kff_decoded == expected);
    }

    {
        gbwtgraph::Key64 encoded = gbwtgraph::Key64::encode(kmer);
        gbwtgraph::Key64::value_type complemented = minimizer_reverse_complement(encoded.get_key(), kmer.length());
        std::string minimizer_decoded = gbwtgraph::Key64(complemented).decode(kmer.length());
        REQUIRE(minimizer_decoded == expected);
    }
}

void recode_block(const std::string& block, const std::vector<std::string>& true_kmers, const uint8_t* encoding) {
    size_t k = true_kmers.front().length();
    std::vector<uint8_t> encoded = kff_encode(block, encoding);
    kff_recoding_t recoding = kff_recoding(encoding);
    std::vector<gbwtgraph::Key64::value_type> recoded = kff_recode(encoded.data(), true_kmers.size(), k, recoding);
    REQUIRE(recoded.size() == true_kmers.size());
    for (size_t i = 0; i < recoded.size(); i++) {
        std::string decoded = gbwtgraph::Key64(recoded[i]).decode(k);
        REQUIRE(decoded == true_kmers[i]);
    }
}

// A raw section of a KFF file.
struct KffSection {
    // Maximum number of kmers in a block.
    size_t max_kmers_per_block;

    // Sequence of each block.
    std::vector<std::string> blocks;
};

// Writes raw sections of `k`-mers with the given encoding to a temporary file.
// Each kmer gets a different count. Returns the file name and the kmers with
// their counts in the minimizer index format.
std::pair<std::string, std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>>> write_kff(
    const std::vector<KffSection>& sections, size_t k, const uint8_t* encoding
) {
    std::string filename = temp_file::create();
    std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> kmers;
    Kff_file file(filename, "w");
    file.write_encoding(encoding[0], encoding[1], encoding[2], encoding[3]);
    for (const KffSection& section : sections) {
        Section_GV variables(&file);
        variables.write_var("k", k);
        variables.write_var("max", section.max_kmers_per_block);
        variables.write_var("data_size", 1);
        variables.close();

        Section_Raw raw(&file);
        for (const std::string& block : section.blocks) {
            std::vector<uint8_t> counts;
            for (size_t start = 0; start + k <= block.length(); start++) {
                counts.push_back(kmers.size() + 1);
                kmers.emplace_back(gbwtgraph::Key64::encode(block.substr(start, k)).get_key(), counts.back());
            }
            std::vector<uint8_t> encoded = kff_encode(block, encoding);
            raw.write_compacted_sequence(encoded.data(), block.length(), counts.data());
        }
        raw.close();
    }
    file.close();
    std::sort(kmers.begin(), kmers.end());
    return { filename, kmers };
}

// Reads all kmers from the file, `n` at a time, using the given number of threads.
std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> read_kff(const std::string& filename, size_t n, size_t threads) {
    ParallelKFFReader reader(filename);
    std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> result;
    #pragma omp parallel num_threads(threads)
    {
        while (true) {
            std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> block = reader.read(n);
            if (block.empty()) {
                break;
            }
            #pragma omp critical
            {
                result.insert(result.end(), block.begin(), block.end());
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

void check_kff_reader(const std::vector<KffSection>& sections, size_t k, const uint8_t* encoding) {
    auto written = write_kff(sections, k, encoding);
    for (size_t n : { 1, 2, 3, 100 }) {
        for (size_t threads : { 1, 4 }) {
            auto read = read_kff(written.first, n, threads);
            REQUIRE(read == written.second);
        }
    }
    temp_file::remove(written.first);
}

} // Anonymous namespace.

//------------------------------------------------------------------------------

TEST_CASE("Encode and decode", "[kff]") {
    SECTION("default encoding") {
        uint8_t encoding[4] = { 0, 1, 2, 3 };
        REQUIRE(kff_is_trivial(encoding));
        for (auto& kmer : kmers) {
            encode_decode(kmer, encoding);
        }
    }

    SECTION("kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        REQUIRE(!kff_is_trivial(encoding));
        for (auto& kmer : kmers) {
            encode_decode(kmer, encoding);
        }
    }
}

TEST_CASE("Recode minimizers", "[kff]") {
    SECTION("default encoding") {
        uint8_t encoding[4] = { 0, 1, 2, 3 };
        for (auto& kmer : kmers) {
            minimizer_recode(kmer, encoding);
        }
    }

    SECTION("trivial recoding") {
        for (auto& kmer : kmers) {
            recode_trivial(kmer);
        }
    }

    SECTION("kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        for (auto& kmer : kmers) {
            minimizer_recode(kmer, encoding);
        }
    }
}

TEST_CASE("Parse integers", "[kff]") {
    SECTION("single byte") {
        uint8_t bytes[1] = { 0x12 };
        uint64_t expected = 0x12;
        uint64_t parsed = kff_parse(bytes, sizeof(bytes));
        REQUIRE(parsed == expected);
    }

    SECTION("multiple bytes") {
        uint8_t bytes[5] = { 0x12, 0x34, 0x56, 0x78, 0x9A };
        uint64_t expected = 0x123456789A;
        uint64_t parsed = kff_parse(bytes, sizeof(bytes));
        REQUIRE(parsed == expected);
    }

    SECTION("maximum length") {
        uint8_t bytes[8] = { 0x12, 0x34, 0x56, 0x78, 0x9A, 0xBC, 0xDE, 0xF0 };
        uint64_t expected = 0x123456789ABCDEF0;
        uint64_t parsed = kff_parse(bytes, sizeof(bytes));
        REQUIRE(parsed == expected);
    }
}

TEST_CASE("Reverse complements", "[kff]") {
    SECTION("default encoding") {
        uint8_t encoding[4] = { 0, 1, 2, 3 };
        for (auto& kmer : kmers) {
            rev_comp(kmer, encoding);
        }
    }

    SECTION("kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        for (auto& kmer : kmers) {
            rev_comp(kmer, encoding);
        }
    }
}

TEST_CASE("Recode blocks", "[kff]") {
    SECTION("single kmer, default encoding") {
        uint8_t encoding[4] = { 0, 1, 2, 3 };
        for (auto& block : kmers) {
            std::vector<std::string> true_kmers { block };
            recode_block(block, true_kmers, encoding);
        }
    }

    SECTION("single kmer, kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        for (auto& block : kmers) {
            std::vector<std::string> true_kmers { block };
            recode_block(block, true_kmers, encoding);
        }
    }

    SECTION("4-mers, default encoding") {
        uint8_t encoding[4] = { 0, 1, 2, 3 };
        size_t k = 4;
        for (auto& block : kmers) {
            std::vector<std::string> true_kmers;
            for (size_t start = 0; start + k <= block.length(); start++) {
                true_kmers.push_back(block.substr(start, k));
            }
            recode_block(block, true_kmers, encoding);
        }
    }

    SECTION("4-mers, kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        size_t k = 4;
        for (auto& block : kmers) {
            std::vector<std::string> true_kmers;
            for (size_t start = 0; start + k <= block.length(); start++) {
                true_kmers.push_back(block.substr(start, k));
            }
            recode_block(block, true_kmers, encoding);
        }
    }
}

TEST_CASE("Parallel KFF reader", "[kff]") {
    // Kmers of length 5 take 2 bytes, with padding in the first byte.
    size_t k = 5;
    KffSection single { 1, { "GATTA", "ATTAC", "TTACA", "CATTA" } };
    KffSection multiple { 4, { "GATTACA", "CATTACAT", "ATTAC" } };

    SECTION("one kmer per block, default encoding") {
        uint8_t encoding[4] = { 0, 1, 2, 3 };
        check_kff_reader({ single }, k, encoding);
    }

    SECTION("one kmer per block, kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        check_kff_reader({ single }, k, encoding);
    }

    SECTION("multiple kmers per block, kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        check_kff_reader({ multiple }, k, encoding);
    }

    SECTION("mixed sections, kff encoding") {
        uint8_t encoding[4] = { 0, 1, 3, 2 };
        check_kff_reader({ single, multiple, single }, k, encoding);
    }
}

//------------------------------------------------------------------------------

}
}
