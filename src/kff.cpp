#include "kff.hpp"

namespace vg {

//------------------------------------------------------------------------------

bool kff_is_trivial(const uint8_t* encoding) {
    for (size_t i = 0; i < 4; i++) {
        if (encoding[i] != i) {
            return false;
        }
    }
    return true;
}

std::string kff_invert(const uint8_t* encoding) {
    std::string result(4, ' ');
    result[encoding[0]] = 'A';
    result[encoding[1]] = 'C';
    result[encoding[2]] = 'G';
    result[encoding[3]] = 'T';
    return result;
}

kff_recoding_t kff_recoding(const uint8_t* encoding) {
    kff_recoding_t result;
    for (size_t i = 0; i < 4; i++) {
        result.data[encoding[i]] = i;
    }
    return result;
}

uint64_t kff_parse(const uint8_t* data, size_t bytes) {
    uint64_t value = 0;
    size_t shift = 8 * bytes;
    for (size_t i = 0; i < bytes; i++) {
        shift -= 8;
        value |= static_cast<uint64_t>(data[i]) << shift;
    }
    return value;
}

//------------------------------------------------------------------------------

// Encode up to 4 characters in one byte.
uint8_t kff_encode(const std::string& kmer, size_t start, size_t limit, const uint8_t* encoding) {
    uint8_t val = 0;
    for (size_t i = start; i < limit; i++) {
        val <<= 2;
        auto packed = gbwtgraph::KmerEncoding::CHAR_TO_PACK[static_cast<uint8_t>(kmer[i])];
        if (packed < 4) {
            val |= encoding[packed];
        }
    }
    return val;
}

std::vector<uint8_t> kff_encode(const std::string& kmer, const uint8_t* encoding) {
    std::vector<uint8_t> result;
    result.reserve(kff_bytes(kmer.length()));

    // If k is not a multiple of 4, KFF adds the padding to the high-order bits
    // of the first byte.
    size_t remainder = kmer.length() & 3;
    if (remainder > 0) {
        result.push_back(kff_encode(kmer, 0, remainder, encoding));
    }
    for (size_t i = remainder; i < kmer.length(); i += 4) {
        result.push_back(kff_encode(kmer, i, i + 4, encoding));
    }

    return result;
}

// Decode up to 4 characters from one byte
void kff_decode(uint8_t byte, size_t chars, const std::string& decoding, std::string& output) {
    size_t offset = 2 * chars;
    for (size_t i = 0; i < chars; i++) {
        offset -= 2;
        output.push_back(decoding[(byte >> offset) & 3]);
    }
}

std::string kff_decode(const uint8_t* kmer, size_t k, const std::string& decoding) {
    std::string result;
    result.reserve(k);

    size_t bytes = kff_bytes(k);
    size_t chars = k & 3;
    if (chars == 0) {
        chars = 4;
    }
    for (size_t i = 0; i < bytes; i++) {
        kff_decode(kmer[i], chars, decoding, result);
        chars = 4;
    }

    return result;
}

//------------------------------------------------------------------------------

// Recode up to 4 characters in one byte.
uint8_t kff_recode(gbwtgraph::Key64::value_type kmer, size_t k, size_t chars, const uint8_t* encoding) {
    size_t offset = 2 * k;
    uint8_t val = 0;
    for (size_t i = 0; i < chars; i++) {
        offset -= 2;
        val = (val << 2) | encoding[(kmer >> offset) & 3];
    }
    return val;
}

std::vector<uint8_t> kff_recode(gbwtgraph::Key64::value_type kmer, size_t k, const uint8_t* encoding) {
    std::vector<uint8_t> result;
    result.reserve(kff_bytes(k));

    size_t remainder = k & 3;
    if (remainder > 0) {
        result.push_back(kff_recode(kmer, k, remainder, encoding));
    }
    for (size_t i = remainder; i < k; i += 4) {
        result.push_back(kff_recode(kmer, k - i, 4, encoding));
    }

    return result;
}

gbwtgraph::Key64::value_type kff_recode(const uint8_t* kmer, size_t k, kff_recoding_t recoding) {
    gbwtgraph::Key64::value_type result = 0;

    size_t bytes = kff_bytes(k);
    size_t chars = k & 3;
    if (chars == 0) {
        chars = 4;
    }
    for (size_t i = 0; i < bytes; i++) {
        size_t offset = 2 * chars;
        for (size_t j = 0; j < chars; j++) {
            offset -= 2;
            result = (result << 2) | recoding.data[(kmer[i] >> offset) & 3];
        }
        chars = 4;
    }

    return result;
}

gbwtgraph::Key64::value_type kff_recode_trivial(const uint8_t* kmer, size_t k, size_t bytes) {
    gbwtgraph::Key64::value_type result = 0;
    for (size_t i = 0; i < bytes; i++) {
        result = (result << 8) | kmer[i];
    }
    return result & sdsl::bits::lo_set[2 * k];
}

std::vector<gbwtgraph::Key64::value_type> kff_recode(const uint8_t* kmers, size_t n, size_t k, kff_recoding_t recoding) {
    std::vector<gbwtgraph::Key64::value_type> result;
    result.reserve(n);

    size_t total_chars = n + k - 1;
    size_t bytes = kff_bytes(total_chars);
    size_t chars = total_chars & 3;
    if (chars == 0) {
        chars = 4;
    }

    gbwtgraph::Key64::value_type curr = 0;
    for (size_t i = 0, processed = 0; i < bytes; i++) {
        size_t offset = 2 * chars;
        for (size_t j = 0; j < chars; j++) {
            offset -= 2;
            curr = (curr << 2) | recoding.data[(kmers[i] >> offset) & 3];
            processed++;
            if (processed >= k) {
                result.push_back(curr & sdsl::bits::lo_set[2 * k]);
            }
        }
        chars = 4;
    }

    return result;
}

//------------------------------------------------------------------------------

uint8_t kff_get(const uint8_t* kmer, size_t i) {
    size_t byte = i / 4;
    size_t offset = 3 - (i & 3);
    return (kmer[byte] >> (2 * offset)) & 3;
}

void kff_set(std::vector<uint8_t>& kmer, size_t i, uint8_t value) {
    size_t byte = i / 4;
    size_t offset = 3 - (i & 3);
    kmer[byte] |= value << (2 * offset);
}

std::vector<uint8_t> kff_reverse_complement(const uint8_t* kmer, size_t k, const uint8_t* encoding) {
    uint8_t complement[4];
    complement[encoding[0]] = encoding[3];
    complement[encoding[1]] = encoding[2];
    complement[encoding[2]] = encoding[1];
    complement[encoding[3]] = encoding[0];

    size_t offset = (4 - (k & 3)) & 3;
    std::vector<uint8_t> result(kff_bytes(k), 0);
    for (size_t i = 0; i < k; i++) {
        kff_set(result, 4 * result.size() - 1 - i, complement[kff_get(kmer, i + offset)]);
    }
    return result;
}

//------------------------------------------------------------------------------

ParallelKFFReader::ParallelKFFReader(const std::string& filename) :
    file(new Kff_file(filename, "r"))
{
    // The variables are defined in the sections before the first block.
    this->find_blocks();
    auto get_var = [&](const std::string& name) -> size_t {
        auto iter = this->file->global_vars.find(name);
        if (iter == this->file->global_vars.end()) {
            throw std::runtime_error("ParallelKFFReader: file " + filename + " does not define variable " + name);
        }
        return iter->second;
    };

    this->k = get_var("k");
    if (this->k > gbwtgraph::Key64::KMER_MAX_LENGTH) {
        throw std::runtime_error("ParallelKFFReader: file " + filename + " contains " + std::to_string(this->k) +
            "-mers; cannot use k > " + std::to_string(gbwtgraph::Key64::KMER_MAX_LENGTH));
    }

    this->max_kmers_per_block = get_var("max");
    this->data_bytes = get_var("data_size");

    for (size_t i = 0; i < 4; i++) {
        this->encoding[i] = this->file->encoding[i];
    }
    this->recoding = kff_recoding(this->encoding);
}

std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> ParallelKFFReader::read(size_t n) {
    std::vector<std::pair<kmer_type, size_t>> result;
    result.reserve(n);

    // Fixed-size records to decode after releasing the file.
    std::vector<std::uint8_t> records;
    size_t kmer_length = 0, count_bytes = 0;
    {
        std::lock_guard<std::mutex> lock(this->mtx);

        while (!this->buffer.empty() && result.size() < n) {
            result.push_back(this->buffer.front());
            this->buffer.pop_front();
        }
        if (result.size() >= n || !this->find_blocks()) {
            return result;
        }

        // A block of a raw section stores the number of kmers in the block only if
        // there can be more than one, so blocks with one kmer are records of a fixed size.
        if (this->section->max == 1 && dynamic_cast<Section_Raw*>(this->section.get()) != nullptr) {
            kmer_length = this->section->k;
            count_bytes = this->section->data_size;
            size_t count = std::min(n - result.size(), static_cast<size_t>(this->section->remaining_blocks));
            records.resize(count * (kff_bytes(kmer_length) + count_bytes));
            this->file->read(records.data(), records.size());
            this->section->remaining_blocks -= count;
        } else {
            this->read_blocks(n, result);
        }
    }

    this->decode_records(records, kmer_length, count_bytes, result);
    return result;
}

bool ParallelKFFReader::find_blocks() {
    while (this->section == nullptr || this->section->remaining_blocks == 0) {
        this->section.reset();
        if (this->file->tellp() >= this->file->end_position) {
            return false;
        }
        char section_type = this->file->read_section_type();
        if (section_type == 'v') {
            // A variable section replaces the variables of the file.
            Section_GV variables(this->file.get());
        } else if (section_type == 'i') {
            Section_Index index(this->file.get());
            index.close();
        } else {
            this->section.reset(Block_section_reader::construct_section(this->file.get()));
            if (this->section == nullptr) {
                throw std::runtime_error(std::string("ParallelKFFReader: unknown section type ") + section_type);
            }
        }
    }
    return true;
}

void ParallelKFFReader::read_blocks(size_t n, std::vector<std::pair<kmer_type, size_t>>& result) {
    // Because we read kmers by blocks, we have to preallocate the buffers.
    size_t kmer_length = this->section->k;
    size_t count_bytes = this->section->data_size;
    std::vector<uint8_t> block(kff_bytes(this->section->max + kmer_length - 1));
    std::vector<uint8_t> data(this->section->max * count_bytes);
    size_t kmer_bytes = kff_bytes(kmer_length);
    bool trivial_encoding(kff_is_trivial(this->encoding));
    while (this->section->remaining_blocks > 0 && result.size() < n) {
        size_t block_size = this->section->read_compacted_sequence(block.data(), data.data());
        if (block_size > 1) {
            std::vector<kmer_type> kmers = kff_recode(block.data(), block_size, kmer_length, this->recoding);
            for (size_t i = 0; i < block_size; i++) {
                std::pair<kmer_type, size_t> kmer(kmers[i], kff_parse(data.data() + i * count_bytes, count_bytes));
                if (result.size() < n) {
                    result.push_back(kmer);
                } else {
                    this->buffer.push_back(kmer);
                }
            }
        } else {
            kmer_type kmer = (trivial_encoding ? kff_recode_trivial(block.data(), kmer_length, kmer_bytes) : kff_recode(block.data(), kmer_length, this->recoding));
            result.push_back({ kmer, kff_parse(data.data(), count_bytes) });
        }
    }
}

void ParallelKFFReader::decode_records(const std::vector<std::uint8_t>& records, size_t kmer_length, size_t count_bytes,
                                       std::vector<std::pair<kmer_type, size_t>>& result) const {
    size_t kmer_bytes = kff_bytes(kmer_length);
    size_t record_bytes = kmer_bytes + count_bytes;
    bool trivial_encoding(kff_is_trivial(this->encoding));
    for (size_t offset = 0; offset < records.size(); offset += record_bytes) {
        const uint8_t* record = records.data() + offset;
        kmer_type kmer = (trivial_encoding ? kff_recode_trivial(record, kmer_length, kmer_bytes) : kff_recode(record, kmer_length, this->recoding));
        result.push_back({ kmer, kff_parse(record + kmer_bytes, count_bytes) });
    }
}

//------------------------------------------------------------------------------

} // namespace vg