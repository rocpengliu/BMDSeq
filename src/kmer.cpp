#include "kmer.h"
#include <sstream>
#include "util.h"

Kmer::Kmer(){
}
Kmer::~Kmer(){

}

void Kmer::processKmerMapPair(Read* & r, std::unordered_map<uint64_t, uint32_t>& kmerMap, int& kmerSize, bool noRevKmer, double& threshold, bool lowComplexFilter){
    processKmerMap(r, kmerMap, kmerSize, threshold, lowComplexFilter);
    if(!noRevKmer){
        Read *rr = r->reverseComplement();
        if(rr){
            processKmerMap(rr, kmerMap, kmerSize, threshold, lowComplexFilter);
            delete rr;
            rr = NULL;
        }
    }
}

void Kmer::processKmerMapPair(Read* & r, std::unordered_map<uint128_t, uint32_t, uint128_hash>& kmerMap, int& kmerSize, bool noRevKmer, double& threshold, bool lowComplexFilter){
    processKmerMap(r, kmerMap, kmerSize, threshold, lowComplexFilter);
    if(!noRevKmer){
        Read *rr = r->reverseComplement();
        if(rr){
            processKmerMap(rr, kmerMap, kmerSize, threshold, lowComplexFilter);
            delete rr;
            rr = NULL;
        }
    }
}

void Kmer::processKmerMap(Read*& r, std::unordered_map<uint64_t, uint32_t>& kmerMap, int& kmerSize, double& threshold, bool lowComplexFilter){
    if(r->length() < kmerSize)
        return;
    std::string dna = *(r->mSeq);
    std::vector<std::string> dna_vec = splitNRead(dna, kmerSize, threshold, lowComplexFilter);
    std::string kmer_str(kmerSize, 'N');
    uint64_t mask = ((uint64_t)1 << (2 * kmerSize)) -1;
    for(const auto & it : dna_vec){
        uint64_t kmer_int = 0;
        const char* subdna = it.c_str();
        for (int i = 0; i < it.length(); ++i){
            kmer_int = (kmer_int <<= 2) | BASE2NUC(&subdna[i]);
            if(i > (kmerSize - 2)){
                if(kmerSize < 32){
                    kmer_int &= mask;
                }
                if(lowComplexFilter){
                    kmer_str = it.substr((i - kmerSize + 1), kmerSize);
                    if(passLowComplexityFilter(kmer_str, threshold)){
                        kmerMap[kmer_int]++;
                    }
                } else {
                    kmerMap[kmer_int]++;
                }
            }
        }
    }
}

void Kmer::processKmerMap(Read*& r, std::unordered_map<uint128_t, uint32_t, uint128_hash>& kmerMap, int& kmerSize, double& threshold, bool lowComplexFilter){
    if(r->length() < kmerSize)
        return;
    std::string dna = *(r->mSeq);
    std::vector<std::string> dna_vec = splitNRead(dna, kmerSize, threshold, lowComplexFilter);
    std::string kmer_str(kmerSize, 'N');
    uint128_t mask = ((uint128_t)1 << (2 * kmerSize)) - 1;
    for(const auto & it : dna_vec){
        uint128_t kmer_int = 0;
        const char* subdna = it.c_str();
        for (int i = 0; i < it.length(); ++i){
            kmer_int = (kmer_int <<= 2) | BASE2NUC(&subdna[i]);
            if(i > (kmerSize - 2)){
                if(kmerSize < 64){
                    kmer_int &= mask;
                }
                if(lowComplexFilter){
                    kmer_str = it.substr((i - kmerSize + 1), kmerSize);
                    if(passLowComplexityFilter(kmer_str, threshold)){
                        kmerMap[kmer_int]++;
                    }
                } else {
                    kmerMap[kmer_int]++;
                }
            }
        }
    }
}

uint64_t Kmer::processSingleKmer(std::string& kmer){
    uint64_t kmer_int = 0;
    const char *dna = kmer.c_str();
    for (int i = 0; i < kmer.length(); ++i){
        kmer_int = (kmer_int << 2) | BASE2NUC(&dna[i]);
    }
    return kmer_int;
}
uint128_t Kmer::processSingleKmer128(std::string& kmer){
    uint128_t kmer_int = 0;
    const char *dna = kmer.c_str();
    for (int i = 0; i < kmer.length(); ++i){
        kmer_int = (kmer_int << 2) | BASE2NUC(&dna[i]);
    }
    return kmer_int;
}

std::vector<std::string> Kmer::splitNRead(string& str, int len, double& threshold, bool lowComplexFilter){
    std::vector<std::string> ret_;
    if (str.empty()){
        return ret_;
    }

    while(!str.empty() && str[0] == 'N'){
        str.erase(0, 1);
    }

    if (str.empty()){
        return ret_;
    }

    while(!str.empty() && str.back() == 'N'){
        str.pop_back();
    }

    if (str.empty()){
        return ret_;
    }

    string tmp;
    string::size_type pos_begin = str.find_first_not_of('N');
    string::size_type n_pos = 0;
    int n_n = 0;
    while (pos_begin != string::npos){
        n_pos = str.find('N', pos_begin);
        if (n_pos != string::npos){
            ++n_n;
            tmp = str.substr(pos_begin, n_pos - pos_begin);
            pos_begin = n_pos + 1;
        } else {
            tmp = str.substr(pos_begin);
            pos_begin = n_pos;
        }
        if(tmp.length() >= len){
            if(lowComplexFilter){
                if(passLowComplexityFilter(tmp, threshold)){
                    ret_.push_back(tmp);
                }
            } else {
                ret_.push_back(tmp);
            }
        }
        tmp.clear();
    }
    if(ret_.empty() && n_n == 0){
        ret_.push_back(str);
    }
    return ret_;
}

bool Kmer::passLowComplexityFilter(const std::string& str, double& threshold) {
    int diff = 0;
    if(str.length() <= 1)
        return false;
    const char* data = str.c_str();
    for(int i=0; i< str.length()-1; i++) {
        if(data[i] != data[i+1])
            diff++;
    }

    if( (double)diff/(double)(str.length()-1) >= threshold ){
        return true;
    } else {
        return false;
    }
}

std::string Kmer::uint2seq(uint64_t key, int& kmerSize){
    std::string dna(kmerSize, 'N');
    for (int i = 0; i != kmerSize; ++i){
        uint8_t nbits = key & 0b11;
        dna[kmerSize -1 - i] = ATCG_BASES[nbits];
        key >>= 2;
    }
    return dna;
}

std::string Kmer::uint2seq128(uint128_t key, int& kmerSize){
    std::string dna(kmerSize, 'N');
    for (int i = 0; i != kmerSize; ++i){
        uint8_t nbits = key & 0b11;
        dna[kmerSize -1 - i] = ATCG_BASES[nbits];
        key >>= 2;
    }
    return dna;
}

std::string Kmer::uint128_to_string(__uint128_t value) {
    if (value == 0) return "0";
    std::string result;
    while (value > 0) {
        result.insert(result.begin(), '0' + value % 10);
        value /= 10;
    }
    return result;
}

uint8_t Kmer::BASE2NUC(const char* base, bool norm){//must be coupled by ATCG_BASES in common.h
    switch(*base){
        case 'A':
            return (norm ? 0b00 : 0b01);
        case 'T':
            return (norm ? 0b01 : 0b00);
        case 'C':
            return (norm ? 0b10 : 0b11);
        case 'G':
            return (norm ? 0b11 : 0b10);
        default:
            cerr << "invalide base is: " << *base << endl;
            throw std::invalid_argument("Invalid nucleotide");
    }
}

std::string Kmer::print_uint_binary(uint128_t value, int& kmerSize) {
    std::string str(kmerSize*2, '0');
    for(int i = 0; i != 2*kmerSize; ++i){
        uint128_t bit = (value >> i) & 1;
        str[2*kmerSize - 1 - i] = static_cast<char>('0' + bit);
    }
    return str;
}

std::string Kmer::print_uint_binary(uint64_t value, int& kmerSize) {
    std::string str(kmerSize*2, '0');
    for(int i = 0; i != 2*kmerSize; ++i){
        uint64_t bit = (value >> i) & 1;
        str[2*kmerSize - 1 - i] = static_cast<char>('0' + bit);
    }
    return str;
}