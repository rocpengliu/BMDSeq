#ifndef KMER_H  // Check if KMER_H is not defined
#define KMER_H  // Define KMER_H to prevent multiple inclusions

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <memory.h>
#include <unordered_map>
#include <vector>
#include "options.h"
#include "read.h"
#include "common.h"

using namespace std;

class Kmer {
public:
    Kmer();
    ~Kmer();
    static void processKmerMapPair(Read*& r, std::unordered_map<uint64_t, uint32_t>& kmerMap, int& kmerSize, bool noRevKmer, double& threshold, bool lowComplexFilter);
    static void processKmerMap(Read*& r, std::unordered_map<uint64_t, uint32_t>& kmerMap, int& kmerSize, double& threshold, bool lowComplexFilter);
    static uint64_t processSingleKmer(std::string& kmer);
    static uint128_t processSingleKmer128(std::string& kmer);
    static std::string uint2seq(uint64_t key, int &kmerSize);
    static std::string uint2seq128(uint128_t key, int& kmerSize);
    static bool passLowComplexityFilter(const std::string &str, double& threshold);
    static std::vector<std::string> splitNRead(string &str, int len, double& threshold, bool lowComplexFilter);
    static void processKmerMapPair(Read*& r, std::unordered_map<uint128_t, uint32_t, uint128_hash>& kmerMap, int& kmerSize, bool noRevKmer, double& threshold, bool lowComplexFilter);
    static void processKmerMap(Read*& r, std::unordered_map<uint128_t, uint32_t, uint128_hash>& kmerMap, int& kmerSize, double& threshold, bool lowComplexFilter);
    static uint8_t BASE2NUC(const char *base, bool norm = true);
    //static std::string print_uint128_binary(uint128_t value);
    static std::string print_uint_binary(uint128_t value, int& kmerSize);
    static std::string print_uint_binary(uint64_t value, int& kmerSize);
    static std::string uint128_to_string(__uint128_t value);
};

#endif