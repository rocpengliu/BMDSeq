#ifndef COMMON_H
#define COMMON_H

#include <functional>
#include <cstdint>
#include <string>

using namespace std;

#define SEQ2KMER_VER "0.0.1"

#define _DEBUG false

#ifndef _WIN32
	typedef long int64;
	typedef unsigned long uint64;
#else
	typedef long long int64;
	typedef unsigned long long uint64;
#endif

typedef int int32;
typedef unsigned int uint32;

typedef short int16;
typedef unsigned short uint16;

typedef char int8;
typedef unsigned char uint8;

typedef __uint128_t uint128_t;

const char ATCG_BASES[] = {'A', 'T', 'C', 'G'};

struct uint128_hash{
	size_t operator()(const uint128_t& key) const{
		size_t high = static_cast<size_t>(key >> 64);
		size_t low = static_cast<size_t>(key);
		size_t hash = std::hash<uint64_t>{}(high);
        hash ^= std::hash<uint64_t>{}(low) + 0x9e3779b9 + (hash << 6) + (hash >> 2); // similar to boost::hash_combine
		return hash;
	}
};

#pragma pack(2)


#pragma pack()


// how many reads one pack has
static const int PACK_SIZE = 256;

// if one pack is produced, but not consumed, it will be kept in the memory
// this number limit the number of in memory packs
// if the number of in memory packs is full, the producer thread should sleep
static const int PACK_IN_MEM_LIMIT = 128;


// different filtering results, bigger number means worse
// if r1 and r2 are both failed, then the bigger one of the two results will be recorded
// we reserve some gaps for future types to be added
static const int PASS_FILTER = 0;
static const int FAIL_POLY_X = 4;
static const int FAIL_OVERLAP = 8;
static const int FAIL_N_BASE = 12;
static const int FAIL_LENGTH = 16;
static const int FAIL_TOO_LONG = 17;
static const int FAIL_QUALITY = 20;
static const int FAIL_COMPLEXITY = 24;

// how many types in total we support
static const int FILTER_RESULT_TYPES = 32;

const static char* FAILED_TYPES[FILTER_RESULT_TYPES] = {
	"passed", "", "", "",
	"failed_polyx_filter", "", "", "",
	"failed_bad_overlap", "", "", "",
	"failed_too_many_n_bases", "", "", "",
	"failed_too_short", "failed_too_long", "", "",
	"failed_quality_filter", "", "", "",
	"failed_low_complexity", "", "", "",
	"", "", "", ""
};

const std::string BMD_MODELS[] = {
	"Lin",
	"Hill",
	"Exp2",
	"Exp3",
	"Exp4",
	"Exp5",
	"Poly2",
	"Poly3",
	"Poly4",
	"Power"
};

#endif /* COMMON_H */
