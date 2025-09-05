#ifndef PE_PROCESSOR_H
#define PE_PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "read.h"
#include <cstdlib>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <unordered_map>
#include "options.h"
#include "threadconfig.h"
#include "filter.h"
#include "umiprocessor.h"
#include "overlapanalysis.h"
#include "writerthread.h"
#include "duplicate.h"
#include "readpool.h"
#include "common.h"

extern mutex kmermtx;
using namespace std;

typedef struct ReadPairRepository ReadPairRepository;

class PairEndProcessor{
public:
    PairEndProcessor(Options* opt);
    ~PairEndProcessor();
    bool process();
public:
    std::unordered_map<uint64_t, uint32_t> *kmerMap;
    std::unordered_map<uint128_t, uint32_t, uint128_hash> *kmerMap128;
    std::map<uint32_t, int> kmerCountFreq;
    //std::map<int, uint64_t> freqMap;

private:
    bool processPairEnd(ReadPack* leftPack, ReadPack* rightPack, ThreadConfig* config);
    void readerTask(bool isLeft);
    void interleavedReaderTask();
    void processorTask(ThreadConfig* config);
    void initConfig(ThreadConfig* config);
    void initOutput();
    void closeOutput();
    void statInsertSize(Read* r1, Read* r2, OverlapResult& ov, int frontTrimmed1 = 0, int frontTrimmed2 = 0);
    int getPeakInsertSize();
    void writerTask(WriterThread* config);
    void recycleToPool1(int tid, Read* r);
    void recycleToPool2(int tid, Read* r);
    void addKmerMap(std::unordered_map<uint64_t, uint32_t>& subKmerMap);
    void addKmerMap(std::unordered_map<uint128_t, uint32_t, uint128_hash>& subKmerMap);
    void writeKmerOutput();

private:
    atomic_bool mLeftReaderFinished;
    atomic_bool mRightReaderFinished;
    atomic_int mFinishedThreads;
    Options* mOptions;
    Filter* mFilter;
    UmiProcessor* mUmiProcessor;
    atomic_long* mInsertSizeHist;
    WriterThread* mLeftWriter;
    WriterThread* mRightWriter;
    WriterThread* mUnpairedLeftWriter;
    WriterThread* mUnpairedRightWriter;
    WriterThread* mMergedWriter;
    WriterThread* mFailedWriter;
    WriterThread* mOverlappedWriter;
    //WriterThread *mKmerOutWriter;
    Duplicate *mDuplicate;
    SingleProducerSingleConsumerList<ReadPack*>** mLeftInputLists;
    SingleProducerSingleConsumerList<ReadPack*>** mRightInputLists;
    size_t mLeftPackReadCounter;
    size_t mRightPackReadCounter;
    atomic_long mPackProcessedCounter;
    ReadPool* mLeftReadPool;
    ReadPool* mRightReadPool;
};


#endif