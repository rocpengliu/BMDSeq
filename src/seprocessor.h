#ifndef SE_PROCESSOR_H
#define SE_PROCESSOR_H

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
#include "writerthread.h"
#include "duplicate.h"
#include "singleproducersingleconsumerlist.h"
#include "readpool.h"
#include "common.h"

extern mutex kmermtx;
using namespace std;

typedef struct ReadRepository ReadRepository;

class SingleEndProcessor{
public:
    SingleEndProcessor(Options* opt);
    ~SingleEndProcessor();
    bool process();
public:
    std::unordered_map<uint64_t, uint32_t> *kmerMap;
    std::unordered_map<uint128_t, uint32_t, uint128_hash> *kmerMap128;
    std::map<uint32_t, int> kmerCountFreq;
    //std::map<int, uint64_t> freqMap;

private:
    bool processSingleEnd(ReadPack* pack, ThreadConfig* config);
    void readerTask();
    void processorTask(ThreadConfig* config);
    void initConfig(ThreadConfig* config);
    void initOutput();
    void closeOutput();
    void writerTask(WriterThread* config);
    void recycleToPool(int tid, Read* r);
    void addKmerMap(std::unordered_map<uint64_t, uint32_t>& subKmerMap);
    void addKmerMap(std::unordered_map<uint128_t, uint32_t, uint128_hash>& subKmerMap128);
    void writeKmerOutput();

private:
    Options* mOptions;
    atomic_bool mReaderFinished;
    atomic_int mFinishedThreads;
    Filter* mFilter;
    UmiProcessor* mUmiProcessor;
    WriterThread* mLeftWriter;
    WriterThread* mFailedWriter;
    Duplicate* mDuplicate;
    SingleProducerSingleConsumerList<ReadPack*>** mInputLists;
    size_t mPackReadCounter;
    atomic_long mPackProcessedCounter;
    ReadPool* mReadPool;
};


#endif