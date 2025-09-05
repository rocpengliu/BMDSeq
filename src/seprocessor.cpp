#include "seprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "adaptertrimmer.h"
#include "polyx.h"
#include "kmer.h"

SingleEndProcessor::SingleEndProcessor(Options* opt){
    mOptions = opt;
    mReaderFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mUmiProcessor = new UmiProcessor(opt);
    mLeftWriter =  NULL;
    mFailedWriter = NULL;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }

    mPackReadCounter = 0;
    mPackProcessedCounter = 0;

    mReadPool = new ReadPool(mOptions);
    kmerMap = new std::unordered_map<uint64_t, uint32_t>();
    kmerMap128 = new std::unordered_map<uint128_t, uint32_t, uint128_hash>();
    //freqMap.clear();
    kmerCountFreq.clear();
}

SingleEndProcessor::~SingleEndProcessor() {
    delete mFilter;
    if(mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
    if(mReadPool) {
        delete mReadPool;
        mReadPool = NULL;
    }
    delete[] mInputLists;
    if(kmerMap){
        delete kmerMap;
    }

    if(kmerMap128){
        delete kmerMap128;
    }
}

void SingleEndProcessor::initOutput() {
    if(!mOptions->failedOut.empty())
        mFailedWriter = new WriterThread(mOptions, mOptions->failedOut);
    if(mOptions->out1.empty())
        return;
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
}

void SingleEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if(mFailedWriter) {
        delete mFailedWriter;
        mFailedWriter = NULL;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;

    if(mOptions->split.enabled) {
        config->initWriterForSplit();
    }
}

bool SingleEndProcessor::process(){
    if(!mOptions->split.enabled)
        initOutput();

    mInputLists = new SingleProducerSingleConsumerList<ReadPack*>*[mOptions->thread];

    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        mInputLists[t] = new SingleProducerSingleConsumerList<ReadPack*>();
        configs[t] = new ThreadConfig(mOptions, t, false);
        configs[t]->setInputList(mInputLists[t]);
        initConfig(configs[t]);
    }

    std::thread readerThread(std::bind(&SingleEndProcessor::readerTask, this));

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::processorTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* failedWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writerTask, this, mLeftWriter));
    if(mFailedWriter)
        failedWriterThread = new std::thread(std::bind(&SingleEndProcessor::writerTask, this, mFailedWriter));

    readerThread.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    std::thread* kmerWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeKmerOutput, this));

    if(!mOptions->split.enabled) {
        if(leftWriterThread)
            leftWriterThread->join();
        if(failedWriterThread)
            failedWriterThread->join();
    }

    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and read filter results
    vector<Stats*> preStats;
    vector<Stats*> postStats;
    vector<FilterResult*> filterResults;
    for(int t=0; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats* finalPreStats = Stats::merge(preStats);
    Stats* finalPostStats = Stats::merge(postStats);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);

    // read filter results to the first thread's
    for(int t=1; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
    }

    cerr << "Read1 before filtering:"<<endl;
    finalPreStats->print();
    cerr << endl;
    cerr << "Read1 after filtering:"<<endl;
    finalPostStats->print();

    cerr << endl;
    cerr << "Filtering result:"<<endl;
    finalFilterResult->print();

    double dupRate = 0.0;
    if(mOptions->duplicate.enabled) {
        dupRate = mDuplicate->getDupRate();
        cerr << endl;
        cerr << "Duplication rate (may be overestimated since this is SE data): " << dupRate * 100.0 << "%" << endl;
    }

    if(kmerWriterThread->joinable()){
        kmerWriterThread->join();
    }

    // make JSON report
    JsonReporter jr(mOptions);
    jr.setDup(dupRate);
    jr.report(finalFilterResult, finalPreStats, finalPostStats);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDup(dupRate);
    if(mOptions->kmerSize > 32){
        hr.setKmer(kmerCountFreq, kmerMap128->size());
    } else {
        hr.setKmer(kmerCountFreq, kmerMap->size());
    }
    hr.report(finalFilterResult, finalPreStats, finalPostStats);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats;
    delete finalPostStats;
    delete finalFilterResult;

    delete[] threads;
    delete[] configs;

    if(leftWriterThread)
        delete leftWriterThread;
    if(failedWriterThread)
        delete failedWriterThread;

    if(!mOptions->split.enabled)
        closeOutput();

    return true;
}

void SingleEndProcessor::recycleToPool(int tid, Read* r) {
    // failed to recycle, then delete it
    if(!mReadPool->input(tid, r))
        delete r;
}

bool SingleEndProcessor::processSingleEnd(ReadPack* pack, ThreadConfig* config){
    string* outstr = new string();
    string* failedOut = new string();
    std::unordered_map<uint64_t, uint32_t> subKmerMap;
    std::unordered_map<uint128_t, uint32_t, uint128_hash> subKmerMap128;
    int tid = config->getThreadId();

    int readPassed = 0;
    for(int p=0;p<pack->count;p++){

        // original read1
        Read* or1 = pack->data[p];

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);

        // handling the duplication profiling
        bool dedupOut = false;
        if(mDuplicate) {
            bool isDup = mDuplicate->checkRead(or1);
            if(mOptions->duplicate.dedup && isDup)
                dedupOut = true;
        }

        // filter by index
        if(mOptions->indexFilter.enabled && mFilter->filterByIndex(or1)) {
            recycleToPool(tid, or1);
            continue;
        }

        // fix MGI
        if(mOptions->fixMGI) {
            or1->fixMGI();
        }
        
        // umi processing
        if(mOptions->umi.enabled)
            mUmiProcessor->process(or1);

        int frontTrimmed = 0;
        // trim in head and tail, and apply quality cut in sliding window
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed);

        if(r1 != NULL) {
            if(mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, config->getFilterResult(), mOptions->polyGTrim.minLen);
        }

        if(r1 != NULL && mOptions->adapter.enabled){
            bool trimmed = false;
            if(mOptions->adapter.hasSeqR1)
                trimmed = AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
            bool incTrimmedCounter = !trimmed;
            if(mOptions->adapter.hasFasta) {
                AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, incTrimmedCounter);
            }
        }

        if(r1 != NULL) {
            if(mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }

        if(r1 != NULL) {
            if( mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
        }

        int result = mFilter->passFilter(r1);

        config->addFilterResult(result, 1);

        if(!dedupOut) {
            if( r1 != NULL &&  result == PASS_FILTER) {
                if(r1->length() >= mOptions->kmerSize){
                    if(mOptions->kmerSize > 32){
                        Kmer::processKmerMapPair(r1, subKmerMap128, mOptions->kmerSize, mOptions->NoRevKmer, mOptions->complexityFilter.threshold, mOptions->complexityFilter.enabled);
                    } else {
                        Kmer::processKmerMapPair(r1, subKmerMap, mOptions->kmerSize, mOptions->NoRevKmer, mOptions->complexityFilter.threshold, mOptions->complexityFilter.enabled);
                    }
                }
                r1->appendToString(outstr);

                // stats the read after filtering
                config->getPostStats1()->statRead(r1);
                readPassed++;
            } else if(mFailedWriter) {
                or1->appendToStringWithTag(failedOut, FAILED_TYPES[result]);
            }
        }

        recycleToPool(tid, or1);
        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL)
            recycleToPool(tid, r1);
    }

    if(mOptions->kmerSize > 32){
        if(!subKmerMap128.empty()){
            addKmerMap(subKmerMap128);
        }
    } else {
        if(!subKmerMap.empty()){
            addKmerMap(subKmerMap);
        }
    }

    if(mOptions->outputToSTDOUT) {
        fwrite(outstr->c_str(), 1, outstr->length(), stdout);
    } else if(mOptions->split.enabled) {
        // split output by each worker thread
        if(!mOptions->out1.empty())
            config->getWriter1()->writeString(outstr);
    }

    if(mLeftWriter) {
        mLeftWriter->input(tid, outstr);
        outstr = NULL;
    }
    if(mFailedWriter) {
        // write failed data
        mFailedWriter->input(tid, failedOut);
        failedOut = NULL;
    }

    if(mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(pack->count);

    if(outstr)
        delete outstr;
    if(failedOut)
        delete failedOut;

    delete pack->data;
    delete pack;

    mPackProcessedCounter++;

    return true;
}

void SingleEndProcessor::readerTask()
{
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    Read** data = new Read*[PACK_SIZE];
    memset(data, 0, sizeof(Read*)*PACK_SIZE);
    FastqReader reader(mOptions->in1, true, mOptions->phred64);
    reader.setReadPool(mReadPool);
    int count=0;
    bool needToBreak = false;
    while(true){
        Read* read = reader.read();
        if(!read || needToBreak){
            // the last pack
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            mInputLists[mPackReadCounter % mOptions->thread]->produce(pack);
            mPackReadCounter++;
            data = NULL;
            if(read) {
                delete read;
                read = NULL;
            }
            break;
        }
        data[count] = read;
        count++;
        // configured to process only first N reads
        if(mOptions->readsToProcess >0 && count + readNum >= mOptions->readsToProcess) {
            needToBreak = true;
        }
        if(mOptions->verbose && count + readNum >= lastReported + 1000000) {
            lastReported = count + readNum;
            string msg = "loaded " + to_string((lastReported/1000000)) + "M reads";
            loginfo(msg);
        }
        // a full pack
        if(count == PACK_SIZE || needToBreak){
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            mInputLists[mPackReadCounter % mOptions->thread]->produce(pack);
            mPackReadCounter++;
            //re-initialize data for next pack
            data = new Read*[PACK_SIZE];
            memset(data, 0, sizeof(Read*)*PACK_SIZE);
            // if the processor is far behind this reader, sleep and wait to limit memory usage
            while( mPackReadCounter - mPackProcessedCounter > PACK_IN_MEM_LIMIT){
                //cerr<<"sleep"<<endl;
                slept++;
                usleep(100);
            }
            readNum += count;
            // if the writer threads are far behind this reader, sleep and wait
            // check this only when necessary
            if(readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
                while(mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT) {
                    slept++;
                    usleep(1000);
                }
            }
            // reset count to 0
            count = 0;
            // re-evaluate split size
            // TODO: following codes are commented since it may cause threading related conflicts in some systems
            /*if(mOptions->split.needEvaluation && !splitSizeReEvaluated && readNum >= mOptions->split.size) {
                splitSizeReEvaluated = true;
                // greater than the initial evaluation
                if(readNum >= 1024*1024) {
                    size_t bytesRead;
                    size_t bytesTotal;
                    reader.getBytes(bytesRead, bytesTotal);
                    mOptions->split.size *=  (double)bytesTotal / ((double)bytesRead * (double) mOptions->split.number);
                    if(mOptions->split.size <= 0)
                        mOptions->split.size = 1;
                }
            }*/
        }
    }

    for(int t=0; t<mOptions->thread; t++)
        mInputLists[t]->setProducerFinished();

    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mReaderFinished = true;
    if(mOptions->verbose) {
        loginfo("Loading completed with " + to_string(mPackReadCounter) + " packs");
    }
    //lock.unlock();

    // if the last data initialized is not used, free it
    if(data != NULL)
        delete[] data;
}

void SingleEndProcessor::processorTask(ThreadConfig* config)
{
    SingleProducerSingleConsumerList<ReadPack*>* input = config->getLeftInput();
    while(true) {
        if(config->canBeStopped()){
            break;
        }
        while(input->canBeConsumed()) {
            ReadPack* data = input->consume();
            processSingleEnd(data, config);
        }
        if(input->isProducerFinished()) {
            if(!input->canBeConsumed()) {
                if(mOptions->verbose) {
                    string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                    loginfo(msg);
                }
                break;
            }
        } else {
            usleep(100);
        }
    }
    input->setConsumerFinished();        

    mFinishedThreads++;
    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
        if(mFailedWriter)
            mFailedWriter->setInputCompleted();
    }

    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void SingleEndProcessor::writerTask(WriterThread* config)
{
    while(true) {
        if(config->isCompleted()){
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if(mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}

void SingleEndProcessor::addKmerMap(std::unordered_map<uint64_t, uint32_t>& subKmerMap){
    kmermtx.lock();
    for(auto & it : subKmerMap){
        (*kmerMap)[it.first] += it.second;
        if(kmerMap->size() % 1000000 == 0){
            loginfo("processed " + std::to_string(kmerMap->size() / 1000000) + " M kmers!");
        }
    }
    kmermtx.unlock();
    subKmerMap.clear();
}

void SingleEndProcessor::addKmerMap(std::unordered_map<uint128_t, uint32_t, uint128_hash>& subKmerMap128){
    kmermtx.lock();
    for(auto & it : subKmerMap128){
        (*kmerMap128)[it.first] += it.second;
        if(kmerMap128->size() % 1000000 == 0){
            loginfo("processed " + std::to_string(kmerMap128->size() / 1000000 ) + " M kmers!");
        }
    }
    kmermtx.unlock();
    subKmerMap128.clear();
}

void SingleEndProcessor::writeKmerOutput(){
    Writer* kmerWriter = new Writer(mOptions, mOptions->kmerOut, mOptions->compression);
    std::ostringstream oss;
    std::string* strData = new std::string();
    uint64_t count = 0;
    int chunk_size = 0;
    oss << "kmer" << "\t" << "count" << "\n";
    if(mOptions->kmerSize > 32){
        for(const auto & it : *kmerMap128){
            kmerCountFreq[it.second]++;
            if (it.second < mOptions->kmerFilter)
                continue;
            oss << Kmer::uint2seq(it.first, mOptions->kmerSize) << "\t" << it.second << "\n";
            if(count % 5000000 == 0 && count > 0){
                loginfo("write " + std::to_string(count / 1000000) + "0 M lines for kmer output file!");
            }
            if (chunk_size == 1000000){
                strData->clear();
                *strData = oss.str();
                kmerWriter->write(strData->data(), strData->length());
                oss.str("");
                oss.clear();
                chunk_size = 0;
            }
            ++chunk_size;
            ++count;
        }
    } else {
        for(const auto & it : *kmerMap){
            kmerCountFreq[it.second]++;
            if(it.second < mOptions->kmerFilter)
                continue;
            oss << Kmer::uint2seq(it.first, mOptions->kmerSize) << "\t" << it.second << "\n";
            if(count % 5000000 == 0 && count > 0){
                loginfo("write " + std::to_string(count / 1000000) + "0 M lines for kmer output file!");
            }
            if (chunk_size == 1000000){
                strData->clear();
                *strData = oss.str();
                kmerWriter->write(strData->data(), strData->length());
                oss.str("");
                oss.clear();
                chunk_size = 0;
            }
            ++chunk_size;
            ++count;
        }
    }
    if(chunk_size > 0){
        strData->clear();
        *strData = oss.str();
        kmerWriter->write(strData->data(), strData->length());
        oss.str("");
        oss.clear();
    }
    if(kmerWriter){
        delete kmerWriter;
        kmerWriter = NULL;
    }
    loginfo(std::to_string(count) + " out of " + std::to_string(mOptions->kmerSize > 32 ? kmerMap128->size() : kmerMap->size()) + 
    "("+ std::to_string(getPer(count, mOptions->kmerSize > 32 ? kmerMap128->size() : kmerMap->size())) + "%) kmers written done!\n");
}
