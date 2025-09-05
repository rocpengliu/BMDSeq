#include <stdio.h>
#include "unittest.h"
#include <time.h>
#include "cmdline.h"
#include "bmdutil.h"
#include "options.h"
#include "common.h"
#include <zlib.h>
#include <unordered_map>
#include <set>
#include "kmer.h"
#include <vector>
#include "writer.h"
#include <sstream>
#include <mutex>
#include <thread>
#include <algorithm>
#include <queue>
#include "bmdseeker.h"

string command;
mutex logmtx;
mutex kmermtxl;
mutex kmermtxr;
using namespace std;
int main(int argc, char* argv[]){
    // display version info if no argument is given
    if(argc == 1) {
        cerr << "kmer2abunt: compile sample kmer count table into a sample kmer table" << endl << "version " << SEQ2KMER_VER << endl;
    }
    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cerr << "kmer2abunt " << SEQ2KMER_VER << endl;
        return 0;
    }
    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in_dir", 'i', "kmer result input dir name", false, "");
    cmd.add<string>("meta_table", 'm', "meta data table", false, "");
    cmd.add<int>("kmer_size", 'k', "kmer size must be between 5 and 64, default: 31", false, 31);
    cmd.add<int>("kmer_filter", 0, "minimum number of kmer with default 0", false, 0);
    cmd.add<string>("prefix", 'f', "prefix of the output file", false, "");
    cmd.add<int>("thread", 'w', "number of threads", false, 2);
    cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4.", false, 4);

    cmd.parse_check(argc, argv);
    if(argc == 1) {
        cerr << cmd.usage() <<endl;
    }
    if(argc == 1) {
        //output citation information
        cerr << "Citation:" <<endl;
        cerr << "Ultrafast all-in-one kmer counter" << endl;
        cerr << endl;
        return 0;
    }
    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();
    time_t t1 = time(NULL);
    Options* opt = new Options();
    opt->in_dir = cmd.get<string>("in_dir");
    opt->metatable = cmd.get<string>("meta_table");
    opt->kmerSize = cmd.get<int>("kmer_size");
    opt->kmerFilter = cmd.get<int>("kmer_filter");
    opt->prefix = cmd.get<string>("prefix");
    opt->thread = cmd.get<int>("thread");
    opt->compression = cmd.get<int>("compression");
    //opt->parseSample();
    opt->readMeta();
    cerr << endl << command << endl;

    std::string pattern = "_kmer_" + std::to_string(opt->kmerSize) + ".tab.gz";
    std::queue<std::string> samQueue;
    for(auto & it : opt->bmd.doseMap){
        for(auto itv : it.second.first){
            samQueue.push(joinpath(opt->in_dir, itv + pattern));
        }
    }

    int numThread = std::min(opt->thread, static_cast<int>(opt->bmd.samMap.size()));
    std::vector<std::thread> threads(numThread);

    for (int i = 0; i < numThread; ++i) {
        threads[i] = std::thread([&samQueue, &opt, &pattern, i]() {
            loginfo("Thread " + std::to_string(i) + " started");
            const int BUFFER_SIZE = 1024;
            char buffer[BUFFER_SIZE];
            uint64_t count = 0;
            uint32_t num = 0;
            std::vector<string> splitVec;
            std::string line = "";
            std::string kmer = "";

            while(true){
                std::unique_lock<std::mutex> lock(kmermtxl);
                if(samQueue.empty()){
                    lock.unlock();
                    break;
                }
                std::string sample = samQueue.front();
                samQueue.pop();
                lock.unlock();
                count = 0;
                num = 0;
                uint128_t key128 = 0;
                uint64_t key = 0;
                std::string sam = remove(basename(sample), pattern);
                int samIndex = opt->bmd.samMap.find(sam)->second;
                loginfo("start to read sample: " + sam);
                gzFile file = gzopen(sample.c_str(), "rb");
                if (!file)
                    error_exit("Error: Failed to open file " + sample);
                while (gzgets(file, buffer, BUFFER_SIZE) != NULL){
                    line = buffer; // Output or process each line
                    if(line.empty())
                        continue;
                    splitVec.clear();
                    split(line, splitVec, "\t");
                    if(splitVec.size() != 2)
                        continue;
                    kmer.clear();
                    kmer = splitVec.at(0);
                    if (kmer == "kmer")
                        continue;
                    if (kmer.length() != opt->kmerSize)
                        error_exit("kmer length is not the kmer_size");
                    num = static_cast<uint32_t>(std::stoul(splitVec.at(1)));
                    if(num < opt->kmerFilter)
                        continue;
                    if (opt->kmerSize > 32){
                        key128 = Kmer::processSingleKmer128(kmer);
                        std::unique_lock<std::mutex> lock2(kmermtxr);
                        auto itm = opt->bmd.kmerMap128.find(key128);
                        if(itm == opt->bmd.kmerMap128.end()){
                            std::vector<uint32_t> v(opt->bmd.samMap.size(), 0);
                            v.at(samIndex) = num;
                            opt->bmd.kmerMap128[key128] = v;
                        } else {
                            opt->bmd.kmerMap128[key128].at(samIndex) = num;
                        }
                        lock2.unlock();
                    } else {
                        key = Kmer::processSingleKmer(kmer);
                        std::unique_lock<std::mutex> lock2(kmermtxr);
                        auto itm = opt->bmd.kmerMap.find(key);
                        if(itm == opt->bmd.kmerMap.end()){
                            std::vector<uint32_t> v(opt->bmd.samMap.size(), 0);
                            v.at(samIndex) = num;
                            opt->bmd.kmerMap[key] = v;
                        } else {
                            opt->bmd.kmerMap[key].at(samIndex) = num;
                        }
                        lock2.unlock();
                    }
                    ++count;
                    if(count % 1000000 == 0 ){
                        loginfo("read " + std::to_string(count / 1000000) + "M lines for " + sam + " file!");
                    }
                }
                gzclose(file);
                loginfo("finish to read sample: " + sam);
            }
        });
    }

    for(auto & it : threads){
        if(it.joinable()){
            it.join();
        }
    }
    loginfo("all sample readings are done!");
    Writer *kmerWriter = new Writer(opt, opt->prefix + "_" + std::to_string(opt->kmerSize) + "_abunt.tab.gz", opt->compression);
    std::string* strData = new std::string();
    int chunk_size = 0;
    uint64_t count = 0;
    loginfo("starting to write kmer table " + std::to_string((opt->kmerSize > 32 ? opt->bmd.kmerMap128.size() : opt->bmd.kmerMap.size())) + "!");
    *strData += "sample\t";
    auto it = opt->bmd.doseMap.begin();
    while(it != opt->bmd.doseMap.end()){
        auto next = std::next(it);
        for(int i = 0; i != it->second.first.size(); ++i){
            *strData += it->second.first.at(i);
            if(next == opt->bmd.doseMap.end()){
                if(i == (it->second.first.size() -1)){
                    *strData += "\n";
                } else {
                    *strData += "\t";
                }
            } else {
                *strData += "\t";
            }
        }
        ++it;
    }
    kmerWriter->write(strData->data(), strData->length());
    strData->clear();
    loginfo("starting to write each kmer !");
    numThread = opt->thread;
    threads.clear();
    threads.resize(numThread);
    for(int i = 0; i < numThread; ++i){
        threads[i] = std::thread([&opt, &strData, &kmerWriter, &chunk_size, &count, i](){
            loginfo("Thread " + std::to_string(i) + " started");
            BmdSeeker* bmdseeker = new BmdSeeker();
            std::ostringstream oss;
            uint32_t sum = 0;
            int start = 0;
            int len = 0;
            bool trend = false;
            std::vector<double> mean_vec;
            mean_vec.reserve(opt->bmd.numDose);
            std::vector<double> y_vec;
            y_vec.reserve(opt->bmd.samMap.size());
            while(true){
                y_vec.clear();
                if(opt->kmerSize > 32){
                    std::unique_lock<std::mutex> lock(kmermtxl);
                    if(opt->bmd.kmerMap128.empty()) break;
                    std::pair<uint128_t, std::vector<uint32_t>> front_pair = *opt->bmd.kmerMap128.begin();
                    opt->bmd.kmerMap128.erase(opt->bmd.kmerMap128.begin());
                    lock.unlock();
                    oss << Kmer::uint2seq128(front_pair.first, opt->kmerSize) << "\t";
                    for(int idx = 0; idx != front_pair.second.size(); ++idx){
                        oss << front_pair.second.at(idx) << (idx == (front_pair.second.size() -1 ) ? "\n" : "\t");
                        y_vec.emplace_back(static_cast<double>(front_pair.second.at(idx)));
                    }
                } else {
                    std::unique_lock<std::mutex> lock(kmermtxl);
                    if(opt->bmd.kmerMap.empty()) break;
                    loginfo("aaaaaaaaaaaaaaaaaaaaaaaaaa");
                    std::pair<uint32_t, std::vector<uint32_t>> front_pair =  *opt->bmd.kmerMap.begin();
                    opt->bmd.kmerMap.erase(opt->bmd.kmerMap.begin());
                    lock.unlock();
                    oss << Kmer::uint2seq(front_pair.first, opt->kmerSize) << "\t";
                    for(int idx = 0; idx != front_pair.second.size(); ++idx){
                        oss << front_pair.second.at(idx) << (idx == (front_pair.second.size() -1 ) ? "\n" : "\t");
                        y_vec.emplace_back(static_cast<double>(front_pair.second.at(idx)));
                    }
                }

                loginfo("bbbbbbbbbbbbbbbbbbbbbbbbbbbb");
                double tot_mean = get_mean(y_vec, 0, y_vec.size());
                if (tot_mean < static_cast<double>(opt->kmerFilter)) continue;
                mean_vec.clear();
                for(const auto & itm : opt->bmd.doseMap){
                    start = itm.second.second;
                    len = itm.second.first.size();
                    mean_vec.emplace_back(get_mean(y_vec, start, len));
                }
                loginfo("cccccccccccccccccccccccccccc");
                trend = (lm(opt->bmd.doseVec, y_vec).slope >= 0) ? true : false;
                
                loginfo("ddddddddddddddddddddddddddd");
                for(const auto & itv : opt->bmd.doseVec){
                    std::cout << itv << "\t";
                }
                std::cout <<"\n";
                for(const auto & itv : y_vec){
                    std::cout << itv << "\t";
                }
                std::cout << "\n";
                try{
                    bmdseeker->runPythonContAnalysis(opt->bmd.doseVec, y_vec, trend);
                    y_vec.clear();
                } catch(const std::exception& ex) {
                    std::cerr << "Exception in thread " << i << ": " << ex.what() << std::endl;
                } catch(...) {
                    std::cerr << "Unknown exception in thread " << i << std::endl;
                }
                
                loginfo("eeeeeeeeeeeeeeeeeeeeeeee");
                
                std::unique_lock<std::mutex> lock2(kmermtxr);
                *strData += oss.str();
                if(count % 1000 == 0){
                    loginfo("write " + std::to_string(count / 1000) + "K lines for kmer output file!");
                }
                if (chunk_size == 1000000){
                    kmerWriter->write(strData->data(), strData->length());
                    strData->clear();
                    chunk_size = 0;
                } else{
                    ++chunk_size;
                }
                ++count;
                lock2.unlock();

                oss.str("");
                oss.clear();
            }

            if(bmdseeker){
                delete bmdseeker;
                bmdseeker = NULL;
            }
        });
    }

    if(!strData->empty()){
        kmerWriter->write(strData->data(), strData->length());
        strData->clear();
    }

    if(strData){
        delete strData;
        strData = NULL;
    }

    if(kmerWriter){
        delete kmerWriter;
        kmerWriter = NULL;
    }
    if(opt)
        delete opt;
    time_t t2 = time(NULL);
    cerr << "kmer2abunt v" << SEQ2KMER_VER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}
