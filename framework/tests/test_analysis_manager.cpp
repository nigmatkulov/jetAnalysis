#include "jetanalysis/AnalysisManager.h"
#include <stdexcept>
#include <string>
#include <vector>
using namespace jetanalysis;
namespace {
void check(bool ok,const char* message) {if(!ok) throw std::runtime_error(message);}
template<class F> void rejects(F f) {try {f();} catch(const std::exception&) {return;} throw std::runtime_error("Expected failure");}
struct Reader final : EventReader {
    int next=0,limit=3;
    std::int64_t entries() const override {return 3;}
    bool readNext(Event& e) override {
        if(next==limit)return false;
        e.info.evt=(std::uint64_t{1}<<63)+next++;
        e.info.vz=2.f;
        return true;
    }
};
struct Probe final : AnalysisTask {
    std::string name;
    std::vector<std::string>& log;
    int count=0; bool fail;
    Probe(std::string n,std::vector<std::string>& l,bool f=false):name(std::move(n)),log(l),fail(f) {}
    void initialize(const RunContext& c) override {check(c.entries==3,"Wrong normalization count");log.push_back(name+" init");}
    void process(const Event& e) override {
        check(e.info.evt==(std::uint64_t{1}<<63)+count,"Task received different event");
        check(e.info.vz==2.f,"Input changed between tasks");
        log.push_back(name+" "+std::to_string(count++));
        if(fail)throw std::runtime_error("task failure");
    }
    void finish() override {log.push_back(name+" finish");}
};
}
int main() {
    Reader reader;std::vector<std::string> log;
    AnalysisManager manager(reader);
    rejects([&]{manager.addTask(nullptr);});
    manager.addTask(std::make_unique<Probe>("a",log));
    manager.addTask(std::make_unique<Probe>("b",log));
    const auto summary=manager.run();
    check(summary.processedEvents==3,"Wrong processed count");
    check(log==std::vector<std::string>{"a init","b init","a 0","b 0","a 1","b 1","a 2","b 2","a finish","b finish"},"Wrong task lifecycle/order");
    rejects([&]{manager.run();});
    rejects([&]{manager.addTask(std::make_unique<Probe>("late",log));});
    Reader shortReader;shortReader.limit=1;std::vector<std::string> partial;
    AnalysisManager shortRun(shortReader);shortRun.addTask(std::make_unique<Probe>("x",partial));
    rejects([&]{shortRun.run();});check(partial.back()=="x 0","Failed reader called finish");
    Reader failureReader;std::vector<std::string> failure;
    AnalysisManager failed(failureReader);failed.addTask(std::make_unique<Probe>("x",failure,true));
    failed.addTask(std::make_unique<Probe>("y",failure));
    rejects([&]{failed.run();});check(failure.back()=="x 0","Failed task dispatched later tasks or called finish");
    Reader emptyReader;AnalysisManager empty(emptyReader);rejects([&]{empty.run();});
}
