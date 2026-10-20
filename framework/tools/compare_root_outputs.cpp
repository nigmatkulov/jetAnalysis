#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TNamed.h>
#include <TArrayD.h>
#include <array>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
void check(bool ok,const std::string& message) {if(!ok) throw std::runtime_error(message);}
void axis(const TAxis& a,const TAxis& b,const std::string& path) {
    check(a.GetNbins()==b.GetNbins(),path+": axis size");
    check(std::string(a.GetTitle())==b.GetTitle(),path+": axis title");
    check(a.GetXbins()->GetSize()==b.GetXbins()->GetSize(),path+": fixed/variable axis");
    for(int i=0;i<a.GetXbins()->GetSize();++i) check(a.GetXbins()->At(i)==b.GetXbins()->At(i),path+": variable edges");
    for(int i=0;i<=a.GetNbins()+1;++i) {
        check(a.GetBinLowEdge(i)==b.GetBinLowEdge(i),path+": bin edge");
        check(std::string(a.GetBinLabel(i))==b.GetBinLabel(i),path+": bin label");
    }
}
std::uint64_t objects=0,bins=0;
void object(TObject& a,TObject& b,const std::string& path) {
    check(std::string(a.ClassName())==b.ClassName(),path+": class");
    check(std::string(a.GetTitle())==b.GetTitle(),path+": title");
    if(auto* h=dynamic_cast<TH1*>(&a)) {
        auto& g=dynamic_cast<TH1&>(b);
        check(h->GetDimension()==g.GetDimension() && h->GetNcells()==g.GetNcells(),path+": dimensions");
        check(h->GetEntries()==g.GetEntries(),path+": entries");
        axis(*h->GetXaxis(),*g.GetXaxis(),path+"/x");
        axis(*h->GetYaxis(),*g.GetYaxis(),path+"/y");
        axis(*h->GetZaxis(),*g.GetZaxis(),path+"/z");
        check(h->GetSumw2N()==g.GetSumw2N(),path+": Sumw2 size");
        for(int i=0;i<h->GetSumw2N();++i) check(h->GetSumw2()->At(i)==g.GetSumw2()->At(i),path+": stored Sumw2");
        std::array<double,13> sa{},sb{}; h->GetStats(sa.data());g.GetStats(sb.data());
        check(sa==sb,path+": histogram statistics");
        for(int i=0;i<h->GetNcells();++i) {
            check(h->GetBinContent(i)==g.GetBinContent(i) && h->GetBinError(i)==g.GetBinError(i),path+": content/error at bin "+std::to_string(i)); ++bins;
        }
    } else if(auto* h=dynamic_cast<THnSparse*>(&a)) {
        auto& g=dynamic_cast<THnSparse&>(b);
        check(h->GetNdimensions()==g.GetNdimensions() && h->GetNbins()==g.GetNbins(),path+": sparse dimensions");
        check(h->GetEntries()==g.GetEntries() && h->GetCalculateErrors()==g.GetCalculateErrors(),path+": sparse entries/errors");
        check(h->GetSumw()==g.GetSumw() && h->GetSumw2()==g.GetSumw2(),path+": sparse statistics");
        for(int i=0;i<h->GetNdimensions();++i) {
            axis(*h->GetAxis(i),*g.GetAxis(i),path+"/axis"+std::to_string(i));
            check(h->GetSumwx(i)==g.GetSumwx(i) && h->GetSumwx2(i)==g.GetSumwx2(i),path+": sparse moments");
        }
        std::vector<int> ca(h->GetNdimensions()),cb(ca.size());
        for(Long64_t i=0;i<h->GetNbins();++i) {
            const auto va=h->GetBinContent(i,ca.data()), vb=g.GetBinContent(i,cb.data());
            check(ca==cb && va==vb && h->GetBinError2(i)==g.GetBinError2(i),path+": sparse coordinate/content/error at bin "+std::to_string(i)); ++bins;
        }
    } else check(a.IsA()==TNamed::Class(),path+": unsupported object type; refusing incomplete comparison");
    ++objects;
}
void directory(TDirectory& a,TDirectory& b,const std::string& path) {
    auto keys=[](TDirectory& d) {
        std::map<std::string,short> result;
        for(auto* item:*d.GetListOfKeys()) {
            auto* k=dynamic_cast<TKey*>(item);
            check(k!=nullptr,"Non-key directory entry");
            check(result.emplace(k->GetName(),k->GetCycle()).second,"Duplicate key cycles are unsupported");
        }
        return result;
    };
    const auto ka=keys(a),kb=keys(b); check(ka==kb,path+": key names/cycles");
    for(const auto& item:ka) {
        std::unique_ptr<TObject> x(a.GetKey(item.first.c_str())->ReadObj()), y(b.GetKey(item.first.c_str())->ReadObj());
        check(x && y,path+": cannot read key");
        const auto name=path+"/"+item.first;
        if(auto* d=dynamic_cast<TDirectory*>(x.get())) {
            auto* e=dynamic_cast<TDirectory*>(y.get());check(e,name+": directory type");directory(*d,*e,name);
        } else object(*x,*y,name);
    }
}
}
int main(int argc,char** argv) {
    if(argc!=3) {std::cerr<<"Usage: compare_root_outputs REFERENCE.root CANDIDATE.root\n";return 2;}
    try {
        TH1::AddDirectory(false);
        TFile a(argv[1],"READ"),b(argv[2],"READ");
        check(!a.IsZombie() && !b.IsZombie() && !a.TestBit(TFile::kRecovered) && !b.TestBit(TFile::kRecovered),"Invalid input ROOT file");
        check(a.GetNkeys()>0,"Empty reference file");directory(a,b,"");
        std::cout<<"EXACT MATCH: "<<objects<<" objects, "<<bins<<" histogram cells/sparse bins; no tolerances.\n";
    } catch(const std::exception& e) {std::cerr<<"MISMATCH: "<<e.what()<<'\n';return 1;}
}
