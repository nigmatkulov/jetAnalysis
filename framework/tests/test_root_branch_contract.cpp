#include "jetanalysis/RootBranchContract.h"

#include <iostream>

template <typename F> void rejects(F check) {
    try { check(); }
    catch (const std::runtime_error&) { return; }
    throw std::runtime_error("Expected schema rejection");
}

int main() {
    try {
        TTree tree("events", "schema checks");
        Float_t pt = 0;
        Double_t wrongPt = 0;
        Int_t count = 2, otherCount = 2;
        UInt_t run = 0;
        ULong64_t event = (ULong64_t{1} << 63) + 17;
        Float_t jets[2] = {1, 2};
        tree.Branch("pt", &pt, "pt/F");
        tree.Branch("wrongPt", &wrongPt, "wrongPt/D");
        tree.Branch("count", &count, "count/I");
        tree.Branch("otherCount", &otherCount, "otherCount/I");
        tree.Branch("run", &run, "run/i");
        tree.Branch("event", &event, "event/l");
        tree.Branch("jets", jets, "jets[count]/F");
        tree.Branch("fixed", jets, "fixed[2]/F");
        tree.Branch("single", jets, "single[1]/F");
        tree.Fill();
        using namespace jetanalysis;
        requireLeaf<Float_t>(tree, "pt");
        requireLeaf<Float_t>(tree, "jets", "count");
        rejects([&] { requireLeaf<Float_t>(tree, "missing"); });
        rejects([&] { requireLeaf<Float_t>(tree, "wrongPt"); });
        rejects([&] { requireLeaf<Int_t>(tree, "run"); });
        rejects([&] { requireLeaf<Float_t>(tree, "jets"); });
        rejects([&] { requireLeaf<Float_t>(tree, "single"); });
        rejects([&] { requireLeaf<Float_t>(tree, "jets", "otherCount"); });
        rejects([&] { requireLeaf<Float_t>(tree, "fixed", "count"); });
        ULong64_t readEvent = 0;
        bindScalar(tree, "event", readEvent);
        tree.GetEntry(0);
        if (readEvent != event || static_cast<std::uint64_t>(readEvent) != event)
            throw std::runtime_error("Event ID lost precision");
        tree.ResetBranchAddresses();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
