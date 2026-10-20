#pragma once

#include <RtypesCore.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TTree.h>

#include <climits>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace jetanalysis {

// Input buffers use ROOT types. Transient event objects may use fixed-width
// integers only after these compile-time representation checks succeed.
static_assert(CHAR_BIT == 8, "Analysis requires eight-bit bytes");
static_assert(sizeof(Int_t) == 4 && sizeof(UInt_t) == 4,
              "Analysis requires 32-bit ROOT integers");
static_assert(sizeof(ULong64_t) == 8, "Analysis requires 64-bit event IDs");
static_assert(sizeof(Float_t) == 4 && std::numeric_limits<Float_t>::is_iec559,
              "Analysis requires IEEE 754 binary32 ROOT Float_t");
static_assert(sizeof(Double_t) == 8 && std::numeric_limits<Double_t>::is_iec559,
              "Analysis requires IEEE 754 binary64 ROOT Double_t");

template <typename T> struct RootLeafType;
template <> struct RootLeafType<Bool_t> { static const char* name() { return "Bool_t"; } };
template <> struct RootLeafType<Float_t> { static const char* name() { return "Float_t"; } };
template <> struct RootLeafType<Double_t> { static const char* name() { return "Double_t"; } };
template <> struct RootLeafType<Int_t> { static const char* name() { return "Int_t"; } };
template <> struct RootLeafType<UInt_t> { static const char* name() { return "UInt_t"; } };
template <> struct RootLeafType<ULong64_t> { static const char* name() { return "ULong64_t"; } };

// Call on each physical input TTree, not just the first tree in a TChain.
// Direct tree lookup deliberately avoids resolving a same-named friend branch.
template <typename T>
TLeaf& requireLeaf(TTree& tree, const char* branchName,
                   const char* countName = nullptr) {
    auto* branch = dynamic_cast<TBranch*>(tree.GetListOfBranches()->FindObject(branchName));
    auto* leaf = branch ? branch->GetLeaf(branchName) : nullptr;
    const std::string context = std::string(tree.GetName()) + "/" + branchName;
    if (!leaf || branch->GetListOfLeaves()->GetEntries() != 1 ||
        std::string(leaf->GetTypeName()) != RootLeafType<T>::name()) {
        throw std::runtime_error(context + ": expected single " + RootLeafType<T>::name() + " leaf");
    }
    const auto* count = leaf->GetLeafCount();
    if (!countName) {
        if (count || leaf->GetLenStatic() != 1 ||
            std::string(leaf->GetTitle()).find('[') != std::string::npos)
            throw std::runtime_error(context + ": expected scalar");
    } else {
        auto& expectedCount = requireLeaf<Int_t>(tree, countName);
        if (count != &expectedCount || leaf->GetLenStatic() != 1)
            throw std::runtime_error(context + ": incorrect array count or shape");
    }
    return *leaf;
}

template <typename T>
void bindScalar(TTree& tree, const char* name, T& buffer) {
    requireLeaf<T>(tree, name);
    tree.SetBranchStatus(name, true);
    // ROOT positive status codes include conversions. Require an exact match.
    if (tree.SetBranchAddress(name, &buffer) != TTree::kMatch)
        throw std::runtime_error(std::string("Cannot bind exact branch type: ") + name);
}

} // namespace jetanalysis
