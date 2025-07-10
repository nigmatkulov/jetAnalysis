/**
 * @file JetUncertainty.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Jet uncertainty class. Original version from Yi Chen
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 * // JetUncertainty
 * // v1.0
 * // Author: Yi Chen
 * // 
 * // This class gives you jet uncertainties
 * //
 * 
 */

#ifndef JetUncertainty_h
#define JetUncertainty_h

// С++ headers
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// ROOT headers
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

//________________
class JetUncertainty {
  public:
    JetUncertainty()                  { Initialized = false; }
    JetUncertainty(std::string File)  { Initialized = false; Initialize(File); }
    virtual ~JetUncertainty()         {}
    void SetJetPT(double value)       { JetPT = value; }
    void SetJetEta(double value)      { JetEta = value; }
    void SetJetPhi(double value)      { JetPhi = value; }
    void SetJetArea(double value)     { JetArea = value; }
    void SetRho(double value)         { Rho = value; }
    void Initialize(std::string FileName);
    std::vector<std::string> BreakIntoParts(std::string Line);
    bool CheckDefinition(std::string Line);
    std::string StripBracket(std::string Line);
    enum Type { TypeNone, TypeJetPT, TypeJetEta, TypeJetPhi, TypeJetArea, TypeRho };
    JetUncertainty::Type ToType(std::string Line);
    std::pair<double, double> GetUncertainty();
    double GetValue(Type T);
  private:
    bool Initialized;
    double JetPT, JetEta, JetPhi, JetArea, Rho;
    std::vector<std::vector<Type>> BinTypes;
    std::vector<std::vector<double>> BinRanges;
    std::vector<std::vector<double>> PTBins;
    std::vector<std::vector<double>> ErrorLow;
    std::vector<std::vector<double>> ErrorHigh;

    ClassDef(JetUncertainty, 0)
};

#endif // #define JetUncertainty_h
