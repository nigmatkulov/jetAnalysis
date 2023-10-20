/**
 * @file JetCorrector.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Jet energy corrector class. Original version from Yi Chen
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 * // SingleJetCorrector
 * // v3.0
 * // Author: Yi Chen
 * // 
 * // This class applies JEC for any given level using TF1 as the workhorse
 * // Supposedly runs faster than v1.0
 * // v3.0: one can add list of text files to apply them one by one
 * 
 */

#ifndef JetCorrector_h
#define JetCorrector_h

// C++ headers
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// ROOT headers
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

// ROOT headers
#include "TObject.h"

//_________________
class SingleJetCorrector {

  public:
    SingleJetCorrector();
    SingleJetCorrector(std::string File);
    virtual ~SingleJetCorrector();
    enum Type { TypeNone, TypeJetPT, TypeJetEta, TypeJetPhi, TypeJetArea, TypeRho };
    void SetJetPT(double value)     { JetPT = value; }
    void SetJetEta(double value)    { JetEta = value; }
    void SetJetPhi(double value)    { JetPhi = value; }
    void SetJetArea(double value)   { JetArea = value; }
    void SetRho(double value)       { Rho = value; }
    void Initialize(std::string FileName);
    std::vector<std::string> BreakIntoParts(std::string Line);
    bool CheckDefinition(std::string Line);
    std::string StripBracket(std::string Line);
    SingleJetCorrector::Type ToType(std::string Line);
    double GetCorrection();
    double GetCorrectedPT();
    double GetValue(Type T);

  private:
    bool Initialized;
    bool IsFunction;
    double JetPT; 
    double JetEta; 
    double JetPhi;
    double JetArea; 
    double Rho;
    std::vector<std::string> Formulas;
    std::vector<std::vector<double>> Parameters;
    std::vector<std::vector<Type>> BinTypes;
    std::vector<std::vector<double>> BinRanges;
    std::vector<std::vector<Type>> Dependencies;
    std::vector<std::vector<double>> DependencyRanges;
    std::vector<TF1 *> Functions;
    std::string Hack4(std::string Formula, char V, int N);

    ClassDef(SingleJetCorrector, 0)
};

//_________________
class JetCorrector {
    
  public:
    JetCorrector();
    JetCorrector(std::string File);
    JetCorrector(std::vector<std::string> Files);
    virtual ~JetCorrector() { /* empty */ }
    void Initialize(std::string File);
    void Initialize(std::vector<std::string> Files);
    void SetJetPT(double value)     { JetPT = value; }
    void SetJetEta(double value)    { JetEta = value; }
    void SetJetPhi(double value)    { JetPhi = value; }
    void SetJetArea(double value)   { JetArea = value; }
    void SetRho(double value)       { Rho = value; }
    double GetCorrection();
    double GetCorrectedPT();

  private:
    std::vector<SingleJetCorrector> JEC;
    double JetPT;
    double JetEta; 
    double JetPhi; 
    double JetArea; 
    double Rho;

    ClassDef(JetCorrector, 0)
};

#endif // #define JetCorrector_h