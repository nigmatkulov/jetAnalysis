/**
 * @file JetCorrector.cc
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

#include "JetCorrector.h"

//_________________
JetCorrector::JetCorrector() : JEC{}, 
    JetPT{-999.}, JetEta{-999.}, JetPhi{-999.}, 
    JetArea{-999.}, Rho{-999.} {
    /* Empty */
}

//________________
JetCorrector::JetCorrector(std::string File) : JEC{}, 
    JetPT{-999.}, JetEta{-999.}, JetPhi{-999.}, 
    JetArea{-999.}, Rho{-999.} { 
    Initialize(File); 
}

//________________
JetCorrector::JetCorrector(std::vector<std::string> Files) : 
    JEC{}, JetPT{-999.}, JetEta{-999.}, JetPhi{-999.}, 
    JetArea{-999.}, Rho{-999.} {
    Initialize(Files); 
}

//________________
void JetCorrector::Initialize(std::string File) { 
    std::vector<std::string> X; 
    X.push_back(File); 
    Initialize(X); 
}

//_________________
void JetCorrector::Initialize(std::vector<std::string> Files) {
    JEC.clear();
    for (auto File : Files)
    {
        JEC.push_back(SingleJetCorrector(File));
    }
}

//_________________
double JetCorrector::GetCorrection() {
    double PT = GetCorrectedPT();
    if (PT < 0) {
        return -1;
    }
    return PT / JetPT;
}

//_________________
double JetCorrector::GetCorrectedPT() {
    double PT = JetPT;

    //std::cout << "GetCorrectedPT first: " << PT << std::endl;

    for (int i = 0; i < (int)JEC.size(); i++) {
        JEC[i].SetJetPT(PT);
        JEC[i].SetJetEta(JetEta);
        JEC[i].SetJetPhi(JetPhi);
        JEC[i].SetRho(Rho);
        JEC[i].SetJetArea(JetArea);

        PT = JEC[i].GetCorrectedPT();

        //std::cout << "GetCorrectedPT mid: " << PT << std::endl;

        if (PT < 0) break;
    }

    //std::cout << "GetCorrectedPT last: " << PT << std::endl;
    return PT;
}

//________________
SingleJetCorrector::SingleJetCorrector() : Initialized{false}, IsFunction{false},
    JetPT{-999.}, JetEta{-999.}, JetPhi{-999.}, JetArea{-999.}, Rho{-999.},
    Formulas{}, Parameters{}, BinTypes{}, BinRanges{}, 
    Dependencies{}, DependencyRanges{}, Functions{} { 
    /* Empty */
}

//________________
SingleJetCorrector::SingleJetCorrector(std::string File) : 
    Initialized{false}, IsFunction{false}, 
    JetPT{-999.}, JetEta{-999.}, JetPhi{-999.}, JetArea{-999.}, Rho{-999.},
    Formulas{}, Parameters{}, BinTypes{}, BinRanges{}, 
    Dependencies{}, DependencyRanges{}, Functions{} { 

    Initialize(File); 
}

//________________
SingleJetCorrector::~SingleJetCorrector() { 
    for (auto P : Functions) {
        if( P ) {
            delete P;
            P = nullptr;
        }
    }
}

//_________________
void SingleJetCorrector::Initialize(std::string FileName) {
    int nvar = 0, npar = 0;
    std::string CurrentFormula = "";
    std::vector<Type> CurrentDependencies;
    std::vector<Type> CurrentBinTypes;

    std::ifstream in(FileName.c_str());

    //std::cout << "SingleJetCorrector::Initialize" << std::endl;
    //std::cout << "Input filename: " << FileName.c_str() << std::endl;

    while (in) {
        char ch[1048576];
        ch[0] = '\0';
        in.getline(ch, 1048577, '\n');
        //std::cout << "line: " << ch << std::endl;

        if (ch[0] == '\0') {
            //std::cout << "Empty string" << std::endl;
            continue;
        }

        bool IsDefinition = CheckDefinition(ch);

        //std::cout << "IsDefinition: " << IsDefinition << std::endl;
        std::vector<std::string> Parts = BreakIntoParts(StripBracket(ch));
        //std::cout << "Parts size: " << Parts.size() << std::endl;
        if (Parts.size() == 0) continue;

        if (IsDefinition == true) {
            // Found a definition line - update current formula

            //std::cout << "IsDefinition = true" << std::endl;

            nvar = atoi(Parts[0].c_str());
            if ( (int)Parts.size() <= (nvar + 1)) continue;

            npar = atoi(Parts[nvar + 1].c_str());
            if ( (int)Parts.size() <= (nvar + 1 + npar + 1)) continue;

            CurrentFormula = Parts[nvar + 1 + npar + 1];

            CurrentBinTypes.clear();
            for (int i = 0; i < nvar; i++) {
                CurrentBinTypes.push_back(ToType(Parts[1 + i]));
            }
                

            CurrentDependencies.clear();
            for (int i = 0; i < npar; i++) {
                CurrentDependencies.push_back(ToType(Parts[nvar + 2 + i]));
            }
        }
        else {
            // Otherwise it's a line with actual JECs, add it to the list

            //std::cout << "IsDefinition = false" << std::endl;

            if ((int)Parts.size() < (nvar * 2 + npar * 2 + 1)) continue;

            std::vector<double> Parameter;
            for (int i = nvar * 2 + npar * 2 + 1; i < (int)Parts.size(); i++) {
                Parameter.push_back(atof(Parts[i].c_str()));
            }
            Parameters.push_back(Parameter);

            Dependencies.push_back(CurrentDependencies);

            if (CurrentDependencies.size() == 4) {
                Formulas.push_back(Hack4(CurrentFormula, 't', Parameter.size()));
            }
            else {
                Formulas.push_back(CurrentFormula);
            }

            std::vector<double> Ranges;
            for (int i = nvar * 2 + 1; i < nvar * 2 + 1 + npar * 2; i++) {
                Ranges.push_back(atof(Parts[i].c_str()));
            }
                
            for (int i = 0; i + 1 < (int)Ranges.size(); i = i + 2) {
                if (Ranges[i] > Ranges[i + 1]) {
                    std::swap(Ranges[i], Ranges[i + 1]);
                }
            }  
            DependencyRanges.push_back(Ranges);
            BinTypes.push_back(CurrentBinTypes);

            Ranges.clear();
            for (int i = 0; i < nvar * 2; i++) {
                Ranges.push_back(atof(Parts[i].c_str()));
            }
            for (int i = 0; i + 1 < (int)Ranges.size(); i = i + 2) {
                if (Ranges[i] > Ranges[i + 1]) {
                    std::swap(Ranges[i], Ranges[i + 1]);
                }
            }
            BinRanges.push_back(Ranges);
            Functions.push_back(nullptr);
        }
    }

    //std::cout << "Closing the input file" << std::endl;
    in.close();
    Initialized = true;
}

//_________________
std::vector<std::string> SingleJetCorrector::BreakIntoParts(std::string Line) {
    std::stringstream str(Line);

    std::vector<std::string> Result;

    while (str) {
        std::string Temp = "";
        str >> Temp;
        if (Temp == "")
            continue;
        Result.push_back(Temp);
    }

    return Result;
}

//_________________
bool SingleJetCorrector::CheckDefinition(std::string Line) {
    for (int i = 0; i < (int)Line.size(); i++) {
        if (Line[i] == '{')
            return true;
        if (Line[i] == ' ')
            continue;
        return false;
    }

    return false;
}

//_________________
std::string SingleJetCorrector::StripBracket(std::string Line) {
    for (int i = 0; i < (int)Line.size(); i++) {
        if (Line[i] == '{' || Line[i] == '}') {
            Line.erase(Line.begin() + i);
            i = i - 1;
        }
    }

    return Line;
}

//_________________
SingleJetCorrector::Type SingleJetCorrector::ToType(std::string Line) {
    if (Line == "JetPt")
        return TypeJetPT;
    if (Line == "JetEta")
        return TypeJetEta;
    if (Line == "JetPhi")
        return TypeJetPhi;
    if (Line == "JetA")
        return TypeJetArea;
    if (Line == "Rho")
        return TypeRho;

    std::cerr << "[SingleJetCorrector] Warning: variable type " << Line << " not found!" << std::endl;

    return TypeNone;
}

//_________________
double SingleJetCorrector::GetCorrection() {
    if (Initialized == false) {
        //std::cout << "Forgot to initialize" << std::endl;
        return -1;
    }

    //std::cout << "Inside GetCorrection" << std::endl;

    int N = Formulas.size();

    //std::cout << "Formulas size: " << N << std::endl;

    for (int iE = 0; iE < N; iE++) {
        //std::cout << "iE: " << iE << std::endl;

        bool InBin = true;

        for (int iB = 0; iB < (int)BinTypes[iE].size(); iB++) {
            //std::cout << "iB: " << iB << std::endl;
            double Value = GetValue(BinTypes[iE][iB]);
            //std::cout << "Value: " << Value << std::endl;
            if (Value < BinRanges[iE][iB * 2] || Value > BinRanges[iE][iB * 2 + 1]) {
                InBin = false;
                //std::cout << "InBin = FALSE" << std::endl;
            }
        }

        if (InBin == false)
            continue;

        if (Dependencies[iE].size() == 0) {
            //std::cout << "Dependencies size: 0" << std::endl; 
            return -1; // huh?
        }
            
        if (Dependencies[iE].size() > 4) {
            std::cerr << "[SingleJetCorrector] There are " << Dependencies[iE].size() << " parameters!" << std::endl;
            return -1; // huh?
        }

        double V[3] = {0, 0, 0};
        for (int i = 0; i < 3; i++) {

            if ((int)Dependencies[iE].size() <= i) continue;

            double Value = GetValue(Dependencies[iE][i]);
            // std::cout << "Value: " << Value << std::endl;
            if (Value < DependencyRanges[iE][i * 2]) {
                Value = DependencyRanges[iE][i * 2];
            }
            if (Value > DependencyRanges[iE][i * 2 + 1]){
                Value = DependencyRanges[iE][i * 2 + 1];
            }
            V[i] = Value;
        }

        TF1 *Function = nullptr;

        if (Functions[iE] == nullptr) {
            if (Dependencies[iE].size() == 1)
                Function = new TF1(Form("Function%d", iE), (Formulas[iE] + "+0*x").c_str());
            if (Dependencies[iE].size() == 2)
                Function = new TF2(Form("Function%d", iE), (Formulas[iE] + "+0*x+0*y").c_str());
            if (Dependencies[iE].size() == 3)
                Function = new TF3(Form("Function%d", iE), (Formulas[iE] + "+0*x+0*y+0*z").c_str());
            if (Dependencies[iE].size() == 4)
                Function = new TF3(Form("Function%d", iE), (Formulas[iE] + "+0*x+0*y+0*z").c_str());

            Functions[iE] = Function;
        }
        else
            Function = Functions[iE];

        for (int i = 0; i < (int)Parameters[iE].size(); i++)
            Function->SetParameter(i, Parameters[iE][i]);
        if (Dependencies[iE].size() == 4)
            Function->SetParameter(Parameters[iE].size(), GetValue(Dependencies[iE][3]));
        double Result = Function->EvalPar(V);


        // std::cout << Formulas[iE] << std::endl;
        // std::cout << "P" << std::endl;
        // for(int i = 0; i < (int)Parameters[iE].size(); i++) {
        //     std::cout << " " << Parameters[iE][i] << std::endl;
        // }
        // std::cout << "V" << std::endl;
        // std::cout << " " << V[0] << std::endl;
        // std::cout << " " << V[1] << std::endl;
        // std::cout << " " << V[2] << std::endl;
        // std::cout << Dependencies[iE].size() << std::endl;
        // std::cout << Function->EvalPar(V) << std::endl;

        return Result;
    } // for (int iE = 0; iE < N; iE++)

    return -1;
}

//_________________
double SingleJetCorrector::GetCorrectedPT() {
    double Correction = GetCorrection();

    //std::cout << "SingleJetCorrector::GetCorrectedPT(): Correction: " << Correction << std::endl;
    if (Correction < 0) {
        return -1.;
    }
    
    //std::cout << "Correction - JetPT*Correction: " << JetPT * Correction << std::endl;

    return JetPT * Correction;
}

//_________________
double SingleJetCorrector::GetValue(Type T) {
    if (T == TypeNone)
        return 0;
    if (T == TypeJetPT)
        return JetPT;
    if (T == TypeJetEta)
        return JetEta;
    if (T == TypeJetPhi)
        return JetPhi;
    if (T == TypeJetArea)
        return JetArea;
    if (T == TypeRho)
        return Rho;

    return -1;
}

//_________________
std::string SingleJetCorrector::Hack4(std::string Formula, char V, int N) {
    int Size = Formula.size();
    for (int i = 0; i < Size; i++) {
        if (Formula[i] != V) continue;

        if (i != 0 && Formula[i - 1] >= 'a' && Formula[i - 1] <= 'z') continue;
        if (i != 0 && Formula[i - 1] >= 'A' && Formula[i - 1] <= 'Z') continue;
        if (i != Size - 1 && Formula[i + 1] >= 'a' && Formula[i + 1] <= 'z') continue;
        if (i != Size - 1 && Formula[i + 1] >= 'A' && Formula[i + 1] <= 'Z') continue;

        Formula.erase(i, 1);
        Formula.insert(i, Form("[%d]", N));
    }

    return Formula;
}
