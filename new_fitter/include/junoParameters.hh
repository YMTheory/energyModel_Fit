#ifndef junoParameters_h
#define junoParameters_h

#include <string>

using namespace std;

enum ScintillatorParameterization {
    kIntegral,
    kSimulation,
    kEmpirical
};

class junoParameters
{
    public:
        static std::string stopPow_File;
        static std::string quenchNL_File;
        static std::string cerenkovNL_File;
        static std::string electronLSNL_File;
        static std::string electronOut_File;
        static std::string quenchNL_outFile;
        static std::string cerenkov_outFile;
    
        static std::string gammaLSNL_File;
        static std::string gammaPdf_File;
        static std::string gammaOut_File;
        static std::string gammaNLOption;

        static double m_gammaError;
        static bool fitGammaSources;

        static ScintillatorParameterization scintillatorParameterization;

        static bool fitB12;
        static std::string B12DataFile;
        static std::string B12PredFile;
        static std::string B12CalcFile;
        static double b12FitMinE;
        static double b12FitMaxE;
        static double b12VertexCut;

        // energy resolution input/output file ...
        static std::string B12_File;
        static std::string B12Out_File;
        static std::string B12Spec_File;

        static bool fitC11;
        static std::string C11DataFile;
        static std::string C11_File;
        static std::string C11Out_File;
        static std::string C11Spec_File;
        static double c11FitMinE;
        static double c11FitMaxE;

        static bool fitC10;
        static std::string C10DataFile;
        static std::string C10Out_File;
        static double c10FitMinE;
        static double c10FitMaxE;

        static double m_gammaScale;


};

#endif
