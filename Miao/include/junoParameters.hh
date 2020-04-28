#ifndef junoParameters_h
#define junoParameters_h

#include <string>

using namespace std;

enum ScintillatorParameterization {
    kIntegral,
    kSimulation
};

class junoParameters
{
    public:
        static std::string stopPow_File;
        static std::string quenchNL_File;
        static std::string cerenkovNL_File;
        static std::string electronLSNL_File;
        static std::string electronOut_File;
    
        static std::string gammaLSNL_File;
        static std::string gammaPdf_File;
        static std::string gammaOut_File;

        static double m_gammaError;

        static ScintillatorParameterization scintillatorParameterization;


};

#endif
