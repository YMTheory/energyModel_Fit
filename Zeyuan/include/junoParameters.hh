#ifndef junoParameters_h
#define junoParameters_h

#include <string>

using namespace std;

class junoParameters
{
    public:
        
        static int nFitParameter;
        static double p0_start;
        static double p1_start;
        static double p2_start;

        static int fitPrintLevel;

        static std::string gammaData_file;
        static std::string gammaPdf_file;

        static bool fitGamma;
        static bool fixScintP0;
        static bool fixScintP1;
        static bool fixScintP2;
        static bool fixGamScale;

        static double gamScale;
        static double gamScaleError;

};

#endif
