#include "electronResponse.hh"
#include "electronQuench.hh"
#include "electronCerenkov.hh"

double electronResponse::getElecNonl(double Etrue)
{
    return electronQuench::ScintillatorNL(Etrue) + electronCerenkov::getCerenkovPE(Etrue);
}
