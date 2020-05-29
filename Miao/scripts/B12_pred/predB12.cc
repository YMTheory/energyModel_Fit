#include "fermiFunc.hh"
#include "finite_size.hh"
#include "screening.hh"
#include "weak_magnetism.hh"
#include "B12prediction.hh"

using namespace std;

int main()
{
    B12prediction* pred = new B12prediction();
    pred->Plot();

    return 1;
}
