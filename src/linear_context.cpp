#include "suitesparse_wrapper/linear_context.h"

#ifdef HAVE_SUITESPARSE
namespace GC {
// global context for linear solvers
LinearContext context;

LinearContext::LinearContext(void)
// constructor
{
  cholmod_l_start(&context);
}

LinearContext::~LinearContext(void)
// destructor
{
  cholmod_l_finish(&context);
}

void LinearContext::setSimplicial(void) {
  context.supernodal = CHOLMOD_SIMPLICIAL;
}

void LinearContext::setSupernodal(void) {
  context.supernodal = CHOLMOD_SUPERNODAL;
}

LinearContext::operator cholmod_common*(void)
// allows LinearContext to be treated as a cholmod_common*
{
  return &context;
}
}
#endif
