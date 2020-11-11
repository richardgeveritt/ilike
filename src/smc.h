#include <RcppArmadillo.h>
using namespace Rcpp;

#ifndef SMC_H
#define SMC_H

class SMC
{
public:

  SMC(void);
  SMC(const SMC &another);
  virtual ~SMC(void);

  void operator=(const SMC &another);
  virtual SMC* duplicate() const=0;

  List do_smc(const List &model,
              const List &algorithm);

protected:

  void make_copy(const SMC &another);
};

#endif
