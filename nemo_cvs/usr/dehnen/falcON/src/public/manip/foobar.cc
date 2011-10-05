#include <public/defman.h>
#include <public/tools.h>
#include <public/io.h>
#include <stdlib.h>

namespace {
  using namespace falcON;
    class foobar : public manipulator {
      public:
      const char* name() const;
      const char* describe const;
      fieldset need() const;
      fieldset change() const;
      fieldset provide() const;
   
      bool manipulate(const snapshot*) const;
      foobar(const double *pars, int npar, const char *file);
      ~foobar();
   };
}
__DEF__MAN(foobar)

