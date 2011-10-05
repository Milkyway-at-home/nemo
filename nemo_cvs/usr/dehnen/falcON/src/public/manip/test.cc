#include <public/defman.h>
#include <public/bodyfunc.h>
#include <sstream>
#include <stdio.h>

namespace falcON { namespace Manipulate {
  //  mutable std::string DS;
  class test : public manipulator {
  public:
    const char*name    () const { return "test"; }
    const char*describe() const {
      return "Just a test'";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::o; }
    fieldset provide() const { return fieldset::f; }
    fieldset change () const { return fieldset::f; }
    //--------------------------------------------------------------------------
    test(const char*p, const char*f) falcON_THROWING  {
	falcON_WarningN("Hi from test!\n");
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*S) const {
	falcON_WarningN("Hi from inside the manipulator!\n");
      /*const BodyFilter*BF = S->pointer<BodyFilter>("filter");
      if(BF && *BF) {
	// make sure flags are supported
	if(!S->have(fieldbit::f))
	  const_cast<snapshot*>(S)->add_field(fieldbit::f);
	// check if all data needed are supported
	if(!S->have_all(BF->need()))
	  falcON_THROW("set_subset::manipulate(): "
		       "filter needs '%s' but snapshot has only '%s'\n",
		       word(BF->need()), word(S->all_data()));
	// loop bodies and run them through the filter
	unsigned sub(0);
	LoopAllBodies(S,b)
	  if((*BF)(b)) {
	    b.into_subset();
	    ++sub;
	  } else
	    b.outof_subset();
	// create description of subset
	std::ostringstream ost;
	ost  << sub << " bodies with \""
	     << BF->expression() << '\"';
	if(BF->npar()) {
	  ost<< " where";
	  for(int n=0; n!=BF->npar(); ++n)
	    ost<<" #"<<n<<'='<<BF->param(n);
	}
	DS = ost.str();
	// put pointer to description into pointer bank
	S->set_pointer(&DS,"subset_description");
      } */
      return false;
    } 
  };
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN__ALT(falcON::Manipulate::test);
