#ifndef STAN_MATH_REV_CORE_START_NESTED_HPP
#define STAN_MATH_REV_CORE_START_NESTED_HPP

//#include <stan/math/rev/core/chainablestack.hpp>
#include <mitsuba/core/stan.h>

namespace stan {
  namespace math {

    /**
     * Record the current position so that <code>recover_memory_nested()</code>
     * can find it.
     */
    static inline void start_nested() {
      mitsuba::StanStack::getInstance()->nested_var_stack_sizes_
        .push_back(mitsuba::StanStack::getInstance()->var_stack_.size());
      mitsuba::StanStack::getInstance()->nested_var_nochain_stack_sizes_
        .push_back(mitsuba::StanStack::getInstance()->var_nochain_stack_.size());
      mitsuba::StanStack::getInstance()->nested_var_alloc_stack_starts_
        .push_back(mitsuba::StanStack::getInstance()->var_alloc_stack_.size());
      mitsuba::StanStack::getInstance()->memalloc_.start_nested();
    }

  }
}
#endif
