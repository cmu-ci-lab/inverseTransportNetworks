/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/stan.h>
#include <mitsuba/core/thread.h>
#include <iostream>
#if defined(MTS_OPENMP)
# include <omp.h>
#else
# error "badness"
#endif

MTS_NAMESPACE_BEGIN


StanStack::StanStack()
{
}


StanStack::~StanStack()
{
}


stan::math::ChainableStack *StanStack::getInstance(int tid) {
    if ( m_instance[tid] == NULL ) m_instance[tid] = new stan::math::ChainableStack();
    return m_instance[tid];
}


void StanStack::resetInstance(int tid) {
    if ( m_instance[tid] ) {
        delete m_instance[tid];
        m_instance[tid] = NULL;
    }
}


void StanStack::record(int tid) {
	m_resetCounters[tid] = __MITSUBA_STAN_RESET_COUNTER;
    if ( m_instance[tid] ) {
        m_stackSizes[6*tid    ] = static_cast<int64_t>(m_instance[tid]->var_stack_.size());
        m_stackSizes[6*tid + 1] = static_cast<int64_t>(m_instance[tid]->var_nochain_stack_.size());
        m_stackSizes[6*tid + 2] = static_cast<int64_t>(m_instance[tid]->var_alloc_stack_.size());

        m_stackSizes[6*tid + 3] = static_cast<int64_t>(m_instance[tid]->nested_var_stack_sizes_.size());
        m_stackSizes[6*tid + 4] = static_cast<int64_t>(m_instance[tid]->nested_var_nochain_stack_sizes_.size());
        m_stackSizes[6*tid + 5] = static_cast<int64_t>(m_instance[tid]->nested_var_alloc_stack_starts_.size());
		
		m_instance[tid]->memalloc_.record();
    }
}


void StanStack::revert(int tid) {
	if ( --m_resetCounters[tid] == 0 ) {
		m_resetCounters[tid] = __MITSUBA_STAN_RESET_COUNTER;
		if ( m_stackSizes[6*tid] >= 0 ) {
			if ( m_instance[tid] == NULL ) m_instance[tid] = new stan::math::ChainableStack();
			m_instance[tid]->var_stack_.resize(m_stackSizes[6*tid]);
			m_instance[tid]->var_nochain_stack_.resize(m_stackSizes[6*tid + 1]);
			m_instance[tid]->var_alloc_stack_.resize(m_stackSizes[6*tid + 2]);
			m_instance[tid]->nested_var_stack_sizes_.resize(m_stackSizes[6*tid + 3]);
			m_instance[tid]->nested_var_nochain_stack_sizes_.resize(m_stackSizes[6*tid + 4]);
			m_instance[tid]->nested_var_alloc_stack_starts_.resize(m_stackSizes[6*tid + 5]);
			
			m_instance[tid]->memalloc_.revert();
		}
	}
}


std::vector<stan::math::ChainableStack *> StanStack::m_instance(NUM_THREADS, NULL);
std::vector<int64_t> StanStack::m_stackSizes(6*NUM_THREADS, -1);
std::vector<int> StanStack::m_resetCounters(NUM_THREADS, __MITSUBA_STAN_RESET_COUNTER);


MTS_NAMESPACE_END
