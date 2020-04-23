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

#pragma once
#if !defined(__MITSUBA_CORE_STAN_H_)
#define __MITSUBA_CORE_STAN_H_

#include <vector>
#include <stan/math/rev/core/chainablestack.hpp>
#include <mitsuba/core/thread.h>

#if !defined(__MITSUBA_STAN_RESET_COUNTER)
#define __MITSUBA_STAN_RESET_COUNTER 1
#endif


MTS_NAMESPACE_BEGIN

class MTS_EXPORT_CORE StanStack {
public:
	static stan::math::ChainableStack *getInstance(int tid = Thread::getID());

    static void resetInstance(int tid = Thread::getID()); // HARD reset

    static void record(int tid = Thread::getID());
    static void revert(int tid = Thread::getID());

protected:
	StanStack();
	~StanStack();

private:
	static std::vector<stan::math::ChainableStack *> m_instance;
    static std::vector<int64_t> m_stackSizes;
	static std::vector<int> m_resetCounters;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_STAN_H_ */
