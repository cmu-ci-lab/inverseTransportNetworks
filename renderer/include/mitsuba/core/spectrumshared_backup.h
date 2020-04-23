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

#ifndef SPECTRUMSHARED_H
#define SPECTRUMSHARED_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/FloatADshared.h>
#include <mitsuba/core/FloatAD.h>
#include <stan/math/rev/core.hpp>


class Spectrum_shared {
public:
	std::vector<FloatAD_shared> m_spec;

	double val() const {
		int tid = mitsuba::Thread::getID();			
		return m_spec[0].m_data[tid].val();
	}
	double adj() const {
		int tid = mitsuba::Thread::getID();			
		return m_spec[0].m_data[tid].adj();
	}	




	Spectrum_shared() : m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp();
			m_spec[i] = temp;
		}
	}

	Spectrum_shared(float x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}

	Spectrum_shared(double x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}

	Spectrum_shared(long double x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}

	Spectrum_shared(bool x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}


	Spectrum_shared(char x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}
	Spectrum_shared(short x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}

	Spectrum_shared(int x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}

	Spectrum_shared(long x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}

	Spectrum_shared(unsigned char x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}	

	Spectrum_shared(unsigned short x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}	

	Spectrum_shared(unsigned int x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}		
	Spectrum_shared(unsigned long x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}



#ifdef _WIN64

	Spectrum_shared(size_t x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}
	Spectrum_shared(ptrdiff_t x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}
	
#endif

#ifdef BOOST_MATH_USE_FLOAT128

	Spectrum_shared(__float128 x): m_spec(SPECTRUM_SAMPLES){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			m_spec[i] = temp;
		}
	}

#endif	
};


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, float y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] *= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, float y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] += y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, float y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] -= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, float y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] /= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}





inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, double y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] *= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, double y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] += y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, double y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] -= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, double y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] /= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, int y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] *= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, int y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] += y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, int y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] -= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, int y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] /= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, FloatAD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] *= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, FloatAD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] += y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, FloatAD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] -= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, FloatAD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] /= y;
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, FloatAD_shared y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] *= y.m_data[tid];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, FloatAD_shared y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] += y.m_data[tid];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, FloatAD_shared y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] -= y.m_data[tid];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, FloatAD_shared y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] /= y.m_data[tid];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] *= y[0];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] += y[0];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] -= y[0];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] /= y[0];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}



inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] *= y[0];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] += y[0];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] -= y[0];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_spec[0].m_data[tid] /= y[0];
	return mitsuba::Spectrum_AD(x.m_spec[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(FloatAD_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_data[tid] += y[0];
	return mitsuba::Spectrum_AD(x.m_data[tid]);
}

												
inline mitsuba::Spectrum_AD operator-(FloatAD_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_data[tid] -= y[0];
	return mitsuba::Spectrum_AD(x.m_data[tid]);
}
inline mitsuba::Spectrum_AD operator*(FloatAD_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_data[tid] *= y[0];
	return mitsuba::Spectrum_AD(x.m_data[tid]);
}
inline mitsuba::Spectrum_AD operator/(FloatAD_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_data[tid] /= y[0];
	return mitsuba::Spectrum_AD(x.m_data[tid]);
}



inline mitsuba::Spectrum_AD operator+(FloatAD_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_data[tid] += y[0];
	return mitsuba::Spectrum_AD(x.m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(FloatAD_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_data[tid] -= y[0];
	return mitsuba::Spectrum_AD(x.m_data[tid]);
}
inline mitsuba::Spectrum_AD operator*(FloatAD_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_data[tid] *= y[0];
	return mitsuba::Spectrum_AD(x.m_data[tid]);
}
inline mitsuba::Spectrum_AD operator/(FloatAD_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x.m_data[tid] /= y[0];
	return mitsuba::Spectrum_AD(x.m_data[tid]);
}



#endif 