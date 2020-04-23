#ifndef FLOATADSHARED_H
#define FLOATADSHARED_H

#include <mitsuba/mitsuba.h>
#include <stan/math/rev/core.hpp>
#include <mitsuba/core/FloatAD.h>
#include <ostream>
#include <vector>


class FloatAD_shared {
public:


	/**
	 * Construct a variable for later assignment.
	 *
	 * This is implemented as a no-op, leaving the underlying implementation
	 * dangling.  Before an assignment, the behavior is thus undefined just
	 * as for a basic double.
	 */

	FloatAD_shared() : m_data(NUM_THREADS) {

		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp();
			m_data[i] = temp;
		}

	}

	FloatAD_shared(float x): m_data(NUM_THREADS){
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}

	FloatAD_shared(double x): m_data(NUM_THREADS) {
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}



	FloatAD_shared(long double x): m_data(NUM_THREADS) {
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}

	FloatAD_shared(bool x) : m_data(NUM_THREADS){
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}


	FloatAD_shared(char x): m_data(NUM_THREADS) {
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}

	FloatAD_shared(short x): m_data(NUM_THREADS) {
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}

	FloatAD_shared(int x) : m_data(NUM_THREADS){
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}

	FloatAD_shared(long x): m_data(NUM_THREADS) {
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}
	FloatAD_shared(unsigned char x) : m_data(NUM_THREADS){
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}

	FloatAD_shared(unsigned short x): m_data(NUM_THREADS) {
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}

	
	FloatAD_shared(unsigned int x) : m_data(NUM_THREADS){
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}
	FloatAD_shared(unsigned long x) : m_data(NUM_THREADS){
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}	


#ifdef _WIN64
	FloatAD_shared(size_t x) : m_data(NUM_THREADS){

		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}
	FloatAD_shared(ptrdiff_t x): m_data(NUM_THREADS) {
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}
	
#endif

#ifdef BOOST_MATH_USE_FLOAT128
	FloatAD_shared(__float128 x) : m_data(NUM_THREADS){
		for (int i = 0; i < NUM_THREADS; i++) {
			stan::math::var temp(x,i);
			m_data[i] = temp;
		}
	}

#endif

	double val() const {
		int tid = mitsuba::Thread::getID();			
		return m_data[tid].val();
	}
	double adj() const {
		int tid = mitsuba::Thread::getID();			
		return m_data[tid].adj();
	}	


	FloatAD& operator+=(double b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] += b;
		return m_data[tid];
	};

	FloatAD& operator*=(double b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] *= b;
		return m_data[tid];
	};

	FloatAD& operator-=(double b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] -= b;
		return m_data[tid];
	};

	FloatAD& operator/=(double b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] /= b;
		return m_data[tid];
	};


	
	FloatAD& operator+=(FloatAD b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] += b;
		return m_data[tid];
	};

	FloatAD& operator*=(FloatAD b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] *= b;
		return m_data[tid];
	};

	FloatAD& operator-=(FloatAD b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] -= b;
		return m_data[tid];
	};

	FloatAD& operator/=(FloatAD b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] /= b;
		return m_data[tid];
	};


	
/////////////////////////////////////
	FloatAD& operator+=(float b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] += b;
		return m_data[tid];
	};

	FloatAD& operator*=(float b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] *= b;
		return m_data[tid];
	};

	FloatAD& operator-=(float b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] -= b;
		return m_data[tid];
	};

	FloatAD& operator/=(float b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] /= b;
		return m_data[tid];
	};

	FloatAD& operator+=(int b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] += b;
		return m_data[tid];
	};

	FloatAD& operator*=(int b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] *= b;
		return m_data[tid];
	};

	FloatAD& operator-=(int b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] -= b;
		return m_data[tid];
	};

	FloatAD& operator/=(int b) {
		int tid = mitsuba::Thread::getID();
		m_data[tid] /= b;
		return m_data[tid];
	};

	FloatAD_shared& operator+=(FloatAD_shared b) {
		for (int i = 0; i < NUM_THREADS; i++) {
			m_data[i] += b.m_data[i];
		}
		return *this;
	};
	FloatAD_shared& operator-=(FloatAD_shared b) {
		for (int i = 0; i < NUM_THREADS; i++) {
			m_data[i] -= b.m_data[i];
		}
		return *this;
	};

	FloatAD_shared& operator*=(FloatAD_shared b) {
		for (int i = 0; i < NUM_THREADS; i++) {
			m_data[i] *= b.m_data[i];
		}
		return *this;
	};	
	FloatAD_shared& operator/=(FloatAD_shared b) {
		for (int i = 0; i < NUM_THREADS; i++) {
			m_data[i] /= b.m_data[i];
		}
		return *this;
	};	



      explicit operator float() {
		  int tid = mitsuba::Thread::getID();
          return (float) m_data[tid];
      }

      explicit operator size_t() {
		int tid = mitsuba::Thread::getID();
          return (size_t) m_data[tid];
      }    



	std::vector<stan::math::var> m_data;
};


	inline FloatAD_shared operator*(FloatAD_shared x, FloatAD_shared b) {
		FloatAD_shared result;
		for (int i = 0; i < NUM_THREADS; i++) {
			result.m_data[i] = x.m_data[i] * b.m_data[i];
		}
		return result;
	};	
	inline FloatAD_shared operator+(FloatAD_shared x, FloatAD_shared b) {
		FloatAD_shared result;
		for (int i = 0; i < NUM_THREADS; i++) {
			result.m_data[i] = x.m_data[i] + b.m_data[i];
		}
		return result;
	};	

	inline FloatAD_shared operator-(FloatAD_shared x, FloatAD_shared b) {
		FloatAD_shared result;
		for (int i = 0; i < NUM_THREADS; i++) {
			result.m_data[i] = x.m_data[i] - b.m_data[i];
		}
		return result;
	};	
	inline FloatAD_shared operator/(FloatAD_shared x, FloatAD_shared b) {
		FloatAD_shared result;
		for (int i = 0; i < NUM_THREADS; i++) {
			result.m_data[i] = x.m_data[i] / b.m_data[i];
		}
		return result;
	};	









inline FloatAD_shared operator+(FloatAD_shared x) {
	for (int i = 0; i < NUM_THREADS; i++) {
		x.m_data[i] = +x.m_data[i];
	}
	return x;
}
inline FloatAD_shared operator-(FloatAD_shared x) {
	for (int i = 0; i < NUM_THREADS; i++) {
		x.m_data[i] = -x.m_data[i];
	}
	return x;
}


inline bool operator<(FloatAD_shared x, FloatAD_shared y) {
	int tid = mitsuba::Thread::getID();
	return x.m_data[tid].val() < y.m_data[tid].val();
}
inline bool operator==(FloatAD_shared x, int y) {
	int tid = mitsuba::Thread::getID();
	return x.m_data[tid].val() == y;
}


inline FloatAD operator*(FloatAD_shared x,double b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] *= b;
	return x.m_data[tid];
};	

inline FloatAD operator+(FloatAD_shared x,double b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] += b;
	return x.m_data[tid];
};	


inline FloatAD operator-(FloatAD_shared x,double b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] -= b;
	return x.m_data[tid];
};	


inline FloatAD operator/(FloatAD_shared x,double b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] /= b;
	return x.m_data[tid];
};		

inline FloatAD operator*(FloatAD_shared x,FloatAD b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] *= b; 
	return x.m_data[tid];
};	

inline FloatAD operator+(FloatAD_shared x,FloatAD b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] += b; 
	return x.m_data[tid];
};	

inline FloatAD operator-(FloatAD_shared x,FloatAD b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] -= b; 
	return x.m_data[tid];
};	
inline FloatAD operator/(FloatAD_shared x,FloatAD b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] /= b; 
	return x.m_data[tid];
};	


inline FloatAD operator*(FloatAD_shared x,float b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] *= b;
	return x.m_data[tid];
};	

inline FloatAD operator+(FloatAD_shared x,float b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] += b;
	return x.m_data[tid];
};	


inline FloatAD operator-(FloatAD_shared x,float b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] -= b;
	return x.m_data[tid];
};	


inline FloatAD operator/(float b,FloatAD_shared x) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] /= b;
	return x.m_data[tid];
};		



inline FloatAD operator*(FloatAD_shared x,int b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] *= b;
	return x.m_data[tid];
};	

inline FloatAD operator+(FloatAD_shared x,int b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] += b;
	return x.m_data[tid];
};	


inline FloatAD operator-(FloatAD_shared x,int b) {
	int tid = mitsuba::Thread::getID();
	x.m_data[tid] -= b;
	return x.m_data[tid];
};	


inline FloatAD operator/(FloatAD_shared x,int b) {

	int tid = mitsuba::Thread::getID();
	x.m_data[tid] /= b;
	return x.m_data[tid];
};		










#endif
