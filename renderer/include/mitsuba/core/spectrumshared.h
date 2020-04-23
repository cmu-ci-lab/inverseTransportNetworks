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
#if !defined(__MITSUBA_CORE_SPECTRUMSHARED_H_)
#define __MITSUBA_CORE_SPECTRUMSHARED_H_

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/FloatADshared.h>
#include <mitsuba/core/FloatAD.h>
#include <stan/math/rev/core.hpp>

#if !defined(SPECTRUM_SAMPLES)
#error The desired number of spectral samples must be \
	specified in the configuration file!
#endif

#define SPECTRUM_MIN_WAVELENGTH   360
#define SPECTRUM_MAX_WAVELENGTH   830
#define SPECTRUM_RANGE                \
	(SPECTRUM_MAX_WAVELENGTH-SPECTRUM_MIN_WAVELENGTH)

MTS_NAMESPACE_BEGIN

template <typename T, int N> struct TSpectrum_shared {
public:
	typedef T          Scalar;
	/// Number of dimensions
	const static int dim = N;

	/// Create a new spectral power distribution, but don't initialize the contents
#if !defined(MTS_DEBUG_UNINITIALIZED)
	inline TSpectrum_shared() { }
#else
	inline TSpectrum_shared() {
		for (int i=0; i<N; i++)
			s[i] = std::numeric_limits<Scalar>::quiet_NaN();
	}
#endif

	/// Create a new spectral power distribution with all samples set to the given value
	explicit inline TSpectrum_shared(Scalar v) {
		for (int i=0; i<N; i++)
			s[i] = (FloatAD_shared) v;
	}

	/// Copy a spectral power distribution
	explicit inline TSpectrum_shared(Scalar spec[N]) {
		memcpy(s, spec, sizeof(Scalar)*N);
	}

	/// Unserialize a spectral power distribution from a binary data stream
	explicit inline TSpectrum_shared(Stream *stream) {
		stream->readArray(s, N);
	}

	/// Initialize with a TSpectrum data type based on a alternate representation
	template <typename AltScalar> explicit TSpectrum_shared(const TSpectrum_shared<AltScalar, N> &v) {
		for (int i=0; i<N; ++i) {
			Scalar temp = (Scalar) v[i].val();
			s[i] = (FloatAD_shared) temp;
		}
			
	}

	/// Add two spectral power distributions
	inline TSpectrum_shared operator+(const TSpectrum_shared &spec) const {
		TSpectrum_shared value = *this;
		for (int i=0; i<N; i++)
			value.s[i] += spec.s[i];
		return value;
	}

	/// Add a spectral power distribution to this instance
	inline TSpectrum_shared& operator+=(const TSpectrum_shared &spec) {
		for (int i=0; i<N; i++)
			s[i] += spec.s[i];
		return *this;
	}

	/// Subtract a spectral power distribution
	inline TSpectrum_shared operator-(const TSpectrum_shared &spec) const {
		TSpectrum_shared value = *this;
		for (int i=0; i<N; i++)
			value.s[i] -= spec.s[i];
		return value;
	}

	/// Subtract a spectral power distribution from this instance
	inline TSpectrum_shared& operator-=(const TSpectrum_shared &spec) {
		for (int i=0; i<N; i++)
			s[i] -= spec.s[i];
		return *this;
	}

	/// Multiply by a scalar
	inline TSpectrum_shared operator*(Scalar f) const {
		TSpectrum_shared value = *this;
		for (int i=0; i<N; i++)
			value.s[i] *= f;
		return value;
	}

	/// Multiply by a scalar
	inline friend TSpectrum_shared operator*(Scalar f, const TSpectrum_shared &spec) {
		return spec * f;
	}

	/// Multiply by a scalar
	inline TSpectrum_shared& operator*=(Scalar f) {
		for (int i=0; i<N; i++)
			s[i] *= f;
		return *this;
	}



	/// Multiply by a scalar
	inline TSpectrum_shared operator*(Float f) const {
		TSpectrum_shared value = *this;
		for (int i=0; i<N; i++)
			value.s[i] *= f;
		return value;
	}

	/// Multiply by a scalar
	inline friend TSpectrum_shared operator*(Float f, const TSpectrum_shared &spec) {
		return spec * f;
	}

	/// Multiply by a scalar
	inline TSpectrum_shared& operator*=(Float f) {
		for (int i=0; i<N; i++)
			s[i] *= f;
		return *this;
	}






	/// Perform a component-wise multiplication by another spectrum_AD
	inline TSpectrum_shared operator*(const TSpectrum_shared &spec) const {
		TSpectrum_shared value = *this;
		for (int i=0; i<N; i++)
			value.s[i] *= spec.s[i];
		return value;
	}

	/// Perform a component-wise multiplication by another spectrum_AD
	inline TSpectrum_shared& operator*=(const TSpectrum_shared &spec) {
		for (int i=0; i<N; i++)
			s[i] *= spec.s[i];
		return *this;
	}



	/// Perform a component-wise division by another spectrum
	inline TSpectrum_shared& operator/=(const TSpectrum_shared &spec) {
		for (int i=0; i<N; i++)
			s[i] /= spec.s[i];
		return *this;
	}

	/// Perform a component-wise division by another spectrum
	inline TSpectrum_shared operator/(const TSpectrum_shared &spec) const {
		TSpectrum_shared value = *this;
		for (int i=0; i<N; i++)
			value.s[i] /= spec.s[i];
		return value;
	}

	/// Divide by a scalar
	inline TSpectrum_shared operator/(Scalar f) const {
		TSpectrum_shared value = *this;
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "TSpectrum_shared: Division by zero!");
#endif
		for (int i=0; i<N; i++)
			value.s[i] *= FloatAD_shared((1.0f / f).val());
		return value;
	}

	/// Equality test
	inline bool operator==(const TSpectrum_shared &spec) const {
		for (int i=0; i<N; i++) {
			if (s[i] != spec.s[i])
				return false;
		}
		return true;
	}

	/// Inequality test
	inline bool operator!=(const TSpectrum_shared &spec) const {
		return !operator==(spec);
	}

	/// Divide by a scalar
	inline friend TSpectrum_shared operator/(Scalar f, TSpectrum_shared &spec) {
		return TSpectrum_shared(f) / spec;
	}

	/// Divide by a scalar
	inline TSpectrum_shared& operator/=(Scalar f) {
#ifdef MTS_DEBUG
		if (f == 0)
			SLog(EWarn, "TTSpectrum_AD: Division by zero!");
#endif
		FloatAD_shared recip = 1.0f / f;
		for (int i=0; i<N; i++)
			s[i] *= recip;
		return *this;
	}

	/// Check for NaNs
	inline bool isNaN() const {
		for (int i=0; i<N; i++)
			if (std::isnan(s[i]))
				return true;
		return false;
	}

	/// Returns whether the spectrum only contains valid (non-NaN, nonnegative) samples
	inline bool isValid() const {
		for (int i=0; i<N; i++)
			if (!std::isfinite(s[i]) || s[i] < 0.0f)
				return false;
		return true;
	}

	/// Multiply-accumulate operation, adds \a weight * \a spec
	inline void addWeighted(Scalar weight, const TSpectrum_shared &spec) {
		for (int i=0; i<N; i++)
			s[i] += weight * spec.s[i];
	}

	/// Return the average over all wavelengths
	inline FloatAD_shared average() const {
		FloatAD_shared result = 0.0f;
		for (int i=0; i<N; i++)
			result += s[i];
		return result * (1.0f / N);
	}

	/// Component-wise absolute value
	inline TSpectrum_shared abs() const {
		TSpectrum_shared value;
		for (int i=0; i<N; i++)
			value.s[i] = std::abs(s[i]);
		return value;
	}

	/// Component-wise square root
	inline TSpectrum_shared sqrt() const {
		TSpectrum_shared value;
		for (int i=0; i<N; i++)
			value.s[i] = std::sqrt(s[i]);
		return value;
	}

	/// Component-wise square root
	inline TSpectrum_shared safe_sqrt() const {
		TSpectrum_shared value;
		for (int i=0; i<N; i++)
			value.s[i] = math::safe_sqrt(s[i]);
		return value;
	}

	/// Component-wise logarithm
	inline TSpectrum_shared log() const {
		TSpectrum_shared value;
		for (int i=0; i<N; i++)
			value.s[i] = math::fastlog(s[i]);
		return value;
	}

	/// Component-wise exponentation
	inline TSpectrum_shared exp() const {
		TSpectrum_shared value;
		for (int i=0; i<N; i++)
			value.s[i] = math::fastexp(s[i]);
		return value;
	}

	/// Component-wise power
	inline TSpectrum_shared pow(Scalar f) const {
		TSpectrum_shared value;
		for (int i=0; i<N; i++)
			value.s[i] = std::pow(s[i], f);
		return value;
	}

	/// Clamp negative values
	inline void clampNegative() {
		for (int i=0; i<N; i++)
			s[i] = std::max((Scalar) 0.0f, s[i]);
	}

	/// Return the highest-valued spectral sample
	inline FloatAD_shared max() const {
		FloatAD_shared result = s[0];
		for (int i=1; i<N; i++)
			result = (FloatAD_shared)std::max(result.val(), s[i].val());
		return result;
	}

	/// Return the lowest-valued spectral sample
	inline FloatAD_shared min() const {
		FloatAD_shared result = s[0];
		for (int i=1; i<N; i++)
			result = (FloatAD_shared)std::min(result.val(), s[i].val());
		return result;
	}

	/// Negate
	inline TSpectrum_shared operator-() const {
		TSpectrum_shared value;
		for (int i=0; i<N; i++)
			value.s[i] = -s[i];
		return value;
	}

	/// Indexing operator
	inline FloatAD_shared &operator[](int entry) {
		return s[entry];
	}

	/// Indexing operator
	inline FloatAD_shared operator[](int entry) const {
		return s[entry];
	}

	/// Check if this spectrum is zero at all wavelengths
	inline bool isZero() const {
		for (int i=0; i<N; i++) {
			if (s[i] != 0.0f)
				return false;
		}
		return true;
	}

	/// Serialize this spectrum to a stream
	inline void serialize(Stream *stream) const {
		stream->writeArray(s, N);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "]";
		return oss.str();
	}

protected:
	FloatAD_shared s[N];
};


struct MTS_EXPORT_CORE Spectrum_shared : public TSpectrum_shared<FloatAD_shared, SPECTRUM_SAMPLES> {
public:
	typedef TSpectrum_shared<FloatAD_shared, SPECTRUM_SAMPLES> Parent;

	/**
	 * \brief When converting from RGB reflectance values to
	 * discretized color spectra, the following `intent' flag
	 * can be provided to improve the results of this highly
	 * under-constrained problem.
	 */
	enum EConversionIntent {
		/// Unitless reflectance data is converted
		EReflectance,

		/// Radiance-valued illumination data is converted
		EIlluminant
	};

	/// Create a new spectral power distribution, but don't initialize the contents

	/// Construct from a TSpectrum instance
	inline Spectrum_shared(const Parent &s) : Parent(s) { }
	/// Initialize with a TSpectrum data type based on a alternate representation
	template <typename AltScalar> explicit Spectrum_shared(const TSpectrum_shared<AltScalar, SPECTRUM_SAMPLES> &v) {
		for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
			Scalar temp = (Scalar) v[i];
			s[i] = (FloatAD_shared) temp;
		}
	}


#if SPECTRUM_SAMPLES == 3
	inline void fromLinearRGB(FloatAD_shared r, FloatAD_shared g, FloatAD_shared b,
		EConversionIntent intent = EReflectance /* unused */) {
	/* Nothing to do -- the renderer is in RGB mode */
	s[0] = r; s[1] = g; s[2] = b;
	}
#else
	
	void fromLinearRGB(FloatAD_shared r, FloatAD_shared g, FloatAD_shared b,
			EConversionIntent intent = EReflectance);
#endif


	double val() const {
		int tid = mitsuba::Thread::getID();			
		return s[0].m_data[tid].val();
	}
	double adj() const {
		int tid = mitsuba::Thread::getID();			
		return s[0].m_data[tid].adj();
	}	



	Spectrum_shared() {
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp();
			s[i] = temp;
		}
	}

	Spectrum_shared(float x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}

	Spectrum_shared(double x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}

	Spectrum_shared(long double x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}

	Spectrum_shared(bool x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}


	Spectrum_shared(char x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}
	Spectrum_shared(short x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}

	Spectrum_shared(int x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}

	Spectrum_shared(long x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}

	Spectrum_shared(unsigned char x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}	

	Spectrum_shared(unsigned short x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}	

	Spectrum_shared(unsigned int x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}		
	Spectrum_shared(unsigned long x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}

	Spectrum_shared(FloatAD_shared x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			
			s[i] = x;
		}
	}

#ifdef _WIN64

	Spectrum_shared(size_t x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}
	Spectrum_shared(ptrdiff_t x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}
	
#endif

#ifdef BOOST_MATH_USE_FLOAT128

	Spectrum_shared(__float128 x){
		for (int i = 0; i < SPECTRUM_SAMPLES; i++) {
			FloatAD_shared temp(x);
			s[i] = temp;
		}
	}

#endif	

	/// Copy a spectral power distribution
	explicit inline Spectrum_shared(FloatAD_shared value[SPECTRUM_SAMPLES]) {
		memcpy(s, value, sizeof(FloatAD_shared)*SPECTRUM_SAMPLES);
	}

	/// Unserialize a spectral power distribution from a binary data stream
	explicit inline Spectrum_shared(Stream *stream) : Parent(stream) { }

};


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, float y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] *= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, float y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] += y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, float y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] -= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, float y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] /= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}





inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, double y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] *= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, double y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] += y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, double y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] -= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, double y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] /= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, int y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] *= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, int y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] += y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, int y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] -= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, int y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] /= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, FloatAD y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] *= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, FloatAD y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] += y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, FloatAD y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] -= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, FloatAD y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] /= y;
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, FloatAD_shared y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] *= y.m_data[tid];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, FloatAD_shared y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] += y.m_data[tid];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, FloatAD_shared y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] -= y.m_data[tid];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, FloatAD_shared y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] /= y.m_data[tid];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] *= y[0];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] += y[0];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] -= y[0];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, mitsuba::Spectrum y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] /= y[0];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}



inline mitsuba::Spectrum_AD operator*(Spectrum_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] *= y[0];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}

inline mitsuba::Spectrum_AD operator+(Spectrum_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] += y[0];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}
inline mitsuba::Spectrum_AD operator-(Spectrum_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] -= y[0];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
}


inline mitsuba::Spectrum_AD operator/(Spectrum_shared x, mitsuba::Spectrum_AD y)  {
	int tid = mitsuba::Thread::getID(); 
	x[0].m_data[tid] /= y[0];
	return mitsuba::Spectrum_AD(x[0].m_data[tid]);
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



inline Spectrum_shared operator+(Spectrum_shared x, Spectrum_shared y)  {

	Spectrum_shared results;
	for (int i = 0; i < NUM_THREADS; i++) {
		x[0].m_data[i] += y[0].m_data[i];
	}
	
	return x;
}
inline Spectrum_shared operator*(Spectrum_shared x, Spectrum_shared y)  {
	
	Spectrum_shared results;
	for (int i = 0; i < NUM_THREADS; i++) {
		x[0].m_data[i] *= y[0].m_data[i];
	}
	
	return x;
}
inline Spectrum_shared operator-(Spectrum_shared x, Spectrum_shared y)  {
	
	Spectrum_shared results;
	for (int i = 0; i < NUM_THREADS; i++) {
		x[0].m_data[i] -= y[0].m_data[i];
	}
	
	return x;
}
inline Spectrum_shared operator/(Spectrum_shared x, Spectrum_shared y)  {
	
	Spectrum_shared results;
	for (int i = 0; i < NUM_THREADS; i++) {
		x[0].m_data[i] /= y[0].m_data[i];
	}
	
	return x;
}

MTS_NAMESPACE_END

#endif 


