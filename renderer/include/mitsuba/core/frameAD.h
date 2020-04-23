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
#if !defined(__MITSUBA_CORE_FRAMEAD_H_)
#define __MITSUBA_CORE_FRAMEAD_H_

#include <mitsuba/mitsuba.h>
#include <stan/math/fwd/scal/fun/sqrt.hpp>
MTS_NAMESPACE_BEGIN

/**
 * \brief Stores a three-dimensional orthonormal coordinate frame
 *
 * This class is mostly used to quickly convert between different
 * cartesian coordinate systems and to efficiently compute certain
 * quantities (e.g. \ref cosTheta(), \ref tanTheta, ..).
 *
 * \ingroup libcore
 * \ingroup libpython
 */
struct FrameAD {
	VectorAD s, t;


	/// Default constructor -- performs no initialization!
	inline FrameAD() { }

	/// Given a normal and tangent vectors, construct a new coordinate frame


	/// Construct a frame from the given orthonormal vectors
	inline FrameAD(const VectorAD &x, const VectorAD &y, const VectorAD &z)
	 : s(x), t(y), n(z) {
	}

	/// Construct a new coordinate frame from a single vector
	inline FrameAD(const VectorAD &n) : n(n) {
		// coordinateSystem(n, s, t);
	}

	operator Frame() {
		Frame result;
		result.s.x = s.x.val();
		result.s.y = s.y.val();
		result.s.z = s.z.val();

		result.t.x = t.x.val();
		result.t.y = t.y.val();
		result.t.z = t.z.val();

		result.n.x = n.x.val();
		result.n.y = n.y.val();
		result.n.z = n.z.val();
		return result;
	}
	/// Unserialize from a binary data stream
	inline FrameAD(Stream *stream) {
		s = VectorAD(stream);
		t = VectorAD(stream);
		n = NormalAD(stream);
	}

	/// Serialize to a binary data stream
	inline void serialize(Stream *stream) const {
		s.serialize(stream);
		t.serialize(stream);
		n.serialize(stream);
	}

	/// Convert from world coordinates to local coordinates
	inline VectorAD toLocal(const VectorAD &v) const {
		return VectorAD(
			dot(v, s),
			dot(v, t),
			dot(v, n)
		);
	}

	inline VectorAD toLocal(const Vector &v) const {
		return VectorAD(
			v.x * s.x + v.y * s.y + v.z * s.z,
			v.x * t.x + v.y * t.y + v.z * t.z,
			v.x * n.x + v.y * n.y + v.z * n.z			
		);
		
	}

	/// Convert from local coordinates to world coordinates
	inline VectorAD toWorld(const VectorAD &v) const {
		return s * v.x + t * v.y + n * v.z;
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared cosine of the angle between the normal and v */
	inline static FloatAD cosTheta2(const VectorAD &v) {
		return v.z * v.z;
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the cosine of the angle between the normal and v */
	inline static FloatAD cosTheta(const VectorAD &v) {
		return v.z;
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the u and v coordinates of the vector 'v' */
	inline static Vector2AD uv(const VectorAD &v) {
		return Vector2AD(v.x, v.y);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared sine of the angle between the normal and v */
	inline static FloatAD sinTheta2(const VectorAD &v) {
		return 1.0f - v.z * v.z;
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the sine of the angle between the normal and v */
	inline static FloatAD sinTheta(const VectorAD &v) {
		FloatAD temp = sinTheta2(v);
		if (temp <= 0.0f)
			return 0.0f;
		return sqrt(temp);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the tangent of the angle between the normal and v */
	inline static FloatAD tanTheta(const VectorAD &v) {
		FloatAD temp = 1 - v.z*v.z;
		if (temp <= 0.0f)
			return 0.0f;
		return sqrt(temp) / v.z;
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared tangent of the angle between the normal and v */
	inline static FloatAD tanTheta2(const VectorAD &v) {
		FloatAD temp = 1 - v.z*v.z;
		if (temp <= 0.0f)
			return 0.0f;
		return temp / (v.z * v.z);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the sine of the phi parameter in spherical coordinates */
	inline static FloatAD sinPhi(const VectorAD &v) {
		FloatAD sinTheta = FrameAD::sinTheta(v);
		if (sinTheta == 0.0f)
			return 1.0f;
		return math::clamp(v.y / sinTheta, (FloatAD) -1.0f, (FloatAD) 1.0f);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the cosine of the phi parameter in spherical coordinates */
	inline static FloatAD cosPhi(const VectorAD &v) {
		FloatAD sinTheta = FrameAD::sinTheta(v);
		if (sinTheta == 0.0f)
			return 1.0f;
		return math::clamp(v.x / sinTheta, (FloatAD) -1.0f, (FloatAD) 1.0f);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared sine of the phi parameter in  spherical
	 * coordinates */
	inline static FloatAD sinPhi2(const VectorAD &v) {
		return math::clamp(v.y * v.y / sinTheta2(v), (FloatAD) 0.0f, (FloatAD) 1.0f);
	}

	/** \brief Assuming that the given direction is in the local coordinate
	 * system, return the squared cosine of the phi parameter in  spherical
	 * coordinates */
	inline static FloatAD cosPhi2(const VectorAD &v) {
		return math::clamp(v.x * v.x / sinTheta2(v), (FloatAD) 0.0f, (FloatAD) 1.0f);
	}

	/// Equality test
	inline bool operator==(const FrameAD &frame) const {
		return frame.s == s && frame.t == t && frame.n == n;
	}

	/// Inequality test
	inline bool operator!=(const FrameAD &frame) const {
		return !operator==(frame);
	}

	/// Return a string representation of this frame
	inline std::string toString() const {
		std::ostringstream oss;
		oss << "Frame[" << std::endl
			<< "]";
		return oss.str();
	}
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_CORE_FRAMEAD_H_ */
