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
#if !defined(__MITSUBA_RENDER_TEXTUREAD_H_)
#define __MITSUBA_RENDER_TEXTUREAD_H_

#include <mitsuba/core/cobject.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shader.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Base class of all textures. Computes values for an arbitrary surface
 * point. \ref Texture2D is a specialization to UV-based textures.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER TextureAD : public ConfigurableObject, public HWResource {
public:
	/**
	 * \brief Return the texture value at \c its
	 * \param filter
	 *    Specifies whether a filtered texture lookup is desired. Note
	 *    that this does not mean that filtering will actually be used.
	 */ 
	virtual Spectrum_AD eval(const Intersection &its, bool filter = true) const;
	virtual Spectrum_AD eval_AD() const;
	virtual void eval3(const Intersection &its, FloatAD &val1, FloatAD &val2, FloatAD &val3, bool filter = true) const;
	virtual void getPixel(const Point2 &uv, FloatAD &val1, FloatAD &val2,bool filter = true) const;	
	/**
	 * \brief Return the texture gradient at \c its
	 *
	 * The parameter \c gradient should point to an array with space for
	 * two \ref Spectrum data structures corresponding to the U and V derivative.
	 *
	 * \remark This function is usually implemented pointwise without any kind of filtering.
	 * The Python signature is <tt>dx, dy = tex.evalGradient(its)</tt>.
	 */
	virtual void evalGradient(const Intersection &its, Spectrum_AD *gradient) const;
	virtual void evalGradient3(const Intersection &its, FloatAD &val00,  FloatAD &val01, FloatAD &val02, FloatAD &val10, FloatAD &val11, FloatAD &val12) const;
	/// Return the component-wise average value of the texture over its domain
	virtual Spectrum_AD getAverage() const;

	/// Return the component-wise minimum of the texture over its domain
	virtual Spectrum_AD getMinimum() const;

	/// Return the component-wise maximum of the texture over its domain
	virtual Spectrum_AD getMaximum() const;

	/// Return the resolution in pixels, if applicable
	virtual Vector3i getResolution() const;

	/// Return whether the texture takes on a constant value everywhere
	virtual bool isConstant() const;

	/**
	 * \brief Return whether the texture is monochromatic / spectrally uniform
	 *
	 * The implementation may conservatively return \c false if it is not sure.
	 */
	virtual bool isMonochromatic() const;

	/**
	 * \brief Does this texture perform any pre-filtering when
	 * ray differentials are available?
	 */
	virtual bool usesRayDifferentials() const;

	/**
	 * \brief Some textures are only proxies for an actual
	 * implementation. This function returns the actual
	 * texture implementation to be used.
	 *
	 * The default implementation returns <tt>this</tt>.
	 */
	virtual ref<TextureAD> expand();

	/// Serialize to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/**
	 * \brief Return a bitmap representation of the texture
	 *
	 * When the class implementing this interface is a bitmap-backed texture,
	 * this function directly returns the underlying bitmap. When it is procedural,
	 * a bitmap version must first be generated. In this case, the parameter
	 * \ref sizeHint is used to control the target size. The default
	 * value <tt>-1, -1</tt> allows the implementation to choose a suitable
	 * size by itself.
	 */
	virtual ref<BitmapAD> getBitmap(const Vector2i &sizeHint = Vector2i(-1, -1)) const;

	MTS_DECLARE_CLASS()
protected:
	TextureAD(const Properties &props);
	TextureAD(Stream *stream, InstanceManager *manager);

	virtual ~TextureAD();
};

/**
 * \brief Base class of all 2D textures
 * \ingroup librender
 */
class MTS_EXPORT_RENDER Texture2DAD : public TextureAD {
public:
	/**
	 * \brief Return the texture value at \c its
	 * \param filter
	 *    Specifies whether a filtered texture lookup is desired. Note
	 *    that this does not mean that filtering will actually be used.
	 */
	Spectrum_AD eval(const Intersection &its, bool filter = true) const;
	void eval3(const Intersection &its, FloatAD &val1, FloatAD &val2, FloatAD &val3, bool filter = true) const;
	void getPixel(const Point2 &uv, FloatAD &val1, FloatAD &val2,bool filter = true) const;	
	/**
	 * \brief Return the texture gradient at \c its
	 *
	 * The parameter \c gradient should point to an array with space for
	 * two \ref Spectrum data structures corresponding to the U and V derivative.
	 *
	 * \remark This function is usually implemented pointwise without any kind of filtering.
	 */
	void evalGradient(const Intersection &its, Spectrum_AD *gradient) const;
	void evalGradient3(const Intersection &its, FloatAD &val00,  FloatAD &val01, FloatAD &val02, FloatAD &val10, FloatAD &val11, FloatAD &val12) const;
	/// Serialize to a binary data stream
	virtual void serialize(Stream *stream, InstanceManager *manager) const;

	/// Unfiltered texture lookup -- Texture2D subclasses must provide this function
	virtual Spectrum_AD eval(const Point2 &uv) const = 0;
	/// Filtered texture lookup -- Texture2D subclasses must provide this function
	virtual Spectrum_AD eval(const Point2 &uv, const Vector2 &d0,
			const Vector2 &d1) const = 0;
	virtual void eval3(const Point2 &uv, const Vector2 &d0,
			const Vector2 &d1, FloatAD &val1, FloatAD &val2, FloatAD &val3) const = 0;

	virtual void getPixel3(const Point2 &uv, FloatAD &val1, FloatAD &val2) const = 0;	
	virtual void eval3(const Point2 &uv, FloatAD &val1, FloatAD &val2, FloatAD &val3) const = 0;


	/// Unfiltered radient lookup lookup -- Texture2D subclasses can optionally provide this function
	virtual void evalGradient(const Point2 &uv, Spectrum_AD *gradient) const;
	virtual void evalGradient3(const Point2 &uv, FloatAD &val00,  FloatAD &val01, FloatAD &val02, FloatAD &val10, FloatAD &val11, FloatAD &val12) const;
	/**
	 * \brief Return a bitmap representation of the texture
	 *
	 * When the class implementing this interface is a bitmap-backed texture,
	 * this function directly returns the underlying bitmap. When it is procedural,
	 * a bitmap version must first be generated. In this case, the parameter
	 * \ref sizeHint is used to control the target size. The default
	 * value <tt>-1, -1</tt> allows the implementation to choose a suitable
	 * size by itself.
	 */
	virtual ref<BitmapAD> getBitmap(const Vector2i &sizeHint = Vector2i(-1, -1)) const;

	MTS_DECLARE_CLASS()
protected:
	Texture2DAD(const Properties &props);
	Texture2DAD(Stream *stream, InstanceManager *manager);

	virtual ~Texture2DAD();
protected:
	Point2 m_uvOffset;
	Vector2 m_uvScale;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_TEXTUREAD_H_ */
