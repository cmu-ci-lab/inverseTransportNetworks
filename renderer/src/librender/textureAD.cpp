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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/mipmapAD.h>

MTS_NAMESPACE_BEGIN

namespace stats {
	StatsCounter mipStorageAD("Texture system", "Cumulative MIP map memory allocations", EByteCount);
	StatsCounter clampedAnisotropyAD("Texture system", "Lookups with clamped anisotropy", EPercentage);
	StatsCounter avgEWASamplesAD("Texture system", "Average EWA samples / lookup", EAverage);
	StatsCounter filteredLookupsAD("Texture system", "Filtered texture lookups", EPercentage);

}

TextureAD::TextureAD(const Properties &props)
 : ConfigurableObject(props) {
}

TextureAD::TextureAD(Stream *stream, InstanceManager *manager)
 : ConfigurableObject(stream, manager) {
}

Vector3i TextureAD::getResolution() const {
	return Vector3i(0);
}
void TextureAD::eval3(const Intersection &its, FloatAD &val1, FloatAD&val2, FloatAD&val3, bool filter) const {NotImplementedError("eval3");}
void TextureAD::getPixel(const Point2 &uv, FloatAD &val1, FloatAD &val2, bool filter) const{NotImplementedError("getpixel");}
Spectrum_AD TextureAD::eval(const Intersection &its, bool filter) const { NotImplementedError("eval"); }
Spectrum_AD TextureAD::eval_AD() const { NotImplementedError("evalAD"); }
Spectrum_AD TextureAD::getAverage() const { NotImplementedError("getAverage"); }
Spectrum_AD TextureAD::getMinimum() const { NotImplementedError("getMinimum"); }
Spectrum_AD TextureAD::getMaximum() const { NotImplementedError("getMaximum"); }
bool TextureAD::isConstant() const { NotImplementedError("isConstant"); }
bool TextureAD::isMonochromatic() const { NotImplementedError("isMonochromatic"); }
bool TextureAD::usesRayDifferentials() const { NotImplementedError("usesRayDifferentials"); }
ref<BitmapAD> TextureAD::getBitmap(const Vector2i &) const { NotImplementedError("getBitmapAD"); }

ref<TextureAD> TextureAD::expand() {
	return this;
}

void TextureAD::evalGradient(const Intersection &_its, Spectrum_AD *gradient) const {
	const Float eps = Epsilon;
	Intersection its(_its);

	Spectrum_AD value = eval(its, false);

	its.p = _its.p + its.dpdu * eps;
	its.uv = _its.uv + Point2(eps, 0);
	Spectrum_AD valueU = eval(its, false);

	its.p = _its.p + its.dpdv * eps;
	its.uv = _its.uv + Point2(0, eps);
	Spectrum_AD valueV = eval(its, false);

	gradient[0] = (valueU - value)*(1/eps);
	gradient[1] = (valueV - value)*(1/eps);
}

void TextureAD::evalGradient3(const Intersection &_its, FloatAD &val00,  FloatAD &val01, FloatAD &val02, FloatAD &val10, FloatAD &val11, FloatAD &val12) const {
	const Float eps = Epsilon;
	Intersection its(_its);
	FloatAD temp1, temp2, temp3, tempU1, tempU2, tempU3, tempV1, tempV2, tempV3;
	eval3(its,temp1, temp2, temp3, false);

	its.p = _its.p + its.dpdu * eps;
	its.uv = _its.uv + Point2(eps, 0);

	eval3(its,tempU1, tempU2, tempU3, false);

	its.p = _its.p + its.dpdv * eps;
	its.uv = _its.uv + Point2(0, eps);

	eval3(its, tempV1, tempV2, tempV3,false);

	val00 = (tempU1 - temp1) * (1/eps);
	val01 = (tempU2 - temp2) * (1/eps);
	val02 = (tempU3 - temp3) * (1/eps);		

	val10 = (tempV1 - temp1) * (1/eps);
	val11 = (tempV2 - temp2) * (1/eps);
	val12 = (tempV3 - temp3) * (1/eps);		
}

TextureAD::~TextureAD() { }

void TextureAD::serialize(Stream *stream, InstanceManager *manager) const {
	ConfigurableObject::serialize(stream, manager);
}

Texture2DAD::Texture2DAD(const Properties &props) : TextureAD(props) {
	if (props.getString("coordinates", "uv") == "uv") {
		m_uvOffset = Point2(
			props.getFloat("uoffset", 0.0f),
			props.getFloat("voffset", 0.0f)
		);
		Float uvscale = props.getFloat("uvscale", 1.0f);
		m_uvScale = Vector2(
			props.getFloat("uscale", uvscale),
			props.getFloat("vscale", uvscale)
		);
	} else {
		Log(EError, "Only UV coordinates are supported at the moment!");
	}
}
void Texture2DAD::eval3(const Intersection &its, FloatAD &val1, FloatAD &val2, FloatAD &val3, bool filter) const {
	Point2 uv = Point2(its.uv.x * m_uvScale.x, its.uv.y * m_uvScale.y) + m_uvOffset;
	if (its.hasUVPartials && filter) {
		eval3(uv,
			Vector2(its.dudx * m_uvScale.x, its.dvdx * m_uvScale.y),
			Vector2(its.dudy * m_uvScale.x, its.dvdy * m_uvScale.y),
			val1, val2, val3);
	} else {		
		eval3(uv, val1, val2, val3);
	}
}





Texture2DAD::Texture2DAD(Stream *stream, InstanceManager *manager)
 : TextureAD(stream, manager) {
	m_uvOffset = Point2(stream);
	m_uvScale = Vector2(stream);
}

Texture2DAD::~Texture2DAD() {
}

void Texture2DAD::serialize(Stream *stream, InstanceManager *manager) const {
	TextureAD::serialize(stream, manager);
	m_uvOffset.serialize(stream);
	m_uvScale.serialize(stream);
}

Spectrum_AD Texture2DAD::eval(const Intersection &its, bool filter) const {
	Point2 uv = Point2(its.uv.x * m_uvScale.x, its.uv.y * m_uvScale.y) + m_uvOffset;
	if (its.hasUVPartials && filter) {
		return eval(uv,
			Vector2(its.dudx * m_uvScale.x, its.dvdx * m_uvScale.y),
			Vector2(its.dudy * m_uvScale.x, its.dvdy * m_uvScale.y));
	} else {
		return eval(uv);
	}
}
void Texture2DAD::getPixel(const Point2 &uv, FloatAD &val1, FloatAD &val2, bool filter) const {
	getPixel3(uv, val1, val2);
}
void Texture2DAD::evalGradient(const Intersection &its, Spectrum_AD *gradient) const {
	Point2 uv = Point2(its.uv.x * m_uvScale.x, its.uv.y * m_uvScale.y) + m_uvOffset;

	evalGradient(uv, gradient);

	gradient[0] *= m_uvScale.x;
	gradient[1] *= m_uvScale.y;
}

void Texture2DAD::evalGradient3(const Intersection &its, FloatAD &val00,  FloatAD &val01, FloatAD &val02, FloatAD &val10, FloatAD &val11, FloatAD &val12) const {
	Point2 uv = Point2(its.uv.x * m_uvScale.x, its.uv.y * m_uvScale.y) + m_uvOffset;

	evalGradient3(uv,val00, val01, val02,val10, val11, val12);

	val00 *= m_uvScale.x;
	val01 *= m_uvScale.x;
	val02 *= m_uvScale.x;

	val10 *= m_uvScale.y;
	val11 *= m_uvScale.y;
	val12 *= m_uvScale.y;
}

void Texture2DAD::evalGradient(const Point2 &uv, Spectrum_AD *gradient) const {
	const Float eps = Epsilon;

	Spectrum_AD value = eval(uv);
	Spectrum_AD valueU = eval(uv + Vector2(eps, 0));
	Spectrum_AD valueV = eval(uv + Vector2(0, eps));

	gradient[0] = (valueU - value)*(1/eps);
	gradient[1] = (valueV - value)*(1/eps);
}
void Texture2DAD::evalGradient3(const Point2 &uv, FloatAD &val00,  FloatAD &val01, FloatAD &val02, FloatAD &val10, FloatAD &val11, FloatAD &val12) const {
	const Float eps = Epsilon;
	FloatAD temp1, temp2, temp3, tempU1, tempU2, tempU3, tempV1, tempV2, tempV3;

	eval3(uv, temp1, temp2, temp3);
	eval3(uv + Vector2(eps, 0),tempU1, tempU2, tempU3);

	eval3(uv + Vector2(0, eps), tempV1, tempV2, tempV3);

	val00 = (tempU1 - temp1) * (1/eps);
	val01 = (tempU2 - temp2) * (1/eps);
	val02 = (tempU3 - temp3) * (1/eps);		

	val10 = (tempV1 - temp1) * (1/eps);
	val11 = (tempV2 - temp2) * (1/eps);
	val12 = (tempV3 - temp3) * (1/eps);		

}
ref<BitmapAD> Texture2DAD::getBitmap(const Vector2i &sizeHint) const {
	Vector2i res(sizeHint);
	if (res.x <= 0 || res.y <= 0)
		res = Vector2i(32);

	Float invX = 1.0f / res.x, invY = 1.0f / res.y;

	ref<BitmapAD> bitmapAD = new BitmapAD(BitmapAD::ESpectrum, BitmapAD::EFloat, res);
	Spectrum_AD *target = (Spectrum_AD *) bitmapAD->getFloatData();
	for (int y=0; y<res.y; ++y)
		for (int x=0; x<res.x; ++x)
			*target++ = eval(Point2((x + 0.5f) * invX, (y + 0.5f) * invY));
	return bitmapAD;
}

MTS_IMPLEMENT_CLASS(TextureAD, true, ConfigurableObject)
MTS_IMPLEMENT_CLASS(Texture2DAD, true, TextureAD)
MTS_NAMESPACE_END
