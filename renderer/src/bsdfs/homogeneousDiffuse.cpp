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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include <stan/math.hpp>


MTS_NAMESPACE_BEGIN

class HomogeneousSmoothDiffuse : public BSDF {
public:
	HomogeneousSmoothDiffuse(const Properties &props)
		: BSDF(props) {
		int idx = 0;
		for (int i = 0; i < props.id_AD_Spectrum_AD.size(); i++) {

			if (props.id_AD_Spectrum_AD[i] == "albedo") {
				idx = i;
				break;
			}
		}
		
		m_reflectance = props.m_theta[idx];
		

	}

	HomogeneousSmoothDiffuse(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_reflectance = this->m_reflectance;

		configure();
	}

	void configure() {
		/* Verify the input parameter and fix them if necessary */
		//m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);

		m_components.clear();
		m_components.push_back(EDiffuseReflection | EFrontSide);
		BSDF::configure();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		Spectrum ret;
		
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) ret[i] = m_reflectance.val();
		return ret;
	}

	/* Just return Spectrum_AD * cos*/
	Spectrum_AD eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		

		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum_AD(0.0f);
		return m_reflectance*(INV_PI * Frame::cosTheta(bRec.wo));
	}


	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;
		
		return warp::squareToCosineHemispherePdf(bRec.wo);
	}


	FloatAD pdfAD(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return FloatAD(0.0f);
		
		return FloatAD(warp::squareToCosineHemispherePdf(bRec.wo));
	}


	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		bRec.wo = warp::squareToCosineHemisphere(sample);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		Spectrum ret;
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) ret[i] = m_reflectance.val();
		return ret;
	}
	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum_AD(0.0f);

		bRec.wo = warp::squareToCosineHemisphere(sample);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;

		Spectrum_AD ret;
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) ret[i] = m_reflectance[i].m_data[Thread::getID()];
		return ret;
	}
	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		bRec.wo = warp::squareToCosineHemisphere(sample);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		pdf = warp::squareToCosineHemispherePdf(bRec.wo);
		
		Spectrum ret;
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) ret[i] = m_reflectance.val();
		return ret;
	}

	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum_AD(0.0f);

		bRec.wo = warp::squareToCosineHemisphere(sample);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EDiffuseReflection;
		pdf = warp::squareToCosineHemispherePdf(bRec.wo);

		Spectrum_AD ret;
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) ret[i] = m_reflectance[i].m_data[Thread::getID()];
		return ret;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
				&& (name == "reflectance" || name == "diffuseReflectance")) {
			m_reflectance = this->m_reflectance;
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);
		cout << "homogenousDiffuse here" << endl;
		
		//manager->serialize(stream, *m_reflectance.val());
	}

	FloatAD getRoughness(const Intersection &its, int component) const {
		return FloatAD(std::numeric_limits<Float>::infinity());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HomogeneousSmoothDiffuse[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  reflectance = [";
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
			if ( i ) oss << ", ";
			oss << m_reflectance.val();
		}
		oss << "]" << endl
		    << "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum_shared m_reflectance;
	Spectrum_shared m_reflectance_new;
};

MTS_IMPLEMENT_CLASS_S(HomogeneousSmoothDiffuse, false, BSDF)
MTS_EXPORT_PLUGIN(HomogeneousSmoothDiffuse, "Smooth homogenousDiffuse BRDF")
MTS_NAMESPACE_END
