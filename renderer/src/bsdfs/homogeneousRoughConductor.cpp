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

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/texture.h>
#include <stan/math.hpp>
#include "microfacet.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

class HomogeneousRoughConductor : public BSDF {
public:
	HomogeneousRoughConductor(const Properties &props) : BSDF(props),
		m_type(MicrofacetDistribution::EGGX), m_sampleVisible(true) {
		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

		m_specularReflectance = props.getSpectrum("specularReflectance", Spectrum(1.0f));
		std::string materialName = props.getString("material", "Cu");

		Spectrum intEta, intK;
		if (boost::to_lower_copy(materialName) == "none") {
			intEta = Spectrum(0.0f);
			intK = Spectrum(1.0f);
		} else {
			intEta.fromContinuousSpectrum(InterpolatedSpectrum(
				fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
			intK.fromContinuousSpectrum(InterpolatedSpectrum(
				fResolver->resolve("data/ior/" + materialName + ".k.spd")));
		}
		Float extEta = lookupIOR(props, "extEta", "air");

		m_eta = props.getSpectrum("eta", intEta) / extEta;
		m_k   = props.getSpectrum("k", intK) / extEta;
		MicrofacetDistribution distr(props);
		m_type = distr.getType();
		m_sampleVisible = distr.getSampleVisible();
		if (props.hasProperty("distribution")) {
			std::string distr = boost::to_lower_copy(props.getString("distribution"));
			if (distr == "ggx")
				m_type = MicrofacetDistribution::EGGX;
			else
				SLog(EError, "Specified an invalid distribution \"%s\", must be "
					"\"ggx\"!", distr.c_str());
		}

		//m_sampleVisible = props.getBoolean("sampleVisible", true);

		int idx;
		for (int i = 0; i < props.id_AD_Spectrum_AD.size(); i++) {
			if (props.id_AD_Spectrum_AD[i] == "alpha") {
				idx = i;
				break;
			}
		}
		m_alpha = props.m_theta[idx];

	}

	HomogeneousRoughConductor(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_alpha = this->m_alpha;
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);
	}

	void configure() {
		unsigned int extraFlags = 0;

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);

		/* Verify the input parameters and fix them if necessary */
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
			m_specularReflectance[i] = (m_specularReflectance[i] < (Float) 1.0f) ? (m_specularReflectance[i]) : ((Float) 1.0f);
		}

		m_usesRayDifferentials = false;

		BSDF::configure();
	}

	/// Helper function: reflect \c wi with respect to a given surface normal
	inline Vector reflect(const Vector &wi, const Normal &m) const {
		return 2 * dot(wi, m) * Vector(m) - wi;
	}

	Spectrum_AD eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		/* Stop if this component was not requested */
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum_AD(0.0f);

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo+bRec.wi);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution_AD distr(
			m_type,
			m_alpha[0].m_data[Thread::getID()],
			m_sampleVisible
		);

		/* Evaluate the microfacet normal distribution */
		FloatAD D = distr.eval(H);

		/* Fresnel factor */
		const Spectrum_AD F = Spectrum_AD(fresnelConductorExact(dot(bRec.wi, H), m_eta, m_k) *
			m_specularReflectance);

		FloatAD G = distr.G(bRec.wi, bRec.wo, H);

		/* Calculate the total amount of reflection */
		FloatAD model = D * G / (4.0f * Frame::cosTheta(bRec.wi));
			
		return F * model;
	}



	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return 0.0f;

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo+bRec.wi);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution sampleDistrdistr(
			m_type,
			(Float) m_alpha.val(),
			(Float) m_alpha.val(),
			m_sampleVisible
		);

		if (m_sampleVisible) {

			return sampleDistrdistr.eval(H) * sampleDistrdistr.smithG1(bRec.wi, H)
				/ (4.0f * Frame::cosTheta(bRec.wi));
		}
		else {

			return sampleDistrdistr.pdf(bRec.wi, H) / (4 * absDot(bRec.wo, H));
		}
	}


	FloatAD pdfAD(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return FloatAD(0.0f);

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo+bRec.wi);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution_AD distr(
			m_type,
			m_alpha[0].m_data[Thread::getID()],
			m_sampleVisible
		);

		if (m_sampleVisible) {
			return distr.eval(H) * distr.smithG1(bRec.wi, H)
				/ (4.0f * Frame::cosTheta(bRec.wi));
		}
		else {

			return FloatAD(distr.pdf(bRec.wi, H) / (4 * absDot(bRec.wo, H)));
		}
	}







	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			(Float) m_alpha.val(),
			(Float) m_alpha.val(),
			m_sampleVisible
		);

		/* Sample M, the microfacet normal */
		Float pdf;
		Normal m = distr.sample(bRec.wi, sample, pdf);

		if (pdf == 0)
			return Spectrum(0.0f);

		/* Perfect specular reflection based on the microfacet normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		Spectrum F = fresnelConductorExact(dot(bRec.wi, m),
			m_eta, m_k) * m_specularReflectance;

		Float weight;
		if (m_sampleVisible) {
			weight = distr.smithG1(bRec.wo, m);
		} else {
			weight = distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (pdf * Frame::cosTheta(bRec.wi));
		}

		return F * weight;
	}

	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum_AD(0.0f);


		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution_AD distr(
			m_type,
			m_alpha[0].m_data[Thread::getID()],
			m_sampleVisible
		);

		/* Trick by Walter et al.: slightly scale the roughness values to
		   reduce importance sampling weights. Not needed for the
		   Heitz and D'Eon sampling technique. */
		MicrofacetDistribution sampleDistr(
				m_type,
				(Float) m_alpha.val(),
				(Float) m_alpha.val(),
				m_sampleVisible);

		/* Sample M, the microfacet normal */
		Float pdf;
		Normal m = sampleDistr.sample(bRec.wi, sample, pdf);

		if (pdf == 0)
			return Spectrum_AD(0.0f);

		/* Perfect specular reflection based on the microfacet normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum_AD(0.0f);

		Spectrum_AD F = Spectrum_AD(fresnelConductorExact(dot(bRec.wi, m),
			m_eta, m_k) * m_specularReflectance);

		Spectrum_AD weight(1.0f);
		if (m_sampleVisible) {
			weight *= distr.smithG1(bRec.wo, m);
		} else {
			// TODO: Check whether we need to make microfacetPDF AD here.
			weight *= distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (pdf * Frame::cosTheta(bRec.wi));
			}


		return F * weight;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			(Float) m_alpha.val(),
			(Float) m_alpha.val(),
			m_sampleVisible
		);

		/* Sample M, the microfacet normal */
		Normal m = distr.sample(bRec.wi, sample, pdf);

		if (pdf == 0)
			return Spectrum(0.0f);

		/* Perfect specular reflection based on the microfacet normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		Spectrum F = fresnelConductorExact(dot(bRec.wi, m),
			m_eta, m_k) * m_specularReflectance;

		Float weight;
		if (m_sampleVisible) {
			weight = distr.smithG1(bRec.wo, m);
		} else {
			weight = distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (pdf * Frame::cosTheta(bRec.wi));
		}

		/* Jacobian of the half-direction mapping */
		pdf /= 4.0f * dot(bRec.wo, m);

		return F * weight;
	}

	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum_AD(0.0f);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution_AD distr(
			m_type,
			m_alpha[0].m_data[Thread::getID()],
			m_sampleVisible
		);

		/* Trick by Walter et al.: slightly scale the roughness values to
		   reduce importance sampling weights. Not needed for the
		   Heitz and D'Eon sampling technique. */
		MicrofacetDistribution sampleDistr(
				m_type,
				(Float) m_alpha.val(),
				(Float) m_alpha.val(),
				m_sampleVisible);

		/* Sample M, the microfacet normal */
		Normal m = sampleDistr.sample(bRec.wi, sample, pdf);

		if (pdf == 0)
			return Spectrum_AD(0.0f);

		/* Perfect specular reflection based on the microfacet normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum_AD(0.0f);

		Spectrum_AD F = Spectrum_AD(fresnelConductorExact(dot(bRec.wi, m),
			m_eta, m_k) * m_specularReflectance);

		FloatAD weight;
		if (m_sampleVisible) {
			weight = distr.smithG1(bRec.wo, m);

		} else {
			// TODO: Check whether we need to make microfacetPDF AD here.
			weight = distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (pdf * Frame::cosTheta(bRec.wi));
					
		}
	
		/* Jacobian of the half-direction mapping */
		pdf /= 4.0f * dot(bRec.wo, m);

		return F * weight;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
				&& (name == "alpha" || name == "roughness")) {
			m_alpha = this->m_alpha;
		} else {
			BSDF::addChild(name, child);
		}
	}

	FloatAD getRoughness(const Intersection &its, int component) const {
		return m_alpha[0].m_data[Thread::getID()];
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughConductor[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
			<< "  sampleVisible = " << m_sampleVisible << "," << endl
			<< "  alpha = [";
			for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
				if ( i ) oss << ", ";
				oss << m_alpha.val();
			}
			oss << "]" << endl;
			oss << "  specularReflectance = [";
			for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
				if ( i ) oss << ", ";
				oss << m_specularReflectance[i];
			}
			oss << "]" << endl;
			oss << "  eta = " << m_eta.toString() << "," << endl
			<< "  k = " << m_k.toString() << endl
			<< "]";
		return oss.str();
	}

	// Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	MicrofacetDistribution::EType m_type;
	Spectrum m_specularReflectance;
	Spectrum_shared m_alpha;
	bool m_sampleVisible;
	Spectrum m_eta, m_k;
};

MTS_IMPLEMENT_CLASS_S(HomogeneousRoughConductor, false, BSDF)
MTS_EXPORT_PLUGIN(HomogeneousRoughConductor, "Homogeneous rough conductor BRDF");
MTS_NAMESPACE_END
