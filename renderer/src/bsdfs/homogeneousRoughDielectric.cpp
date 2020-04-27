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
#include <mitsuba/render/sampler.h>
#include <mitsuba/render/texture.h>
#include <stan/math.hpp>
#include "microfacet.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

class HomogeneousRoughDielectric : public BSDF {
public:
	HomogeneousRoughDielectric(const Properties &props) : BSDF(props),
		m_type(MicrofacetDistribution::EGGX), m_sampleVisible(true) {
		m_specularReflectance = props.getSpectrum("specularReflectance", Spectrum(1.0f));
		m_specularTransmittance = props.getSpectrum("specularTransmittance", Spectrum(1.0f));

		/* Specifies the internal index of refraction at the interface */
		Float intIOR = lookupIOR(props, "intIOR", "bk7");
		//FloatAD intIOR = props.m_theta_FloatAD[0]
		/* Specifies the external index of refraction at the interface */
		Float extIOR = lookupIOR(props, "extIOR", "air");

		if (intIOR < 0 || extIOR < 0 || intIOR == extIOR)
			Log(EError, "The interior and exterior indices of "
				"refraction must be positive and differ!");

		m_eta = intIOR / extIOR;
		m_invEta = 1 / m_eta;

		if (props.hasProperty("distribution")) {
			std::string distr = boost::to_lower_copy(props.getString("distribution"));
			if (distr == "ggx")
				m_type = MicrofacetDistribution::EGGX;
			else
				SLog(EError, "Specified an invalid distribution \"%s\", must be "
					"\"ggx\"!", distr.c_str());
		}

		m_sampleVisible = props.getBoolean("sampleVisible", true);
		
		int idx = 0;
		for (int i = 0; i < props.id_AD_Spectrum_AD.size(); i++) {
			if (props.id_AD_Spectrum_AD[i] == "alpha") {
				idx = i;
				break;
			}
		}
		
		
		m_alpha = props.m_theta[idx];
	}

	HomogeneousRoughDielectric(Stream *stream, InstanceManager *manager)
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
		m_components.push_back(EGlossyReflection | EFrontSide
			| EBackSide | EUsesSampler | extraFlags | 0);
		m_components.push_back(EGlossyTransmission | EFrontSide
			| EBackSide | EUsesSampler | ENonSymmetric | extraFlags | 0);

		/* Verify the input parameters and fix them if necessary */
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
			m_specularReflectance[i] = (m_specularReflectance[i] < (Float) 1.0f) ? (m_specularReflectance[i]) : ((Float) 1.0f);
			m_specularTransmittance[i] = (m_specularTransmittance[i] < (Float) 1.0f) ? (m_specularTransmittance[i]) : ((Float) 1.0f);
		}

		m_usesRayDifferentials = false;

		BSDF::configure();
	}

	Spectrum_AD eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle || Frame::cosTheta(bRec.wi) == 0)
			return Spectrum_AD(0.0f);

		/* Determine the type of interaction */
		bool reflect = Frame::cosTheta(bRec.wi)
			* Frame::cosTheta(bRec.wo) > 0;

		Vector H;
		if (reflect) {
			/* Stop if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 0)
				|| !(bRec.typeMask & EGlossyReflection))
				return Spectrum_AD(0.0f);

			/* Calculate the reflection half-vector */
			H = normalize(bRec.wo+bRec.wi);
		} else {
			/* Stop if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 1)
				|| !(bRec.typeMask & EGlossyTransmission))
				return Spectrum_AD(0.0f);

			/* Calculate the transmission half-vector */
			Float eta = Frame::cosTheta(bRec.wi) > 0
				? m_eta : m_invEta;

			H = normalize(bRec.wi + bRec.wo*eta);
		}

		/* Ensure that the half-vector points into the
		   same hemisphere as the macrosurface normal */
		H *= math::signum(Frame::cosTheta(H));

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution_AD distr(
			m_type,
			m_alpha[0].m_data[Thread::getID()],
			m_sampleVisible
		);

		/* Evaluate the microfacet normal distribution */
		const FloatAD D = distr.eval(H);

		/* Fresnel factor */
		const Float F = fresnelDielectricExt(dot(bRec.wi, H), m_eta);

		/* Smith's shadow-masking function */
		const FloatAD G = distr.G(bRec.wi, bRec.wo, H);

		if (reflect) {
			/* Calculate the total amount of reflection */
			FloatAD value = F * D * G /
				(4.0f * std::abs(Frame::cosTheta(bRec.wi)));

			return Spectrum_AD(m_specularReflectance) * value;
		} else {
			Float eta = Frame::cosTheta(bRec.wi) > 0.0f ? m_eta : m_invEta;

			/* Calculate the total amount of transmission */
			Float sqrtDenom = dot(bRec.wi, H) + eta * dot(bRec.wo, H);
			FloatAD value = ((1 - F) * D * G * eta * eta
				* dot(bRec.wi, H) * dot(bRec.wo, H)) /
				(Frame::cosTheta(bRec.wi) * sqrtDenom * sqrtDenom);

			/* Missing term in the original paper: account for the solid angle
			   compression when tracing radiance -- this is necessary for
			   bidirectional methods */
			Float factor = (bRec.mode == ERadiance)
				? (Frame::cosTheta(bRec.wi) > 0 ? m_invEta : m_eta) : 1.0f;

			return Spectrum_AD(m_specularTransmittance)
				* stan::math::abs(value * factor * factor);
		}
	}


	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle)
			return 0.0f;

		/* Determine the type of interaction */
		bool hasReflection   = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     reflect         = Frame::cosTheta(bRec.wi)
				             * Frame::cosTheta(bRec.wo) > 0;

		Vector H;
		Float dwh_dwo;

		if (reflect) {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 0)
				|| !(bRec.typeMask & EGlossyReflection))
				return 0.0f;

			/* Calculate the reflection half-vector */
			H = normalize(bRec.wo+bRec.wi);

			/* Jacobian of the half-direction mapping */
			dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, H));
		} else {
			/* Zero probability if this component was not requested */
			if ((bRec.component != -1 && bRec.component != 1)
				|| !(bRec.typeMask & EGlossyTransmission))
				return 0.0f;

			/* Calculate the transmission half-vector */
			Float eta = Frame::cosTheta(bRec.wi) > 0
				? m_eta : m_invEta;

			H = normalize(bRec.wi + bRec.wo*eta);

			/* Jacobian of the half-direction mapping */
			Float sqrtDenom = dot(bRec.wi, H) + eta * dot(bRec.wo, H);
			dwh_dwo = (eta*eta * dot(bRec.wo, H)) / (sqrtDenom*sqrtDenom);
		}

		/* Ensure that the half-vector points into the
		   same hemisphere as the macrosurface normal */
		H *= math::signum(Frame::cosTheta(H));

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution sampleDistr(
			m_type,
			(Float) m_alpha.val(),
			(Float) m_alpha.val(),
			m_sampleVisible
		);

		/* Trick by Walter et al.: slightly scale the roughness values to
		   reduce importance sampling weights. Not needed for the
		   Heitz and D'Eon sampling technique. */
		if (!m_sampleVisible)
			sampleDistr.scaleAlpha(1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		/* Evaluate the microfacet model sampling density function */
		Float prob = sampleDistr.pdf(math::signum(Frame::cosTheta(bRec.wi)) * bRec.wi, H);

		if (hasTransmission && hasReflection) {
			Float F = fresnelDielectricExt(dot(bRec.wi, H), m_eta);
			prob *= reflect ? F : (1-F);
		}

		return std::abs(prob * dwh_dwo);
	}

	FloatAD pdfAD(const BSDFSamplingRecord &bRec, EMeasure measure) const {

		return FloatAD(1.0);
	}





	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     sampleReflection = hasReflection;

		if (!hasReflection && !hasTransmission)
			return Spectrum(0.0f);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			(Float) m_alpha.val(),
			(Float) m_alpha.val(),
			m_sampleVisible
		);

		/* Trick by Walter et al.: slightly scale the roughness values to
		   reduce importance sampling weights. Not needed for the
		   Heitz and D'Eon sampling technique. */
		MicrofacetDistribution sampleDistr(distr);
		if (!m_sampleVisible)
			sampleDistr.scaleAlpha(1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		/* Sample M, the microfacet normal */
		Float microfacetPDF;
		const Normal m = sampleDistr.sample(math::signum(Frame::cosTheta(bRec.wi)) * bRec.wi, sample, microfacetPDF);
		if (microfacetPDF == 0)
			return Spectrum(0.0f);

		Float cosThetaT;
		Float F = fresnelDielectricExt(dot(bRec.wi, m), cosThetaT, m_eta);
		Spectrum weight(1.0f);

		if (hasReflection && hasTransmission) {
			if (bRec.sampler->next1D() > F)
				sampleReflection = false;
		} else {
			weight = Spectrum(hasReflection ? F : (1-F));
		}

		if (sampleReflection) {
			/* Perfect specular reflection based on the microfacet normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

			weight *= m_specularReflectance;
		} else {
			if (cosThetaT == 0)
				return Spectrum(0.0f);

			/* Perfect specular transmission based on the microfacet normal */
			bRec.wo = refract(bRec.wi, m, m_eta, cosThetaT);
			bRec.eta = cosThetaT < 0 ? m_eta : m_invEta;
			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum(0.0f);

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? m_invEta : m_eta) : 1.0f;

			weight *= m_specularTransmittance * (factor * factor);
		}

		if (m_sampleVisible)
			weight *= distr.smithG1(bRec.wo, m);
		else
			weight *= std::abs(distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (microfacetPDF * Frame::cosTheta(bRec.wi)));

		return weight;
	}

	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
			 hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
			 sampleReflection = hasReflection;

		if (!hasReflection && !hasTransmission)
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
		if (!m_sampleVisible)
			sampleDistr.scaleAlpha(1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		/* Sample M, the microfacet normal */
		Float microfacetPDF;
		const Normal m = sampleDistr.sample(math::signum(Frame::cosTheta(bRec.wi)) * bRec.wi, sample, microfacetPDF);
		if (microfacetPDF == 0)
			return Spectrum_AD(0.0f);

		Float cosThetaT;
		Float F = fresnelDielectricExt(dot(bRec.wi, m), cosThetaT, m_eta);
		Spectrum_AD weight(1.0f);

		if (hasReflection && hasTransmission) {
			if (bRec.sampler->next1D() > F)
				sampleReflection = false;
		} else {
			weight = Spectrum_AD(hasReflection ? F : (1-F));
		}

		if (sampleReflection) {
			/* Perfect specular reflection based on the microfacet normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum_AD(0.0f);

			weight *= m_specularReflectance;
		} else {
			if (cosThetaT == 0)
				return Spectrum_AD(0.0f);

			/* Perfect specular transmission based on the microfacet normal */
			bRec.wo = refract(bRec.wi, m, m_eta, cosThetaT);
			bRec.eta = cosThetaT < 0 ? m_eta : m_invEta;
			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum_AD(0.0f);

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? m_invEta : m_eta) : 1.0f;

			weight *= m_specularTransmittance * (factor * factor);
		}

		if (m_sampleVisible)
			weight *= distr.smithG1(bRec.wo, m);
		else
			// TODO: Check whether we need to make microfacetPDF AD here.
			weight *= stan::math::abs(distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (microfacetPDF * Frame::cosTheta(bRec.wi)));

		return weight;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
		     hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
		     sampleReflection = hasReflection;

		if (!hasReflection && !hasTransmission)
			return Spectrum(0.0f);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			(Float) m_alpha.val(),
			(Float) m_alpha.val(),
			m_sampleVisible
		);

		/* Trick by Walter et al.: slightly scale the roughness values to
		   reduce importance sampling weights. Not needed for the
		   Heitz and D'Eon sampling technique. */
		MicrofacetDistribution sampleDistr(distr);
		if (!m_sampleVisible)
			sampleDistr.scaleAlpha(1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		/* Sample M, the microfacet normal */
		Float microfacetPDF;
		const Normal m = sampleDistr.sample(math::signum(Frame::cosTheta(bRec.wi)) * bRec.wi, sample, microfacetPDF);
		if (microfacetPDF == 0)
			return Spectrum(0.0f);
		pdf = microfacetPDF;

		Float cosThetaT;
		Float F = fresnelDielectricExt(dot(bRec.wi, m), cosThetaT, m_eta);
		Spectrum weight(1.0f);

		if (hasReflection && hasTransmission) {
			if (bRec.sampler->next1D() > F) {
				sampleReflection = false;
				pdf *= 1-F;
			} else {
				pdf *= F;
			}
		} else {
			weight *= hasReflection ? F : (1-F);
		}

		Float dwh_dwo;
		if (sampleReflection) {
			/* Perfect specular reflection based on the microfacet normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

			weight *= m_specularReflectance;

			/* Jacobian of the half-direction mapping */
			dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, m));
		} else {
			if (cosThetaT == 0)
				return Spectrum(0.0f);

			/* Perfect specular transmission based on the microfacet normal */
			bRec.wo = refract(bRec.wi, m, m_eta, cosThetaT);
			bRec.eta = cosThetaT < 0 ? m_eta : m_invEta;
			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum(0.0f);

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? m_invEta : m_eta) : 1.0f;

			weight *= m_specularTransmittance * (factor * factor);

			/* Jacobian of the half-direction mapping */
			Float sqrtDenom = dot(bRec.wi, m) + bRec.eta * dot(bRec.wo, m);
			dwh_dwo = (bRec.eta*bRec.eta * dot(bRec.wo, m)) / (sqrtDenom*sqrtDenom);
		}

		if (m_sampleVisible)
			weight *= distr.smithG1(bRec.wo, m);
		else
			weight *= std::abs(distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (microfacetPDF * Frame::cosTheta(bRec.wi)));

		pdf *= std::abs(dwh_dwo);

		return weight;
	}

	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
		Point2 sample(_sample);

		bool hasReflection = ((bRec.component == -1 || bRec.component == 0)
							  && (bRec.typeMask & EGlossyReflection)),
			 hasTransmission = ((bRec.component == -1 || bRec.component == 1)
							  && (bRec.typeMask & EGlossyTransmission)),
			 sampleReflection = hasReflection;

		if (!hasReflection && !hasTransmission)
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
		if (!m_sampleVisible)
			sampleDistr.scaleAlpha(1.2f - 0.2f * std::sqrt(
				std::abs(Frame::cosTheta(bRec.wi))));

		/* Sample M, the microfacet normal */
		Float microfacetPDF;
		const Normal m = sampleDistr.sample(math::signum(Frame::cosTheta(bRec.wi)) * bRec.wi, sample, microfacetPDF);
		if (microfacetPDF == 0)
			return Spectrum_AD(0.0f);
		pdf = microfacetPDF;

		Float cosThetaT;
		Float F = fresnelDielectricExt(dot(bRec.wi, m), cosThetaT, m_eta);
		Spectrum_AD weight(1.0f);

		if (hasReflection && hasTransmission) {
			if (bRec.sampler->next1D() > F) {
				sampleReflection = false;
				pdf *= 1-F;
			} else {
				pdf *= F;
			}
		} else {
			weight *= hasReflection ? F : (1-F);
		}

		Float dwh_dwo;
		if (sampleReflection) {
			/* Perfect specular reflection based on the microfacet normal */
			bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum_AD(0.0f);

			weight *= m_specularReflectance;

			/* Jacobian of the half-direction mapping */
			dwh_dwo = 1.0f / (4.0f * dot(bRec.wo, m));
		} else {
			if (cosThetaT == 0)
				return Spectrum_AD(0.0f);

			/* Perfect specular transmission based on the microfacet normal */
			bRec.wo = refract(bRec.wi, m, m_eta, cosThetaT);
			bRec.eta = cosThetaT < 0 ? m_eta : m_invEta;
			bRec.sampledComponent = 1;
			bRec.sampledType = EGlossyTransmission;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0)
				return Spectrum_AD(0.0f);

			/* Radiance must be scaled to account for the solid angle compression
			   that occurs when crossing the interface. */
			Float factor = (bRec.mode == ERadiance)
				? (cosThetaT < 0 ? m_invEta : m_eta) : 1.0f;

			weight *= m_specularTransmittance * (factor * factor);

			/* Jacobian of the half-direction mapping */
			Float sqrtDenom = dot(bRec.wi, m) + bRec.eta * dot(bRec.wo, m);
			dwh_dwo = (bRec.eta*bRec.eta * dot(bRec.wo, m)) / (sqrtDenom*sqrtDenom);
		}

		if (m_sampleVisible)
			weight *= distr.smithG1(bRec.wo, m);
		else
			// TODO: Check whether we need to make microfacetPDF AD here.
			weight *= stan::math::abs(distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (microfacetPDF * Frame::cosTheta(bRec.wi)));

		pdf *= std::abs(dwh_dwo);

		return weight;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
				&& (name == "alpha" || name == "roughness")) {
			m_alpha = this->m_alpha;
		} else {
			BSDF::addChild(name, child);
		}
	}

	Float getEta() const {
		return m_eta;
	}

	FloatAD getRoughness(const Intersection &its, int component) const {
		return m_alpha[0].m_data[Thread::getID()];
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "homogemesouRoughDielectric[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
			<< "  sampleVisible = " << m_sampleVisible << "," << endl
			<< "  eta = " << m_eta << "," << endl
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
		oss << "  specularTransmittance = [";
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
			if ( i ) oss << ", ";
			oss << m_specularTransmittance[i];
		}
		oss << "]" << endl;
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	MicrofacetDistribution::EType m_type;
	Spectrum m_specularTransmittance;
	Spectrum m_specularReflectance;
	Spectrum_shared m_alpha;
	Float m_eta, m_invEta;
	bool m_sampleVisible;
};

MTS_IMPLEMENT_CLASS_S(HomogeneousRoughDielectric, false, BSDF)
MTS_EXPORT_PLUGIN(HomogeneousRoughDielectric, "Rough homogeneousDielectric BSDF");
MTS_NAMESPACE_END
