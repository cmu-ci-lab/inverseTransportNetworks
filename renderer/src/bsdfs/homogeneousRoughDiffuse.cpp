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
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{roughdiffuse}{Rough diffuse material}
 * \order{2}
 * \icon{bsdf_roughdiffuse}
 * \parameters{
 *     \parameter{reflectance}{\Spectrum\Or\Texture}{
 *         Specifies the diffuse albedo of the
 *         material. \default{0.5}
 *     }
 *     \parameter{alpha}{\Spectrum\Or\Texture}{
 *         Specifies the roughness of the unresolved surface micro-geometry
 *         using the \emph{root mean square} (RMS) slope of the
 *         microfacets. \default{0.2}
 *     }
 *     \parameter{useFastApprox}{\Boolean}{
 *         This parameter selects between the full version of the model
 *         or a fast approximation that still retains most qualitative features.
 *         \default{\texttt{false}, i.e. use the high-quality version}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Smooth diffuse surface ($\alpha=0$)}
 *        {bsdf_roughdiffuse_0}
 *     \rendering{Very rough diffuse surface ($\alpha=0.7$)}
 *        {bsdf_roughdiffuse_0_7}
 *   \vspace{-3mm}
 *   \caption{The effect of switching from smooth to rough diffuse scattering
 *   is fairly subtle on this model---generally, there will be higher
 *   reflectance at grazing angles, as well as an overall reduced contrast.}\vspace{3mm}
 * }
 *
 * This reflectance model describes the interaction of light with a \emph{rough}
 * diffuse material, such as plaster, sand, clay, or concrete, or ``powdery''
 * surfaces. The underlying theory was developed by Oren and Nayar
 * \cite{Oren1994Generalization}, who model the microscopic surface structure as
 * unresolved planar facets arranged in V-shaped grooves, where each facet
 * is an ideal diffuse reflector. The model takes into account shadowing,
 * masking, as well as interreflections between the facets.
 *
 * Since the original publication, this approach has been shown to
 * be a good match for many real-world materials, particularly compared
 * to Lambertian scattering, which does not take surface roughness into account.
 *
 * The implementation in Mitsuba uses a surface roughness parameter $\alpha$ that
 * is slightly different from the slope-area variance in the original 1994 paper.
 * The reason for this change is to make the parameter $\alpha$ portable
 * across different models (i.e. \pluginref{roughdielectric}, \pluginref{roughplastic},
 * \pluginref{roughconductor}).
 *
 * To get an intuition about the effect of the
 * parameter $\alpha$, consider the following approximate classification:
 * a value of $\alpha=0.001-0.01$ corresponds to a material
 * with slight imperfections on an otherwise smooth surface (for such small
 * values, the model will behave identically to \pluginref{diffuse}), $\alpha=0.1$
 * is relatively rough, and $\alpha=0.3-0.7$ is \emph{extremely} rough
 * (e.g. an etched or ground surface).
 *
 * Note that this material is one-sided---that is, observed from the
 * back side, it will be completely black. If this is undesirable,
 * consider using the \pluginref{twosided} BRDF adapter plugin.
 */
class homogeneousRoughDiffuse : public BSDF {
public:
	homogeneousRoughDiffuse(const Properties &props)
		: BSDF(props) {
		/* For better compatibility with other models, support both
		   'reflectance' and 'diffuseReflectance' as parameter names */
		m_reflectance = props.m_theta[0];
		m_useFastApprox = props.getBoolean("useFastApprox", false);
		m_alpha = props.m_theta_FloatAD[0];

	}

	homogeneousRoughDiffuse(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_reflectance = this->m_reflectance;
		m_alpha = this->m_alpha;
		m_useFastApprox = stream->readBool();

		configure();
	}

	void configure() {
		/* Verify the input parameter and fix them if necessary */
		//m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide
			| 0);
		BSDF::configure();
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		Spectrum ret;
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) ret[i] = m_reflectance[i].val();
		return ret;
	}
	Spectrum_AD getRef() const {
		return m_reflectance;
	}
	Spectrum_AD eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum_AD(0.0f);

		/* Conversion from Beckmann-style RMS roughness to
		   Oren-Nayar-style slope-area variance. The factor
		   of 1/sqrt(2) was found to be a perfect fit up
		   to extreme roughness values (>.5), after which
		   the match is not as good anymore */
		const Float conversionFactor = 1 / std::sqrt((Float) 2);

		FloatAD sigma = m_alpha * conversionFactor;
		const FloatAD sigma2 = sigma*sigma;

		Float sinThetaI = Frame::sinTheta(bRec.wi),
			  sinThetaO = Frame::sinTheta(bRec.wo);

		Float cosPhiDiff = 0;
		if (sinThetaI > Epsilon && sinThetaO > Epsilon) {
			/* Compute cos(phiO-phiI) using the half-angle formulae */
			Float sinPhiI = Frame::sinPhi(bRec.wi),
				  cosPhiI = Frame::cosPhi(bRec.wi),
				  sinPhiO = Frame::sinPhi(bRec.wo),
				  cosPhiO = Frame::cosPhi(bRec.wo);
			cosPhiDiff = cosPhiI * cosPhiO + sinPhiI * sinPhiO;
		}

		if (m_useFastApprox) {
			FloatAD A = 1.0f - 0.5f * sigma2 / (sigma2 + 0.33f),
				    B = 0.45f * sigma2 / (sigma2 + 0.09f),
				    sinAlpha, tanBeta;

			if (Frame::cosTheta(bRec.wi) > Frame::cosTheta(bRec.wo)) {
				sinAlpha = sinThetaO;
				tanBeta = sinThetaI / Frame::cosTheta(bRec.wi);
			} else {
				sinAlpha = sinThetaI;
				tanBeta = sinThetaO / Frame::cosTheta(bRec.wo);
			}

			return m_reflectance * (INV_PI * Frame::cosTheta(bRec.wo) * (A + B
				* std::max(cosPhiDiff, (Float) 0.0f) * sinAlpha * tanBeta));
		} else {
			Float sinThetaI = Frame::sinTheta(bRec.wi),
				  sinThetaO = Frame::sinTheta(bRec.wo),
				  thetaI = math::safe_acos(Frame::cosTheta(bRec.wi)),
				  thetaO = math::safe_acos(Frame::cosTheta(bRec.wo)),
				  alpha = std::max(thetaI, thetaO),
				  beta = std::min(thetaI, thetaO);

			Float sinAlpha, sinBeta, tanBeta;
			if (Frame::cosTheta(bRec.wi) > Frame::cosTheta(bRec.wo)) {
				sinAlpha = sinThetaO; sinBeta = sinThetaI;
				tanBeta = sinThetaI / Frame::cosTheta(bRec.wi);
			} else {
				sinAlpha = sinThetaI; sinBeta = sinThetaO;
				tanBeta = sinThetaO / Frame::cosTheta(bRec.wo);
			}

			FloatAD tmp = sigma2 / (sigma2 + 0.09f),
				    tmp2 = (4*INV_PI*INV_PI) * alpha * beta,
				    tmp3 = 2*beta*INV_PI;

			FloatAD C1 = 1.0f - 0.5f * sigma2 / (sigma2 + 0.33f),
				    C2 = 0.45f * tmp,
				    C3 = 0.125f * tmp * tmp2 * tmp2,
				    C4 = 0.17f * sigma2 / (sigma2 + 0.13f);

			if (cosPhiDiff > 0)
				C2 *= sinAlpha;
			else
				C2 *= sinAlpha - tmp3*tmp3*tmp3;

			/* Compute tan(0.5 * (alpha+beta)) using the half-angle formulae */
			Float tanHalf = (sinAlpha + sinBeta) / (
					math::safe_sqrt(1.0f - sinAlpha * sinAlpha) +
					math::safe_sqrt(1.0f - sinBeta  * sinBeta));

			Spectrum_AD rho = m_reflectance,
					    snglScat = rho * (C1 + cosPhiDiff * C2 * tanBeta +
					    (1.0f - std::abs(cosPhiDiff)) * C3 * tanHalf),
					    dblScat = rho * rho * (C4 * (1.0f - cosPhiDiff*tmp3*tmp3));

			return  (snglScat + dblScat) * (INV_PI * Frame::cosTheta(bRec.wo));
		}
	}
	std::vector<stan::math::var> getDiff(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		std::vector<stan::math::var> theta;
		return theta;
	}	
	stan::math::var eval_test() const {
		stan::math::var a = 1.0f;
		return  a;
	}

	std::pair<stan::math::var, stan::math::var> eval_test2() const {
		stan::math::var a = 1.0f;
		std::pair<stan::math::var, stan::math::var> ret;
		ret.first = a;
		ret.second = a;

		return  ret;
	}	
	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (!(bRec.typeMask & EGlossyReflection) || measure != ESolidAngle
			|| Frame::cosTheta(bRec.wi) <= 0
			|| Frame::cosTheta(bRec.wo) <= 0)
			return 0.0f;

		return warp::squareToCosineHemispherePdf(bRec.wo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (!(bRec.typeMask & EGlossyReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		bRec.wo = warp::squareToCosineHemisphere(sample);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;
		return Spectrum(eval(bRec, ESolidAngle)[0].val()) /
			warp::squareToCosineHemispherePdf(bRec.wo);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (!(bRec.typeMask & EGlossyReflection) || Frame::cosTheta(bRec.wi) <= 0)
			return Spectrum(0.0f);

		bRec.wo = warp::squareToCosineHemisphere(sample);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;
		pdf = warp::squareToCosineHemispherePdf(bRec.wo);
		return Spectrum(eval(bRec, ESolidAngle)[0].val()) / pdf;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "reflectance" || name == "diffuseReflectance")
				m_reflectance = this->m_reflectance;
			else if (name == "alpha")
				m_alpha = this->m_alpha;
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		//manager->serialize(stream, m_reflectance.get());
		//manager->serialize(stream, m_alpha.get());
		stream->writeBool(m_useFastApprox);
	}

	Float getRoughness(const Intersection &its, int component) const {
		return std::numeric_limits<Float>::infinity();
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "RoughDiffuse[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  reflectance = " << m_reflectance[0].val() << "," << endl
			<< "  alpha = " << m_alpha << "," << endl
			<< "  useFastApprox = " << m_useFastApprox << endl
			<< "]";
		return oss.str();
	}



	MTS_DECLARE_CLASS()
private:
	Spectrum_AD m_reflectance;
	FloatAD m_alpha;
	bool m_useFastApprox;
};



MTS_IMPLEMENT_CLASS_S(homogeneousRoughDiffuse, false, BSDF)
MTS_EXPORT_PLUGIN(homogeneousRoughDiffuse, "Rough Homogeneous diffuse BRDF")
MTS_NAMESPACE_END
