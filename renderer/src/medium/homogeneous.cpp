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
#include <stan/math.hpp>
#include "maxexp.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{homogeneous}{Homogeneous participating medium}
 * \order{1}
 * \parameters{
 *     \parameter{material}{\String}{
 *         Name of a material preset, see
 *         \tblref{medium-coefficients}. \default{\texttt{skin1}}
 *     }
 *     \parameter{sigmaA, sigmaS}{\Spectrum}{
 *         Absorption and scattering
 *         coefficients of the medium in inverse scene units.
 *         These parameters are mutually exclusive with \code{sigmaT} and \code{albedo}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{sigmaT, albedo}{\Spectrum}{
 *         Extinction coefficient in inverse scene units
 *         and a (unitless) single-scattering albedo.
 *         These parameters are mutually exclusive with \code{sigmaA} and \code{sigmaS}
 *         \default{configured based on \code{material}}
 *     }
 *     \parameter{\footnotesize{scale}}{\Float}{
 *         Optional scale factor that will be applied to the \code{sigma*} parameters.
 *         It is provided for convenience when accomodating data based on different units,
 *         or to simply tweak the density of the medium. \default{1}
 *     }
 *     \parameter{\Unnamed}{\Phase}{
 *          A nested phase function that describes the directional
 *          scattering properties of the medium. When none is specified,
 *          the renderer will automatically use an instance of
 *          \pluginref{isotropic}.
 *     }
 * }
 *
 * This class implements a flexible homogeneous participating
 * medium with support for arbitrary phase functions and various
 * medium sampling methods. It provides two different ways of configuring
 * the medium properties. One possibility is to load a material preset
 * using the \code{material} parameter---see \tblref{medium-coefficients}
 * for details. Alternatively, when specifying parameters by hand, they can either
 * be provided using the scattering and absorption coefficients, or
 * by declaring the extinction coefficient and single scattering
 * albedo (whichever is more convenient). Mixing these parameter
 * initialization methods is not allowed.
 *
 * All scattering parameters (named \code{sigma*}) should
 * be provided in inverse scene units. For instance, when a world-space
 * distance of 1 unit corresponds to a meter, the scattering coefficents should
 * have units of inverse meters. For convenience, the \code{scale}
 * parameter can be used to correct the units. For instance, when the scene is
 * in meters and the coefficients are in inverse millimeters, set
 * \code{scale} to \code{1000}.
 *
 * \renderings{
 *    \rendering{A squishy ball rendered with subsurface scattering and
 *    a dielectric BSDF (courtesy of Chanxi Zheng)}{medium_homogeneous_squishy.jpg}
 * }
 *
 * \begin{xml}[caption=Declaration of a forward scattering medium with high albedo]
 * <medium id="myMedium" type="homogeneous">
 *     <spectrum name="sigmaS" value="1"/>
 *     <spectrum name="sigmaA" value="0.05"/>
 *
 *     <phase type="hg">
 *         <float name="g" value="0.7"/>
 *     </phase>
 * </medium>
 * \end{xml}
 *
 * \textbf{Note}: Rendering media that have a spectrally
 * varying extinction coefficient can be tricky, since all
 * commonly used medium sampling methods suffer from high
 * variance in that case. Here, it may often make more sense to render
 * several monochromatic images separately (using only the coefficients for
 * a single channel) and then merge them back into a RGB image. There
 * is a \code{mtsutil} (\secref{mtsutil}) plugin named \code{joinrgb}
 * that will perform this RGB merging process.
 *
 * \begin{table}[h!]
 *     \centering
 *     \vspace{3mm}
 *     {\footnotesize
 *     \begin{tabular}{>{\ttfamily}p{3.8cm}p{.4cm}>{\ttfamily}p{3.8cm}p{.4cm}>{\ttfamily}p{3.8cm}}
 *         \toprule
 *         \rmfamily \small\textbf{Name} &&
 *         \rmfamily \small\textbf{Name} &&
 *         \rmfamily \small\textbf{Name} \\
 *         \cmidrule{1-1} \cmidrule{3-3} \cmidrule{5-5}
 *         Apple && Chicken1 && Chicken2 \\
 *         Cream && Ketchup  && Potato \\
 *         Skimmilk && Skin1 && Skin2 \\
 *         Spectralon && Wholemilk && \\
 *         \cmidrule{1-1} \cmidrule{3-3} \cmidrule{5-5}
 *         Lowfat Milk              &&  Gatorade                &&    White Grapefruit Juice     \\
 *         Reduced Milk             &&  Chardonnay              &&    Shampoo                    \\
 *         Regular Milk             &&  White Zinfandel         &&    Strawberry Shampoo         \\
 *         Espresso                 &&  Merlot                  &&    \mbox{Head \& Shoulders Shampoo}  \\
 *         Mint Mocha Coffee        &&  Budweiser Beer          &&    Lemon Tea Powder           \\
 *         Lowfat Soy Milk          &&  Coors Light Beer        &&    Orange Juice Powder        \\
 *         Regular Soy Milk         &&  Clorox                  &&    Pink Lemonade Powder       \\
 *         Lowfat Chocolate Milk    &&  Apple Juice             &&    Cappuccino Powder          \\
 *         Regular Chocolate Milk   &&  Cranberry Juice         &&    Salt Powder                \\
 *         Coke                     &&  Grape Juice             &&    Sugar Powder               \\
 *         Pepsi                    &&  Ruby Grapefruit Juice   &&    Suisse Mocha               \\
 *         Sprite                   &&                          &&                               \\
 *         \bottomrule
 *     \end{tabular}}
 *     \caption{\label{tbl:medium-coefficients}This
 *          table lists all supported medium material presets. The
 *          top entries are from Jensen et al. \cite{Jensen2001Practical}, and the
 *          bottom ones are from Narasimhan et al. \cite{Narasimhan2006Acquiring}.
 *          They all use units of $\frac{1}{mm}$, so remember to set
 *          \code{scale} appropriately when your scene is not
 *          in units of millimeters.
 *          These material presets can be used with the plugins
 *          \pluginref{homogeneous},\
 *          \pluginref{dipole}, and \
 *          \pluginref{hk}
 *     }
 * \end{table}
 */



class HomogeneousMedium : public Medium {
public:
	/**
	 * This class supports the following sampling strategies for choosing
	 * a suitable scattering location when sampling the RTE
	 */
	enum ESamplingStrategy {
		EBalance,  /// Exponential distrib.; pick a random channel each time
		ESingle,   /// Exponential distrib.; pick a specified channel
		EManual,   /// Exponential distrib.; manually specify the falloff
		EMaximum   /// Maximum-of-exponential distribution
	};

	HomogeneousMedium(const Properties &props)
			: Medium(props), m_samplingDensity(0.0f), m_maxExpDist(NULL) {
		std::string strategy = props.getString("strategy", "balance");

		/**
		 * The goal of the medium sampling weight is to be able to
		 * sample medium intarctions according to
		 *    sigma_s(t) * tau(0 <-> t)
		 * as opposed to
		 *    sigma_t(t) * tau(0 <-> t)
		 * See the separate writeup for more details.
		 */
		m_mediumSamplingWeight = props.getFloat("mediumSamplingWeight", -1);
		if (m_mediumSamplingWeight == -1) {
			for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
				/// Record the highest albedo values across channels
				Float albedo = m_albedo[i].m_data[Thread::getID()].val();
				if (albedo > m_mediumSamplingWeight && m_sigmaT[i].m_data[Thread::getID()].val() != 0)
					m_mediumSamplingWeight = albedo;
			}
			if (m_mediumSamplingWeight > 0) {
				/* The medium scatters some light -> place at least half
				   of the samples in it, otherwise we will render lots
				   of spatially varying noise where one pixel has a
				   medium interaction and the neighbors don't */
				m_mediumSamplingWeight = std::max(m_mediumSamplingWeight,
					(Float) 0.5f);
			}
		}

		if (strategy == "balance") {
			m_strategy = EBalance;
		} else if (strategy == "single") {
			m_strategy = ESingle;

			/* By default, choose the lowest-variance channel
			   (the one with the smallest sigma_t, that is) */
			int channel = 0;
			Float smallest = std::numeric_limits<Float>::infinity();
			for (int i=0; i<SPECTRUM_SAMPLES; ++i) {
				if (m_sigmaT[i].m_data[Thread::getID()].val() < smallest) {
					smallest = m_sigmaT[i].m_data[Thread::getID()].val();
					channel = i;
				}
			}

			channel = props.getInteger("channel", channel);
			Assert(channel >= 0 && channel < SPECTRUM_SAMPLES);
			m_samplingDensity = m_sigmaT[channel].m_data[Thread::getID()].val();
		
//			if (props.getBoolean("monochromatic", false)) {
//				/* Optionally turn this into a monochromatic medium
//				   based on the chosen color channel. This is useful
//				   when the whole scene is rendered once per channel
//				   and then recombined to create a color image */
//				m_sigmaA = Spectrum(m_sigmaA[channel]);
//				m_sigmaS = Spectrum(m_sigmaS[channel]);
//				m_sigmaT = m_sigmaA + m_sigmaS;
//			}
		} else if (strategy == "maximum") {
			m_strategy = EMaximum;
			std::vector<Float> coeffs(SPECTRUM_SAMPLES);
			for (int i=0; i<SPECTRUM_SAMPLES; ++i)
				coeffs[i] = m_sigmaT[i].m_data[Thread::getID()].val();
			m_maxExpDist = new MaxExpDist(coeffs);
		} else if (strategy == "manual") {
			m_strategy = EManual;
			m_samplingDensity = props.getFloat("samplingDensity");
		} else {
			Log(EError, "Specified an unknown sampling strategy");
		}
	}

	HomogeneousMedium(Stream *stream, InstanceManager *manager)
		: Medium(stream, manager), m_maxExpDist(NULL) {
	}

	virtual ~HomogeneousMedium() {
		if (m_maxExpDist)
			delete m_maxExpDist;
	}

	void configure() {
		Medium::configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Medium::serialize(stream, manager);
	}

	Spectrum_AD evalTransmittance(const Ray &ray, Sampler *) const {
		Float negLength = ray.mint - ray.maxt;
		Spectrum_AD transmittance;
		if (std::isinf(negLength)) {
			transmittance[0] = FloatAD(0.0f);
		} else {
			transmittance[0] = stan::math::exp(m_sigmaT[0].m_data[Thread::getID()] * negLength);
		}
//			transmittance[i] = m_sigmaT[i] != 0
//				? math::fastexp(m_sigmaT[i] * negLength) : (Float) 1.0f;
		return transmittance;
	}

	bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec,
			Sampler *sampler) const {
		Float rand = sampler->next1D(), sampledDistance;

		Float samplingDensity = m_samplingDensity;

		if (rand < m_mediumSamplingWeight) {
			rand /= m_mediumSamplingWeight;
			if (m_strategy == EBalance) {
				
				int channel = std::min((int) (sampler->next1D()
						* SPECTRUM_SAMPLES), SPECTRUM_SAMPLES-1);
				samplingDensity = Float(m_sigmaT[channel].m_data[Thread::getID()].val());
			}
			sampledDistance = -math::fastlog(1-rand) / samplingDensity;
		} else {
			/* Don't generate a medium interaction */
			sampledDistance = std::numeric_limits<Float>::infinity();
		}
		Float distSurf = ray.maxt - ray.mint;
		bool success = true;
		
		if (sampledDistance < distSurf) {
			mRec.t = sampledDistance + ray.mint;
			mRec.p = ray(mRec.t);
			
			mRec.sigmaA = getSigmaA();
			mRec.sigmaS = getSigmaS();
//			for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
//				mRec.sigmaA[i] = m_sigmaA[i].m_data[Thread::getID()];
//				mRec.sigmaS[i] = m_sigmaS[i].m_data[Thread::getID()];
//			}
			mRec.time = ray.time;
			mRec.medium = this;

			/* Fail if there is no forward progress
			   (e.g. due to roundoff errors) */
			if (mRec.p == ray.o)
				success = false;
		} else {
			sampledDistance = distSurf;
			success = false;
		}

		switch (m_strategy) {
			case EBalance:
				mRec.pdfFailure = 0;
				mRec.pdfSuccess = 0;
				Float tmp;
				if (std::isinf(sampledDistance)) {
					tmp = Float(0.0f);
				} else {
					tmp = math::fastexp(-Float(m_sigmaT[0].m_data[Thread::getID()].val()) * sampledDistance);
				}
				mRec.pdfFailure += tmp;
				mRec.pdfSuccess += Float(m_sigmaT[0].m_data[Thread::getID()].val()) * tmp;
				break;

			default:
				Log(EError, "Unknown sampling strategy!");
		}

		if (std::isinf(sampledDistance)) {
			mRec.transmittance[0] = FloatAD(0.0f);
		} else {
			mRec.transmittance[0] = stan::math::exp(m_sigmaT[0].m_data[Thread::getID()] * (-sampledDistance));
		}
		mRec.pdfSuccessRev = mRec.pdfSuccess = mRec.pdfSuccess * m_mediumSamplingWeight;
		mRec.pdfFailure = m_mediumSamplingWeight * mRec.pdfFailure + (1-m_mediumSamplingWeight);
		mRec.medium = this;
//		if (mRec.transmittance.max() < 1e-20)
//			mRec.transmittance = Spectrum(0.0f);

		return success;
	}

	void eval(const Ray &ray, MediumSamplingRecord &mRec) const {
		Float distance = ray.maxt - ray.mint;
		switch (m_strategy) {
			case EBalance: {
					mRec.pdfSuccess = 0;
					mRec.pdfFailure = 0;
					Float temp;
					if (std::isinf(distance)) {
						temp = Float(0.0f);
					} else {
						temp = math::fastexp(-Float(m_sigmaT[0].m_data[Thread::getID()].val()) * distance);
					}
					mRec.pdfSuccess += Float(m_sigmaT[0].m_data[Thread::getID()].val()) * temp;
					mRec.pdfFailure += temp;
				}
				break;

			default:
				Log(EError, "Unknown sampling strategy!");
		}

		if (std::isinf(distance)) {
			mRec.transmittance[0] = FloatAD(0.0f);
		} else {
			mRec.transmittance[0] = stan::math::exp(m_sigmaT[0].m_data[Thread::getID()] * (-distance));
		}
		mRec.pdfSuccess = mRec.pdfSuccessRev = mRec.pdfSuccess * m_mediumSamplingWeight;
		mRec.pdfFailure = mRec.pdfFailure * m_mediumSamplingWeight + (1-m_mediumSamplingWeight);
		mRec.sigmaA = getSigmaA();
		mRec.sigmaS = getSigmaS();
//		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
//			mRec.sigmaA[i] = m_sigmaA[i].m_data[Thread::getID()];
//			mRec.sigmaS[i] = m_sigmaS[i].m_data[Thread::getID()];
//		}

		mRec.time = ray.time;
		mRec.medium = this;
//		if (mRec.transmittance.max() < 1e-20)
//			mRec.transmittance = Spectrum(0.0f);
	}

	bool isHomogeneous() const {
		return true;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HomogeneousMedium[" << endl
			<< "  albedo = " << m_albedo.toString() << "," << endl
			<< "  sigmaT = " << m_sigmaT.toString() << "," << endl
			<< "  mediumSamplingWeight = " << m_mediumSamplingWeight << "," << endl
			<< "  samplingDensity = " << m_samplingDensity << "," << endl
			<< "  strategy = ";

		switch (m_strategy) {
			case ESingle: oss << "single," << endl; break;
			case EManual: oss << "manual," << endl; break;
			case EBalance: oss << "balance," << endl; break;
			case EMaximum: oss << "maximum," << endl; break;
		}

		oss << "  phase = " << indent(m_phaseFunction.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Float m_samplingDensity, m_mediumSamplingWeight;
	ESamplingStrategy m_strategy;
	MaxExpDist *m_maxExpDist;
};


MTS_IMPLEMENT_CLASS_S(HomogeneousMedium, false, Medium)
MTS_EXPORT_PLUGIN(HomogeneousMedium, "Homogeneous medium");
MTS_NAMESPACE_END
