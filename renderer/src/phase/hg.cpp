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

#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/frame.h>

#include <stan/math.hpp>
#include <random>
#include <iostream>
#include <fstream>
#include <time.h>  
MTS_NAMESPACE_BEGIN

/*!\plugin{hg}{Henyey-Greenstein phase function}
 * \order{2}
 * \parameters{
 *     \parameter{g}{\Float}{
 *       This parameter must be somewhere in the range $-1$ to $1$
 *       (but not equal to $-1$ or $1$). It denotes the \emph{mean cosine}
 *       of scattering interactions. A value greater than zero indicates that
 *       medium interactions predominantly scatter incident light into a similar
 *       direction (i.e. the medium is \emph{forward-scattering}), whereas
 *       values smaller than zero cause the medium to be
 *       scatter more light in the opposite direction.
 *     }
 * }
 * This plugin implements the phase function model proposed by
 * Henyey and Greenstein \cite{Henyey1941Diffuse}. It is
 * parameterizable from backward- ($g<0$) through
 * isotropic- ($g=0$) to forward ($g>0$) scattering.
 */
class HGPhaseFunction : public PhaseFunction {
public:
	HGPhaseFunction(const Properties &props)
		: PhaseFunction(props) {

		int idx = 0;
		for (int i = 0; i < props.id_AD_Spectrum_AD.size(); i++) {
			if (props.id_AD_Spectrum_AD[i] == "g") {
				idx = i;
				break;
			}
		}
		m_g = props.m_theta[idx];
		if (m_g.val() >= 1 || m_g.val() <= -1)
			Log(EError, "The asymmetry parameter must lie in the interval (-1, 1)!");
	
	}

	HGPhaseFunction(Stream *stream, InstanceManager *manager)
		: PhaseFunction(stream, manager) {
		configure();
	}

	virtual ~HGPhaseFunction() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
	}

	void configure() {

		PhaseFunction::configure();
		m_type = EAngleDependence;
	}

	inline Float sample(PhaseFunctionSamplingRecord &pRec,
			Sampler *sampler) const {
		Point2 sample(sampler->next2D());

		Float cosTheta;
		if (std::abs(Float(m_g.val())) < Epsilon) {
			cosTheta = 1 - 2*sample.x;
		} else {
			Float sqrTerm = (1 - Float(m_g.val()) * Float(m_g.val())) / (1 - Float(m_g.val()) + 2 * Float(m_g.val()) * sample.x);
			cosTheta = (1 + Float(m_g.val()) * Float(m_g.val()) - sqrTerm * sqrTerm) / (2 * Float(m_g.val()));
		}

		Float sinTheta = math::safe_sqrt(1.0f-cosTheta*cosTheta),
			  sinPhi, cosPhi;

		math::sincos(2*M_PI*sample.y, &sinPhi, &cosPhi);

		pRec.wo = Frame(-pRec.wi).toWorld(Vector(
			sinTheta * cosPhi,
			sinTheta * sinPhi,
			cosTheta
		));

		return 1.0f;
	}

	inline Spectrum_AD sampleAD(PhaseFunctionSamplingRecord &pRec,
			Sampler *sampler) const {
		Point2 sample(sampler->next2D());

		Float cosTheta;
		if (std::abs(Float(m_g.val())) < Epsilon) {
			cosTheta = 1 - 2*sample.x;
		} else {
			Float sqrTerm = (1 - Float(m_g.val()) * Float(m_g.val())) / (1 - Float(m_g.val()) + 2 * Float(m_g.val()) * sample.x);
			cosTheta = (1 + Float(m_g.val()) * Float(m_g.val()) - sqrTerm * sqrTerm) / (2 * Float(m_g.val()));
		}

		Float sinTheta = math::safe_sqrt(1.0f-cosTheta*cosTheta),
			  sinPhi, cosPhi;

		math::sincos(2*M_PI*sample.y, &sinPhi, &cosPhi);

		pRec.wo = Frame(-pRec.wi).toWorld(Vector(
			sinTheta * cosPhi,
			sinTheta * sinPhi,
			cosTheta
		));

		Spectrum_AD val = HGPhaseFunction::eval(pRec);
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
			val[i] = val[i] / Float(val[i].val());
		}
		return val;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
			Float &pdf, Sampler *sampler) const {

		HGPhaseFunction::sample(pRec, sampler);
		pdf = Float(HGPhaseFunction::eval(pRec)[0].val());
		return 1.0f;
	}

	Spectrum_AD sampleAD(PhaseFunctionSamplingRecord &pRec,
			Float &pdf, Sampler *sampler) const {

		HGPhaseFunction::sample(pRec, sampler);
		Spectrum_AD val = HGPhaseFunction::eval(pRec);
		pdf = Float(val[0].val());
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
			val[i] = val[i] / Float(val[i].val());
		}
		
		return val;
	}

	Spectrum_AD eval(const PhaseFunctionSamplingRecord &pRec) const {
		Spectrum_AD temp, ret;
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
			temp[i] = 1.0f + m_g[i].m_data[Thread::getID()] * m_g[i].m_data[Thread::getID()]
						+ 2.0f * m_g[i].m_data[Thread::getID()] * dot(pRec.wi, pRec.wo);
			ret[i] = INV_FOURPI * (1.0f - m_g[i].m_data[Thread::getID()]*m_g[i].m_data[Thread::getID()])
								/ (temp[i] * stan::math::sqrt(temp[i]));
		}
		return ret;
	}

	FloatAD getMeanCosine() const {
		return m_g[0].m_data[Thread::getID()];
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HGPhaseFunction[g= [";
		for ( int i = 0; i < SPECTRUM_SAMPLES; ++i ) {
			if ( i ) oss << ", ";
			oss << m_g.val();
		}
		oss << "]" << endl;
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	Spectrum_shared m_g;
};

MTS_IMPLEMENT_CLASS_S(HGPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(HGPhaseFunction, "Henyey-Greenstein phase function");
MTS_NAMESPACE_END
