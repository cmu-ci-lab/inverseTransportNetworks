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
#include <mitsuba/hw/gpuprogram.h>

MTS_NAMESPACE_BEGIN

/*! \plugin{mixturebsdf}{Mixture material}
 * \order{16}
 * \parameters{
 *     \parameter{weights}{\String}{A comma-separated list of BSDF weights}
 *     \parameter{\Unnamed}{\BSDF}{Multiple BSDF instances that should be
 *     mixed according to the specified weights}
 * }
 * \renderings{
 *     \medrendering{Smooth glass}{bsdf_mixturebsdf_smooth}
 *     \medrendering{Rough glass}{bsdf_mixturebsdf_rough}
 *     \medrendering{An mixture of 70% smooth glass and 30% rough glass
 *     results in a more realistic smooth material with imperfections
 *     (\lstref{mixture-example})}{bsdf_mixturebsdf_result}
 * }
 *
 * This plugin implements a ``mixture'' material, which represents
 * linear combinations of multiple BSDF instances. Any surface scattering
 * model in Mitsuba (be it smooth, rough, reflecting, or transmitting) can
 * be mixed with others in this manner to synthesize new models. There
 * is no limit on how many models can be mixed, but their combination
 * weights must be non-negative and sum to a value of one or less to ensure
 * energy balance. When they sum to less than one, the material will
 * absorb a proportional amount of the incident illlumination.
 *
 * \vspace{4mm}
 * \begin{xml}[caption={A material definition for a mixture of 70% smooth
 *     and 30% rough glass},
 *     label=lst:mixture-example]
 * <bsdf type="mixturebsdf">
 *     <string name="weights" value="0.7, 0.3"/>
 *
 *     <bsdf type="dielectric"/>
 *
 *     <bsdf type="roughdielectric">
 *         <float name="alpha" value="0.3"/>
 *     </bsdf>
 * </bsdf>
 * \end{xml}
 */

class MixtureBSDFAD : public BSDF {
public:
	MixtureBSDFAD(const Properties &props)
		: BSDF(props) {
		/* Parse the weight parameter */
		for (int iIdx = 0; iIdx < props.id_AD_FloatAD.size(); iIdx++) {
			if (props.id_AD_FloatAD[iIdx].find("w") != std::string::npos) {
				m_weights.push_back(props.m_theta_FloatAD[iIdx]);

			}
		}
		testalpha = props.m_theta[0];

	}

	MixtureBSDFAD(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {

		configure();
	}

	virtual ~MixtureBSDFAD() {
		for (size_t i=0; i<m_bsdfs.size(); ++i)
			m_bsdfs[i]->decRef();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);


	}

	void configure() {
		m_usesRayDifferentials = false;
		size_t componentCount = 0;
		if (m_bsdfs.size() != m_weights.size())
			Log(EError, "BSDF count mismatch: " SIZE_T_FMT " bsdfs, but specified " SIZE_T_FMT " weights",
				m_bsdfs.size(), m_bsdfs.size());

		Float totalWeight = 0;
		for (size_t i=0; i<m_weights.size(); ++i)
			totalWeight += m_weights[i].val();

		if (totalWeight <= 0)
			Log(EError, "The weights must sum to a value greater than zero!");

		if (m_ensureEnergyConservation && totalWeight > 1) {
			std::ostringstream oss;
			Float scale = 1.0f / totalWeight;
			oss << "The BSDF" << endl << toString() << endl
				<< "potentially violates energy conservation, since the weights "
				<< "sum to " << totalWeight << ", which is greater than one! "
				<< "They will be re-scaled to avoid potential issues. Specify "
				<< "the parameter ensureEnergyConservation=false to prevent "
				<< "this from happening.";
			Log(EWarn, "%s", oss.str().c_str());
			for (size_t i=0; i<m_weights.size(); ++i)
				m_weights[i] *= scale;
		}

		for (size_t i=0; i<m_bsdfs.size(); ++i) {
				componentCount += m_bsdfs[i]->getComponentCount();
		}


		m_pdf = DiscreteDistribution(m_bsdfs.size());
		
		m_components.reserve(componentCount);
		m_components.clear();
		m_indices.reserve(componentCount);
		m_indices.clear();
		m_offsets.reserve(m_bsdfs.size());
		m_offsets.clear();

		int offset = 0;
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			const BSDF *bsdf = m_bsdfs[i];
			m_offsets.push_back(offset);

			for (int j=0; j<bsdf->getComponentCount(); ++j) {
				int componentType = bsdf->getType(j);
				m_components.push_back(componentType);
				m_indices.push_back(std::make_pair((int) i, j));
			}

			offset += bsdf->getComponentCount();
			m_usesRayDifferentials |= bsdf->usesRayDifferentials();
			m_pdf.append(m_weights[i].val());
		}
		
		m_pdf.normalize();
		BSDF::configure();
	}

	Spectrum_AD eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		Spectrum_AD result(0.0f);

		for (size_t i=0; i<m_bsdfs.size(); ++i)
			result += m_weights[i] * m_bsdfs[i]->eval(bRec, measure);
		return result;


	}


	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		Float result = 0.0f;
		for (size_t i=0; i<m_bsdfs.size(); ++i)
			result += m_bsdfs[i]->pdf(bRec, measure) * m_weights[i].val();
		return result;

	}


	FloatAD pdfAD(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		cout << "Unimplemented funmction" << endl;
		return FloatAD(1.0f);
	}


	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
		cout << "Unimplemented funmction" << endl;
		return Spectrum(1.0f);
	}


	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
		cout << "Unimplemented funmction" << endl;
		Point2 sample(_sample);
			/* Choose a component based on the normalized weights */
		size_t entry = m_pdf.sampleReuse(sample.x);
		Float pdf;
		Spectrum_AD result = m_bsdfs[entry]->sampleAD(bRec, pdf, sample);
		if (result.isZero()) // sampling failed
			return result;

		result *= m_weights[entry].m_data[Thread::getID()] * pdf;
		pdf *= m_pdf[entry];

		EMeasure measure = BSDF::getMeasure(bRec.sampledType);
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			if (entry == i)
				continue;
			pdf += m_bsdfs[i]->pdf(bRec, measure) * m_pdf[i];
			result += m_bsdfs[i]->eval(bRec, measure) * m_weights[i].m_data[Thread::getID()];
		}

		bRec.sampledComponent += m_offsets[entry];
		return result / pdf;
	}


	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
		cout << "Unimplemented funmction" << endl;
		return Spectrum(1.0f);
	}


	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
		
		Point2 newSample(_sample);
	
		size_t entry = m_pdf.sampleReuse(newSample.x);
		Spectrum_AD temp = m_bsdfs[entry]->sampleAD(bRec, pdf, newSample);

		if (temp.isZero()) // sampling failed
			return temp;

		Spectrum_AD result = Spectrum_AD(0.0f);
		FloatAD t = FloatAD(0.0f);
		EMeasure measure = BSDF::getMeasure(bRec.sampledType);
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			t += m_bsdfs[i]->pdfAD(bRec, measure) * m_pdf[i];
			result += m_weights[i] * m_bsdfs[i]->eval(bRec, measure);					

		}
		pdf = t;

		return result/t;
	}
	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			BSDF *bsdf = static_cast<BSDF *>(child);
			m_bsdfs.push_back(bsdf);
			bsdf->incRef();
		} else {
			BSDF::addChild(name, child);
		}
	}

	FloatAD getRoughness(const Intersection &its, int component) const {
		int bsdfIndex = m_indices[component].first;
		component = m_indices[component].second;
		return m_bsdfs[bsdfIndex]->getRoughness(its, component);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MixtureBSDF[" << endl
   			<< "  id = \"" << getID() << "\"," << endl
			<< "  weights = {";
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			oss << " " << m_weights[i].val();
			if (i + 1 < m_bsdfs.size())
				oss << ",";
		}
		oss << " }," << endl
			<< "  bsdfs = {" << endl;
		for (size_t i=0; i<m_bsdfs.size(); ++i)
			oss << "    " << indent(m_bsdfs[i]->toString(), 2) << "," << endl;
		oss << "  }" << endl
			<< "]";
		return oss.str();
	}


	MTS_DECLARE_CLASS()
private:
	Spectrum_shared testalpha;
	std::vector<FloatAD_shared> m_weights;
	std::vector<std::pair<int, int> > m_indices;
	std::vector<int> m_offsets;
	std::vector<BSDF *> m_bsdfs;
	DiscreteDistribution m_pdf;
};

MTS_IMPLEMENT_CLASS_S(MixtureBSDFAD, false, BSDF)
MTS_EXPORT_PLUGIN(MixtureBSDFAD, "Mixture BSDFAD")
MTS_NAMESPACE_END
