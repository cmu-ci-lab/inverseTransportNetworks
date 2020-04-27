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
#include <mitsuba/hw/basicshader.h>

MTS_NAMESPACE_BEGIN

//testing here

class BlendBSDFAD : public BSDF {
public:
	BlendBSDFAD(const Properties &props)
		: BSDF(props) {
		int idx;
		for (int i = 0; i < props.id_AD_FloatAD.size(); i++) {
			if (props.id_AD_FloatAD[i] == "weight") {
				idx = i;
				break;
			}
		}
		m_weight = props.m_theta_FloatAD[idx];

	}

	BlendBSDFAD(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {

		m_weight = this->m_weight;
		m_bsdfs.push_back(static_cast<BSDF *>(manager->getInstance(stream)));
		m_bsdfs.push_back(static_cast<BSDF *>(manager->getInstance(stream)));
		m_bsdfs[0]->incRef();
		m_bsdfs[1]->incRef();
		configure();
	}

	virtual ~BlendBSDFAD() {
		for (size_t i=0; i<m_bsdfs.size(); ++i)
			m_bsdfs[i]->decRef();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		Assert(m_bsdfs.size() == 2);
		//manager->serialize(stream, m_weight);
		manager->serialize(stream, m_bsdfs[0]);
		manager->serialize(stream, m_bsdfs[1]);
	}

	void configure() {
		m_usesRayDifferentials = false;
		size_t componentCount = 0;

		if (m_bsdfs.size() != 2)
			Log(EError, "BSDF count mismatch: expected two nested BSDF instances!");

		for (size_t i=0; i<m_bsdfs.size(); ++i)
			componentCount += m_bsdfs[i]->getComponentCount();

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
		}
		BSDF::configure();
	}

	Spectrum_AD eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		std::cout << "start eval" << std::endl;
		FloatAD weight = m_weight.m_data[Thread::getID()];
		if (bRec.component == -1) {
			Spectrum_AD test = (1.0-weight) * m_bsdfs[0]->eval(bRec, measure)+
                                weight * m_bsdfs[1]->eval(bRec, measure);
			std::cout << "gg" << std::endl;			
			return
				(1.0-weight) * m_bsdfs[0]->eval(bRec, measure)+
				weight * m_bsdfs[1]->eval(bRec, measure);

		} else {
			/* Pick out an individual component */
			int idx = m_indices[bRec.component].first;
			if (idx == 0)
				weight = 1.0-weight;
			BSDFSamplingRecord bRec2(bRec);
			bRec2.component = m_indices[bRec.component].second;
			return m_bsdfs[idx]->eval(bRec2, measure) * weight;
		}




	}


	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		Float weight = m_weight.m_data[Thread::getID()].val();
		Float tt =  m_bsdfs[0]->pdf(bRec, measure) * (1-weight) +
                                m_bsdfs[1]->pdf(bRec, measure) * weight;
		std::cout << "start pdf" << std::endl;
		if (bRec.component == -1) {
			std::cout << "endpdf" << std::endl;
			return
				m_bsdfs[0]->pdf(bRec, measure) * (1-weight) +
				m_bsdfs[1]->pdf(bRec, measure) * weight;
		} else {
			/* Pick out an individual component */
			int idx = m_indices[bRec.component].first;
			if (idx == 0)
				weight = 1-weight;
			BSDFSamplingRecord bRec2(bRec);
			bRec2.component = m_indices[bRec.component].second;
			return m_bsdfs[idx]->pdf(bRec2, measure) * weight;
		}		

	}

	FloatAD pdfAD(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		cout << "Unimplemented function" << endl;
		return FloatAD(1.0f);
	}
	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
		cout << "Unimplemented function" << endl;
		Point2 sample(_sample);

		Float weights[2];
		weights[1] = m_weight.m_data[Thread::getID()].val();
		weights[0] = 1-weights[1];

		if (bRec.component == -1) {
			size_t entry;
			if (sample.x < weights[0]) {
				entry = 0; sample.x /= weights[0];
			} else {
				entry = 1; sample.x = (sample.x - weights[0]) / weights[1];
			}

			Float pdf;
			Spectrum result = m_bsdfs[entry]->sample(bRec, pdf, sample);
			if (result.isZero()) // sampling failed
				return result;

			result *= weights[entry] * pdf;
			pdf *= weights[entry];

			EMeasure measure = BSDF::getMeasure(bRec.sampledType);
			for (size_t i=0; i<m_bsdfs.size(); ++i) {
				if (entry == i)
					continue;
				pdf += m_bsdfs[i]->pdf(bRec, measure) * weights[i];
				Spectrum tmp = Spectrum(m_bsdfs[i]->eval(bRec, measure)[0]);
				result += tmp * weights[i];
			}

			bRec.sampledComponent += m_offsets[entry];
			return result/pdf;
		} else {
			/* Pick out an individual component */
			int requestedComponent = bRec.component;
			int bsdfIndex = m_indices[requestedComponent].first;
			bRec.component = m_indices[requestedComponent].second;
			Spectrum result = m_bsdfs[bsdfIndex]->sample(bRec, sample)
				* weights[bsdfIndex];
			bRec.component = bRec.sampledComponent = requestedComponent;
			return result;
		}
	}

	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, const Point2 &_sample) const {
		cout << "Unimplemented function" << endl;
		Point2 sample(_sample);

		FloatAD weights[2];
		weights[1] = m_weight.m_data[Thread::getID()];
		weights[0] = 1.0 -weights[1];


		size_t entry;
		if (sample.x < weights[0].val()) {
			entry = 0; sample.x /= weights[0].val();
		} else {
			entry = 1; sample.x = (sample.x - weights[0].val()) / weights[1].val();
		}

		Float pdf;
		Spectrum_AD result = m_bsdfs[entry]->sampleAD(bRec, pdf, sample);
		if (result.isZero()) // sampling failed
			return result;

		result *= weights[entry] * double(pdf);
		pdf *= weights[entry].val();

		EMeasure measure = BSDF::getMeasure(bRec.sampledType);
		for (size_t i=0; i<m_bsdfs.size(); ++i) {
			if (entry == i)
				continue;
			pdf += m_bsdfs[i]->pdf(bRec, measure) * weights[i].val();
			result += weights[i] * m_bsdfs[i]->eval(bRec, measure);
		}

		bRec.sampledComponent += m_offsets[entry];
		return result/pdf;
		
	}


	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
		cout << "Unimplemented function" << endl;
		Point2 sample(_sample);

		Float weights[2];
		weights[1] = m_weight.m_data[Thread::getID()].val();
		weights[0] = 1-weights[1];

		if (bRec.component == -1) {
			size_t entry;
			if (sample.x < weights[0]) {
				entry = 0; sample.x /= weights[0];
			} else {
				entry = 1; sample.x = (sample.x - weights[0]) / weights[1];
			}

			Spectrum result = m_bsdfs[entry]->sample(bRec, pdf, sample);
			if (result.isZero()) // sampling failed
				return result;

			result *= weights[entry] * pdf;
			pdf *= weights[entry];

			EMeasure measure = BSDF::getMeasure(bRec.sampledType);
			for (size_t i=0; i<m_bsdfs.size(); ++i) {
				if (entry == i)
					continue;
				pdf += m_bsdfs[i]->pdf(bRec, measure) * weights[i];
				Spectrum tmp = Spectrum(m_bsdfs[i]->eval(bRec, measure)[0]);
				result +=  tmp * weights[i];
			}

			bRec.sampledComponent += m_offsets[entry];
			return result/pdf;
		} else {
			/* Pick out an individual component */
			int requestedComponent = bRec.component;
			int bsdfIndex = m_indices[requestedComponent].first;
			bRec.component = m_indices[requestedComponent].second;
			Spectrum result = m_bsdfs[bsdfIndex]->sample(bRec, pdf, sample)
				* weights[bsdfIndex];
			bRec.component = bRec.sampledComponent = requestedComponent;
			return result;
		}
	}


	Spectrum_AD sampleAD(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const {
		std::cout << "start sampleAD" << std::endl;
		Point2 sample(_sample);

		FloatAD weights[2];
		weights[1] = m_weight.m_data[Thread::getID()];
		weights[0] = 1.0-weights[1];
		if (bRec.component == -1) {
			size_t entry;
			if (sample.x < weights[0].val()) {
				entry = 0; sample.x /= weights[0].val();
			} else {
				entry = 1; sample.x = (sample.x - weights[0].val()) / weights[1].val();
			}
			Spectrum_AD tmp = m_bsdfs[entry]->sampleAD(bRec, pdf, sample);
			if (tmp.isZero()) // sampling failed
				return tmp;
			EMeasure measure = BSDF::getMeasure(bRec.sampledType);

			Spectrum_AD result = Spectrum_AD(0.0f);
			Float t = 0.0f;

			for (size_t i=0; i<m_bsdfs.size(); ++i) {
				result += weights[i] * m_bsdfs[i]->eval(bRec, measure);				
				t += m_bsdfs[i]->pdf(bRec, measure) * weights[i].val();

			}



			bRec.sampledComponent += m_offsets[entry];
			pdf = t;
			std::cout << "end" << std::endl;
			return result/t;
		} else {
			/* Pick out an individual component */
			int requestedComponent = bRec.component;
			int bsdfIndex = m_indices[requestedComponent].first;
			bRec.component = m_indices[requestedComponent].second;
			Spectrum_AD result = m_bsdfs[bsdfIndex]->sampleAD(bRec, pdf, sample)
				* weights[bsdfIndex];
			bRec.component = bRec.sampledComponent = requestedComponent;
			return result;
		}
	}






	FloatAD getRoughness(const Intersection &its, int component) const {
		int bsdfIndex = m_indices[component].first;
		component = m_indices[component].second;
		return m_bsdfs[bsdfIndex]->getRoughness(its, component);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			BSDF *bsdf = static_cast<BSDF *>(child);
			m_bsdfs.push_back(bsdf);
			bsdf->incRef();
		} else if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "weight") {
			m_weight = static_cast<Texture *>(child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "blended bsdf" << endl;
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	std::vector<BSDF *> m_bsdfs;
	FloatAD_shared m_weight;
	std::vector<std::pair<int, int> > m_indices;
	std::vector<int> m_offsets;
};


MTS_IMPLEMENT_CLASS_S(BlendBSDFAD, false, BSDF)
MTS_EXPORT_PLUGIN(BlendBSDFAD, "BlendAD BSDF")
MTS_NAMESPACE_END
