#include <mitsuba/render/scene.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/render/mipmap.h>
#include <mitsuba/hw/gpuprogram.h>
#include <mitsuba/hw/gputexture.h>

MTS_NAMESPACE_BEGIN
class mixenv : public Emitter {
public:
	mixenv(const Properties &props) : Emitter(props){
		for (int iIdx = 0; iIdx < props.id_AD_FloatAD.size(); iIdx++) {
			if (props.id_AD_FloatAD[iIdx].find("w") != std::string::npos) {
				m_weights.push_back(props.m_theta_FloatAD[iIdx]);
				m_weights_float.push_back(props.m_theta_FloatAD[iIdx].val());
			}
		}
		m_type |= EOnSurface | EEnvironmentEmitter;
	}

	mixenv(Stream *stream, InstanceManager *manager) : Emitter(stream, manager) {

		size_t envCount = stream->readSize();
		m_weights.resize(envCount);
		for (size_t i=0; i<envCount; ++i) {
			//m_weights[i] = stream->readFloat();
			Emitter *envmap = static_cast<Emitter *>(manager->getInstance(stream));
			envmap->incRef();
			m_envs.push_back(envmap);
		}


	}

	virtual ~mixenv() {
		for (size_t i=0; i<m_envs.size(); ++i)
			m_envs[i]->decRef();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Emitter::serialize(stream, manager);

		stream->writeSize(m_envs.size());
		for (size_t i=0; i<m_envs.size(); ++i) {
			//stream->writeFloat(m_weights[i]);
			manager->serialize(stream, m_envs[i]);
		}
	}

	void configure() {
		
		size_t componentCount = 0;

		if (m_envs.size() != m_weights.size())
			Log(EError, "Emitter count mismatch: " SIZE_T_FMT " bsdfs, but specified " SIZE_T_FMT " weights",
				m_envs.size(), m_weights.size());

		Float totalWeight = 0;
		for (size_t i=0; i<m_weights.size(); ++i)
			totalWeight += m_weights[i].val();

		if (totalWeight <= 0)
			Log(EError, "The weights must sum to a value greater than zero!");

		for (size_t i=0; i<m_envs.size(); ++i)
			componentCount += m_envs[i]->getComponentCount();
		m_pdf = DiscreteDistribution(m_envs.size());
		m_components.reserve(componentCount);
		m_components.clear();
		m_indices.reserve(componentCount);
		m_indices.clear();
		m_offsets.reserve(m_envs.size());
		m_offsets.clear();

		int offset = 0;
		for (size_t i=0; i<m_envs.size(); ++i) {
		 	const Emitter *envmap = m_envs[i];
		 	m_offsets.push_back(offset);

		 	for (int j=0; j<envmap->getComponentCount(); ++j) {
		 		m_indices.push_back(std::make_pair((int) i, j));
		 	}

		 	offset += envmap->getComponentCount();
		 	m_pdf.append(m_weights[i].val());
		}
		m_pdf.normalize();
		Emitter::configure();
	}

	bool fillDirectSamplingRecord(DirectSamplingRecord &dRec, const Ray &ray) const {
		bool result = true;
		for (size_t i=0; i<m_envs.size(); ++i) {
			result &= m_envs[i]->fillDirectSamplingRecord(dRec, ray);
		}

		return result;
	}


	std::vector<ref<Shape>> createShapes(const Scene *scene) {
		if (m_envs.size() != m_shapes.size()) {
			for (size_t i=0; i<m_envs.size(); ++i) {
				m_shapes.push_back(m_envs[i]->createShape(scene));
			}
		} else {
			for (size_t i=0; i<m_envs.size(); ++i) {
				m_shapes[i] = m_envs[i]->createShape(scene);
			}
		}

		return m_shapes;
	}

	bool isRegularEnv() const {
		return false;
	}
	Spectrum eval(const Intersection &its, const Vector &d) const {
		Spectrum result;
		const Intersection temp_its = its;  
		const Vector temp_d = d;  
		for (size_t i=0; i<m_envs.size(); ++i) {
			result +=  m_weights[i].val() * m_envs[i]->eval(temp_its, temp_d);
		}

		return result;
	}
	Spectrum_AD evalAD(const Intersection &its, const Vector &d) const {
		Spectrum_AD result;
		const Intersection temp_its = its;  
		const Vector temp_d = d;  
		for (size_t i=0; i<m_envs.size(); ++i) {
			result += m_weights[i] * m_envs[i]->evalAD(temp_its, temp_d);
		}

		return result;
	}

	Spectrum_AD evalEnvironment(const RayDifferential &ray) const {
		Spectrum_AD result = Spectrum_AD(0.0f);
		const RayDifferential temp_ray = ray;  
		for (size_t i=0; i<m_envs.size(); ++i) {
			result +=  m_weights[i] * m_envs[i]->evalEnvironment(temp_ray);
		}
		
		return result;
	}
	Spectrum evalEnvironment_sky(const RayDifferential &ray) const {
		cout << " Unimplemented funmction " << endl;
	}

	Spectrum samplePosition(PositionSamplingRecord &pRec, const Point2 &sample,
			const Point2 *extra) const {

		return Spectrum(1.0f);
	}



	Spectrum evalPosition(const PositionSamplingRecord &pRec) const {

		return Spectrum(1.0f);
	}

	Float pdfPosition(const PositionSamplingRecord &pRec) const {

		return Float(1.0f);
	}

	Spectrum sampleDirection(DirectionSamplingRecord &dRec,
			PositionSamplingRecord &pRec,
			const Point2 &sample,
			const Point2 *extra) const {

		return Spectrum(1.0f);

	}

	Float pdfDirection(const DirectionSamplingRecord &dRec,
			const PositionSamplingRecord &pRec) const {


		return Float(1.0f);
	}

	Spectrum evalDirection(const DirectionSamplingRecord &dRec,
			const PositionSamplingRecord &pRec) const {

		return Spectrum(1.0f);

	}

	Spectrum sampleRay(Ray &ray,
			const Point2 &spatialSample,
			const Point2 &directionalSample,
			Float time) const {

		return Spectrum(1.0f);
	}

	Spectrum_AD sampleDirectAD(DirectSamplingRecord &dRec, const Point2 &sample) const {

		//first randomly samples a direction using one of the envmaps
		Point2 newSample(sample);
		
		size_t entry = m_pdf.sampleReuse(newSample.x);
		Spectrum temp = m_envs[entry]->sampleDirect(dRec,newSample);
		

		Spectrum_AD result = Spectrum_AD(0.0f);
		Float t = 0.0f;
		const DirectSamplingRecord &newdRec = dRec;
		for (size_t i=0; i<m_envs.size(); ++i) {
			Float temppdf;
			Spectrum tem_result = m_envs[i]->evalDirect(newdRec, temppdf);

			result += m_weights[i] * tem_result * temppdf;
			t += temppdf * m_weights[i].val();
		}
		dRec.pdf = t;
		return result/t;
			

	}

	Float pdfDirect(const DirectSamplingRecord &dRec) const {
				
		Float result = Float(0.0f);
		for (size_t i=0; i<m_envs.size(); ++i) {

			result += m_weights[i].val() * m_envs[i]->pdfDirect(dRec);
		}
		return result;
	}

	std::vector<AABB> getAABBs(){
		/* The scene sets its bounding box so that it contains all shapes and
		   emitters, but this particular emitter always wants to be *a little*
		   bigger than the scene. To avoid a silly recursion, just return a
		   point here. */

		if (m_envs.size() != m_AABBs.size()) {
			for (size_t i=0; i<m_envs.size(); ++i) {
				m_AABBs.push_back(m_envs[i]->getAABB());
			}
		} else {
			for (size_t i=0; i<m_envs.size(); ++i) {
				m_AABBs[i] = m_envs[i]->getAABB();
			}
		}

		return m_AABBs;
	}

	AABB getAABB() const {
		/* The scene sets its bounding box so that it contains all shapes and
		   emitters, but this particular emitter always wants to be *a little*
		   bigger than the scene. To avoid a silly recursion, just return a
		   point here. */
		return AABB();
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		
		if (child->getClass()->derivesFrom(MTS_CLASS(Emitter))) {
			Emitter *envmap = static_cast<Emitter *>(child);
			m_envs.push_back(envmap);
			envmap->incRef();
		} else {
			Emitter::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;

		return oss.str();
	}

	

	MTS_DECLARE_CLASS()
private:
	
	std::vector<FloatAD_shared> m_weights;
	std::vector<Float> m_weights_float;
	std::vector<Emitter *> m_envs;
	std::vector<ref<Shape>> m_shapes;
	std::vector<AABB> m_AABBs;
	std::vector<std::pair<int, int> > m_indices;
	std::vector<int> m_offsets;
	DiscreteDistribution m_pdf;
	
};

MTS_IMPLEMENT_CLASS_S(mixenv, false, Emitter)
MTS_EXPORT_PLUGIN(mixenv, "mixture Environment map");
MTS_NAMESPACE_END
