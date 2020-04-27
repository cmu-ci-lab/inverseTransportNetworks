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

#include <mitsuba/core/statistics.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/renderproc.h>
#include <stan/math.hpp>

MTS_NAMESPACE_BEGIN

Integrator::Integrator(const Properties &props)
 : NetworkedObject(props) { }

Integrator::Integrator(Stream *stream, InstanceManager *manager)
 : NetworkedObject(stream, manager) { }

bool Integrator::preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) { return true; }
void Integrator::postprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) { }
void Integrator::serialize(Stream *stream, InstanceManager *manager) const {
	NetworkedObject::serialize(stream, manager);
}
void Integrator::configureSampler(const Scene *scene, Sampler *sampler) {
	/* Prepare the sampler for bucket-based rendering */
	sampler->setFilmResolution(scene->getFilm()->getCropSize(),
		getClass()->derivesFrom(MTS_CLASS(SamplingIntegrator)));
}
const Integrator *Integrator::getSubIntegrator(int idx) const { return NULL; }

SamplingIntegrator::SamplingIntegrator(const Properties &props)
 : Integrator(props) {

	m_theta_integrator = props.m_theta;
	m_theta_integrator_FloatAD = props.m_theta_FloatAD;
	m_id_AD_integrator_Spectrum_AD = props.id_AD_Spectrum_AD;
	m_id_AD_integrator_FloatAD = props.id_AD_FloatAD;
	m_theta_integrator_TextureAD = props.m_theta_TextureAD;
	m_id_AD_integrator_TextureAD = props.id_AD_TextureAD;

 }

SamplingIntegrator::SamplingIntegrator(Stream *stream, InstanceManager *manager)
 : Integrator(stream, manager) { 
	 cout << "calling stream" << endl;
 }

void SamplingIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
		 cout << "calling seralization" << endl;
	Integrator::serialize(stream, manager);
}

Spectrum SamplingIntegrator::E(const Scene *scene, const Intersection &its,
		const Medium *medium, Sampler *sampler, int nSamples, bool handleIndirect) const {
	Spectrum E(0.0f);
	RadianceQueryRecord query(scene, sampler);
	DirectSamplingRecord dRec(its);
	Frame frame(its.shFrame.n);

	sampler->generate(Point2i(0));
	for (int i=0; i<nSamples; i++) {
		/* Sample the direct illumination component */
		int maxIntermediateInteractions = -1;
		Spectrum directRadiance = Spectrum(scene->sampleAttenuatedEmitterDirect(
			dRec, its, medium, maxIntermediateInteractions, query.nextSample2D())[0].val());

		if (!directRadiance.isZero()) {
			Float dp = dot(dRec.d, its.shFrame.n);
			if (dp > 0)
				E += directRadiance * dp;
		}

		/* Sample the indirect illumination component */
		if (handleIndirect) {
			query.newQuery(RadianceQueryRecord::ERadianceNoEmission, medium);
			Vector d = frame.toWorld(warp::squareToCosineHemisphere(query.nextSample2D()));
			++query.depth;
			query.medium = medium;
			E += Spectrum(Li(RayDifferential(its.p, d, its.time), query)[0].val()) * M_PI;
		}

		sampler->advance();
	}

	return E / (Float) nSamples;
}

void SamplingIntegrator::cancel() {
	if (m_process)
		Scheduler::getInstance()->cancel(m_process);
}

bool SamplingIntegrator::render(Scene *scene,
		RenderQueue *queue, const RenderJob *job,
		int sceneResID, int sensorResID, int samplerResID) {
	ref<Scheduler> sched = Scheduler::getInstance();
	ref<Sensor> sensor = static_cast<Sensor *>(sched->getResource(sensorResID));
	ref<Film> film = sensor->getFilm();

	size_t nCores = sched->getCoreCount();
	const Sampler *sampler = static_cast<const Sampler *>(sched->getResource(samplerResID, 0));
	size_t sampleCount = sampler->getSampleCount();

	Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SIZE_T_FMT
		" %s, " SSE_STR ") ..", film->getCropSize().x, film->getCropSize().y,
		sampleCount, sampleCount == 1 ? "sample" : "samples", nCores,
		nCores == 1 ? "core" : "cores");

	/* This is a sampling-based integrator - parallelize */
	ref<ParallelProcess> proc = new BlockedRenderProcess(job,
		queue, scene->getBlockSize());
	int integratorResID = sched->registerResource(this);
	proc->bindResource("integrator", integratorResID);
	proc->bindResource("scene", sceneResID);
	proc->bindResource("sensor", sensorResID);
	proc->bindResource("sampler", samplerResID);
	scene->bindUsedResources(proc);
	bindUsedResources(proc);
	sched->schedule(proc);
	m_process = proc;
	sched->wait(proc);
	m_process = NULL;
	sched->unregisterResource(integratorResID);

	return proc->getReturnStatus() == ParallelProcess::ESuccess;
}

void SamplingIntegrator::bindUsedResources(ParallelProcess *) const {
	/* Do nothing by default */
}

void SamplingIntegrator::wakeup(ConfigurableObject *parent,
	std::map<std::string, SerializableObject *> &) {
	/* Do nothing by default */
}

void SamplingIntegrator::renderBlock(const Scene *scene,
		const Sensor *sensor, Sampler *sampler, ImageBlock *block,
		const bool &stop, const std::vector< TPoint2<uint8_t> > &points, int counter) const {

	Float diffScaleFactor = 1.0f /
		std::sqrt((Float) sampler->getSampleCount());
	
	int nDerivatives = m_theta_integrator.size() + m_theta_integrator_FloatAD.size() + m_theta_integrator_TextureAD.size();
	bool needsApertureSample = sensor->needsApertureSample();
	bool needsTimeSample = sensor->needsTimeSample();

	RadianceQueryRecord rRec(scene, sampler);
	Point2 apertureSample(0.5f);
	Float timeSample = 0.5f;
	RayDifferential sensorRay;
	block->clear();
	uint32_t queryType = RadianceQueryRecord::ESensorRay;
	
	if (!sensor->getFilm()->hasAlpha()) /* Don't compute an alpha channel if we don't have to */
		queryType &= ~RadianceQueryRecord::EOpacity;

	for (size_t i = 0; i<points.size(); ++i) {
		Point2i offset = Point2i(points[i]) + Vector2i(block->getOffset());
		if (stop)
			break;
		sampler->generate(offset);
		for (size_t j = 0; j<sampler->getSampleCount(); j++) {
			rRec.newQuery(queryType, sensor->getMedium());
			Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));

			if (needsApertureSample)
				apertureSample = rRec.nextSample2D();
			if (needsTimeSample)
				timeSample = rRec.nextSample1D();
			Spectrum spec = sensor->sampleRayDifferential(
				sensorRay, samplePos, apertureSample, timeSample);	
			sensorRay.scaleDifferential(diffScaleFactor);
			
			std::vector<Float> specVec = std::vector<Float>(nDerivatives + 1);

			Spectrum_AD specNew = Li(sensorRay, rRec);

			specVec[0] = spec[0];
			for (int num = 1; num < nDerivatives + 1; num++) {
				specVec[num] = 1.0f;
			}
			specVec[0] *= specNew[0].val();
			stan::math::set_zero_all_adjoints();
			specNew[0].grad();
			for (int i_theta = 0; i_theta < m_theta_integrator.size(); i_theta++) {
				specVec[i_theta + 1] *= m_theta_integrator[i_theta].adj();
			}
			for (int i_theta_FloatAD = 0; i_theta_FloatAD < m_theta_integrator_FloatAD.size(); i_theta_FloatAD++) {
				specVec[i_theta_FloatAD + 1 + m_theta_integrator.size()] *= m_theta_integrator_FloatAD[i_theta_FloatAD].adj();
			}



			block->put(samplePos, specVec, rRec.alpha, nDerivatives);
			StanStack::revert();


			sampler->advance();
		}

	}

}

MonteCarloIntegrator::MonteCarloIntegrator(const Properties &props) : SamplingIntegrator(props) {

	/* Depth to begin using russian roulette */
	m_rrDepth = props.getInteger("rrDepth", 5);

	/* Longest visualized path depth (\c -1 = infinite).
	   A value of \c 1 will visualize only directly visible light sources.
	   \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
	m_maxDepth = props.getInteger("maxDepth", -1);

	/**
	 * This parameter specifies the action to be taken when the geometric
	 * and shading normals of a surface don't agree on whether a ray is on
	 * the front or back-side of a surface.
	 *
	 * When \c strictNormals is set to \c false, the shading normal has
	 * precedence, and rendering proceeds normally at the risk of
	 * introducing small light leaks (this is the default).
	 *
	 * When \c strictNormals is set to \c true, the random walk is
	 * terminated when encountering such a situation. This may
	 * lead to silhouette darkening on badly tesselated meshes.
	 */
	m_strictNormals = props.getBoolean("strictNormals", false);

	/**value
	 * When this flag is set to true, contributions from directly
	 * visible emitters will not be included in the rendered image
	 */
	m_hideEmitters = props.getBoolean("hideEmitters", false);

	if (m_rrDepth <= 0)
		Log(EError, "'rrDepth' must be set to a value greater than zero!");

	if (m_maxDepth <= 0 && m_maxDepth != -1)
		Log(EError, "'maxDepth' must be set to -1 (infinite) or a value greater than zero!");
}

MonteCarloIntegrator::MonteCarloIntegrator(Stream *stream, InstanceManager *manager)
	: SamplingIntegrator(stream, manager) {
	m_rrDepth = stream->readInt();
	m_maxDepth = stream->readInt();
	m_strictNormals = stream->readBool();
	m_hideEmitters = stream->readBool();
}

void MonteCarloIntegrator::serialize(Stream *stream, InstanceManager *manager) const {
	SamplingIntegrator::serialize(stream, manager);
	stream->writeInt(m_rrDepth);
	stream->writeInt(m_maxDepth);
	stream->writeBool(m_strictNormals);
	stream->writeBool(m_hideEmitters);
}




std::string RadianceQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "RadianceQueryRecord[" << endl
		<< "  type = { ";
	if (type & EEmittedRadiance) oss << "emitted ";
	if (type & ESubsurfaceRadiance) oss << "subsurface ";
	if (type & EDirectSurfaceRadiance) oss << "direct ";
	if (type & EIndirectSurfaceRadiance) oss << "indirect ";
	if (type & ECausticRadiance) oss << "caustic ";
	if (type & EDirectMediumRadiance) oss << "inscatteredDirect ";
	if (type & EIndirectMediumRadiance) oss << "inscatteredIndirect ";
	if (type & EDistance) oss << "distance ";
	if (type & EOpacity) oss << "opacity ";
	if (type & EIntersection) oss << "intersectiLion ";
	oss << "}," << endl
		<< "  depth = " << depth << "," << endl
		<< "  its = " << indent(its.toString()) << endl
		<< "  alpha = " << alpha << "," << endl
		<< "  extra = " << extra << "," << endl
		<< "]" << endl;
	return oss.str();
}


MTS_IMPLEMENT_CLASS(Integrator, true, NetworkedObject)
MTS_IMPLEMENT_CLASS(SamplingIntegrator, true, Integrator)
MTS_IMPLEMENT_CLASS(MonteCarloIntegrator, true, SamplingIntegrator)
MTS_NAMESPACE_END
