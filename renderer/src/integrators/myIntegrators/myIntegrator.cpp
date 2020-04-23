#include  <mitsuba/render/scene.h>
MTS_NAMESPACE_BEGIN
class  MyIntegrator  :  public  SamplingIntegrator  {
public:
        MTS_DECLARE_CLASS()
        ///  Initialize  the  integrator  with  the  specified  properties
        MyIntegrator(const  Properties  &props)  :  SamplingIntegrator(props)  {
                Spectrum  defaultColor;
                defaultColor.fromLinearRGB(0.2f,  0.5f,  0.2f);
                m_color  =  props.getSpectrum("color",  defaultColor);
     }
             ///  Unserialize  from  a  binary  data  stream
        MyIntegrator(Stream  *stream,  InstanceManager  *manager)
                :  SamplingIntegrator(stream,  manager)  {
                m_color  =  Spectrum(stream);
        }
           ///  Serialize  to  a  binary  data  stream
        void  serialize(Stream  *stream,  InstanceManager  *manager)  const  {
                SamplingIntegrator::serialize(stream,  manager);
                m_color.serialize(stream);
        }
                ///  Query  for  an  unbiased  estimate  of  the  radiance  along  <tt>r</tt>
        Spectrum  Li(const  RayDifferential  &r,  RadianceQueryRecord  &rRec)  const  {
                return  m_color;
        }
private:
        Spectrum  m_color;
};
MTS_IMPLEMENT_CLASS_S(MyIntegrator,  false,  SamplingIntegrator)
MTS_EXPORT_PLUGIN(MyIntegrator,  "A  contrived  integrator");
MTS_NAMESPACE_END
