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

#include <mitsuba/core/bitmapAD.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/sched.h>
#include <mitsuba/render/textureAD.h>
#include <mitsuba/render/mipmapAD.h>
#include <mitsuba/hw/renderer.h>
#include <mitsuba/hw/gputexture.h>
#include <mitsuba/hw/gpuprogram.h>
#include <boost/algorithm/string.hpp>
#include <mitsuba/core/barray.h>
MTS_NAMESPACE_BEGIN


class BitmapTextureAD : public Texture2DAD {
public:
	/* Store texture data using half precision, but perform computations in
	   single/double precision based on compilation flags. The following
	   generates efficient implementations for both luminance and RGB data */
	typedef TSpectrum_AD<FloatAD, 1> Color1;
	typedef TSpectrum_AD<FloatAD, 3> Color3;
	typedef TSpectrum_AD<half, 1>  Color1h;
	typedef TSpectrum_AD<half, 3>  Color3h;
	typedef TMIPMapAD<Color1, Color1h> MIPMap1;
	typedef TMIPMapAD<Color3, Color3h> MIPMap3;
	BitmapTextureAD(const Properties &props) : Texture2DAD(props) {
		uint64_t timestamp = 0;
		bool tryReuseCache = false;
		fs::path cacheFile;
		ref<Bitmap> bitmap;

		m_channel = boost::to_lower_copy(props.getString("channel", ""));

		if (props.hasProperty("bitmapAD")) {
			/* Support initialization via raw data passed from another plugin */
			bitmap = reinterpret_cast<Bitmap *>(props.getData("bitmapAD").ptr);
		} else {
			m_filename = Thread::getThread()->getFileResolver()->resolve(
				props.getString("filename"));

			Log(EInfo, "Loading texture \"%s\"", m_filename.filename().string().c_str());
			if (!fs::exists(m_filename))
				Log(EError, "Texture file \"%s\" could not be found!", m_filename.string().c_str());

			boost::system::error_code ec;
			timestamp = (uint64_t) fs::last_write_time(m_filename, ec);
			if (ec.value())
				Log(EError, "Could not determine modification time of \"%s\"!", m_filename.string().c_str());

			cacheFile = m_filename;

			if (m_channel.empty())
				cacheFile.replace_extension(".mip");
			else
				cacheFile.replace_extension(formatString(".%s.mip", m_channel.c_str()));

			tryReuseCache = fs::exists(cacheFile) && props.getBoolean("cache", true);
		}

		std::string filterType = boost::to_lower_copy(props.getString("filterType", "ewa"));
		std::string wrapMode = props.getString("wrapMode", "repeat");
		m_wrapModeU = parseWrapMode(props.getString("wrapModeU", wrapMode));
		m_wrapModeV = parseWrapMode(props.getString("wrapModeV", wrapMode));

		m_gamma = props.getFloat("gamma", 0);

		if (filterType == "ewa")
			m_filterType = EEWA;
		else if (filterType == "bilinear")
			m_filterType = EBilinear;
		else if (filterType == "trilinear")
			m_filterType = ETrilinear;
		else if (filterType == "nearest")
			m_filterType = ENearest;
		else
			Log(EError, "Unknown filter type '%s' -- must be "
				"'ewa', 'trilinear', or 'nearest'!", filterType.c_str());

		m_maxAnisotropy = props.getFloat("maxAnisotropy", 20);

		if (m_filterType != EEWA)
			m_maxAnisotropy = 1.0f;

		if (tryReuseCache && MIPMap3::validateCacheFile(cacheFile, timestamp,
				Bitmap::ERGB, m_wrapModeU, m_wrapModeV, m_filterType, m_gamma)) {
			/* Reuse an existing MIP map cache file */
			m_mipmap3 = new MIPMap3(cacheFile, m_maxAnisotropy);
		} else if (tryReuseCache && MIPMap1::validateCacheFile(cacheFile, timestamp,
				Bitmap::ELuminance, m_wrapModeU, m_wrapModeV, m_filterType, m_gamma)) {
			/* Reuse an existing MIP map cache file */
			m_mipmap1 = new MIPMap1(cacheFile, m_maxAnisotropy);
		} else {
			if (bitmap == NULL) {
				/* Load the input image if necessary */
				ref<Timer> timer = new Timer();
				ref<FileStream> fs = new FileStream(m_filename, FileStream::EReadOnly);
				bitmap = new Bitmap(Bitmap::EAuto, fs);
				if (m_gamma != 0)
					bitmap->setGamma(m_gamma);
				Log(EDebug, "Loaded \"%s\" in %i ms", m_filename.filename().string().c_str(),
					timer->getMilliseconds());
			}

			Bitmap::EPixelFormat pixelFormat;
			if (!m_channel.empty()) {
				/* Create a texture from a certain channel of an image */
				pixelFormat = Bitmap::ELuminance;
				bitmap = bitmap->extractChannel(findChannel(bitmap, m_channel));
				if (m_channel == "a")
					bitmap->setGamma(1.0f);
			} else {
				switch (bitmap->getPixelFormat()) {
					case Bitmap::ELuminance:
					case Bitmap::ELuminanceAlpha:
						pixelFormat = Bitmap::ELuminance;
						break;
					case Bitmap::ERGB:
					case Bitmap::ERGBA:
						pixelFormat = Bitmap::ERGB;
						break;
					default:
						Log(EError, "The input image has an unsupported pixel format!");
						return;
				}
			}

			/* (Re)generate the MIP map hierarchy; downsample using a
			    2-lobed Lanczos reconstruction filter */
			Properties rfilterProps("lanczos");
			rfilterProps.setInteger("lobes", 2);
			ref<ReconstructionFilter> rfilter = static_cast<ReconstructionFilter *> (
				PluginManager::getInstance()->createObject(
				MTS_CLASS(ReconstructionFilter), rfilterProps));
			rfilter->configure();
			
			/* Potentially create a new MIP map cache file */
			bool createCache = !cacheFile.empty() && props.getBoolean("cache",
				bitmap->getSize().x * bitmap->getSize().y > 1024*1024);
			if (pixelFormat == Bitmap::ELuminance)
				m_mipmap1 = new MIPMap1(bitmap, pixelFormat, Bitmap::EFloatAD,
					rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
					createCache ? cacheFile : fs::path(), timestamp);
			else
				m_mipmap3 = new MIPMap3(bitmap, pixelFormat, Bitmap::EFloatAD,
					rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
					createCache ? cacheFile : fs::path(), timestamp);
		}
	}

	static int findChannel(const Bitmap *bitmap, const std::string channel) {
		int found = -1;
		std::string channelNames;
		for (int i=0; i<bitmap->getChannelCount(); ++i) {
			std::string name = boost::to_lower_copy(bitmap->getChannelName(i));
			if (name == channel)
				found = i;
			channelNames += name;
			if (i + 1 < bitmap->getChannelCount())
				channelNames += std::string(", ");
		}

		if (found == -1) {
			Log(EError, "Channel \"%s\" not found! Must be one of: [%s]",
				channel.c_str(), channelNames.c_str());
		}

		return found;
	}

	inline ReconstructionFilter::EBoundaryCondition parseWrapMode(const std::string &wrapMode) {
		if (wrapMode == "repeat")
			return ReconstructionFilter::ERepeat;
		else if (wrapMode == "clamp")
			return ReconstructionFilter::EClamp;
		else if (wrapMode == "mirror")
			return ReconstructionFilter::EMirror;
		else if (wrapMode == "zero" || wrapMode == "black")
			return ReconstructionFilter::EZero;
		else if (wrapMode == "one" || wrapMode == "white")
			return ReconstructionFilter::EOne;
		else
			Log(EError, "Unknown wrap mode '%s' -- must be "
				"'repeat', 'clamp', 'black', or 'white'!", wrapMode.c_str());
		return ReconstructionFilter::EZero; // make gcc happy
	}

	BitmapTextureAD(Stream *stream, InstanceManager *manager)
	 : Texture2DAD(stream, manager) {
		m_filename = stream->readString();
		Log(EDebug, "Unserializing texture \"%s\"", m_filename.filename().string().c_str());
		m_filterType = (EMIPFilterType) stream->readUInt();
		m_wrapModeU = (ReconstructionFilter::EBoundaryCondition) stream->readUInt();
		m_wrapModeV = (ReconstructionFilter::EBoundaryCondition) stream->readUInt();
		m_gamma = stream->readFloat();
		m_maxAnisotropy = stream->readFloat();
		m_channel = stream->readString();

		size_t size = stream->readSize();
		ref<MemoryStream> mStream = new MemoryStream(size);
		stream->copyTo(mStream, size);
		mStream->seek(0);
		ref<Bitmap> bitmap = new Bitmap(Bitmap::EAuto, mStream);
		if (m_gamma != 0)
			bitmap->setGamma(m_gamma);

		/* Downsample using a 2-lobed Lanczos reconstruction filter */
		Properties rfilterProps("lanczos");
		rfilterProps.setInteger("lobes", 2);
		ref<ReconstructionFilter> rfilter = static_cast<ReconstructionFilter *> (
			PluginManager::getInstance()->createObject(
			MTS_CLASS(ReconstructionFilter), rfilterProps));
		rfilter->configure();

		Bitmap::EPixelFormat pixelFormat;
		if (!m_channel.empty()) {
			/* Create a texture from a certain channel of an image */
			pixelFormat = Bitmap::ELuminance;
			bitmap = bitmap->extractChannel(findChannel(bitmap, m_channel));
			if (m_channel == "a")
				bitmap->setGamma(1.0f);
		} else {
			switch (bitmap->getPixelFormat()) {
				case Bitmap::ELuminance:
				case Bitmap::ELuminanceAlpha:
					pixelFormat = Bitmap::ELuminance;
					break;
				case Bitmap::ERGB:
				case Bitmap::ERGBA:
					pixelFormat = Bitmap::ERGB;
					break;
				default:
					Log(EError, "The input image has an unsupported pixel format!");
					return;
			}
		}
		if (pixelFormat == Bitmap::ELuminance)
			m_mipmap1 = new MIPMap1(bitmap, pixelFormat, Bitmap::EFloatAD,
				rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
				fs::path(), 0);
		else
			m_mipmap3 = new MIPMap3(bitmap, pixelFormat, Bitmap::EFloatAD,
				rfilter, m_wrapModeU, m_wrapModeV, m_filterType, m_maxAnisotropy,
				fs::path(), 0);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Texture2DAD::serialize(stream, manager);
		stream->writeString(m_filename.string());
		stream->writeUInt(m_filterType);
		stream->writeUInt(m_wrapModeU);
		stream->writeUInt(m_wrapModeV);
		stream->writeFloat(m_gamma);
		stream->writeFloat(m_maxAnisotropy);

		if (!m_filename.empty() && fs::exists(m_filename)) {
			/* We still have access to the original image -- use that, since
			   it is probably much smaller than the in-memory representation */
			ref<Stream> is = new FileStream(m_filename, FileStream::EReadOnly);
			stream->writeString(m_channel);
			stream->writeSize(is->getSize());
			is->copyTo(stream);
		} else {
			/* No access to the original image anymore. Create an EXR image
			   from the top MIP map level and serialize that */
			ref<MemoryStream> mStream = new MemoryStream();
			ref<BitmapAD> bitmapAD = m_mipmap1.get() ?
				m_mipmap1->toBitmapAD() : m_mipmap3->toBitmapAD();
			bitmapAD->write(BitmapAD::EOpenEXR, mStream);

			stream->writeString("");
			stream->writeSize(mStream->getSize());
			stream->write(mStream->getData(), mStream->getSize());
		}
	}

	Spectrum_AD eval(const Point2 &uv) const {
		/* There are no ray differentials to do any kind of
		   prefiltering. Evaluate the full-resolution texture */

		Spectrum_AD result;
		
		
		return result;
	}

	void evalGradient3(const Point2 &uv, FloatAD &val00,  FloatAD &val01, FloatAD &val02, FloatAD &val10, FloatAD &val11, FloatAD &val12) const {
		/* There are no ray differentials to do any kind of
		   prefiltering. Evaluate the full-resolution texture */

	
		Color3 result[2];
		if (m_mipmap3->getFilterType() != ENearest) {
			m_mipmap3->evalGradientBilinear(0, uv, result);
			val00 = result[0][0];
			val01 = result[0][1];
			val02 = result[0][2];

			val10 = result[1][0];
			val11 = result[1][1];
			val12 = result[1][2];		
		} else {
			val00 = FloatAD(0.0f);
			val01 = FloatAD(0.0f);
			val02 = FloatAD(0.0f);

			val10 = FloatAD(0.0f);
			val11 = FloatAD(0.0f);
			val12 = FloatAD(0.0f);
		}
		
		stats::filteredLookupsAD.incrementBase();
	}

	ref<BitmapAD> getBitmap(const Vector2i &/* unused */) const {
		return m_mipmap1.get() ? m_mipmap1->toBitmapAD() : m_mipmap3->toBitmapAD();
	}

	Spectrum_AD eval(const Point2 &uv, const Vector2 &d0, const Vector2 &d1) const {
		stats::filteredLookupsAD.incrementBase();
		++stats::filteredLookupsAD;

		Spectrum_AD result(0.0f);
		// if (m_mipmap3.get()) {
		// 	Color3 value = m_mipmap3->eval(uv, d0, d1);
		// 	result.fromLinearRGB(value[0], value[1], value[2]);
		// } else {
		// 	Color1 value = m_mipmap1->eval(uv, d0, d1);
		// 	result = Spectrum(value[0]);
		// }
		return result;
	}

	void eval3(const Point2 &uv, const Vector2 &d0, const Vector2 &d1, FloatAD &val1, FloatAD &val2, FloatAD &val3) const {
		stats::filteredLookupsAD.incrementBase();
		++stats::filteredLookupsAD;


		Color3 value = m_mipmap3->eval(uv, d0, d1);
		val1 = value[0];
		val2 = value[1];
		val3 = value[2];


	}


	void getPixel3(const Point2 &uv, FloatAD &val1, FloatAD &val2) const {
		Color3 value;

		value = m_mipmap3->evalBox3(0, uv);

		val1 = value[0];
		val2 = value[1];	
		stats::filteredLookupsAD.incrementBase();
	}

	void eval3(const Point2 &uv, FloatAD &val1, FloatAD &val2, FloatAD &val3) const {
		/* There are no ray differentials to do any kind of
		   prefiltering. Evaluate the full-resolution texture */
		Color3 value;
		if (m_mipmap3->getFilterType() != ENearest) {
			value = m_mipmap3->evalBilinear(0, uv);
			
		} else {
			value = m_mipmap3->evalBox(0, uv);
		}
		val1 = value[0];
		val2 = value[1];
		val3 = value[2];		
		stats::filteredLookupsAD.incrementBase();
		
	}	


	Spectrum_AD getAverage() const {
		Spectrum_AD result(0.0f);
		// if (m_mipmap3.get()) {
		// 	Color3 value = m_mipmap3->getAverage();
		// 	result.fromLinearRGB(value[0], value[1], value[2]);
		// } else {
		// 	Color1 value = m_mipmap1->getAverage();
		// 	result = Spectrum(value[0]);
		// }
		return result;
	}

	Spectrum_AD getMaximum() const {
		Spectrum_AD result(0.0f);
		// if (m_mipmap3.get()) {
		// 	Color3 value = m_mipmap3->getMaximum();
		// 	result.fromLinearRGB(value[0], value[1], value[2]);
		// } else {
		// 	Color1 value = m_mipmap1->getMaximum();
		// 	result = Spectrum(value[0]);
		// }
		return result;
	}

	Spectrum_AD getMinimum() const {
		Spectrum_AD result(0.0f);
		// if (m_mipmap3.get()) {
		// 	Color3 value = m_mipmap3->getMinimum();
		// 	result.fromLinearRGB(value[0], value[1], value[2]);
		// } else {
		// 	Color1 value = m_mipmap1->getMinimum();
		// 	result = Spectrum(value[0]);
		// }
		return result;
	}

	bool isConstant() const {
		return false;
	}

	bool usesRayDifferentials() const {
		return true;
	}

	bool isMonochromatic() const {
		return m_mipmap1.get() != NULL;
	}

	Vector3i getResolution() const {
		if (m_mipmap3.get()) {
			return Vector3i(
				m_mipmap3->getWidth(),
				m_mipmap3->getHeight(),
				1
			);
		} else {
			return Vector3i(
				m_mipmap1->getWidth(),
				m_mipmap1->getHeight(),
				1
			);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "BitmapTexture[" << endl
			<< "  filename = \"" << m_filename.string() << "\"," << endl;

		if (m_mipmap3.get())
			oss << "  mipmap = " << indent(m_mipmap3.toString()) << endl;
		else
			oss << "  mipmap = " << indent(m_mipmap1.toString()) << endl;

		oss << "]";
		return oss.str();
	}


	MTS_DECLARE_CLASS()
protected:
	ref<MIPMap1> m_mipmap1;
	ref<MIPMap3> m_mipmap3;
	EMIPFilterType m_filterType;
	ReconstructionFilter::EBoundaryCondition m_wrapModeU;
	ReconstructionFilter::EBoundaryCondition m_wrapModeV;
	Float m_gamma, m_maxAnisotropy;
	std::string m_channel;
	fs::path m_filename;
};


MTS_IMPLEMENT_CLASS_S(BitmapTextureAD, false, Texture2DAD)
MTS_EXPORT_PLUGIN(BitmapTextureAD, "Bitmap textureAD");
MTS_NAMESPACE_END
