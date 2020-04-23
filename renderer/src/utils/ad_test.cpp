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

#include <iostream>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/stan.h>
#include <mitsuba/render/util.h>
#include <stan/math/rev/core.hpp>
#include <mitsuba/render/textureAD.h>
#include <mitsuba/core/bitmapAD.h>
MTS_NAMESPACE_BEGIN

class AutoDiffTesting : public Utility {
public:
	int run(int argc, char **argv) {
	// 	// Properties props("bitmap");
	// 	// props.setString("filename", "/usr0/home/cche/CMU/Research/PhD/Project/mitsuba_autodiff/dist/normalmaptest/mitsuba.exr");
	// 	// //PluginManager::getInstance()->createObject(MTS_CLASS(TextureAD), props);
	// 	// Texture2D *m_texture = static_cast<Texture2D *> (PluginManager::getInstance()->createObject(MTS_CLASS(Texture2D), props));



	// 	Properties props("bitmapAD");
	// 	//props.setString("filename", "myplugin.exr");
	// 	props.setString("filterType", "nearest");		
	// 	props.setString("filename", "/usr0/home/cche/CMU/Research/PhD/Project/mitsuba_autodiff/dist/normalmaptest/myplugin.exr");
	// 	Texture2DAD *m_texture = static_cast<Texture2DAD *> (PluginManager::getInstance()->createObject(MTS_CLASS(Texture2DAD), props));
		
		
	// 	int x = m_texture->getResolution().x;
	// 	int y = m_texture->getResolution().y;
	// 	int i = 0;
	// 	int j = 0;
	// 	Point2 uv(i, j);
		
		
	// 	FloatAD theta;
	// 	FloatAD phi;
	// 	FloatAD r;
	// 	FloatAD theta1;
	// 	FloatAD phi1;
	// 	FloatAD r1;
	// 	m_texture->eval3(uv,theta, phi, r);
	// 	m_texture->eval3(uv,theta1, phi1, r1);

	// 	stan::math::set_zero_all_adjoints();
	// 	theta1.grad();
	// 	cout << theta.val() << endl;
	// 	cout << theta.adj() << endl;
	// 	// az = (az/255)*2*3.1415926f;
	// 	// ele = (ele/255)*3.1415926f;

        // // FloatAD newVal = r * stan::math::sin(ele) * stan::math::cos(az);
        // // // n.y = stan::math::sin(theta) * stan::math::sin(phi);
	// 	// // n.z = stan::math::cos(theta);



	// 	// m_texture->eval(uv);
	// 	// cout <<xx << endl;
	// 	return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(AutoDiffTesting, "Auto diff. tester")
MTS_NAMESPACE_END
