# Inverse Transport Networks

This repository is an implementation of the method described in the following paper: 

["Towards Learning-based Inverse Subsurface Scattering" [project website]](http://imaging.cs.cmu.edu/inverse_transport_networks/)\
[Chengqian Che](https://brucect2.github.io/), [Fujun Luan](https://www.cs.cornell.edu/~fujun/), [Shuang Zhao](https://shuangz.com/), [Kavita Bala](http://www.cs.cornell.edu/~kb/), and [Ioannis Gkioulekas](https://www.cs.cmu.edu/~igkioule/)\
IEEE International Conference on Computational Photography  (ICCP), 2020

## Getting Started

These instruction constains three parts:
* ITNSceneFiles - scripts and additional files needed for rendering our dataset images and derivatives
* Renderer - a Monte-Carlo Differentiable renderer based on [Mitsuba 0.5.0](https://www.mitsuba-renderer.org/releases/current/documentation.pdf)
* Learning - codes for training and evaluating our networks

### Rendering Scripts and Dataset

We used Mitsuba to generate our [dataset](http://http://imaging.cs.cmu.edu/inverse_transport_networks/). Our image file name follows the convention as:
[shape]_e[sunlight_direction]_d[sigmaT]_a[albedo]_g[g]_q[sampleCount].exr

For example, one can render the following scenes:
* shape: cube
* sunDirection: azimuth angle 30 and elevation angle 60
* materials: sigmaT 100, volumetric albedo 0.8 and g 0.2
* number of samples: 4096

```
mitsuba scenes/cube_sunsky.xml -Dmeshmodel=cube -DsigmaT=100 -Dalbedo=0.8 -Dg=0.2 -DnumSamples=4096 -Dx=0.433 -Dy=0.866 -Dz=0.25 -o cube_e30_d100_a0.8_g0.2_q4096.exr
```
One can also render a class of images using the following bash scripts with Sun Grid Engine:

```
./create_jobs_sge.sh
```
### Differentiable Renderer

We developed our differentiable renderer based on Mitsuba 0.5.0 and it compiles the same way as compiling Mistuba. To compile and render a signle image with derivatives:
```
cd renderer
scons
mitsuba scenesAD/cube_sunsky.xml -Dmeshmodel=cube -DsigmaT=100 -Dalbedo=0.8 -Dg=0.2 -DnumSamples=4096 -Dx=0.433 -Dy=0.866 -Dz=0.25 -o cube_e30_d100_a0.8_g0.2_q4096.exr
```
The current renderer supports computing derivatives with respect to the following parameters:
* scattering parameters: extinction coefficient, volumetric albedo and g in Henyey-Greenstein phase function
* bsdf parameters:roughness, reflectance and weights in a mixure of bsdfs
* lighitng: weights in a mixture of environment maps


In the sceneAD files, one can define 
The output image contains multiple channels with corrsponding channel names:
| channel name                      | Description                                                       |
| ----------------------------------| ------------------------------------------------------------      |
| forward                           | forward rendering                                                 |
| sigmaT                            | derivatives with respect to extinction coefficient                |
| albedo                            | derivatives with respect to volemetric albedo                                |
| g                                 | derivatives with respect to averaged cosine Henyey-Greenstein phase function |
| reflectance                       | derivatives with respect to surface albedo                                |
| alpha                             | derivatives with respect to surface roughness                                |
| weight                            | derivatives with respect to weights in a mixture models                                |

## Learning code

Codes used to train and evaluate our approach is inside folder learning/. Pre-trained models with 5 different networks can be downloaded here.

## Built With

* [Pytorch](https://pytorch.cn/previous-versions/)
* [Pandas](https://pandas.pydata.org/pandas-docs/version/0.15/tutorials.html)
* [scikit-image](https://scikit-image.org/docs/dev/api/skimage.html)
* [scikit-learn](https://scikit-learn.org/stable/)
* [pyEXR](https://pypi.org/project/PyEXR/)
* [scipy](https://www.scipy.org/)
* [EdgeBox](https://github.com/pdollar/edges) for creating weight maps

The networks were trained using Amazon EC2 clusters. All image names are in ITNSceneFiles/imgNames/. One can evaluate our model by doing:

```
python eval.py
``` 
And to use our models to initialize analysis by synthesis, one can run:
```
 python eval_calibrated.py
```
