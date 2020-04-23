from __future__ import print_function, division

import torch
import torchvision
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.autograd import Variable
from torch.utils.data import Dataset, DataLoader
from torchvision import transforms, utils

import pandas as pd
from skimage import io, transform
import numpy as np 
import scipy
import scipy.misc
import time
import os
import math

import pyexr
import sys 
from random import shuffle
import scipy.io
from sklearn.feature_extraction.image import extract_patches


class MatNet(nn.Module):
	def __init__(self):
		super(MatNet, self).__init__()
		# Material network based on http://openaccess.thecvf.com/content_ICCV_2017/supplemental/Liu_Material_Editing_Using_ICCV_2017_supplemental.pdf
		self.conv1 = nn.Conv2d(2   , 16  , 3, padding=1)
		self.conv2 = nn.Conv2d(16  , 32  , 3, padding=1)
		self.conv3 = nn.Conv2d(32  , 64  , 3, padding=1) 
		self.conv4 = nn.Conv2d(64  , 128 , 3, padding=1)
		self.conv5 = nn.Conv2d(128 , 256 , 3, padding=1)
		self.conv6 = nn.Conv2d(256 , 512 , 3, padding=1)  
		self.conv7 = nn.Conv2d(512 , 1024, 3, padding=1)
		self.fc8   = nn.Linear(1024, 3)

	
	def forward(self, x):
    		# x -> 3 x 256 x 256
		x, self.pool_idx1 = F.max_pool2d(F.relu(self.conv1(x)), kernel_size=2, return_indices=True) 
		# x -> 16 x 128 x 128
		x, self.pool_idx2 = F.max_pool2d(F.relu(self.conv2(x)), kernel_size=2, return_indices=True)
		# x -> 32 x 64 x 64
		x, self.pool_idx3 = F.max_pool2d(F.relu(self.conv3(x)), kernel_size=2, return_indices=True)
		# x -> 64 x 32 x 32
		x, self.pool_idx4 = F.max_pool2d(F.relu(self.conv4(x)), kernel_size=2, return_indices=True)
		# # x -> 128 x 16 x 16
		x, self.pool_idx5 = F.max_pool2d(F.relu(self.conv5(x)), kernel_size=2, return_indices=True)
		# # x -> 256 x 8 x 8
		x, self.pool_idx6 = F.max_pool2d(F.relu(self.conv6(x)), kernel_size=2, return_indices=True)
		# # x -> 512 x 4 x 4
		x, self.pool_idx7 = F.max_pool2d(F.relu(self.conv7(x)), kernel_size=2, return_indices=True) 
		# x -> 1024 x 2 x 2
		x = x.view(-1, self.num_flat_features(x))
		x = self.fc8(x)

		return x


	def num_flat_features(self, x):
		size = x.size()[1:]  # all dimensions except the batch dimension
		num_features = 1
		for s in size:
			num_features *= s
		return num_features

class AllShapeDataset(Dataset):
	def __init__(self, csv_file, root_dir, mask_dir, env_dir, wm_dir, transform=None):
				self.data_frame = pd.read_csv(csv_file)
				self.root_dir   = root_dir
				self.mask_dir   = mask_dir
				self.env_dir   = env_dir				
				self.transform  = transform
				self.wm_dir     = wm_dir
	def __len__(self):
		return len(self.data_frame)
	def __getitem__(self, idx):
  		
			input     = np.zeros(shape = (128, 128, 2))

			img_fn     = os.path.join(self.root_dir, self.data_frame.ix[idx, 0])
			luminance  = pyexr.open(img_fn).get() # luminance 
			wm_fn      = os.path.join(self.wm_dir, self.data_frame.ix[idx, 0])
			wm         = pyexr.open(wm_fn).get() # luminance 
			mask_fn    =  os.path.join(self.mask_dir, self.data_frame.ix[idx, 1] + '.exr')
			mask       = pyexr.open(mask_fn).get()
			luminance  = 28.2486 * luminance * mask

			shape 	  = self.data_frame.ix[idx, 1]
			light 	  = self.data_frame.ix[idx, 2]

			sigmaT 	  = self.data_frame.ix[idx, 3].astype('float')
			albedo    = self.data_frame.ix[idx, 4].astype('float')
			hg        = self.data_frame.ix[idx, 5].astype('float') 


			patches_luminance = extract_patches(luminance[:,:,0], patch_shape=(4, 4), extraction_step=(4, 4))
			input[:,:,0] = patches_luminance.mean(-1).mean(-1)
			patches_wm = extract_patches(wm[:,:,0], patch_shape=(4, 4), extraction_step=(4, 4))
                	
			input[:,:,1] = patches_wm.mean(-1).mean(-1)
			weights = input[:,:,1]
			target = torch.from_numpy(patches_luminance.mean(-1).mean(-1)).unsqueeze(0)

			input = input.transpose((2, 0, 1))
			input = torch.from_numpy(input)

			params    = torch.Tensor(3)
			params[0] = sigmaT
			params[1] = albedo 
			params[2] = hg 

			sample = {
					'input' : input,
					'params': params,
					'target' : target,
					'shape' : shape,
					'light' : light,
					'sigmaT' : sigmaT,
					'albedo' : albedo,
					'hg' : hg,
					'fn' : self.data_frame.ix[idx, 0],
					'weight':weights
			}

			return sample

test = torch.load('ITNNetworks/ITN_similar_map.bin', map_location=lambda storage, loc: storage)
net = test['model']



class MitsubaLossFunction(torch.autograd.Function):
	def __init__(self):
		super(MitsubaLossFunction, self).__init__()

	def forward(self, input, tmp_img, tmp_grad_sigmaT, tmp_grad_albedo, tmp_grad_hg, gt_img, wm_img):

		diff = 28.2486*tmp_img - gt_img

		self.save_for_backward(input, tmp_img.float(), tmp_grad_sigmaT.float(), tmp_grad_albedo.float(), tmp_grad_hg.float(), gt_img.float(), wm_img.float())
		loss = torch.Tensor(1, 1)
		loss_tmp = torch.pow(diff,2)*wm_img
		loss[0][0] = loss_tmp.mean()
		

		return loss

	def backward(self, grad_out):
		input, tmp_img,tmp_grad_sigmaT, tmp_grad_albedo, tmp_grad_hg, gt_img, wm_img= self.saved_tensors
		diff = 28.2486*tmp_img - gt_img
		grad_in = torch.Tensor(input.size()).float()

		#print( ((diff * tmp_grad_sigmaT).mean(2).mean(1) * grad_out.clone())[0][0])	
		#print(grad_in[0])
		grad_in[0] = ((diff * tmp_grad_sigmaT*wm_img).mean(2).mean(1) * grad_out.clone())[0][0]
		grad_in[1] = ((diff * tmp_grad_albedo*wm_img).mean(2).mean(1) * grad_out.clone())[0][0]
		grad_in[2] = ((diff * tmp_grad_hg*wm_img).mean(2).mean(1) * grad_out.clone())[0][0]


		return grad_in, None, None, None, None, None, None


class MitsubaLossModule(torch.nn.Module):
	def __init__(self):
		super(MitsubaLossModule, self).__init__()

	def forward(self, input, tmp_img,tmp_grad_sigmaT, tmp_grad_albedo, tmp_grad_hg,  gt_img, wm_img):
		return MitsubaLossFunction()(input, tmp_img, tmp_grad_sigmaT, tmp_grad_albedo, tmp_grad_hg, gt_img, wm_img)




def loadFigure(filename):
	try:
		tmp = pyexr.open(filename)
		return tmp

	except Exception as e:
		return None




batch_size = 1



five_shape_ten_light_train_dataset = AllShapeDataset(
						csv_file='calibrated.csv',
						root_dir='testSet/',
						mask_dir='masks/',
						env_dir='xx/',
						wm_dir='maps/'

					)


data_train_loader = DataLoader(five_shape_ten_light_train_dataset, 
							batch_size=batch_size, 
							shuffle=False, 
							num_workers=16)



criterion_Para = nn.MSELoss(size_average=False)




totalST = 0

totalABD = 0

totalG = 0


for i_batch, sample in enumerate(data_train_loader):
	print('batch strats now----------',i_batch)    
	#print(sample['fn'])
	
	target_params = Variable(sample['params'].float())#.cuda()
	input = Variable(sample['input'].float())#.cuda()
	weightMap_batch = sample['weight'].float()

	output = net.forward(input)
	sigmaT_pred_ = torch.clamp(output[:, 0] + output[:, 1], 0, 300)
	albedo_pred_ = torch.clamp(output[:, 0] / (output[:, 0] + output[:, 1]), 0.35, 0.97)
	g_pred_      = torch.clamp(1 - output[:, 2] / output[:, 0], -0.1, 0.9)

	#print(output)
	#sigmaT_pred_ = torch.clamp(torch.exp(output[:, 0]), 0, 300)
	#albedo_pred_ = torch.clamp(output[:, 1], 0.35, 0.97)
	#g_pred_      = torch.clamp(output[:, 2], -0.1, 0.9)
	#rest_pred_ = output[:, 1:]


	params = torch.Tensor(3)
	#print(target_params[0].data[0])
	params[0] =  np.log(sigmaT_pred_)
	params[1] =  albedo_pred_
	params[2] =  g_pred_

	shape = sample['shape'][0]
	light = sample['light'][0]
	xml_fn = shape + '_sunsky.xml'
	
	fi = int(light)
	xValue = math.sin(math.radians(theta))*math.cos(math.radians(fi))
	zValue = math.sin(math.radians(theta))*math.sin(math.radians(fi))
	yValue = math.cos(math.radians(theta))

	target_batch = sample['target'].float()
	params_new = Variable(params, requires_grad = True)
	optimizer = optim.Adam([params_new], lr = 0.01)
	allSigmaT = []
	allAlbedo = []
	allG = []
	allloss = []
	print(params)		
	for j in range(num_epochs):
	
		sigmaT_ac = np.clip(params_new.data[0], np.log(20.0),np.log(350.0)) #torch.clamp(params[0], 10.0, 300.0)
		albedo_ac = np.clip(params_new.data[1], 0.35, 0.97)  #torch.clamp(params[1], 0.05, 0.95)
		g_ac = np.clip(params_new.data[2], -0.1, 0.9)  #torch.clamp(params[2], -0.2, 0.9)
		optimizer.param_groups[0]['params'][0].data[0] = sigmaT_ac
		optimizer.param_groups[0]['params'][0].data[1] = albedo_ac
		optimizer.param_groups[0]['params'][0].data[2] = g_ac


		sigmT_afterexp = torch.exp(params_new[0])
		rest_pred = params_new[1:]
		adterEXP_pred_clamped = torch.cat((sigmT_afterexp, rest_pred))

		
		tmpf = shape
		cmd = 'cd ~/GfxMLForInverse/mitsuba_autodiff/dist && ./mitsuba_AD ~/GfxMLForInverse/GfxMLForInverse_datasets/scenesAD_sunsky_inv/'+xml_fn+' -Dmeshmodel='+shape +' -Dx='+str(xValue) +' -Dy='+str(yValue)+' -Dz='+str(zValue) + ' -DsigmaT='+str(adterEXP_pred_clamped.data[0])+' -Dalbedo='+str(adterEXP_pred_clamped.data[1])+' -Dg=' + str(adterEXP_pred_clamped.data[2])+' -DnumSamples='+str(nsamples)+' -q -b 4 -o ~/GfxMLForInverse/Invtmp/RGnaive/'+tmpf+'.exr'

		#print(cmd)
		os.system(cmd)

		tmp  = pyexr.open('tmp/'+tmpf+'.exr')
		batch_result_fwd = torch.from_numpy(tmp.get('Forward'))[:,:,0].unsqueeze(0)
		batch_result_sigmaT = torch.from_numpy(tmp.get('sigmaT'))[:,:,0].unsqueeze(0)
		batch_result_albedo = torch.from_numpy(tmp.get('albedo'))[:,:,0].unsqueeze(0)
		batch_result_hg = torch.from_numpy(tmp.get('g'))[:,:,0].unsqueeze(0)

		gt_img_var = Variable(target_batch[:, 0, :, :])
		wm_img_var = Variable(weightMap_batch[:, 0, :, :])
		tmp_img_var = Variable(batch_result_fwd)
		tmp_grad_sigmaT_var = Variable(batch_result_sigmaT)
		tmp_grad_albedo_var = Variable(batch_result_albedo)
		tmp_grad_hg_var = Variable(batch_result_hg)
		
		loss = criterion_MTS(adterEXP_pred_clamped, tmp_img_var, tmp_grad_sigmaT_var, tmp_grad_albedo_var, tmp_grad_hg_var, gt_img_var, wm_img_var)
		#print(loss)
	
		print('loss = ', loss.data[0][0], ' pred sigmaT = ', adterEXP_pred_clamped.data[0], ' albedo = ', adterEXP_pred_clamped.data[1], ' hg = ', adterEXP_pred_clamped.data[2])

		allSigmaT += [adterEXP_pred_clamped.data[0]]
		allAlbedo += [adterEXP_pred_clamped.data[1]]
		allG += [adterEXP_pred_clamped.data[2]]
		allloss += [loss.data[0][0]]
		optimizer.zero_grad()
		loss.backward()
		optimizer.step()
	totalloss += [allloss]
	totalSigmaT += [allSigmaT]
	totalAlbedo += [allAlbedo]
	totalG += [allG]
	print(totalloss)
	# np.save('inv/itn_similar_wm/SigmaT_001_scheduled.npy', np.array(totalSigmaT))
	# np.save('inv/itn_similar_wm/albedo_001_scheduled.npy', np.array(totalAlbedo))
	# np.save('inv/itn_similar_wm/g_001_scheduled.npy', np.array(totalG))
	# np.save('inv/itn_similar_wm/loss_001_scheduled.npy', np.array(totalloss)) = torch.cat((sigmaT_pred_.unsqueeze(1), albedo_pred_.unsqueeze(1), g_pred_.unsqueeze(1)), 1)






