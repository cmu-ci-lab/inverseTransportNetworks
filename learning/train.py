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

from sklearn.feature_extraction.image import extract_patches
angle_train = '0_30_60_90_train_train'


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
			#x, self.pool_idx8 = F.max_pool2d(F.relu(self.conv8(x)), kernel_size=2, return_indices=True)
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
	def __init__(self, csv_file, root_dir, info_dir, env_dir, wm_dir, transform=None):
				self.data_frame = pd.read_csv(csv_file)
				self.root_dir   = root_dir				
				self.transform  = transform
				self.wm_dir    = wm_dir
	def __len__(self):
		return len(self.data_frame)
	def __getitem__(self, idx):
  		
			input     = np.zeros(shape = (128, 128, 2))

			img_fn    = os.path.join(self.root_dir, self.data_frame.ix[idx, 0])
			luminance = pyexr.open(img_fn).get() # luminance 
			wm_fn     = os.path.join(self.wm_dir, self.data_frame.ix[idx,0])
			wm        = pyexr.open(wm_fn).get()

			shape 	  = self.data_frame.ix[idx, 1]
			light 	  = self.data_frame.ix[idx, 2]

			sigmaT 	  = self.data_frame.ix[idx, 3].astype('float')
			albedo    = self.data_frame.ix[idx, 4].astype('float')
			hg        = self.data_frame.ix[idx, 5].astype('float') 

			sigmaS = sigmaT * albedo
			sigmaA = sigmaT - sigmaS
			sigmaS_reduced = (1 - hg) * sigmaS

			patches_luminance = extract_patches(luminance[:,:,0], patch_shape=(4, 4), extraction_step=(4, 4))
			input[:,:,0] = patches_luminance.mean(-1).mean(-1)
			patches_wm = extract_patches(wm[:,:,0], patch_shape=(4, 4), extraction_step=(4, 4))
			input[:,:,1] = patches_wm.mean(-1).mean(-1)

			target = torch.from_numpy(patches_luminance.mean(-1).mean(-1)).unsqueeze(0)
			weight = torch.from_numpy(patches_wm.mean(-1).mean(-1)).unsqueeze(0)
			input = input.transpose((2, 0, 1))
			input = torch.from_numpy(input)

			params    = torch.Tensor(3)
			params[0] = sigmaS
			params[1] = sigmaA 
			params[2] = sigmaS_reduced 

			sample = {
					'input' : input,
					'params': params,
					'target' : target,
					'shape' : shape,
					'light' : light,
					'sigmaS' : sigmaS,
					'sigmaA' : sigmaA,
					'sigmaS_reduced' : sigmaS_reduced,
					'fn' : self.data_frame.ix[idx, 0],
					'weightMap': weight
			}

			return sample

 

five_shape_ten_light_train_dataset = AllShapeDataset(
																		csv_file='/shared/e0_30_60_90_train_train.csv',
																		root_dir='/shared/dataset_s2b_blur_train_mean1_nb/',
																		wm_dir  = '/shared/dataset_s2b_blur_train_wm/'
																)

net = MatNet()
lr_rate =  0.0001
optimizer = optim.Adam(net.parameters(), lr=lr_rate)
# can also initialize with a pre-trained regressor


#preTrain = torch.load('/phoenix/S6/chec/results/data_allshregressor9.bin')
#net = preTrain['model']
#optimizer = preTrain['optimizer']



class MitsubaLossFunction(torch.autograd.Function):

	normalization = 28.2486

	def __init__(self):
		super(MitsubaLossFunction, self).__init__()

	def forward(self, input, tmp_img, tmp_grad_sigmaT, tmp_grad_albedo, tmp_grad_hg, gt_img, wm_img):
 
		diff = self.normalization * tmp_img - gt_img
		self.save_for_backward(input, tmp_img.float(), tmp_grad_sigmaT.float(), tmp_grad_albedo.float(), tmp_grad_hg.float(), gt_img.float(), wm_img.float())
		loss = torch.Tensor(1, 1)
		wm_square_loss = torch.pow(diff, 2)*wm_img
		sums = wm_square_loss.sum(2).sum(1)
		sums_wm = wm_img.sum(2).sum(1)
		loss[0][0] = (sums/sums_wm).mean()
		
		return loss

	def backward(self, grad_out):
		input, tmp_img,tmp_grad_sigmaT, tmp_grad_albedo, tmp_grad_hg, gt_img, wm_img = self.saved_tensors
		diff = self.normalization*tmp_img - gt_img
		grad_in = torch.Tensor(input.size()).float()
		tmp = wm_img.sum(2).sum(1).unsqueeze(1).unsqueeze(2)
		wm_by_const = 2*wm_img / (tmp.repeat(1, 1, 128, 128).squeeze(0))

		grad_in[:, 0] = (diff * self.normalization * tmp_grad_sigmaT * wm_by_const).mean(2).mean(1) * grad_out.clone()[0][0]
		grad_in[:, 1] = (diff * self.normalization * tmp_grad_albedo * wm_by_const).mean(2).mean(1) * grad_out.clone()[0][0]
		grad_in[:, 2] = (diff * self.normalization * tmp_grad_hg * wm_by_const).mean(2).mean(1) * grad_out.clone()[0][0]

		
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


batch_size = 60
num_epochs = 50

nsamples = 256	
theta = 30
data_train_loader = DataLoader(five_shape_ten_light_train_dataset, 
							   batch_size=batch_size, 
							   shuffle=True, 
							   num_workers=16)

momentum= 0.9

criterion_MTS = MitsubaLossModule()
criterion_Para = nn.MSELoss()



all_training_results = []
all_training_results_para = []
all_training_results_mts = []


min_testing_loss = 1e10
IPs = np.loadtxt('ips.txt', dtype = str)
ip_add = IPs.tolist()


for epoch in range(num_epochs):
	epoch_str = 'Epoch is ' + str(epoch)
	cur_results = []
	cur_results_para = []
	cur_results_mts = []		
	print('Running epoch i = ', epoch)
	epoch_loss = 0
	cnt = 0

	for i_batch, sample in enumerate(data_train_loader):

		print('batch strats now----------',i_batch,epoch)

		batch_str = 'batch is ' + str(i_batch)

		input_batch = sample['input'].float()
		params_batch = sample['params'].float()
		target_batch = sample['target'].float()
		weightMap_batch = sample['weightMap'].float()

		shape_batch = sample['shape']
		light_batch = sample['light']
		fn_batch = sample['fn']

		target_params = Variable(sample['params'].float())
		input = Variable(sample['input'].float())
		output = net.forward(input)
		loss_Para = criterion_Para(output, target_params)

		print('PARAM [GT, PD]:\n sigmaS = [{:.4f},{:.4f}]'.format(target_params.data[0][0], output.data[0][0]), 
							 '\n sigmaA = [{:.4f},{:.4f}]'.format(target_params.data[0][1], output.data[0][1]),  
							 '\n sigmaS_red = [{:.4f},{:.4f}]'.format(target_params.data[0][2], output.data[0][2]))

		
		sigmaS_pred = output[:,0]
		sigmaA_pred = output[:,1]
		sigmaS_r_pred = output[:,2]
		
		sigmaT_pred = sigmaS_pred + sigmaA_pred
		albedo_pred = sigmaS_pred / sigmaT_pred
		hg_pred = 1 - (sigmaS_r_pred / sigmaS_pred)

		sigmaT_pred_clamped = torch.clamp(sigmaT_pred, 20.0, 300.0)
		albedo_pred_clamped = torch.clamp(albedo_pred, 0.35, 0.97)
		hg_pred_clamped = torch.clamp(hg_pred, -0.1, 0.9)


		output_MTSready = torch.cat((sigmaT_pred_clamped.unsqueeze(1), albedo_pred_clamped.unsqueeze(1), hg_pred_clamped.unsqueeze(1)), 1)

		os.system('rm -r -f tmp/*.exr')
		for i in range(len(input_batch)):
			cur_paras = params_batch[i]
			cur_output = output[i]

			sigmaT_pred = output_MTSready[i, 0].data.cpu().double()
			albedo_pred = output_MTSready[i, 1].data.cpu().double()
			hg_pred = output_MTSready[i, 2].data.cpu().double()

			shape = shape_batch[i]
			light = light_batch[i]
			xml_fn = shape + '_sunsky.xml' 

			fi = int(light)
			xValue = math.sin(math.radians(theta))*math.cos(math.radians(fi))
			zValue = math.sin(math.radians(theta))*math.sin(math.radians(fi))
			yValue = math.cos(math.radians(theta))

			cmd = 'nohup ssh '+ip_add[i]+' "cd ~/mitsuba/dist && ./mitsuba_AD ~/GfxMLForInverse/GfxMLForInverse_datasets_full/scenesAD_sunsky/'+xml_fn+' -Dmeshmodel='+shape +' -Dx='+str(xValue) +' -Dy='+str(yValue)+' -Dz='+str(zValue) + ' -DsigmaT='+str(sigmaT_pred[0])+' -Dalbedo='+str(albedo_pred[0])+' -Dg=' + str(hg_pred[0])+' -DnumSamples='+str(nsamples)+' -b 4 -p 18 -o ~/GfxMLForInverse/tmp/tmp'+str(i)+'.exr" > /shared/outputs/output'+str(i)+'.out 2> /shared/errors/error'+str(i)+'.txt &'
			
			os.system(cmd)
	
		print('------------------------current input length------------------',len(input_batch))

		while(len(os.listdir('tmp/'))!=len(input_batch)):
			pass

		flag = 0
		already_init = False
		batch_result_fwd = None
		batch_result_sigmaT = None
		batch_result_albedo = None
		batch_result_hg = None
	

		while(flag != len(input_batch)):
			fileName = 'tmp/tmp' + str(flag) + '.exr'
			tmp = loadFigure(fileName)
			if (tmp != None):
				flag += 1


				if (already_init == False):
					batch_result_fwd = torch.from_numpy(tmp.get('Forward'))[:,:,0].unsqueeze(0)
					batch_result_sigmaT = torch.from_numpy(tmp.get('sigmaT'))[:,:,0].unsqueeze(0)
					batch_result_albedo = torch.from_numpy(tmp.get('albedo'))[:,:,0].unsqueeze(0)
					batch_result_hg = torch.from_numpy(tmp.get('g'))[:,:,0].unsqueeze(0)
					already_init = True
				else:
					batch_result_fwd = torch.cat((batch_result_fwd, torch.from_numpy(tmp.get('Forward'))[:,:,0].unsqueeze(0)), 0)
					batch_result_sigmaT = torch.cat((batch_result_sigmaT, torch.from_numpy(tmp.get('sigmaT'))[:,:,0].unsqueeze(0)), 0)
					batch_result_albedo = torch.cat((batch_result_albedo, torch.from_numpy(tmp.get('albedo'))[:,:,0].unsqueeze(0)), 0)
					batch_result_hg = torch.cat((batch_result_hg, torch.from_numpy(tmp.get('g'))[:,:,0].unsqueeze(0)), 0)
												
			else:
				print('wait')
				time.sleep(1)

		gt_img_var = Variable(target_batch[:, 0, :, :])
		wm_img_var = Variable(weightMap_batch[:, 0, :, :])
		tmp_img_var = Variable(batch_result_fwd)
		tmp_grad_sigmaT_var = Variable(batch_result_sigmaT)
		tmp_grad_albedo_var = Variable(batch_result_albedo)
		tmp_grad_hg_var = Variable(batch_result_hg)
		
		lossMTS = criterion_MTS(output_MTSready, tmp_img_var, tmp_grad_sigmaT_var, tmp_grad_albedo_var, tmp_grad_hg_var, gt_img_var, wm_img_var)
		
		loss = 0.5 * (1.5 * lossMTS + loss_Para)  
		print('mts loss = ', 1.5 * lossMTS.data[0][0])
		print('param loss = ', loss_Para.data[0])
		cur_results += [loss.data[0]]
		cur_results_para += [loss_Para.data[0]]
		cur_results_mts += [1.5 * lossMTS.data[0][0]]

		epoch_loss += loss.data[0]
		optimizer.zero_grad()
		loss.backward()
		optimizer.step()
		cnt += 1

	all_training_results += [cur_results]
	all_training_results_para += [cur_results_para]
	all_training_results_mts += [cur_results_mts]

	np.save('/shared/results/training_error.npy', np.array(all_training_results))
	np.save('/shared/results/training_error_para.npy', np.array(all_training_results_para))
	np.save('/shared/results_/training_error_mts.npy', np.array(all_training_results_mts))
	epoch_loss /= cnt
	netToSave = {'model': net, 'optimizer': optimizer}
	torch.save(netToSave, '/nets/ITN_'+str(epoch)+'.bin')
	print('Epoch i =', epoch, 'done. Training Loss =', epoch_loss)	

