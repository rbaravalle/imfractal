#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 12:12:40 2017

@author: deby
"""
import argparse
from sklearn.externals import joblib

from imfractal import MFS
from PIL import Image
import time
import csv
import sys
import os
from subprocess import *

from sklearn.ensemble import RandomForestClassifier
from sklearn import cross_validation
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib

import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import misc
from scipy.stats import mode


sea_label=1
dolphin_label=2
true_values = ['t', 'T', '1', 1, 'true', 'True', 'TRUE']

parser = argparse.ArgumentParser(description='Binarize an image using classifiers models')
parser.add_argument("-ws", dest="windowSizes", type=int, required=True, nargs='+', help="Sizes in pixels for sliding windows. Must be corresponding to the ws used in the classification model given")
parser.add_argument("-pmodel", dest="path_models", type=str, required=True, nargs='+', help="Path to the models files .pkl")   #eg.: sample-point-sea,sample-point-sky,sample-point-dolphin             
parser.add_argument("-c", dest="classifier", type=str, required=True, nargs='+', help="Type of classifier: RF (for Random Forest) or SVC (for Support Vector Classification)")
parser.add_argument("-t", dest="training", type=str, required=True, nargs='+', help="Type of classification training: tt (training and test) or cv (cross validation)")
parser.add_argument("-i", dest="images", type=str, required=True, nargs='+', help=".jpg files image to be binarized")
parser.add_argument("-core", dest="core", type=float, required=True, nargs='+', help="Size in pixels or percentage expressed as decimal for the sliding window core where is putted the result of the classification.")
parser.add_argument("-dfs", dest="dfs", type=int, required=True, nargs=1, help="Amount of fractal dimensions in the MFS")
parser.add_argument("-yiq", dest="yiq", type=str, required=True, nargs=1, help="Transform to YIQ or not")

args = parser.parse_args()

yiq_str = "yiq" if args.yiq[0] in true_values else "no_yiq"
dfs_str = str(args.dfs[0])


def sliding_window(image, stepSize, windowSize):
    # slide a window across the image
    high, width = image.shape
    for y in xrange(0, high, stepSize): #filas
        for x in xrange(0, width, stepSize): #columnas
            # yield the current window           
            yield (x, y, image[y:y + windowSize[0], x:x + windowSize[1]])


def binarize(filename_img_to_binarize, filename_model, windowSizes, coreSize):

    #levanto imagen a binarizar
    img = []
    if filename_img_to_binarize != '':
        img = Image.open(filename_img_to_binarize)
        # Preprocessing: if IM is a color image convert it to a gray image
        img = img.convert("L")        
        img_np= np.array(img)        

    else:
         print "Must specify a filename_img_to_binarize"
         exit()    

    high, width = img_np.shape
    
    #creo nueva img del tamano de la original para el resultado de la binarizacion       
    # con cantidad de dimensiones igual a  la cantidad de tamanos de las ventanas 
    #en cada capa del axis 2 (dim), guardo los resultados de cada tamano de ventana, 
    #en la capa adicional guardo el resultado del voting de las capas de cada tamanio de ventana
    binarize_img=np.zeros((len(windowSizes) + 1,high,width),dtype=np.int)

    #seteo MFS igual al que use cuando entrené los clasificadores
    ins = MFS()
    ins.setDef(1,args.dfs[0],3)
    
    #string con los tamaños de ventanas usadas. lo uso para concatenarlo al nombre del archivo de imagen "binarizado"
    wsrealized=""

    
    for s in range(len(windowSizes)):
        
        winSize = windowSizes[s]
        wsrealized+="-" + str(winSize)    
        (winH,winW) = (winSize,winSize)

        #levanto modelo del clasificador segun tamano de grilla 
        if filename_model != '':
            filename_model_size=filename_model.format(winSize, yiq_str, dfs_str)
            clf = joblib.load(filename_model_size)
        else:
            print "Must specify a filename_model"
            exit()
        
        # si es un porcentaje de la ventana
        if (coreSize > 0 and coreSize < 1):
            core_size = int (coreSize*winSize)
        else:
            core_size= int (coreSize)
        
        if(core_size%2!=0):
                core_size=core_size+1            

        stepSize=core_size

        # loop over the sliding window for each layer of the pyramid
        for (col, row, window) in sliding_window(img_np, stepSize, windowSize=(winH,winW)):


        # if the window does not meet our desired window size, ignore it
            if window.shape[0] != winH or window.shape[1] != winW:
                continue
            mfs_window = ins.getFDs("",window)
            mfs_window = mfs_window.reshape(1, mfs_window.shape[0])
            predicted_value= clf.predict(mfs_window)
#            print "row {} , col {} ".format(row,col)

#            print predicted_value.astype(np.int64)[0]           
            row_win_center = row + (winSize/2)
            col_win_center = col + (winSize/2)
            
            row_from = row_win_center - (core_size/2)
            row_to = row_win_center + (core_size/2)
            
            col_from = col_win_center - (core_size/2)
            col_to = col_win_center + (core_size/2)
            
#            print "row: {} : {} - col {}: {}".format(row_from,row_to,col_from,col_to)
            
            binarize_img[s,row_from:row_to,col_from:col_to]= predicted_value.astype(np.int64)[0]

#            print "binarize_img [{},{}:{},{}:{}]".format(s,row_from,row_to,col_from,col_to)
        # los pixeles quedarian 0=si no tiene etiqueta 1=etiqueta1 2=etiqueta2 and so on
    
    # aplico mode(voting) sobre cada pixel en todas las capas
    #guardo el resultado del voting en la capa correspondiente de la imgen binarizada binarize_img
    for row in range(high):
        for column in range(width):
            #en la ultima capa de binarize_img guardo el voting pixel a pixel de todas las capas de cada winSize
      
            label_voted =mode(binarize_img[:-1,row,column])
            binarize_img[-1,row,column]= label_voted[0]
   
    # aplico mascara sobre la imagen original para segmentar  (pinto de blanco segun la etiqueta) 
        
    dolphin_masked_img = img_np.copy()
    dolphin_mask= binarize_img[-1,:,:] == dolphin_label
    dolphin_masked_img [dolphin_mask] = 255
    
    sea_masked_img = img_np.copy()    
    sea_mask= binarize_img[-1,:,:] == sea_label   
    sea_masked_img [sea_mask] = 255

#    f, (ax0, ax1, ax2) = plt.subplots(1, 2, 3)
#    ax0.imshow(img_np, cmap='gray')
#    ax1.imshow(dolphin_masked_img, cmap='gray')
#    ax1.imshow(sea_masked_img, cmap='gray')
#    plt.show()

    #path
    outdir=os.path.dirname(filename_img_to_binarize)
    #make new file name to save the binarized images
    orginalfn=os.path.basename(filename_img_to_binarize)
            
    binarized_dolphinfn="{}_dolphin_{}_{}_n{}{}_{}_{}{}".format(os.path.splitext(orginalfn)[0],args.classifier[0],args.training[0], str(coreSize),wsrealized,yiq_str,dfs_str,os.path.splitext(orginalfn)[1])
    outdirfn="{}/{}".format(outdir,binarized_dolphinfn)
    
    im_binarized= Image.fromarray(dolphin_masked_img)
    im_binarized.save(outdirfn)
    print 'saving '+ outdirfn
    
    
    binarized_seafn="{}_sea_{}_{}_n{}{}_{}_{}{}".format(os.path.splitext(orginalfn)[0],args.classifier[0],args.training[0], str(coreSize),wsrealized,yiq_str,dfs_str,os.path.splitext(orginalfn)[1])
    outdirfn="{}/{}".format(outdir,binarized_seafn)
    
    im_binarized= Image.fromarray(sea_masked_img)
    im_binarized.save(outdirfn)    
    print 'saving '+ outdirfn

def do_test():
   
    print args.path_models[0]
    print args.classifier[0]
    print args.training[0]
    
    filename_model = args.path_models[0] + "/"+ args.classifier[0] +"_"+ args.training[0] +"_{}_{}_{}.pkl"
        
    for filename_img_to_binarize in args.images:
            binarize(filename_img_to_binarize, filename_model, args.windowSizes, args.core[0])            

do_test()
