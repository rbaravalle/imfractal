import os
import argparse
import subprocess

def main(args):

    print "Training........."

    subprocess.call("python test_sea.py -imgs "+args.train_images[0]+" -data "+args.data[0]+" -dfs "+str(args.dfs[0])+" -tr "+args.tr[0]+" -p "+str(args.ptrain[0])+" -eq "+str(args.equalize[0]), shell=True)

    print "Binarizing images......"

    subprocess.call("python binarize.py -ws "+str(args.ws[0])+" -pmodel "+args.model[0]+" -c "+args.classifier[0]+" -core "+str(args.core[0])+" -t "+args.training[0]+" -i "+args.images[0]+" -dfs "+str(args.dfs[0])+" -tr "+args.tr[0]+" -o "+str(args.output_folder[0])+" -eq "+str(args.equalize[0]), shell=True)


    print "Comparing with Ground Truth........"
    subprocess.call("python compare_with_gt.py -g "+args.gt_dir[0]+" -i "+args.output_folder[0], shell=True)
   


if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Compare a folder with GT')
    parser.add_argument("-ws", dest="ws", type=int, required=True, nargs='+', help="Sizes in pixels for sliding windows. Must be corresponding to the ws used in the classification model given")
    parser.add_argument("-pmodel", dest="model", type=str, required=True, nargs='+', help="Path to the models files .pkl")   #eg.: sample-point-sea,sample-point-sky,sample-point-dolphin             
    parser.add_argument("-c", dest="classifier", type=str, required=True, nargs='+', help="Type of classifier: RF (for Random Forest) or SVC (for Support Vector Classification)")
    parser.add_argument("-t", dest="training", type=str, required=True, nargs='+', help="Type of classification training: tt (training and test) or cv (cross validation)")
    parser.add_argument("-i", dest="images", type=str, required=True, nargs='+', help=".jpg files image to be binarized")
    parser.add_argument("-core", dest="core", type=float, required=True, nargs='+', help="Size in pixels or percentage expressed as decimal for the sliding window core where is putted the result of the classification.")
    parser.add_argument("-dfs", dest="dfs", type=int, required=True, nargs=1, help="Amount of fractal dimensions in the MFS")
    parser.add_argument("-tr", dest="tr", type=str, required=True, nargs=1, help="Transform or not")

    parser.add_argument("-data", dest="data", type=str, required=True, nargs=1, help="Path where data will be saved")

    parser.add_argument("-ptrain", dest="ptrain", type=float, required=True, nargs=1, help="Percentage of dataset to be used for training (float)")

    parser.add_argument("-timgs", dest="train_images", type=str, required=True, nargs=1, help="Path with images for training")
    parser.add_argument("-g", dest="gt_dir", type=str, required=True, nargs=1, help="Dir with GT binarized images")
    parser.add_argument("-o", dest="output_folder", type=str, required=True, nargs=1, help="Output folder for binarized images")
    parser.add_argument("-eq", dest="equalize", type=str, required=True, nargs=1, help="Equalize image or not for training")



    args = parser.parse_args()

    main(args)
