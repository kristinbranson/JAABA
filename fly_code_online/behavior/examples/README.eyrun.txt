
**Sample train command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/synth_flies/Params/LungeWalkSVMBehaviorParams.txt -F Data/synth_flies/Params/LungeWalkSVMFeatureParams.txt -d Data/synth_flies/TrainingData/traindata.txt -o Data/synth_flies/Model/LungeWalkSVMStructLearntModel.txt

NOTE1: Run with ADD_DUMMY_BEHAVIORS=0 and ALLOW_SAME_TRANSITIONS=1.
NOTE2: This has to be run with DEBUG=0. This is because in debug mode there are checks to see whether the predicted label has a higher score than the ground truth label. Since this version allows for self-transitions and only goes 50 frames back, a 500 frame long segment of type "other" would have a score of say 20, but ten 50 frame long segments of type "other" would have a score of 200. I wanted to fix this by weighing the bout score with the number of frames it spans, but that currently causes a bug in the code. I will be looking into that. The predicted output looks good though.

**Sample test command:
bin/debug_static/svm_fly_behavior_sequence.out -B Data/synth_flies/Params/LungeWalkSVMBehaviorParams.txt -F Data/synth_flies/Params/LungeWalkSVMFeatureParams.txt -t Data/synth_flies/TrainingData/testdata.txt Data/synth_flies/TrainingData/testdata_prediction.txt -i Data/synth_flies/Model/LungeWalkSVMStructLearntModel.txt

**To compare ground truth to predicted output:
In Matlab, run compare_pred_gt.m
