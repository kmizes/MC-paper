Raw kinematic data for "The role of motor cortex in motor sequence execution depends on demands for flexibility"


Data is preprocessed and saved in .mat files, as vectors of trials.

Summary of variables:
 - joint: label of the tracked body part
 - lever: sequence performed
 - seqchange: indicator of lesion condition. First digit is either mock break or a placeholder, second indicates whhen the unilateral lesion occured, second is the bilateral.
 - sessID: index of which session the trial came from, for behavioral structures.
 - tSess: matlab datetime of when the trial session occured.
 - tapon: tap x trial vector of the frame when the lever press occured
 - traj: frame x 2 vector of the joint position in pixels. all videos were captured at 40hz framerate
 - traj_prob: frame x 1 vector of the confidence of the tracking, provided from DeeperCut. 