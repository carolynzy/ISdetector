"""
Module: src/config.py
Project: ISdetector v1.0
Author: Yang Zhou
Purpose: Store default values for parameters and flags used across pipeline stages.
"""

# Stage 1 & 3: 
MIN_ANCHOR_LEN = 20

# Stage 3: Clustering
CLUSTER_DIST = 50
MIN_CLUSTER_SUPPORT = 20
MIN_PEAK_SUPPORT = 5
PEAK_DISTANCE = 30
DEPTH_DIFF_THRESHOLD = 0.3
ZERO_THRESH = 3 

# Stage 4: Detection
PAIR_GAP = 20
MAX_SV_SIZE = 30000
