"""
Module: src/config.py
Project: ISdetector
Version: 1.1.0
Author: Yang Zhou
Purpose: Store default values for parameters and flags used across pipeline stages.
"""

# Stage 1 & 3: 
MIN_ANCHOR_LEN = 20

# Stage 2: Reference IS removal (BLAST filters; CLI --identity / --fraction)
BLAST_IDENTITY = 95.0   # Minimum percent identity (0-100)
BLAST_FRACTION = 0.90   # Minimum aligned length / IS length (0-1)

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
