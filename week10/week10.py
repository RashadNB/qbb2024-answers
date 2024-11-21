#!/usr/bin/env python

import numpy as np
import scipy
import matplotlib.pyplot as plt
import imageio
import pandas as pd
import plotly.express as px
import plotly

# Set initial constants
genes = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
fields = ["field0", "field1"]
channels = ["DAPI", "nascentRNA", "PCNA"]
min_pixels = 100

# Load images into RGB format
image_arrays = []
for gene in genes:
    for field in fields:
        img_array = None
        for i, channel in enumerate(channels):
            filename = f"{gene}_{field}_{channel}.tif"
            img = imageio.v3.imread(filename).astype(np.uint16)
            if img_array is None:
                img_array = np.zeros((img.shape[0], img.shape[1], 3), dtype=np.uint16)
            img_array[:, :, i] = img
        image_arrays.append(img_array)

# Function to find labels in binary mask
def find_labels(mask):
    l = 0
    labels = np.zeros(mask.shape, np.int32)
    equivalence = [0]

    if mask[0, 0]:
        l += 1
        equivalence.append(l)
        labels[0, 0] = l

    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                l += 1
                equivalence.append(l)
                labels[0, y] = l

    for x in range(1, mask.shape[0]):
        if mask[x, 0]:
            if mask[x - 1, 0]:
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                l += 1
                equivalence.append(l)
                labels[x, 0] = l

        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    if mask[x - 1, y - 1]:
                        labels[x, y] = min(equivalence[labels[x - 1, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x - 1, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    elif mask[x, y - 1]:
                        labels[x, y] = min(equivalence[labels[x, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    else:
                        labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l

        if mask[x, -1]:
            if mask[x - 1, -1]:
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                labels[x, -1] = equivalence[labels[x, -2]]
            else:
                l += 1
                equivalence.append(l)
                labels[x, -1] = l

    equivalence = np.array(equivalence)
    for i in range(1, len(equivalence))[::-1]:
        labels[np.where(labels == i)] = equivalence[i]
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        labels[np.where(labels == j)] = i
    return labels

# Function to filter labels by size
def filter_by_size(labels, minsize, maxsize):
    sizes = np.bincount(labels.ravel())
    for i in range(1, sizes.shape[0]):
        if sizes[i] < minsize or sizes[i] > maxsize:
            where = np.where(labels == i)
            labels[where] = 0
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        labels[np.where(labels == j)] = i
    return labels

# Process images to create masks and extract data
data = []
for img_idx, image in enumerate(image_arrays):
    dapi_channel = image[:, :, 0]
    mask = dapi_channel > np.mean(dapi_channel)
    labels = find_labels(mask)
    labels = filter_by_size(labels, min_pixels, 1000000)

    sizes = np.bincount(labels.ravel())[1:]
    mean_size = np.mean(sizes)
    std_size = np.std(sizes)
    labels = filter_by_size(labels, mean_size - std_size, mean_size + std_size)

    num_nuclei = np.max(labels) + 1
    for nuc in range(1, num_nuclei):
        where = np.where(labels == nuc)
        pcna_signal = np.mean(image[where][1])
        nrna_signal = np.mean(image[where][2])
        log2_ratio = np.log2(nrna_signal / pcna_signal)
        data.append({
            "gene": genes[img_idx // 2],
            "field": fields[img_idx % 2],
            "nucleiNumber": nuc,
            "nascentRNA": nrna_signal,
            "PCNA": pcna_signal,
            "log2_ratio": log2_ratio
        })

# Save results to CSV
pd.DataFrame(data).to_csv("Nuclei.csv", sep=",", header=True, index=False, mode="w")


